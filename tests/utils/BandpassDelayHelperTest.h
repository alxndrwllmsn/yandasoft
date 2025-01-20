/// @file
///
/// Unit test for BandpassDelayHelper class which is intended to solve for and modify delays in bandpass solutions
///
///
/// @copyright (c) 2025 CSIRO
/// Australia Telescope National Facility (ATNF)
/// Commonwealth Scientific and Industrial Research Organisation (CSIRO)
/// PO Box 76, Epping NSW 1710, Australia
/// atnf-enquiries@csiro.au
///
/// This file is part of the ASKAP software distribution.
///
/// The ASKAP software distribution is free software: you can redistribute it
/// and/or modify it under the terms of the GNU General Public License as
/// published by the Free Software Foundation; either version 2 of the License,
/// or (at your option) any later version.
///
/// This program is distributed in the hope that it will be useful,
/// but WITHOUT ANY WARRANTY; without even the implied warranty of
/// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/// GNU General Public License for more details.
///
/// You should have received a copy of the GNU General Public License
/// along with this program; if not, write to the Free Software
/// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
///
/// @author Max Voronkov <maxim.voronkov@csiro.au>

#include <askap/utils/BandpassDelayHelper.h>
#include <askap/calibaccess/CachedCalSolutionAccessor.h>
#include <askap/scimath/fitting/CalParamNameHelper.h>
#include <askap/scimath/fitting/Params.h>
#include <askap/askap/AskapUtil.h>

#include <cppunit/extensions/HelperMacros.h>

namespace askap {

namespace utils {

class BandpassDelayHelperTest : public CppUnit::TestFixture 
{
   CPPUNIT_TEST_SUITE(BandpassDelayHelperTest);
   CPPUNIT_TEST(testLoad);
   CPPUNIT_TEST(testStore);
   CPPUNIT_TEST(testIdeal);
   CPPUNIT_TEST(testDelayCalcAndApply);
   CPPUNIT_TEST_SUITE_END();
public:

   void testLoad() {
      const casacore::uInt nAnt = 6u;
      const casacore::uInt nBeam = 9u;
      const casacore::uInt nChan = 16u;
      BandpassDelayHelper bdh(nAnt, nBeam, nChan);

      // instead of table-based solution accessor we use CachedCalSolutionAccessor adapter which stores info in Params class
      // first fill it with data and then try to load it into the helper class.
      scimath::Params buffer;
      for (casacore::uInt ant = 0; ant < nAnt; ++ant) {
           for (casacore::uInt beam = 0; beam < nBeam; ++beam) {
                // we'll have only one polarisation defined for testing depending on whether ant == beam or otherwise
                const std::string parName = scimath::CalParamNameHelper::paramName(ant, beam, ant == beam ? casacore::Stokes::YY : casacore::Stokes::XX);
                for (casacore::uInt chan = 0; chan < nChan; ++chan) {
                     // ensure we have different parameters for each antenna / beam / channel combination
                     buffer.add(scimath::CalParamNameHelper::addChannelInfo(parName, chan), casacore::Complex(ant + 1, beam + 1) * (static_cast<int>(chan) + 1));
                }
           }
      }

      // setup params-based solution accessor
      boost::shared_ptr<scimath::Params> bufferPtr(&buffer, utility::NullDeleter());
      accessors::CachedCalSolutionAccessor acc(bufferPtr);

      // load bandpass into the helper class
      bdh.loadBandpass(acc, 1e6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1e6, bdh.itsResolution, 1.);

      // now check the content of the helper class to verify that bandpass has been loaded correctly along with the validity flags
      for (casacore::uInt ant = 0; ant < nAnt; ++ant) {
           for (casacore::uInt beam = 0; beam < nBeam; ++beam) {
                for (casacore::uInt chan = 0; chan < nChan; ++chan) {
                     const casacore::IPosition validPoint(4, chan, beam == ant ? 1 : 0, beam, ant);
                     const casacore::IPosition invalidPoint(4, chan, beam == ant ? 0 : 1, beam, ant);
                     CPPUNIT_ASSERT(bdh.itsBandpassValid(validPoint));
                     CPPUNIT_ASSERT(!bdh.itsBandpassValid(invalidPoint));
                     const casacore::Complex expectedVal = casacore::Complex(ant + 1, beam + 1) * (static_cast<int>(chan) + 1);
                     CPPUNIT_ASSERT_DOUBLES_EQUAL(casacore::real(expectedVal), casacore::real(bdh.itsBandpass(validPoint)), 1e-6);
                     CPPUNIT_ASSERT_DOUBLES_EQUAL(casacore::imag(expectedVal), casacore::imag(bdh.itsBandpass(validPoint)), 1e-6);
                }
           }
      }
   }

   void testStore() {
      const casacore::uInt nAnt = 6u;
      const casacore::uInt nBeam = 9u;
      const casacore::uInt nChan = 16u;
      BandpassDelayHelper bdh(nAnt, nBeam, nChan);

      // fill the values in the buffer directly
      for (casacore::uInt ant = 0; ant < nAnt; ++ant) {
           for (casacore::uInt beam = 0; beam < nBeam; ++beam) {
                for (casacore::uInt chan = 0; chan < nChan; ++chan) {
                     const casacore::IPosition validPoint(4, chan, beam == ant ? 1 : 0, beam, ant);
                     const casacore::IPosition invalidPoint(4, chan, beam == ant ? 0 : 1, beam, ant);
                     bdh.itsBandpassValid(validPoint) = true;
                     // check that invalid point is invalid by default instead of assigning
                     CPPUNIT_ASSERT(!bdh.itsBandpassValid(invalidPoint));
                     bdh.itsBandpass(validPoint) = casacore::Complex(ant + 1, beam + 1) * (static_cast<int>(chan) + 1);
                }
           }
      }
      // params-based accessor to store the bandpass, Params object is held internally
      accessors::CachedCalSolutionAccessor acc;
      bdh.storeBandpass(acc);

      // check the results
      for (casacore::uInt ant = 0; ant < nAnt; ++ant) {
           for (casacore::uInt beam = 0; beam < nBeam; ++beam) {
                // we'll have only one polarisation defined for testing depending on whether ant == beam or otherwise
                const std::string parName = scimath::CalParamNameHelper::paramName(ant, beam, ant == beam ? casacore::Stokes::YY : casacore::Stokes::XX);
                const std::string otherPolParName = scimath::CalParamNameHelper::paramName(ant, beam, ant == beam ? casacore::Stokes::XX : casacore::Stokes::YY);
                // the actual parameter names have channel number encoded, so these ones should be missing
                CPPUNIT_ASSERT(!acc.cache().has(otherPolParName));
                CPPUNIT_ASSERT(!acc.cache().has(parName));
                for (casacore::uInt chan = 0; chan < nChan; ++chan) {
                     const std::string fullParName = scimath::CalParamNameHelper::addChannelInfo(parName, chan);
                     CPPUNIT_ASSERT(acc.cache().has(fullParName));
                     CPPUNIT_ASSERT(!acc.cache().has(scimath::CalParamNameHelper::addChannelInfo(otherPolParName, chan)));
                     const casacore::Complex expected = casacore::Complex(ant + 1, beam + 1) * (static_cast<int>(chan) + 1);
                     const casacore::Complex obtained = acc.cache().complexValue(scimath::CalParamNameHelper::addChannelInfo(parName, chan));
                     CPPUNIT_ASSERT_DOUBLES_EQUAL(casacore::real(expected), casacore::real(obtained), 1e-6);
                     CPPUNIT_ASSERT_DOUBLES_EQUAL(casacore::imag(expected), casacore::imag(obtained), 1e-6);
                }
           }
      }
   }

   void testIdeal() {
      const casacore::uInt nAnt = 6u;
      const casacore::uInt nBeam = 3u;
      const casacore::uInt nChan = 10u;
      BandpassDelayHelper bdh(nAnt, nBeam, nChan);

      bdh.setIdealBandpass(1e6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1e6, bdh.itsResolution, 1.);

      // check the results
      for (casacore::uInt ant = 0; ant < nAnt; ++ant) {
           for (casacore::uInt beam = 0; beam < nBeam; ++beam) {
                for (casacore::uInt pol = 0; pol < 2u; ++pol) {
                     for (casacore::uInt chan = 0; chan < nChan; ++chan) {
                          const casacore::IPosition index(4, chan, pol, beam, ant);
                          CPPUNIT_ASSERT(bdh.itsBandpassValid(index));
                          CPPUNIT_ASSERT_DOUBLES_EQUAL(1., casacore::real(bdh.itsBandpass(index)), 1e-6);
                          CPPUNIT_ASSERT_DOUBLES_EQUAL(0., casacore::imag(bdh.itsBandpass(index)), 1e-6);
                     }
                }
           }
      }
   }

   void testDelayCalcAndApply() {
      const casacore::uInt nAnt = 6u;
      const casacore::uInt nBeam = 9u;
      const casacore::uInt nChan = 64u;
      BandpassDelayHelper bdh(nAnt, nBeam, nChan);
      CPPUNIT_ASSERT(bdh.itsDelay.shape() == casacore::IPosition(3, 2, nBeam, nAnt));

      bdh.setIdealBandpass(1e6);
      CPPUNIT_ASSERT_DOUBLES_EQUAL(1e6, bdh.itsResolution, 1.);

      // fill in delays only with some validity flags
      for (casacore::uInt ant = 0; ant < nAnt; ++ant) {
           for (casacore::uInt beam = 0; beam < nBeam; ++beam) {
                const casacore::IPosition validPoint(3, beam == ant ? 1 : 0, beam, ant);
                const casacore::IPosition invalidPoint(3, beam == ant ? 0 : 1, beam, ant);
                // check that invalid point is invalid by default instead of assigning it
                CPPUNIT_ASSERT(!bdh.itsDelayValid(invalidPoint));
                bdh.itsDelayValid(validPoint) = true;
                // some unique number of nanoseconds for each antenna and beam combination
                bdh.itsDelay(validPoint) = (ant * nBeam + beam) * 1e-9;
           }
      } 

      // now apply the delays set above to the ideal bandpass
      bdh.applyDelays();

      // reset delays and validity flags
      bdh.itsDelay.set(0.f);
      bdh.itsDelayValid.set(false);

      // compute delays and we should get figures close to what we put in originally
      bdh.calcDelays();

      // check the results
      for (casacore::uInt ant = 0; ant < nAnt; ++ant) {
           for (casacore::uInt beam = 0; beam < nBeam; ++beam) {
                const casacore::IPosition validPoint(3, beam == ant ? 1 : 0, beam, ant);
                const casacore::IPosition invalidPoint(3, beam == ant ? 0 : 1, beam, ant);
                CPPUNIT_ASSERT(!bdh.itsDelayValid(invalidPoint));
                CPPUNIT_ASSERT(bdh.itsDelayValid(validPoint));
                const float expected = (ant * nBeam + beam) * 1e-9;
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, bdh.itsDelay(validPoint), 1e-6);
           }
      } 
      
   }
};
    
} // namespace utils

} // namespace askap

