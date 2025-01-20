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

      // now check the content of the helper class to verify that bandpass has been loaded correctly along with the validity flags
      for (casacore::uInt ant = 0; ant < nAnt; ++ant) {
           for (casacore::uInt beam = 0; beam < nBeam; ++beam) {
                for (casacore::uInt chan = 0; chan < nChan; ++chan) {
                     const casacore::IPosition validPoint(4, chan, beam == ant ? 1 : 0, beam, ant);
                     const casacore::IPosition invalidPoint(4, chan, beam == ant ? 0 : 1, beam, ant);
                     //std::cout<<ant<<" "<<beam<<" "<<chan<<" : "<<bdh.itsBandpassValid(casacore::IPosition(4, chan, 0, beam, ant)) << " "<< bdh.itsBandpassValid(casacore::IPosition(4, chan, 1, beam, ant))<<" values: "<<bdh.itsBandpass(casacore::IPosition(4, chan, 0, beam, ant)) << " "<< bdh.itsBandpass(casacore::IPosition(4, chan, 1, beam, ant)) << " "<<std::endl;
                     CPPUNIT_ASSERT(bdh.itsBandpassValid(validPoint));
                     CPPUNIT_ASSERT(!bdh.itsBandpassValid(invalidPoint));
                     const casacore::Complex expectedVal = casacore::Complex(ant + 1, beam + 1) * (static_cast<int>(chan) + 1);
                     CPPUNIT_ASSERT_DOUBLES_EQUAL(casacore::real(expectedVal), casacore::real(bdh.itsBandpass(validPoint)), 1e-6);
                     CPPUNIT_ASSERT_DOUBLES_EQUAL(casacore::imag(expectedVal), casacore::imag(bdh.itsBandpass(validPoint)), 1e-6);
                }
           }
      }
   }
};
    
} // namespace utils

} // namespace askap

