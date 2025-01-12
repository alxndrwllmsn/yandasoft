/// @file
///
/// Helper class to manage bandpass importing/exporting as well as 
/// delay estimates and some operations to aggregate several solutions.
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
///

// casa

// own
#include <askap/utils/BandpassDelayHelper.h>
#include <askap/scimath/utils/DelayEstimator.h>
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>

ASKAP_LOGGER(logger, ".BandpassDelayHelper");

namespace askap {

namespace utils {

/// @brief constructor
/// @param[in] nAnt number of antennas
/// @param[in] nBeam number of beams
/// @param[in] nChan number of spectral channels
/// @note We need these parameters to setup the buffer appropriately. In the current calibration accessor interfaces
/// there is no way to get this information (because it was intended that the setup won't change). For table-based accessor
/// which we use, the intention is to set them up from the same parset values as for the output table. In principle, this tool can
/// be extended in the future to expand the data (e.g. casting to more channels) or to select a subset of the data. But currently it is out of
/// scope.
BandpassDelayHelper::BandpassDelayHelper(casacore::uInt nAnt, casacore::uInt nBeam, casacore::uInt nChan) :
     itsBandpass(casacore::IPosition(4, nChan, 2, nBeam, nAnt)), itsBandpassValid(casacore::IPosition(4, nChan, 2, nBeam, nAnt), false),
     itsDelay(2, nBeam, nAnt), itsDelayValid(2, nBeam, nAnt, false) {}
    
/// @brief calculate delays
/// @details This method estimates delays for all antennas, beams and polarisations for the currently
/// buffered bandpass solution. Note, it must be read or obtained some other way. For now, the simplest apporach
/// of phase-slope fitting is used. It is possible to extend it into two stage approach (FFT-spectral averaging-phase slope)
/// similar to the one used by the delay solver later on. 
void BandpassDelayHelper::calcDelays() 
{
   scimath::DelayEstimator de(itsResolution);
}
    
/// @brief add delays from another object
/// @details This method adds delays stored in another object of this type (which maps to another calibration table). It is a bit messy to do it
/// this way but handy to encapsulate all operations with validity flags (this logic can be extracted into a separate class which we can return instead of
/// having such data as data members of this class).
/// @param[in] other other instance of the class with calculated delays (which are added to the delays stored in this class)
void BandpassDelayHelper::addDelays(const BandpassDelayHelper &other)
{
   const casacore::IPosition shape = itsDelay.shape();
   const casacore::IPosition otherShape = other.itsDelay.shape();
   ASKAPDEBUGASSERT(shape.nelements() == 3u);
   ASKAPDEBUGASSERT(otherShape.nelements() == 3u);
   ASKAPDEBUGASSERT(shape[0] == otherShape[0]);
   ASKAPDEBUGASSERT(shape[0] == 2);
   if (otherShape > shape) {
       ASKAPLOG_INFO_STR(logger, "adding delays computed for a bandpass table with larger number of antennas or beams - some data are going to be ignored");
   }
   if (otherShape < shape) {
       ASKAPLOG_INFO_STR(logger, "adding delays computed for a bandpass table with less antennas or beams - some good data may be flagged as invalid");
   }
   for (int ant = 0; ant < shape[2]; ++ant) {
        for (int beam = 0; beam < shape[1]; ++beam) {
             const bool otherIsOutOfBounds = (ant >= otherShape[2]) || (beam >= otherShape[1]);
             for (int pol = 0; pol < shape[0]; ++pol) {
                  itsDelayValid(pol, beam, ant) &= !otherIsOutOfBounds && other.itsDelayValid(pol, beam, ant);
                  if (itsDelayValid(pol, beam, ant)) {
                      itsDelay(pol, beam, ant) += other.itsDelay(pol, beam, ant);
                  }
             }
        }
   }
}

/// @brief flip the sign of stored delays
/// @details This method essentially multiplies all stored delays by -1. It can be used together with addDelays to implement the subtraction or,
/// together with calcDelays and applyDelays, to compensate the delays present in the bandpass solution
void BandpassDelayHelper::negateDelays()
{
   // just flip the sign for everything, validity flags remain the same
   itsDelay *= -1.f;
}

/// @brief apply currently stored delays to the currently stored bandpass
/// @details This method modifies currently stored bandpass by applying the phase slope implied by the appropriate delay term. So, combined with
/// calcDelays and negateDelays calls, this can be used to remove the best fit delay from the loaded bandpass.
void BandpassDelayHelper::applyDelays()
{
   const casacore::IPosition shape = itsDelay.shape();
   ASKAPDEBUGASSERT(itsBandpass.shape().nelements() == shape.nelements() + 1);
   // keep the type signed because IPosition class has signed coordinates
   const int nChan = itsBandpass.shape()[0];
   ASKAPDEBUGASSERT(shape == itsBandpass.shape().getLast(shape.nelements()));

   for (int ant = 0; ant < shape[2]; ++ant) {
        for (int beam = 0; beam < shape[1]; ++beam) {
             for (int pol = 0; pol < shape[0]; ++pol) {
                  const bool delayValid = itsDelay(pol, beam, ant);
                  // taking the slice to get bandpass for the given pol, beam, ant and benefit from reference semantics
                  // (note, due to the order of axes the resulting vector should be contiguous in memory and hence, it should be efficient)
                  const casacore::IPosition start(4, 0, pol, beam, ant);
                  const casacore::IPosition length(4, nChan, 1, 1, 1);
                  casacore::Vector<casacore::Complex> bpVec = itsBandpass(casacore::Slicer(start, length));
                  casacore::Vector<bool> bpValidVec = itsBandpassValid(casacore::Slicer(start, length));
                  if (delayValid) {
                      // apply the slope, keep the original bandpass validity flags
                      // note - check the sign
                      const float phaseStepPerChannel = 2.*casacore::C::pi * itsResolution * itsDelay(pol, beam, ant);
                      for (int chan = 0; chan < nChan; ++chan) {
                           bpVec[chan] *= casacore::polar(1.f, phaseStepPerChannel * chan);
                      }
                  } else {
                      // invalidate all bandpass points as the corresponding delay term is invalid
                      bpValidVec.set(false);
                  }
             }
        }
   }  
}
    
/// @brief load bandpass via calibration accessor
/// @details This method loads bandpass parameters from the calibration solution accessor for the number of
/// antennas, beams and channels (two polarisations are always expected) this class is initialised with. We
/// also need spectral resolution (which is not stored with the calibration table as the design assumed a fixed configuration)
/// to be able to compute (or apply) delays later on in physical units. Working with physical units of delay (rather than per channel) allows
/// us to combine results obtained with different resolution (e.g. one bandpass table is in full resolution and another is averaged to 1 MHz)
/// @param[in] acc calibration solution accessor to read bandpass from
/// @param[in] resolution spectral resolution in Hz
void BandpassDelayHelper::loadBandpass(const accessors::ICalSolutionConstAccessor &acc, double resolution)
{
   itsResolution = resolution;
}

/// @brief store bandpass into a calibration accessor
/// @details This method does the opposite of loadBandpass and stores the current bandpass parameters using the supplied calibration 
/// solution accessor. Only valid parameters are stored (which ensures that they remain flagged according to the general calibration accessor logic),
/// but an exception is expected if, say, the number of antennas or beams in the current bandpass with valid data is larger than the given accessor 
/// can handle. Note, this cannot occur if the output accessor and this class are initialised with the same antennas, beams, etc numbers or if 
/// the bandpass of the same array has been loaded via loadBandpass first before any manipulation.
/// @param[in] acc calibration accessor to store the bandpass into
void BandpassDelayHelper::storeBandpass(accessors::ICalSolutionAccessor &acc) const
{
}

/// @brief initialise the class with ideal bandpass
/// @details This method initialises the class with the amplitude of 1. and zero phase. It doesn't affect the delays if they're calculated earlier.
/// Thus, this method can be used to export a fake bandpass containing just the delay slope or the inverted one providing the correction term.
/// This option exists largely for experiments. All antennas, beams, channels and polarisations this class has been setup to handle are marked as valid.
/// However, the validity flags of delays will be respected if delays are applied to this ideal bandpass.
/// @param[in] resolution spectral resolution in Hz to be used in subsequent calcDelays and applyDelays calls (as for the bandpass loaded via the accessor)
void BandpassDelayHelper::setIdealBandpass(double resolution)
{
   itsResolution = resolution;
   itsBandpass.set(casacore::Complex(1.f, 0.f));
   itsBandpassValid.set(true);
}

} // namespace utils

} // namespace askap
