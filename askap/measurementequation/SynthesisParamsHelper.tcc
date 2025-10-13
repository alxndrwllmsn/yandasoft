/// @copyright (c) 2007 CSIRO
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

#ifndef SYNSYNTHESISPARAMSHELPER_TCC
#define SYNSYNTHESISPARAMSHELPER_TCC

#include <askap/imagemath/utils/MultiDimArrayPlaneIter.h>
#include <askap/scimath/utils/PaddingUtils.h>
#include <askap/gridding/SupportSearcher.h>


#include <casacore/lattices/LatticeMath/Fit2D.h>


//#include <askap/askap/AskapLogging.h>
//ASKAP_LOGGER(logger, ".measurementequation.synthesisparamshelper");

//#include <askap/askap/AskapError.h>
//#include <casacore/measures/Measures/Stokes.h>

using namespace askap::scimath;
using namespace casa;

namespace askap
{
  namespace synthesis
  {
    template<typename T>
    inline casacore::Vector<casacore::Quantum<double> > SynthesisParamsHelper::fitBeam(const casacore::Array<T> &psfArray,
       const scimath::Axes &axes, const double cutoff, const int maxsupport)
    {
        ASKAPCHECK(axes.hasDirection(), "Direction axes are missing from the PSF parameter, unable to convert pixels to angular units");
        const casacore::Vector<casacore::Double> increments = axes.directionAxis().increment();
        ASKAPCHECK(increments.nelements() == 2, "Expect just two elements for increments of the direction axis, you have "<<
                   increments);
        ASKAPCHECK(increments[1]>0, "Expect positive increment on the declination axis. increments="<<increments);
        ASKAPCHECK(fabs(fabs(increments[0])-fabs(increments[1]))<1e-6,
                   "Different cell sizes mean that the current beam fitting code would give a wrong position angle. increments="
                   <<increments);

       casacore::Vector<double> result = fitBeam(psfArray, cutoff, maxsupport);

       casacore::Vector<casacore::Quantum<double> > beam(3);
       beam[0] = casacore::Quantum<double>(fabs(increments[0])*result[0],"rad");
       beam[1] = casacore::Quantum<double>(fabs(increments[1])*result[1],"rad");
       // position angle in radians
       // fitBeam call above already does the pi/2 offset
       //double pa = increments[0]<0 ? result[2] - casacore::C::pi/2 : casacore::C::pi/2 - result[2];
       double pa = increments[0]<0 ? result[2] : -result[2];
       if (pa < -casacore::C::pi/2) {
           pa += casacore::C::pi;
       }
       beam[2] = casacore::Quantum<double>(pa,"rad");
       return beam;
    }

    template<typename T>
    inline casacore::Vector<double> SynthesisParamsHelper::fitBeam(const casacore::Array<T> &psfArray,
        const double cutoff, const int maxsupport) {

       const casacore::IPosition shape = psfArray.shape();
       ASKAPCHECK(shape.nelements()>=2,"PSF image is supposed to be at least 2-dimensional, shape="<<psfArray.shape());
       ASKAPCHECK(cutoff>0 && cutoff<1,"beam cutoff level should be between 0 and 1");
       // make maxsupport a valid odd value
       int maxSupport = maxsupport;
       maxSupport = casacore::min(maxSupport,shape(0)-1);
       maxSupport = casacore::min(maxSupport,shape(1)-1);
       maxSupport = casacore::max(3,maxSupport);
       maxSupport = 2 * (maxSupport/2) + 1;

       if (shape.product() != shape[0]*shape[1]) {
           //ASKAPLOG_WARN_STR(logger, "Multi-dimensional PSF is present (shape="<<shape<<
           //                  "), using the first 2D plane only to fit the beam");
       }
       ASKAPCHECK(shape[0] >= 3 && shape[1] >= 3, "Expect at least 3x3 pixel images, you have shape = "<<shape);

       // we need a non-const reference to use the plane iterator but make sure we don't change it
       casacore::Array<T> localPSF(psfArray);
       casacore::Matrix<T> psfSlice = imagemath::MultiDimArrayPlaneIter::getFirstPlane(localPSF).nonDegenerate();

       // search for support to speed up beam fitting
       //ASKAPLOG_INFO_STR(logger, "Searching for support with the relative cutoff of "<<cutoff<<" to speed fitting up");
       SupportSearcher ss(cutoff);
       ss.search(psfSlice);
       // search only looks in x and y direction, now extend the support to diagonals
       ss.extendedSupport(psfSlice);
       casacore::uInt support = ss.symmetricalSupport(psfSlice.shape());

       // limit support size to roughly 100 pixels per beam
       // some mostly flagged channels can have very large support, fit will take too long
       if (support > maxSupport) {
           // maybe we should just fail here and abort processing of the channel
           //ASKAPLOG_DEBUG_STR(logger, "Reducing support size for beam fit from "<<support<<" to "<<maxSupport);
           support = maxSupport;
       }

       if (support < 3) {
           support = 3;
       } else {
          if (support % 2 == 0) {
              // if even, move up to next odd.
              ++support;
          }
       }
       ASKAPDEBUGASSERT(support % 2 == 1);

       //ASKAPLOG_INFO_STR(logger, "Extracting support of "<<support<<" pixels for 2D gaussian fitting");
       const casa::IPosition newShape(2,support,support);
       for (int dim=0; dim<2; ++dim) {
            ASKAPCHECK(psfSlice.shape()[dim] >= int(support), "Support is greater than the original size, shape="<<
                       psfSlice.shape());
       }

       // make a copy since we rescale in next step and we're not allowed to change the psfArray argument
       casa::Array<T> tempPSFSlice;
       tempPSFSlice = scimath::PaddingUtils::centeredSubArray(psfSlice,newShape); 

       // normalise to 1 - technically unnecessary as we should have it already normalised
       const T maxPSF = casa::max(tempPSFSlice);
       if (fabs(maxPSF-1.)>1e-6) {
            tempPSFSlice /= maxPSF;
       }

       // actual fitting
       // the beam fitter fails at times - producing no or a bad solution
       // We've tried a few approaches:
       //  - retry with different support - often works, but still some failures
       //  - use setIncludeRange to exclude pixels below the cutoff - works sometimes, but worse at times
       //  - use estimate to set the initial guess - this seems the best solution so far
       casa::LogIO os;
       casa::Fit2D fitter(os);
       // do not set cutoff for small support, i.e. all pixels will be used for the fit
       // note, we can't get support less than 3 given the code above
       if (support > 3) {
           // for large support case ensure that we get at least 4 neighbouring pixels near the peak above the cutoff
           double newCutoff = cutoff;
           const int centrePixel = static_cast<int>(support) / 2;
           for (int pix = 0; pix < 4; ++pix) {
                // this gives +/- 1 pixel from centre on each axis and all combinations, assumes integer math.
                const IPosition cursor(2, centrePixel + (pix % 2) * 2 - 1, centrePixel + (pix / 2) * 2 - 1);
                const float curVal = tempPSFSlice(cursor);
                if (curVal < newCutoff) {
                    newCutoff = curVal;
                }
           }
           if (newCutoff < cutoff) {
               //ASKAPLOG_DEBUG_STR(logger, "Reducing cutoff for fitter to ensure enough pixels are taken into account, new cutoff "<<newCutoff<<
               //          ", it probably means the beam is bad");
           }
           fitter.setIncludeRange(newCutoff,1.0);
       }
       // Using estimate to set the initial guess seems to make the fit more robust
       casa::Vector<casa::Double> initialEstimate = fitter.estimate(casa::Fit2D::GAUSSIAN,tempPSFSlice);
       //ASKAPLOG_DEBUG_STR(logger,"Initial beam fit: "<<initialEstimate);
       initialEstimate[0]=1.; // PSF peak is always 1
       initialEstimate[1]=newShape[0]/2; // centre
       initialEstimate[2]=newShape[1]/2; // centre
       casa::Vector<casa::Bool> parameterMask(6,casa::False);
       parameterMask[3] = casa::True; // fit maj
       parameterMask[4] = casa::True; // fit min
       parameterMask[5] = casa::True; // fit pa

       fitter.addModel(casa::Fit2D::GAUSSIAN,initialEstimate,parameterMask);
       const casa::Array<T> sigma(tempPSFSlice.shape(),T(1.0));
       const casa::Fit2D::ErrorTypes fitError = fitter.fit(tempPSFSlice,sigma);
       ASKAPCHECK(fitError == casa::Fit2D::OK, "Error fitting the beam. fitError="<<fitError<<
                  " message: "<<fitter.errorMessage());
       ASKAPCHECK(fitter.numberPoints() > 4, "The number of points available for fitting the beam ("<<fitter.numberPoints()<<
                  ") is too small. Either cell size or cutoff are too large");

       const casa::Vector<casa::Double> result = fitter.availableSolution();
       const casa::Vector<casa::Double> errors = fitter.availableErrors();

       //ASKAPLOG_DEBUG_STR(logger, "Got fit result (in pixels) "<<result<<" and uncertainties "<< errors<<
       //                   " number of points used for the fit: "<<fitter.numberPoints());
       ASKAPCHECK(result.nelements() == 6, "Expect 6 parameters for 2D gaussian, result vector has "<<result.nelements());
       // Check if fit is reasonable
       for (casacore::uInt i = 0; i < result.nelements(); ++i) {
            ASKAPCHECK(!isnan(result[i]), "Beam fitting produces NaN in the result ("<<result<<"), most likely beam is not sampled properly, decrease the cell size");
            ASKAPCHECK(!isinf(result[i]), "Beam fitting produces infinity in the result ("<<result<<"), most likely beam is not sampled properly, decrease the cell size");
       }
       ASKAPCHECK(result[3] <= 2.0 * support && result[4] <= 2.0 * support,
                  "Error fitting the beam. Gaussian width >2x support");
       casacore::Vector<double> beam(3);
       beam[0] = result[3];
       beam[1] = result[4];
       // position angle on the pixel grid in radians
       double pa = result[5] - casacore::C::pi/2;
       if (pa < -casacore::C::pi/2) {
           pa += casacore::C::pi;
       }
       beam[2] = pa;
       return beam;
    }
  }
}
#endif // SYNSYNTHESISPARAMSHELPER_TCC
