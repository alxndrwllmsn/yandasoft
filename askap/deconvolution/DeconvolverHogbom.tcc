/// @file DeconvolverHogbom.tcc
/// @brief Class for a deconvolver based on the Hogbom CLEAN
/// @details This concrete class defines a deconvolver used to estimate an
/// image from a dirty image, psf optionally using a mask and a weights image.
/// @ingroup Deconvolver
///
///
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
/// @author Tim Cornwell <tim.cornwell@csiro.au>
///

#include <string>

#include <casacore/casa/aips.h>
#include <boost/shared_ptr.hpp>
#include <casacore/casa/Arrays/Array.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(dechogbomlogger, ".deconvolution.hogbom");

#include <askap/deconvolution/DeconvolverHogbom.h>

namespace askap {

    namespace synthesis {

        /// @brief Class for a deconvolver based on the Hogbom Clean
        /// @details This base class defines a deconvolver used to estimate an
        /// image from a dirty image, psf optionally using a mask and a weights image.
        /// The template argument T is the type, and FT is the transform
        /// e.g. DeconvolverHogbom<Double, DComplex>
        /// @ingroup Deconvolver

        template<class T, class FT>
        DeconvolverHogbom<T, FT>::~DeconvolverHogbom()
        {
        };

        template<class T, class FT>
        DeconvolverHogbom<T, FT>::DeconvolverHogbom(Vector<Array<T> >& dirty, Vector<Array<T> >& psf)
                : DeconvolverBase<T, FT>::DeconvolverBase(dirty, psf)
        {
            if (this->itsNumberDirtyTerms > 1) {
                throw(AskapError("Hogbom CLEAN cannot perform multi-term deconvolutions"));
            }
        };

        template<class T, class FT>
        DeconvolverHogbom<T, FT>::DeconvolverHogbom(Array<T>& dirty, Array<T>& psf)
                : DeconvolverBase<T, FT>::DeconvolverBase(dirty, psf)
        {
        };

        template<class T, class FT>
        void DeconvolverHogbom<T, FT>::initialise()
        {
            DeconvolverBase<T, FT>::initialise();
        }

        template<class T, class FT>
        bool DeconvolverHogbom<T, FT>::deconvolve()
        {
            this->initialise();

            ASKAPLOG_INFO_STR(dechogbomlogger, "Performing Hogbom CLEAN for " << this->control()->targetIter() << " iterations");
            do {
                this->oneIteration();
                this->monitor()->monitor(*(this->state()));
                this->state()->incIter();
            } while (!this->control()->terminate(*(this->state())));

            ASKAPLOG_INFO_STR(dechogbomlogger, "Performed Hogbom CLEAN for " << this->state()->currentIter() << " iterations");

            ASKAPLOG_INFO_STR(dechogbomlogger, this->control()->terminationString());

            this->finalise();

            return True;
        }

        template<class T, class FT>
        void DeconvolverHogbom<T, FT>::configure(const LOFAR::ParameterSet& parset)
        {
            DeconvolverBase<T, FT>::configure(parset);
        }

        // This contains the heart of the Hogbom Clean algorithm
        template<class T, class FT>
        bool DeconvolverHogbom<T, FT>::oneIteration()
        {
            const bool isMasked(this->weight(0).shape().conform(this->dirty(0).shape()));

            // Find peak in residual image
            casacore::IPosition minPos;
            casacore::IPosition maxPos;
            T minVal, maxVal;
            if (isMasked) {
                casacore::minMaxMasked(minVal, maxVal, minPos, maxPos, this->dirty(0), this->weight(0));
                minVal = this->dirty(0)(minPos);
                maxVal = this->dirty(0)(maxPos);
            } else {
                casacore::minMax(minVal, maxVal, minPos, maxPos, this->dirty(0));
            }
            //
            ASKAPLOG_INFO_STR(dechogbomlogger, "Maximum = " << maxVal << " at location " << maxPos);
            ASKAPLOG_INFO_STR(dechogbomlogger, "Minimum = " << minVal << " at location " << minPos);

            T absPeakVal = 0.0;
            casacore::IPosition absPeakPos;
            if (abs(minVal) < abs(maxVal)) {
                absPeakVal = maxVal;
                absPeakPos = maxPos;
            } else {
                absPeakVal = minVal;
                absPeakPos = minPos;
            }

            this->state()->setPeakResidual(absPeakVal);
            this->state()->setObjectiveFunction(absPeakVal);
            this->state()->setTotalFlux(sum(this->model()));

            // Has this terminated for any reason?
            if (this->control()->terminate(*(this->state()))) {
                return True;
            }

            const uInt nx(this->psf(0).shape()(0));
            const uInt ny(this->psf(0).shape()(1));

            const IPosition subPsfShape = this->findSubPsfShape();

            // Now we adjust model and residual for this component
            const casacore::IPosition residualShape(this->dirty(0).shape().nonDegenerate());
            const IPosition subPsfStart(2, nx / 2 - subPsfShape(0) / 2, ny / 2 - subPsfShape(1) / 2);
            const IPosition subPsfEnd(2, nx / 2 + subPsfShape(0) / 2 - 1, ny / 2 + subPsfShape(1) / 2 - 1);
            const IPosition subPsfStride(2, 1, 1);

            Slicer subPsfSlicer(subPsfStart, subPsfEnd, subPsfStride, Slicer::endIsLast);

            const casacore::IPosition psfShape(2, nx, ny);

            casacore::IPosition residualStart(2, 0), residualEnd(2, 0), residualStride(2, 1);
            casacore::IPosition psfStart(2, 0), psfEnd(2, 0), psfStride(2, 1);

            const casacore::IPosition modelShape(this->model(0).shape().nonDegenerate());
            casacore::IPosition modelStart(2, 0), modelEnd(2, 0), modelStride(2, 1);

            const casacore::IPosition peakPSFPos = this->getPeakPSFPosition();
            ASKAPDEBUGASSERT(peakPSFPos.nelements() >= 2);
            // Wrangle the start, end, and shape into consistent form.
            for (uInt dim = 0; dim < 2; dim++) {
                residualStart(dim) = max(0, Int(absPeakPos(dim) - psfShape(dim) / 2));
                residualEnd(dim) = min(Int(absPeakPos(dim) + psfShape(dim) / 2 - 1), Int(residualShape(dim) - 1));
                // Now we have to deal with the PSF. Here we want to use enough of the
                // PSF to clean the residual image.
                psfStart(dim) = max(0, Int(peakPSFPos(dim) - (absPeakPos(dim) - residualStart(dim))));
                psfEnd(dim) = min(Int(peakPSFPos(dim) - (absPeakPos(dim) - residualEnd(dim))),
                                  Int(psfShape(dim) - 1));

                modelStart(dim) = residualStart(dim);
                modelEnd(dim) = residualEnd(dim);
            }

            casacore::Slicer psfSlicer(psfStart, psfEnd, psfStride, Slicer::endIsLast);
            casacore::Slicer residualSlicer(residualStart, residualEnd, residualStride, Slicer::endIsLast);
            casacore::Slicer modelSlicer(modelStart, modelEnd, modelStride, Slicer::endIsLast);

            if (!(residualSlicer.length() == psfSlicer.length()) || !(residualSlicer.stride() == psfSlicer.stride())) {
                ASKAPLOG_INFO_STR(dechogbomlogger, "Peak of PSF  : " << peakPSFPos);
                ASKAPLOG_INFO_STR(dechogbomlogger, "Peak of residual: " << absPeakPos);
                ASKAPLOG_INFO_STR(dechogbomlogger, "Residual start  : " << residualStart << " end: " << residualEnd);
                ASKAPLOG_INFO_STR(dechogbomlogger, "PSF   start  : " << psfStart << " end: " << psfEnd);
                ASKAPLOG_INFO_STR(dechogbomlogger, "Residual slicer : " << residualSlicer);
                ASKAPLOG_INFO_STR(dechogbomlogger, "PSF slicer   : " << psfSlicer);
                throw AskapError("Mismatch in slicers for residual and psf images");
            }

            // Add to model
            this->model()(absPeakPos) = this->model()(absPeakPos) + this->control()->gain() * absPeakVal;

            // Subtract entire PSF from residual image

            this->dirty()(residualSlicer) = this->dirty()(residualSlicer)
                                            - this->control()->gain() * absPeakVal * this->psf()(psfSlicer);

            return True;
        }

    } // namespace synthesis

} // namespace askap
