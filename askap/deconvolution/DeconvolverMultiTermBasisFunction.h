/// @file DeconvolverMultiTermBasisFunction.h
/// @brief Class for a BasisFunction-Clean-based deconvolver
/// @details This interface class defines a deconvolver used to estimate an
/// image from a dirty image, psf optionally using a mask and a weights image.
/// This version can deal with multiple terms
/// D = B(0)*I(0) + B(1)*I(1) + B(2)*I(2)
/// The most common example is where the B's are the spectral dirty beams
/// and the I(0) are the Taylor series approximation to the frequency
/// dependent frequencies.
/// @ingroup Deconvolver
///
///
/// @copyright (c) 2007,2024 CSIRO
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
/// @author Mark Wieringa <Mark.Wieringa@csiro.au>
///

#ifndef ASKAP_SYNTHESIS_DECONVOLVERMULTITERMBASISFUNCTION_H
#define ASKAP_SYNTHESIS_DECONVOLVERMULTITERMBASISFUNCTION_H

#include <string>

#include <boost/shared_ptr.hpp>

#include <askap/deconvolution/DeconvolverBase.h>
#include <askap/deconvolution/BasisFunction.h>


namespace askap {

    namespace synthesis {

        /// @brief Class for a deconvolver using the BasisFunction Clean algorithm
        /// @details This base class defines a deconvolver used to estimate an
        /// image from a dirty image, psf optionally using a mask and a weights image.
        /// This algorithm is similar to the MultiScale Clean (Cornwell 2009) with changes
        /// to improve performance and flexibility.
        ///
        /// The template argument T is the type, and FT is the transform
        /// e.g. DeconvolverBasisFunction<double, DComplex>
        /// @ingroup Deconvolver
        template<class T, class FT>
        class DeconvolverMultiTermBasisFunction : public DeconvolverBase<T, FT> {

            public:

                typedef boost::shared_ptr<DeconvolverMultiTermBasisFunction<T, FT>> ShPtr;

                /// @brief Construct from dirty image and psf
                /// @detail Construct a deconvolver from a dirty image and
                /// the corresponding PSF. Note that both dirty image
                /// and psf can have more than 2 dimensions. We use a vector
                /// here to allow multiple dirty images and PSFs for the
                /// same model (e.g. as in MFS)
                /// @param[in] dirty Dirty image (array)
                /// @param[in] psf Point Spread Function (array)
                /// @param[in] psf Point Spread Function containing 2*nTaylor-1 terms (array)
                DeconvolverMultiTermBasisFunction(Vector<Array<T>>& dirty,
                                                  Vector<Array<T>>& psf,
                                                  Vector<Array<T>>& psfLong);

                /// @brief Construct from dirty image and psf
                /// @detail Construct a deconvolver from a dirty image and
                /// the corresponding PSF. Note that both dirty image
                /// and psf can have more than 2 dimensions. We keep this
                /// version for compatibility
                /// @param[in] dirty Dirty image (array)
                /// @param[in] psf Point Spread Function (array)
                DeconvolverMultiTermBasisFunction(Array<T>& dirty, Array<T>& psf);

                virtual ~DeconvolverMultiTermBasisFunction();

                /// @brief Set the basis function to be used
                /// @details The algorithm can work with different basis functions
                /// PointBasisFunction, MultiScaleBasisFunction.
                /// @param[in] bf Shared pointer to basisfunction instance
                void setBasisFunction(boost::shared_ptr<BasisFunction<T>> bf);

                /// @brief Return the basis function to be used
                /// @details The algorithm can work with different basis functions
                /// PointBasisFunction, MultiScaleBasisFunction
                boost::shared_ptr<BasisFunction<T>> basisFunction();

                /// @brief Set the type of solution used in finding the optimum component
                void setSolutionType(const std::string& solutionType);

                // Perform many iterations using OpenMP
                void ManyIterations();

                /// @brief Perform the deconvolution
                /// @detail This is the main deconvolution method.
                virtual bool deconvolve();

                /// @brief Initialize the deconvolution
                /// @detail Initialise e.g. set weighted mask
                virtual void initialise();

                /// @brief Finalise the deconvolution
                /// @detail Finalise the deconvolution
                virtual void finalise();

                /// @brief configure basic parameters of the solver
                /// @details This method encapsulates extraction of basic solver parameters from the parset.
                /// @param[in] parset parset
                virtual void configure(const LOFAR::ParameterSet &parset);

                /// @brief Update only the dirty image
                /// @detail Update an existing deconvolver for a changed dirty image
                /// @param[in] dirty Dirty image (array)
                /// @param[in] term term to update
                virtual void updateDirty(Array<T>& dirty, uInt term = 0);

                /// @brief Update only the dirty images
                /// @detail Update an existing deconvolver for a changed dirty images.
                /// @param[in] dirty Dirty image (vector of arrays)
                virtual void updateDirty(Vector<Array<T>>& dirty);

                /// @brief export the scale mask
                /// @detail Give access to the scale mask used in the deconvolution
                /// @return Matrix<T> scale mask (bitmask of scales for each pixel)
                const Matrix<T> scaleMask();

                /// @brief import initial scale mask
                /// @detail Load an initial scale mask to use in the deconvolution
                /// @param[in]scaleMask a Matrix<T> with bitmask of scales for each pixel
                void setScaleMask(const Matrix<T>& scaleMask);

                /// @brief release working memory not needed between major cycles
                /// @details Deconvolvers can use a lot of memory, try to release as
                /// much as possible without losing essential state
                virtual void releaseMemory();

            private:

                // Initialise the PSFs - only need to do this once per change in basis functions
                void initialisePSF();

                // Initialise the residual bases if required
                void initialiseResidual();

                // Initialise the deep clean masks if required
                void initialiseMask();

                // Initialise the PSFs - only need to do this once per change in basis functions
                virtual void initialiseForBasisFunction(bool force);

                // Not currently implemented, extract from ManyIterations and reinstate after timer update
                void chooseComponent(uInt& optimumBase, IPosition& absPeakPos, T& absPeakVal, Vector<T>& peakValues);

                /// Fill the vector of high (or active) pixels for each base
                void fillHighPixelList(std::vector<std::vector<uInt>>& highPixels, const Matrix<T>& weight);

                /// Work out the slice of the residual and psf to use
                void getResidualAndPSFSlice(const IPosition& absPeakPos,
                    IPosition& shape, IPosition& residualStart, IPosition& psfStart);

                /// Add the optimum basisfunction component to the model
                void addComponentToModel(const Vector<T>& peakValues,
                        const IPosition& shape, const IPosition& resStart, const IPosition& psfStart,
                        const uInt optimumBase, Matrix<T>& model);

                /// Subtract the basisfunction component from all residual bases and terms
                void subtractPSF(const Vector<T>& peakValues,
                    const IPosition& shape, const IPosition& resStart, const IPosition& psfStart,
                    uInt optimumBase);

                /// Subtract the basisfunction component from all residual bases and terms
                /// This version uses the list of active pixels instead of the full arrays
                void subtractPSFPixels(const Vector<T>& peakValues, const IPosition& peakPos,
                    uInt optimumBase, bool useHighPixels, const std::vector<std::vector<uInt>>& highPixels);

                /// Long vector of PSFs
                Vector<Matrix<T>> itsPsfLongVec;

                /// Residual images convolved with basis functions, [nx,ny][nterms][nbases]
                Vector<Vector<Matrix<T>>> itsResidualBasis;

                /// Number of scales is limited to 24 as we write the scale mask out as a float image
                /// we could use std::numeric_limits<T>::digits, but for double
                /// the uInt used doesn't have enough range - anyway 24 scales seems plenty
                const uInt itsMaxScales = 24;

                /// Bitmask listing active pixels for each scale
                /// Would need to use ulong/size_t instead of uInt for images > 64k^2
                std::vector<std::set<uInt>> itsScalePixels;

                /// Point spread functions convolved with cross terms
                // [nxsub,nysub][nterms,nterms][nbases,nbases]
                Matrix<Matrix<Matrix<T>>> itsPSFCrossTerms;

                /// The coupling between different terms for each basis [nterms,nterms][nbases]
                Vector<Matrix<double>> itsCouplingMatrix;

                /// Inverse of the coupling matrix [nterms,nterms][nbases]
                Vector<Matrix<double>> itsInverseCouplingMatrix;

                /// Basis function used in the deconvolution
                boost::shared_ptr<BasisFunction<T>> itsBasisFunction;

                /// The flux subtracted on each term and scale [nterms][nbases]
                Vector< Vector<T>> itsTermBaseFlux;

                /// Flag to indicate the dirty / residual images have been updated
                bool itsDirtyChanged;

                /// FLag to indicate the BasisFunction has been updated
                bool itsBasisFunctionChanged;

                /// The Clean solution type - MAXBASE or MAXCHISQ
                std::string itsSolutionType;

                /// Use a list of high / active pixels
                bool itsUsePixelLists;

                /// @brief pixel list tolerance
                /// @details Multiply the clean threshold by (1-tolerance) to decide which
                /// pixels go in the pixel list for a major cycle. Suggested value 0.1
                float itsPixelListTolerance;

                /// @brief pixel list n sigma limit
                /// @detail Don't put pixels in the pixellist with amplitude < limit*noise
                /// The noise is determined separately for each scale / residual basis
                float itsPixelListNSigma;

                /// Read a pre-existing scale mask with the name given if not empty
                std::string itsScaleMaskName;
        };

    } // namespace synthesis

} // namespace askap

#include <askap/deconvolution/DeconvolverMultiTermBasisFunction.tcc>

#endif
