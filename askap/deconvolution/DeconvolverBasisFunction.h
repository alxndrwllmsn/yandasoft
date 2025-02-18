/// @file DeconvolverBasisFunction.h
/// @brief Class for a BasisFunction-Clean-based deconvolver
/// @details This interface class defines a deconvolver used to estimate an
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

#ifndef ASKAP_SYNTHESIS_DECONVOLVERBASISFUNCTION_H
#define ASKAP_SYNTHESIS_DECONVOLVERBASISFUNCTION_H

#include <string>

#include <boost/shared_ptr.hpp>
#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/Array.h>

#include <askap/deconvolution/DeconvolverBase.h>
#include <askap/deconvolution/DeconvolverState.h>
#include <askap/deconvolution/DeconvolverControl.h>
#include <askap/deconvolution/DeconvolverMonitor.h>
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
        /// e.g. DeconvolverBasisFunction<Double, DComplex>
        /// @ingroup Deconvolver
        template<class T, class FT>
        class DeconvolverBasisFunction : public DeconvolverBase<T, FT> {

            public:
                typedef boost::shared_ptr<DeconvolverBasisFunction<T, FT> > ShPtr;

                /// @brief Construct from dirty image and psf
                /// @detail Construct a deconvolver from a dirty image and
                /// the corresponding PSF. Note that both dirty image
                /// and psf can have more than 2 dimensions. We use a vector
                /// here to allow multiple dirty images and PSFs for the
                /// same model (e.g. as in MFS)
                /// @param[in] dirty Dirty image (array)
                /// @param[in] psf Point Spread Function (array)
                DeconvolverBasisFunction(casacore::Vector<casacore::Array<T> >& dirty,
                                         casacore::Vector<casacore::Array<T> >& psf);

                /// @brief Construct from dirty image and psf
                /// @detail Construct a deconvolver from a dirty image and
                /// the corresponding PSF. Note that both dirty image
                /// and psf can have more than 2 dimensions. We keep this
                /// version for compatibility
                /// @param[in] dirty Dirty image (array)
                /// @param[in] psf Point Spread Function (array)
                DeconvolverBasisFunction(casacore::Array<T>& dirty, casacore::Array<T>& psf);

                virtual ~DeconvolverBasisFunction();

                /// @brief Set the basis function to be used
                /// @details The algorithm can work with different basis functions
                /// PointBasisFunction, MultiScaleBasisFunction.
                /// @param[in] bf Shared pointer to basisfunction instance
                void setBasisFunction(const boost::shared_ptr<BasisFunction<T> >& bf);

                /// @brief Return the basis function to be used
                /// @details The algorithm can work with different basis functions
                /// PointBasisFunction, MultiScaleBasisFunction
                const boost::shared_ptr<BasisFunction<T> >& basisFunction() const;

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
            protected:

                /// @brief report in the log on the product of coupling matrix and its inverse
                void reportOnCouplingMatrix() const;

                /// @brief report in the log on the level of coupling between adjacent scales
                void reportOnAdjacentScaleCoupling() const;

            private:

                /// @brief Perform the deconvolution
                /// @detail This is the main deconvolution method.
                bool oneIteration();

                void initialisePSF();

                void initialiseResidual();

                void minMaxMaskedScales(T& minVal, T& maxVal,
                                        casacore::IPosition& minPos, casacore::IPosition& maxPos,
                                        const casacore::Array<T>& dataArray,
                                        const casacore::Array<T>& maskArray);

                // Find the coefficients for each scale by applying the
                // inverse of the coupling matrix
                casacore::Vector<T> findCoefficients(const casacore::Matrix<casacore::Double>& invCoupling,
                                                 const casacore::Vector<T>& peakValues);

                /// @brief basically matrix multiplication across the basis function domain
                /// @details This method applies inverse coupling matrix to every pixel of the data array
                /// @param[in] invCoupling inverse coupling matrix (i.e. the matrix multiplied by the data array)
                /// @param[in] dataArray data array (with basis function decomposition)
                /// @return new data array after the inverse coupling matrix is applied
                static casacore::Array<T> applyInverse(const casacore::Matrix<casacore::Double>& invCoupling,
                                            const casacore::Array<T>& dataArray);

                static casacore::Array<T> apply(const casacore::Matrix<casacore::Double>& invCoupling,
                                     const casacore::Vector<T> dataVector);

                /// Residual images convolved with basis functions
                casacore::Array<T> itsResidualBasisFunction;

                /// Point spread functions convolved with basis functions
                casacore::Array<T> itsPSFBasisFunction;

                /// Use cross terms in the source removal step?
                casacore::Bool itsUseCrossTerms;

                /// The coupling between different scales.
                casacore::Matrix<casacore::Double> itsCouplingMatrix;

                /// Inverse of the coupling matrix
                casacore::Matrix<casacore::Double> itsInverseCouplingMatrix;

                /// Determinant of the coupling Matrix
                casacore::Double itsDetCouplingMatrix;

                /// Point spread functions convolved with cross terms of basis functions
                casacore::Array<T> itsPSFCrossTerms;

                casacore::Bool itsDecouple;

                casacore::String itsDecouplingAlgorithm;

                /// Basis function used in the deconvolution
                boost::shared_ptr<BasisFunction<T> > itsBasisFunction;

                /// We keep track of the strength and location of each component
                /// identified. This allows calculation of the L1 norm of the
                /// model. We use clean to minimise the L1 so this is a good
                /// check to make. Ideally for stokes I, this should be equal to the flux.
                casacore::Vector<casacore::Array<T> > itsL1image;

                /// The flux subtracted on each scale
                casacore::Vector<T> itsScaleFlux;

                /// The peak of the convolved PSF as a function of scale
                casacore::Vector<T> itsPSFScales;
        };

    } // namespace synthesis

} // namespace askap

#include <askap/deconvolution/DeconvolverBasisFunction.tcc>

#endif  // #ifndef I_DECONVOLVERBASISFUNCTION_H
