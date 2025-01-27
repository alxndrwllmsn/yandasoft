/// @file SphFuncVisGridder.tcc
///
/// @copyright (c) 2016 CSIRO
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
/// @author Daniel Mitchell <daniel.mitchell@csiro.au>
///
#ifndef SPHVISGRIDDER_TCC_
#define SPHVISGRIDDER_TCC_

#include <askap/scimath/utils/PaddingUtils.h>

namespace askap
{
  namespace synthesis
  {

    /// @brief estimate the spheroidal function at nu=1
    /// @param[in] func vector to be interpolated
    template<typename T>
    void SphFuncVisGridder::interpolateEdgeValues(casacore::Vector<T> &func)
    {
      const int length = func.shape()[0];
      ASKAPASSERT(length>3);

      //func(0) = func(1) + (func(1)-func(2)) + (func(1)-2.0*func(2)+func(3));
      func(0) = T(3.) * (func(1) - func(2)) + func(3);
      if (length%2==1) { // nu=1 for the last element as well
        func(length-1) = func(0);
      }

    }


    template<typename T>
    inline void SphFuncVisGridder::correctConvolution(casacore::Array<T>& grid,
                    scimath::SpheroidalFunction& sf, int support,
                    bool interpolate)
    {
        casacore::IPosition shape = grid.shape();
        ASKAPDEBUGASSERT(shape.nelements()>=2);
        ASKAPDEBUGASSERT(shape(0)>1);
        ASKAPDEBUGASSERT(shape(1)>1);

        const casacore::Int xHalfSize = shape(0)/2;
        const casacore::Int yHalfSize = shape(1)/2;
        casacore::Vector<double> ccfx(shape(0));
        casacore::Vector<double> ccfy(shape(1));

        // initialise buffers to enable a filtering of the correction
        // function in Fourier space.
        casacore::Vector<typename askap::scimath::ComplexTypeTrait<T>::type> bufx(shape(0));
        casacore::Vector<typename askap::scimath::ComplexTypeTrait<T>::type> bufy(shape(1));

        // note grdsf(1)=0.
        for (int ix=0; ix<shape(0); ++ix)
        {
            const double nux=std::abs(double(ix-xHalfSize))/double(xHalfSize);
            const double val = sf(nux);
            bufx(ix) = typename askap::scimath::ComplexTypeTrait<T>::type(val,0.0);
        }

        for (int iy=0; iy<shape(1); ++iy)
        {
            const double nuy=std::abs(double(iy-yHalfSize))/double(yHalfSize);
            const double val = sf(nuy);
            bufy(iy) = typename askap::scimath::ComplexTypeTrait<T>::type(val,0.0);
        }

        if (interpolate) {
            // The spheroidal is undefined and set to zero at nu=1, but that
            // is not the numerical limit. Estimate it from its neighbours.
            interpolateEdgeValues(bufx);
            interpolateEdgeValues(bufy);
        }

        // Fourier filter the spheroidal (crop in Fourier space in line with
        // gridding kernel support size)
        const bool doFiltering = true;
        if (doFiltering) {
            // Some more advanced gridders have support>3 (e.g. w-proj).
            //
            int support = 3;
            const typename askap::scimath::ComplexTypeTrait<T>::type maxBefore = bufx(shape(0)/2);
            scimath::fft(bufx, true);
            scimath::fft(bufy, true);
            for (int ix=0; ix<shape(0)/2-support; ++ix) {
                bufx(ix) = 0.0;
            }
            for (int ix=shape(0)/2+support+1; ix<shape(0); ++ix) {
                bufx(ix) = 0.0;
            }
            for (int iy=0; iy<shape(1)/2-support; ++iy) {
                bufy(iy) = 0.0;
            }
            for (int iy=shape(1)/2+support+1; iy<shape(1); ++iy) {
                bufy(iy) = 0.0;
            }
            scimath::fft(bufx, false);
            scimath::fft(bufy, false);
            // Normalise after filtering.
            const typename askap::scimath::ComplexTypeTrait<T>::type normalisation = maxBefore / bufx(shape(0)/2);
            bufx *= normalisation;
            bufy *= normalisation;
        }

        for (int ix=0; ix<shape(0); ++ix) {
            double val = real(bufx(ix));
            ccfx(ix) = casacore::abs(val) > 1e-10 ? 1.0/val : 0.;
        }
        for (int iy=0; iy<shape(1); ++iy) {
            double val = real(bufy(iy));
            ccfy(iy) = casacore::abs(val) > 1e-10 ? 1.0/val : 0.;
        }

        casacore::ArrayIterator<T> it(grid, 2);
        while (!it.pastEnd())
        {
            casacore::Matrix<T> mat(it.array());
            ASKAPDEBUGASSERT(int(mat.nrow()) <= shape(0));
            ASKAPDEBUGASSERT(int(mat.ncolumn()) <= shape(1));
            for (int iy=0; iy<shape(1); iy++)
            {
                for (int ix=0; ix<shape(0); ix++)
                {
                    mat(ix, iy)*=ccfx(ix)*ccfy(iy);
                }
            }
            it.next();
        }
    }

  }
}

#endif
