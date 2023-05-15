/// @file
/// @brief Specialised gridder just for uv-weight construction
/// @details We don't need everything from the gridder (i.e. the actual gridding, CF generation, etc) for
/// the weight construction. This is essentially a cut down version of the Box gridder trimmed specifically
/// so it can only construct weight.
///
/// @copyright (c) 2023 CSIRO
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

// own includes
#include <askap/gridding/UVWeightGridder.h>

namespace askap {

namespace synthesis {

/// @brief default constructor
/// @note this class constructed via the default constructor will be useless without the builder set (via setUVWeightBuilder call)
UVWeightGridder::UVWeightGridder() : itsPaddingFactor(1.f), itsUCellSize(0.), itsVCellSize(0.), itsMaxPointingSeparation(-1.)
{}

/// @brief Initialise the gridding and the associated builder class
/// @details This method is supposed to be called before gridding first data. For convenience parameters resemble those
/// the proper gridders from the IVisGridder class hierarchy are using. In particular, the shape parameter is 4-dimensional
/// (as used for the gridders) with uSize, vSize, nPol and nChan as opposed to the 3-dimensional shape used for weight grids
/// (uSize, vSize, nChan - i.e. it is assumed that we always have the same weight for all polarisation products).
/// @param axes axes specifications
/// @param shape desired shape of the weight grid, same as passed to the proper gridder for image creation, i.e. u, v, pol, chan
/// @note this method plays the role of initialiseGrid in the gridder hierarchy
void UVWeightGridder::initialise(const scimath::Axes& axes, const casacore::IPosition& shape)
{
}

/// @brief process the visibility data.
/// @param acc const data accessor to work with
/// @note this method plays the role of 'generic' or 'grid' methods in the gridder hierarchy. I (MV) not sure at this stage whether
/// we need some selection methods to control what actually contributes to weights or should use the accessor selector instead 
/// (as this would be a separate iteration over the data anyway). The method is 'const' because the actual accumulation is done
/// by the builder and this class is unchanged except for various caches (like frequency mapper)
void UVWeightGridder::accumulate(accessors::IConstDataAccessor& acc) const
{
}

/// @brief obtain the tangent point
/// @details  This method extracts the tangent point (reference position) from the
/// coordinate system.
/// @return direction measure corresponding to the tangent point
casacore::MVDirection UVWeightGridder::getTangentPoint() const
{
   ASKAPCHECK(itsAxes.hasDirection(),"Direction axis is missing. axes="<<itsAxes);
   const casacore::Vector<casacore::Double> refVal(itsAxes.directionAxis().referenceValue());
   ASKAPDEBUGASSERT(refVal.nelements() == 2);
   const casacore::Quantum<double> refLon(refVal[0], "rad");
   const casacore::Quantum<double> refLat(refVal[1], "rad");
   const casacore::MVDirection out(refLon, refLat);
   return out;
}

// (MV): code duplication with the gridder classes tells me that we probably need to move it to Axes (and getTangentPoint too)
// there is a complication with shape and padding, though

/// @brief obtain the centre of the image
/// @details This method extracts RA and DEC axes from itsAxes and
/// forms a direction measure corresponding to the middle of each axis.
/// @return direction measure corresponding to the image centre
casacore::MVDirection UVWeightGridder::getImageCentre() const
{
   ASKAPCHECK(itsAxes.hasDirection(),"Direction axis is missing. axes="<<itsAxes);
   casacore::MDirection out;
   casacore::Vector<casacore::Double> centrePixel(2);
   ASKAPDEBUGASSERT(itsShape.nelements()>=2);
   ASKAPDEBUGASSERT(paddingFactor()>0);
   for (size_t dim=0; dim<2; ++dim) {
        centrePixel[dim] = double(itsShape[dim])/2./double(paddingFactor());
   }
   // MV: note, there were experiments running gridder code under OpenMP which would hit the issue with the lack of
   // thread safety for casacore routines. To abstract this out, syncHelper was written but it is a bit ugly to
   // include it outside of TableVisGridder. If thread safety is required here, one would need to change this code
   // to the way similar to TableVisGridder and probably made the access to syncHelper via a proper singleton pattern.
   // At this stage, I don't think we ever run this part from multiple threads, so I use casacore pixel to world conversion
   // directly.
   ASKAPCHECK(itsAxes.directionAxis().toWorld(out, centrePixel),
        "Unable to obtain world coordinates for the centre of the image. Something is wrong with the coordinate system");
   return out.getValue();
}

} // namespace synthesis

} // namespace askap

