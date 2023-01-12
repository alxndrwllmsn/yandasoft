/// @file
/// @brief Grid of weights in the uv-domain
/// @details This class is related to traditional weighting. It represents a
/// single grid with weight in the uv-domain. Conceptually, this is just a cube
/// (the 3rd dimension is necessary because our gridders can be setup with non-degenerate
/// third dimension). The role of this class is to restrict the access interface and enable
/// a possible extension in the future with an abstract interface and multiple possible
/// implementations.
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

#ifndef ASKAP_SYNTHESIS_GRIDDING_UV_WEIGHT_H
#define ASKAP_SYNTHESIS_GRIDDING_UV_WEIGHT_H

// casa includes
#include <casacore/casa/Arrays/Cube.h>
#include <casacore/casa/aipstype.h>

// boost includes
#include <boost/noncopyable.hpp>

namespace askap {

namespace synthesis {

/// @brief Grid of weights in the uv-domain
/// @details This class is related to traditional weighting. It represents a
/// single grid with weight in the uv-domain. Conceptually, this is just a cube
/// (the 3rd dimension is necessary because our gridders can be setup with non-degenerate
/// third dimension). The role of this class is to restrict the access interface and enable
/// a possible extension in the future with an abstract interface and multiple possible
/// implementations.
/// @ingroup gridding
struct UVWeight { //: public boost::noncopyable {

   /// @brief construct a weight class with the given shape
   /// @details This constructor initialises the weight class and fills the grid with zeros.
   /// @param[in] uSize size in the direction of u-coordinate
   /// @param[in] vSize size in the direction of v-coordinate
   /// @param[in] nPlanes number of planes for the given weight (3rd dimension)
   /// @note This class is agnostic about the physical pixel sizes (e.g. uv cell size) -
   /// we assume it is always the same as the setup of the gridder
   UVWeight(casacore::uInt uSize, casacore::uInt vSize, casacore::uInt nPlanes) :
        itsWeightCube(uSize, vSize, nPlanes, 0.f) {}

   // do we need a default constructor? Not needed if push reference semantics everywhere and don't 
   // have issues with containers like with early C++ standards

   // also a method to setup from Cube may be handy. In the reverse direction, we may need to rely on 
   // friend relation between say weight builder and this class. An alternative is to have a proper 
   // interface with abstract virtual functions but this seems to be an overkill at this stage.

   /// @brief wrap an exsiting cube to make a UVWeight object
   /// @details This constructor makes UVWeight object from an existing cube benefiting from
   /// reference semantics of the casacore arrays. We want it to be a possible type cast, so no "explicit".
   /// @param[in] wt an object to setup from (weight in the form of casacore cube; note - reference semantics)
   UVWeight(const casacore::Cube<float> &wt) : itsWeightCube(wt) {}

   /// @brief obtain weight for the given coordinates 
   /// @brief[in] u first coordinate
   /// @brief[in] v second coordinate
   /// @brief[in] plane 3rd coordinate
   /// @return weight for the given point
   /// @note there is no cross-check that the coordinate is within the bounds
   /// (as it is expected anyway on the user's side and we'll use this code in
   /// the inner-most loop)
   inline float operator()(casa::uInt u, casa::uInt v, casa::uInt plane) const
   { return itsWeightCube(u,v,plane); }

   /// @brief read-write access to weight at the given coordinates
   /// @brief[in] u first coordinate
   /// @brief[in] v second coordinate
   /// @brief[in] plane 3rd coordinate
   /// @return reference to weight for the given point
   /// @note there is no cross-check that the coordinate is within the bounds
   /// (as it is expected anyway on the user's side and we'll use this code in
   /// the inner-most loop)
   inline float& operator()(casa::uInt u, casa::uInt v, casa::uInt plane)
   { return itsWeightCube(u,v,plane); }

   /// @brief obtain a size of the grid in the u-direction
   /// @return a size of the grid in the u-direction
   inline casacore::uInt uSize() const { return itsWeightCube.nrow(); }

   /// @brief obtain a size of the grid in the v-direction
   /// @return a size of the grid in the v-direction
   inline casacore::uInt vSize() const { return itsWeightCube.ncolumn(); }

   /// @brief obtain a size of the grid on the 3rd axis
   /// @return a size of the grid along the 3rd dimension
   inline casacore::uInt nPlane() const { return itsWeightCube.nplane(); }


private:
   /// @brief array of values describing the weight
   casacore::Cube<float> itsWeightCube;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef ASKAP_SYNTHESIS_GRIDDING_UV_WEIGHT_H

