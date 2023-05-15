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

#ifndef ASKAP_SYNTHESIS_GRIDDING_UV_WEIGHT_GRIDDER_H
#define ASKAP_SYNTHESIS_GRIDDING_UV_WEIGHT_GRIDDER_H

// casa includes
#include <casacore/casa/aipstype.h>
#include <casacore/casa/BasicSL/Complex.h>

// own includes
#include <askap/gridding/IUVWeightBuilder.h>
#include <askap/dataaccess/IConstDataAccessor.h>
#include <askap/gridding/FrequencyMapper.h>
#include <askap/scimath/fitting/Axes.h>

// boost includes (although it would be included through interfaces)
#include <boost/shared_ptr.hpp>

namespace askap {

namespace synthesis {

/// @brief Specialised gridder just for uv-weight construction
/// @details We don't need everything from the gridder (i.e. the actual gridding, CF generation, etc) for
/// the weight construction. This is essentially a cut down version of the Box gridder trimmed specifically
/// so it can only construct weight.
/// @ingroup gridding
struct UVWeightGridder  {

   // need to think whether we should provide a constructor which sets the builder to the generic version up front. At least it would be a handy short-cut for now.

   /// @brief assign uv weight builder
   /// @details The UV weight gridder is to be used together with one of the builder classes doing actual
   /// accumulation. This separation of the roles allows us to be able to use the same builder/index translation
   /// with either generic gridders or with this class to avoid gridder overhead where the gridding of visibilities
   /// is not necessary at the time of the weight construction. This method is used to set the builder object to
   /// work with.
   /// @param[in] wtBuilder shared pointer to the weight builder
   /// @note Unlike the proper gridder, this class is not useful if the builder class is not set. So perhaps, this needs
   /// to be promoted to the constructor parameter. 
   inline void setUVWeightBuilder(const boost::shared_ptr<IUVWeightBuilder> &wtBuilder) { itsUVWeightBuilder = wtBuilder;}

   /// @brief Initialise the gridding and the associated builder class
   /// @details This method is supposed to be called before gridding first data. For convenience parameters resemble those
   /// the proper gridders from the IVisGridder class hierarchy are using. In particular, the shape parameter is 4-dimensional
   /// (as used for the gridders) with uSize, vSize, nPol and nChan as opposed to the 3-dimensional shape used for weight grids
   /// (uSize, vSize, nChan - i.e. it is assumed that we always have the same weight for all polarisation products).
   /// @param axes axes specifications
   /// @param shape desired shape of the weight grid, same as passed to the proper gridder for image creation, i.e. u, v, pol, chan
   /// @note this method plays the role of initialiseGrid in the gridder hierarchy
   void initialise(const scimath::Axes& axes, const casacore::IPosition& shape);

   /// @brief process the visibility data.
   /// @param acc const data accessor to work with
   /// @note this method plays the role of 'generic' or 'grid' methods in the gridder hierarchy. I (MV) not sure at this stage whether
   /// we need some selection methods to control what actually contributes to weights or should use the accessor selector instead 
   /// (as this would be a separate iteration over the data anyway). The method is 'const' because the actual accumulation is done
   /// by the builder and this class is unchanged except for various caches (like frequency mapper)
   void accumulate(accessors::IConstDataAccessor& acc) const;

   /*
   // check whether we need to keep padding factor like the gridder does, may be it is sufficient to pass it to initialise method
   // and later use the shape or cell size  

   /// @brief obtain padding factor
   /// @details To mimic the behaviour of proper gridders, this class also implements optional padding where the actual grid used
   /// by this class is slightly larger than what the user has requested.
   /// @return current padding factor
   float inline paddingFactor() const { return itsPaddingFactor;}
   */

   /// @brief obtain the tangent point
   /// @details  This method extracts the tangent point (reference position) from the
   /// coordinate system.
   /// @return direction measure corresponding to the tangent point
   casacore::MVDirection getTangentPoint() const;


private:

   /// @brief mapping class between image planes and accessor channels
   /// @details Correspondence between planes of the image cube and accessor channels may be
   /// non-trivial. This class takes care of the mapping.
   /// @note (MV) this approach was just copied from the gridder, ideally a proper reprojection of the spectral axis is required
   mutable FrequencyMapper itsFreqMapper;

   /// @brief Axes definition for the image associated to this weight grid
   /// @details It is needed for proper cell sizes, etc and set via the initialise method
   askap::scimath::Axes itsAxes;

   /// @brief cell sizes in wavelengths along the U direction
   double itsUCellSize;

   /// @brief cell sizes in wavelengths along the U direction
   double itsVCellSize;

   /// @brief uv weight builder
   /// @details The builder object is responsible for the actual weight book-keeping. This cutdown
   /// gridder class doesn't make much sense without the builder set. Therefore, accumulate method
   /// throws an exception if this is the case at that stage.
   boost::shared_ptr<IUVWeightBuilder> itsUVWeightBuilder;

};

} // namespace synthesis

} // namespace askap

#endif // #ifndef ASKAP_SYNTHESIS_GRIDDING_UV_WEIGHT_GRIDDER_H

