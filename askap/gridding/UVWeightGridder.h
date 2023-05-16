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

   /// @brief default constructor
   /// @note this class constructed via the default constructor will be useless without the builder set (via setUVWeightBuilder call)
   UVWeightGridder();

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

   // as mentioned in the notes for the accumulate method, it may be more correct (from design-purist point of view)
   // to delegate all data selection to the accessor level. However, gridders implement some selection (and some is
   // even implicit like wmax rejection which would be very difficult to take into account in a generic way). So
   // we lack this functionality in the accessor selection. To move foward faster, I (MV) will copy some of this
   // gridder functionality here. It can be removed later on, if we ever had a cleaner redesign of gridder classes.

   /// @brief set the largest angular separation between the pointing centre and the image centre
   /// @details If the threshold is positive, it is interpreted as the largest allowed angular
   /// separation between the beam (feed in the accessor terminology) pointing centre and the
   /// image centre. This option matches the gridder option - We need it to allow imaging of a subset of data 
   /// (i.e. a smaller field of view) and reject all pointings located outside this smaller image. All accessor rows with
   /// pointingDir1 separated from the image centre by more than this threshold will be ignored in accumulate.
   /// If the threshold is negative (default), no data rejection based on the pointing direction is done.
   /// The class is initialised by default with a negative threshold, i.e. all data are used by default.
   /// @param[in] threshold largest allowed angular separation in radians, use negative value to select all data
   void inline maxPointingSeparation(double threshold = -1.) { itsMaxPointingSeparation = threshold; }

   /// @brief set or reset flag controlling selection of the representative beam and pointing
   /// @details Change itsDoBeamAndFieldSelection
   /// @param[in] doSelection new value of the flag
   void inline useAllDataForPSF(const bool doSelection) { itsDoBeamAndFieldSelection = doSelection;}

   /// @brief assign source index to be used for all future accumulated weights
   /// @details This is essentially an arbitrary index (zero by default) which is passed to the
   /// builder class as the 3rd parameter. There is no particular need to have it in our current use cases,
   /// it has been added because gridder classes support it for various research-related experiments.
   void setSourceIndex(casacore::uInt index) { itsSourceIndex = index; }

protected:

   /// @brief obtain current field index
   /// @details Although it is not great, we use the fact that only one field (i.e. dish pointing)
   /// can be represented by a single accessor. It is the case in the current implementation, but
   /// is not, strictly speaking, required by the interface or MS standard. In principle, only potentially
   /// bad performance stops us doing it per row rather than per accessor, so the limitation is not fundamental.
   /// This method returns the field 
   /// corresponding to the accessor passed during the last call to indexField.
   /// @return current field index
   inline casacore::uInt currentField() const { return itsCurrentField; }
  
   /// @brief checks whether the current field has been updated
   /// @details See currentField for the description of limitations. This method detects field changes in the field pointing (and numbers them in the 
   /// order they are encountered). If at a later stage we find that the fields need to be numbered in a particular way, this can be implemented.
   /// @note To match implementation of the gridder classes, we detect changes in the pointing of the first encountered beam. It has implications if
   /// either 3rd axis is operated in a non-tracking way or accessor row structure is different from one iteration to another. I (MV) suspect it was done
   /// this way because in early days we're trying to simulate equatorial vs. alt-az mounts and, technically, physical beam pointing matters.
   /// @param[in] acc input const accessor to analyse
   void indexField(const accessors::IConstDataAccessor &acc);

   /// @brief obtain the tangent point
   /// @details This method extracts the tangent point (reference position) from the
   /// coordinate system.
   /// @return direction measure corresponding to the tangent point
   casacore::MVDirection getTangentPoint() const;

   /// @brief obtain the centre of the image
   /// @details This method extracts RA and DEC axes from itsAxes and
   /// forms a direction measure corresponding to the middle of each axis.
   /// @return direction measure corresponding to the image centre
   casacore::MVDirection getImageCentre() const;

   
   // check whether we need to keep padding factor like the gridder does, may be it is sufficient to pass it to initialise method
   // and later use the shape or cell size  

   /// @brief obtain padding factor
   /// @details To mimic the behaviour of proper gridders, this class also implements optional padding where the actual grid used
   /// by this class is slightly larger than what the user has requested.
   /// @return current padding factor
   float inline paddingFactor() const { return itsPaddingFactor;}

private:

   /// @brief mapping class between image planes and accessor channels
   /// @details Correspondence between planes of the image cube and accessor channels may be
   /// non-trivial. This class takes care of the mapping.
   /// @note (MV) this approach was just copied from the gridder, ideally a proper reprojection of the spectral axis is required
   mutable FrequencyMapper itsFreqMapper;

   /// @brief Axes definition for the image associated to this weight grid
   /// @details It is needed for proper cell sizes, etc and set via the initialise method
   askap::scimath::Axes itsAxes;

   /// @brief shape of the associated image
   /// @note it includes extra dimensions the gridder would deal with
   casacore::IPosition itsShape;

   /// @brief internal padding factor, 1 by default
   float itsPaddingFactor;    

   /// @brief cell sizes in wavelengths along the U direction
   double itsUCellSize;

   /// @brief cell sizes in wavelengths along the U direction
   double itsVCellSize;

   /// @brief uv weight builder
   /// @details The builder object is responsible for the actual weight book-keeping. This cutdown
   /// gridder class doesn't make much sense without the builder set. Therefore, accumulate method
   /// throws an exception if this is the case at that stage.
   boost::shared_ptr<IUVWeightBuilder> itsUVWeightBuilder;

   /// @brief largest angular separation between the pointing centre and the image centre
   /// @details If the value is positive, it is interpreted as the largest allowed angular
   /// separation between the beam (feed in the accessor terminology) pointing centre and the
   /// image centre. It is intended to allow imaging of a subset of data (i.e. smaller field of view)
   /// and reject all pointings located outside this smaller image. All accessor rows with
   /// pointingDir1 separated from the image centre by more than this threshold are ignored.
   /// If the value is negative, no data rejection based on the pointing direction is done.
   /// Values are in radians.
   double itsMaxPointingSeparation;

   // the following fields are used to implement functionality similar to representative 
   // feed and field selection for PSF in the ordinary gridders (i.e. if multiple 
   // beams and pointings are gridded onto the same grid, we don't want them all to contribute
   // to the uv-weight).

   /// @brief true if no visibilities have been accumulated since the last initialise
   /// @details By default, we only take the first encountered beam (feed in the accessor 
   /// terminology) and field (it is an approximation that they all are the same). This flag
   /// is reset in the initialise call enabling reuse of the object with potentially different 
   /// representative feed and field.
   mutable bool itsFirstAccumulatedVis;

   /// @brief an index of the beam (feed in accessor) which is accepted
   /// @details This data member is initialized when the first visibility is accumulated,
   /// only this beam (feed in the accessor terminology) contributes to the weight
   /// @note the value only makes sense if itsFirstAccumulatedVis is false and itsDoBeamAndFieldSelection is true
   mutable casacore::uInt itsSelectedBeam;

   /// @brief pointing direction of the beam which is accepted
   /// @details This data member is initialized when the first visibility is accumulated
   /// enabling selection of the particular field contributing to the weight
   /// @note the value only makes sense if itsFirstAccumulatedVis is false and itsDoBeamAndFieldSelection is true
   mutable casacore::MVDirection itsSelectedPointing;

   /// @brief flag controlling data selection for accumulation
   /// @details By default we only accumulate a representative beam (feed in the accessor terminology)
   /// and field to construct the weight. For research purposes we need an option which allows us to take all 
   /// available data into account. Setting this is flag to false will result in itsSelectedBeam and itsSelectedPointing 
   /// being ignored. The default value is true.
   /// @note The flag fulfilling the similar role in the ordinary gridders is defined in the opposite 
   /// sense (false to do the selection).
   bool itsDoBeamAndFieldSelection;

   /// @brief current "source index" to be passed to the builder class
   /// @details This is essentially an arbitrary index (zero by default) which is passed to the
   /// builder class as the 3rd parameter. There is no particular need to have it in our current use cases,
   /// it has been added because gridder classes support it for various research-related experiments.
   casacore::uInt itsSourceIndex;

   /// @brief cache for the current field index
   /// @details See indexField for more info
   mutable casacore::uInt itsCurrentField;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef ASKAP_SYNTHESIS_GRIDDING_UV_WEIGHT_GRIDDER_H

