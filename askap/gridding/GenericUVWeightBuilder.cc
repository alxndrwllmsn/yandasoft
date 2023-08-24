/// @file
/// @brief Simple implementation of IUVWeightBuilder interface
/// @details It holds the collection of weights with flat index and converts the indices
/// required by the interace by applying some custom index translator (a linear translation is
/// setup by default via GenericUVWeightIndexTranslator) in a way similar to GenericUVWeightAccessor.
/// This is sufficient for ignoring some or all indices (via multiplying them by zero) or translate
/// into a flat index with the desired order. 
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
#include <askap/gridding/GenericUVWeightBuilder.h>

namespace askap {

namespace synthesis {

/// @brief construct weight builder from the given collection
/// @details The default translation of indices is set up by this method. Also, by default, all 
/// coefficients are zero which means that only index 0 will be used from the collection. 
/// A different index translation class can be assigned via the appropriate setter method of UVWeightIndexTranslationHelper.
/// Unlike with weight accessor class, the builder owes the collection. It is setup on demand using parameters passed in 
/// the initialised method. The translated index is coeffBeam * beam + coeffField * field + coeffSource * source
/// @param[in] coeffBeam beam index coefficient
/// @param[in] coeffField field index coefficient
/// @param[in] coeffSource source index coefficient
GenericUVWeightBuilder::GenericUVWeightBuilder(casacore::uInt coeffBeam,
                           casacore::uInt coeffField, casacore::uInt coeffSource) : itsUSize(0u), itsVSize(0u), itsNPlanes(0u)
{
   // could've used the constructor directly, but it is safer this way and more clear that we won't have a dangling pointer
   const boost::shared_ptr<GenericUVWeightIndexTranslator> translator(new GenericUVWeightIndexTranslator(coeffBeam, coeffField, coeffSource));
   setTranslator(translator);
}

/// @brief obtain weight grid for a given index
/// @details index interpretation is left for the implementation which can be different
/// for different gridders. Trying to make it more general, we pass beam, field and source indices
/// as they are currently used in various gridders. This can change in the future, if necessary.
/// @note Strictly speaking, detection of field changes will be subject to some assumptions. And the 
/// numbering will be per gridder (figured out from the accessor supplied), other gridders may have
/// different indexing. The field index is expected to match currentField() output in AProjectGridderBase.
/// @param[in] beam beam index (from accessor for the given row, it is assumed that we don't have cross-beam correlations)
/// @param[in] field field index if the gridder is a mosaicing one, zero otherwise
/// @param[in] source source index used to form direction-dependent offset index (not sure if it is needed, 
/// add it here just to keep things general as it is used in the gridder code)
/// @return UVWeight object with the selected grid of weights
UVWeight GenericUVWeightBuilder::getWeight(casacore::uInt beam, casacore::uInt field, casacore::uInt source) const
{
   const casacore::uInt index = indexOf(beam, field, source);
   return itsWeights.get(index);
}


/// ensure uv weights are collected on the grid with the same dimensions. At this stage we assume
/// that the gridder has the same grid size for all potential weight grids (i.e. for different beams).
/// This is the only case we support anyway inside the gridder.
/// @param[in] uSize size in the direction of u-coordinate
/// @param[in] vSize size in the direction of v-coordinate
/// @param[in] nPlanes number of planes for the given weight (3rd dimension)
/// @note This class is agnostic about the physical pixel sizes (e.g. uv cell size) -
/// we assume it is always the same as in the setup of the gridder
void GenericUVWeightBuilder::initialise(casacore::uInt uSize, casacore::uInt vSize, casacore::uInt nPlanes) 
{
   itsUSize = uSize;
   itsVSize = vSize;
   itsNPlanes = nPlanes;
}

/// @brief obtain weight grid for writing for the given metadata (indices)
/// @details index interpretation is left for the implementation which can be different
/// for different gridders. Trying to make it more general, we pass beam, field and source indices
/// as they are currently used in various gridders. This can change in the future, if necessary.
/// @note Strictly speaking, detection of field changes will be subject to some assumptions. And the 
/// numbering will be per gridder (figured out from the accessor supplied), other gridders may have
/// different indexing. The field index is expected to match currentField() output in AProjectGridderBase.
/// In this class, the index translator (initially setup to the generic version with parameters passed in the constructor)
/// is used to convert generic indices passed to this method to the flat index of the collection.
/// @param[in] beam beam index (from accessor for the given row, it is assumed that we don't have cross-beam correlations)
/// @param[in] field field index if the gridder is a mosaicing one, zero otherwise
/// @param[in] source source index used to form direction-dependent offset index (not sure if it is needed, 
/// add it here just to keep things general as it is used in the gridder code)
/// @return UVWeight object with the selected grid of weights
/// @note This method can be called multiple times for the same indices, if necessary. The new grid will be created on demand.
UVWeight GenericUVWeightBuilder::addWeight(casacore::uInt beam, casacore::uInt field, casacore::uInt source) 
{
   const casacore::uInt index = indexOf(beam, field, source);
   if (itsWeights.exists(index)) {
       return itsWeights.get(index);
   } 
   // have to initialise the new weight grid using default parameters
   // this can probably be checked only in debug mode as any potential issue triggering this would be a logic bug rather than user parameter settings
   ASKAPCHECK(itsUSize > 0u && itsVSize > 0u && itsNPlanes > 0u, "New weight grid is requested (beam = "<<beam<<", field = "<<field<<", source = "<<source<<" (index = "<<index<<
              ") but default shape is not set, perhaps initialise method has not been called");
   itsWeights.add(index, itsUSize, itsVSize, itsNPlanes);
   // the reference to the cube returned by get will be wrapped into a UVWeight object implicitly
   return itsWeights.get(index);
}

/// @brief merge weight information from the other builder
/// @details This method is expected to be used in conjunction with the EstimatorAdapter and normal equation tree reduction
/// routines. It enables merging multiple weight grids (corresponding to the same indices) together in distributed data
/// processing.
/// @param[in] src other builder to merge from
/// @note The source should be the object which can be cast to the same type
void GenericUVWeightBuilder::merge(const IUVWeightBuilder &src) {
   try {
      const UVWeightCollection &otherCollection = dynamic_cast<const GenericUVWeightBuilder&>(src).itsWeights;
      itsWeights.merge(otherCollection);
   }
   catch (const std::bad_cast&) {
      ASKAPTHROW(AskapError, "Logic error! An attempt to merge an incompatible type in GenericUVWeightBuilder::merge");
   }
}

/// @brief finalise uv weight computation by applying some post-processing to accumulated weight grids
/// @details This method is expected to call the given object function for each weight in the collection.
/// Examples are calculation of robust weights, ensure Hermitian symmetry, uv-taper or a combination of those.
/// It is envisaged that we'll work with a chain of algorithms, but keep this complexity outside of this interface.
/// @param[in] calc processing algorithm/object function 
/// @return a reference to the uv-weight collection which is the final result (may be useful if the weight collection is
/// owed by implementation of this interface instead of being supplied by reference).
UVWeightCollection& GenericUVWeightBuilder::finalise(const IUVWeightCalculator &calc) 
{
   // in principle, we can add a safeguard here to ensure finalise method is called only once
   std::set<casacore::uInt> indices = itsWeights.indices();
   // we can add OpenMP parallelism later for either index or plane loop (or both), if required
   // (although care must be taken with slice creation as this part is not thread safe)
   for (std::set<casacore::uInt>::const_iterator ci = indices.begin(); ci != indices.end(); ++ci) {
        ASKAPDEBUGASSERT(itsWeights.exists(*ci));
        casacore::Cube<float> &cube = itsWeights.get(*ci);
        for (casacore::uInt plane = 0; plane < cube.nplane(); ++plane) {
             casacore::Matrix<float> slice = cube.xyPlane(plane);
             calc.process(slice);
        }
   }

   return itsWeights;
}

} // namespace synthesis

} // namespace askap

