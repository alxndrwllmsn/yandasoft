/// @file
/// @brief Simple implementation of IUVWeightAccessor interface
/// @details It holds the collection of weights with flat index and converts the indices
/// required by the interace by applying some custom index translator (a linear translation is
/// setup by default via GenericUVWeightIndexTranslator). 
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
#include <askap/gridding/GenericUVWeightAccessor.h>

namespace askap {

namespace synthesis {

/// @brief construct weight accessor from the given collection
/// @details The default translation of indices is set up by this method. Also, by default, all 
/// coefficients are zero which means that only index 0 will be used from the collection. 
/// A different index translation class can be assigned via the appropriate setter method of UVWeightIndexTranslationHelper.
/// The constructor is not "explicit" by intention. This way we can cast it from the collection type directly.
/// Otherwise, the resulting index is coeffBeam * beam + coeffField * field + coeffSource * source
/// @param[in] wts weight collection
/// @param[in] coeffBeam beam index coefficient
/// @param[in] coeffField field index coefficient
/// @param[in] coeffSource source index coefficient
GenericUVWeightAccessor::GenericUVWeightAccessor(const UVWeightCollection &wts, casacore::uInt coeffBeam,
                           casacore::uInt coeffField, casacore::uInt coeffSource) : itsWeights(wts)
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
UVWeight GenericUVWeightAccessor::getWeight(casacore::uInt beam, casacore::uInt field, casacore::uInt source) const
{
   const casacore::uInt index = indexOf(beam, field, source);
   return itsWeights.get(index);
}

} // namespace synthesis

} // namespace askap

