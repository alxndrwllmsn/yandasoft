/// @file
///
/// @brief Calibration effect: frequency-dependent antenna leakages
/// @details This effect represents bandpass leakage calibration parameters
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
/// @author Max Voronkov <maxim.voronkov@csiro.au>
/// @author Mark Wieringa <mark.wieringa@csiro.au>

#ifndef FREQ_DEPENDENT_LEAKAGE_H
#define FREQ_DEPENDENT_LEAKAGE_H

// own includes
#include <askap/scimath/fitting/ComplexDiffMatrix.h>
#include <askap/scimath/fitting/ComplexDiff.h>
#include <askap/scimath/fitting/Params.h>
#include <askap/dataaccess/IConstDataAccessor.h>
#include <askap/askap/AskapError.h>
#include <askap/measurementequation/ParameterizedMEComponent.h>
#include <askap/scimath/utils/PolConverter.h>

// std includes
#include <string>
#include <utility>

namespace askap {

namespace synthesis {

/// @brief Calibration effect: frequency-dependent antenna leakages
/// @details This effect represents bandpass leakage calibration parameters. Note, an
/// exception is thrown if the number of spectral channels in the accessor doesn't match
/// the dimension of the bandpass parameter
/// @ingroup measurementequation
struct FreqDependentLeakage : public ParameterizedMEComponent<true> {

   /// @brief constructor, store reference to paramters
   /// @param[in] par shared pointer to parameters
   inline explicit FreqDependentLeakage(const scimath::Params::ShPtr &par) :
                              ParameterizedMEComponent<true>(par) {}

   /// @brief main method returning Mueller matrix and derivatives
   /// @details This method has to be overloaded (in the template sense) for
   /// all classes representing various calibration effects. CalibrationME
   /// template will call it when necessary. It returns
   /// @param[in] chunk accessor to work with
   /// @param[in] row row of the chunk to work with
   /// @return ComplexDiffMatrix filled with Mueller matrix corresponding to
   /// this effect
   inline scimath::ComplexDiffMatrix get(const accessors::IConstDataAccessor &chunk,
                                casacore::uInt row) const;

};

} // namespace synthesis

} // namespace askap

#include <askap/measurementequation/FreqDependentLeakage.tcc>

#endif // #ifndef FREQ_DEPENDENT_LEAKAGE_H
