/// @file TilingUtils.h
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
/// @author Mark Wieringa <mark.wieringa@csiro.au>
///

#ifndef ASKAP_YANDASOFT_TILING_UTILS_H
#define ASKAP_YANDASOFT_TILING_UTILS_H

#include <askap/dataaccess/IConstDataSource.h>

#include <string>


namespace askap {
namespace utils {

/// @brief set selection to distribute work over ranks by tiles
/// @param[in] IDataSelectorPtr sel : shared pointer to data selector
/// @param[in] string column : name of column to inspect (usually DATA or FLAG)
/// @param[in] uint nRanks: number of MPI ranks
/// @param[in] uint rank : MPI rank of current process
void distributeByTile(accessors::IDataSelectorPtr& sel, const std::string& column, uint nRanks, uint rank);

}
}
#endif
