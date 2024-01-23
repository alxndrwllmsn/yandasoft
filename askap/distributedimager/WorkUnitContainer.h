/// @file
/// @brief Class managing a collection of work units
/// @details This class is intended to encapsulate management of a vector of work units
/// (adding elements, various preprocessing like squashing channels, iteration and breaking
/// iteration in a certain way). It allows us to separate iteration logic from the rest of
/// the algorithm and thus to have a cleaner code. 
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

#ifndef ASKAP_SYNTHESIS_DISTRIBUTEDIMAGER_WORK_UNIT_CONTAINER_H
#define ASKAP_SYNTHESIS_DISTRIBUTEDIMAGER_WORK_UNIT_CONTAINER_H

// own includes
#include <askap/messages/ContinuumWorkUnit.h>

// std includes
#include <vector>

// boost includes
#include <boost/noncopyable.hpp>

namespace askap {

namespace synthesis {

/// @brief Class managing a collection of work units
/// @details This class is intended to encapsulate management of a vector of work units
/// (adding elements, various preprocessing like squashing channels, iteration and breaking
/// iteration in a certain way). It allows us to separate iteration logic from the rest of
/// the algorithm and thus to have a cleaner code. 
/// @ingroup distributedimager
class WorkUnitContainer : public boost::noncopyable { 
public:
   /// @brief add a work unit to the container
   /// @details The new unit is prepended to the existing vector of work units
   /// (to match the code behaviour prior to refactoring)
   /// @param[in] wu work unit to add
   void add(const cp::ContinuumWorkUnit &wu);

   /// @brief squash work units with adjacent channels into one work unit
   /// @details This method was copied pretty much as it was from ContinuumWorker as part of
   /// the refactoring. It modifies the container in situ by merging work units corresponding to 
   /// adjacent channels. This allows us to save on processing in the continuum case.
   void compressWorkUnits(); 

   /// @brief obtain the number of work units in the container
   /// @return the current number of stored workunits
   inline size_t size() const { return itsWorkUnits.size(); }

private:
   /// @brief vector of work units 
   /// @note It may be more appropriate to use some other container type here based on the usage, but at this stage 
   /// mimic the original code prior to refactoring to minimise work
   std::vector<cp::ContinuumWorkUnit> itsWorkUnits;
};

} // namespace synthesis

} // namespace askap

#endif // #ifndef ASKAP_SYNTHESIS_DISTRIBUTEDIMAGER_WORK_UNIT_CONTAINER_H

