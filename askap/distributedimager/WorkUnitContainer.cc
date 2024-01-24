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

// own includes
#include "askap/distributedimager/WorkUnitContainer.h"
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>
#include <askap/askap/AskapUtil.h>

ASKAP_LOGGER(logger, ".WorkUnitContainer");

namespace askap {

namespace synthesis {

/// @brief add a work unit to the container
/// @details The new unit is prepended to the existing vector of work units
/// (to match the code behaviour prior to refactoring)
/// @param[in] wu work unit to add
void WorkUnitContainer::add(const cp::ContinuumWorkUnit &wu) {
   itsWorkUnits.insert(itsWorkUnits.begin(),wu); 
}

/// @brief squash work units with adjacent channels into one work unit
/// @details This method was copied pretty much as it was from ContinuumWorker as part of
/// the refactoring. It modifies the container in situ by merging work units corresponding to 
/// adjacent channels. This allows us to save on processing in the continuum case.
void WorkUnitContainer::mergeAdjacentChannels() 
{
   // current working element, everything before that one has already been processed
   std::vector<cp::ContinuumWorkUnit>::iterator cursorIt = itsWorkUnits.begin();
   for (; cursorIt != itsWorkUnits.end(); ++cursorIt) {
        size_t contiguousCount = 1u;
        int sign = 1;
        std::vector<cp::ContinuumWorkUnit>::const_iterator testIt = cursorIt;
        const unsigned int cursorChan = cursorIt->get_localChannel();
        for (++testIt; testIt != itsWorkUnits.end(); ++testIt, ++contiguousCount) {
             // break if beam or dataset changes (can add more conditions like that later on)
             if ((cursorIt->get_dataset() != testIt->get_dataset()) || (cursorIt->get_beam() != testIt->get_beam())) {
                 break;
             }
             const unsigned int testChan = testIt->get_localChannel();
             ASKAPDEBUGASSERT(cursorChan != testChan);
             if (contiguousCount == 1u) {
                 sign = testChan > cursorChan ? 1 : -1;
             }
             // gap condition
             if (testChan != cursorChan + sign * contiguousCount) {
                 break;
             }
        }
        if (contiguousCount > 1u) {
            // can merge all work units up to (but not including) testIt into the one pointed by cursorIt
            // For now leave frequencies, width, etc untouched and just update the number of channels. This matches the 
            // old behaviour prior to refactoring
            cursorIt->set_nchan(contiguousCount);
            if (sign < 0) {
                // always keep the lowest channel number for the group. By default the unit pointed to by cursorIt will have the first
                ASKAPDEBUGASSERT(cursorChan + 1 >= contiguousCount);
                cursorIt->set_localChannel(cursorChan + 1 - contiguousCount);
            }
            contiguousCount = 1u;
            std::vector<cp::ContinuumWorkUnit>::const_iterator obsoletePartIt = cursorIt;
            // we keep the element pointed by cursorIt, but delete the ones after that
            ++obsoletePartIt;
            ASKAPDEBUGASSERT(obsoletePartIt != itsWorkUnits.end());
            // the following invalidates iterators from obsoletePartIt onwards, including itsWorkUnits.end(), but we don't 
            // cache it (so the call as part of the loop condition is ok). The cursorIt iterator remains valid.
            itsWorkUnits.erase(obsoletePartIt, testIt);
        }
   }
  
   // the original code prior to refactoring updated parset as well from inside this procedure (to select a group of channels rather than one)
   // This needs to be performed elsewhere. I leave this original code here commented out as a reminder. 
   // string ChannelParam = "["+toString(contiguousCount)+","+
   //       toString(compressedWorkUnit.get_localChannel())+"]";
   //       ASKAPLOG_DEBUG_STR(logger, "compressWorkUnit: ChannelParam = "<<ChannelParam);
   //       itsParset.replace("Channels",ChannelParam);
}

} // namespace synthesis

} // namespace askap
