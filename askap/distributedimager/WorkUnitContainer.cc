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
void WorkUnitContainer::compressWorkUnits() {
    // This takes the list of workunits and reprocesses them so that all the contiguous
    // channels are compressed into single workUnits for multiple channels
    // this is not applicable for the spectral line experiment but can markedly reduce
    // the number of FFT required for the continuum processesing mode

    // In preProcessWorkUnit we made a list of all the channels in the allocation
    // but the workunit may contain different measurement sets so I suppose it is
    // globalChannel that is more important for the sake of allocation ... but
    // the selector only works on one measurement set.

    // So the upshot is this simple scheme cannot combine channels from different
    // measurement sets into the same grid as we are using the MS accessor as the vehicle to
    // provide the integration.

    // So we need to loop through our workunit list and make a new list that just contains a
    // single workunit for each contiguous group of channels.

    // First lets loop through our workunits

    std::vector<cp::ContinuumWorkUnit> compressedList; // probably easier to generate a new list

    ASKAPDEBUGASSERT(itsWorkUnits.size() > 0u);
    cp::ContinuumWorkUnit startUnit = itsWorkUnits[0];

    unsigned int contiguousCount = 1;
    int sign = 1;
    if (itsWorkUnits.size() == 1) {
        ASKAPLOG_WARN_STR(logger,"Asked to compress channels but workunit count 1");
    }
    cp::ContinuumWorkUnit compressedWorkUnit = startUnit;


    for ( size_t count = 1; count < itsWorkUnits.size(); ++count) {

        cp::ContinuumWorkUnit nextUnit = itsWorkUnits[count];

        const std::string startDataset = startUnit.get_dataset();
        const int startChannel = startUnit.get_localChannel();
        const std::string nextDataset = nextUnit.get_dataset();
        const int nextChannel = nextUnit.get_localChannel();
        if (contiguousCount == 1 && nextChannel == startChannel - 1) {
            // channels are running backwards
            sign = -1;
        }


        if (startDataset == nextDataset) { // same dataset
            ASKAPLOG_DEBUG_STR(logger,"nextChannel "<<nextChannel<<" startChannel "<<startChannel<<" contiguousCount"<< contiguousCount);
            if (nextChannel == (startChannel + sign * contiguousCount)) { // next channel is contiguous to previous
                ++contiguousCount;
                ASKAPLOG_DEBUG_STR(logger, "contiguous channel detected: count " << contiguousCount);
                if (sign < 0) {
                    compressedWorkUnit = nextUnit;
                }
                compressedWorkUnit.set_nchan(contiguousCount); // update the nchan count for this workunit
                /*
                // MV: do we really need to update the parset here? It doesn't make sense to me because it will be overwritten by successive
                // blocks. There is some technical debt here. When this class is used we need to ensure that parset is updated in the loop before
                // the actual processing takes place for the given work unit with squashed channels
                // ah, the comment below refers to this and AXA-1004 ticket
 

                // Now need to update the parset details
                string ChannelParam = "["+toString(contiguousCount)+","+
                    toString(compressedWorkUnit.get_localChannel())+"]";
                ASKAPLOG_DEBUG_STR(logger, "compressWorkUnit: ChannelParam = "<<ChannelParam);
                itsParset.replace("Channels",ChannelParam);
                */
            } else { // no longer contiguous channels reset the count
                contiguousCount = 0;
            }
        }
        else { // different dataset reset the count
            ASKAPLOG_DEBUG_STR(logger, "Datasets differ resetting count");
            contiguousCount = 0;
        }
        if (count == (itsWorkUnits.size()-1) || contiguousCount == 0) { // last unit
            ASKAPLOG_DEBUG_STR(logger, "Adding unit to compressed list");
            compressedList.insert(compressedList.end(),compressedWorkUnit);
            startUnit = nextUnit;
            compressedWorkUnit = startUnit;
        }

    }
    if (compressedList.size() > 0) {
        ASKAPLOG_INFO_STR(logger, "Replacing workUnit list of size " << itsWorkUnits.size() << " with compressed list of size " << compressedList.size());
        ASKAPLOG_INFO_STR(logger,"A corresponding change has been made to the parset");
        itsWorkUnits = compressedList;
    }
    else {
        ASKAPLOG_WARN_STR(logger,"No compression performed");
    }
    // MV: this check should probably go on the user's side as this class can implement the general case with no issues
    ASKAPCHECK(compressedList.size() < 2, "The number of compressed workunits is greater than one. Channel parameters may be incorrect - see AXA-1004 and associated technical debt tickets");
}

} // namespace synthesis

} // namespace askap
