/// @file DeconvolveTimerUTils.cc
///
/// @copyright (c) 2016 CSIRO
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
/// @author Minh Vuong <minh.vuong@csiro.au>
///

#include <askap/utils/DeconvolveTimerUtils.h>

#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapUtil.h>

#include <fstream>
#include <iomanip>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

ASKAP_LOGGER(logger, ".deconvolveTimerUtils");

using namespace askap;
using namespace askap::utils;

Timer::Timer()
{
#ifdef HAVE_MPI
    int flag = 0;
    int r = MPI_Initialized(&flag);
    if (r == MPI_SUCCESS && flag == 1) {
        itsTimer.reset(new MPITimer {});
    } else {
        ASKAPLOG_INFO_STR(logger,"MPI is installed but MPI_Initialized() is not yet called. Using StdTimer instead");
        itsTimer.reset(new StdTimer {});
    }
#else
    itsTimer.reset(new StdTimer {});
#endif
}

double Timer::elapsedTime() const
{
    return itsTimer->elapsedTime();
}

void Timer::start()
{
    itsTimer->start();
}

void Timer::stop()
{
    itsTimer->stop();
}

std::string
Timer::summary() const
{
    return itsTimer->summary();
}

/////////////////////////////////

StdTimer::StdTimer()
: itsElapsedTime(0)
{
    itsState = State::STOP;
}

void StdTimer::start()
{
    if ( itsState == State::STOP ) {
        auto now = std::chrono::system_clock::now();
        itsStartTime = std::chrono::system_clock::to_time_t(now);
        itsState = State::START;
    }
}

void StdTimer::stop()
{
    if ( itsState == State::START ) {
        auto now = std::chrono::system_clock::now();
        itsStopTime = std::chrono::system_clock::to_time_t(now);
        if ( itsElapsedTime == 0 ) {
            itsElapsedTime = itsStopTime - itsStartTime;
        } else {
            itsElapsedTime += itsStopTime - itsStartTime;
        }
        itsState = State::STOP;
    }
}

std::string StdTimer::summary() const
{
    ASKAPASSERT(itsState == State::STOP);
    std::string summary = std::string("Elapsed Time: ") + std::to_string(itsElapsedTime);
    return summary;
}

double StdTimer::elapsedTime() const
{
    ASKAPASSERT(itsState == State::STOP);
    double et =  static_cast<double>(itsElapsedTime);
    return et;
}

//////////////////////////
#ifdef HAVE_MPI
MPITimer::MPITimer()
:  itsElapsedTime(0.0)
{
    itsState = State::STOP;
}


void MPITimer::start()
{
    if ( itsState == State::STOP ) {
        itsStartTime = MPI_Wtime();
        itsState = State::START;
    }
}

void MPITimer::stop()
{
    if ( itsState == State::START ) {
        itsStopTime = MPI_Wtime();;
        if ( itsElapsedTime == 0.0 ) {
            itsElapsedTime = itsStopTime - itsStartTime;
        } else {
            itsElapsedTime += itsStopTime - itsStartTime;
        }
        itsState = State::STOP;
    }
}

std::string MPITimer::summary() const
{
    ASKAPASSERT(itsState == State::STOP);
    std::string summary = std::string("Elapsed Time: ") + std::to_string(itsElapsedTime);
    return summary;
}

double MPITimer::elapsedTime() const
{
    ASKAPASSERT(itsState == State::STOP);
    return itsElapsedTime;
}
#endif

//////////////////////////////
SectionTimer::SectionTimer(unsigned int sections)
{
    for(unsigned int section = 0; section < sections; section++) {
        itsTimers.insert(std::make_pair(section, new Timer{}));
    }
}

void SectionTimer::start(unsigned int section)
{
    auto timerIter = itsTimers.find(section);
    if ( timerIter != itsTimers.end() ) {
        auto timer = timerIter->second;
        timer->start();
    } else {
        ASKAPLOG_WARN_STR(logger,"SectionTimer::start - section = " 
                            << section << " is not in the map");
    }
}

void SectionTimer::stop(unsigned int section)
{
    auto timerIter = itsTimers.find(section);
    if ( timerIter != itsTimers.end() ) {
        auto timer = timerIter->second;
        timer->stop();
    } else {
        ASKAPLOG_WARN_STR(logger,"SectionTimer::start - section = "
                            << section << " is not in the map");
    }
}

void SectionTimer::summary() const
{
    for ( const auto& kvp : itsTimers ) {
        const auto& timer = kvp.second;
        ASKAPLOG_INFO_STR(logger,"timer " << kvp.first << ", " << timer->summary());
    }
}

double SectionTimer::totalElapsedTime() const
{
    double total = 0.0;
    for ( const auto& kvp : itsTimers ) {
        const auto& timer = kvp.second;
        total += static_cast<double> (timer->elapsedTime());
    }
    return total;
}
