/// @file DeconvolveTimerUtils.h
///
/// @copyright (c) 2022 CSIRO
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
/// @author Minh Vuong <minh.vvuong@csiro.au>
///

#ifndef ASKAP_DECONVOLVE_TIMER_UTILS_H
#define ASKAP_DECONVOLVE_TIMER_UTILS_H

#include <map>
#include <string>
#include <memory>
#include <ctime>
#include <new>


namespace askap {
namespace utils {

enum class State {START, STOP};

//class StdTimer;
//#ifdef HAVE_MPI
//class MPITimer;
//#endif

class ITimer {
    public:
        virtual ~ITimer() {}
        virtual void start() = 0;
        virtual void stop() = 0;
        virtual std::string summary() const = 0;
        virtual double elapsedTime() const = 0;

    protected:
        askap::utils::State itsState;
};

class Timer final : public ITimer {
    public:
        Timer();
        Timer(const Timer&) = default;
        Timer& operator=(const Timer&) = default;
        ~Timer() {}
        void start() override;
        void stop() override;
        std::string summary() const override;
        double elapsedTime() const override;

    private:
        std::shared_ptr<ITimer> itsTimer;
};

class StdTimer final : public ITimer
{
  public:
    friend class Timer;
    ~StdTimer() {}
  private:
    StdTimer();
    StdTimer(const StdTimer&) = default;
    StdTimer& operator=(const StdTimer&) = default;

    void start() override;
    void stop() override;
    double elapsedTime() const override;
    std::string summary() const override;

    std::time_t itsElapsedTime;
    std::time_t itsStartTime;
    std::time_t itsStopTime;
};

#ifdef HAVE_MPI
class MPITimer final : public ITimer
{
  public:
    friend class Timer;
    ~MPITimer() {};
  private:
    MPITimer();
    MPITimer(const MPITimer&) = default;
    MPITimer& operator=(const MPITimer&) = default;
    void start() override;
    void stop() override;
    double elapsedTime() const override;
    std::string summary() const override;
    double itsElapsedTime;
    double itsStartTime;
    double itsStopTime;
    askap::utils::State itsState;
};
#endif

class SectionTimer final
{
  public:
    SectionTimer(unsigned int sections);
    SectionTimer(const SectionTimer&) = delete;
    SectionTimer& operator=(const SectionTimer&) = delete;
    ~SectionTimer() {}

    void start(unsigned int section);
    void stop(unsigned int section);
    void summary() const;
    double totalElapsedTime() const;
  private:
    std::map<unsigned int, std::shared_ptr<Timer>> itsTimers;

};
}
}
#endif
