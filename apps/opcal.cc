/// @file
///
/// @brief Generic operations-specific calibration
/// @details This application is intended for the types of calibration which 
/// cannot follow ASKAP's predict-forward approach, i.e. which require 
/// observations of various fields done in some special way. It is intended
/// for experimentation with calibration as well as some operation-specific
/// tasks like baseline and pointing refinements.  
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
/// @author Max Voronkov <Maxim.Voronkov@csiro.au>
// Package level header file
#include <askap/askap_synthesis.h>

// ASKAPsoft includes
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>
#include <askap/askap/Application.h>
#include <askap/askap/StatReporter.h>
#include <askap/askapparallel/AskapParallel.h>
#include <Common/ParameterSet.h>
#include <askap/opcal/OpCalImpl.h>
#include <askap/opcal/BaselineSolver.h>
#include <askap/opcal/PointingSolver.h>


ASKAP_LOGGER(logger, ".opcal");

using namespace std;
using namespace askap;
using namespace askap::synthesis;

class OpCalApp : public askap::Application
{
    public:
        virtual int run(int argc, char* argv[])
        {
            StatReporter stats;

            // This class must have scope outside the main try/catch block
            askap::askapparallel::AskapParallel comms(argc, const_cast<const char**>(argv));

            try {

                 ASKAPLOG_INFO_STR(logger, "ASKAP synthesis opcal application " << ASKAP_PACKAGE_VERSION);
               
                 if (comms.isMaster()) {
                   ASKAPLOG_INFO_STR(logger, "Parset file contents:\n" << config());
                 }
                 
                 OpCalImpl impl(comms,config());
                 // we'll eventually get a proper factory method here; for now create explicitly
                 const std::string hlSolver = config().getString("solver","");
                 if (hlSolver != "") {
                     ASKAPLOG_INFO_STR(logger, "Calibration results will be passed to '"<<hlSolver<<"' solver");
                     if (hlSolver == "baseline") {
                         impl.setHighLevelSolver(boost::shared_ptr<BaselineSolver>(new BaselineSolver(config())));
                     } else if (hlSolver == "pointing") {
                         impl.setHighLevelSolver(boost::shared_ptr<PointingSolver>(new PointingSolver(config())));
                     } else {
                       ASKAPCHECK(false, "Unknown high-level solver: "<<hlSolver);
                     }
                 } else {
                     ASKAPLOG_INFO_STR(logger, "High-level solver has not been defined, raw calibration results will be printed in the log");
                 }
                 
                 impl.run();                 
                 
            } catch (const askap::AskapError& e) {
                ASKAPLOG_FATAL_STR(logger, "Askap error in " << argv[0] << ": " << e.what());
                std::cerr << "Askap error in " << argv[0] << ": " << e.what() << std::endl;
                exit(1);
            } catch (const std::exception& e) {
                ASKAPLOG_FATAL_STR(logger, "Unexpected exception in " << argv[0] << ": " << e.what());
                std::cerr << "Unexpected exception in " << argv[0] << ": " << e.what() << std::endl;
                exit(1);
            }
            stats.logSummary();
            return 0;
        }

    private:
        std::string getVersion() const override {
            const std::string pkgVersion = std::string("yandasoft:") + ASKAP_PACKAGE_VERSION;
            return pkgVersion;
        }
};

// Main function
int main(int argc, char* argv[])
{
    OpCalApp app;
    return app.main(argc, argv);
}
            
