/// @file cflagger.cc
///
/// Application to flag data
/// Control parameters are passed in from a LOFAR ParameterSet file.
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
/// @author Mark Wieringa

// Package level header file
#include <askap/askap_synthesis.h>

// ASKAPsoft includes
#include <askap/askap/AskapLogging.h>
#include <askap/askap/AskapError.h>
#include <askap/askap/Application.h>
#include <askap/askap/StatReporter.h>
#include <askap/askapparallel/AskapParallel.h>
#include <Common/ParameterSet.h>

// Local package includes
#include <askap/parallel/FlagParallel.h>

ASKAP_LOGGER(logger, ".cflagger");

using namespace std;
using namespace askap;
using namespace askap::synthesis;

class CFlaggerApp : public askap::Application
{
    public:
        virtual int run(int argc, char* argv[]) override
        {
            StatReporter stats;

            // This class must have scope outside the main try/catch block
            askap::askapparallel::AskapParallel comms(argc, const_cast<const char**>(argv));

            try {
                LOFAR::ParameterSet subset(config().makeSubset("cflagger."));
                if (subset.empty()) {
                    subset = config().makeSubset("Cflag.");
                }

                // Perform %w or %r substitutions for all keys.
                if (comms.isParallel()) {
                    for (LOFAR::ParameterSet::iterator it = subset.begin(); it != subset.end(); ++it) {
                        it->second = LOFAR::ParameterValue(comms.substitute(it->second));
                    }
                }
                // Make the master do work
                subset.replace(LOFAR::KVpair("masterDoesWork",true));
                
                // We cannot issue log messages until MPI is initialized!
                FlagParallel flag(comms, subset);

                ASKAPLOG_INFO_STR(logger, "ASKAP synthesis flagging application " << ASKAP_PACKAGE_VERSION);

                if (comms.isMaster()) {
                    ASKAPLOG_INFO_STR(logger, "Parset file contents:\n" << config());
                }

                flag.doFlag();
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
    CFlaggerApp app;
    return app.main(argc, argv);
}
