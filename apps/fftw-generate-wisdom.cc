/// @file fftw-generate-wisdom.cc
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

#include <casacore/casa/Arrays/Matrix.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/askap/AskapError.h>
#include <askap/askap/AskapLogging.h>
#include <askap/scimath/fft/FFT2DWrapper.h>
#include <askap/utils/CommandLineParser.h>
#include <Common/OpenMP.h>
#include <Common/ParameterSet.h>


#include <iostream>
#include <fftw3.h>

using namespace std;
using namespace askap;

ASKAP_LOGGER(logger, "synthesis.fftwWisdom");

// Main function
// Generate fftw3 wisdom for 2d complex to complex single or double precision transforms
// fft sizes used in pipelines: 1024, 1536, 2048, 2560, 4000, 4096, 5120, 6144 + nyquist grids
// Code below generates wisdom for all sizes from minSize up to maxSize, for sizes that factorise
// into powers of 2, 3, 5 and / or 7.
int main(int argc, const char** argv) {
    // Now we have to initialize the logger before we use it
    // If a log configuration exists in the current directory then
    // use it, otherwise try to use the programs default one
    std::ifstream config("askap.log_cfg", std::ifstream::in);
    if (config) {
        ASKAPLOG_INIT("askap.log_cfg");
    } else {
        std::ostringstream ss;
        ss << argv[0] << ".log_cfg";
        ASKAPLOG_INIT(ss.str().c_str());
    }

    utils::Parser parser;
    utils::FlaggedParameter<std::string> parsetPar("-c","");
    // this parameter is optional
    parser.add(parsetPar, utils::Parser::return_default);
    parser.process(argc, argv);
    const std::string parsetFile = parsetPar;
    const LOFAR::ParameterSet parset = (parsetPar.defined() ? LOFAR::ParameterSet(parsetPar) : LOFAR::ParameterSet());

    const int minSize = parset.getInt("minsize",1024);
    const int maxSize = parset.getInt("maxsize",6144);
    const int numThreads = parset.getInt("numthreads",1);
    const bool doDouble = parset.getBool("double",false);
    // plan options are FFTW_ESTIMATE, FFTW_MEASURE, FFTW_PATIENT, or FFTW_EXHAUSTIVE
    // ESTIMATE is the default without wisdom
    // time for planning single 5000^2 FFT: Measure: ~45s, Patient: ~45m, exhaustive: ?
    const string plan = parset.getString("plan","measure");
    uint planOption = FFTW_MEASURE;
    if (plan == "patient") {
        planOption = FFTW_PATIENT;
    } else if (plan == "exhaustive") {
        planOption = FFTW_EXHAUSTIVE;
    } else if (plan != "measure") {
        ASKAPTHROW(AskapError, "Unknown plan type (use measure, patient or exhaustive) :"<< plan);
    }
    ASKAPCHECK(numThreads <= LOFAR::OpenMP::maxThreads(),"not enough threads available");

    if (doDouble) {
        scimath::FFT2DWrapper<casacore::DComplex> fft2d;
        fft2d.generateWisdom(minSize, maxSize, numThreads, planOption);
    } else {
        scimath::FFT2DWrapper<casacore::Complex> fft2d;
        fft2d.generateWisdom(minSize, maxSize, numThreads, planOption);
    }
}
