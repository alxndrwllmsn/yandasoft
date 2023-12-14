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

/// @brief generate FFTW3 wisdom for range of sizes of 2D Complex transform
/// @param[in] int minSize : smallest size FFT to plan
/// @param[in] int maxSize : largest size FFT to plan
/// @param[in] int numThreads : number of threads to plan for
/// @param[in] uint planOption : FFTW plan option (FFTW_MEASURE, FFTW_PATIENT, or FFTW_EXHAUSTIVE)
void generateFloatWisdom(int minSize, int maxSize, int numThreads, uint planOption)
{
    casacore::Matrix<casacore::Complex> data(maxSize,maxSize);
    fftwf_complex* buf = reinterpret_cast<fftwf_complex*>(data.data());
    #ifdef HAVE_FFTW_OPENMP
    ASKAPCHECK(fftwf_init_threads(),"Failure initialising fftw threads");
    fftwf_plan_with_nthreads(numThreads);
    #endif
    int n = minSize;
    ASKAPLOG_INFO_STR(logger,"Planning 2D Complex FFT sizes: "<<n <<" to "<<maxSize<<" with "<<numThreads<<" threads");
    const char* env = std::getenv("FFTW_WISDOM");
    ASKAPCHECK(env,"Environment variable FFTW_WISDOM pointing to wisdom directory needs to be defined");
    string file(env);
    file += "/fftwf-threads-"+utility::toString(numThreads)+".wisdom";
    int status = fftwf_import_wisdom_from_filename(file.c_str());
    if (status == 0) {
        ASKAPLOG_WARN_STR(logger,"Error reading FFTW wisdom file: " << file);
    }
    while (n <= maxSize) {
        fftwf_plan plan1 = fftwf_plan_dft_2d(n, n, buf, buf, FFTW_FORWARD, planOption);
        fftwf_plan plan2 = fftwf_plan_dft_2d(n, n, buf, buf, FFTW_BACKWARD, planOption);
        n = askap::synthesis::SynthesisParamsHelper::nextFactor2357(n + 1);
        // save results regularly
        if (n % 1024 == 0) {
            ASKAPCHECK(fftwf_export_wisdom_to_filename(file.c_str()),"Error writing wisdom to file "<< file);
            ASKAPLOG_INFO_STR(logger, "Planning complete up to "<< n);
        }
    }
    ASKAPCHECK(fftwf_export_wisdom_to_filename(file.c_str()),"Error writing wisdom to file "<< file);
    ASKAPLOG_INFO_STR(logger,"Done - wisdom exported to "<< file);
}

/// @brief generate FFTW3 wisdom for range of sizes of 2D Double Complex transform
/// @param[in] int minSize : smallest size FFT to plan
/// @param[in] int maxSize : largest size FFT to plan
/// @param[in] int numThreads : number of threads to plan for
/// @param[in] uint planOption : FFTW plan option (FFTW_MEASURE, FFTW_PATIENT, or FFTW_EXHAUSTIVE)
void generateDoubleWisdom(int minSize, int maxSize, int numThreads, uint planOption)
{
    casacore::Matrix<casacore::DComplex> data(maxSize,maxSize);
    fftw_complex* buf = reinterpret_cast<fftw_complex*>(data.data());
    #ifdef HAVE_FFTW_OPENMP
    ASKAPCHECK(fftw_init_threads(),"Failure initialising fftw threads");
    fftw_plan_with_nthreads(numThreads);
    #endif
    int n = minSize;
    ASKAPLOG_INFO_STR(logger,"Planning 2D DComplex FFT sizes: "<<n <<" to "<<maxSize<<" with "<<numThreads<<" threads");
    const char* env = std::getenv("FFTW_WISDOM");
    ASKAPCHECK(env,"Environment variable FFTW_WISDOM pointing to wisdom diectory needs to be defined");
    string file(env);
    file += "/fftw-double-"+utility::toString(numThreads)+"threads.wisdom";
    const int status = fftw_import_wisdom_from_filename(file.c_str());
    if (status == 0) {
        ASKAPLOG_WARN_STR(logger,"Error reading FFTW wisdom file: " << file);
    }
    while (n <= maxSize) {
        fftw_plan plan1 = fftw_plan_dft_2d(n, n, buf, buf, FFTW_FORWARD, planOption);
        fftw_plan plan2 = fftw_plan_dft_2d(n, n, buf, buf, FFTW_BACKWARD, planOption);
        n = askap::synthesis::SynthesisParamsHelper::nextFactor2357(n + 1);
        // save results regularly
        if (n % 1024 == 0) {
            ASKAPCHECK(fftw_export_wisdom_to_filename(file.c_str()),"Error writing wisdom to file "<< file);
            ASKAPLOG_INFO_STR(logger, "Planning complete up to "<< n);
        }
    }
    ASKAPCHECK(fftw_export_wisdom_to_filename(file.c_str()),"Error writing wisdom to file "<< file);
    ASKAPLOG_INFO_STR(logger,"Done - wisdom exported to "<< file);
}


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
