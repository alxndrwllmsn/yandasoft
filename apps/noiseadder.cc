//
// @file noiseadder.cc : Small utility to add random noise to the existing measurement set
//
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


#include <askap/dataaccess/TableDataSource.h>
#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".noiseadder");

#include <askap/askap/AskapError.h>
#include <askap/dataaccess/SharedIter.h>
#include <askap/dataaccess/ParsetInterface.h>

#include <askap/dataaccess/TableManager.h>
#include <askap/dataaccess/IDataConverterImpl.h>

// casa
#include <casacore/measures/Measures/MFrequency.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/casa/OS/Timer.h>
#include <casacore/casa/Arrays/ArrayMath.h>

#include <askap/utils/CommandLineParser.h>

#include <askap/measurementequation/GaussianNoiseME.h>
#include <askap/measurementequation/IMeasurementEquation.h>
#include <askap/dataaccess/MemBufferDataAccessor.h>
#include <askap/askapparallel/MPIComms.h>

// std
#include <stdexcept>
#include <iostream>

using std::cout;
using std::cerr;
using std::endl;

using namespace askap;
using namespace synthesis;
using namespace askap::accessors;


void addNoise(const IDataSource &ds, const IMeasurementEquation &ime) {
  IDataSelectorPtr sel=ds.createSelector();
  //sel->chooseFeed(1);  
  IDataConverterPtr conv=ds.createConverter();
  conv->setFrequencyFrame(casa::MFrequency::Ref(casa::MFrequency::TOPO),"MHz");
  conv->setEpochFrame(casa::MEpoch(casa::Quantity(53635.5,"d"),
                      casa::MEpoch::Ref(casa::MEpoch::UTC)),"s");
  IDataSharedIter it=ds.createIterator(sel,conv);
  //for (size_t run=0;run<10;++run)
  for (it.init();it!=it.end();it.next()) {
       MemBufferDataAccessor buffer(*it);
       ime.predict(buffer);
       it->rwVisibility() += buffer.visibility();
  }
}

int main(int argc, char **argv) {
  try {
     casa::Timer timer;

     timer.mark();
     
     utils::Parser parser; // a command line parser
     
     // command line parameters
     utils::GenericParameter<std::string> msName;
     utils::GenericParameter<double> noiseVariance;
     parser.add(msName, utils::Parser::throw_exception);
     parser.add(noiseVariance, utils::Parser::throw_exception);
     parser.process(argc, argv);
     
     // Initialize MPI (also succeeds if no MPI available).
     askap::askapparallel::MPIComms comms(argc, argv);

     casa::Int seed1 = casa::Int(time(0));
     casa::Int seed2 = casa::Int(comms.rank());
     std::cerr<<"Using seeds: "<<seed1<<" "<<seed2<<std::endl;
     GaussianNoiseME noiseME(noiseVariance,seed1,seed2);
     
     TableDataSource ds(msName.getValue(),TableDataSource::MEMORY_BUFFERS | TableDataSource::WRITE_PERMITTED);     
     
     std::cerr<<"Initialization: "<<timer.real()<<std::endl;

     timer.mark();
     addNoise(ds,noiseME);
     std::cerr<<"Job: "<<timer.real()<<std::endl;
     
  }
  catch(const utils::XParser &) {
     cerr<<"Usage "<<argv[0]<<" measurement_set noise_variance_in_Jy_squared"<<endl;
	 return -2;    
  }
  catch(const AskapError &ce) {
     cerr<<"AskapError has been caught. "<<ce.what()<<endl;
     return -1;
  }
  catch(const std::exception &ex) {
     cerr<<"std::exception has been caught. "<<ex.what()<<endl;
     return -1;
  }
  catch(...) {
     cerr<<"An unexpected exception has been caught"<<endl;
     return -1;
  }
  return 0;
}
