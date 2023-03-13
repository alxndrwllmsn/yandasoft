/// @file linmos-mpi.cc
///
/// @brief combine a number of images as a linear mosaic
/// @details This is a utility to merge images into a mosaic. Images can be set
/// explicitly or found automatically based on input tags.
///
/// @copyright (c) 2012,2014,2021 CSIRO
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

// Package level header file
#include <askap/askap_synthesis.h>
#include <askap/utils/StatsAndMask.h>

/// ASKAP includes
#include <askap/askap/AskapUtil.h>
#include <askap/utils/LinmosUtils.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/askapparallel/AskapParallel.h>
#include <askap/imageaccess/WeightsLog.h>
#include <askap/imagemath/linmos/LinmosAccumulator.h>
#include <askap/imagemath/utils/MultiDimArrayPlaneIter.h>

/// 3rd party
#include <Common/ParameterSet.h>

#include <map>

ASKAP_LOGGER(logger, ".linmos");

using namespace askap::synthesis;
using namespace std;
using namespace casacore;

namespace askap {

using ImageBlcTrcMapT =  std::map<std::string,std::pair<casacore::IPosition,casacore::IPosition>>;
// the key of  the map is channel no
// the pair is <channel weight pixel, max channel weight pixel>
using MaxWgtPerChannelMapT = std::map<int,float>;

/// @brief This function returns the weight pixel of a channel
/// @param[in] accumulator - object that provides functions which are heavily utilised by this task
/// @param[in] iacc - image access object
/// @param[in] inImgName - name of the input image
/// @param[in] inWgtName - name of the weight image
/// @param[in] channel - the channel of the weight pixel
/// @param[out] inWgtPix - the weight pixel
static void loadWgtImage(imagemath::LinmosAccumulator<float>& accumulator,
                         const accessors::IImageAccess<casacore::Float>& iacc,
                         const std::string& inImgName,
                         const std::string& inWgtName,
                         const int channel,
                         float& inWgtPix)
{

  if (accumulator.weightType() == FROM_WEIGHT_IMAGES || accumulator.weightType() == COMBINED) {
    if (accumulator.useWeightsLog()) {
      //ASKAPLOG_INFO_STR(logger,"Reading weights log file :"<< inWgtName);
      accessors::WeightsLog wtlog(inWgtName);
      wtlog.read();
      //inWgtPix.resize(inPix.shape());
      inWgtPix = wtlog.weight(channel);
    } else if (accumulator.useWtFromHdr()) {
      //ASKAPLOG_INFO_STR(logger,"Reading weights from input image file :"<< inImgName);
      const casacore::Vector<casacore::Float> wts = readWeightsTable(inImgName);
      const size_t size = wts.size();
      if (size > 1) {
        //inWgtPix.resize(inPix.shape());
        ASKAPCHECK(channel < size,"Not enough channels in weights table");
        inWgtPix = wts(channel);
      } else if (size == 1) {
        //inWgtPix.resize(inPix.shape());
        inWgtPix = wts(0);
      }
      ASKAPCHECK(size>0,"No weights found in image header or extension for image: "<<inImgName);
    } else {
      //inWgtPix = iacc.read(inWgtName);
      ASKAPCHECK(false,"Neither useweightslog or useWtFromHdr is set");
    }
  }
}
/// @brief This function returns the blc and trc of the input image
/// @detail if trimming is true, this function returns blc and trc of the 
///         trimmed input image otherwise it returns the full input image size.
///         Also, if trimming is true, it uses the imgName argument to look up
///         the blc and trc in the map otherwise it uses the imgInputOrWgtOrSenName
///         argument to obtain the full image size.
/// @param[in] trimming - flag to indicate if trimming is used or not
/// @param[in] iacc - image access object
/// @param[in] imgName - name of input image
/// @param[in] imageBlcTrcMap - a map which contains the input images' trc and blc
/// @param[out] blc - bottom left corner of the output image
/// @param[out] trc - top right corner of the output image
static void 
getBlcTrc(const bool trimming,
          const accessors::IImageAccess<casacore::Float>& iacc,
          const std::string& imgName,
          const ImageBlcTrcMapT& imageBlcTrcMap,
          casacore::IPosition& blc,
          casacore::IPosition& trc)
{
  if ( ! trimming ) {
    const casa::IPosition shape = iacc.shape(imgName);
    blc = casa::IPosition(shape.nelements(),0);
    trc = casa::IPosition(shape-1);
  } else {
    const auto imageBlcTrcMapIter = imageBlcTrcMap.find(imgName);
    const auto& blcTrcPair = imageBlcTrcMapIter->second;
    blc = blcTrcPair.first;
    trc = blcTrcPair.second;;
  }
}

/// @brief This function calculates the maximum weight pixel of each channel
/// @details This function iterates over the weight files and determines the
///          maximum weight pixel for each channel. The result is stored in
///          a map (maxWgtPerChannelMap). The key of the map is the channel
///          and the value is the max weight pixel.
/// @param[in] comms - mpi communicator
/// @param[in] accumulator - object that provides functions which are heavily utilised by this task
/// @param[in] iacc - image access object
/// @param[in] inImgNames - a list of input images
/// @param[in] inWgtNames - a list of weight images
/// @param[in] channel - the channel of the weight pixel
/// @param[out] maxWgtPerChannelMap - contains a map of channel and max weight pixel.
static void calcMaxWgtPerChannels(askap::askapparallel::AskapParallel &comms,
                               imagemath::LinmosAccumulator<float>& accumulator,
                               const accessors::IImageAccess<casacore::Float>& iacc,
                               const std::vector<std::string>& inImgNames,
                               const std::vector<std::string>& inWgtNames,
                               const int channel,
                               MaxWgtPerChannelMapT& maxWgtPerChannelMap)
{

  int inWgtNamesIndex = 0;

  // contains a collection of weight pixel values
  std::vector<float> wgtPixVect;
  for (vector<string>::const_iterator it = inImgNames.begin(); it != inImgNames.end(); ++it) {
    // find the max weight pixel for this each input image
    float inWgtPix = -1.0;
    loadWgtImage(accumulator,iacc,*it,inWgtNames[inWgtNamesIndex],channel,inWgtPix);
    wgtPixVect.push_back(inWgtPix);
    if ( comms.isMaster() ) {
        ASKAPLOG_INFO_STR(logger,"weightlog file: " << inWgtNames[inWgtNamesIndex] << ": channel" << channel << " , weight value: " << inWgtPix);
    }
    inWgtNamesIndex += 1;
  }    
  // get the max weight pixel for this rank
  auto iter = std::max_element(wgtPixVect.begin(),wgtPixVect.end());
  float thisRankMaxWgt = *iter;
  
  if ( comms.isMaster() ) {
    ASKAPLOG_INFO_STR(logger,"Master rank max weight pixel = "  << thisRankMaxWgt);
  }

  comms.barrier();

  float thisChannelMaxPix = 0;
  // get the max weight pixel for this channel
  wgtPixVect.resize(0); // just reuse this vector to store max weight pixel value for each rank
  int nProcs = comms.nProcs();

  if ( comms.isMaster() ) {
    wgtPixVect.push_back(thisRankMaxWgt);
    float buffer;
    for (int sender = 1; sender < nProcs; sender++) {
      comms.receive(&buffer,sizeof(float),sender);
      wgtPixVect.push_back(buffer);
    }
    auto maxIter = std::max_element(wgtPixVect.begin(),wgtPixVect.end());
    thisChannelMaxPix = *maxIter;
    // now send the max weight pixel of this image to all the workers
    for (int rcv = 1; rcv < nProcs; rcv++) {
      comms.send(&thisChannelMaxPix,sizeof(float),rcv);
      wgtPixVect.push_back(buffer);
    }
  } else {
    // send this rank max weight pixel to master
    comms.send((void *) &thisRankMaxWgt,sizeof(int),0);
    // wait and receive the max weight pixel of this weight image from master
    comms.receive((void *) &thisChannelMaxPix,sizeof(int),0);
  }
  if ( comms.isMaster() ) {
    ASKAPLOG_INFO_STR(logger,"calcMaxWgtPerChannels: channel" << channel << ", thisChannelMaxPix: " << thisChannelMaxPix);
  }
  maxWgtPerChannelMap.insert(std::make_pair(channel,thisChannelMaxPix));
  comms.barrier();
}
                                            
/// @brief This function calculates the blc and trc  of the input image
/// @details This function iterates over the channels and calculates the blc and trc
///          of the input image. The imageBlcTrcMap map stores the blc and trc of the images.
///          The key of the map is the input name and the value is an std::pair consisting of
///          blc and trc values.
/// @param[in] parset - the task configuration settings
/// @param[in] comms - mpi communicator
/// @param[in] accumulator - object that provides functions which are heavily utilised by this task
/// @param[in] iacc - image access object
/// @param[in] inImgName - name of input image
/// @param[in] firstChannel - first channel allocated to this rank
/// @param[in] lastChannel - last channel allocated to this rank
/// @param[in] channelInc - channel increment
/// @param[in] nchanCube - numer of channels in the cube
/// @param[in] beanmCentreIndex - index of the input images
/// @param[out] imageBlcTrcMap - contains the blc and trc of the input images
/// @param[in] useWgtLog - true indicate if the weight log file is used
static void calcMinMaxXYInputImagePlanes(const LOFAR::ParameterSet &parset,
                                         askap::askapparallel::AskapParallel &comms,
                                         imagemath::LinmosAccumulator<float>& accumulator,
                                         const accessors::IImageAccess<casacore::Float>& iacc,
                                         const std::string& inImgName,
                                         const std::string& inWgtName,
                                         const int firstChannel, const int lastChannel,
                                         const int channelInc, const int nchanCube,
                                         const int numberOfInImg,
                                         const int beamCentreIndex,
                                         MaxWgtPerChannelMapT& maxWgtPerChannelMap,
                                         ImageBlcTrcMapT& imageBlcTrcMap,
                                         const bool useWgtLog)
{
  const float cutoff = parset.getFloat("cutoff",0.01);
  float beamPackingFactor = 1.0;

  if ( parset.isDefined("beampackingfactor") ) {
    beamPackingFactor = parset.getFloat("beampackingfactor");
  } else if (numberOfInImg >= 36) {
    beamPackingFactor = 1.4;
  }
  
  // used to store the min and max of the x and y dimension of the input image planes/channels
  std::vector<int> xMinMaxVect;
  std::vector<int> yMinMaxVect;

  // calculates the maximum weight pixel for each channel
  for (int channel = firstChannel; channel <= lastChannel; channel += channelInc) {
    float scaledCutoff = cutoff;
    if (useWgtLog) {
      float inWgtPix = 0.0;
      loadWgtImage(accumulator,iacc,inImgName,inWgtName,channel,inWgtPix);
      ASKAPCHECK(maxWgtPerChannelMap.find(channel) != maxWgtPerChannelMap.end() || maxWgtPerChannelMap[channel] != 0.0,
               "Channel: " << channel << " is not found or it has 0 max weight pixel");
      const float maxWgt = maxWgtPerChannelMap[channel];
      scaledCutoff = (inWgtPix == 0 ? 0 : cutoff * beamPackingFactor * sqrt(maxWgt/inWgtPix));
      if ( comms.isMaster() ) {
        ASKAPLOG_INFO_STR(logger,"channel: " << channel 
                        << "; scaledCutoff: " << scaledCutoff
                        << "; inWgtPix: << " << inWgtPix
                        << "; maxWgt: " << maxWgt);
      }
    }

    int xmin = -1;
    int xmax = -1;
    int ymin = -1;
    int ymax = -1;
    if ( scaledCutoff != 0 ) {
      // scaledCutoff is 0 if the weight pixel of the channel is 0 (ie no valid data).
      // Hence, we dont need to calculate min and max for this channel
      CoordinateSystem coordSys = iacc.coordSys(inImgName);
      casacore::IPosition shape = iacc.shape(inImgName);
      accumulator.calcWeightInputShape(coordSys,shape,beamCentreIndex,
                                       channel,scaledCutoff,xmin,xmax,ymin,ymax);
      // Dont reset these vectors to 0. Just add xmin, xmax, ymin and ymax to them
      xMinMaxVect.push_back(xmin);
      xMinMaxVect.push_back(xmax);
      yMinMaxVect.push_back(ymin);
      yMinMaxVect.push_back(ymax);
    }
  }

  // xMinMaxVect and yMinMaxVect now contains the min and max of all the image planes
  // of the input image
  int smallestX = -1.0; int largestX = -1.0; int smallestY = -1.0; int largestY = -1.0;

  if (xMinMaxVect.size() != 0 && yMinMaxVect.size() != 0) {
    auto xMinMaxIter = std::minmax_element(xMinMaxVect.begin(),xMinMaxVect.end());
    smallestX = *xMinMaxIter.first;
    largestX = *xMinMaxIter.second;

    auto yMinMaxIter = std::minmax_element(yMinMaxVect.begin(),yMinMaxVect.end());
    smallestY = *yMinMaxIter.first;
    largestY = *yMinMaxIter.second;
  }

  // smallestX, largestX, smallestY, largestY are the smallest and largest min max
  // values of the x and y dimension of the image planes/channels allocated to this rank
  // What we need to do now is to collect these values from all the ranks and work out
  // the smallestX, largestX, smallestY, largestY from those collected values and this will
  // give us the smallestX, largestX, smallestY, largestY of the input image

  // reset xMinMaxVect and yMinMaxVect
  xMinMaxVect.resize(0);
  yMinMaxVect.resize(0);
  // xMinMaxVect and yMinMaxVect contains the smallest and largest min and max values belonged
  // to this rank
  if ( smallestX != -1 && largestX != -1 && smallestY != -1 && largestY != -1 ) {
    // only interested in channels that have valid smallest and largest x and y 
    xMinMaxVect.push_back(smallestX);
    xMinMaxVect.push_back(largestX);
    yMinMaxVect.push_back(smallestY);
    yMinMaxVect.push_back(largestY);
  }

  int xyMinMax[4];
  int nProcs = comms.nProcs();

  // now find the smallest and largest min max x y values of all the ranks.
  // This gives us the smallest blc and largest trc of the input image
  int thisImgSmallestX = 0;
  int thisImgSmallestY = 0;
  int thisImgLargestX = 0;
  int thisImgLargestY = 0;

  if ( comms.isMaster() ) {
    // then collect min and max of x and y of workers to master
    for (int sender = 1; sender < nProcs; sender++) {
      comms.receive(xyMinMax,4*sizeof(int),sender);
      if ( xyMinMax[0] != -1.0 && xyMinMax[1] != -1.0 && 
           xyMinMax[2] != -1.0 && xyMinMax[3] != -1.0 ) {
        // now copy from worker ranks
        xMinMaxVect.push_back(xyMinMax[0]);
        xMinMaxVect.push_back(xyMinMax[2]);
        yMinMaxVect.push_back(xyMinMax[1]);
        yMinMaxVect.push_back(xyMinMax[3]);
      }
    }

    // The master now has collected all the min and max from all the workers.
    // Next, work out the smallest and largest from all the ranks
    auto xMinMaxIter = std::minmax_element(xMinMaxVect.begin(),xMinMaxVect.end());
    thisImgSmallestX = *xMinMaxIter.first;
    thisImgLargestX = *xMinMaxIter.second;

    auto yMinMaxIter = std::minmax_element(yMinMaxVect.begin(),yMinMaxVect.end());
    thisImgSmallestY = *yMinMaxIter.first;
    thisImgLargestY = *yMinMaxIter.second;   

    // thisImgSmallestX, thisImgSmallestY, thisImgLargestX and thisImgLargestY
    // are the smallest and largest min and max of this image.
    xyMinMax[0] = thisImgSmallestX;
    xyMinMax[1] = thisImgSmallestY;
    xyMinMax[2] = thisImgLargestX;
    xyMinMax[3] = thisImgLargestY;
    // now send the smallest and largest x and y of this image all the workers
    for (int receiver = 1; receiver < nProcs; receiver++) {
      comms.send((void *) xyMinMax,4*sizeof(int),receiver);
    }
  } else {
    // send min and max of x and y to the master
    xyMinMax[0] = smallestX;
    xyMinMax[1] = smallestY;
    xyMinMax[2] = largestX;
    xyMinMax[3] = largestY;
    comms.send((void *) xyMinMax, 4*sizeof(int),0);
    //  wait and receive the smallest and largest x and y of this image from master
    comms.receive(xyMinMax,4*sizeof(int),0);

    thisImgSmallestX = xyMinMax[0];
    thisImgSmallestY = xyMinMax[1];
    thisImgLargestX = xyMinMax[2];
    thisImgLargestY = xyMinMax[3];
  }    
  comms.barrier();
  
  // by the time the code gets here, each rank knows the smallest and largest min and max
  // of the image. Dedue the blc and trc from these values and store them in the map
  casacore::IPosition blc(4,thisImgSmallestX,thisImgSmallestY,0,nchanCube-1);
  casacore::IPosition trc(4,thisImgLargestX,thisImgLargestY,0,nchanCube-1);
  std::pair<casacore::IPosition,casacore::IPosition> p;
  p.first = blc;
  p.second = trc;
  imageBlcTrcMap.insert(std::make_pair(inImgName,p));
}

/// @brief This function calculates the bounding box of pixels above the beam cutoff for all 
///        the image cubes.
/// @detail This function calculates required blc and trc for the output image. This is 
///         achieved by each rank first determines its own smallest and largest 
///         min and max for the x and y dimension. The master (rank 0) then collects
///         these values from other ranks and deduces the required blc and trc for the cube.
///         Finally, the blc and trc are sent to other ranks
/// @param[in] parset - the task configuration settings
/// @param[in] comms - mpi communicator object
/// @param[in] accumulator - object that provides functions which are heavily utilised by this task
/// @param[in] iacc - image access object
/// @param[in] inImgNames - input images
/// @param[in] nchanCube - how many channels in the cube
/// @param[in] firstChannel - first channel of this rank
/// @param[in] lastChannel - last channel of this rank
/// @param[in] channelInc - channel increment
/// @param[out] imageBlcTrcMap - contains the blc and trc of the output limos images.
///                              It is used to set the size of the output images.
/// @param[out] imageBlcTrcLimosShapeMap - contains the blc and trc of the input images
///                              It is used to set the size of the input images.
static void findBoundingBoxes(const LOFAR::ParameterSet &parset,
                            askap::askapparallel::AskapParallel &comms,
                            imagemath::LinmosAccumulator<float>& accumulator,
                            accessors::IImageAccess<casacore::Float>& iacc,
                            const std::vector<string>& inImgNames,
                            const std::vector<string>& inWgtNames,
                            const int nchanCube, const int firstChannel, 
                            const int lastChannel, const int channelInc,
                            ImageBlcTrcMapT& imageBlcTrcMap,
                            ImageBlcTrcMapT& imageBlcTrcLimosShapeMap)
{
  ASKAPLOG_INFO_STR(logger,"findBoundingBoxes");

  // find the max weight pixel for each channel
  // const bool useWgtLog = parset.getBool("useweightslog", false);
  const bool useWgtLog = (accumulator.useWeightsLog() || accumulator.useWtFromHdr()); 
  MaxWgtPerChannelMapT maxWgtPerChannelMap;
  if ( useWgtLog ) {
    ASKAPLOG_INFO_STR(logger,"useWgtLog: true");
    for (int channel = firstChannel; channel <= lastChannel; channel += channelInc) {
      calcMaxWgtPerChannels(comms,accumulator,iacc,inImgNames,inWgtNames,channel,maxWgtPerChannelMap);
    }
    comms.barrier();
  } else {
    ASKAPLOG_INFO_STR(logger,"useWgtLog: false");
  }

  int beamCentreIndex = 0;
  int inWgtNamesIndex = 0;
  const int numberOfInImg = inImgNames.size();
  for (vector<string>::const_iterator it = inImgNames.begin(); it != inImgNames.end(); ++it) {
    calcMinMaxXYInputImagePlanes(parset,comms,accumulator,iacc,*it,inWgtNames[inWgtNamesIndex],
                                 firstChannel,lastChannel,channelInc,nchanCube,numberOfInImg,
                                 beamCentreIndex,maxWgtPerChannelMap,imageBlcTrcMap,useWgtLog);
    beamCentreIndex += 1;
    inWgtNamesIndex += 1;
  }

  // wait for all ranks get to here
  comms.barrier();

  if ( useWgtLog ) {
    // if weight log file is used, we need to get the blc and trc for the input images
    beamCentreIndex = 0;
    inWgtNamesIndex = 0;
    MaxWgtPerChannelMapT unusedMap;
    for (vector<string>::const_iterator it = inImgNames.begin(); it != inImgNames.end(); ++it) {
      calcMinMaxXYInputImagePlanes(parset,comms,accumulator,iacc,*it,inWgtNames[inWgtNamesIndex],
                                 firstChannel,lastChannel,channelInc,nchanCube,numberOfInImg,
                                 beamCentreIndex,unusedMap,imageBlcTrcLimosShapeMap,false);
      beamCentreIndex += 1;
      inWgtNamesIndex += 1;
    }

    comms.barrier();
  }
}
// @brief do the merge
/// @param[in] parset subset with parameters
static void mergeMPI(const LOFAR::ParameterSet &parset, askap::askapparallel::AskapParallel &comms) {

  ASKAPLOG_INFO_STR(logger, "ASKAP linear (parallel) mosaic task (MPI+LOWMEM)" << ASKAP_PACKAGE_VERSION);
  if (comms.isMaster()) {
      ASKAPLOG_INFO_STR(logger, "Parset parameters:\n" << parset);
  }

  imagemath::LinmosAccumulator<float> accumulator;

  // get the image history keyword if it is defined
  const std::vector<std::string> historyLines = parset.getStringVector("imageHistory",{},false);

  // get the calcstats flag from the parset. if it is true, then this task also calculates the image statistics
  const bool calcstats = parset.getBool("calcstats", false);
  // file to store the statistics
  const std::string outputStats = parset.getString("outputStats","");

  const bool useWgtLog = parset.getBool("useweightslog", false);

  // get the list of keywords to copy from the input
  // Set some defaults for simple beam mosaic, won't be correct for more complex cases
  const std::vector<std::string> keywordsToCopy = parset.getStringVector("keywordsToCopy",
      {"TELESCOP","PROJECT","SBID","DATE-OBS","DURATION"});

  // Original shape
  int nchanCube = -1;
  // load the parset
  if ( !accumulator.loadParset(parset) ) return;
  // initialise an image accessor
  accessors::IImageAccess<casacore::Float>& iacc = SynthesisParamsHelper::imageHandler();
  // Do we want to use collective I/O or individual/independent I/O
  bool collective = parset.getString("imageaccess","individual")=="collective";
  // CASA images need to be written one process at a time, for fits we have a choice
  // options "serial", "parallel"
  bool serialWrite = parset.getString("imageaccess.write","serial")=="serial" ||
                     parset.getString("imagetype")=="casa";


  // if we have Taylor terms and we need to correct them for the beam spectral
  // index - do it now ...


  // loop over the mosaics, reading each in and adding to the output pixel arrays
  vector<string> inImgNames, inWgtNames, inSenNames, inStokesINames;
  string outImgName, outWgtName, outSenName;
  map<string,string> outWgtNames = accumulator.outWgtNames();

  // Ahh this loops over the output mosaicks first
  // This will be just a single image or all taylor terms for an image (if findmosaics=false)

  for(map<string,string>::iterator ii=outWgtNames.begin(); ii!=outWgtNames.end(); ++ii) {

    // get output files for this mosaic
    outImgName = (*ii).first;
    outWgtName = accumulator.outWgtNames()[outImgName];
    ASKAPLOG_INFO_STR(logger, "++++++++++++++++++++++++++++++++++++++++++");
    ASKAPLOG_INFO_STR(logger, "Preparing mosaic " << outImgName);
    if (!accumulator.outWgtDuplicates()[outImgName]) {
      ASKAPLOG_INFO_STR(logger, " - also weights image " << outWgtName);
    }
    accumulator.doSensitivity(false);
    if (accumulator.genSensitivityImage()[outImgName]) {
      outSenName = accumulator.outSenNames()[outImgName];
      accumulator.doSensitivity(true);
      ASKAPLOG_INFO_STR(logger, " - also sensitivity image " << outSenName);
    }

    // get input files for this mosaic
    inImgNames = accumulator.inImgNameVecs()[outImgName];
    ASKAPLOG_INFO_STR(logger, "output mosaic " <<outImgName << " has input images: "<<inImgNames);

    if (accumulator.weightType() == FROM_WEIGHT_IMAGES || accumulator.weightType() == COMBINED ) {
      inWgtNames = accumulator.inWgtNameVecs()[outImgName];
      if (inWgtNames.size()==0) {
        // want weights, but have not specified any files - try to read from image headers
        bool useWtFromHdr = true;
        for (auto i : inImgNames) {
          useWtFromHdr = readWeightsTable(i).size()>0;
          if (!useWtFromHdr) {
            break;
          }
        }
        ASKAPCHECK(useWtFromHdr, "Image weighting requested, but no weights found in headers");
        ASKAPLOG_INFO_STR(logger,"   (reading weights from input images)");
        accumulator.setUseWtFromHdr(useWtFromHdr);
      } else if (accumulator.useWeightsLog()) {
          ASKAPLOG_INFO_STR(logger, " - input weightslog files: " << inWgtNames);
      } else {
          ASKAPLOG_INFO_STR(logger, " - input weights images: " << inWgtNames);
      }
    }

    if (accumulator.weightType() == FROM_BP_MODEL|| accumulator.weightType() == COMBINED) {
      accumulator.setBeamCentres(loadBeamCentres(parset,iacc,inImgNames));
    }

    if (accumulator.doSensitivity()) {
      inSenNames = accumulator.inSenNameVecs()[outImgName];
      ASKAPLOG_INFO_STR(logger, " - input sensitivity images: " << inSenNames);
    }

    if (accumulator.doLeakage()) {
        inStokesINames = accumulator.inStokesINameVecs()[outImgName];
        if (inStokesINames.size()>0) {
            ASKAPLOG_INFO_STR(logger, " - and using input Stokes I images: " << inStokesINames);
        }
    }

    // set the output coordinate system and shape, based on the overlap of input images
    // for this output image

    vector<IPosition> inShapeVec;
    vector<CoordinateSystem> inCoordSysVec;
    // if the useWgtLog is true, the vectors below store the coordinate systems and shapes of the
    // linmos input shape. This is needed because if the useWgtLog is true, the output shape and
    // the input shape (see setOutputParameters() and setInputParameters()) are not the same i.e
    // in this case, we want the input shape to be at beam cutoff whereas the output shape to be
    // at a scaled beam cutoff (1.4 * cutoff * sqrt(maxWgtPix/WgtPix) 
    vector<IPosition> inLinmosShapeVec;
    vector<CoordinateSystem> inLinmosCoordSysVec;
    // What fraction of the full problem does this rank
    // THe units here are channels - not polarisations.
    int numChannelsLocal = 0;
    int firstChannel = 0;
    int channelInc = 1;
    int lastChannel = 0;
    // Where a rank is in its allocation
    int channel = 0;
    //
    // this calculates the allocations for the images but stupidly does it for every image ...
    // anyway worry about that later
    //

    // I need to change this so that this is done for the current channel ....
    //
    // Lets get an example infile

    string testImage = inImgNames[0];

    //
    // lets get its shape

    const casa::IPosition shape = iacc.shape(testImage);

    ASKAPCHECK(shape.nelements()==4,"linmos-mpi can only work with 4D cubes: (RA,Dec,Pol,Freq)");
    ASKAPLOG_INFO_STR(logger," - ImageAccess Shape " << shape);

    // lets calculate the allocations ...

    casa::IPosition blc(shape.nelements(),0);
    casa::IPosition trc(shape);
    nchanCube = trc[3];

    if (comms.rank() >= nchanCube) {
      ASKAPLOG_WARN_STR(logger,"Rank " << comms.rank() << " has no work to merge");
      return;
    }
    if (nchanCube % comms.nProcs() != 0) {
      ASKAPLOG_WARN_STR(logger,"Unbalanced allocation: num of ranks:" << comms.nProcs() << " not a factor of number of channels: "<< trc[3]);
      if (collective)  ASKAPLOG_WARN_STR(logger,"Turning off collective I/O due to unbalanced allocation");
      collective = false;
    }
    if (comms.nProcs() > nchanCube) {
      collective = false; // can't do it easily
      numChannelsLocal = 1;
    }
    else {
      numChannelsLocal = nchanCube/comms.nProcs();
      if (nchanCube % comms.nProcs() != 0) numChannelsLocal +=1;
    }

    // Two schemes for reading/writing - not sure what is best for independent I/O
    // 1. process channels evenly distributed over file, each rank has its own section
    // 2. process contiguous sections of file, all ranks process each section
    bool contiguous = parset.getString("imageaccess.order","distributed")!="distributed";
    // For collective I/O only contiguous sections are currently supported - so enforce that
    contiguous |= collective;

    if (contiguous) {
        if (collective) ASKAPLOG_INFO_STR(logger,"Using collective I/O and contiguous channel access");
        else ASKAPLOG_INFO_STR(logger,"Using contiguous channel access");
        firstChannel = comms.rank();
        channelInc = comms.nProcs();
        lastChannel = firstChannel + (numChannelsLocal - 1) * channelInc;
        // unless beyond end
        if (lastChannel >= nchanCube) {
           numChannelsLocal -= 1;
           lastChannel = firstChannel + (numChannelsLocal - 1) * channelInc;
        }
    } else {
        ASKAPLOG_INFO_STR(logger,"Using distributed channel access");
        firstChannel = comms.rank() * numChannelsLocal;
        if (firstChannel >=nchanCube) {
            ASKAPLOG_INFO_STR(logger,"Rank " << comms.rank() << " has no work to merge");
            return;
        }
        channelInc = 1;
        lastChannel = firstChannel + numChannelsLocal - 1;
        if (lastChannel >= nchanCube) {
           numChannelsLocal = nchanCube - firstChannel;
           lastChannel = nchanCube - 1;
        }
    }

    ASKAPLOG_INFO_STR(logger,"Channel allocation starts at " << firstChannel << " and is " << numChannelsLocal << " in size"
        << " with increment "<< channelInc << " and stops at "<< lastChannel);

    // So the plan is to iterate over each channel ...
    // Calculate the inShapes for each channel and file ....

    // create a stats object for this image. It is a bit messy. The StatsAndMask class requires a shared_ptr for the
    // IImageAccess so we have to create a shared pointer from iacc. Also, we dont want this shared pointer to delete
    // the object (iacc) when it goes out of scope
    boost::shared_ptr<askap::accessors::IImageAccess<>> iaccPtr;
    boost::shared_ptr<askap::utils::StatsAndMask> statsAndMask;
    if ( calcstats ) {
      iaccPtr.reset(&iacc,askap::utility::NullDeleter{});
      statsAndMask.reset(new askap::utils::StatsAndMask{comms,outImgName,iaccPtr});
    }

    ImageBlcTrcMapT imageBlcTrcMap;
    ImageBlcTrcMapT imageBlcTrcLimosShapeMap;
    const bool trimming = parset.getBool("trimming",false);
    if ( trimming ) {
      if (accumulator.weightType() == FROM_BP_MODEL|| accumulator.weightType() == COMBINED) {
        findBoundingBoxes(parset,comms,accumulator,iacc,inImgNames,inWgtNames,nchanCube,
                               firstChannel,lastChannel,channelInc,imageBlcTrcMap,imageBlcTrcLimosShapeMap);
      } else {
        ASKAPCHECK(false,"trimming is only supported for weighttype FromPrimaryBeamModel or Combined");
      }
    }

    for (channel = firstChannel; channel <= lastChannel; channel += channelInc) {

      // clear the lists of input coordinates and shapes
      inCoordSysVec.clear();
      inShapeVec.clear();
      inLinmosCoordSysVec.clear();
      inLinmosShapeVec.clear();

      for (vector<string>::iterator it = inImgNames.begin(); it != inImgNames.end(); ++it) {
        ASKAPLOG_INFO_STR(logger,"Processing Channel " << channel << " of input image " << *it << " which is part of output mosaick " << outImgName);
        if ( !trimming ) {
          const casa::IPosition shape = iacc.shape(*it);

          ASKAPCHECK(shape.nelements()==4,"Work with 4D cubes!");
          ASKAPLOG_INFO_STR(logger," - ImageAccess Shape " << shape);

          casa::IPosition inblc(shape.nelements(),0); // input bottom left corner of this allocation
          casa::IPosition intrc(shape-1);
          nchanCube = shape(3);
          inblc[3] = channel;
          intrc[3] = channel;


          ASKAPCHECK(inblc[3]>=0 && inblc[3]<shape[3], "Start channel is outside the number of channels or negative, shape: "<<shape);
          ASKAPCHECK(trc[3]<=shape[3], "Subcube extends beyond the original cube, shape:"<<shape);

          ASKAPLOG_INFO_STR(logger, " - Corners " << "input bottom lc  = " << inblc << ", input top rc = " << intrc << "\n");
          inCoordSysVec.push_back(iacc.coordSysSlice(*it,inblc,intrc));
          // reset the shape to be the size ...
          intrc = shape;
          intrc[3] = 1;
          const casa::IPosition shape3(intrc);
          ASKAPLOG_INFO_STR(logger, " - Calculated Shape for this accumulator and this image is" << shape3);
          inShapeVec.push_back(shape3);
        } else {
          const auto imageBlcTrcMapIter = imageBlcTrcMap.find(*it);
          ASKAPCHECK(imageBlcTrcMapIter != imageBlcTrcMap.end(),"input image is not in imageBlcTrcMap");
          const auto& blcTrcPair = imageBlcTrcMapIter->second;
          casa::IPosition tempblc = blcTrcPair.first;
          casa::IPosition temptrc = blcTrcPair.second;
          ASKAPLOG_INFO_STR(logger,"trimmed output shape: blc = " << tempblc << "; trc = " << temptrc);
          tempblc[3] = channel;
          temptrc[3] = channel;
          inCoordSysVec.push_back(iacc.coordSysSlice(*it,tempblc,temptrc));
          //casa::IPosition trimmedShape = blcTrcPair.second - blcTrcPair.first + 1;
          casa::IPosition trimmedShape = blcTrcPair.second - blcTrcPair.first + 1;
          trimmedShape[3] = 1;
          inShapeVec.push_back(trimmedShape);    

          if ( useWgtLog ) {
            const auto imageBlcTrcMapIter2 = imageBlcTrcLimosShapeMap.find(*it);
            ASKAPCHECK(imageBlcTrcMapIter2 != imageBlcTrcLimosShapeMap.end(),"input image is not in imageBlcTrcLimosShapeMap");
            const auto& blcTrcPair2 = imageBlcTrcMapIter2->second;
            casa::IPosition tempblc2 = blcTrcPair2.first;
            casa::IPosition temptrc2 = blcTrcPair2.second;
            ASKAPLOG_INFO_STR(logger,"trimmed output shape: blc2 = " << tempblc2 << "; trc2 = " << temptrc2);
            tempblc2[3] = channel;
            temptrc2[3] = channel;
            inLinmosCoordSysVec.push_back(iacc.coordSysSlice(*it,tempblc2,temptrc2));
            //casa::IPosition trimmedShape = blcTrcPair.second - blcTrcPair.first + 1;
            casa::IPosition trimmedShape2 = blcTrcPair2.second - blcTrcPair2.first + 1;
            trimmedShape2[3] = 1;
            inLinmosShapeVec.push_back(trimmedShape2);
          }
        }
      } // got the input shapes for this output image


      // I wonder if we can re-use the accumulator ... lets find out

      casa::IPosition example = inShapeVec[0];
      ASKAPLOG_INFO_STR(logger, "Number of channels in allocations is " << example[3]);

      accumulator.setOutputParameters(inShapeVec, inCoordSysVec);
      ASKAPLOG_INFO_STR(logger, " - Output Shape " << accumulator.outShape());

      casa::IPosition sliceShape = accumulator.outShape();



      // Build the full output cube here:
      // test for master .... (and first channel)
      // this loop uses the accumulator.outShape() method - but only the spatial dimensions
      //


      if (channel == firstChannel) { // build this cube - does not need to loop over the outWgtNames
        if (comms.isMaster()) {
            ASKAPLOG_INFO_STR(logger, "++++++++++++++++++++++++++++++++++++++++++");
            ASKAPLOG_INFO_STR(logger, "Building output mosaic " << outImgName);
            ASKAPLOG_INFO_STR(logger, "++++++++++++++++++++++++++++++++++++++++++");
            casa::IPosition outShape  = accumulator.outShape();
            // has the channel dimension of the allocation - so lets fix that.
            outShape[3] = nchanCube;

            // set one of the input images as a reference for metadata (the first by default)
            const uint psfref = parset.getUint("psfref",0);

            ASKAPLOG_INFO_STR(logger, "Getting brightness info for the output image from input number " << psfref);
            // get pixel units from the selected reference image
            const string units = iacc.getUnits(inImgNames[psfref]);
            ASKAPLOG_INFO_STR(logger, "Got units as " << units);

            ASKAPLOG_INFO_STR(logger, "Getting PSF beam info for the output image from input number " << psfref);
            // get psf beam information from the selected reference image
            Vector<Quantum<double> > psf = iacc.beamInfo(inImgNames[psfref]);
            bool psfValid = (psf.nelements()==3) && (psf[0].getValue("rad")>0) && (psf[1].getValue("rad")>0);
            accessors::BeamList refBeamList = iacc.beamList(inImgNames[psfref]);

            ASKAPLOG_INFO_STR(logger, " Creating output file - Shape " << outShape << " nchanCube " << nchanCube);
            iacc.create(outImgName, outShape, accumulator.outCoordSys());
            if (psfValid) {
                iacc.setBeamInfo(outImgName, psf[0].getValue("rad"), psf[1].getValue("rad"), psf[2].getValue("rad"));
            }
            if (!refBeamList.empty()) {
                iacc.setBeamInfo(outImgName,refBeamList);
            }
            copyKeywords(outImgName, accumulator.getReference(outImgName), keywordsToCopy);
            iacc.setUnits(outImgName,units);
            iacc.addHistory(outImgName, historyLines);
            iacc.makeDefaultMask(outImgName);
            // Save table of mosaic pointing centres & beamsizes?
            saveMosaicTable(outImgName,inImgNames,accumulator.getBeamCentres());

            if (accumulator.outWgtDuplicates()[outImgName]) {
              ASKAPLOG_INFO_STR(logger, "Accumulated weight image " << outWgtName << " already written");
            } else {
              outWgtName = accumulator.outWgtNames()[outImgName];
              ASKAPLOG_INFO_STR(logger, "Writing accumulated weight image to " << outWgtName);
              iacc.create(outWgtName, outShape, accumulator.outCoordSys());
              copyKeywords(outWgtName, accumulator.getReference(outImgName), keywordsToCopy);
              iacc.setUnits(outWgtName,units);
              iacc.addHistory(outWgtName, historyLines);
              iacc.makeDefaultMask(outWgtName);
            }

            if (accumulator.doSensitivity()) {
              outSenName = accumulator.outSenNames()[outImgName];
              ASKAPLOG_INFO_STR(logger, "Writing accumulated sensitivity image to " << outSenName);
              iacc.create(outSenName, outShape, accumulator.outCoordSys());
              copyKeywords(outSenName, accumulator.getReference(outImgName), keywordsToCopy);
              iacc.setUnits(outSenName,units);
              iacc.addHistory(outSenName, historyLines);
              iacc.makeDefaultMask(outSenName);
            }

        }
        // built the output cube for this image, all wait until done (for parallel write)
        if (!serialWrite) {
            comms.barrier();
        }
      }

      // here we have to loop over the channels and let everything else take care of the
      // other planes ...
      //
      // set up the output pixel arrays
      // lets do this one channel at a time


      Array<float> outPix(sliceShape,0.f);
      Array<float> outWgtPix(sliceShape,0.f);
      Array<bool> outMask(sliceShape,casacore::False);
      Array<float> outSenPix;
      //
      if (accumulator.doSensitivity()) {
        outSenPix = Array<float>(sliceShape,0.f);
      }

      // set up an indexing vector for the arrays
      casa::IPosition curpos(outPix.ndim(),0);
      ASKAPASSERT(curpos.nelements()>=2);

      // iterator over planes (e.g. freq & polarisation), regridding and accumulating weights and weighted images

      imagemath::MultiDimArrayPlaneIter planeIter(accumulator.inShape());
      // loop over the input images, reading each in and adding to the output pixel arrays
      // remember this is for the current output mosaick

      for (; planeIter.hasMore(); planeIter.next()) { // this is a loop over the polarisations as well as channels
        for (uInt img = 0; img < inImgNames.size(); ++img ) {
        // set up an iterator for all directionCoordinate planes in the input images

        // short cuts
          string inImgName = inImgNames[img];
          string inWgtName, inSenName, inStokesIName;

          ASKAPLOG_INFO_STR(logger, "Processing input image " << inImgName);
          if (accumulator.weightType() == FROM_WEIGHT_IMAGES || accumulator.weightType() == COMBINED) {
            if (!accumulator.useWtFromHdr()) {
              inWgtName = inWgtNames[img];
              if (accumulator.useWeightsLog()) {
                  ASKAPLOG_INFO_STR(logger, " - and input weightslog " << inWgtName);
              } else {
                  ASKAPLOG_INFO_STR(logger, " - and input weight image " << inWgtName);
              }
            }
          }
          if (accumulator.doSensitivity()) {
            inSenName = inSenNames[img];
            ASKAPLOG_INFO_STR(logger, " - and input sensitivity image " << inSenName);
          }

          if (accumulator.doLeakage() && inStokesINames.size()>img) {
              inStokesIName = inStokesINames[img];
              ASKAPLOG_INFO_STR(logger, " - and input Stokes I image " << inStokesIName);
          }

          casa::IPosition blc; 
          casa::IPosition trc;
          if ( useWgtLog ) {
            getBlcTrc(trimming,iacc,inImgName,imageBlcTrcLimosShapeMap,blc,trc);
          } else {
            getBlcTrc(trimming,iacc,inImgName,imageBlcTrcMap,blc,trc);
          }
          ASKAPLOG_INFO_STR(logger,"trimmed Limos input shape: blc = " << blc << ", trc = " << trc);

          if (nchanCube < 0) {
            nchanCube = shape(3);
          }
          else {
            ASKAPCHECK(nchanCube == shape(3),"Nchan missmatch in merge" );
          }
          // this assumes all allocations
          blc[3] = channel;
          trc[3] = channel;

          if ( ! trimming ) {
            accumulator.setInputParameters(inShapeVec[img], inCoordSysVec[img], img);
          } else {
            if ( useWgtLog ) {
                ASKAPLOG_INFO_STR(logger,"Trimming using weight log file");
              accumulator.setInputParameters(inLinmosShapeVec[img], inLinmosCoordSysVec[img], img);
            } else {
                ASKAPLOG_INFO_STR(logger,"Trimming but not using weight log file");
                accumulator.setInputParameters(inShapeVec[img], inCoordSysVec[img], img);
            }
          }
          Array<float> inPix = iacc.read(inImgName,blc,trc);


          ASKAPLOG_INFO_STR(logger, "Shapes " << shape << " blc " << blc << " trc " << trc << " inpix " << inPix.shape());

          if (parset.getBool("removebeam",false)) {

              Array<float> taylor0;
              Array<float> taylor1;
              Array<float> taylor2;

              ASKAPLOG_INFO_STR(logger, "Scaling Taylor terms -- inImage = " << inImgNames[img]);
              // need to get all the taylor terms for this image
              string ImgName = inImgName;
              int inPixIsTaylor = 0;
              for (int n = 0; n < accumulator.numTaylorTerms(); ++n) {
                  const string taylorN = "taylor." + boost::lexical_cast<string>(n);
                  // find the taylor.0 image for this image
                  size_t pos0 = ImgName.find(taylorN);
                  if (pos0!=string::npos) {
                      ImgName.replace(pos0, taylorN.length(), accumulator.taylorTag());
                      ASKAPLOG_INFO_STR(logger, "This is a Taylor " << n << " image");
                      inPixIsTaylor = n;
                      break;
                  }

                  // now go through each taylor term

              }


              ASKAPLOG_INFO_STR(logger, "To avoid altering images on disk re-reading the Taylor terms");
              for (int n = 0; n < accumulator.numTaylorTerms(); ++n) {

                  size_t pos0 = ImgName.find(accumulator.taylorTag());
                  if (pos0!=string::npos) {
                      const string taylorN = "taylor." + boost::lexical_cast<string>(n);
                      ImgName.replace(pos0, taylorN.length(), taylorN);


                      switch (n)
                      {
                      case 0:
                          ASKAPLOG_INFO_STR(logger, "Reading -- Taylor0");
                          ASKAPLOG_INFO_STR(logger, "Reading -- inImage = " << ImgName);
                          taylor0 = iacc.read(ImgName,blc,trc);
                          ASKAPLOG_INFO_STR(logger, "Shape -- " << taylor0.shape());
                          break;
                      case 1:
                          ASKAPLOG_INFO_STR(logger, "Reading -- Taylor1");
                          ASKAPLOG_INFO_STR(logger, "Reading -- inImage = " << ImgName);
                          taylor1 = iacc.read(ImgName,blc,trc);
                          ASKAPLOG_INFO_STR(logger, "Shape -- " << taylor1.shape());
                          break;
                      case 2:
                          ASKAPLOG_INFO_STR(logger, "Reading -- Taylor2");
                          ASKAPLOG_INFO_STR(logger, "Reading -- inImage = " << ImgName);
                          taylor2 = iacc.read(ImgName,blc,trc);
                          ASKAPLOG_INFO_STR(logger, "Shape -- " << taylor2.shape());
                          break;

                      }
                      ImgName.replace(pos0, accumulator.taylorTag().length(), accumulator.taylorTag());
                  }
                  // now go through each taylor term
              }


              casa::IPosition thispos(taylor0.shape().nelements(),0);
              ASKAPLOG_INFO_STR(logger, " removing Beam for Taylor terms - slice " << thispos);
              accumulator.removeBeamFromTaylorTerms(taylor0,taylor1,taylor2,thispos,iacc.coordSys(inImgName));


              // now we need to set the inPix to be the scaled version
              // Note this means we are reading the Taylor terms 3 times for every
              // read. But I'm not sure this matters.
              // The .copy is not needed, Array assignment doesn't reference, only the (copy)constructor does

              switch (inPixIsTaylor)
              {
              case 0:
                  inPix = taylor0;
                  break;
              case 1:
                  inPix = taylor1;
                  break;
              case 2:
                  inPix = taylor2;
                  break;
              }

          }

          if (parset.getBool("removeleakage",false)) {
              ASKAPCHECK(inPix.shape()[2]==1,"Pol axis should have size 1 for removeleakage");
              // only do this if we're processing a Q, U or V image
              int pol = 0;
              size_t pos = inImgName.find(".q.");
              bool found = false;
              if (pos != string::npos) {
                  found = true;
                  pol = 1;
              }
              if (!found) {
                  pos = inImgName.find(".u.");
                  if (pos != string::npos) {
                      found = true;
                      pol = 2;
                  }
              }
              if (!found) {
                  pos = inImgName.find(".v.");
                  if (pos !=string::npos) {
                      found = true;
                      pol = 3;
                  }
              }
              if (found) {
                  // find corresponding Stokes I image
                  if (inStokesIName=="") {
                      //try to find it
                      inStokesIName = inImgName;
                      inStokesIName.replace(pos,3,".i.");
                  }
                  Array<float> stokesI = iacc.read(inStokesIName,blc,trc);
                  ASKAPCHECK(stokesI.shape()==inPix.shape(),"Stokes I and Pol image shapes don't match");

                  // do leakage correction
                  casa::IPosition thispos(blc);
                  ASKAPLOG_INFO_STR(logger," removing Stokes I leakage using "<<inStokesIName<<" for channel " << thispos(3));
                  accumulator.removeLeakage(inPix,stokesI,pol,thispos,iacc.coordSys(inImgName));
              } else {
                  ASKAPLOG_WARN_STR(logger,"Skipping removeLeakage - cannot determine polarisation of input");
              }
          }

          Array<float> inWgtPix;
          Array<float> inSenPix;

          if (accumulator.weightType() == FROM_WEIGHT_IMAGES || accumulator.weightType() == COMBINED) {

            if (accumulator.useWeightsLog()) {
                ASKAPLOG_INFO_STR(logger,"Reading weights log file :"<< inWgtName);
                accessors::WeightsLog wtlog(inWgtName);
                wtlog.read();

                inWgtPix.resize(inPix.shape());
                inWgtPix = wtlog.weight(channel);

            } else if (accumulator.useWtFromHdr()) {
              ASKAPLOG_INFO_STR(logger,"Reading weights from input image file :"<< inImgName);
              const casacore::Vector<casacore::Float> wts = readWeightsTable(inImgName);
              const size_t size = wts.size();
              if (size > 1) {
                inWgtPix.resize(inPix.shape());
                ASKAPCHECK(channel < size,"Not enough channels in weights table");
                inWgtPix = wts(channel);
              } else if (size == 1) {
                inWgtPix.resize(inPix.shape());
                inWgtPix = wts(0);
              }
              ASKAPCHECK(size>0,"No weights found in image header or extension for image: "<<inImgName);

            } else {
                // use the same blc and trc obtained from the previous getBlcTrc() call for 
                // input image since input/weight/sensitive message (should) have the same size
                inWgtPix = iacc.read(inWgtName,blc,trc);
            }
            ASKAPASSERT(inPix.shape() == inWgtPix.shape());
          }
          if (accumulator.doSensitivity()) {

            // use the same blc and trc obtained from the previous getBlcTrc() call for
            // input image since input/weight/sensitive message (should) have the same size
            inSenPix = iacc.read(inSenName,blc,trc);

            ASKAPASSERT(inPix.shape() == inSenPix.shape());
          }

          // test whether to simply add weighted pixels, or whether a regrid is required
          bool regridRequired = (!accumulator.coordinatesAreEqual()) ;

          // if regridding is required, set up buffer some images
          if ( regridRequired ) {
            ASKAPLOG_INFO_STR(logger, " - regridding -- input pixel grid is different from the output");
            // currently all output planes have full-size, so only initialise once
            // would be faster if this was reduced to the size of the current input image
            if ( accumulator.outputBufferSetupRequired() ) {
              ASKAPLOG_INFO_STR(logger, " - initialising output buffers and the regridder");
              // set up temp images required for regridding
              //accumulator.initialiseOutputBuffers();
              // set up regridder
              accumulator.initialiseRegridder();
            }
            // set up temp images required for regridding
            // need to do this here if some do and some do not have sensitivity images
            // are those of the previous iteration correctly freed?
            accumulator.initialiseOutputBuffers();
            accumulator.initialiseInputBuffers();
          }
          else {
            ASKAPLOG_INFO_STR(logger, " - not regridding -- input pixel grid is the same as the output");
            // not regridding so point output image buffers at the input buffers
            accumulator.initialiseInputBuffers();
            accumulator.redirectOutputBuffers();
          }

          // set the indices of any higher-order dimensions for this slice
          curpos = planeIter.position();
          // i think we have to edit the plane Iter ...
          ASKAPLOG_INFO_STR(logger, " - input slice " << curpos);

          // load input buffer for the current plane
          // Since the plane iterator is set up outside the linmos loop and
          // is shared by all of the input image cubes, when the planes of
          // the input images have different shapes a new temporary plane
          // accessor is needed with a unique shape. To do this, send the
          // current iterator position, rather than the iterator itself.
          accumulator.loadAndWeightInputBuffers(curpos, inPix, inWgtPix, inSenPix);
          
          if ( regridRequired ) {
            // call regrid for any buffered images
            accumulator.regrid();
          }

          // update the accumulation arrays for this plane
          accumulator.accumulatePlane(outPix, outWgtPix, outSenPix, curpos);

        } // over the input images for this
      } // iterated over the polarisation - the accumulator is FULL for this CHANNEL

      //build the mask
      //use the outWgtPix to define the mask

      float itsCutoff = 0.01;

      if (parset.isDefined("cutoff")) itsCutoff = parset.getFloat("cutoff");

      /// This logic is in addition to the mask in the accumulator
      /// which works on an individual beam weight
      /// this is masking the output mosaick - it therefore has slightly
      /// different criteria. In this case the final weight has to be equal to
      /// or bigger than the cutoff.
      /// There is a possible failure mode where due to rounding a pixel maybe
      /// have been masked by the accumulator but missed here.
      /// The first time I implemented this I just used a looser conditional:
      /// > instead of >= - but in the second attempt I decided to replace all
      /// masked pixels by NaN - which has the nice secondary effect of implementing
      /// the FITS mask.

      /// Since I added the scheme to incorporate a variance weight in the mosaicking
      /// I can no longer assume max weight is 1.

      float minVal, maxVal;
      IPosition minPos, maxPos;

      minMax(minVal,maxVal,minPos,maxPos,outWgtPix);
      ASKAPLOG_INFO_STR(logger, "Maximum pixel weight is " << maxVal);
      ASKAPLOG_INFO_STR(logger, "Power fraction cutoff is " << itsCutoff*itsCutoff);

      float wgtCutoff = itsCutoff * itsCutoff * maxVal;

      for(size_t i=0;i<outMask.size();i++){
          if (outWgtPix.data()[i] >= wgtCutoff) {
              outMask.data()[i] = casa::True;
          } else {
              outMask.data()[i] = casa::False;
              setNaN(outWgtPix.data()[i]);
          }
      }

      // deweight the image pixels
      // use another iterator to loop over planes
      ASKAPLOG_INFO_STR(logger, "Deweighting accumulated images");
      imagemath::MultiDimArrayPlaneIter deweightIter(accumulator.outShape());
      for (; deweightIter.hasMore(); deweightIter.next()) {
        curpos = deweightIter.position();
        accumulator.deweightPlane(outPix, outWgtPix, outSenPix, curpos);
      }
      // write accumulated images and weight images
      ASKAPLOG_INFO_STR(logger, "Writing accumulated image to " << outImgName);
      casa::IPosition outShape = accumulator.outShape();

      if (!collective && serialWrite) {
          if (comms.isMaster()) {

            ASKAPLOG_INFO_STR(logger, "Ensuring serial access to cubes");

          }
          else { // this is essentially a serializer - it is required for CASA image types
                 // but not FITS
            int buf;
            int from = comms.rank() - 1;
            comms.receive((void *) &buf,sizeof(int),from);
          }
      }

      casa::IPosition loc(outShape.nelements(),0);
      loc[3] = channel;
      ASKAPLOG_INFO_STR(logger, " - location " << loc);
      iacc.write(outImgName,outPix,outMask,loc);
      // calculate the statistics for the array slice
      if ( calcstats ) {
        statsAndMask->calculate(outImgName,channel,outPix);
      }
      if (accumulator.outWgtDuplicates()[outImgName]) {
          ASKAPLOG_INFO_STR(logger, "Accumulated weight image " << outWgtName << " already written");
      } else {
          iacc.write(outWgtName,outWgtPix,outMask,loc);
      }
      if (accumulator.doSensitivity()) {
          iacc.write(outSenName,outSenPix,outMask,loc);
      }

      if (!collective && serialWrite) {
          if (comms.rank() < comms.nProcs()-1) { // last rank does not use this method
            int buf;
            int to = comms.rank()+1;
            comms.send((void *) &buf,sizeof(int),to);
          }
      }
    }
    // Update stats when all the writing is done
    ASKAPLOG_INFO_STR(logger,"rank " << comms.rank() << " got to barrier");
    comms.barrier();

    if ( calcstats ) {
      // Since the processing of the image channels is distributed among the MPI ranks,
      // the master has to collect all the stats from the worker ranks prior to writing
      // the stats to the image table
      if (comms.isMaster()) {
        statsAndMask->receiveStats();
        statsAndMask->writeStatsToImageTable(outImgName);
        if ( outputStats != "" ) {
            statsAndMask->writeStatsToFile(outputStats);
        }
      } else {
        statsAndMask->sendStats();
      }
    }
  }
}

class linmosMPIApp : public askap::Application
{
    public:

        virtual int run(int argc, char* argv[]) {

            // This class must have scope outside the main try/catch block
            askap::askapparallel::AskapParallel comms(argc, const_cast<const char**>(argv));

            try {
                StatReporter stats;
                LOFAR::ParameterSet subset(config().makeSubset("linmos."));
                SynthesisParamsHelper::setUpImageHandler(subset,comms);
                mergeMPI(subset, comms);
                stats.logSummary();
                return 0;
            }
            catch (const askap::AskapError& e) {
                ASKAPLOG_FATAL_STR(logger, "Askap error in " << argv[0] << ": " << e.what());
                std::cerr << "Askap error in " << argv[0] << ": " << e.what() << std::endl;
                return 1;
            } catch (const std::exception& e) {
                ASKAPLOG_FATAL_STR(logger, "Unexpected exception in " << argv[0] << ": " << e.what());
                std::cerr << "Unexpected exception in " << argv[0] << ": " << e.what()
                    << std::endl;
                return 1;
            }


    };

    private:
        std::string getVersion() const override {
            const std::string pkgVersion = std::string("yandasoft:") + ASKAP_PACKAGE_VERSION;
            return pkgVersion;
        }
};

} // end namespace askap

int main(int argc, char *argv[])
{
    askap::linmosMPIApp app;
    return app.main(argc, argv);
}
