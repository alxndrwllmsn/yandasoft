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


ASKAP_LOGGER(logger, ".linmos");

using namespace askap::synthesis;

namespace askap {

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


    for (channel = firstChannel; channel <= lastChannel; channel += channelInc) {

      // clear the lists of input coordinates and shapes
      inCoordSysVec.clear();
      inShapeVec.clear();
      for (vector<string>::iterator it = inImgNames.begin(); it != inImgNames.end(); ++it) {

        ASKAPLOG_INFO_STR(logger,"Processing Channel " << channel << " of input image " << *it << " which is part of output mosaick " << outImgName);


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

          const casa::IPosition shape = iacc.shape(inImgName);
          casa::IPosition blc(shape.nelements(),0);
          casa::IPosition trc(shape-1);

          if (nchanCube < 0) {
            nchanCube = shape(3);
          }
          else {
            ASKAPCHECK(nchanCube == shape(3),"Nchan missmatch in merge" );
          }
          // this assumes all allocations
          blc[3] = channel;
          trc[3] = channel;

          accumulator.setInputParameters(inShapeVec[img], inCoordSysVec[img], img);
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
                // use weights images
                const casa::IPosition shape = iacc.shape(inWgtName);
                casa::IPosition blc(shape.nelements(),0);
                casa::IPosition trc(shape-1);

                blc[3] = channel;
                trc[3] = channel;

                inWgtPix = iacc.read(inWgtName,blc,trc);
            }
            ASKAPASSERT(inPix.shape() == inWgtPix.shape());
          }
          if (accumulator.doSensitivity()) {

            const casa::IPosition shape = iacc.shape(inSenName);
            casa::IPosition blc(shape.nelements(),0);
            casa::IPosition trc(shape-1);

            blc[3] = channel;
            trc[3] = channel;

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
