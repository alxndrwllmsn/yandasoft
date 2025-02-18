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

#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".measurementequation.imagefftequation");

#include <askap/askap/AskapError.h>
//#include <fft/FFTWrapper.h>

#include <askap/dataaccess/SharedIter.h>
#include <askap/dataaccess/MemBufferDataAccessor.h>
#include <askap/dataaccess/DDCalBufferDataAccessor.h>
#include <askap/scimath/fitting/Params.h>
#include <askap/measurementequation/ImageFFTEquation.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/gridding/BoxVisGridder.h>
#include <askap/gridding/SphFuncVisGridder.h>
#include <askap/scimath/fitting/ImagingNormalEquations.h>
#include <askap/scimath/fitting/DesignMatrix.h>
#include <askap/scimath/fitting/Axes.h>
#include <askap/profile/AskapProfiler.h>

#include <casacore/scimath/Mathematics/RigidVector.h>

#include <casacore/casa/BasicSL/Constants.h>
#include <casacore/casa/BasicSL/Complex.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/ArrayMath.h>

#include <stdexcept>

using askap::scimath::Params;
using askap::scimath::Axes;
using askap::scimath::ImagingNormalEquations;
using askap::scimath::DesignMatrix;
using namespace askap::accessors;

namespace askap
{
  namespace synthesis
  {

    ImageFFTEquation::ImageFFTEquation(const askap::scimath::Params& ip,
        IDataSharedIter& idi) : scimath::Equation(ip),
      askap::scimath::ImagingEquation(ip), itsIdi(idi),
      itsSphFuncPSFGridder(false), itsBoxPSFGridder(false),
      itsUsePreconGridder(false), itsNDir(1)
    {
      itsGridder = IVisGridder::ShPtr(new SphFuncVisGridder());
      init();
    }


    ImageFFTEquation::ImageFFTEquation(IDataSharedIter& idi) :
      itsIdi(idi), itsSphFuncPSFGridder(false), itsBoxPSFGridder(false),
      itsUsePreconGridder(false), itsNDir(1)
    {
      itsGridder = IVisGridder::ShPtr(new SphFuncVisGridder());
      reference(defaultParameters().clone());
      init();
    }

    ImageFFTEquation::ImageFFTEquation(const askap::scimath::Params& ip,
        IDataSharedIter& idi, IVisGridder::ShPtr gridder) :
      scimath::Equation(ip), askap::scimath::ImagingEquation(ip),
      itsGridder(gridder), itsIdi(idi), itsSphFuncPSFGridder(false),
      itsBoxPSFGridder(false), itsUsePreconGridder(false), itsNDir(1)
    {
      init();
    }
    ;

    ImageFFTEquation::ImageFFTEquation(const askap::scimath::Params& ip,
        IDataSharedIter& idi, IVisGridder::ShPtr gridder,
        const LOFAR::ParameterSet& parset) : scimath::Equation(ip),
      askap::scimath::ImagingEquation(ip), itsGridder(gridder), itsIdi(idi),
      itsSphFuncPSFGridder(false), itsBoxPSFGridder(false),
      itsUsePreconGridder(false), itsNDir(1)
    {
      useAlternativePSF(parset);
      init();
    }

    ImageFFTEquation::ImageFFTEquation(IDataSharedIter& idi,
        IVisGridder::ShPtr gridder) :
      itsGridder(gridder), itsIdi(idi), itsSphFuncPSFGridder(false),
      itsBoxPSFGridder(false), itsUsePreconGridder(false), itsNDir(1)
    {
      reference(defaultParameters().clone());
      init();
    }

    ImageFFTEquation::~ImageFFTEquation()
    {
    }

    /// @brief define whether to use an alternative gridder for the PSF
    /// and/or the preconditioner function
    void ImageFFTEquation::useAlternativePSF(const LOFAR::ParameterSet& parset)
    {
      const bool useGentlePCF = parset.getBool("preconditioner.preservecf", parset.isDefined("preconditioner.Names"));
      const bool useBoxPSF = parset.getBool("boxforpsf", false);
      const bool useSphPSF = parset.getBool("sphfuncforpsf", false);

      if (useGentlePCF) {
         itsUsePreconGridder = true;
         ASKAPLOG_INFO_STR(logger,
             "A separate tophat-style gridder will be used to calculate the preconditioner function");
      }

      if (useBoxPSF) {
         itsBoxPSFGridder = true;
         ASKAPLOG_INFO_STR(logger,
             "The box (nearest neighbour) gridder will be used to calculate the PSF");
         if (useSphPSF) {
             ASKAPLOG_WARN_STR(logger,
                 "The spheroidal function will not be used to calculate the PSF");
         }
      } else if (useSphPSF) {
         itsSphFuncPSFGridder = true;
         ASKAPLOG_INFO_STR(logger,
             "The default spheroidal function gridder will be used to calculate the PSF");
      } else {
         ASKAPLOG_INFO_STR(logger,
             "The PSF will be calculated by the same gridder type as used for the model and residuals");
      }
    }

    askap::scimath::Params ImageFFTEquation::defaultParameters()
    {
      Params ip(true);
      ip.add("image");
      return ip;
    }

    ImageFFTEquation::ImageFFTEquation(const ImageFFTEquation& other) :
          Equation(), ImagingEquation()
    {
      operator=(other);
    }

    ImageFFTEquation& ImageFFTEquation::operator=(const ImageFFTEquation& other)
    {
      if(this!=&other)
      {
        static_cast<askap::scimath::Equation*>(this)->operator=(other);
        itsIdi=other.itsIdi;
        itsGridder = other.itsGridder;
        itsSphFuncPSFGridder = other.itsSphFuncPSFGridder;
        itsBoxPSFGridder = other.itsBoxPSFGridder;
        itsUsePreconGridder = other.itsUsePreconGridder;
        itsVisUpdateObject = other.itsVisUpdateObject;
      }
      return *this;
    }

    void ImageFFTEquation::init()
    {
    }

    /// Clone this into a shared pointer
    /// @return shared pointer to a copy
    ImageFFTEquation::ShPtr ImageFFTEquation::clone() const
    {
      return ImageFFTEquation::ShPtr(new ImageFFTEquation(*this));
    }

    /// @brief setup object function to update degridded visibilities
    /// @details For the parallel implementation of the measurement equation we need
    /// inter-rank communication. To avoid introducing cross-dependency of the measurement
    /// equation and the MPI one can use polymorphic object function to sum degridded visibilities
    /// across all required ranks in the distributed case and do nothing otherwise.
    /// By default, this class doesn't alter degridded visibilities.
    /// @param[in] obj new object function (or an empty shared pointer to turn this option off)
    void ImageFFTEquation::setVisUpdateObject(const boost::shared_ptr<IVisCubeUpdate> &obj)
    {
      itsVisUpdateObject = obj;
    }

    /// @brief helper method to verify whether a parameter had been changed
    /// @details This method checks whether a particular parameter is tracked. If
    /// yes, its change monitor is used to verify the status since the last call of
    /// the method, otherwise new tracking begins and true is returned (i.e. to
    /// update all dependent cache).
    /// @param[in] name name of the parameter
    /// @return true if parameter has been updated since the previous call
    bool ImageFFTEquation::notYetDegridded(const std::string &name) const
    {
      std::map<std::string,scimath::ChangeMonitor>::iterator it = itsImageChangeMonitors.find(name);
      bool result = true;
      if (it != itsImageChangeMonitors.end()) {
          result = parameters().isChanged(name,it->second);
          it->second = parameters().monitorChanges(name);
      } else {
          itsImageChangeMonitors[name] = parameters().monitorChanges(name);
      }
      return result;
    }


    void ImageFFTEquation::predict() const
    {
      ASKAPTRACE("ImageFFTEquation::predict");
      const std::vector<std::string> completions(parameters().completions("image"));

      // To minimize the number of data passes, we keep copies of the gridders in memory, and
      // switch between these. This optimization may not be sufficient in the long run.

      itsIdi.chooseOriginal();

      // DDCALTAG -- set increased buffer size
      if (itsNDir > 1) {
         DDCalBufferDataAccessor accBuffer(*itsIdi);
         ASKAPLOG_DEBUG_STR(logger, "calling accBuffer.setNDir("<<itsNDir<<")");
         accBuffer.setNDir(itsNDir);
      }
      int dirIndex = itsNDir > 1 ? -1 : 0;

      ASKAPLOG_DEBUG_STR(logger, "Initialising for model degridding");
      for (std::vector<std::string>::const_iterator it=completions.begin();it!=completions.end();it++)
      {
        string imageName("image"+(*it));
        SynthesisParamsHelper::clipImage(parameters(),imageName);

        if(itsModelGridders.count(imageName)==0) {
          itsModelGridders[imageName]=itsGridder->clone();
        }
	    itsModelGridders[imageName]->customiseForContext(*it);

        if (notYetDegridded(imageName)) {
            ASKAPLOG_DEBUG_STR(logger, "Degridding image "<<imageName);
            const Axes axes(parameters().axes(imageName));
            casacore::Array<imtype> imagePixels;
            #ifdef ASKAP_FLOAT_IMAGE_PARAMS
                imagePixels = parameters().valueF(imageName);
            #else
                imagePixels = parameters().value(imageName);
            #endif
            const casacore::IPosition imageShape(imagePixels.shape());
            itsModelGridders[imageName]->initialiseDegrid(axes, imagePixels);
        }

        // DDCALTAG -- set parameter for increased buffer size
        // could alternatively pass dirIndex * itsIdi->nRow(). See if one is better than the other
        if (itsNDir > 1) {
            if ( imageName.find(".taylor.") != std::string::npos ) {
                // Taylor terms are used, so only increment dirIndex if this is taylor.0
                if ( imageName.find(".taylor.0") != std::string::npos ) {
                    ASKAPLOG_DEBUG_STR(logger, imageName<<" is taylor 0, so incrementing DD-cal buffer index");
                    dirIndex++;
                } else {
                    ASKAPLOG_DEBUG_STR(logger, imageName<<" is not taylor 0 -- not incrementing DD-cal buffer index");
                }
            } else {
                // Taylor terms are not used, so increment dirIndex
                dirIndex++;
            }
            ASKAPLOG_DEBUG_STR(logger, "degridding "<<imageName<<" into buffer "<<dirIndex);
            try {
                // Only possible in DDCalBufferDataAccessor, so cast first
                const boost::shared_ptr<TableVisGridder const> &tGridder =
                       boost::dynamic_pointer_cast<TableVisGridder const>(itsModelGridders[imageName]);
                tGridder->setSourceIndex(dirIndex);
            }
            catch (std::bad_cast&) {}
            ASKAPCHECK(dirIndex <= itsNDir,
                "The number of specified calibration directions is less than the number calibration image fields");
        }

      }
      // Loop through degridding the data
      ASKAPLOG_DEBUG_STR(logger, "Starting to degrid model" );

      // report every 5000000 degridded rows into log in the debug mode
      #ifdef ASKAP_DEBUG
      unsigned long report_every = 5000000;
      unsigned long total_rows = 0;
      unsigned long current_rows = 0;
      #endif // #ifdef ASKAP_DEBUG

      for (itsIdi.init();itsIdi.hasMore();itsIdi.next())
      {
        itsIdi->rwVisibility().set(0.0);
        /*
        ASKAPDEBUGASSERT(itsIdi->nPol() == 4);
        for (casacore::uInt p=1;p<3;++p) {
             itsIdi->rwVisibility().xyPlane(p).set(0.);
        }
        */
        for (std::vector<std::string>::const_iterator it=completions.begin();it!=completions.end();it++)
        {
            string imageName("image"+(*it));
            itsModelGridders[imageName]->degrid(*itsIdi);
        }
        #ifdef ASKAP_DEBUG
        const casacore::uInt nRow = itsIdi->nRow();
        total_rows += nRow;
        current_rows += nRow;
        if (current_rows > report_every) {
            current_rows = 0;
            ASKAPLOG_DEBUG_STR(logger, "Degridded "<<total_rows<<" rows of data");
        }
        #endif // #ifdef ASKAP_DEBUG
      }
      ASKAPLOG_DEBUG_STR(logger, "Finished degridding model" );
    };

    /// @brief assign a different iterator
    /// @details This is a temporary method to assign a different iterator.
    /// All this business is a bit ugly, but should go away when all
    /// measurement equations are converted to work with accessors.
    /// @param idi shared pointer to a new iterator
    void ImageFFTEquation::setIterator(IDataSharedIter& idi)
    {
      itsIdi = idi;
    }


    // Calculate the residual visibility and image. We transform the model on the fly
    // so that we only have to read (and write) the data once. This uses more memory
    // but cuts down on IO
    void ImageFFTEquation::calcImagingEquations(askap::scimath::ImagingNormalEquations& ne) const
    {
      ASKAPTRACE("ImageFFTEquation::calcImagingEquations");


      // We will need to loop over all completions i.e. all sources
      const std::vector<std::string> completions(parameters().completions("image"));

      // To minimize the number of data passes, we keep copies of the gridders in memory, and
      // switch between these. This optimization may not be sufficient in the long run.
      // Set up initial gridders for model and for the residuals. This enables us to
      // do both at the same time.


      for (std::vector<std::string>::const_iterator it=completions.begin();it!=completions.end();it++)
      {
        const string imageName("image"+(*it));
        SynthesisParamsHelper::clipImage(parameters(),imageName);
        if(itsModelGridders.count(imageName)==0) {
           itsModelGridders[imageName]=itsGridder->clone();
        }
        if(itsResidualGridders.count(imageName)==0) {
          itsResidualGridders[imageName]=itsGridder->clone();
        }
        if(itsPSFGridders.count(imageName)==0) {
          if (itsBoxPSFGridder) {
             boost::shared_ptr<BoxVisGridder> psfGridder(new BoxVisGridder);
             itsPSFGridders[imageName] = psfGridder;
          } else if (itsSphFuncPSFGridder) {
             boost::shared_ptr<SphFuncVisGridder> psfGridder(new SphFuncVisGridder);
             itsPSFGridders[imageName] = psfGridder;
          } else {
             itsPSFGridders[imageName] = itsGridder->clone();
          }
        }
        if(itsUsePreconGridder && itsPreconGridders.count(imageName)==0) {
           // preconditioning of higher order terms is set from term 0
           bool isMFS = (imageName.find(".taylor.") != std::string::npos);
           bool isTT0 = (imageName.find(".taylor.0") != std::string::npos);
           if (isTT0 || !isMFS) {
             // Should this be a clone of the psf or the image?
             //itsPreconGridders[imageName] = itsGridder->clone();
             itsPreconGridders[imageName] = itsPSFGridders[imageName]->clone();
           }
        }
        if (itsCoordSystems.count(imageName) == 0) {
          itsCoordSystems[imageName] = SynthesisParamsHelper::coordinateSystem(parameters(),imageName);
        }
      }
      // Now we initialise appropriately
      ASKAPLOG_DEBUG_STR(logger, "Initialising for model degridding and residual gridding" );
      if (completions.size() == 0) {
          ASKAPLOG_WARN_STR(logger,
              "Found no free image parameters, this rank will not contribute usefully to normal equations");
      }
      bool somethingHasToBeDegridded = false;
      for (std::vector<std::string>::const_iterator it=completions.begin();it!=completions.end();it++)
      {
        string imageName("image"+(*it));
        ASKAPLOG_DEBUG_STR(logger, "Initialising for " << imageName);
        const Axes axes(parameters().axes(imageName));
        casacore::Array<imtype> imagePixels;
        #ifdef ASKAP_FLOAT_IMAGE_PARAMS
            imagePixels = parameters().valueF(imageName);
        #else
            imagePixels = parameters().valueT(imageName);
        #endif
        const casacore::IPosition imageShape(imagePixels.shape());
        /// First the model
        itsModelGridders[imageName]->customiseForContext(*it);
        itsModelGridders[imageName]->initialiseDegrid(axes, imagePixels);
        if (!itsModelGridders[imageName]->isModelEmpty()) {
            somethingHasToBeDegridded = true;
        }
        /// Now the residual images, dopsf=false, dopcf=false
        itsResidualGridders[imageName]->customiseForContext(*it);
        itsResidualGridders[imageName]->initialiseGrid(axes, imageShape, false);
        // and PSF gridders, dopsf=true, dopcf=false
        itsPSFGridders[imageName]->customiseForContext(*it);
        itsPSFGridders[imageName]->initialiseGrid(axes, imageShape, true);
        // and PCF gridders, dopsf=false, dopcf=true
        if (itsUsePreconGridder && (itsPreconGridders.count(imageName)>0)) {
            itsPreconGridders[imageName]->customiseForContext(*it);
            itsPreconGridders[imageName]->initialiseGrid(axes, imageShape, false, true);
        }
      }
      // synchronise emtpy flag across multiple ranks if necessary
      if (itsVisUpdateObject) {
          itsVisUpdateObject->aggregateFlag(somethingHasToBeDegridded);
      }
      // Now we loop through all the data
      ASKAPLOG_DEBUG_STR(logger, "Starting degridding model and gridding residuals" );
      size_t counterGrid = 0, counterDegrid = 0;
      for (itsIdi.init();itsIdi.hasMore();itsIdi.next())
      {
        // buffer-accessor, used as a replacement for proper buffers held in the subtable
        // effectively, an array with the same shape as the visibility cube is held by this class
        MemBufferDataAccessor accBuffer(*itsIdi);

        // Accumulate model visibility for all models
        accBuffer.rwVisibility().set(0.0);
        if (somethingHasToBeDegridded) {
            for (std::vector<std::string>::const_iterator it=completions.begin();it!=completions.end();++it) {
                 const std::string imageName("image"+(*it));
                 const std::map<std::string, IVisGridder::ShPtr>::iterator grdIt = itsModelGridders.find(imageName);
                 ASKAPDEBUGASSERT(grdIt != itsModelGridders.end());
                 const IVisGridder::ShPtr degridder = grdIt->second;
                 ASKAPDEBUGASSERT(degridder);
                 if (!degridder->isModelEmpty()) {
                     degridder->degrid(accBuffer);
                     counterDegrid+=accBuffer.nRow();
                 }
            }
            // optional aggregation of visibilities in the case of distributed model
            // somethingHasToBeDegridded is supposed to have consistent value across all participating ranks
            if (itsVisUpdateObject) {
                itsVisUpdateObject->update(accBuffer.rwVisibility());
            }
            //
        }
        accBuffer.rwVisibility() -= itsIdi->visibility();
        accBuffer.rwVisibility() *= float(-1.);
        /// Now we can calculate the residual visibility and image
        size_t tempCounter = 0;
        for (size_t i = 0; i<completions.size(); ++i) {
            const string imageName("image"+completions[i]);
            if (parameters().isFree(imageName)) {
                itsResidualGridders[imageName]->grid(accBuffer);
                itsPSFGridders[imageName]->grid(accBuffer);
                if (itsUsePreconGridder && (itsPreconGridders.count(imageName)>0)) {
                    itsPreconGridders[imageName]->grid(accBuffer);
                }
                tempCounter += accBuffer.nRow();
            }
        }
        counterGrid += tempCounter;
      }
      ASKAPLOG_DEBUG_STR(logger, "Finished degridding model and gridding residuals" );
      ASKAPLOG_DEBUG_STR(logger, "Number of accessor rows iterated through is "<<counterGrid<<" (gridding) and "<<
                        counterDegrid<<" (degridding)");

      // We have looped over all the data, so now we have to complete the
      // transforms and fill in the normal equations with the results from the
      // residual gridders
      if (itsUsePreconGridder) {
        ASKAPLOG_DEBUG_STR(logger,
            "Adding residual image, PSF, preconditioner function and weights image to the normal equations" );
      } else {
        ASKAPLOG_DEBUG_STR(logger,
            "Adding residual image, PSF and weights image to the normal equations" );
      }
      for (std::vector<std::string>::const_iterator it=completions.begin();it!=completions.end();++it)
      {
        const string imageName("image"+(*it));
        const casacore::IPosition imageShape(parameters().shape(imageName));
        ASKAPLOG_INFO_STR(logger,"Name: " << imageName << " Shape: " << imageShape);
        casacore::Array<imtype> imageDeriv(imageShape);
        casacore::Array<imtype> imagePSF(imageShape);
        casacore::Array<imtype> imageWeight(imageShape);
        itsResidualGridders[imageName]->finaliseGrid(imageDeriv);

        // for debugging/research, store grid prior to FFT
        // boost::shared_ptr<TableVisGridder> tvg = boost::dynamic_pointer_cast<TableVisGridder>(itsResidualGridders[imageName]);
        // if (tvg) {
        //    tvg->storeGrid("uvcoverage"+(*it),0);
        //}
        // end debugging code

        itsPSFGridders[imageName]->finaliseGrid(imagePSF);
        itsResidualGridders[imageName]->finaliseWeights(imageWeight);

        // just to deallocate the grid memory
        itsModelGridders[imageName]->finaliseDegrid();

        /*{
          casacore::Array<imtype> imagePSFWeight(imageShape);
          itsPSFGridders[imageName]->finaliseWeights(imagePSFWeight);
          const double maxPSFWeight = casacore::max(imagePSFWeight);
          if (maxPSFWeight > 0) {
              const double psfScalingFactor = casacore::max(imageWeight)/maxPSFWeight;
              //std::cout<<"psf peak = "<<casacore::max(imagePSF)<<" maxPSFWeight = "<<maxPSFWeight<<" factor="<<
              //     psfScalingFactor<<" maxImageWeight="<<casacore::max(imageWeight)<<std::endl;
              imagePSF *= psfScalingFactor;
              // now psf has the same peak as the weight image
          } else {
             if (maxPSFWeight < 0) {
                 ASKAPTHROW(AskapError, "PSF weight is supposed to be non-negative, you have "<<maxPSFWeight);
             }
             // do nothing for zero weight, it just means that this part of normal equations is empty. However,
             // we may still be able to have data after summing all parts of the NE
          }
        }*/

        casacore::IPosition vecShape(1, imagePSF.nelements());

        casacore::Vector<imtype> imagePreconVec;
        if (itsUsePreconGridder && (itsPreconGridders.count(imageName)>0)) {
          casacore::Array<imtype> imagePrecon(imageShape);
          itsPreconGridders[imageName]->finaliseGrid(imagePrecon);
          imagePreconVec.reference(imagePrecon.reform(vecShape));
        }

        {
          casacore::IPosition reference(4, imageShape(0)/2, imageShape(1)/2, 0, 0);
          casacore::Vector<imtype> imagePSFVec(imagePSF.reform(vecShape));
          casacore::Vector<imtype> imageWeightVec(imageWeight.reform(vecShape));
          casacore::Vector<imtype> imageDerivVec(imageDeriv.reform(vecShape));
          ne.addSlice(imageName, imagePSFVec, imageWeightVec, imagePreconVec,
              imageDerivVec, imageShape, reference,itsCoordSystems[imageName]);
        }
      }
    }

  }

}
