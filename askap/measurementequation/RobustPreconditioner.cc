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

#include <askap/measurementequation/RobustPreconditioner.h>

#include <askap/askap_synthesis.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".measurementequation.wienerpreconditioner");

#include <askap/askap/AskapError.h>
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/scimath/utils/PaddingUtils.h>
#include <askap/profile/AskapProfiler.h>

#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/MatrixMath.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/lattices/Lattices/SubLattice.h>
#include <casacore/lattices/Lattices/ArrayLattice.h>
#include <casacore/lattices/LatticeMath/LatticeFFT.h>
#include <casacore/lattices/LEL/LatticeExpr.h>
using namespace casa;

#include <iostream>
#include <cmath>
using std::abs;

namespace askap
{
  namespace synthesis
  {

    RobustPreconditioner::RobustPreconditioner() :
	    itsRobust(0.0)
    {
    }
    
    RobustPreconditioner::RobustPreconditioner(const float& robust) :
	    itsRobust(robust)
    {
    }
        
    IImagePreconditioner::ShPtr RobustPreconditioner::clone()
    {
	    return IImagePreconditioner::ShPtr(new RobustPreconditioner(*this));
    }
    
    bool RobustPreconditioner::doPreconditioning(casacore::Array<float>& psf,
                                                 casacore::Array<float>& dirty,
                                                 casacore::Array<float>& pcf) const
    {
      ASKAPTRACE("RobustPreconditioner::doPreconditioning");
      ASKAPLOG_INFO_STR(logger, "Applying Robust filter with robustness parameter " << itsRobust);
      
      const float maxPSFBefore=casacore::max(psf);
      ASKAPLOG_INFO_STR(logger, "Peak of PSF before Robust filtering = " << maxPSFBefore);
      casacore::ArrayLattice<float> lpsf(psf);
      casacore::ArrayLattice<float> ldirty(dirty);
      
      // we need to pad to twice the size in the image plane in order to avoid wraparound
      casacore::IPosition shape = lpsf.shape();
      casacore::ArrayLattice<casacore::Complex> scratch(shape);
      scratch.copyData(casacore::LatticeExpr<casacore::Complex>(toComplex(lpsf)));
      
      LatticeFFT::cfft2d(scratch, True);

      // Construct a Robust filter
      
      casacore::ArrayLattice<casacore::Complex> robustfilter(shape);
      // Normalize relative to the average weight
      const double noisepower(pow(10.0, 2*itsRobust));
      const double rnp(1.0/(noisepower*maxPSFBefore));
      robustfilter.copyData(casacore::LatticeExpr<casacore::Complex>(1.0/(sqrt(real(scratch*conj(scratch)))*rnp+1.0)));
            
      // Apply the filter to the lpsf
      // (reuse the ft(lpsf) currently held in 'scratch')
      scratch.copyData(casacore::LatticeExpr<casacore::Complex> (robustfilter * scratch));
      
      /*
	SynthesisParamsHelper::saveAsCasaImage("dbg.img",casacore::amplitude(scratch.asArray()));       
	//SynthesisParamsHelper::saveAsCasaImage("dbg.img",lpsf.asArray());
	throw AskapError("This is a debug exception");
      */
      
      LatticeFFT::cfft2d(scratch, False);       
      lpsf.copyData(casacore::LatticeExpr<float>(real(scratch)));
      const float maxPSFAfter = casacore::max(psf);
      ASKAPLOG_INFO_STR(logger, "Peak of PSF after Robust filtering  = " << maxPSFAfter);
      psf *= maxPSFBefore/maxPSFAfter;
 
      ASKAPLOG_INFO_STR(logger, "Normalized to unit peak");
     
      // Apply the filter to the dirty image
      scratch.copyData(casacore::LatticeExpr<casacore::Complex>(toComplex(ldirty)));       
      
      LatticeFFT::cfft2d(scratch, True);
      
      scratch.copyData(casacore::LatticeExpr<casacore::Complex> (robustfilter * scratch));
      LatticeFFT::cfft2d(scratch, False);

      ldirty.copyData(casacore::LatticeExpr<float>(real(scratch)));      
      dirty *= maxPSFBefore/maxPSFAfter;
      
      return true;
    }
  }
}


