/// @file
/// @brief Simple disk illumination model
/// @details This class represents a simple illumination model,
/// which is just a disk of a certain radius with a hole in the centre.
/// Optionally a phase slope can be applied to simulate offset pointing.
///
/// @copyright (c) 2008 CSIRO
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
/// @author Max Voronkov <maxim.voronkov@csiro.au>

#include <casacore/casa/BasicSL/Constants.h>
#include <casacore/casa/Arrays/ArrayMath.h>


#include <askap/gridding/DiskIllumination.h>
#include <askap/askap/AskapError.h>
// temporary
#include <askap/measurementequation/SynthesisParamsHelper.h>
#include <askap/askap/AskapLogging.h>
ASKAP_LOGGER(logger, ".gridding.diskillumination");
//

#include <askap/profile/AskapProfiler.h>

using namespace askap;
using namespace askap::synthesis;

/// @brief construct the model
/// @param[in] diam disk diameter in metres
/// @param[in] blockage a diameter of the central hole in metres
DiskIllumination::DiskIllumination(double diam, double blockage) :
   itsDiameter(diam), itsBlockage(blockage)
{
  ASKAPDEBUGASSERT(diam>0);
  ASKAPDEBUGASSERT(blockage>=0);
  ASKAPDEBUGASSERT(diam > blockage);
}

/// @brief obtain illumination pattern
/// @details This is the main method which populates the
/// supplied uv-pattern with the values corresponding to the model
/// represented by this object. It has to be overridden in the
/// derived classes. An optional phase slope can be applied to
/// simulate offset pointing.
/// @param[in] freq frequency in Hz for which an illumination pattern is required
/// @param[in] pattern a UVPattern object to fill
/// @param[in] l angular offset in the u-direction (in radians)
/// @param[in] m angular offset in the v-direction (in radians)
/// @param[in] pa parallactic angle, or strictly speaking the angle between
/// uv-coordinate system and the system where the pattern is defined (unused)
void DiskIllumination::getPattern(double freq, UVPattern &pattern, double l,
                          double m, double pa) const
{
    // zero value of the pattern by default
    pattern.pattern().set(0.);
}

/// @param[in] freq frequency in Hz for which an illumination pattern is required
/// @param[in] pattern a UVPattern object to fill
/// @param[in] imageCentre ra & dec of the image centre
/// @param[in] beamCentre ra & dec of the beam pointing centre
/// @param[in] pa polarisation position angle (in radians)
/// @param[in] isPSF bool indicting if this gridder is for a PSF
/// @param[in] feed  feed number for case where pattern differs between feeds
void DiskIllumination::getPattern(double freq, UVPattern &pattern,
                          const casacore::MVDirection &imageCentre,
                          const casacore::MVDirection &beamCentre,
                          const double pa, const bool isPSF, const int feed) const
{
    ASKAPTRACE("DiskIllumination::getPattern");
    const casacore::uInt oversample = pattern.overSample();
    const double cellU = pattern.uCellSize()/oversample;
    const double cellV = pattern.vCellSize()/oversample;

    // calculate beam offset
    double l = 0.0;
    double m = 0.0;
    if (!isPSF) {
        l = sin(beamCentre.getLong() - imageCentre.getLong()) * cos(beamCentre.getLat());
        m = sin(beamCentre.getLat()) * cos(imageCentre.getLat()) -
            cos(beamCentre.getLat()) * sin(imageCentre.getLat()) *
            cos(beamCentre.getLong() - imageCentre.getLong());
    }

    // scaled l and m to take the calculations out of the loop
    // these quantities are effectively dimensionless
    const double lScaled = 2.*casacore::C::pi*cellU *l;
    const double mScaled = 2.*casacore::C::pi*cellV *m;

    // zero value of the pattern by default
    pattern.pattern().set(0.);

    ASKAPCHECK(std::abs(std::abs(cellU/cellV)-1.)<1e-7,
               "Rectangular cells are not supported at the moment, you have ("<<cellU<<" , "<<cellV<<")");

    const double cell = std::abs(cellU*(casacore::C::c/freq));

    const double dishRadiusInCells = itsDiameter/(2.0*cell);

    // squares of the disk and blockage area radii
    const double rMaxSquared = casacore::square(dishRadiusInCells);
    double rMinSquared = casacore::square(itsBlockage/(2.0*cell));

    // sizes of the grid to fill with pattern values
    const casacore::uInt nU = pattern.uSize();
    const casacore::uInt nV = pattern.vSize();

    ASKAPCHECK(rMaxSquared>rMinSquared, "Disk hole is supposed to be less than the disk size, you have diameter="<<
               itsDiameter<<" (rMaxSquared="<<rMaxSquared<<") blockage="<<itsBlockage<<" (rMinSquared="<<
               rMinSquared<<")");

    ASKAPCHECK((casacore::square(double(nU)) > rMaxSquared) &&
               (casacore::square(double(nV)) > rMaxSquared),
               "The pattern buffer passed to DiskIllumination::getPattern is too small for the given model. "
               "Sizes should be greater than "<<sqrt(rMaxSquared)<<" on each axis, you have "
                <<nU<<" x "<<nV);

    // maximum possible support for this class corresponds to the dish size
    pattern.setMaxSupport(1+2*casacore::uInt(dishRadiusInCells)/oversample);

    if (rMinSquared + 1 >= rMaxSquared) {
        // we don't have enough resolution in the uv-plane to represent the blockage
        // (image selected is probably too small to get the sidelobes affected by the blockage)
        // we need at least 1 pixel wide annulus.
        --rMinSquared;
    }

    double sum=0.; // normalisation factor
    #ifdef _OPENMP_WORKING
    #pragma omp parallel default(shared)
    {
        #pragma omp for reduction(+:sum)
    #endif
        for (casacore::uInt iU=0; iU<nU; ++iU) {
             const double offsetU = double(iU)-double(nU)/2.;
             const double offsetUSquared = casacore::square(offsetU);
             for (casacore::uInt iV=0; iV<nV; ++iV) {
                  const double offsetV = double(iV)-double(nV)/2.;
                  const double offsetVSquared = casacore::square(offsetV);
                  const double radiusSquared = offsetUSquared + offsetVSquared;
                  if ( (radiusSquared >= rMinSquared) && (radiusSquared <= rMaxSquared)) {
                       // don't need to multiply by wavelength here because we
                       // divided the radius (i.e. the illumination pattern is given
                       // in a relative coordinates in frequency
                       const double phase = lScaled*offsetU + mScaled*offsetV;
                       pattern(iU, iV) = imtypeComplex(cos(phase), -sin(phase));
                       sum += 1.;
                  }
             }
	    }
    #ifdef _OPENMP_WORKING
    }
    #endif

    ASKAPCHECK(sum > 0., "Integral of the disk should be non-zero");
    pattern.pattern() *= imtypeComplex(float(nU)*float(nV)/float(sum),0.);
}

/// @brief check whether the pattern is symmetric
/// @details Some illumination patterns like this one are trivial and known a priori to
/// be symmetric. This method always returns true to reflect this
/// @return always true
bool DiskIllumination::isSymmetric() const
{
  return true;
}

/// @brief check whether the output pattern is image-based, rather than an illumination pattern.
/// @details Some illumination patterns need to be generated in the image domain, and given
/// the standard usage (FFT to image-domain for combination with other functions) any image
/// domain function may as well stay in the image domain. So check the state before doing the FFT.
/// @return false
bool DiskIllumination::isImageBased() const
{
  return false;
}

/// @brief check whether the output pattern is feed dependent
/// @details Some illumination patterns vary with feed (number) and no shortcuts can
/// be taken
/// @return false
bool DiskIllumination::isFeedDependent() const
{
    return false;
}
