/// @file
///
/// This file contains utilities to search for a support of the convolution function.
/// They could, in principle, be moved to a higher level (to Base), but left here
/// for now as they are not logically a part of fitting, but would introduce
/// a casacore dependence if moved to askap.
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
/// @author Max Voronkov <maxim.voronkov@csiro.au>

#ifndef SUPPORT_SEARCHER_TCC
#define SUPPORT_SEARCHER_TCC

#include <askap/askap/AskapError.h>
#include <askap/profile/AskapProfiler.h>
#include <complex>
//#include <askap/askap/AskapLogging.h>
//ASKAP_LOGGER(gsslogger, ".gridding.supportsearcher");

namespace askap {

namespace synthesis {

// std::norm (and therefore casacore::norm) in gcc g++ library appear to
// be implemented optimised for polar encoding of complex number
// (abs(x)**2) ! This is faster in theory and practice on my
// workstation
// Overload to avoid unnecessary squaring of real data
template<typename Type>
Type ampFunction(const std::complex<Type>& c)
{
  return real(c)*real(c) + imag(c)*imag(c);
}

// could use template<typename Type>, but want to exclude types like casacore::Complex
inline double ampFunction(const double& c)
{
  return fabs(c);
}
inline float ampFunction(const float& c)
{
  return fabs(c);
}

// Include const reference to matrix c to ensure that the correct amplitude function is used
template<typename Type>
double processAmplitudeThreshold(const casacore::Matrix<std::complex<Type> >& c, const double& cutoff)
{
  return fabs(cutoff)*fabs(cutoff);
}

inline double processAmplitudeThreshold(const casacore::Matrix<double>& c, const double& cutoff)
{
  return fabs(cutoff);
}
inline float processAmplitudeThreshold(const casacore::Matrix<float>& c, const double& cutoff)
{
  return fabs(cutoff);
}

/// @brief determine the peak and its position
/// @details This method fills only itsPeakPos and itsPeakVal. It is
/// normally called from one of the search methods, but could be called
/// separately.
/// @param[in] in input 2D matrix with an image
template<typename T>
void SupportSearcher::findPeak(const casacore::Matrix<T> &in)
{
  ASKAPDEBUGTRACE("SupportSearcher::findPeak");

  itsPeakPos.resize(in.shape().nelements(),casacore::False);
  itsPeakPos = 0;
  itsPeakVal = -1;
  #ifdef _OPENMP_WORKING
  #pragma omp parallel default(shared)
  {
  #pragma omp for
  #endif
  for (int iy=0;iy<int(in.ncolumn());++iy) {
       double tempPeakNorm = -1;
       int tempPeakX = 0, tempPeakY = 0;
       for (int ix=0;ix<int(in.nrow());++ix) {
       const double curNorm = ampFunction(casacore::DComplex(in(ix,iy)));
            if(tempPeakNorm< curNorm) {
               tempPeakX = ix;
               tempPeakY = iy;
               tempPeakNorm = curNorm;
            }
       }
       #ifdef _OPENMP_WORKING
       #pragma omp critical
       {
       #endif
       if (itsPeakVal < 0 || (itsPeakVal*itsPeakVal < tempPeakNorm)) {
           itsPeakPos(0) = tempPeakX;
           itsPeakPos(1) = tempPeakY;
           itsPeakVal = sqrt(tempPeakNorm);
       }
       #ifdef _OPENMP_WORKING
       }
       #endif
  }
  #ifdef _OPENMP_WORKING
  }
  #endif
#ifdef ASKAP_DEBUG
  if (itsPeakVal<0) {
      ASKAPTHROW(CheckError, "An empty matrix has been passed to SupportSearcher::findPeak, please investigate. Shape="<<
                 in.shape());
  }
  if (std::isinf(itsPeakVal) || std::isnan(itsPeakVal)) {
      debugStoreImage(in);
  }
  ASKAPCHECK(!std::isnan(itsPeakVal), "Peak value is not a number, please investigate. itsPeakPos="<<itsPeakPos);
  ASKAPCHECK(!std::isinf(itsPeakVal), "Peak value is infinite, please investigate. itsPeakPos="<<itsPeakPos);
#endif // #ifdef ASKAP_DEBUG

  ASKAPCHECK(itsPeakVal>0.0, "Unable to find peak in the support searcher (in either making a convolution function or fitting the PSF), all values appear to be zero. itsPeakVal="
             << itsPeakVal);
}

/// @brief debug method to save the matrix
/// @details This method is used for debugging only and stores given
/// complex matrix (i.e. if NaN appears in the CF calculation). It does
/// nothing for the generic value type.
/// @param[in] in input 2D matrix with an image
template<>
void SupportSearcher::debugStoreImage(const casacore::Matrix<casacore::Complex> &in);


/// @brief full search which determines the peak
/// @details This search method doesn't assume anything about the peak and
/// searches for its position and peak beforehand. The search starts at the
/// edges and progresses towards the peak. The edge of the support region
/// is where the value first time exceeds the cutoff*peakVal, or cutoff*value
/// if value given as a second parameter is positive
/// @param[in] in input 2D matrix with an image
/// @param[in] value optional peak value, if a positive value is given it will be used
/// instead of the peak amplitude (although the positon of the peak will still be searched for)
template<typename T>
void SupportSearcher::search(const casacore::Matrix<T> &in, const double value)
{
  findPeak(in);
  if (value > 0) {
      itsPeakVal = value;
  }
  doSupportSearch(in);
}


/// @brief do actual support search
/// @details This method assumes that peak has already been found and
/// implements the actual search of blc and trc of the support region.
/// @param[in] in input 2D matrix with an image
template<typename T>
void SupportSearcher::doSupportSearch(const casacore::Matrix<T> &in)
{
  ASKAPDEBUGTRACE("SupportSearcher::doSupportSearch");

  ASKAPDEBUGASSERT(in.shape().nelements() == 2);
  ASKAPDEBUGASSERT(itsPeakPos.nelements() == 2);
  itsBLC.resize(2,casacore::False);
  itsTRC.resize(2,casacore::False);
  itsBLC = -1;
  itsTRC = -1;
  ASKAPCHECK(itsPeakVal>0.0,
             "A positive peak value of the convolution function is expected inside doSupportSearch, itsPeakVal=" <<
             itsPeakVal);
  ASKAPCHECK(!std::isinf(itsPeakVal), "Peak value is infinite, this shouldn't happen. itsPeakPos="<<itsPeakPos);
  ASKAPCHECK(itsPeakPos[0]>0 && itsPeakPos[1]>0, "Peak position of the convolution function "<<itsPeakPos<<
             " is too close to the edge, increase maxsupport");
  ASKAPCHECK(itsPeakPos[0] + 1 < int(in.nrow()) && itsPeakPos[1] + 1 < int(in.ncolumn()),
             "Peak position of the convolution function "<<itsPeakPos<<" is too close to the edge, increase maxsupport");

  const double absCutoff = processAmplitudeThreshold(in, itsCutoff*itsPeakVal);

  #ifdef _OPENMP_WORKING
  #pragma omp parallel sections
  {
  #pragma omp section
  {
  #endif
  for (int ix = 0; ix<=itsPeakPos(0); ++ix) {
       if (ampFunction(in(ix, itsPeakPos(1))) > absCutoff) {
           itsBLC(0) = ix;
           break;
       }
  }

  #ifdef _OPENMP_WORKING
  }
  #pragma omp section
  {
  #endif

  for (int iy = 0; iy<=itsPeakPos(1); ++iy) {
       if (ampFunction(in(itsPeakPos(0),iy)) > absCutoff) {
           itsBLC(1) = iy;
           break;
       }
  }

  #ifdef _OPENMP_WORKING
  }
  #pragma omp section
  {
  #endif

  for (int ix = int(in.nrow())-1; ix>=itsPeakPos(0); --ix) {
       if (ampFunction(in(ix, itsPeakPos(1))) > absCutoff) {
           itsTRC(0) = ix;
           break;
       }
  }

  #ifdef _OPENMP_WORKING
  }
  #pragma omp section
  {
  #endif

  for (int iy = int(in.ncolumn())-1; iy>=itsPeakPos(1); --iy) {
       if (ampFunction(in(itsPeakPos(0),iy)) > absCutoff) {
           itsTRC(1) = iy;
           break;
       }
  }

  #ifdef _OPENMP_WORKING
  }
  }
  #endif

  ASKAPCHECK((itsBLC(0)>=0) && (itsBLC(1)>=0) && (itsTRC(0)>=0) &&
             (itsTRC(1)>=0), "Unable to find the support on one of the coordinates (try decreasing the value of .gridder.cutoff) Effective support is 0. itsBLC="<<itsBLC<<" itsTRC="<<itsTRC<<" itsPeakPos="<<itsPeakPos<<" in.shape()="<<in.shape()<<" absCutoff="<<absCutoff<<" itsPeakVal="<<itsPeakVal);
}

/// @brief extend support area to include diagonals
/// @details This search method assumes the peak has been found
/// and the support area has been set.
/// It then tries to extend this area by determining the pixels furthest
/// in x or y from the peak that are above the cutoff and connected to the peak
/// @param[in] in input 2D matrix with an image
template<typename T>
void SupportSearcher::extendedSupport(const casacore::Matrix<T> &in)
{
  const double absCutoff = processAmplitudeThreshold(in, itsCutoff*itsPeakVal);

  // zig-zag to furthest point above cutoff, first in x, then find extent in y, at extreme find extent in x, etc
  // assumes a roughly elliptical beam shape at the cutoff level
  int x0 = itsPeakPos(0);
  int y0 = itsPeakPos(1);
  int x=itsTRC(0);
  int y=y0;
  int lastx = 0;
  int lasty = 0;
  int xmax = in.nrow() - 1;
  int ymax = in.ncolumn() - 1;
  bool pos = true;
  while (lastx != x || lasty != y){
      lastx = x;
      lasty = y;
      // find highest y still above cutoff
      while (y < ymax && ampFunction(in(x,y+1)) > absCutoff) y++;
      int yposstep = y - y0;
      y = lasty;
      // find lowest y above cutoff
      while (y > 0 && ampFunction(in(x,y-1)) > absCutoff) y--;
      int ynegstep = y0 - y;
      pos = (yposstep >= ynegstep);
      // pick largest step
      if (pos) {
          y = y0 + yposstep;
      } else {
          y = y0 - ynegstep;
      }
      // now find highest x above cutoff at new y value
      while (x < xmax && ampFunction(in(x+1,y)) > absCutoff) x++;
  }
  //ASKAPLOG_INFO_STR(gsslogger, "Original support area: "<<itsBLC<<", "<<itsTRC);
  // Found furthest point
  itsTRC(0) = x;
  if (pos) itsTRC(1) = casacore::max(y,itsTRC(1));
  if (!pos) itsBLC(1) = casacore::min(y,itsBLC(1));
  // Do we need to do the same on other side? For now just assume symmetry
  itsBLC(0) = x0 - (itsTRC(0) - x0);
  if (pos) itsBLC(1) = y0 - (itsTRC(1) - y0);
  if (!pos) itsTRC(1) = y0 + (y0 - itsBLC(1));
  //ASKAPLOG_INFO_STR(gsslogger, "Extended support area: "<<itsBLC<<", "<<itsTRC);
}

/// @brief search assuming the peak is in the centre
/// @details This search method assumes the peak is in the centre of the
/// image and has a given value. The search starts at the edges and
/// terminated as soon as the absolute value higher than
/// cutoff*value has been found. Giving value of 1. effectively means
/// that the cutoff is an absolute cutoff (default).
/// @param[in] in input 2D matrix with an image
/// @param[in] value assumed peak value
template<typename T>
void SupportSearcher::searchCentered(const casacore::Matrix<T> &in, double value)
{
  itsPeakVal = value;
  itsPeakPos.resize(in.shape().nelements(), casacore::False);
  ASKAPDEBUGASSERT(itsPeakPos.nelements() == 2);
  itsPeakPos = in.shape();
  itsPeakPos(0)/=2;
  itsPeakPos(1)/=2;
  doSupportSearch(in);
}


} // namespace synthesis

} // namespace askap


#endif // #ifndef SUPPORT_SEARCHER_TCC
