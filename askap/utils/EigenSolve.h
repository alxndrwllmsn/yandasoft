/// @file
/// 
/// @brief interface to eigen problem solver
/// @details
///
/// An interface module for low-level routines to solve eigenproblem or
/// a singular value decomposition problem
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

#ifndef EIGENSOLVE_H
#define EIGENSOLVE_H

#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Exceptions/Error.h>
#include <casacore/casa/BasicSL/Complex.h>

/// @brief interface to eigen problem solver
/// @details
///
/// An interface module for low-level routines to solve eigenproblem or
/// a singular value decomposition problem
class EigenSolver {
   casa::Matrix<casa::Complex> vec; // eigenvectors
   casa::Vector<casa::Double> val;  // eigenvalues
   casa::Matrix<casa::Complex> vecV; // SVD's vector V. (not V^t)
public:
   // solve eigenproblem and fill vec and val data members
   void solveEigen(const casa::Matrix<casa::Complex> &in);
   // return results
   const casa::Matrix<casa::Complex>& getEigenVectors() const;
   const casa::Vector<casa::Double>& getEigenValues() const;
   const casa::Matrix<casa::Complex>& getV() const;
};

#endif // EIGENSOLVE_H
