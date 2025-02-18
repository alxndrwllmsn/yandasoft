/// @file
///
/// IVisWeights: Abstract base class for
///              Interface definition for visibility weight calculators
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
/// @author Urvashi Rau <rurvashi@aoc.nrao.edu>
///
#ifndef ASKAP_SYNTHESIS_IVISWEIGHTS_H_
#define ASKAP_SYNTHESIS_IVISWEIGHTS_H_

#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/Cube.h>

#include <askap/scimath/fitting/Axes.h>

#include <boost/shared_ptr.hpp>
#include <askap/dataaccess/SharedIter.h>

#include <boost/shared_ptr.hpp>

namespace askap
{
	namespace synthesis
	{

		/// @brief Abstract Base Class for all visibility weight calculators.
		/// Visibilities can be weighted in various ways while implenenting
		/// matched filtering methods that search for the presence or absence
		/// of a particular pattern.
		/// 
		/// @ingroup gridding
		class IVisWeights
		{
			public:

			/// Shared pointer definition
			typedef boost::shared_ptr<IVisWeights> ShPtr;

                        /// @brief empty destructor to keep the compiler happy
			virtual ~IVisWeights();
			
			/// Clone a copy of this Gridder
			virtual ShPtr clone() const = 0;

			/// @brief Set the context
			/// @param order The index of the enumerated expansion
			virtual void setParameters(int order)=0;
			
                        /// @brief Calculate the visibility weight.
                        /// @param freq Frequency
                        /// @note we used to pass sample index and polarisation here which has not been used so far, so I (MV) took it out
                        virtual float getWeight(double freq) const = 0;

		};
	}
}
#endif                                            /*VISWEIGHTS_H_*/
