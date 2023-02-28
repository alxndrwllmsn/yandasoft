/// @file
///
/// Unit test for various procedures computing uv-weight or doing other operations with it
///
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

#include <askap/gridding/ConjugatesAdderFFT.h>
#include <cppunit/extensions/HelperMacros.h>


namespace askap {

namespace synthesis {

class UVWeightCalculatorTest : public CppUnit::TestFixture 
{
   CPPUNIT_TEST_SUITE(UVWeightCalculatorTest);
   CPPUNIT_TEST(testConjugatesAdderFFT);
   CPPUNIT_TEST_SUITE_END();
public:

   void testConjugatesAdderFFT() {
        casacore::Int pos[4][2] = {{200, 233}, {242, 324}, {270, 300}, {300, 222}};
        casacore::Matrix<float> buf(512, 512, 0.f);
        for (casacore::uInt pt = 0; pt < 4; ++pt) {
             buf(pos[pt][0], pos[pt][1]) = 2.f * pt / 6;
        }
 
        ConjugatesAdderFFT calc;
        calc.process(buf);
        
        for (casacore::uInt pt = 0; pt < 4; ++pt) {
             // factor of 2 still remains because we added conjugates as opposed to just taking the real part, so this way of doing it preserves the values
             CPPUNIT_ASSERT_DOUBLES_EQUAL(2.f*pt/6, buf(pos[pt][0], pos[pt][1]), 1e-6);
             // check conjugates as well (note, we don't have 0th row or column which do not have conjugates with even-sized grids)
             CPPUNIT_ASSERT_DOUBLES_EQUAL(2.f*pt/6, buf(512 - pos[pt][0], 512 - pos[pt][1]), 1e-6);
        }
   }

};
    
} // namespace synthesis

} // namespace askap

