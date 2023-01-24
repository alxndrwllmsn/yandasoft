/// @file
///
/// Unit test for UVWeight and related classes (used in traditional weighting)
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

#include <askap/gridding/UVWeightCollection.h>
#include <cppunit/extensions/HelperMacros.h>

#include <casacore/casa/Arrays/IPosition.h>


namespace askap {

namespace synthesis {

class UVWeightTest : public CppUnit::TestFixture 
{
   CPPUNIT_TEST_SUITE(UVWeightTest);
   CPPUNIT_TEST(testUVWeight);
   CPPUNIT_TEST(testUVWeightFromCube);
   CPPUNIT_TEST_SUITE_END();
public:
   
   void testUVWeight() {
      UVWeight wtGrid1(10u,15u,1u);
      CPPUNIT_ASSERT_EQUAL(10u, wtGrid1.uSize());
      CPPUNIT_ASSERT_EQUAL(15u, wtGrid1.vSize());
      CPPUNIT_ASSERT_EQUAL(1u, wtGrid1.nPlane());
      CPPUNIT_ASSERT(!wtGrid1.empty());

      // reference semantics in copy constructor
      UVWeight wtGrid2(wtGrid1);
      CPPUNIT_ASSERT_EQUAL(10u, wtGrid2.uSize());
      CPPUNIT_ASSERT_EQUAL(15u, wtGrid2.vSize());
      CPPUNIT_ASSERT_EQUAL(1u, wtGrid2.nPlane());
      CPPUNIT_ASSERT(!wtGrid2.empty());

      UVWeight wtGrid3;
      CPPUNIT_ASSERT(wtGrid3.empty());
      // reference semantics in assignment operator
      wtGrid3 = wtGrid1;
      CPPUNIT_ASSERT_EQUAL(10u, wtGrid3.uSize());
      CPPUNIT_ASSERT_EQUAL(15u, wtGrid3.vSize());
      CPPUNIT_ASSERT_EQUAL(1u, wtGrid3.nPlane());
      CPPUNIT_ASSERT(!wtGrid3.empty());
      
      // the following assignment should propagate to both weight grids (reference semantics)
      wtGrid1(5,7,0) = 1.5;
      for (casacore::uInt u = 0; u < wtGrid1.uSize(); ++u) {
           for (casacore::uInt v = 0; v < wtGrid1.vSize(); ++v) {
                const float expected = ( (u == 5) && (v == 7) ) ? 1.5 : 0.;
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, wtGrid1(u,v,0u), 1e-6);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, wtGrid2(u,v,0u), 1e-6);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, wtGrid3(u,v,0u), 1e-6);
           }
      }
   }

   void testUVWeightFromCube() {
      casacore::Cube<float> wtCube(10,15,2, 1.);
      UVWeight wtGrid(wtCube);
      CPPUNIT_ASSERT_EQUAL(10u, wtGrid.uSize());
      CPPUNIT_ASSERT_EQUAL(15u, wtGrid.vSize());
      CPPUNIT_ASSERT_EQUAL(2u, wtGrid.nPlane());
      CPPUNIT_ASSERT(!wtGrid.empty());
 
      // check reference semantics
      wtCube(5,7,0) = 0.5;
      for (casacore::uInt u = 0; u < wtGrid.uSize(); ++u) {
           for (casacore::uInt v = 0; v < wtGrid.vSize(); ++v) {
                for (casacore::uInt plane = 0; plane < wtGrid.nPlane(); ++plane) {
                     const float expected = ( (u == 5) && (v == 7)  && (plane == 0) ) ? 0.5 : 1.;
                     CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, wtGrid(u,v,plane), 1e-6);
                }
           }
      }
   }
};
    
} // namespace synthesis

} // namespace askap

