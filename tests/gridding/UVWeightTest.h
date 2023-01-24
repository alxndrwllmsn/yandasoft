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
#include <askap/gridding/GenericUVWeightIndexTranslator.h>
#include <askap/gridding/GenericUVWeightAccessor.h>
#include <askap/gridding/UVWeightIndexTranslationHelper.h>

#include <cppunit/extensions/HelperMacros.h>

#include <casacore/casa/Arrays/IPosition.h>
#include <askap/askap/AskapUtil.h>


namespace askap {

namespace synthesis {

class UVWeightTest : public CppUnit::TestFixture 
{
   CPPUNIT_TEST_SUITE(UVWeightTest);
   CPPUNIT_TEST(testUVWeight);
   CPPUNIT_TEST(testUVWeightFromCube);
   CPPUNIT_TEST(testUVWeightCollection);
   CPPUNIT_TEST_EXCEPTION(testUVWeightCollectionBadIndex, AskapError);
   CPPUNIT_TEST_EXCEPTION(testUVWeightCollectionBadIndexConstAccess, AskapError);
   CPPUNIT_TEST(testIndexTranslation);
   CPPUNIT_TEST(testGenericUVWeightAccessor);
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
   
   
   void testUVWeightCollection() {
       UVWeightCollection collection;
       collection.add(3u, 5u,7u, 1u);

       casacore::Cube<float> wtCube(10,15,1, 1.);
       collection.add(1u, wtCube);

       // check reference semantics
       wtCube(5,7,0) = 0.5;
       // test sparse indices
       for (casa::uInt index = 0; index < 4u; ++index) {
            CPPUNIT_ASSERT(collection.exists(index) == (index % 2 == 1));
       }
       // test shapes (this relies on the non-const interface)
       CPPUNIT_ASSERT(collection.get(1u).shape() == casacore::IPosition(3, 10, 15, 1));
       CPPUNIT_ASSERT(collection.get(3u).shape() == casacore::IPosition(3, 5, 7, 1));
  
       // test values
       for (casacore::uInt u = 0; u < wtCube.nrow(); ++u) {
           for (casacore::uInt v = 0; v < wtCube.ncolumn(); ++v) {
                if ((u < 5u) && (v<7u)) {
                    const float value = collection.get(3u)(u, v, 0u);
                    CPPUNIT_ASSERT_DOUBLES_EQUAL(0., value, 1e-6);
                    CPPUNIT_ASSERT_DOUBLES_EQUAL(0., getWeightViaConstInterface(collection, 3u, u,v, 0u), 1e-6);
                }
                const float expected = (u == 5) && (v == 7) ? 0.5 : 1.;
                const float actual = collection.get(1u)(u, v, 0u);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, actual, 1e-6);
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, getWeightViaConstInterface(collection, 1u, u,v, 0u), 1e-6);
           }
       }
   }

   void testUVWeightCollectionBadIndex() {
       UVWeightCollection collection;
       collection.get(0u);
   }

   void testUVWeightCollectionBadIndexConstAccess() {
       UVWeightCollection collection;
       collection.add(3u, 5u,7u, 1u);
       getWeightViaConstInterface(collection, 1u, 0u, 0u, 0u);
   }

   void testIndexTranslation() {
       GenericUVWeightIndexTranslator ttor(1u, 36u, 36u*100u); 
       boost::shared_ptr<GenericUVWeightIndexTranslator> pTtor(&ttor, utility::NullDeleter());
       UVWeightIndexTranslationHelper<VoidClass> hlp(pTtor);

       for (casacore::uInt beam = 0; beam < 36u; ++beam) {
            for (casacore::uInt field = 0; field < 2u; ++field) {
                 for (casacore::uInt facet = 0; facet < 36u; ++facet) {
                      const casacore::uInt expected = beam + 36u*(field + 100u*facet);
                      CPPUNIT_ASSERT_EQUAL(expected, ttor.indexOf(beam, field, facet));
                      CPPUNIT_ASSERT_EQUAL(expected, hlp.indexOf(beam, field, facet));
                 }
            }
       }

       boost::shared_ptr<GenericUVWeightIndexTranslator> pOtherTtor(new GenericUVWeightIndexTranslator(1u));
       hlp.setTranslator(pOtherTtor);

       for (casacore::uInt beam = 0; beam < 36u; ++beam) {
            for (casacore::uInt field = 0; field < 2u; ++field) {
                 for (casacore::uInt facet = 0; facet < 36u; ++facet) {
                      const casacore::uInt oldExpected = beam + 36u*(field + 100u*facet);
                      CPPUNIT_ASSERT_EQUAL(oldExpected, ttor.indexOf(beam, field, facet));
                      CPPUNIT_ASSERT_EQUAL(beam, pOtherTtor->indexOf(beam, field, facet));
                      CPPUNIT_ASSERT_EQUAL(beam, hlp.indexOf(beam, field, facet));
                 }
            }
       }
   }

   void testGenericUVWeightAccessor() {
       casacore::Cube<float> wtCube(10,15,1, 1.);
       UVWeightCollection collection;
       collection.add(0u, wtCube);

       // check reference semantics
       wtCube(5,7,0) = 0.5;
       const boost::shared_ptr<GenericUVWeightAccessor> genericAcc(new GenericUVWeightAccessor(collection, 0u, 0u, 0u));
       // cast to the interface to exercise polymorphic behaviour
       const boost::shared_ptr<IUVWeightAccessor> acc = genericAcc;
       CPPUNIT_ASSERT(acc);

       // get UVWeight object via the accessor interface for some random beam, field and source/facet (to check translation)
       const UVWeight wt = acc->getWeight(35u, 2u, 1u);
       // check shape
       CPPUNIT_ASSERT_EQUAL(static_cast<casacore::uInt>(wtCube.nrow()), wt.uSize());
       CPPUNIT_ASSERT_EQUAL(static_cast<casacore::uInt>(wtCube.ncolumn()), wt.vSize());
       CPPUNIT_ASSERT_EQUAL(static_cast<casacore::uInt>(wtCube.nplane()), wt.nPlane());
       // test values
       for (casacore::uInt u = 0; u < wtCube.nrow(); ++u) {
           for (casacore::uInt v = 0; v < wtCube.ncolumn(); ++v) {
                const float expected = (u == 5) && (v == 7) ? 0.5 : 1.;
                CPPUNIT_ASSERT_DOUBLES_EQUAL(expected, wt(u, v, 0u), 1e-6);
           }
       }
       
   }

protected:

   float getWeightViaConstInterface(const UVWeightCollection &collection, casacore::uInt index, casacore::uInt u, casacore::uInt v, casacore::uInt plane) {
      const UVWeight wtGrid = collection.get(index);
      return wtGrid(u,v,plane);
   }

   struct VoidClass {
   };

};
    
} // namespace synthesis

} // namespace askap

