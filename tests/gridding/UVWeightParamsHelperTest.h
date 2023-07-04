/// @file
///
/// Unit test for UVWeightParamsHelper class 
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

#include <askap/gridding/UVWeightParamsHelper.h>
#include <askap/scimath/fitting/Params.h>
#include <cppunit/extensions/HelperMacros.h>

namespace askap {

namespace synthesis {

class UVWeightParamsHelperTest : public CppUnit::TestFixture 
{
   CPPUNIT_TEST_SUITE(UVWeightParamsHelperTest);
   CPPUNIT_TEST(testAdd);
   CPPUNIT_TEST_SUITE_END();
public:
   void setUp() {
      itsParams.reset(new scimath::Params());
   }
  
   void tearDown() {
      itsParams.reset();
   }

   void testAdd() {
      UVWeightParamsHelper hlp(itsParams);
      CPPUNIT_ASSERT(!hlp.exists("gc"));
      CPPUNIT_ASSERT(itsParams);

      UVWeightCollection collection;
      collection.add(3u, 5u,7u, 1u);

      casacore::Cube<float> wtCube(10,15,1, 1.);
      collection.add(1u, wtCube);

      boost::shared_ptr<GenericUVWeightIndexTranslator> ttor(new GenericUVWeightIndexTranslator(1u,0u,0u));

      hlp.addUVWeights("gc", collection, ttor);

      CPPUNIT_ASSERT(hlp.exists("gc"));
      CPPUNIT_ASSERT(itsParams->has("uvweight.gc.3"));
      CPPUNIT_ASSERT(itsParams->has("uvweight.gc.1"));
      CPPUNIT_ASSERT(itsParams->has("uvweight_indices.gc"));
   }
private:
   boost::shared_ptr<scimath::Params> itsParams;
};
    
} // namespace synthesis

} // namespace askap

