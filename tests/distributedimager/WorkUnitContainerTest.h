/// @file
///
/// Unit test for WorkUnitContainer class (used in distributed imager)
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

#include <askap/distributedimager/WorkUnitContainer.h>

#include <cppunit/extensions/HelperMacros.h>

namespace askap {

namespace synthesis {

class WorkUnitContainerTest : public CppUnit::TestFixture 
{
   CPPUNIT_TEST_SUITE(WorkUnitContainerTest);
   CPPUNIT_TEST(testAdd);
   CPPUNIT_TEST(testChannelMerge);
   CPPUNIT_TEST_SUITE_END();
public:
   
   void testAdd() {
      WorkUnitContainer wuc;
      cp::ContinuumWorkUnit wu;
      wu.set_dataset("test1.ms");
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(0u), wuc.size());
      CPPUNIT_ASSERT(wuc.begin() == wuc.end());
      wuc.add(wu);
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1u), wuc.size());
      CPPUNIT_ASSERT(wuc.begin() != wuc.end());
      addWorkUnits(wuc, wu, 5u, 1);
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(6u), wuc.size());
      CPPUNIT_ASSERT(wuc.begin() != wuc.end());
   }

   void testChannelMerge() {
      WorkUnitContainer wuc;
      cp::ContinuumWorkUnit wu;
      wu.set_localChannel(1u);
      wu.set_dataset("test1.ms");
      wu.set_beam(0u);
      CPPUNIT_ASSERT_EQUAL(1u, wu.get_nchan());
      addWorkUnits(wuc, wu, 5u, 1);
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(5u), wuc.size());
      wuc.compressWorkUnits();
      CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1u), wuc.size());
      CPPUNIT_ASSERT(wuc.begin() != wuc.end());
      const cp::ContinuumWorkUnit newWu = *wuc.begin();
      CPPUNIT_ASSERT_EQUAL(wu.get_dataset(), newWu.get_dataset());
      CPPUNIT_ASSERT_EQUAL(wu.get_beam(), newWu.get_beam());
      CPPUNIT_ASSERT_EQUAL(5u, newWu.get_nchan());
   }
protected:
 
   /// @brief helper method to add work units for a range of local channels
   /// @details 
   /// @param[in] wuc container to work with
   /// @param[in] wuTemplate  template work unit for the units to be added (only local channel 
   ///                        is altered by adding increments)
   /// @param[in] nchan number of channels 
   /// @param[in] increment local channel increment for the added work units
   static void addWorkUnits(WorkUnitContainer& wuc, const cp::ContinuumWorkUnit& wuTemplate, unsigned int nchan, int increment) {
      for (unsigned int chan = 0; chan < nchan; ++chan) {
           const int adjustment = increment * chan;
           const unsigned int startChan = wuTemplate.get_localChannel();
           if (adjustment < 0) {
               CPPUNIT_ASSERT(startChan >= -adjustment);
           }
           cp::ContinuumWorkUnit newWu(wuTemplate);
           newWu.set_localChannel(startChan + adjustment);
           wuc.add(newWu);
      }
   }
   

};
    
} // namespace synthesis

} // namespace askap

