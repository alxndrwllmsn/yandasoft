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

#include <askap/measurementequation/ComponentEquation.h>
#include <askap/scimath/fitting/LinearSolver.h>
#include <askap/dataaccess/DataIteratorStub.h>
#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/BasicSL/Constants.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/casa/Quanta/MVPosition.h>

#include <cppunit/extensions/HelperMacros.h>

#include <askap/askap/AskapError.h>

#include <cmath>

using std::abs;

#include <boost/shared_ptr.hpp>

using namespace askap;
using namespace askap::scimath;

namespace askap
{
  namespace synthesis
  {

    class ComponentEquationSpectralTest : public CppUnit::TestFixture
    {

      CPPUNIT_TEST_SUITE(ComponentEquationSpectralTest);
      CPPUNIT_TEST(testConstructNormalEquations);
      CPPUNIT_TEST(testSolveNormalEquationsLSQR);
      CPPUNIT_TEST_EXCEPTION(testNoFree, AskapError);
      CPPUNIT_TEST_SUITE_END();

      private:
        ComponentEquation *p1, *p2;
        Params *params1, *params2;
        accessors::IDataSharedIter idi;

      public:
        void setUp()
        {
          idi= accessors::IDataSharedIter(new accessors::DataIteratorStub(1));

          params1 = new Params;
          params1->add("flux.i.cena", 100.0);
          params1->add("direction.ra.cena", 0.5);
          params1->add("direction.dec.cena", -0.3);
          params1->add("shape.bmaj.cena", 30.0*casacore::C::arcsec);
          params1->add("shape.bmin.cena", 20.0*casacore::C::arcsec);
          params1->add("shape.bpa.cena", -55*casacore::C::degree);
          params1->add("flux.spectral_index.cena", -0.7);
          params1->add("flux.ref_freq.cena", 1e9);
          params1->fix("flux.ref_freq.cena");

          p1 = new ComponentEquation(*params1, idi);

          params2 = new Params;
          params2->add("flux.i.cena", 100.0);
          params2->add("direction.ra.cena", 0.500005);
          params2->add("direction.dec.cena", -0.300003);
          params2->add("shape.bmaj.cena", 33.0*casacore::C::arcsec);
          params2->add("shape.bmin.cena", 22.0*casacore::C::arcsec);
          params2->add("shape.bpa.cena", -57*casacore::C::degree);
          params2->add("flux.spectral_index.cena", -0.71);
          params2->add("flux.ref_freq.cena", 1e9);
          params2->fix("flux.ref_freq.cena");

          p2 = new ComponentEquation(*params2, idi);

        }

        void tearDown()
        {
          delete p1;
          delete p2;
        }

        void testCopy()
        {
          Params ip;
          ip.add("Value0");
          ip.add("Value1");
          ip.add("Value2");
          delete p1;
          p1 = new ComponentEquation(ip, idi);
          delete p2;
          p2 = new ComponentEquation(*p1);
          CPPUNIT_ASSERT(p2->parameters().names().size()==3);
          CPPUNIT_ASSERT(p2->parameters().names()[0]=="Value0");
          CPPUNIT_ASSERT(p2->parameters().names()[1]=="Value1");
          CPPUNIT_ASSERT(p2->parameters().names()[2]=="Value2");
        }
        
        void testConstructNormalEquations()
        {
          GenericNormalEquations ne; //(*params1);
          p2->calcEquations(ne);
          std::vector<std::string> names(params1->freeNames());
          for (uint row=0;row<names.size();row++) {
            for (uint col=0;col<names.size();col++) {
                if (ne.normalMatrix(names[row], names[col]).nelements() != 0) {
                    const casacore::Matrix<double> nm = ne.normalMatrix(names[row], names[col]);
                    casacore::IPosition ip(nm.shape());
                    CPPUNIT_ASSERT(ip(0)==1);
                    CPPUNIT_ASSERT(ip(1)==1);
                }
            }
          }
        }

        void testSolveNormalEquationsLSQR()
        {
            testSolveNormalEquations("LSQR");
        }

        void testNoFree()
        {
          GenericNormalEquations ne; //(*params1);
          p1->predict();
          //p2->calcEquations(ne);
          Quality q;
          params2->fix("flux.i.cena");
          params2->fix("direction.ra.cena");
          params2->fix("direction.dec.cena");
          params2->fix("shape.bmaj.cena");
          params2->fix("shape.bmin.cena");
          params2->fix("shape.bpa.cena");
          params2->fix("flux.spectral_index.cena");
          params2->fix("flux.ref_freq.cena");
          LinearSolver solver1(LinearSolver::KeepAllSingularValues);
          solver1.addNormalEquations(ne);
          // the following line is not supposed to throw exceptions any more
          solver1.solveNormalEquations(*params2,q);
        }

      private:
        void testSolveNormalEquations(const std::string& solverType)
        {
// Predict with the "perfect" parameters"
          p1->predict();
// Calculate gradients using "imperfect" parameters"
          GenericNormalEquations ne; //(*params2);
          p2->calcEquations(ne);
          Quality q;
          LinearSolver solver1(LinearSolver::KeepAllSingularValues);
          solver1.addNormalEquations(ne);
          solver1.setAlgorithm(solverType);
          solver1.solveNormalEquations(*params2,q);
          //CPPUNIT_ASSERT(abs(q.cond()/4.99482e+12-1.0)<0.001);
          // check that it converged to a correct result
          const Params &solution = *params2;
          CPPUNIT_ASSERT(abs(solution.scalarValue("direction.dec.cena")+0.3)<1e-5);
          CPPUNIT_ASSERT(abs(solution.scalarValue("direction.ra.cena")-0.5)<1e-5);
          CPPUNIT_ASSERT(abs(solution.scalarValue("flux.i.cena")-100)<3e-1);
          CPPUNIT_ASSERT(abs(solution.scalarValue("flux.spectral_index.cena")+0.7)<1e-3);
          CPPUNIT_ASSERT(abs(solution.scalarValue("shape.bmaj.cena")-0.00014859)<1e-6);
          CPPUNIT_ASSERT(abs(solution.scalarValue("shape.bmin.cena")-9.66813e-5)<1e-6);
          CPPUNIT_ASSERT(abs(solution.scalarValue("shape.bpa.cena")+0.952876)<1e-2);
        }

    };

  }
}
