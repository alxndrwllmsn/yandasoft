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

#include <askap/gridding/BoxVisGridder.h>
#include <askap/gridding/SphFuncVisGridder.h>
#include <askap/gridding/AWProjectVisGridder.h>
#include <askap/gridding/AProjectWStackVisGridder.h>
#include <askap/gridding/WStackVisGridder.h>
#include <askap/gridding/AWProjectVisGridder.h>
#include <askap/gridding/WProjectVisGridder.h>
#include <askap/gridding/ATCAIllumination.h>
#include <askap/gridding/DiskIllumination.h>
#include <askap/scimath/fitting/Params.h>
#include <askap/measurementequation/ComponentEquation.h>
#include <askap/dataaccess/DataIteratorStub.h>
#include <casacore/casa/aips.h>
#include <casacore/casa/Arrays/Matrix.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/casa/Quanta/MVPosition.h>
#include <casacore/casa/BasicSL/Constants.h>
#include <askap/askap/AskapError.h>
#include <askap/gridding/VisGridderFactory.h>
#include <askap/gridding/IUVWeightAccessor.h>
#include <askap/scimath/utils/ImageUtils.h>
#include <askap/scimath/fft/FFT2DWrapper.h>
#include <askap/gridding/SupportSearcher.h>
#include <askap/gridding/GenericUVWeightBuilder.h>
#include <askap/gridding/UVWeightGridder.h>
#include <askap/gridding/ConjugatesAdderFFT.h>

#include <cppunit/extensions/HelperMacros.h>

#include <stdexcept>
#include <boost/shared_ptr.hpp>

using namespace askap::scimath;

namespace askap
{
  namespace synthesis
  {

    class TableVisGridderTest : public CppUnit::TestFixture
    {

      CPPUNIT_TEST_SUITE(TableVisGridderTest);
      CPPUNIT_TEST(testForwardBox);
      CPPUNIT_TEST(testReverseBox);
      CPPUNIT_TEST(testCreateAbstract);
      CPPUNIT_TEST(testUnknownGridder);
      CPPUNIT_TEST(testForwardSph);
      CPPUNIT_TEST(testReverseSph);
      CPPUNIT_TEST(testUVWeightApplication);
      CPPUNIT_TEST(testBuildingUVWeight);
      CPPUNIT_TEST(testUVWeightBuilderInit);
      CPPUNIT_TEST(testUVWeightGridder);
      CPPUNIT_TEST(testUVWeightGridderBeamAndFieldSelection);
      CPPUNIT_TEST(testForwardAWProject);
      CPPUNIT_TEST(testReverseAWProject);
      CPPUNIT_TEST(testForwardWProject);
      CPPUNIT_TEST(testReverseWProject);
      CPPUNIT_TEST(testForwardWStack);
      CPPUNIT_TEST(testReverseWStack);
      CPPUNIT_TEST(testForwardAWProject);
      CPPUNIT_TEST(testReverseAWProject);
      CPPUNIT_TEST(testForwardAProjectWStack);
      CPPUNIT_TEST(testReverseAProjectWStack);
      CPPUNIT_TEST(testForwardATCAIllumination);
      CPPUNIT_TEST(testReverseATCAIllumination);
      CPPUNIT_TEST_SUITE_END();

  private:
      boost::shared_ptr<BoxVisGridder> itsBox;
      boost::shared_ptr<SphFuncVisGridder> itsSphFunc;
      boost::shared_ptr<AWProjectVisGridder> itsAWProject;
      boost::shared_ptr<WProjectVisGridder> itsWProject;
      boost::shared_ptr<WStackVisGridder> itsWStack;
      boost::shared_ptr<AProjectWStackVisGridder> itsAProjectWStack;

      accessors::IDataSharedIter idi;
      boost::shared_ptr<Axes> itsAxes;
      boost::shared_ptr<casa::Array<imtype> > itsModel;
      boost::shared_ptr<casa::Array<imtype> > itsModelPSF;
      boost::shared_ptr<casa::Array<imtype> > itsModelWeights;

      // for a test of uv-weight application
      struct UVWeightAccessorTester : virtual public IUVWeightAccessor {
           explicit UVWeightAccessorTester(const UVWeight &wt) : itsWt(wt) {}
           virtual UVWeight getWeight(casacore::uInt beam, casacore::uInt field, casacore::uInt source) const {
               itsBeams.insert(beam);
               itsFields.insert(field);
               itsSourceIndices.insert(source);
               return itsWt;
           }
           // fields to keep track the parameters getWeight is called with
           mutable std::set<casacore::uInt> itsBeams;
           mutable std::set<casacore::uInt> itsFields;
           mutable std::set<casacore::uInt> itsSourceIndices;
        private:
           UVWeight itsWt;
      };

      // for a test of gridder-based weight building
      struct UVWeightBuilderTester : virtual public IUVWeightBuilder {
           // mocked up weight access method
           virtual UVWeight getWeight(casacore::uInt, casacore::uInt, casacore::uInt) const {
               ASKAPTHROW(AskapError, "This method is not supposed to be called in the test we have");
           }

           // mocked up initialise callback method
           virtual void initialise(casacore::uInt uSize, casacore::uInt vSize, casacore::uInt nPlanes) {
               itsWt.resize(uSize, vSize, nPlanes);
               itsWt.set(0.f);
           }

           // weight addition method, ignores all indices returning (by reference) the same object
           virtual UVWeight addWeight(casacore::uInt beam, casacore::uInt field, casacore::uInt source) {
               itsBeams.insert(beam);
               itsFields.insert(field);
               itsSourceIndices.insert(source);
               
               return itsWt;
           }

           // mocked up merge method
           virtual void merge(const IUVWeightBuilder &) {
               ASKAPTHROW(AskapError, "This method is not supposed to be called in the current tests");
           }

           // mocked up finialise method, calls calculator for the first plane of the cube
           virtual UVWeightCollection& finalise(const IUVWeightCalculator &calc) {
               CPPUNIT_ASSERT(itsWt.nplane() > 0u);
               casacore::Matrix<float> slice = itsWt.xyPlane(0);
               calc.process(slice);
               itsWtCollection.reset(new UVWeightCollection());
               itsWtCollection->add(0, itsWt);
               return *itsWtCollection;
           }

           // fields to keep track the parameters addWeight is called with
           mutable std::set<casacore::uInt> itsBeams;
           mutable std::set<casacore::uInt> itsFields;
           mutable std::set<casacore::uInt> itsSourceIndices;
 
           // helper method to get the shape of the cube (to be able to test it earlier before addWeight is called)
           casacore::IPosition wtShape() const { return itsWt.shape(); }

        private:
           casacore::Cube<float> itsWt;
           boost::shared_ptr<UVWeightCollection> itsWtCollection;
      };

      // mocked up weight calculator class, the only thing it does is to check the shape
      // and write a weight of 100 to the central pixel (to be able to check that this method
      // has been called
      struct TestUVWeightCalculator : virtual public IUVWeightCalculator {
          TestUVWeightCalculator(const casacore::IPosition& expectedShape) : itsExpectedShape(expectedShape) {}

          // mocked up weight calculation method
          virtual void process(casacore::Matrix<float> &wt) const  {
             CPPUNIT_ASSERT_EQUAL(itsExpectedShape, wt.shape());
             CPPUNIT_ASSERT(itsExpectedShape.nelements() > 0u);
             wt(wt.nrow() / 2, wt.ncolumn() / 2) = 100.;
          }
        private:
          /// expected shape, note that process is called with one plane only, so it's 2D shape
          casacore::IPosition itsExpectedShape;
      };

  public:
      void setUp()
      {
        idi = accessors::IDataSharedIter(new accessors::DataIteratorStub(1));

        Params ip;
        ip.add("flux.i.cena", 100.0);
        ip.add("direction.ra.cena", 0.5*casacore::C::degree);
        ip.add("direction.dec.cena", -0.3*casacore::C::degree);
        ip.add("shape.bmaj.cena", 0.0);
        ip.add("shape.bmin.cena", 0.0);
        ip.add("shape.bpa.cena", 0.0);

        ComponentEquation ce(ip, idi);
        ce.predict();

	itsBox.reset(new BoxVisGridder());
        itsSphFunc.reset(new SphFuncVisGridder());
        boost::shared_ptr<IBasicIllumination> illum(new DiskIllumination(120.0, 10.0));
        // illumination models can easily be shared between gridders if parameters,
        // like dish size and blockage are the same
        itsAWProject.reset(new AWProjectVisGridder(illum, 10000.0, 9, 1e-3, 1, 128, 1));
        itsAProjectWStack.reset(new AProjectWStackVisGridder(illum, 10000.0, 9, 0., 1, 128, 1));
        itsWProject.reset(new WProjectVisGridder(10000.0, 9, 1e-3, 1, 128, 0, ""));
        itsWStack.reset(new WStackVisGridder(10000.0, 9));

        double cellSize=10*casa::C::arcsec;

        casa::Matrix<double> xform(2,2,0.);
        xform.diagonal().set(1.);

        itsAxes.reset(new Axes());
        itsAxes->addDirectionAxis(casa::DirectionCoordinate(casa::MDirection::J2000,
                     casa::Projection(casa::Projection::SIN), 0.,0.,cellSize,cellSize,xform,256.,256.));

        itsModel.reset(new casa::Array<imtype>(casa::IPosition(4,512,512,1,1)));
        itsModel->set(0.0);
        itsModelPSF.reset(new casa::Array<imtype>(casa::IPosition(4,512,512,1,1)));
        itsModelPSF->set(0.0);
        itsModelWeights.reset(new casa::Array<imtype>(casa::IPosition(4,512,512,1,1)));
        itsModelWeights->set(0.0);
      }

      void tearDown()
      {
      }

      void testReverseBox()
      {
         itsBox->initialiseGrid(*itsAxes, itsModel->shape(), false);
         itsBox->grid(*idi);
         itsBox->finaliseGrid(*itsModel);
         itsBox->finaliseWeights(*itsModelWeights);
         itsBox->initialiseGrid(*itsAxes, itsModel->shape(), true);
         itsBox->grid(*idi);
         itsBox->finaliseGrid(*itsModelPSF);
      }
      void testForwardBox()
      {
         itsBox->initialiseDegrid(*itsAxes, *itsModel);
         itsBox->degrid(*idi);
      }

      void testCreateAbstract()
      {
         // TableVisGridder is still an abstract class
         // calling createGridder static method should raise an
         // exception
         LOFAR::ParameterSet parset;
         CPPUNIT_ASSERT_THROW(TableVisGridder::createGridder(parset),AskapError);
      }

      void testUnknownGridder()
      {
         // here we attempt to create a gridder, which is not recognized by a factory
         // The code should try to load a properly named dynamic library and because
         // it is not there an exception will be thrown
         LOFAR::ParameterSet parset;
         parset.add("gridder","AGridderWithSuchNameShouldDefinitelyBeUnknownToTheSystsem");
         CPPUNIT_ASSERT_THROW(VisGridderFactory::make(parset),AskapError);
      }

      void testReverseSph()
      {
        itsSphFunc->initialiseGrid(*itsAxes, itsModel->shape(), true);
        itsSphFunc->grid(*idi);
        itsSphFunc->finaliseGrid(*itsModelPSF);
        itsSphFunc->initialiseGrid(*itsAxes, itsModel->shape(), false);
        itsSphFunc->grid(*idi);
        itsSphFunc->finaliseGrid(*itsModel);
        itsSphFunc->finaliseWeights(*itsModelWeights);
      }
      void testForwardSph()
      {
        itsSphFunc->initialiseDegrid(*itsAxes, *itsModel);
        itsSphFunc->degrid(*idi);
      }

      void testUVWeightApplication() {
        // setup weighting with some square window
        UVWeight wt(512,512,1);
        for (casacore::uInt iu = 240; iu <= 272; ++iu) {
             for (casacore::uInt iv = 240; iv <= 272; ++iv) {
                  wt(iu,iv,0) = 1.;
             }
        }
        boost::shared_ptr<UVWeightAccessorTester> wtAcc(new UVWeightAccessorTester(wt));
        itsSphFunc->setUVWeightAccessor(wtAcc);

        casacore::Array<imtype> newImg(itsModel->shape(), static_cast<imtype>(0.));
        itsSphFunc->initialiseGrid(*itsAxes, newImg.shape(), false);
        itsSphFunc->grid(*idi);
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1u), wtAcc->itsBeams.size());
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1u), wtAcc->itsFields.size());
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1u), wtAcc->itsSourceIndices.size());
        CPPUNIT_ASSERT_EQUAL(0u, *wtAcc->itsBeams.begin());
        CPPUNIT_ASSERT_EQUAL(0u, *wtAcc->itsFields.begin());
        CPPUNIT_ASSERT_EQUAL(0u, *wtAcc->itsSourceIndices.begin());
        itsSphFunc->finaliseGrid(newImg);
        //scimath::saveAsCasaImage("tst.img", newImg);
        // assess the result now
        casacore::Matrix<imtypeComplex> buf(newImg.shape().getFirst(2), imtypeComplex(static_cast<imtype>(0.)));
        casacore::setReal(buf, newImg.nonDegenerate());
        FFT2DWrapper<imtypeComplex> fftWrapper;
        fftWrapper(buf, false);
        // use the support searcher to get rough bounds of the unweighted part
        // simulated source is 100 Jy, the result should be windowed gridded visibility largely concentrated in
        // the inner 32 pixels. Leave the cutoff sufficiently high (i.e. 20% of the simulated flux - here it is
        // the absolute cutoff) to ignore low-level rumble and a bit of aliasing
        SupportSearcher ss(20);
        ss.searchCentered(buf);
        CPPUNIT_ASSERT(ss.support() < 32u);
        //scimath::saveAsCasaImage("tst.img", casacore::real(buf));
      }

      void testBuildingUVWeight() {
        // set up the mocked up builder
        boost::shared_ptr<UVWeightBuilderTester> wtBuilder(new UVWeightBuilderTester());
        itsSphFunc->setUVWeightBuilder(wtBuilder);
        CPPUNIT_ASSERT_EQUAL(casacore::IPosition(3,0,0,0), wtBuilder->wtShape());       

        // now, a normal gridding job as for the individual gridder tests
        casacore::Array<imtype> newImg(itsModel->shape(), static_cast<imtype>(0.));
        itsSphFunc->initialiseGrid(*itsAxes, newImg.shape(), false);
        CPPUNIT_ASSERT_EQUAL(newImg.shape().getFirst(3), wtBuilder->wtShape());       

        itsSphFunc->grid(*idi);
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1u), wtBuilder->itsBeams.size());
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1u), wtBuilder->itsFields.size());
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1u), wtBuilder->itsSourceIndices.size());
        CPPUNIT_ASSERT_EQUAL(0u, *wtBuilder->itsBeams.begin());
        CPPUNIT_ASSERT_EQUAL(0u, *wtBuilder->itsFields.begin());
        CPPUNIT_ASSERT_EQUAL(0u, *wtBuilder->itsSourceIndices.begin());
        itsSphFunc->finaliseGrid(newImg);

        TestUVWeightCalculator calc(newImg.shape().getFirst(2));
        UVWeightCollection& wtCollection = wtBuilder->finalise(calc);
        CPPUNIT_ASSERT(wtCollection.exists(0u));
        casacore::Cube<float> &wt = wtCollection.get(0u);
        CPPUNIT_ASSERT_EQUAL(newImg.shape().getFirst(3), wt.shape());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(100., wt(wt.nrow() / 2, wt.ncolumn() / 2, 0), 1e-6);
        // scimath::saveAsCasaImage("tst.img", wt.xyPlane(0));

        // again, use the support searcher to assess the weight grid. It is pretty much the
        // uv-coverage of mocked up data accessor. There is no low-level rumble here, so set the
        // absolute cutoff just below 1.0 to catch the bounding box of the uv coverage and then test blc, trc
        SupportSearcher ss(0.0009);
        ss.search(wt.xyPlane(0));
        CPPUNIT_ASSERT_EQUAL(casacore::IPosition(2,wt.nrow() / 2, wt.ncolumn() / 2), ss.peakPos());
        CPPUNIT_ASSERT_DOUBLES_EQUAL(100., ss.peakVal(), 1e-6);
        // MV: something isn't right with the support searcher in the asymmetric case. 
        // it is worth checking at some point to ensure no bug / secret assuption is lurking 
        // allow a generous buffer in the test condition below
        ss.searchCentered(wt.xyPlane(0));
        casacore::uInt sz = ss.symmetricalSupport(wt.shape().getFirst(2));
        // need to figure out why the number doesn't match manual calculation from the image
        CPPUNIT_ASSERT(sz > 160u && sz < 190u);
      }
    
      void testUVWeightBuilderInit() {
        // we rely on correct override of gridder initialisation behaviour for the weight builder initialisation to work
        // it is worth testing it for other gridders (tests above are only done for the SphFunc gridder)

        // set up the mocked up builder
        boost::shared_ptr<UVWeightBuilderTester> wtBuilder(new UVWeightBuilderTester());
        itsAWProject->setUVWeightBuilder(wtBuilder);
        CPPUNIT_ASSERT_EQUAL(casacore::IPosition(3,0,0,0), wtBuilder->wtShape());       

        itsAWProject->initialiseGrid(*itsAxes, itsModel->shape(), false);
        CPPUNIT_ASSERT_EQUAL(itsModel->shape().getFirst(3), wtBuilder->wtShape());       

        // reset to a new object to ensure the shape is not set
        wtBuilder.reset(new UVWeightBuilderTester());
        itsWProject->setUVWeightBuilder(wtBuilder);
        CPPUNIT_ASSERT_EQUAL(casacore::IPosition(3,0,0,0), wtBuilder->wtShape());       

        itsWProject->initialiseGrid(*itsAxes, itsModel->shape(), false);
        CPPUNIT_ASSERT_EQUAL(itsModel->shape().getFirst(3), wtBuilder->wtShape());       

        // reset to a new object to ensure the shape is not set
        wtBuilder.reset(new UVWeightBuilderTester());
        itsWStack->setUVWeightBuilder(wtBuilder);
        CPPUNIT_ASSERT_EQUAL(casacore::IPosition(3,0,0,0), wtBuilder->wtShape());       

        itsWStack->initialiseGrid(*itsAxes, itsModel->shape(), false);
        CPPUNIT_ASSERT_EQUAL(itsModel->shape().getFirst(3), wtBuilder->wtShape());       

        // reset to a new object to ensure the shape is not set
        wtBuilder.reset(new UVWeightBuilderTester());
        itsAProjectWStack->setUVWeightBuilder(wtBuilder);
        CPPUNIT_ASSERT_EQUAL(casacore::IPosition(3,0,0,0), wtBuilder->wtShape());       

        itsAProjectWStack->initialiseGrid(*itsAxes, itsModel->shape(), false);
        CPPUNIT_ASSERT_EQUAL(itsModel->shape().getFirst(3), wtBuilder->wtShape());       
      }

      void testUVWeightGridder() {
        // this test could really be in a separate file as it is not about TableVisGridder, but it is handy to have it here
        // to be able to reuse setup done for other tests

        // set up two identical builders, one to be used with normal gridder and another with specialised UVWeightGridder
        const boost::shared_ptr<GenericUVWeightBuilder> genericBuilder1(new GenericUVWeightBuilder(0u, 0u, 0u));
        const boost::shared_ptr<GenericUVWeightBuilder> genericBuilder2(new GenericUVWeightBuilder(0u, 0u, 0u));

        itsBox->setUVWeightBuilder(genericBuilder1);

        // now, a normal gridding job as for the individual gridder tests
        itsBox->initialiseGrid(*itsAxes, itsModel->shape(), false);
        itsBox->grid(*idi);

        // and the same via the specialised weight griddder which doesn't have overheads of the normal gridder
        UVWeightGridder wtg;
        wtg.setUVWeightBuilder(genericBuilder2);
        wtg.initialise(*itsAxes, itsModel->shape());
        wtg.accumulate(*idi);
        // simple calculator (instead of writing a stubbed version we can run the real code and test it as well)
        const ConjugatesAdderFFT calc;
        const UVWeightCollection& wts1 = genericBuilder1->finalise(calc);
        const UVWeightCollection& wts2 = genericBuilder2->finalise(calc);
        // get weight for the 0th plane, this should be the only plane present in either collection
        CPPUNIT_ASSERT(wts1.exists(0u));
        CPPUNIT_ASSERT(wts2.exists(0u));
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1u), wts1.indices().size());
        CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1u), wts2.indices().size());

        UVWeight wt1 = wts1.get(0u);
        UVWeight wt2 = wts2.get(0u);
        // compare the content of two weight grids
        CPPUNIT_ASSERT_EQUAL(wt1.uSize(), wt2.uSize());
        CPPUNIT_ASSERT_EQUAL(wt1.vSize(), wt2.vSize());
        CPPUNIT_ASSERT_EQUAL(wt1.nPlane(), wt2.nPlane());
        CPPUNIT_ASSERT_EQUAL(1u, wt1.nPlane());
 
        for (casacore::uInt iu = 0u; iu < wt1.uSize(); ++iu) {
             for (casacore::uInt iv = 0u; iv < wt1.vSize(); ++iv) {
                  // note, need to use the box gridder above, otherwise, it looks like the condition for which I submitted AXA-2485 
                  // is triggered often and as a result weight grids don't match. The other option may be to use
                  // a mosaicing gridder and set oversampling factor to 1.
                  CPPUNIT_ASSERT_DOUBLES_EQUAL(wt1(iu,iv,0u), wt2(iu,iv,0u), 1e-6);
             }
        }
      }

      // test selection of representative beam and field for the uv-weight gridder 
      // (specialised gridder used in traditional weighting). Note, it may be more correct to have uv-weight gridder tested in 
      // its own file, but it is handy to reuse infrastructure related to tests of gridder-based approaches (like simulated accessors, etc)
      void testUVWeightGridderBeamAndFieldSelection() {
         // this test needs to modify the content of the test accessor, so get the appropriate reference and work with it
         boost::shared_ptr<accessors::DataIteratorStub> di = idi.dynamicCast<accessors::DataIteratorStub>();
         CPPUNIT_ASSERT(di);
         // it looks like stubbed accessor field is public in the iterator stub. Probably ok, given this class exists largely for tests. 
         // But it would be better to have a separate method to get stubbed accessor and keep the field itself private 
         // (and, in principle, the appropriate reference can be obtained via standard interface and dynamic cast).
         accessors::DataAccessorStub& accStub = di->itsAccessor;

         // now setup two uv-weight gridders, one which ignores all indices and the other which takes them into account

         // this builder ingores indices (because all coefficients are zero)
         const boost::shared_ptr<GenericUVWeightBuilder> genericBuilder1(new GenericUVWeightBuilder(0u, 0u, 0u));

         // this builder maps beam, field, source to the grid index
         const boost::shared_ptr<GenericUVWeightBuilder> genericBuilder2(new GenericUVWeightBuilder(1u, 10u, 100u));

         UVWeightGridder wtg1;
         wtg1.setUVWeightBuilder(genericBuilder1);

         UVWeightGridder wtg2;
         wtg2.setUVWeightBuilder(genericBuilder2);

         // we don't really use the 3rd index (named "source") now, but this is a nice opportunity to test it to ensure the code works as expected
         // (the first one should be ignored, the second one will cause all indices to be incremented by 100u).
         wtg1.setSourceIndex(2u);
         wtg2.setSourceIndex(1u);

         // initialise both gridders (and builders indirectly)
         wtg1.initialise(*itsAxes, itsModel->shape());
         wtg2.initialise(*itsAxes, itsModel->shape());

         // first, accumulate the default accessor by both gridders
         wtg1.accumulate(accStub);
         wtg2.accumulate(accStub);
        
         // test the current state. Note, although getWeight method does not really have to be implemented for builders (nor,
         // the builder is required to be derived from the accessor interface, we have it this way for the GenericUVWeightBuilder class. 
         // This helps with the test because we don't have to finalise to check the state. And if index is wrong an exception will be thrown.
         // The alternative is to call finalise multiple times (which we can, technically), although this is, strictly speaking, not the
         // expected use case either. But the appropriate calculator class can be written for the test which can leave the weight intact
         // (and therefore, the next accumulate call easy to understand)
         {
            UVWeight wt1 = genericBuilder1->getWeight(0u,0u,0u);
            UVWeight wt2 = genericBuilder2->getWeight(0u,0u,1u);
            CPPUNIT_ASSERT_EQUAL(wt1.uSize(), wt2.uSize());
            CPPUNIT_ASSERT_EQUAL(wt1.vSize(), wt2.vSize());
            CPPUNIT_ASSERT_EQUAL(wt1.nPlane(), wt2.nPlane());
            for (casacore::uInt plane = 0; plane < wt1.nPlane(); ++plane) {
                 for (casacore::uInt u = 0; u < wt1.uSize(); ++u) {
                      for (casacore::uInt v = 0; v < wt1.vSize(); ++v) {
                           CPPUNIT_ASSERT_DOUBLES_EQUAL(wt1(u,v,plane), wt2(u,v,plane), 1e-6);
                      }
                 }
            }
         }
         
         // all results have been checked, so run finalise and explore the content of each collection
         // (to ensure that all indices are accounted for). For this particular test we don't
         // care about the calculator class type or actual weights (they're already checked above). 
         // Therefore, use TestUVWeightCalculator which just checks that the shape matches the 
         // expected shape passed at construction and sets the central pixel to 100.
         TestUVWeightCalculator calc(itsModel->shape().getFirst(2));
         const UVWeightCollection& wts1 = genericBuilder1->finalise(calc);
         const UVWeightCollection& wts2 = genericBuilder2->finalise(calc);

         // check existance of the appropriate planes and the correct size of the collection
         CPPUNIT_ASSERT(wts1.exists(0u));
         CPPUNIT_ASSERT(wts2.exists(100u));
         CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1u), wts1.indices().size());
         CPPUNIT_ASSERT_EQUAL(static_cast<size_t>(1u), wts2.indices().size());
      }

      void testReverseAWProject()
      {
        itsAWProject->initialiseGrid(*itsAxes, itsModel->shape(), false);
        itsAWProject->grid(*idi);
        itsAWProject->finaliseGrid(*itsModel);
        itsAWProject->finaliseWeights(*itsModelWeights);
        itsAWProject->initialiseGrid(*itsAxes, itsModel->shape(), true);
        itsAWProject->grid(*idi);
        itsAWProject->finaliseGrid(*itsModelPSF);
      }
      void testForwardAWProject()
      {
        itsAWProject->initialiseDegrid(*itsAxes, *itsModel);
        itsAWProject->degrid(*idi);
      }
      void testReverseWProject()
      {
        itsWProject->initialiseGrid(*itsAxes, itsModel->shape(), false);
        itsWProject->grid(*idi);
        itsWProject->finaliseGrid(*itsModel);
        itsWProject->finaliseWeights(*itsModelWeights);
        itsWProject->initialiseGrid(*itsAxes, itsModel->shape(), true);
        itsWProject->grid(*idi);
        itsWProject->finaliseGrid(*itsModelPSF);
      }
      void testForwardWProject()
      {
        itsWProject->initialiseDegrid(*itsAxes, *itsModel);
        itsWProject->degrid(*idi);
      }
      void testReverseWStack()
      {
        itsWStack->initialiseGrid(*itsAxes, itsModel->shape(), false);
        itsWStack->grid(*idi);
        itsWStack->finaliseGrid(*itsModel);
        itsWStack->finaliseWeights(*itsModelWeights);
        itsWStack->initialiseGrid(*itsAxes, itsModel->shape(), true);
        itsWStack->grid(*idi);
        itsWStack->finaliseGrid(*itsModelPSF);
      }
      void testForwardWStack()
      {
        itsWStack->initialiseDegrid(*itsAxes, *itsModel);
        itsWStack->degrid(*idi);
      }
      void testReverseAProjectWStack()
      {
        itsAProjectWStack->initialiseGrid(*itsAxes, itsModel->shape(), false);
        itsAProjectWStack->grid(*idi);
        itsAProjectWStack->finaliseGrid(*itsModel);
        itsAProjectWStack->finaliseWeights(*itsModelWeights);
        itsAProjectWStack->initialiseGrid(*itsAxes, itsModel->shape(), true);
        itsAProjectWStack->grid(*idi);
        itsAProjectWStack->finaliseGrid(*itsModelPSF);
      }
      void testForwardAProjectWStack()
      {
        itsAProjectWStack->initialiseDegrid(*itsAxes, *itsModel);
        itsAProjectWStack->degrid(*idi);
      }
      void testReverseATCAIllumination()
      {
        boost::shared_ptr<ATCAIllumination> illum(new ATCAIllumination(120.0, 10.0));
        illum->simulateTapering(1.0);
        itsAProjectWStack.reset(new AProjectWStackVisGridder(illum, 10000.0, 9, 0., 1, 128, 1));
        itsAProjectWStack->initialiseGrid(*itsAxes, itsModel->shape(), false);
        itsAProjectWStack->grid(*idi);
        itsAProjectWStack->finaliseGrid(*itsModel);
        itsAProjectWStack->finaliseWeights(*itsModelWeights);
        itsAProjectWStack->initialiseGrid(*itsAxes, itsModel->shape(), true);
        itsAProjectWStack->grid(*idi);
        itsAProjectWStack->finaliseGrid(*itsModelPSF);
      }
      void testForwardATCAIllumination()
      {
        boost::shared_ptr<ATCAIllumination> illum(new ATCAIllumination(12.0, 2.0));
        illum->simulateTapering(1.0);
        itsAProjectWStack.reset(new AProjectWStackVisGridder(illum, 10000.0, 9, 0., 1, 128, 1));
        itsAProjectWStack->initialiseDegrid(*itsAxes, *itsModel);
        itsAProjectWStack->degrid(*idi);
      }
    };

  }
}
