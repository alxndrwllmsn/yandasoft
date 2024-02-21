/// @file CalcCore.h
///
/// @copyright (c) 2016 CSIRO
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
/// @author Stephen Ord <stephen.ord@csiro.au>
/// based upon SolverCore by
/// Ben Humphreys <ben.humphreys@csiro.au>

#ifndef ASKAP_CP_ASKAP_IMAGER_CALCCORE_H
#define ASKAP_CP_ASKAP_IMAGER_CALCCORE_H

// System includes
#include <string>

// ASKAPsoft includes
#include <Common/ParameterSet.h>
#include <askap/scimath/fitting/INormalEquations.h>
#include <askap/scimath/fitting/Params.h>
#include <askap/scimath/fitting/Solver.h>
#include <askap/scimath/fitting/Quality.h>
#include <askap/parallel/ImagerParallel.h>
#include <askap/gridding/IVisGridder.h>
#include <askap/gridding/VisGridderFactory.h>

#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/casa/Arrays/Vector.h>
#include <askap/dataaccess/IDataSource.h>
#include <askap/dataaccess/SharedIter.h>
#include <askap/measurementequation/ImageFFTEquation.h>

// boost includes
#include <boost/noncopyable.hpp>

// Local includes


namespace askap {
namespace cp {

/// @brief Core Normal Equation Calculation functionality required by the imager.
    class CalcCore : public synthesis::ImagerParallel,
                     public boost::noncopyable

    {
    public:
        /// @brief Constructor
        CalcCore(LOFAR::ParameterSet& parset,
                   askap::askapparallel::AskapParallel& comms,
                   accessors::IDataSource& ds, int localChannel=1, double frequency=0);

        /// @brief Constructor that maintains the gridder
        CalcCore(LOFAR::ParameterSet& parset,
                askap::askapparallel::AskapParallel& comms,
                accessors::IDataSource& ds, askap::synthesis::IVisGridder::ShPtr gdr,
                 int localChannel=1, double frequency=0);

        /// @brief Calc the normal equations
        /// @detail Overrides the virtual function in the ImagerParallel base
        virtual void calcNE() override;

        /// @brief Solve the normal equations (runs in the solver)
        /// @details Runs minor cycle
        virtual void solveNE() override;

        void doCalc();

        void init();

        void updateSolver() const;

        void reset() const;

        void zero() const;

        void check() const;

        void restoreImage() const;

        void writeLocalModel(const std::string& postfix) const;

        /// @brief obtain the current gridder template
        /// @return shared pointer to the gridder which can be cloned
        askap::synthesis::IVisGridder::ShPtr gridder() const { return itsGridder;};

        /// @brief return the residual grid
        casacore::Array<casacore::Complex> getGrid() const;
        /// @brief return the PCF grid
        casacore::Array<casacore::Complex> getPCFGrid() const;
        /// @brief return the PSF grid
        casacore::Array<casacore::Complex> getPSFGrid() const;

        /// @brief store all complex grids in the model object for future writing
        /// @details This method calls getGrid, getPCFGrid and getPSFGrid and stores
        /// returned arrays in the model so they can be exported later.
        void addGridsToModel();

        /// @brief iterate over data and accumulate samples for uv weights
        /// @details This method is used to build the sample density in the uv-plane via the appropriate gridder
        /// and weight builder class. It expects the builder already setup and accessible via the normal equations 
        /// shared pointer. Unlike the variant from the base class which works with the iterator supplied as a parameter,
        /// this version uses the iterator returned by makeDataIterator (wrapped into the calibration adapter, if needed)
        void accumulateUVWeights() const;

        /// @brief configure normal equation for linear mosaicing
        /// @details When linmos is expected to happen during merge of normal equations we need to configure
        /// NEs appropriately to interpret weight correctly. This helper method does it. 
        /// @note Normal equations should already be setup (although could be empty) before this method is called.
        /// Otherwise, an exception will be thrown. Also, we could've do this setup automatically based on the 
        /// gridder type. But at the moment the same approach is followed as we had prior to refactoring.
        void configureNormalEquationsForMosaicing() const;

        /// @brief merge normal equations from another CalcCore
        /// @details This is a convenience method to merge in normal equations held by other CalcCore
        /// object. In principle, we can have this method in one of the base classes (and require 
        /// broader type rather than CalcCore as the input) because all of the required functionality is
        /// in the base classes. But we only use it with CalcCore, so keep it in this class as well.
        /// @note Normal equations should be initialised (and with the consistent type) in both
        /// this and other CalcCore instances, but could be empty. The method is const because it doesn't change
        /// this class (only changes normal equations held by pointer).
        /// @param[in] other an instance of CalcCore to merge from
        void mergeNormalEquations(const CalcCore &other) const;

    protected:
        /// @brief keep the base class' version accessible here
        /// @note for quick reference, it has the following signature:
        /// void accumulateUVWeights(const boost::shared_ptr<accessors::IConstDataIterator> &iter) const;
        using ImagerParallel::accumulateUVWeights; 

        /// @brief make data iterator
        /// @details This helper method makes an iterator based on the configuration in the current parset and
        /// data fields of this class such as itsChannel and itsFrequency
        /// @return shared pointer to the iterator over original data
        accessors::IDataSharedIter makeDataIterator() const;

        /// @brief create measurement equation 
        /// @details This method creates measurement equation as appropriate (with calibration application or without) using
        /// internal state of this class and the parset
        void createMeasurementEquation();

        /// @brief first image name in the model
        /// @details This is a helper method to obtain the name of the first encountered image parameter in the model. 
        /// @note It is written as part of the refactoring of various getGrid methods. However, in principle we could have multiple
        /// image parameters simultaneously. The original approach getting the first one won't work in this case.
        /// @return name of the first encountered image parameter in the model
        std::string getFirstImageName() const;

    private:

        // Communications class
        askap::askapparallel::AskapParallel& itsComms;

        /// @brief shared pointer to the solver
        askap::scimath::Solver::ShPtr itsSolver;

        /// @brief run restore solver?
        bool itsRestore;

        /// @brief data source to work with (essentially a measurement set)
        accessors::IDataSource& itsDataSource;

        /// @brief shared pointer to the gridder prototype
        /// @details WARNING this is cloned by the Equation - so you get little from specifying it
        askap::synthesis::IVisGridder::ShPtr itsGridder;

        // Its channel in the dataset
        int itsChannel;

        // Its frequency in the output cube
        double itsFrequency;

};
};
};

#endif
