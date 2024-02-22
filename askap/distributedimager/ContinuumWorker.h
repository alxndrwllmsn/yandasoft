/// @file ContinuumWorker.h
///
/// @copyright (c) 2009 CSIRO
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
/// @author Ben Humphreys <ben.humphreys@csiro.au>
/// @author Stephen Ord <stephen.ord@csiro.au>

#ifndef ASKAP_CP_SIMAGER_CONTINUUMWORKER_H
#define ASKAP_CP_SIMAGER_CONTINUUMWORKER_H

// System includes
#include <string>

// boost includes
#include <boost/noncopyable.hpp>

// ASKAPsoft includes
#include "boost/shared_ptr.hpp"
#include <Common/ParameterSet.h>
#include <askap/scimath/fitting/INormalEquations.h>
#include <askap/scimath/fitting/Params.h>
#include <askap/dataaccess/TableDataSource.h>
#include <askap/gridding/IVisGridder.h>
#include <askap/askap/StatReporter.h>

// Local includes
#include "askap/distributedimager/AdviseDI.h"
#include "askap/distributedimager/MSSplitter.h"
#include "askap/distributedimager/CalcCore.h"
#include "askap/messages/ContinuumWorkUnit.h"
#include "askap/distributedimager/WorkUnitContainer.h"
#include "askap/distributedimager/CubeBuilder.h"
#include "askap/distributedimager/CubeComms.h"
#include <askap/utils/StatsAndMask.h>

namespace askap {
namespace cp {

class ContinuumWorker : public boost::noncopyable

    {
    public:
        ContinuumWorker(LOFAR::ParameterSet& parset,
                           CubeComms& comms, StatReporter& stats);

        void run(void);

        void writeCubeStatistics();

    private:
        /// @brief configure allocation in channels
        /// @details This method sets up channel allocation for cube writing (in local solver mode)
        /// or combines channels in the global solver mode if configured in the parset.
        void configureChannelAllocation();

        /// @brief configure reference channel used for the restoring beam
        /// @details This method populates itsBeamReferenceChannel based on the parset and
        /// the number of channels allocated to the cube handled by this rank
        void configureReferenceChannel();

        /// @brief initialise cube writing if in local solver mode
        /// @details This method encapsulates the code which handles cube writing in the local solver mode
        /// (i.e. when it is done from the worker). Safe to call in continuum mode too (as itsComms.isWriter
        /// would return false in this case)
        void initialiseCubeWritingIfNecessary();

        /// @brief check stopping thresholds in the model 
        /// @details This method is used at the end of minor cycle deconvolution to check whether to continue iterations.
        /// @param[in] model shared pointer to the scimath::Params object with the model
        /// @return true if stopping is required
        /// @note We do similar checks in both the master and in workers. So this method can be moved somewhere else to be shared.
        bool checkStoppingThresholds(const boost::shared_ptr<scimath::Params> &model) const;


        /// @brief add current image as a model
        /// @details This method adds fullres (if present) or ordinary image as model.slice in the given params object.
        /// It is expected that model.slice will be absent and either fullres or Nyquist resolution image should be
        /// present.
        /// @param[in] params shared pointer to the params object to work with (should be non-empty)
        /// @note I (MV) think there could be untidy design here - we probably make an extra copy which could be avoided
        static void addImageAsModel(const boost::shared_ptr<scimath::Params> &params);

        /// @brief perform one write job for a remote client
        /// @details This method is expected to be used for cube writing ranks only. It receives a single
        /// write job and performs it.
        void performSingleWriteJob(); 

        /// @brief cleanup outstanding write jobs
        /// @details This method is used for cube writing ranks (does nothing for non-writers), it loops over all outstanding
        /// work units and performs write operation assigned to this rank.
        /// @param[in] targetOutstanding desired number of outstanding write jobs at the end of execution
        ///                              (to spread writing across the iteration), default is all jobs
        /// @param[in] minOutstanding    minimal number of outstanding jobs to remain (default - none)
        /// @note (MV) I didn't fully understand the logic behind targetOutstanding and minOutstanding (one should be
        /// sufficient), the same behaviour as we had prior to refactoring has been implemented.
        void performOutstandingWriteJobs(int targetOutstanding = 0, int minOutstanding = -1);

        /// @brief perform write job allocated to this rank
        /// @details unlike performOutstandingWriteJobs or performSingleWriteJob this method deals with the write
        /// job handled entirely by this rank (i.e. its own write job) and, hence, provided explicitly rather than
        /// received from another rank.
        /// @param[in] globalChannel global channel (i.e. channel in the whole cube) to write
        /// @param[in] params shared pointer to the model with required info (should not be empty)
        void performOwnWriteJob(unsigned int globalChannel, const boost::shared_ptr<scimath::Params> &params);

        /// @brief send blank image to writer
        /// @details This method is expected to be used when calculation of a spectral plane is failed for some reason,
        /// but some other rank is responsible for writing it. Essentially it sends a blank image with parameters 
        /// (like frequency and channel) filled from the work unit. 
        /// @param[in] wu work unit to take the information from
        /// @note (MV:) This doesn't seem like a good design, but the behaviour is left the same as it was prior to
        /// the refactoring.
        void sendBlankImageToWriter(const cp::ContinuumWorkUnit &wu) const;

        /// @brief figure out if preconditioning is to be done
        /// @details This method encapsulates checks of the parset indicating that preconditioning is going to be done.
        /// It is necessary to configure writing of additional data products, although preconditioning itself 
        /// is enabled and done by the appropriate solver class.
        /// @param[in] parset parset to use (this method is expected to be used in the constructor, so it is handy not
        /// to rely on the itsParset data field).
        /// @return true if preconditioning is to be done, false otherwise
        static bool doingPreconditioning(const LOFAR::ParameterSet &parset);

        /// @brief helper method to obtain the number of writers for the cube
        /// @details This method obtains the number of writers from the parset and adjusts it if necessary.
        /// It is intended to be used in the constructor to fill itsNumWriters data field and requires 
        /// itsGridType and itsParset to be valid.
        /// @return number of writer ranks for grid export
        int configureNumberOfWriters();

        // My Advisor
        boost::shared_ptr<synthesis::AdviseDI> itsAdvisor;

        /// @brief the work units
        WorkUnitContainer itsWorkUnits;

        /// @brief whether preconditioning has been requested
        /// @details It is populated in the constructor based on the parset using doingPreconditioning method
        const bool itsDoingPreconditioning;

        // Process a workunit
        void preProcessWorkUnit(ContinuumWorkUnit& wu);

        // Compress all continuous channel allocations into individual workunits
        void compressWorkUnits();

        //For all workunits .... process
        void processChannels();

        // Setup the image specified in parset and add it to the Params instance.
        void setupImage(const askap::scimath::Params::ShPtr& params,
                    double channelFrequency, bool shapeOveride = false) const;

        void buildSpectralCube();

        // Root Parameter set good for information common to all workUnits
        LOFAR::ParameterSet& itsParset;

        // Communications class
        CubeComms& itsComms;

        // statistics
        StatReporter& itsStats;

        // ID of the master process
        static const int itsMaster = 0;

        // the basechannel number assigned to this worker
        unsigned int itsBaseChannel;

        // the baseFrequency associated with this channel
        double itsBaseFrequency;

        // the baseFrequency associated with the cube if being built
        double itsBaseCubeFrequency;

        // the global channel associated with this part of the cube
        int itsBaseCubeGlobalChannel;

        // the number of channels in this cube (if writer)
        int itsNChanCube;

        /// @brief true if solver is run locally (spectral line mode), 
        /// false for central solver (continuum). 
        const bool itsLocalSolver;

        boost::shared_ptr<CubeBuilder<casacore::Float> > itsImageCube;
        boost::shared_ptr<CubeBuilder<casacore::Float> > itsPSFCube;
        boost::shared_ptr<CubeBuilder<casacore::Complex> > itsPCFGridCube;
        boost::shared_ptr<CubeBuilder<casacore::Complex> > itsPSFGridCube;
        boost::shared_ptr<CubeBuilder<casacore::Complex> > itsVisGridCube;
        boost::shared_ptr<CubeBuilder<casacore::Float> > itsPCFGridCubeReal;
        boost::shared_ptr<CubeBuilder<casacore::Float> > itsPSFGridCubeReal;
        boost::shared_ptr<CubeBuilder<casacore::Float> > itsVisGridCubeReal;
        boost::shared_ptr<CubeBuilder<casacore::Float> > itsPCFGridCubeImag;
        boost::shared_ptr<CubeBuilder<casacore::Float> > itsPSFGridCubeImag;
        boost::shared_ptr<CubeBuilder<casacore::Float> > itsVisGridCubeImag;
        boost::shared_ptr<CubeBuilder<casacore::Float> > itsResidualCube;
        boost::shared_ptr<CubeBuilder<casacore::Float> > itsWeightsCube;
        boost::shared_ptr<CubeBuilder<casacore::Float> > itsPSFimageCube;
        boost::shared_ptr<CubeBuilder<casacore::Float> > itsRestoredCube;
        boost::shared_ptr<askap::utils::StatsAndMask> itsRestoredStatsAndMask;
        boost::shared_ptr<askap::utils::StatsAndMask> itsResidualStatsAndMask;
        std::string itsWeightsName;

        void handleImageParams(askap::scimath::Params::ShPtr params, unsigned int chan);

        void copyModel(askap::scimath::Params::ShPtr SourceParams, askap::scimath::Params::ShPtr SinkParams) const;

        void initialiseBeamLog(const unsigned int numChannels);
        void recordBeam(const askap::scimath::Axes &axes, const unsigned int globalChannel);

        std::map<unsigned int, casacore::Vector<casacore::Quantum<double> > > itsBeamList;
        unsigned int itsBeamReferenceChannel;
        void logBeamInfo() const;

        void initialiseWeightsLog(const unsigned int numChannels);
        void recordWeight(float wt, const unsigned int globalChannel);
        std::map<unsigned int, float> itsWeightsList;
        void logWeightsInfo() const;

        /// @brief Do we want a restored image?
        const bool itsRestore;

        /// @brief Do we want a residual image
        const bool itsWriteResidual;

        /// @brief write 'raw', unnormalised, natural weight psf
        const bool itsWritePsfRaw;

        /// @brief write normalised, preconditioned psf
        const bool itsWritePsfImage;

        /// @brief write weights image
        const bool itsWriteWtImage;

        /// @brief write a weights log
        const bool itsWriteWtLog;

        /// @brief write out the (clean) model image
        const bool itsWriteModelImage;

        /// @brief write out the gridded data, pcf and psf
        const bool itsWriteGrids;

        /// @brief grid image type
        const std::string itsGridType;

        /// @brief write out the grids with UV coordinate grid
        const bool itsGridCoordUV;

        /// @brief write out the FFT of the grids
        const bool itsGridFFT;

        /// @brief the number of rank that can write to the cube
        const int itsNumWriters;

};

};
};

#endif
