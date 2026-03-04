/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MonteCarloSimulation.hpp"
#include "Log.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PhotonPackets.hpp"
#include "SourceSystem.hpp"
#include "StringUtils.hpp"
#include "TimeLogger.hpp"
#include <cstdio>

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::setupSimulation()
{
    // perform regular setup for the hierarchy and wait for all processes to finish
    {
        TimeLogger logger(log(), "setup");
        _config->setup();  // first of all perform setup for the configuration object
        SimulationItem::setup();
    }

    // write setup output
    {
        TimeLogger logger(log(), "setup output");

        // notify the probe system
        probeSystem()->probeSetup();
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::setupSelfAfter()
{
    Simulation::setupSelfAfter();
}

////////////////////////////////////////////////////////////////////

Configuration* MonteCarloSimulation::config() const
{
    return _config;
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::runSimulation()
{
    // run the simulation
    {
        TimeLogger logger(log(), "the run");

        runPrimaryEmission();
    }

    // write final output
    {
        TimeLogger logger(log(), "final output");

        // notify the probe system
        probeSystem()->probeRun();

        // write instrument output
        instrumentSystem()->write();
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::runPrimaryEmission()
{
    string segment = "primary emission";
    TimeLogger logger(log(), segment);

    // shoot photons from primary sources, if needed
    size_t Npp = _config->numPrimaryPackets();
    size_t Nbp = _config->numBatchPackets();
    if (!Npp)
    {
        log()->warning("Skipping primary emission because no photon packets were requested");
    }
    else if (!sourceSystem()->luminosity())
    {
        log()->warning("Skipping primary emission because the total luminosity of primary sources is zero");
    }
    else
    {
        initProgress(segment, Npp);
        sourceSystem()->prepareForLaunch(Npp);

        auto parallel = find<ParallelFactory>()->parallel();

        PhotonPackets pp;
        size_t currentBatch;
        for (size_t firstIndex = 0; firstIndex < Npp; firstIndex += currentBatch)
        {
            // batch is at most the total photons left
            currentBatch = std::min(Nbp, Npp - firstIndex);
            pp.setBatchSize(currentBatch);

            // launch photon packets in parallel
            // this prepares the data-oriented structure, PhotonPackets, for the GPU kernel
            parallel->call(currentBatch,
                           [this, &pp, &firstIndex](size_t i, size_t n) { launch(pp, firstIndex + i, n); });

            // GPU kernel with prepared photon packets
        }

        // flush intruments
    }
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::initProgress(string segment, size_t numTotal)
{
    _segment = segment;

    log()->info("Launching " + StringUtils::toString(static_cast<double>(numTotal)) + " " + _segment
                + " photon packets");
    log()->infoSetElapsed(numTotal);
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::logProgress(size_t numDone)
{
    // log message if the minimum time has elapsed
    log()->infoIfElapsed("Launched " + _segment + " photon packets: ", numDone);
}

////////////////////////////////////////////////////////////////////

void MonteCarloSimulation::launch(PhotonPackets& pp, size_t firstIndex, size_t numIndices)
{
    auto ss = sourceSystem();
    for (size_t i = 0; i != numIndices; ++i)
    {
        ss->launch(pp, firstIndex + i, i);
    }
}

////////////////////////////////////////////////////////////////////
