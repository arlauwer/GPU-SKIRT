/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "MonteCarloSimulation.hpp"
#include "CartesianSpatialGrid.hpp"
#include "Log.hpp"
#include "MediumSystem.hpp"
#include "Parallel.hpp"
#include "ParallelFactory.hpp"
#include "PhotonPackets.hpp"
#include "Position.hpp"
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
        auto cart = find<CartesianSpatialGrid>();
        _kernel = new SimulationKernel(cart);

        PhotonPackets& photons = _kernel->photons();
        size_t currentBatch;
        for (size_t firstIndex = 0; firstIndex < Npp; firstIndex += currentBatch)
        {
            // batch is at most the total photons left
            currentBatch = std::min(Nbp, Npp - firstIndex);
            photons.setBatchSize(currentBatch);

            // launch photon packets in parallel
            // this prepares the data-oriented structure, PhotonPackets, for the GPU kernel
            parallel->call(currentBatch,
                           [this, &photons, firstIndex](size_t i, size_t n) { launch(photons, firstIndex + i, n); });

            // GPU kernel with prepared photon packets
            _kernel->runBatch();
        }

        delete _kernel;  // why is _config not deleted?
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
    auto cart = dynamic_cast<CartesianSpatialGrid*>(mediumSystem()->grid());
    for (size_t b = 0; b != numIndices; ++b)
    {
        ss->launch(pp, firstIndex + b, b);

        // set initial cell index
        double x = pp.rxv[b];
        double y = pp.ryv[b];
        double z = pp.rzv[b];
        Position r(x, y, z);

        int i, j, k;
        cart->cellIndices(i, j, k, r);
        int m = cart->index(i, j, k);

        // set cell indices
        pp.iv[b] = i;
        pp.jv[b] = j;
        pp.kv[b] = k;
        pp.mv[b] = m;
    }
}

////////////////////////////////////////////////////////////////////
