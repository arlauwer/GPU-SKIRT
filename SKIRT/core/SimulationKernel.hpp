/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SIMULATIONKERNEL_HPP
#define SIMULATIONKERNEL_HPP

#include "CartesianSpatialGrid.hpp"
#include "PhotonPackets.hpp"
#include "SourceSystem.hpp"

class MediumSystem;

class SimulationKernel
{
public:
    SimulationKernel(SourceSystem* ss, MediumSystem* ms);

    virtual ~SimulationKernel();

    void runBatch();

    PhotonPackets& photons() { return _photons; }

private:
    SourceSystem* _ss;
    MediumSystem* _ms;

    PhotonPackets _photons;
    PhotonPackets _pphotons;

    // cartesian hard coded grid
    int _Nx{0};
    int _Ny{0};
    int _Nz{0};
    Array _gxv;
    Array _gyv;
    Array _gzv;

    // source wavelength grid
    int _Nss{0};
    double* _ss_lambdav{nullptr};
    double* _crossv{nullptr};
};

#endif
