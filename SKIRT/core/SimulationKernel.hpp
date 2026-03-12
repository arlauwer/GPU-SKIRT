/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SIMULATIONKERNEL_HPP
#define SIMULATIONKERNEL_HPP

#include "CartesianSpatialGrid.hpp"
#include "Random.hpp"
#include "PhotonPackets.hpp"
#include "SourceSystem.hpp"

class MediumSystem;

class SimulationKernel
{
public:
    SimulationKernel(SourceSystem* ss, MediumSystem* ms);

    virtual ~SimulationKernel();

    void runBatch();

	void traverse(PhotonPackets& pp);

    PhotonPackets& photons() { return _photons; }

private:
    SourceSystem* _ss;
    MediumSystem* _ms;

    PhotonPackets _photons;
    PhotonPackets _pphotons;

    // cartesian grid geometry
    size_t _Nx{0};
    size_t _Ny{0};
    size_t _Nz{0};
    Array _gxv;  // x cell borders [Nx+1]
    Array _gyv;  // y cell borders [Ny+1]
    Array _gzv;  // z cell borders [Nz+1]

    // number densities per cell, indexed by [m]
    size_t _Ncell{0};
    double* _nv{nullptr};

    // tabulated radiation field, indexed by [m * Nrad + ell]
    size_t _Nrad{0};
    double* _rad{nullptr};
    // log-spaced wavelength grid for radiation field lookups
    double _rad_logLambdaMin{0.0};
    double _rad_logInvBinWidth{0.0};

    // tabulated extinction cross sections
    size_t _Nsec{0};
    double* _crossv{nullptr};
    // log-spaced wavelength grid for cross section lookups
    double _sec_logLambdaMin{0.0};
    double _sec_logInvBinWidth{0.0};

	Random* _random;
};

#endif
