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

    // cartesian grid geometry
    int _Nx{0};
    int _Ny{0};
    int _Nz{0};
    Array _gxv;  // x cell borders [Nx+1]
    Array _gyv;  // y cell borders [Ny+1]
    Array _gzv;  // z cell borders [Nz+1]

    // number densities per cell, indexed by cell index m [Ncell]
    size_t _Ncell{0};
    double* _nv{nullptr};

    // radiation field accumulation buffer, indexed as [m * Nrad + ell]
    size_t _Nrad{0};
    double* _rf1{nullptr};

    // log-spaced wavelength grid for cross section lookup on GPU
    // index lookup: ell = (log(lambda) - _xsec_logLambdaMin) * _xsec_logInvBinWidth
    // cross section at index ell: _crossv[ell]
    int _Nxsec{0};
    double _xsec_logLambdaMin{0.0};    // log(lambdaMin), cached for GPU lookup
    double _xsec_logInvBinWidth{0.0};  // (Nxsec-1) / log(lambdaMax/lambdaMin), cached for GPU lookup
    double* _crossv{nullptr};          // extinction cross section per wavelength bin [Nxsec]
};

#endif
