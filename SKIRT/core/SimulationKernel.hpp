/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SIMULATIONKERNEL_HPP
#define SIMULATIONKERNEL_HPP

#include "CartesianSpatialGrid.hpp"
#include "PhotonPackets.hpp"

class SimulationKernel
{
public:
    SimulationKernel(const CartesianSpatialGrid* grid);

    void runBatch();

    PhotonPackets& photons() { return _photons; }

private:
    PhotonPackets _photons;

    // cartesian hard coded grid
    int _Nx{0};
    int _Ny{0};
    int _Nz{0};
    std::vector<double> _gxv;
    std::vector<double> _gyv;
    std::vector<double> _gzv;
};

#endif
