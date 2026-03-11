#include "SimulationKernel.hpp"
#include "CartesianSpatialGrid.hpp"
#include "MediumSystem.hpp"
#include "NR.hpp"
#include <omp.h>
#include <iostream>

namespace
{
    inline double cartesianNext(const double* gxv, const double* gyv, const double* gzv, int Nx, int Ny, int Nz, int& i,
                                int& j, int& k, int& m, double& rx, double& ry, double& rz, double& kx, double& ky,
                                double& kz)
    {
        int xdir = kx < 0.0;
        int ydir = ky < 0.0;
        int zdir = kz < 0.0;

        double xE = (xdir) ? gxv[i] : gxv[i + 1];
        double yE = (ydir) ? gyv[j] : gyv[j + 1];
        double zE = (zdir) ? gzv[k] : gzv[k + 1];
        double dsx = (fabs(kx) > 1e-15) ? (xE - rx) / kx : DBL_MAX;
        double dsy = (fabs(ky) > 1e-15) ? (yE - ry) / ky : DBL_MAX;
        double dsz = (fabs(kz) > 1e-15) ? (zE - rz) / kz : DBL_MAX;

        // dsx
        if (dsx <= dsy && dsx <= dsz)
        {
            rx = xE;
            ry += ky * dsx;
            rz += kz * dsx;

            i += xdir ? -1 : 1;
            m += xdir * Nz * Ny;

            return dsx;
        }
        // dsy
        else if (dsy < dsx && dsy <= dsz)
        {
            ry = yE;
            rx += kx * dsy;
            rz += kz * dsy;

            j += ydir ? -1 : 1;
            m += ydir * Nz;

            return dsy;
        }
        // dsz
        else  // if (dsz < dsx && dsz < dsy)
        {
            rz = zE;
            rx += kx * dsz;
            ry += ky * dsz;

            k += zdir ? -1 : 1;
            m += zdir;

            return dsz;
        }
    }
}

// hard coded table index for MediumSystem::_rf1
inline size_t radIndex(size_t m, size_t ell, size_t num_ell)
{
    return m * num_ell + ell;
}

SimulationKernel::SimulationKernel(SourceSystem* ss, MediumSystem* ms) : _ss(ss), _ms(ms)
{
    const auto grid = dynamic_cast<CartesianSpatialGrid*>(_ms->grid());
    _Nx = grid->_Nx;
    _Ny = grid->_Ny;
    _Nz = grid->_Nz;
    _gxv = grid->_xv;
    _gyv = grid->_yv;
    _gzv = grid->_zv;

    // source wavelength grid
    _Nss = 1000;
    double lambdaMax = _ss->maxWavelength();
    double lambdaMin = _ss->minWavelength();
    Array ss_lambdaa;
    NR::buildLogGrid(ss_lambdaa, lambdaMin, lambdaMax, _Nss);
    _ss_lambdav = new double[ss_lambdaa.size()];
    std::copy(begin(ss_lambdaa), end(ss_lambdaa), _ss_lambdav);

    _crossv = new double[ss_lambdaa.size()];
    for (int ell = 0; ell < _Nss; ++ell)
    {
        double lam = _ss_lambdav[ell];
        _crossv[ell] = _ms->mix(0, 0)->sectionExt(lam);
    }
}

SimulationKernel::~SimulationKernel()
{
    delete[] _ss_lambdav;
    delete[] _crossv;
}

void SimulationKernel::runBatch()
{
    // grid
    size_t Nx = _Nx;
    size_t Ny = _Ny;
    size_t Nz = _Nz;
    double* gxv = std::begin(_gxv);
    double* gyv = std::begin(_gyv);
    double* gzv = std::begin(_gzv);
    // (_numCells, _wavelengthGrid->numBins())
    size_t Ncell = _ms->numCells();
    size_t Nrad = _ms->_wavelengthGrid->numBins();
    double* rf1 = std::begin(_ms->_rf1.data());

    int Nss = _Nss;
    double* ss_lambdav = _ss_lambdav;
    double* crossv = _crossv;

    // photon packets
    size_t Nb = _photons.batchSize();

    // photon packet
    double* lambdav = _photons.lambdav.data();
    double* weightv = _photons.weightv.data();
    int* iv = _photons.iv.data();
    int* jv = _photons.jv.data();
    int* kv = _photons.kv.data();
    int* mv = _photons.mv.data();
    double* rxv = _photons.rxv.data();
    double* ryv = _photons.ryv.data();
    double* rzv = _photons.rzv.data();
    double* kxv = _photons.kxv.data();
    double* kyv = _photons.kyv.data();
    double* kzv = _photons.kzv.data();

    // peel-off packet
    double* p_lambdav = _pphotons.lambdav.data();
    double* p_weightv = _pphotons.weightv.data();
    int* p_iv = _pphotons.iv.data();
    int* p_jv = _pphotons.jv.data();
    int* p_kv = _pphotons.kv.data();
    int* p_mv = _pphotons.mv.data();
    double* p_rxv = _pphotons.rxv.data();
    double* p_ryv = _pphotons.ryv.data();
    double* p_rzv = _pphotons.rzv.data();
    double* p_kxv = _pphotons.kxv.data();
    double* p_kyv = _pphotons.kyv.data();
    double* p_kzv = _pphotons.kzv.data();

#define NEXT(b, p) \
    cartesianNext(gxv, gyv, gzv, Nx, Ny, Nz, iv[b], p##jv[b], p##kv[b], p##mv[b], p##rxv[b], p##ryv[b], p##rzv[b], \
                  p##kxv[b], p##kyv[b], p##kzv[b]);

#pragma omp target data map(to : gxv[0 : Nx + 1], gyv[0 : Ny + 1], gzv[0 : Nz + 1]) \
    map(to : ss_lambdav[0 : Nrad], crossv[0 : Nrad]) MAP_PACKETS(Nb, ) MAP_PACKETS(Nb, p_) \
    map(tofrom : rf1[0 : Ncell * Nrad])
    {
#pragma omp target teams distribute parallel for firstprivate(Nx, Ny, Nz, Nb, Ncell, Nrad, Nss)
        for (size_t b = 0; b != Nb; ++b)
        {
            int Ntot = Nrad * Ncell;

            for (int r = 0; r < Nrad; ++r)
            {
                rf1[r] = 3.14;
            }

            // if (mv[b] < 0) continue;

            // int m = mv[b];

            // first attempt: propagate and store rad field
            // double ds = NEXT(b, );

            // peel-off emission
            //// propagate
            //// peel-off scatter
            //// scatter
        }
    }

    // debug test
    int Ntot = Nrad * Ncell;
    for (int r = 0; r < Nrad; ++r)
    {
        std::cout << r << " " << rf1[r] << std::endl;
    }
}
