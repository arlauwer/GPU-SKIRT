#include "SimulationKernel.hpp"
#include <omp.h>
#include <iterator>

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

SimulationKernel::SimulationKernel(const CartesianSpatialGrid* grid)
    : _Nx(grid->_Nx), _Ny(grid->_Ny), _Nz(grid->_Nz),    //
      _gxv(std::begin(grid->_xv), std::end(grid->_xv)),  //
      _gyv(std::begin(grid->_yv), std::end(grid->_yv)),  //
      _gzv(std::begin(grid->_zv), std::end(grid->_zv))
{}

void SimulationKernel::runBatch()
{
    // grid
    double* gxv = _gxv.data();
    double* gyv = _gyv.data();
    double* gzv = _gzv.data();
    int Nx = _Nx;
    int Ny = _Ny;
    int Nz = _Nz;

    // photon packets
    size_t Nb = _photons.batchSize();
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

#pragma omp target data map(to : gxv[0 : Nx + 1], gyv[0 : Ny + 1], gzv[0 : Nz + 1]) \
    map(to : iv[0 : Nb], jv[0 : Nb], kv[0 : Nb], mv[0 : Nb]) map(to : lambdav[0 : Nb], weightv[0 : Nb]) \
    map(to : kxv[0 : Nb], kyv[0 : Nb], kzv[0 : Nb]) map(to : rxv[0 : Nb], ryv[0 : Nb], rzv[0 : Nb])
    {
#pragma omp target teams distribute parallel for firstprivate(Nx, Ny, Nz, Nb)
        for (size_t b = 0; b != Nb; ++b)
        {
            double ds = cartesianNext(gxv, gyv, gzv, Nx, Ny, Nz, iv[b], jv[b], kv[b], mv[b], rxv[b], ryv[b], rzv[b],
                                      kxv[b], kyv[b], kzv[b]);

            if (mv[b] < 0) continue;

            // launch peel-off
            //// propagate
            //// peel-off scatter
            //// scatter
        }
    }
}
