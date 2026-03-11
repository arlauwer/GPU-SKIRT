#include "SimulationKernel.hpp"
#include "CartesianSpatialGrid.hpp"
#include "MediumSystem.hpp"
#include "NR.hpp"
#include <omp.h>
#include <cmath>

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

        if (dsx <= dsy && dsx <= dsz)
        {
            rx = xE;
            ry += ky * dsx;
            rz += kz * dsx;
            i += xdir ? -1 : 1;
            m += xdir * Nz * Ny;
            return dsx;
        }
        else if (dsy < dsx && dsy <= dsz)
        {
            ry = yE;
            rx += kx * dsy;
            rz += kz * dsy;
            j += ydir ? -1 : 1;
            m += ydir * Nz;
            return dsy;
        }
        else
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

// Flat index into radiation field buffer: rf1[m * Nrad + ell]
inline size_t radIndex(size_t m, size_t ell, size_t Nrad)
{
    return m * Nrad + ell;
}

SimulationKernel::SimulationKernel(SourceSystem* ss, MediumSystem* ms) : _ss(ss), _ms(ms)
{
    // --- grid ---
    const auto grid = dynamic_cast<CartesianSpatialGrid*>(_ms->grid());
    _Nx = grid->_Nx;
    _Ny = grid->_Ny;
    _Nz = grid->_Nz;
    _gxv = grid->_xv;
    _gyv = grid->_yv;
    _gzv = grid->_zv;

    // --- number densities [Ncell] ---
    _Ncell = _ms->numCells();
    _nv = new double[_Ncell];
    for (size_t m = 0; m < _Ncell; ++m) _nv[m] = _ms->numberDensity(m, 0);

    // --- radiation field buffer [Ncell * Nrad] ---
    _Nrad = _ms->_wavelengthGrid->numBins();
    _rf1 = std::begin(_ms->_rf1.data());

    // --- cross section lookup grid ---
    // Build a log-spaced wavelength grid and sample the extinction cross section.
    // On the GPU, the bin index for a given wavelength lambda is:
    //   ell = (log(lambda) - _xsec_logLambdaMin) * _xsec_logInvBinWidth
    _Nxsec = 1000;
    double lambdaMin = _ss->minWavelength();
    double lambdaMax = _ss->maxWavelength();
    _xsec_logLambdaMin = std::log(lambdaMin);
    _xsec_logInvBinWidth = (_Nxsec - 1) / std::log(lambdaMax / lambdaMin);

    Array xsec_lambdaa;
    NR::buildLogGrid(xsec_lambdaa, lambdaMin, lambdaMax, _Nxsec);

    _crossv = new double[_Nxsec];
    for (int ell = 0; ell < _Nxsec; ++ell) _crossv[ell] = _ms->mix(0, 0)->sectionExt(xsec_lambdaa[ell]);
}

SimulationKernel::~SimulationKernel()
{
    delete[] _nv;
    delete[] _crossv;
}

void SimulationKernel::runBatch()
{
    // --- grid ---
    int Nx = _Nx;
    int Ny = _Ny;
    int Nz = _Nz;
    double* gxv = std::begin(_gxv);
    double* gyv = std::begin(_gyv);
    double* gzv = std::begin(_gzv);

    // --- number densities ---
    size_t Ncell = _Ncell;
    double* nv = _nv;

    // --- radiation field ---
    size_t Nrad = _Nrad;
    double* rf1 = _rf1;

    // --- cross section lookup ---
    int Nxsec = _Nxsec;
    double xsec_logLambdaMin = _xsec_logLambdaMin;
    double xsec_logInvBinWidth = _xsec_logInvBinWidth;
    double* crossv = _crossv;

    // --- photon packets ---
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

    // --- peel-off packets ---
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

// Advance photon packet b (with optional prefix p_ for peel-off) one cell step.
#define NEXT(b, p) \
    cartesianNext(gxv, gyv, gzv, Nx, Ny, Nz, p##iv[b], p##jv[b], p##kv[b], p##mv[b], p##rxv[b], p##ryv[b], p##rzv[b], \
                  p##kxv[b], p##kyv[b], p##kzv[b])

// Compute cross section bin index from wavelength (log-space lookup, O(1)).
#define XSEC_INDEX(lambda) static_cast<int>((std::log(lambda) - xsec_logLambdaMin) * xsec_logInvBinWidth)

#pragma omp target data map(to : gxv[0 : Nx + 1], gyv[0 : Ny + 1], gzv[0 : Nz + 1]) map(to : nv[0 : Ncell]) \
    map(to : crossv[0 : Nxsec]) MAP_PACKETS(Nb, ) MAP_PACKETS(Nb, p_) map(tofrom : rf1[0 : Ncell * Nrad])
    {
#pragma omp target teams distribute parallel for firstprivate(Nx, Ny, Nz, Nb, Ncell, Nrad, Nxsec, xsec_logLambdaMin, \
                                                                  xsec_logInvBinWidth)
        for (size_t b = 0; b != Nb; ++b)
        {
            if (mv[b] < 0) continue;

            // TODO: ray-trace packet b through grid, accumulating optical depth
            // Example usage:
            //   double ds = NEXT(b, );
            //   int ell = XSEC_INDEX(lambdav[b]);
            //   double dtau = nv[mv[b]] * crossv[ell] * ds;
            //   rf1[radIndex(mv[b], ell, Nrad)] += weightv[b] * ds;  // needs atomic
        }
    }
}
