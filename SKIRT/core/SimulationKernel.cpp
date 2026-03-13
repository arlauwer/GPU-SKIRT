#include "SimulationKernel.hpp"
#include "CartesianSpatialGrid.hpp"
#include "MediumSystem.hpp"
#include "NR.hpp"
#include "PhotonPackets.hpp"
#include <omp.h>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstddef>

namespace
{
    inline double cartesianNext(const double* gxv, const double* gyv, const double* gzv, int Nx, int Ny, int Nz, int& i,
                                int& j, int& k, int& m, double& rx, double& ry, double& rz, double& kx, double& ky,
                                double& kz)
    {
        int xdir = std::copysign(1, kx);
        int ydir = std::copysign(1, ky);
        int zdir = std::copysign(1, kz);

        assert(i >= 0 && i < Nx + 1);
        double xE = (xdir > 0) ? gxv[i + 1] : gxv[i];
        double yE = (ydir > 0) ? gyv[j + 1] : gyv[j];
        double zE = (zdir > 0) ? gzv[k + 1] : gzv[k];
        double dsx = (fabs(kx) > 1e-15) ? (xE - rx) / kx : DBL_MAX;
        double dsy = (fabs(ky) > 1e-15) ? (yE - ry) / ky : DBL_MAX;
        double dsz = (fabs(kz) > 1e-15) ? (zE - rz) / kz : DBL_MAX;

        if (dsx <= dsy && dsx <= dsz)
        {
            rx = xE;
            ry += ky * dsx;
            rz += kz * dsx;
            i += xdir;
            m += xdir * Nz * Ny;
            return dsx;
        }
        else if (dsy < dsx && dsy <= dsz)
        {
            ry = yE;
            rx += kx * dsy;
            rz += kz * dsy;
            j += ydir;
            m += ydir * Nz;
            return dsy;
        }
        else
        {
            rz = zE;
            rx += kx * dsz;
            ry += ky * dsz;
            k += zdir;
            m += zdir;
            return dsz;
        }
    }

    inline bool inside(int Nx, int Ny, int Nz, int i, int j, int k)
    {
        return i >= 0 && i < Nx && j >= 0 && j < Ny && k >= 0 && k < Nz;
    }

    inline size_t logIndex(double llambda, double logMin, double logInvBinWidth)
    {
        return (llambda - logMin) * logInvBinWidth;
    }

    inline double lnmean(double x1, double x2, double lnx1, double lnx2)
    {
        if (x1 > x2)
        {
            std::swap(x1, x2);
            std::swap(lnx1, lnx2);
        }
        if (x1 <= 0) return 0.;

        double x = x2 / x1 - 1.;
        if (x < 1e-3)
        {
            return x1
                   / (1. - 1. / 2. * x + 1. / 3. * x * x - 1. / 4. * x * x * x + 1. / 5. * x * x * x * x
                      - 1. / 6. * x * x * x * x * x);
        }
        else
        {
            return (x2 - x1) / (lnx2 - lnx1);
        }
    }
}

inline size_t radIndex(size_t m, size_t ell, size_t Nrad)
{
    return m * Nrad + ell;
}

SimulationKernel::SimulationKernel(SourceSystem* ss, MediumSystem* ms) : _ss(ss), _ms(ms)
{
    _random = ss->find<Random>();
    const auto grid = dynamic_cast<CartesianSpatialGrid*>(_ms->grid());
    _Nx = grid->_Nx;
    _Ny = grid->_Ny;
    _Nz = grid->_Nz;
    _gxv = grid->_xv;
    _gyv = grid->_yv;
    _gzv = grid->_zv;

    _Ncell = _ms->numCells();
    _nv = new double[_Ncell];
    for (size_t m = 0; m < _Ncell; ++m) _nv[m] = _ms->numberDensity(m, 0);

    _rad = std::begin(_ms->_rf1.data());
    // radiation field wavelength grid log lookup
    _Nrad = _ms->_wavelengthGrid->numBins();
    double radLambdaMin = _ms->_wavelengthGrid->leftBorder(0);
    double radLambdaMax = _ms->_wavelengthGrid->rightBorder(_Nrad - 1);
    _rad_logLambdaMin = std::log10(radLambdaMin);
    _rad_logInvBinWidth = (_Nrad - 1) / std::log10(radLambdaMax / radLambdaMin);

    // cross section wavelength grid log lookup
    _Nsec = 1000;
    double lambdaMin = _ss->minWavelength();
    double lambdaMax = _ss->maxWavelength();
    _sec_logLambdaMin = std::log10(lambdaMin);
    _sec_logInvBinWidth = (_Nsec - 1) / std::log10(lambdaMax / lambdaMin);

    // cross section wavelength grid
    Array sec_lambdaa;
    NR::buildLogGrid(sec_lambdaa, lambdaMin, lambdaMax, _Nsec);

    // tabulated extinction cross sections
    _crossv = new double[_Nsec];
    for (size_t ell = 0; ell < _Nsec; ++ell) _crossv[ell] = _ms->mix(0, 0)->sectionExt(sec_lambdaa[ell]);
}

SimulationKernel::~SimulationKernel()
{
    delete[] _nv;
    delete[] _crossv;
}

void SimulationKernel::runBatch()
{

    // Peelof emission

    // Loop
	bool done = false;
	size_t Nb = _photons.batchSize();
	while (!done)
    {
        // Calculate tauinteractv
        for (size_t b = 0; b != Nb; ++b)
        {
            double tauinteract = -log(_random->uniform());
            _photons.tauinteractv[b] = tauinteract;
        }
        traverse(_photons);
        // Remove photons that left the system
		done = true;
        for (size_t b = 0; b < Nb; ++b)
        {
			if (inside(_Nx, _Ny, _Nz, _photons.iv[b], _photons.jv[b], _photons.kv[b]))
			{
				done = false;
			}
			else
			{
				_photons.mv[b] = -1;
			}
        }
        // Not used?
        vector<int> mv = _photons.mv;
        vector<double> sv = _photons.sv;

        // Scatter
		if (!done) _ms->simulateScattering(_random, _photons);
        // Peelof scatter
    }
}

void SimulationKernel::traverse(PhotonPackets& pp)
{
    // cartesian grid geometry
    int Nx = _Nx;
    int Ny = _Ny;
    int Nz = _Nz;
    double* gxv = std::begin(_gxv);
    double* gyv = std::begin(_gyv);
    double* gzv = std::begin(_gzv);

    // number densities per cell, indexed by cell index m [Ncell]
    size_t Ncell = _Ncell;
    double* nv = _nv;

    // tabulated radiation field, indexed by [m * Nrad + ell]
    size_t Nrad = _Nrad;
    double* rad = _rad;
    // log-spaced wavelength grid for radiation field lookups// log-spaced lookup for radiation field wavelength grid
    double rad_logLambdaMin = _rad_logLambdaMin;
    double rad_logInvBinWidth = _rad_logInvBinWidth;

    // tabulated extinction cross sections
    int Nsec = _Nsec;
    double* crossv = _crossv;
    // log-spaced wavelength grid for cross section lookups
    double sec_logLambdaMin = _sec_logLambdaMin;
    double sec_logInvBinWidth = _sec_logInvBinWidth;

    size_t Nb = pp.batchSize();

    // photon packets
    double* lambdav = pp.lambdav.data();
    double* weightv = pp.weightv.data();
    int* iv = pp.iv.data();
    int* jv = pp.jv.data();
    int* kv = pp.kv.data();
    int* mv = pp.mv.data();
    double* rxv = pp.rxv.data();
    double* ryv = pp.ryv.data();
    double* rzv = pp.rzv.data();
    double* kxv = pp.kxv.data();
    double* kyv = pp.kyv.data();
    double* kzv = pp.kzv.data();
    double* sv = pp.sv.data();
    double* tauinteractv = pp.tauinteractv.data();

#define NEXT(b, p) \
    cartesianNext(gxv, gyv, gzv, Nx, Ny, Nz, p##iv[b], p##jv[b], p##kv[b], p##mv[b], p##rxv[b], p##ryv[b], p##rzv[b], \
                  p##kxv[b], p##kyv[b], p##kzv[b])

#pragma omp target data map(to : gxv[0 : Nx + 1], gyv[0 : Ny + 1], gzv[0 : Nz + 1], nv[0 : Ncell], crossv[0 : Nsec]) \
    MAP_PACKETS(Nb, ) map(tofrom : rad[0 : Ncell * Nrad])
    {
#pragma omp target teams distribute parallel for firstprivate( \
        Nx, Ny, Nz, Nb, Ncell, Nrad, sec_logLambdaMin, sec_logInvBinWidth, rad_logLambdaMin, rad_logInvBinWidth)
        for (size_t b = 0; b != Nb; ++b)
        {
            int m = mv[b];
            if (m < 0) continue;

            double lambda = lambdav[b];
            double llambda = std::log10(lambda);
            double weight = weightv[b];
            size_t sec_l = logIndex(llambda, sec_logLambdaMin, sec_logInvBinWidth);
            size_t rad_l = logIndex(llambda, rad_logLambdaMin, rad_logInvBinWidth);

            // Not clamping sec or rad!

            double lnExtBeg = 0.0;
            double extBeg = 1.0;

            double tau = 0;

            // Traverse all cells until the packet exits the grid
            while (inside(Nx, Ny, Nz, iv[b], jv[b], kv[b]) && tau < tauinteractv[b])
            {
                m = mv[b];

                double ds = NEXT(b, );

                sv[b] += ds;

                // Optical depth contribution from this segment
                double kappa = nv[m] * crossv[sec_l];
                double dtau = kappa * ds;
				tau += dtau;
                double lnExtEnd = lnExtBeg - dtau;
                double extEnd = exp(lnExtEnd);

                // Logarithmic mean extinction over the segment (lnmean of extBeg, extEnd)
                double extMean = lnmean(extEnd, extBeg, lnExtEnd, lnExtBeg);

                double Lds = weight * extMean * ds;

// Atomic update of the radiation field
#pragma omp atomic
                rad[radIndex(m, rad_l, Nrad)] += Lds;

                lnExtBeg = lnExtEnd;
                extBeg = extEnd;
            }
        }
    }
}
