#include "PhotonPackets.hpp"

void PhotonPackets::launch(size_t batchIndex, double lambda, double L, Position bfr, Direction bfk)
{
    lambdav[batchIndex] = lambda;
    weightv[batchIndex] = L;
    // does not set cell index
    rxv[batchIndex] = bfr.x();
    ryv[batchIndex] = bfr.y();
    rzv[batchIndex] = bfr.z();
    kxv[batchIndex] = bfk.x();
    kyv[batchIndex] = bfk.y();
    kzv[batchIndex] = bfk.z();
}
