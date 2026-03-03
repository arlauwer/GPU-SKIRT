#include "PhotonPackets.hpp"

void PhotonPackets::launch(size_t batchIndex, double lambda, double L, Position bfr, Direction bfk)
{
    _lambdav[batchIndex] = lambda;
    _weightv[batchIndex] = L;
    _bfrv[batchIndex] = bfr;
    _bfkv[batchIndex] = bfk;
}
