/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PHOTONPACKETS_HPP
#define PHOTONPACKETS_HPP

#include "Basics.hpp"
#include "Position.hpp"
#include <cstddef>

class PhotonPackets
{
public:
    void setBatchSize(size_t batchSize)
    {
        _b = batchSize;
        _lambdav.resize(_b, 0);
        _weightv.resize(_b, 0);
        _bfrv.resize(_b, Position());
        _bfkv.resize(_b, Direction());
    }

    /** This function launches a photon packet with the specified wavelength, luminosity, and
        position and direction. The batch index is the index of the photon packet in this batch. */
    void launch(size_t batchIndex, double lambda, double L, Position bfr, Direction bfk);

public:            // debug public
    size_t _b{0};  // batch size
    vector<double> _lambdav;
    vector<double> _weightv;
    vector<Position> _bfrv;
    vector<Direction> _bfkv;
};

#endif
