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

        lambdav.resize(_b, 0);
        weightv.resize(_b, 0);

        iv.resize(_b, 0);
        jv.resize(_b, 0);
        kv.resize(_b, 0);
        mv.resize(_b, 0);

        // packet
        rxv.resize(_b, 0);
        ryv.resize(_b, 0);
        rzv.resize(_b, 0);

        kxv.resize(_b, 0);
        kyv.resize(_b, 0);
        kzv.resize(_b, 0);

		sv.resize(_b, 0);
		tauinteractv.resize(_b, 0);
    }

    /** This function launches a photon packet with the specified wavelength, luminosity, and
        position and direction. The batch index is the index of the photon packet in this batch. */
    void launch(size_t batchIndex, double lambda, double L, Position bfr, Direction bfk);

    size_t batchSize() const { return _b; }

public:
    size_t _b{0};  // batch size

    vector<double> lambdav;
    vector<double> weightv;

    vector<int> iv;
    vector<int> jv;
    vector<int> kv;
    vector<int> mv;

    // packet
    vector<double> rxv;
    vector<double> ryv;
    vector<double> rzv;

    vector<double> kxv;
    vector<double> kyv;
    vector<double> kzv;

	vector<double> sv;
	vector<double> tauinteractv;
};

#define MAP_PACKETS(Nb, p) \
    map(to : p##lambdav [0:Nb], p##weightv [0:Nb]) map(to : p##iv [0:Nb], p##jv [0:Nb], p##kv [0:Nb], p##mv [0:Nb]) \
        map(to : p##rxv [0:Nb], p##ryv [0:Nb], p##rzv [0:Nb]) map(to : p##kxv [0:Nb], p##kyv [0:Nb], p##kzv [0:Nb]) \
			map(to : p##sv [0:Nb], p##tauinteractv [0:Nb])

#endif
