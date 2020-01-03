/*
Created by Maziar Ghorbani - Brunel University on 12/06/19.
*/

#ifndef __LRHLS_H__
#define __LRHLS_H__

#include "L1Trigger/TrackFindingTMTT/interface/HLS/LRHLS_types.h"
#include "L1Trigger/TrackFindingTMTT/interface/HLS/LRHLS_top.h"

namespace TMTT {

class LRHLS {

public:

    LRHLS(Track *trackIn, Track *trackOut);

    ~LRHLS() {}

    void produce();

    void initFit();

    uint1_t checkValidity();

public:
    Track *trackIn_;
    Track *trackOut_;

    array_t<LRStub> stubs_;
    LRStub stub_;

    uint3_t layerPopulation_[7];

    uint4_t nStubs_;
    uint3_t nLayers_;
    uint3_t nLayersPS_;
    uint1_t valid_;

};

}

#endif