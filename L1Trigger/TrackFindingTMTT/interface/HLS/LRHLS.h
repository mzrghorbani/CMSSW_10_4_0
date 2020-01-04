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
    void calcHelix();
    void calcResidual();
    void killResidual();
    void create();

public:

    Track *trackIn_;
    Track *trackOut_;
    LRTrack LRParameter_;
    residData residual_;
    array_t<LRStub> stubs_;
    array_t<residData> residuals_;
    int layerPopulation_[7];
    stubData layerPos_[7];
    int nStubs_;
    uint3_t nLayers_;
    uint3_t nLayersPS_;
    uint1_t valid_;

};

}

#endif