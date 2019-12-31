/*
Created by Maziar Ghorbani - Brunel University on 12/06/19.
*/

#ifndef __LRHLS_H__
#define __LRHLS_H__

#include "L1Trigger/TrackFindingTMTT/interface/HLS/LRHLS_top.h"

namespace TMTT {

class LRHLS {

public:

    LRHLS(Track* in, Track* out);

    ~LRHLS() {}

    void produce();

public:

    void initFit();
    uint1_t checkValidity();
    void calcHelix();
    void calcResiduals();
    void killLargestResidual();

    Track* trackIn_;
    Track* trackOut_;
    LRTrack LRParameter_;

    LRStub stubs_[12];
    LRStub stub_;
    residData residuals_[12];
    residData residual_;
    uint4_t layerPopulation_[7];
    stubData layerPos_[7];
    uint3_t nLayers_;
    uint4_t nStubs_;
    uint3_t nIterations_;
    uint1_t valid_;
};

}

#endif
