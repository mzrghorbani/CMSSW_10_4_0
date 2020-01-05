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
    bool checkValidity();
    void calcHelix();
    void calcResidual();
    bool killLargestResidual();
    void findLargestResidual();
    void create();

public:

    Track *trackIn_;
    Track *trackOut_;
    array_t<LRStub> stubs_;
    array_t<residData> residuals_;
    LRTrack HTParameter_;
    LRTrack LRParameter_;
    residData largestResid_;
    stubData layerPos_[7];
    uint4_t layerPopulation_[7]{};
    uint4_t nIterations_;
    uint4_t nLayers_;
    uint4_t nLayersPS_;
    uint4_t nStubs_;
    bool valid_;
};

}

#endif