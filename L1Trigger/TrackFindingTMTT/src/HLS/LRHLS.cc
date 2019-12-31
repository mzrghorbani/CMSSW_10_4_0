/*
Created by Maziar Ghorbani - Brunel University on 12/06/19.
*/

#include "L1Trigger/TrackFindingTMTT/interface/HLS/LRHLS.h"

namespace TMTT {

LRHLS::LRHLS(Track* in, Track* out) : trackIn_(in), trackOut_(out) {

    uint4_t i;

    for(auto stub : trackIn_->stubs()) {

        stubs_[i].r = stub->r();
        stubs_[i].phi = stub->phi();
        stubs_[i].z = stub->z();
        stubs_[i].layerId = stub->layerId();
        stubs_[i].valid = 1;
        i++;
    }

    for(i = 0; i < 7; ++i) {
        layerPopulation_[i] = 0;
    }

    nLayers_ = 0;
    nStubs_ = 0;
    valid_ = 0;
}

void LRHLS::produce() {

    initFit();

    valid_ = checkValidity();

    if(not valid_) {
        trackOut_ = trackIn_;
        return;;
    }

    calcHelix();
    calcResiduals();
    calcResiduals();
    killLargestResidual();
}

void LRHLS::initFit() {

    uint4_t i;

    uint3_t layers[7] = {0,0,0,0,0,0,0};

    for(i = 0; i < 12; ++i) {
        stub_ = stubs_[i];

        if (stub_.layerId) {
            switch (stub_.layerId) {
                case 1:
                    layerPopulation_[1]++;
                    layers[1] = 1;
                    nStubs_++;
                    break;
                case 2:
                    layerPopulation_[2]++;
                    layers[2] = 1;
                    nStubs_++;
                    break;
                case 3:
                    layerPopulation_[3]++;
                    layers[3] = 1;
                    nStubs_++;
                    break;
                case 4:
                    layerPopulation_[4]++;
                    layers[4] = 1;
                    nStubs_++;
                    break;
                case 5:
                    layerPopulation_[5]++;
                    layers[5] = 1;
                    nStubs_++;
                    break;
                case 6:
                    layerPopulation_[6]++;
                    layers[6] = 1;
                    nStubs_++;
                    break;
                default:
                    layerPopulation_[0]++;
                    layers[0] = 1;
                    nStubs_++;
            }
        }
    }

    for(i = 0; i < 7; ++i) {
        nLayers_ += layers[i];
    }
}

uint1_t LRHLS::checkValidity() {

    uint1_t valid = 0;

    if(nLayers_ < 4)
        valid = 0;
    if(nStubs_ < 4)
        valid = 0;

    return valid;
}

void LRHLS::calcHelix() {

    uint4_t i = 0;
    uint4_t j = 0;

    sumData phiSums;
    sumData zSums;

    stubData layerMinPos[7];

    for (i = 0; i < 7; i++) {
        layerMinPos[i].RPhi = 4095;
        layerMinPos[i].Phi = 4095;
        layerMinPos[i].RZ = 4095;
        layerMinPos[i].Z = 4095;
    }

    stubData layerMaxPos[7];

    for (i = 0; i < 7; i++) {
        layerMaxPos[i].RPhi = -4096;
        layerMaxPos[i].Phi = -4096;
        layerMaxPos[i].RZ = -4096;
        layerMaxPos[i].Z = -4096;
    }

    for(i = 0; i < 12; ++i) {
        stub_ = stubs_[i];

        stubData pos(stub_.r, stub_.phi, stub_.r, stub_.z);

        switch (stub_.layerId) {
            case 0:
                layerMinPos[0] <= pos;
                layerMaxPos[0] >= pos;
                layerPos_[0] = layerMinPos[0] + layerMaxPos[0];
                layerPos_[0] /= 2;
                break;

            case 1:
                layerMinPos[1] <= pos;
                layerMaxPos[1] >= pos;
                layerPos_[1] = layerMinPos[1] + layerMaxPos[1];
                layerPos_[1] /= 2;
                break;

            case 2:
                layerMinPos[2] <= pos;
                layerMaxPos[2] >= pos;
                layerPos_[2] = layerMinPos[2] + layerMaxPos[2];
                layerPos_[2] /= 2;
                break;

            case 3:
                layerMinPos[3] <= pos;
                layerMaxPos[3] >= pos;
                layerPos_[3] = layerMinPos[3] + layerMaxPos[3];
                layerPos_[3] /= 2;
                break;

            case 4:
                layerMinPos[4] <= pos;
                layerMaxPos[4] >= pos;
                layerPos_[4] = layerMinPos[4] + layerMaxPos[4];
                layerPos_[4] /= 2;
                break;

            case 5:
                layerMinPos[5] <= pos;
                layerMaxPos[5] >= pos;
                layerPos_[5] = layerMinPos[5] + layerMaxPos[5];
                layerPos_[5] /= 2;
                break;

            case 6:
                layerMinPos[6] <= pos;
                layerMaxPos[6] >= pos;
                layerPos_[6] = layerMinPos[6] + layerMaxPos[6];
                layerPos_[6] /= 2;
                break;
        }
    }

    for(i = 0; i < 7; i++) {
        if(layerPopulation_[i] > 0) {
            phiSums += make_pair_t(layerPos_[i].RPhi, layerPos_[i].Phi);
            zSums += make_pair_t(layerPos_[i].RZ, layerPos_[i].Z);
        }
    }

    const pair_t<dtf_t, dtf_t> &phiParameter = phiSums.calcLinearParameter();
    const pair_t<dtf_t, dtf_t> &zParameter = zSums.calcLinearParameter();

    LRParameter_ = LRTrack(phiParameter.first, phiParameter.second, zParameter.first, zParameter.second);
}

void LRHLS::calcResiduals() {

    uint4_t i = 0;
    uint4_t j = 0;

    dtf_t phiPredicted[7] = {0, 0, 0, 0, 0, 0, 0};

    dtf_t zPredicted[7] = {0, 0, 0, 0, 0, 0, 0};

    for (i = 0; i < 7; i++) {
        if (layerPopulation_[i] > 0) {
            phiPredicted[i] = dtf_t(LRParameter_.phiT + dtf_t(LRParameter_.qOverPt * layerPos_[i].Phi));
            zPredicted[i] = dtf_t(LRParameter_.zT + dtf_t(LRParameter_.cotTheta * layerPos_[i].Z));
        }
    }

    for(i = 0; i < 10; i++) {
        residuals_[i].phi = 0;
        residuals_[i].z = 0;
        residuals_[i].layerId = 0;
        residuals_[i].stubId = 0;
        residuals_[i].valid = 0;
    }

    for (i = 0; i < 10; i++) {
        stub_ = stubs_[i];

        if (stub_.valid) {
            for (j = 0; j < 7; j++) {
                if (stub_.layerId == j) {
                    residuals_[i].phi = abs_t(stub_.phi - phiPredicted[j]);
                    residuals_[i].z = abs_t(stub_.z - zPredicted[j]);
                    residuals_[i].layerId = j;
                    residuals_[i].stubId = i;
                    residuals_[i].valid = 1;
                    break;
                }
            }
        }
    }

    residual_.phi = 0;
    residual_.z = 0;
    residual_.layerId = 0;
    residual_.stubId = 0;

    dtf_t residComb = 0;
    dtf_t largest = 0;

    for (i = 0; i < 12; i++) {

        if(residuals_[i].valid) {

            residComb = residuals_[i].combined();
            largest = residual_.combined();

            if (residComb > largest) {
                residual_ = residuals_[i];
            }
        }
    }
}

void LRHLS::killLargestResidual() {

#ifndef __SYNTHESIS__
    std::cout << "residual_.phi: " << residual_.phi << std::endl;
    std::cout << "residual_.z: " << residual_.z << std::endl;
    std::cout << "residual_.layerId: " << residual_.layerId << std::endl;
    std::cout << "residual_.stubId: " << residual_.stubId << std::endl;
#endif

    stubs_[residual_.stubId].valid = 0;
    layerPopulation_[residual_.layerId]--;
    nStubs_--;
}

}
