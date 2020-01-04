/*
Created by Maziar Ghorbani - Brunel University on 12/06/19.
*/

#include "L1Trigger/TrackFindingTMTT/interface/HLS/LRHLS.h"

namespace TMTT {


LRHLS::LRHLS(Track *trackIn, Track *trackOut)
        : trackIn_(trackIn), trackOut_(trackOut), nStubs_(0), nLayers_(0), nLayersPS_(0), valid_(0) {

    uint4_t i;

    trackOut_ = trackIn_;

    LRStub stub;

    for(auto s : trackIn_->stubs()) {
        stub.r = s->r();
        stub.phi = s->phi();
        stub.z = s->z();
        stub.layerId = s->layerId();
        stub.psModule = s->psModule();
        stub.barrel = s->barrel();
        stub.valid = 1;
        stubs_.push_back(stub);
    }

    nStubs_ = int(stubs_.size());
}

void LRHLS::produce() {

    uint4_t i, j;

    initFit();

    valid_ = checkValidity();

    if(valid_ == 0) {
        return;
    }

    for(i = 0; i < 7; i++) {
        if ((nStubs_ > 4) || (layerPopulation_[i] > 1)) {

            calcHelix();
            calcResidual();
            killResidual();
        }
    }

    for(i = 0; i < 7; i++) {
        cout << "  " << layerPopulation_[i];
    }
    cout << "        " << nStubs_;
    cout << endl;

    create();
}

void LRHLS::initFit() {

    uint4_t i;

    uint3_t layers[7];
    uint3_t layersPS[7];

    for(i = 0; i < 7; ++i) {
        layers[i] = 0;
        layersPS[i] = 0;
        layerPopulation_[i] = 0;
    }

    for(auto stub : stubs_) {

        if (stub.valid) {
            switch (stub.layerId) {
                case 1:
                    layerPopulation_[1]++;
                    layers[1]++;
                    if(stub.psModule)
                        layersPS[1]++;
                    break;
                case 2:
                    layerPopulation_[2]++;
                    layers[2]++;
                    if(stub.psModule)
                        layersPS[2]++;
                    break;
                case 3:
                    layerPopulation_[3]++;
                    layers[3]++;
                    if(stub.psModule)
                        layersPS[3]++;
                    break;
                case 4:
                    layerPopulation_[4]++;
                    layers[4]++;
                    if(stub.psModule)
                        layersPS[4]++;
                    break;
                case 5:
                    layerPopulation_[5]++;
                    layers[5]++;
                    if(stub.psModule)
                        layersPS[5]++;
                    break;
                case 6:
                    layerPopulation_[6]++;
                    layers[6]++;
                    if(stub.psModule)
                        layersPS[6]++;
                    break;
                default:
                    layerPopulation_[0]++;
                    layers[0]++;
                    if(stub.psModule)
                        layersPS[0]++;
            }
        }
    }

    for(i = 0; i < 7; ++i) {
        if(layers[i])
            nLayers_ += 1;
        if(layersPS[i])
            nLayersPS_ += 1;
    }
}

uint1_t LRHLS::checkValidity() {

    uint1_t valid = 1;

    if(nStubs_ < 3)
        valid = 0;
    if(nLayers_ < 4)
        valid = 0;
    if(nLayersPS_ < 2)
        valid = 0;

    return valid;
}

void LRHLS::calcHelix() {

    uint4_t i;

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

    uint1_t ps[7];

    for(i = 0; i < 7; i++) {
        ps[i] = 0;
    }

    for(auto stub : stubs_) {

        if(stub.valid) {
            stubData pos(stub.r, stub.phi, stub.r, stub.z);

            switch (stub.layerId) {
                case 1:
                    if(stub.psModule) {
                        ps[1] = 1;
                        layerMinPos[1] <= pos;
                        layerMaxPos[1] >= pos;
                    } else {
                        layerMinPos[1].RPhi = min_t(layerMinPos[1].RPhi, pos.RPhi);
                        layerMinPos[1].RPhi = min_t(layerMinPos[1].Phi, pos.Phi);
                        layerMaxPos[1].RPhi = max_t(layerMaxPos[1].RPhi, pos.RPhi);
                        layerMaxPos[1].RPhi = max_t(layerMaxPos[1].Phi, pos.Phi);
                    }
                    layerPos_[1] = layerMinPos[1] + layerMaxPos[1];
                    layerPos_[1] /= 2;
                    break;

                case 2:
                    if(stub.psModule) {
                        ps[2] = 1;
                        layerMinPos[2] <= pos;
                        layerMaxPos[2] >= pos;
                    } else {
                        layerMinPos[2].RPhi = min_t(layerMinPos[2].RPhi, pos.RPhi);
                        layerMinPos[2].RPhi = min_t(layerMinPos[2].Phi, pos.Phi);
                        layerMaxPos[2].RPhi = max_t(layerMaxPos[2].RPhi, pos.RPhi);
                        layerMaxPos[2].RPhi = max_t(layerMaxPos[2].Phi, pos.Phi);
                    }
                    layerPos_[2] = layerMinPos[2] + layerMaxPos[2];
                    layerPos_[2] /= 2;
                    break;

                case 3:
                    if(stub.psModule) {
                        ps[3] = 1;
                        layerMinPos[3] <= pos;
                        layerMaxPos[3] >= pos;
                    } else {
                        layerMinPos[3].RPhi = min_t(layerMinPos[3].RPhi, pos.RPhi);
                        layerMinPos[3].RPhi = min_t(layerMinPos[3].Phi, pos.Phi);
                        layerMaxPos[3].RPhi = max_t(layerMaxPos[3].RPhi, pos.RPhi);
                        layerMaxPos[3].RPhi = max_t(layerMaxPos[3].Phi, pos.Phi);
                    }
                    layerPos_[3] = layerMinPos[3] + layerMaxPos[3];
                    layerPos_[3] /= 2;
                    break;

                case 4:
                    if(stub.psModule) {
                        ps[4] = 1;
                        layerMinPos[4] <= pos;
                        layerMaxPos[4] >= pos;
                    } else {
                        layerMinPos[4].RPhi = min_t(layerMinPos[4].RPhi, pos.RPhi);
                        layerMinPos[4].RPhi = min_t(layerMinPos[4].Phi, pos.Phi);
                        layerMaxPos[4].RPhi = max_t(layerMaxPos[4].RPhi, pos.RPhi);
                        layerMaxPos[4].RPhi = max_t(layerMaxPos[4].Phi, pos.Phi);
                    }
                    layerPos_[4] = layerMinPos[4] + layerMaxPos[4];
                    layerPos_[4] /= 2;
                    break;

                case 5:
                    if(stub.psModule) {
                        ps[5] = 1;
                        layerMinPos[5] <= pos;
                        layerMaxPos[5] >= pos;
                    } else {
                        layerMinPos[5].RPhi = min_t(layerMinPos[5].RPhi, pos.RPhi);
                        layerMinPos[5].RPhi = min_t(layerMinPos[5].Phi, pos.Phi);
                        layerMaxPos[5].RPhi = max_t(layerMaxPos[5].RPhi, pos.RPhi);
                        layerMaxPos[5].RPhi = max_t(layerMaxPos[5].Phi, pos.Phi);
                    }
                    layerPos_[5] = layerMinPos[5] + layerMaxPos[5];
                    layerPos_[5] /= 2;
                    break;

                case 6:
                    if(stub.psModule) {
                        ps[6] = 1;
                        layerMinPos[6] <= pos;
                        layerMaxPos[6] >= pos;
                    } else {
                        layerMinPos[6].RPhi = min_t(layerMinPos[6].RPhi, pos.RPhi);
                        layerMinPos[6].RPhi = min_t(layerMinPos[6].Phi, pos.Phi);
                        layerMaxPos[6].RPhi = max_t(layerMaxPos[6].RPhi, pos.RPhi);
                        layerMaxPos[6].RPhi = max_t(layerMaxPos[6].Phi, pos.Phi);
                    }
                    layerPos_[6] = layerMinPos[6] + layerMaxPos[6];
                    layerPos_[6] /= 2;
                    break;

                default:
                    if(stub.psModule) {
                        ps[0] = 1;
                        layerMinPos[0] <= pos;
                        layerMaxPos[0] >= pos;
                    } else {
                        layerMinPos[0].RPhi = min_t(layerMinPos[0].RPhi, pos.RPhi);
                        layerMinPos[0].RPhi = min_t(layerMinPos[0].Phi, pos.Phi);
                        layerMaxPos[0].RPhi = max_t(layerMaxPos[0].RPhi, pos.RPhi);
                        layerMaxPos[0].RPhi = max_t(layerMaxPos[0].Phi, pos.Phi);
                    }
                    layerPos_[0] = layerMinPos[0] + layerMaxPos[0];
                    layerPos_[0] /= 2;
                    break;
            }
        }
    }
    for(i = 0; i < 7; i++) {
        if(layerPopulation_[i] > 0) {
            phiSums += make_pair_t(layerPos_[i].RPhi, layerPos_[i].Phi);
            if(ps[i] == 1)
                zSums += make_pair_t(layerPos_[i].RZ, layerPos_[i].Z);
        }
    }

    const pair_t<dtf_t, dtf_t> &phiParameter = phiSums.calcLinearParameter();
    const pair_t<dtf_t, dtf_t> &zParameter = zSums.calcLinearParameter();

    LRParameter_ = LRTrack(phiParameter.first, phiParameter.second, zParameter.first, zParameter.second);

}

void LRHLS::calcResidual() {

    uint4_t i;
    uint4_t j;

    dtf_t phiPredicted[7];
    dtf_t zPredicted[7];

    for (i = 0; i < 7; ++i) {
        phiPredicted[i] = 0;
        zPredicted[i] = 0;
    }

    for (i = 0; i < 7; i++) {
        if (layerPopulation_[i] > 0) {
            phiPredicted[i] = dtf_t(LRParameter_.phiT + dtf_t(LRParameter_.qOverPt * layerPos_[i].Phi));
            zPredicted[i] = dtf_t(LRParameter_.zT + dtf_t(LRParameter_.cotTheta * layerPos_[i].Z));
        }
    }

    for(i = 0; i < residuals_.size(); i++) {
        residuals_[i].phi = 0;
        residuals_[i].z = 0;
        residuals_[i].layerId = 0;
        residuals_[i].stubId = 0;
        residuals_[i].valid = 0;
    }

    for(i = 0; i < stubs_.size(); i++) {

        if (stubs_[i].valid) {
            for (j = 0; j < 7; j++) {
                if (stubs_[i].layerId == j) {
                    residual_.phi = abs_t(stubs_[i].phi - phiPredicted[j]);
                    residual_.z = abs_t(stubs_[i].z - zPredicted[j]);
                    residual_.layerId = j;
                    residual_.stubId = i;
                    residual_.valid = stubs_[i].valid;
                    residuals_.push_back(residual_);
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

    for (auto residual : residuals_) {

        if(residual.valid) {

            residComb = residual.combined();
            largest = residual_.combined();

            if (residComb > largest) {
                residual_ = residual;
            }
        }
    }

}

void LRHLS::killResidual() {

    stubs_[residual_.stubId].valid = 0;
    layerPopulation_[residual_.layerId]--;
    nStubs_--;

}

void LRHLS::create() {

    uint4_t i;

    trackOut_->qOverPt_ = LRParameter_.qOverPt;
    trackOut_->phi_ = LRParameter_.phiT;
    trackOut_->cot_ = LRParameter_.cotTheta;
    trackOut_->z_ = LRParameter_.zT;

    for(i = 0; i < stubs_.size(); i++) {
        trackOut_->stubs_[i]->r_ = stubs_[i].r;
        trackOut_->stubs_[i]->phi_ = stubs_[i].phi;
        trackOut_->stubs_[i]->z_ = stubs_[i].z;
        trackOut_->stubs_[i]->valid_ = stubs_[i].valid;
        trackOut_->stubs_[i]->module_->layerId_ = stubs_[i].layerId;
        trackOut_->stubs_[i]->module_->psModule_ = stubs_[i].psModule;
        trackOut_->stubs_[i]->module_->barrel_ = stubs_[i].barrel;
    }

}


}