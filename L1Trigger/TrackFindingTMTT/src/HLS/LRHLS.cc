/*
Created by Maziar Ghorbani - Brunel University on 12/06/19.
*/

#include "L1Trigger/TrackFindingTMTT/interface/HLS/LRHLS.h"

namespace TMTT {


LRHLS::LRHLS(Track *trackIn, Track *trackOut)
        : trackIn_(trackIn), trackOut_(trackOut), nStubs_(0), nLayers_(0), nLayersPS_(0), valid_(0) {

    LRStub stub;

    for(auto i : trackIn_->stubs()) {
        stub.r = i->r();
        stub.phi = i->phi();
        stub.z = i->z();
        stub.layerId = i->layerId();
        stub.psModule = i->psModule();
        stub.barrel = i->barrel();
        stub.valid = 1;
        stubs_.push_back(stub);
    }

    nStubs_ = stubs_.size();
}

void LRHLS::produce() {

    initFit();

    valid_ = checkValidity();

    if(not valid_) {
        trackIn_->valid_ = false;
    }

    trackOut_ = trackIn_;

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

    if(nLayers_ < 6)
        valid = 0;
    if(nLayersPS_ < 4)
        valid = 0;

    return valid;
}


}