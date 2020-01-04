/*
Created by Maziar Ghorbani - Brunel University on 12/06/19.
*/

#include "L1Trigger/TrackFindingTMTT/interface/HLS/LRHLS_top.h"

namespace TMTT {

LRHLS_top::LRHLS_top(const Settings *settings, Data *data) : settings_(settings), data_(data) {

}

void LRHLS_top::produce() {

    uint32_t i, j;

    data_->tracksLR().clear();

    for(i = 0; i < data_->tracksMHT().size(); i++) {
        LRHLS lrhls(data_->tracksMHT_[i], data_->tracksLR_[i]);
        lrhls.produce();
        data_->tracksLR_.push_back(lrhls.trackOut_);
    }

//    for(auto track : data_->tracksLR()) {
//        for (auto stub : track->stubs()) {
//            cout << stub->valid() << "  ";
//        }
//        cout << endl;
//    }
}

}