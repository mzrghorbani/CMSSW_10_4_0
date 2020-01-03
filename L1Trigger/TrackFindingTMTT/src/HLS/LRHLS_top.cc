/*
Created by Maziar Ghorbani - Brunel University on 12/06/19.
*/

#include "L1Trigger/TrackFindingTMTT/interface/HLS/LRHLS_top.h"

namespace TMTT {

LRHLS_top::LRHLS_top(const Settings *settings, Data *data) : settings_(settings), data_(data) {

}

void LRHLS_top::produce() {

    data_->tracksLR().clear();
    cout << "trackIn cleared : " << data_->tracksLR().size() << endl;

    cout << "trackIn size : " << data_->tracksMHT().size() << endl;


    for(unsigned int i = 0; i < data_->tracksMHT().size(); i++) {
        LRHLS lrhls(data_->tracksMHT_[i], data_->tracksLR_[i]);
        lrhls.produce();
        if(lrhls.trackOut_->valid())
            data_->tracksLR_.push_back(lrhls.trackOut_);
    }

    for(auto & i : data_->tracksLR()) {
        cout << i->valid() << " ";
    }

    cout << endl;

    cout << "trackOut size : " << data_->tracksLR().size() << endl;

//    TPs tps;
//    for (const Track *track : data_->tracksLR()) {
//        for (TP *tp : track->tps()) {
//            if (tp->useForAlgEff()) {
//                tps->push_back(tp);
//            }
//        }
//    }

}

}