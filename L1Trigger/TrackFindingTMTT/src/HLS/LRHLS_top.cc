/*
Created by Maziar Ghorbani - Brunel University on 12/06/19.
*/

#include "L1Trigger/TrackFindingTMTT/interface/HLS/LRHLS_top.h"

namespace TMTT {

LRHLS_top::LRHLS_top(const Settings* settings, Data *data) : settings_(settings), data_(data) {

}

void LRHLS_top::produce() {

    for (unsigned int i = 0; i < data_->tracksMHT().size(); ++i) {

        LRHLS lrhls(data_->tracksMHT_[i], data_->tracksLR_[i]);
        lrhls.produce();

    }

//    TPs tps;
//    for (const Track *track : tracks) {
//        for (TP *tp : track->tps()) {
//            if (tp.useForAlgEff()) {
//                tps.push_back(tp);
//            }
//        }
//    }

}

}