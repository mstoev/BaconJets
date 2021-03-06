// system headers
#include<iostream>
#include<fstream>
#include <sys/stat.h>

// ROOT headers
#include<TCollection.h>
#include<TMath.h>
#include<TFile.h>
#include<TRandom3.h>
#include<TSystem.h>
#include<TVector3.h>
#include <vector>

#include "UHH2/BaconJets/include/mc_weight.h"

#include "UHH2/bacondataformats/interface/TJet.hh"
using namespace std;

namespace uhh2bacon {


McWeight::McWeight(uhh2::Context & ctx) :
    context(ctx),
    event(0),
    fPuReweighting_histo(NULL)
{
    h_jets = context.declare_event_input<TClonesArray>("Jet05");
    h_eventInfo = context.declare_event_input<baconhep::TEventInfo>("Info");

    fPuReweighting_histoname.push_back(hPuReweighting_histo40);
    fPuReweighting_histoname.push_back(hPuReweighting_histo80);
    fPuReweighting_histoname.push_back(hPuReweighting_histo140);
    fPuReweighting_histoname.push_back(hPuReweighting_histo200);
    fPuReweighting_histoname.push_back(hPuReweighting_histo260);
    fPuReweighting_histoname.push_back(hPuReweighting_histo320);
    fPuReweighting_histoname.push_back(hPuReweighting_histo400);

    TString DATABASE_PATH = "/nfs/dust/cms/user/kovalch/DataPileup/PuWeights";
    TString fPuReweighting_filename[] = {"nPu_reweighting_40.root", "nPu_reweighting_80.root", "nPu_reweighting_140.root", "nPu_reweighting_200.root", "nPu_reweighting_260.root", "nPu_reweighting_320.root", "nPu_reweighting_400.root"};

    for (int j = 0; j < fPuReweighting_histoname.size(); j++) {
        file = new TFile (DATABASE_PATH+"/"+fPuReweighting_filename[j]);
        fPuReweighting_histoname[j] = (TH1F*) file -> Get("histo_substr");

        if (fPuReweighting_histoname[j] == NULL) {
            cout << "was not possible to retrieve a histogram from "<< fPuReweighting_filename[j]<< endl;
            abort();
        } else {
            cout << "successfully got the PU reweighting histogram from "<< fPuReweighting_filename[j]<< endl;
        }

    }
}

void McWeight::SetEvent(uhh2::Event& evt)
{
    event = &evt;
    assert(event);

}

//gets a weighting factor for true-PU reweighting 
float  McWeight::getPuReweighting() {

    assert(event);
    const baconhep::TEventInfo & info = event->get(h_eventInfo);
    baconhep::TEventInfo* eventInfo= new baconhep::TEventInfo(info);
    assert(eventInfo);

    Float_t   Pu_true = eventInfo->nPUmean;
    // to check which triggers were fired
//     for (int j = 0; j < 7; j++) {
//         cout << " trig pass "<< eventInfo->triggerBits[j] << endl;
//     }

    // declare variable to store a weighting factor
    Double_t    weighting_factor = 1;
    unsigned    bin = 0;

    // sanity check: a pointer to the histogram should be valid
    for (int i = 0; i < fPuReweighting_histoname.size(); i++) {
        if (fPuReweighting_histoname[i] == NULL) {
            cout << "was not possible to retrieve a histogram " << endl;
            abort();
        }
    }

    // weights used for PU reweighting is normalized
    // using the weight corresponds to the highest fired trigger
    if (info.triggerBits[6]==1){ //PFJetAve400
        bin = fPuReweighting_histoname[6] -> FindBin(Pu_true);
        weighting_factor = fPuReweighting_histoname[6] -> GetBinContent(bin);
//        cout << "PU distribution for trigger PFJetAve400 was used"<<endl;
    } else if (info.triggerBits[5]==1) { //PFJetAve320
        bin = fPuReweighting_histoname[5] -> FindBin(Pu_true);
        weighting_factor = fPuReweighting_histoname[5] -> GetBinContent(bin);
//         cout << "PU distribution for trigger PFJetAve320 was used"<<endl;
    } else if (info.triggerBits[4]==1) { //PFJetAve260
        bin = fPuReweighting_histoname[4] -> FindBin(Pu_true);
        weighting_factor = fPuReweighting_histoname[4] -> GetBinContent(bin);
//         cout << "PU distribution for trigger PFJetAve260 was used"<<endl;
    } else if (info.triggerBits[3]==1) { //PFJetAve200
        bin = fPuReweighting_histoname[3] -> FindBin(Pu_true);
        weighting_factor = fPuReweighting_histoname[3] -> GetBinContent(bin);
//         cout << "PU distribution for trigger PFJetAve200 was used"<<endl;
    } else if (info.triggerBits[2]==1) { //PFJetAve140
        bin = fPuReweighting_histoname[2] -> FindBin(Pu_true);
        weighting_factor = fPuReweighting_histoname[2] -> GetBinContent(bin);
//         cout << "PU distribution for trigger PFJetAve140 was used"<<endl;
    } else if (info.triggerBits[1]==1) {//PFJetAve80
        bin = fPuReweighting_histoname[1] -> FindBin(Pu_true);
        weighting_factor = fPuReweighting_histoname[1] -> GetBinContent(bin);
//         cout << "PU distribution for trigger PFJetAve80 was used"<<endl;
    } else if (info.triggerBits[0]==1) {//PFJetAve40
        bin = fPuReweighting_histoname[0] -> FindBin(Pu_true);
        weighting_factor = fPuReweighting_histoname[0] -> GetBinContent(bin);
//         cout << "PU distribution for trigger PFJetAve40 was used"<<endl;
    }

    // finished successfully, return the weighting factor
//       cout <<" ev nu: "<<eventInfo->evtNum << " nPU = "<< Pu_true<< " weighting_factor = "<< weighting_factor<< " bin # = "<< bin <<endl;

    if (weighting_factor!=0) return      weighting_factor;
    else return 1;
}


McWeight::~McWeight()
{
}

}