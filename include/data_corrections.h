//#include <iostream>
// ROOT includes
#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TString.h>
#include <TList.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TMinuit.h>
#include <TLorentzVector.h>
#include <TProfile.h>

#include "TClonesArray.h"

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"

#include "UHH2/bacondataformats/interface/TEventInfo.hh"
//#include "UHH2/bacondataformats/interface/BaconAnaDefs.hh"
using namespace std;

namespace uhh2bacon {

class DataCorr {

    private:
    uhh2::Context& context;
    uhh2::Event* event;



    public:
    DataCorr(uhh2::Context & ctx);
    ~DataCorr();
    // histogram name and pointer to it for Pu reweighting

    TFile *             file;
    float j2L1corr;
    float j3L1corr;

    void        SetEvent(uhh2::Event& evt);
     float       jet12getL1correction(float Eta1, float Eta2, float & j2L1corr);
     float       jet2getL1correction(float Eta);
     float       jet3getL1correction(float Eta);
//     float  getL1correction(float j1Eta, float j2Eta, float j3Eta, float & j2L1corr, float & j3L1corr);

};

}
