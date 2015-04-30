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

class McWeight {

    private:
    uhh2::Context& context;
    uhh2::Event* event;

    uhh2::Event::Handle<TClonesArray> h_jets;
    uhh2::Event::Handle<baconhep::TEventInfo> h_eventInfo;

    bool        fPuReweighting;

    public:
    McWeight(uhh2::Context & ctx);
    ~McWeight();
    // histogram name and pointer to it for Pu reweighting
    TH1F *              fPuReweighting_histo;
    std::vector<TH1F*>  fPuReweighting_histoname;
    TH1F *              hPuReweighting_histo40; 
    TH1F *              hPuReweighting_histo80;
    TH1F *              hPuReweighting_histo140;
    TH1F *              hPuReweighting_histo200;
    TH1F *              hPuReweighting_histo260;
    TH1F *              hPuReweighting_histo320;
    TH1F *              hPuReweighting_histo400;
    TFile *             file;

    void        SetEvent(uhh2::Event& evt);
    float       getPuReweighting();


};

}