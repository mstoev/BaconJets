#include "UHH2/BaconJets/include/JECAnalysisHists.h"
#include "UHH2/BaconJets/include/constants.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/core/include/Jet.h"

#include "UHH2/bacondataformats/interface/TJet.hh"
#include "UHH2/bacondataformats/interface/TEventInfo.hh"
#include "UHH2/bacondataformats/interface/BaconAnaDefs.hh"
#include "UHH2/bacondataformats/interface/TVertex.hh"



#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <getopt.h>
using namespace std;
using namespace uhh2;
using namespace baconhep;
    uhh2::Event::Handle<TClonesArray> h_pv;
JECAnalysisHists::JECAnalysisHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
    // book all histograms here
    // jets
    TH1::SetDefaultSumw2();

    book<TH1F>("N_jets", "N_{jets}", 20, -0.5, 19.5);
    book<TH1F>("pt","p_{T} all jets; p_{T} (GeV)",100,0,1500);
    book<TH1F>("eta","#eta all jets; #eta",100,-5,5);
    book<TH1F>("phi","#phi all jets; #phi",50,-M_PI,M_PI);
    book<TH1F>("MET","MET all jets; MET",400,0,400);

    book<TH1F>("nPu","Number of PU events",60,0,60);
    book<TH1F>("N_PV","Number of PVtx",60,0,60);
    book<TH1F>("N_PV_sel","Number of PVtx after cuts on z, y, x",60,0,60);

    book<TH1F>("weight_histo","weight_histo ",20,0,2);

    book<TH1F>("pt_1","p_{T} jet 1",100,0,600);
    book<TH1F>("eta_1","#eta jet 1",100,-5,5);

    book<TH1F>("pt_2","p_{T} jet 2",100,0,600);
    book<TH1F>("eta_2","#eta jet 2",100,-5,5);

    book<TH1F>("pt_3","p_{T} jet 3",100,0,600);
    book<TH1F>("eta_3","#eta jet 3",100,-5,5);

    book<TH1F>("ptRaw_barrel","p^{Raw}_{T} barrel jet; p_{T}^{Raw,barrel} (GeV)",100,0,600);
    book<TH1F>("ptRaw_probe","p^{Raw}_{T} probe jet; p_{T}^{Raw,probe} (GeV)",100,0,600);
    book<TH1F>("pt_barrel","p_{T} barrel jet; p_{T}^{barrel} (GeV)",100,0,600);
    book<TH1F>("pt_probe","p_{T} probe jet; p_{T}^{probe} (GeV)",100,0,600);
    book<TH1F>("eta_barrel","#eta barrel jet; #eta^{barrel}",100,-5,5);
    book<TH1F>("eta_probe","#eta probe jet #eta^{probe}",100,-5,5);

    book<TH1F>("pt_ave","p_{T} ave jet; p_{T}^{ave} (GeV)",100,0,600);

    book<TH1F>("pt_rel","p_{T}^{jet3} / p_{T}^{ave}; #alpha ", 50, 0, 1);
    book<TH1F>("generic_pt_rel","generic p_{T}^{jet3} / p_{T}^{ave}; #alpha ", 50, 0, 1);

    book<TH1F>("asym","asymmetrie jet 1 and jet 2 in loop; Asymmetry",100,-1,1);
    book<TH1F>("generic_asym","generic asymmetrie jet 1 and jet 2; Asymmetry",100,-1,1);
    book<TH1F>("mpf","MPF response; MPF response",100,0.5,1.5);
    book<TH1F>("generic_mpf","generic MPF response; MPF response",100,0.5,1.5);
    book<TH1F>("r_rel","R_{rel}; R_{rel}; Relative response",100,0.5,1.5);
    book<TH1F>("generic_r_rel","generic R_{rel}; Relative response",100,0.5,1.5);

    book<TH1F>("DeltaPhi_Jet1_Jet2", "#Delta#Phi(first jet, second jet); #Delta #Phi", 100, 0, 7);
    book<TH2F>("ptrel_vs_deltaphi","delta phi vs pt_rel", 50, 0, 1 ,50, 0, 3.14);

    book<TH1F>("pt_ave_hltDiPFJetAve40","p_{T} ave40 jet; p_{T}^{ave} (GeV)",120,0,600);
    book<TH1F>("pt_ave_hltDiPFJetAve60","p_{T} ave60 jet; p_{T}^{ave} (GeV)",120,0,600);
    book<TH1F>("pt_ave_hltDiPFJetAve80","p_{T} ave80 jet; p_{T}^{ave} (GeV)",120,0,600);
    book<TH1F>("pt_ave_hltDiPFJetAve140","p_{T} ave140 jet; p_{T}^{ave} (GeV)",120,0,600);
    book<TH1F>("pt_ave_hltDiPFJetAve200","p_{T} ave200 jet; p_{T}^{ave} (GeV)",120,0,600);
    book<TH1F>("pt_ave_hltDiPFJetAve260","p_{T} ave260 jet; p_{T}^{ave} (GeV)",120,0,600);
    book<TH1F>("pt_ave_hltDiPFJetAve320","p_{T} ave320 jet; p_{T}^{ave} (GeV)",120,0,600);
    book<TH1F>("pt_ave_hltDiPFJetAve400","p_{T} ave400 jet; p_{T}^{ave} (GeV)",120,0,600);
    book<TH1F>("pt_ave_hltDiPFJetAve500","p_{T} ave500 jet; p_{T}^{ave} (GeV)",120,0,600);

    book<TH2F>("Rrel_vs_Npv","Rrel vs. Npv ", 40, 0, 40 ,100, 0.1,4.5);
//     book<TProfile>("prof_Rrel_vs_Npv","Rrel vs. Npv ", 50, 0, 50 ,100, 0.5,1.5);
//     TProfile * prof_Rrel_vs_Npv  = new TProfile("prof_Rrel_vs_Npv","Rrel vs. Npv",100, 0, 50, 0.5, 1.5);


/*    uhh2::Event::Handle<TClonesArray> h_pv;*/
    h_jets = ctx.get_handle<TClonesArray>("AK4PFCHS");
    h_eventInfo = ctx.get_handle<baconhep::TEventInfo>("Info");
    h_pv = ctx.get_handle<TClonesArray>("PV");

}

void JECAnalysisHists::fill(const uhh2::Event & ev){
    fill(ev, 0);
}
void JECAnalysisHists::fill(const uhh2::Event & ev, const int rand){
    // fill the histograms. Please note the comments in the header file:
    // 'hist' is used here a lot for simplicity, but it will be rather
    // slow when you have many histograms; therefore, better
    // use histogram pointers as members as in 'UHH2/common/include/ElectronHists.h'
    // Don't forget to always use the weight when filling.

    const TClonesArray & js = ev.get(h_jets);
    const baconhep::TEventInfo & info = ev.get(h_eventInfo);
    baconhep::TEventInfo* eventInfo= new baconhep::TEventInfo(info);
    const TClonesArray & pvs = ev.get(h_pv);

    double weight = ev.weight;
//    //cout << "weight_histo = " <<weight<<endl;
    Int_t njets = js.GetEntries();
    hist("N_jets")->Fill(njets, weight);

    for (int i=0; i<njets; i++){
        baconhep::TJet* jets = (baconhep::TJet*)js[i];
        hist("pt")->Fill(jets->pt, weight);
        hist("eta")->Fill(jets->eta, weight);
        hist("phi")->Fill(jets->phi, weight);
        hist("MET")->Fill(eventInfo->pfMET, weight);
       // hist("nPu")->Fill(eventInfo->nPU, weight);
        hist("weight_histo")->Fill(weight, 1);
    }
        hist("nPu")->Fill(ev.nPU, weight);//for data only

    float nPrVer = 0;
    Int_t nvertices = pvs.GetEntries();
    hist("N_PV")->Fill(nvertices, weight);

    // require in the event that there is at least one reconstructed vertex
    if(nvertices>0) {
        // pick the first (i.e. highest sum pt) verte
        for (int i=0; i<nvertices; i++){
            baconhep::TVertex* vertices = (baconhep::TVertex*)pvs[i];
            // require that the vertex meets certain criteria
            if((fabs(vertices->z) < s_n_Pv_z) && (fabs(vertices->y) < s_n_Pv_xy) && (fabs(vertices->x) < s_n_Pv_xy) ){
                nPrVer++;
            }
        }
    }
    hist("N_PV_sel")->Fill(nPrVer, weight);

    baconhep::TJet* jet1 = (baconhep::TJet*)js[0];
    hist("pt_1")->Fill(ev.jet1_pt, weight);
    hist("eta_1")->Fill(jet1->eta, weight);

    baconhep::TJet* jet2 = (baconhep::TJet*)js[1];
    hist("pt_2")->Fill(ev.jet2_pt, weight);
    hist("eta_2")->Fill(jet2->eta, weight);


    hist("pt_ave")          ->Fill(ev.pt_ave, weight);
    hist("ptRaw_barrel")    ->Fill(ev.barreljet_ptRaw, weight);
    hist("ptRaw_probe")     ->Fill(ev.probejet_ptRaw , weight);
    hist("pt_barrel")   ->Fill(ev.barreljet_pt, weight);
    hist("pt_probe")    ->Fill(ev.probejet_pt , weight);
    hist("eta_barrel")  ->Fill(ev.barreljet_eta, weight);
    hist("eta_probe")   ->Fill(ev.probejet_eta , weight);
    hist("mpf")         ->Fill(ev.mpf_r, weight);
    hist("asym")        ->Fill(ev.asymmetry, weight);
    hist("r_rel")       ->Fill(ev.rel_r, weight);

    double deltaPhi = abs(jet1->phi - jet2->phi);
    hist("DeltaPhi_Jet1_Jet2")->Fill(deltaPhi, weight);

//     TString FileName[7] = {"pt_ave_hltDiPFJetAve40","pt_ave_hltDiPFJetAve80", "pt_ave_hltDiPFJetAve140", "pt_ave_hltDiPFJetAve200", "pt_ave_hltDiPFJetAve260", "pt_ave_hltDiPFJetAve320","pt_ave_hltDiPFJetAve400"};

//     HLT_DiPFJetAve140 = triggerBits[1]
//     HLT_DiPFJetAve200 = triggerBits[3]
//     HLT_DiPFJetAve260 = triggerBits[5]
//     HLT_DiPFJetAve320 = triggerBits[7]
//     HLT_DiPFJetAve40  = triggerBits[8]
//     HLT_DiPFJetAve400 = triggerBits[9]
//     HLT_DiPFJetAve500 = triggerBits[10]
//     HLT_DiPFJetAve60  = triggerBits[11]
//     HLT_DiPFJetAve80  = triggerBits[13]

    TString FileNameRun2[9] = {"pt_ave_hltDiPFJetAve140","pt_ave_hltDiPFJetAve200", "pt_ave_hltDiPFJetAve260", "pt_ave_hltDiPFJetAve320","pt_ave_hltDiPFJetAve40",  "pt_ave_hltDiPFJetAve400","pt_ave_hltDiPFJetAve500","pt_ave_hltDiPFJetAve60","pt_ave_hltDiPFJetAve80"};
    if(eventInfo->triggerBits[1]==1) hist(FileNameRun2[0])->Fill(ev.pt_ave, weight);
    if(eventInfo->triggerBits[3]==1) hist(FileNameRun2[1])->Fill(ev.pt_ave, weight);
    if(eventInfo->triggerBits[5]==1) hist(FileNameRun2[2])->Fill(ev.pt_ave, weight);
    if(eventInfo->triggerBits[7]==1) hist(FileNameRun2[3])->Fill(ev.pt_ave, weight);
    if(eventInfo->triggerBits[8]==1) hist(FileNameRun2[4])->Fill(ev.pt_ave, weight);
    if(eventInfo->triggerBits[9]==1) hist(FileNameRun2[5])->Fill(ev.pt_ave, weight);
    if(eventInfo->triggerBits[10]==1) hist(FileNameRun2[6])->Fill(ev.pt_ave, weight);
    if(eventInfo->triggerBits[11]==1) hist(FileNameRun2[7])->Fill(ev.pt_ave, weight);
    if(eventInfo->triggerBits[13]==1) hist(FileNameRun2[8])->Fill(ev.pt_ave, weight);


//     TString FileNameRun2HF[6] = {"HLT_DiPFJetAve100_HFJEC","HLT_DiPFJetAve160_HFJEC", "HLT_DiPFJetAve220_HFJEC", "HLT_DiPFJetAve300_HFJEC","HLT_DiPFJetAve60_HFJEC",  "HLT_DiPFJetAve80_HFJEC"};
//     if(eventInfo->triggerBits[0]==1) hist(FileNameRun2HF[0])->Fill(ev.pt_ave, weight);
//     if(eventInfo->triggerBits[2]==1) hist(FileNameRun2HF[1])->Fill(ev.pt_ave, weight);
//     if(eventInfo->triggerBits[4]==1) hist(FileNameRun2HF[2])->Fill(ev.pt_ave, weight);
//     if(eventInfo->triggerBits[6]==1) hist(FileNameRun2HF[3])->Fill(ev.pt_ave, weight);
//     if(eventInfo->triggerBits[12]==1) hist(FileNameRun2HF[4])->Fill(ev.pt_ave, weight);
//     if(eventInfo->triggerBits[14]==1) hist(FileNameRun2HF[5])->Fill(ev.pt_ave, weight);


    TVector2 pt, met;
    if(fabs(jet1->eta) < s_eta_barr){
        pt.Set(ev.jet1_pt * cos(jet1->phi),ev.jet1_pt * sin(jet1->phi));
        hist("generic_asym")->Fill((ev.jet2_pt - ev.jet1_pt) / (ev.jet2_pt + ev.jet1_pt), weight);
        //j(E_{jet}) = 1 + \frac{ {E^{\gamma}_{T}} \cdot { \slashed{E}_{T} } } { (E^{\gamma}_{T})^{2} } 
        hist("generic_mpf")->Fill(1 + (met.Px()*pt.Px() + met.Py()*pt.Py())/(pt.Px()*pt.Px() + pt.Py()*pt.Py()), weight);
        hist("generic_r_rel")->Fill(ev.jet2_pt / ev.jet1_pt, weight);
        hist("Rrel_vs_Npv")->Fill(nvertices,ev.jet2_pt / ev.jet1_pt);
//         hist("prof_Rrel_vs_Npv")->Fill(nvertices,ev.jet2_pt / ev.jet1_pt);

    }
    if(fabs(jet2->eta) < s_eta_barr){
        pt.Set(ev.jet2_pt * cos(jet2->phi),ev.jet2_pt * sin(jet2->phi));
        hist("generic_asym")->Fill((ev.jet1_pt - ev.jet2_pt) / (ev.jet2_pt + ev.jet1_pt), weight);
        //j(E_{jet}) = 1 + \frac{ {E^{\gamma}_{T}} \cdot { \slashed{E}_{T} } } { (E^{\gamma}_{T})^{2} }
        hist("generic_mpf")->Fill(1 + (met.Px()*pt.Px() + met.Py()*pt.Py())/(pt.Px()*pt.Px() + pt.Py()*pt.Py()), weight);
        hist("generic_r_rel")->Fill(ev.jet1_pt / ev.jet2_pt, weight);
        hist("Rrel_vs_Npv")->Fill(nvertices,ev.jet1_pt / ev.jet2_pt);
//         hist("prof_Rrel_vs_Npv")->Fill(nvertices,ev.jet1_pt / ev.jet2_pt);

    }

    baconhep::TJet* jet3 = (baconhep::TJet*)js[2];
    if (njets > 2){
        hist("pt_3")->Fill(ev.jet3_pt, weight);
        hist("eta_3")->Fill(jet3->eta, weight);
        hist("pt_rel")->Fill(ev.jet3_pt/(0.5*(ev.barreljet_pt + ev.probejet_pt )),weight);
        hist("generic_pt_rel")->Fill(ev.jet3_pt/(0.5*(ev.jet1_pt + ev.jet2_pt)),weight);
        hist("ptrel_vs_deltaphi")->Fill(ev.jet3_pt/(0.5*(ev.barreljet_pt + ev.probejet_pt)),deltaPhi);
    }
}
JECAnalysisHists::~JECAnalysisHists(){}
