#include "UHH2/BaconJets/include/JECAnalysisHists.h"
#include "UHH2/core/include/Event.h"
#include "UHH2/bacondataformats/interface/TJet.hh"
#include "UHH2/bacondataformats/interface/TEventInfo.hh"
#include "UHH2/bacondataformats/interface/BaconAnaDefs.hh"
#include "TH1F.h"
#include "TH2F.h"
#include <iostream>
using namespace std;
using namespace uhh2;
using namespace baconhep;
JECAnalysisHists::JECAnalysisHists(Context & ctx, const string & dirname): Hists(ctx, dirname){
    // book all histograms here
    // jets
    TH1::SetDefaultSumw2();

    book<TH1F>("N_jets", "N_{jets}", 20, -0.5, 19.5);
    book<TH1F>("pt","p_{T} all jets",100,0,1500);
    book<TH1F>("eta","#eta all jets",100,-5,5);
    book<TH1F>("phi","#phi all jets",50,-M_PI,M_PI);
    book<TH1F>("MET","MET all jets",400,0,400);

    book<TH1F>("pt_1","p_{T} jet 1",100,0,1500);
    book<TH1F>("eta_1","#eta jet 1",100,-5,5);

    book<TH1F>("pt_2","p_{T} jet 2",100,0,1500);
    book<TH1F>("eta_2","#eta jet 2",100,-5,5);

    book<TH1F>("pt_3","p_{T} jet 3",100,0,1500);
    book<TH1F>("eta_3","#eta jet 3",100,-5,5);

    book<TH1F>("pt_barrel","p_{T} barrel jet",100,0,1500);
    book<TH1F>("pt_probe","p_{T} probe jet",100,0,1500);
    book<TH1F>("pt_ave","p_{T} ave jet",100,0,1500);
    book<TH1F>("pt_rel","p_{T} jet 3 / #overline{P}_{T}", 50, 0, 1);

    book<TH1F>("asym","asymmetrie jet 1 and jet 2 in loop",100,-1,1);
    book<TH1F>("generic_asym","generic asymmetrie jet 1 and jet 2",100,-1,1);
    book<TH1F>("mpf","MPF response",100,0.5,1.5);
    book<TH1F>("generic_mpf","generic MPF response",100,0.5,1.5);
    book<TH1F>("r_rel","R_{rel}",100,0.5,1.5);
    book<TH1F>("generic_r_rel","generic R_{rel}",100,0.5,1.5);

    book<TH1F>("DeltaPhi_Jet1_Jet2", "#Delta#Phi(first jet, second jet)", 100, 0, 7);
    book<TH1F>("generic_pt_rel","generic p_{T} jet 3 / #overline{P}_{T}", 50, 0, 1);
    book<TH2F>("ptrel_vs_deltaphi","delta phi vs pt_rel", 50, 0, 1 ,50, 0, 3.14);

    book<TH1F>("pt_ave_hltDiPFJetAve40","p_{T} ave40 jet",500,0,500);
    book<TH1F>("pt_ave_hltDiPFJetAve80","p_{T} ave80 jet",500,0,500);
    book<TH1F>("pt_ave_hltDiPFJetAve140","p_{T} ave140 jet",500,0,500);
    book<TH1F>("pt_ave_hltDiPFJetAve200","p_{T} ave200 jet",500,0,500);
    book<TH1F>("pt_ave_hltDiPFJetAve260","p_{T} ave260 jet",500,0,500);
    book<TH1F>("pt_ave_hltDiPFJetAve320","p_{T} ave320 jet",500,0,500);
    book<TH1F>("pt_ave_hltDiPFJetAve400","p_{T} ave400 jet",500,0,500);

    h_jets = ctx.get_handle<TClonesArray>("Jet05");
    h_eventInfo = ctx.get_handle<baconhep::TEventInfo>("Info");
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

    double weight = ev.weight;

    Int_t njets = js.GetEntries();
    hist("N_jets")->Fill(njets, weight);

    for (int i=0; i<njets; i++){
        baconhep::TJet* jets = (baconhep::TJet*)js[i];
        hist("pt")->Fill(jets->pt, weight);
        hist("eta")->Fill(jets->eta, weight);
        hist("phi")->Fill(jets->phi, weight);
        hist("MET")->Fill(eventInfo->pfMET, weight);
    }

    baconhep::TJet* jet1 = (baconhep::TJet*)js[0];
  //  std::cout <<"histo: corr jet1->pt " << jet1->pt<< std::endl;
    hist("pt_1")->Fill(jet1->pt, weight);
    hist("eta_1")->Fill(jet1->eta, weight);

    baconhep::TJet* jet2 = (baconhep::TJet*)js[1];
    hist("pt_2")->Fill(jet2->pt, weight);
    hist("eta_2")->Fill(jet2->eta, weight);

    double pt_ave = (jet1->pt + jet2->pt)/2;
    hist("pt_ave")->Fill(pt_ave, weight);

    double deltaPhi = abs(jet1->phi - jet2->phi);
    hist("DeltaPhi_Jet1_Jet2")->Fill(deltaPhi, weight);

    TString FileName[7] = {"pt_ave_hltDiPFJetAve40","pt_ave_hltDiPFJetAve80", "pt_ave_hltDiPFJetAve140", "pt_ave_hltDiPFJetAve200", "pt_ave_hltDiPFJetAve260", "pt_ave_hltDiPFJetAve320","pt_ave_hltDiPFJetAve400"};
    for (int j = 0; j < 7; j++) {
        if(eventInfo->triggerBits[j]==1) {
            hist(FileName[j])->Fill(pt_ave, weight);
        }
    }// end loop over pt for diff. triggers

    double barreljet = 0.0;
    double barreljet_eta = 0.0;
    double probejet = 0.0;
    double asymmetry = 0.0;

    TVector2 pt, met;
    met.Set(eventInfo->pfMET * cos(eventInfo->pfMETphi),eventInfo->pfMET * sin(eventInfo->pfMETphi));

    LorentzVector ETmiss(0,0,0,0);
    ETmiss.SetPt(eventInfo->pfMET);
    ETmiss.SetPhi(eventInfo->pfMETphi);

    int numb = rand % 2 + 1;
    // barrel region |eta|<1.3
    if ((fabs(jet1->eta)<1.3)&&(fabs(jet2->eta)<1.3)) {
        if(numb==1){
            if(fabs(jet1->eta)<1.3){
                barreljet += jet1->pt;
		barreljet_eta += jet1->eta;
                probejet += jet2->pt;
                asymmetry += (jet2->pt - jet1->pt)/(jet2->pt + jet1->pt);
                pt.Set(jet1->pt * cos(jet1->phi),jet1->pt * sin(jet1->phi));
                //hist("mpf")->Fill(1 + (met.Px()*pt.Px() + met.Py()*pt.Py())/(pt.Px()*pt.Px() + pt.Py()*pt.Py()), weight);
		//hist("mpf")->Fill(1+(cos(ETmiss.phi()-barreljet_eta)* ETmiss.pt()/barreljet  ), weight);
                hist("pt_barrel")->Fill(jet1->pt, weight);
                hist("pt_probe")->Fill(jet2->pt, weight);
                hist("asym")->Fill(((jet2->pt - jet1->pt)/(jet2->pt + jet1->pt)), weight);
                //hist("mpf")->Fill(1 + ((eventInfo->pfMET * jet1->pt)/pow(jet1->pt,2)), weight);
                hist("r_rel")->Fill(jet2->pt / jet1->pt, weight);

            }
        } if(numb==2){
            if(fabs(jet2->eta)<1.3){
                barreljet += jet2->pt;
		barreljet_eta += jet2->eta;
                probejet += jet1->pt;
                asymmetry += (jet1->pt - jet2->pt)/(jet1->pt + jet2->pt);
                pt.Set(jet2->pt * cos(jet2->phi),jet2->pt * sin(jet2->phi));
                //hist("mpf")->Fill(1 + (met.Px()*pt.Px() + met.Py()*pt.Py())/(pt.Px()*pt.Px() + pt.Py()*pt.Py()), weight);
		//hist("mpf")->Fill(1+(cos(ETmiss.phi()-barreljet_eta)* ETmiss.pt()/barreljet  ), weight);
                hist("pt_barrel")->Fill(jet2->pt, weight);
                hist("pt_probe")->Fill(jet1->pt, weight);
                hist("asym")->Fill(((jet1->pt - jet2->pt)/(jet1->pt + jet2->pt)), weight);
                //hist("mpf")->Fill(1 + ((eventInfo->pfMET * jet2->pt)/pow(jet2->pt,2)), weight);
                hist("r_rel")->Fill(jet1->pt / jet2->pt, weight);

            }
        }
    } else if ((fabs(jet1->eta)<1.3)||(fabs(jet2->eta)<1.3)){
        if(fabs(jet1->eta)<1.3){
            barreljet += jet1->pt;
	    barreljet_eta += jet1->eta;
            probejet += jet2->pt;
            asymmetry += (jet2->pt - jet1->pt)/(jet2->pt + jet1->pt);
            pt.Set(jet1->pt * cos(jet1->phi),jet1->pt * sin(jet1->phi));
            //hist("mpf")->Fill(1 + (met.Px()*pt.Px() + met.Py()*pt.Py())/(pt.Px()*pt.Px() + pt.Py()*pt.Py()), weight);
	    //hist("mpf")->Fill(1+(cos(ETmiss.phi()-barreljet_eta)* ETmiss.pt()/barreljet  ), weight);
            hist("pt_barrel")->Fill(jet1->pt, weight);
            hist("pt_probe")->Fill(jet2->pt, weight);
            hist("asym")->Fill(((jet2->pt - jet1->pt)/(jet2->pt + jet1->pt)), weight);
            //hist("mpf")->Fill(1 + ((eventInfo->pfMET * jet1->pt)/pow(jet1->pt,2)), weight);
            hist("r_rel")->Fill(jet2->pt / jet1->pt, weight);

        } if(fabs(jet2->eta)<1.3){
            barreljet += jet2->pt;
	    barreljet_eta += jet2->eta;
            probejet += jet1->pt;
            asymmetry += (jet1->pt - jet2->pt)/(jet1->pt + jet2->pt);
            pt.Set(jet2->pt * cos(jet2->phi),jet2->pt * sin(jet2->phi));
            //hist("mpf")->Fill(1 + (met.Px()*pt.Px() + met.Py()*pt.Py())/(pt.Px()*pt.Px() + pt.Py()*pt.Py()), weight);
	    //hist("mpf")->Fill(1+(cos(ETmiss.phi()-barreljet_eta)* ETmiss.pt()/barreljet  ), weight);
            hist("pt_barrel")->Fill(jet2->pt, weight);
            hist("pt_probe")->Fill(jet1->pt, weight);
            hist("asym")->Fill(((jet1->pt - jet2->pt)/(jet1->pt + jet2->pt)), weight);
            //hist("mpf")->Fill(1 + ((eventInfo->pfMET * jet2->pt)/pow(jet2->pt,2)), weight);
            hist("r_rel")->Fill(jet1->pt / jet2->pt, weight);

        }
    }

    hist("mpf")->Fill(1+(cos(ETmiss.phi()-barreljet_eta)* ETmiss.pt()/barreljet  ), weight);



    if(fabs(jet1->eta)<1.3){
        pt.Set(jet1->pt * cos(jet1->phi),jet1->pt * sin(jet1->phi));
        hist("generic_asym")->Fill((jet2->pt - jet1->pt) / (jet2->pt + jet1->pt), weight);
        //j(E_{jet}) = 1 + \frac{ {E^{\gamma}_{T}} \cdot { \slashed{E}_{T} } } { (E^{\gamma}_{T})^{2} } 
        hist("generic_mpf")->Fill(1 + (met.Px()*pt.Px() + met.Py()*pt.Py())/(pt.Px()*pt.Px() + pt.Py()*pt.Py()), weight);
        hist("generic_r_rel")->Fill(jet2->pt / jet1->pt, weight);

    }
    if(fabs(jet2->eta)<1.3){
        pt.Set(jet2->pt * cos(jet2->phi),jet2->pt * sin(jet2->phi));
        hist("generic_asym")->Fill((jet1->pt - jet2->pt) / (jet2->pt + jet1->pt), weight);
        //j(E_{jet}) = 1 + \frac{ {E^{\gamma}_{T}} \cdot { \slashed{E}_{T} } } { (E^{\gamma}_{T})^{2} }
        hist("generic_mpf")->Fill(1 + (met.Px()*pt.Px() + met.Py()*pt.Py())/(pt.Px()*pt.Px() + pt.Py()*pt.Py()), weight);
        hist("generic_r_rel")->Fill(jet1->pt / jet2->pt, weight);

    }

    baconhep::TJet* jet3 = (baconhep::TJet*)js[2];
    if (njets>2){
        hist("pt_3")->Fill(jet3->pt, weight);
        hist("eta_3")->Fill(jet3->eta, weight);
        hist("pt_rel")->Fill(jet3->pt/(0.5*(barreljet + probejet)),weight);
        hist("generic_pt_rel")->Fill(jet3->pt/(0.5*(jet1->pt + jet2->pt)),weight);
        hist("ptrel_vs_deltaphi")->Fill(jet3->pt/(0.5*(barreljet + probejet)),deltaPhi);
    }
}
JECAnalysisHists::~JECAnalysisHists(){}
