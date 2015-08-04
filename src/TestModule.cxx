#include <iostream>
#include <memory>
#include <stdlib.h>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "../include/JECAnalysisHists.h"

#include "UHH2/common/include/ObjectIdUtils.h"
#include "UHH2/common/include/Utils.h"
#include "UHH2/common/include/JetCorrections.h"

#include "UHH2/BaconJets/include/selection.h"
#include "UHH2/BaconJets/include/jet_corrections.h"
#include "UHH2/BaconJets/include/mc_weight.h"
#include "UHH2/BaconJets/include/constants.h"
//#include "UHH2/BaconJets/include/TSetTree.h"

#include "UHH2/bacondataformats/interface/TGenEventInfo.hh"
#include "UHH2/bacondataformats/interface/TJet.hh"
#include "UHH2/bacondataformats/interface/TEventInfo.hh"
#include "UHH2/bacondataformats/interface/BaconAnaDefs.hh"



#include "TClonesArray.h"
#include "TString.h"


TTree   *fCurrentTree;
Int_t   Runnr;
Int_t   Eventnr;
TFile *fCurrentTreeFile;
using namespace std;
using namespace uhh2;

namespace uhh2bacon {

  class TestModule: public AnalysisModule {
  public:
    explicit TestModule(Context & ctx);
    virtual bool process(Event & event) override;
    ~TestModule();

  private:

    Event::Handle<TClonesArray> h_jets  ;
    Event::Handle<baconhep::TEventInfo> h_eventInfo;
    Event::Handle<baconhep::TGenEventInfo> h_genInfo;

    std::unique_ptr<JetCorrector> jet_corrector;



    std::unique_ptr<Hists> h_nocuts, h_sel, h_dijet, h_match;
    //   std::vector<double> eta_range, pt_range, alpha_range;
    std::vector<JECAnalysisHists> h_pt_bins, h_eta_bins, h_pt_bins_a01, h_eta_bins_a01, h_eta_bins_a03, h_eta_bins_a04, h_noalpha_bins, h_mikko_barrel;

    Selection sel;
    JetCorrections jetcorr;
    McWeight mcweight;
    bool is_mc;


    //TSetTree cSetTree;

  };


  TestModule::TestModule(Context & ctx) :
    sel(ctx),
    jetcorr(ctx),
    mcweight(ctx)//,
    //cSetTree()
  {
    auto dataset_type = ctx.get("dataset_type");
    is_mc = dataset_type  == "MC";
    h_jets = ctx.declare_event_input<TClonesArray>("AK4PFCHS");
    h_eventInfo = ctx.declare_event_input<baconhep::TEventInfo>("Info");
    if(is_mc){ /// apply for MC only
      h_genInfo = ctx.declare_event_input<baconhep::TGenEventInfo>("GenEvtInfo");
    }

    h_nocuts.reset(new JECAnalysisHists(ctx,"noCuts"));
    h_dijet.reset(new JECAnalysisHists(ctx,"diJet"));
    h_match.reset(new JECAnalysisHists(ctx,"JetMatching"));
    h_sel.reset(new JECAnalysisHists(ctx,"Selection"));

    const std::vector<std::string> PHYS14_L123_MC = {"JetMETObjects/data/PHYS14_25_V2_L1FastJet_AK4PFchs.txt",
						     "JetMETObjects/data/PHYS14_25_V2_L2Relative_AK4PFchs.txt", "JetMETObjects/data/PHYS14_25_V2_L3Absolute_AK4PFchs.txt"};
    const std::vector<std::string> SUMMER15_L1 = {"UHH2/BaconJets/src/Summer15_50nsV2_DATA_L1FastJet_AK4PFchs.txt"};
    //jet_corrector.reset(new JetCorrector(PHYS14_L123_MC));
    //jet_corrector.reset(new JetCorrector(SUMMER15_L1));

    // int size = sizeof(eta_range)/sizeof(double); // to get size of string object
    // eta_range.size() // to get size of vector

    std::vector<std::string> pt_range_name;
    for( unsigned int i=0; i < pt_range.size(); ++i ){
      char pt_buffer [50];
      sprintf (pt_buffer, "%5.3f", pt_range[i]);
      pt_range_name.push_back(pt_buffer);
    }
    std::vector<std::string> eta_range_name;
    for( unsigned int j=0; j < eta_range.size(); ++j ){
      char eta_buffer [50];
      sprintf (eta_buffer, "%5.3f", eta_range[j]);
      eta_range_name.push_back(eta_buffer);
    }
    std::vector<std::string> alpha_range_name;
    for( unsigned int k=0; k < alpha_range.size(); ++k ){
      char alpha_buffer [50];
      sprintf (alpha_buffer, "%5.3f", alpha_range[k]);
      alpha_range_name.push_back(alpha_buffer);
    }
    std::vector<std::string> eta_range_mikko_barrel_name;
    for( unsigned int j=0; j < eta_range_mikko_barrel.size(); ++j ){
      char eta_buffer_mikko_barrel [50];
      sprintf (eta_buffer_mikko_barrel, "%5.3f", eta_range_mikko_barrel[j]);
      eta_range_mikko_barrel_name.push_back(eta_buffer_mikko_barrel);
    }

    cout << "alpha range "<<alpha_range_name[3]<<"eta range "<<eta_range_name[3]<< "pt range "<<pt_range_name[3]<< endl;
    // for( unsigned int k=0; k < alpha_range.size()-1; ++k ){
    //     for( unsigned int j=0; j < eta_range.size()-1; ++j ){
    for( unsigned int i=0; i < pt_range.size()-1; ++i ){
      h_pt_bins.push_back(JECAnalysisHists(ctx,(std::string)("a02/pt_"+pt_range_name[i]+"_"+pt_range_name[i+1])));
      h_pt_bins_a01.push_back(JECAnalysisHists(ctx,(std::string)("a01/pt_"+pt_range_name[i]+"_"+pt_range_name[i+1])));

    }
    //     }
    //   }
    for( unsigned int i=0; i < eta_range.size()-1; ++i ){
      h_eta_bins.push_back(JECAnalysisHists(ctx,(std::string)("a02/eta_"+eta_range_name[i]+"_"+eta_range_name[i+1])));
      h_eta_bins_a01.push_back(JECAnalysisHists(ctx,(std::string)("a01/eta_"+eta_range_name[i]+"_"+eta_range_name[i+1])));
      h_eta_bins_a03.push_back(JECAnalysisHists(ctx,(std::string)("a03/eta_"+eta_range_name[i]+"_"+eta_range_name[i+1])));
      h_eta_bins_a04.push_back(JECAnalysisHists(ctx,(std::string)("a04/eta_"+eta_range_name[i]+"_"+eta_range_name[i+1])));
    }

    // histos for the pT extrapolations
    for( unsigned int j=0; j < eta_range.size()-1; ++j ){
      for( unsigned int i=0; i < pt_range.size()-1; ++i ){
	h_noalpha_bins.push_back(JECAnalysisHists(ctx,(std::string)("/eta_"+eta_range_name[j]+"_"+eta_range_name[j+1]+"/pt_"+pt_range_name[i]+"_"+pt_range_name[i+1])));
      }
    }

    for( unsigned int j=0; j < eta_range_mikko_barrel.size()-1; ++j ){
      for( unsigned int i=0; i < pt_range.size()-1; ++i ){
	h_mikko_barrel.push_back(JECAnalysisHists(ctx,(std::string)("eta_"+eta_range_mikko_barrel_name[j]+"_"+eta_range_mikko_barrel_name[j+1]+"/pt_"+pt_range_name[i]+"_"+pt_range_name[i+1])));
      }
    }

    //h_eta_bin1.reset(new JECAnalysisHists(ctx,"eta_bin1"));
    //h_eta_bin2.reset(new JECAnalysisHists(ctx,"eta_bin2"));

  }

  TestModule::~TestModule() {
    //cSetTree.general();
  }

  bool TestModule::process(Event & event) {

    sel.SetEvent(event);
    jetcorr.SetEvent(event);
    mcweight.SetEvent(event);

    //jet_corrector->process(event);

    const TClonesArray & js = event.get(h_jets);
    baconhep::TJet* jet1 = (baconhep::TJet*)js[0];
    baconhep::TJet* jet2 = (baconhep::TJet*)js[1];
    Int_t njets = js.GetEntries();

    const baconhep::TEventInfo & info = event.get(h_eventInfo);
    baconhep::TEventInfo* eventInfo= new baconhep::TEventInfo(info);


    if(is_mc){ /// apply for MC only
      const baconhep::TGenEventInfo & geninfo = event.get(h_genInfo);
      baconhep::TGenEventInfo* genInfo= new baconhep::TGenEventInfo(geninfo);
      // event.weight = event.weight * genInfo->weight * mcweight.getPuReweighting() * mcweight.getEvReweighting();
      event.weight = event.weight * genInfo->weight * mcweight.getEvReweighting();

      //! matching from GEN to RECO
      if(!jetcorr.JetMatching()) return false;

      //! JER smearing
      if(!jetcorr.JetResolutionSmearer()) return false;
    }



    float probejet_eta = -99.;
    float probejet_pt = 0;
    float probejet_phi = -99.;
    float probejet_ptRaw = -99.;

    float barrel_eta = -99.;
    float barrel_pt = 0;
    float barrel_phi = -99.;
    float barrel_ptRaw = -99.;

    float rel_r = -99.;
    float mpf_r = -99.;
    float asymmetry = -99.;

    TVector2 pt, met;
    TVector2* MET = new TVector2(1,1);
    MET->SetMagPhi(eventInfo->pfMET ,eventInfo->pfMETphi);
    met.Set(eventInfo->pfMET * cos(eventInfo->pfMETphi),eventInfo->pfMET * sin(eventInfo->pfMETphi));

    int ran = rand();
    int numb = ran % 2 + 1;
    if ((fabs(jet1->eta)<s_eta_barr)&&(fabs(jet2->eta)<s_eta_barr)) {
      if(numb==1){
        probejet_eta = jet2->eta;
        probejet_pt = jet2->pt;
        probejet_phi = jet2->phi;
        probejet_ptRaw = jet2->ptRaw;

        barrel_eta = jet1->eta;
        barrel_pt = jet1->pt;
        barrel_phi = jet1->phi;
        barrel_ptRaw = jet1->ptRaw;

        asymmetry = (jet2->pt - jet1->pt)/(jet2->pt + jet1->pt);
        rel_r = jet2->pt / jet1->pt;

        pt.Set(jet1->pt * cos(jet1->phi),jet1->pt * sin(jet1->phi));
        mpf_r = 1 + (met.Px()*pt.Px() + met.Py()*pt.Py())/(pt.Px()*pt.Px() + pt.Py()*pt.Py());


      }
      if(numb==2){
        probejet_eta = jet1->eta;
        probejet_pt = jet1->pt;
        probejet_phi = jet1->phi;
        probejet_ptRaw = jet1->ptRaw;

        barrel_eta = jet2->eta;
        barrel_pt = jet2->pt;
        barrel_phi = jet2->phi;
        barrel_ptRaw = jet2->ptRaw;

        asymmetry = (jet1->pt - jet2->pt)/(jet1->pt + jet2->pt);
        rel_r = jet1->pt / jet2->pt;

        pt.Set(jet2->pt * cos(jet2->phi),jet2->pt * sin(jet2->phi));
        mpf_r = 1 + (met.Px()*pt.Px() + met.Py()*pt.Py())/(pt.Px()*pt.Px() + pt.Py()*pt.Py());

      }
    } else if ((fabs(jet1->eta)<s_eta_barr)||(fabs(jet2->eta)<s_eta_barr)){
      if(fabs(jet1->eta)<s_eta_barr){
        probejet_eta = jet2->eta;
        probejet_pt = jet2->pt;
        probejet_phi = jet2->phi;
        probejet_ptRaw = jet2->ptRaw;

        barrel_eta = jet1->eta;
        barrel_pt = jet1->pt;
        barrel_phi = jet1->phi;
        barrel_ptRaw = jet1->ptRaw;

        asymmetry = (jet2->pt - jet1->pt)/(jet2->pt + jet1->pt);
        rel_r = jet2->pt / jet1->pt;

        pt.Set(jet1->pt * cos(jet1->phi),jet1->pt * sin(jet1->phi));
        mpf_r = 1 + (met.Px()*pt.Px() + met.Py()*pt.Py())/(pt.Px()*pt.Px() + pt.Py()*pt.Py());


      }
      else{
        probejet_eta = jet1->eta;
        probejet_pt = jet1->pt;
        probejet_phi = jet1->phi;
        probejet_ptRaw = jet1->ptRaw;

        barrel_eta = jet2->eta;
        barrel_pt = jet2->pt;
        barrel_phi = jet2->phi;
        barrel_ptRaw = jet2->ptRaw;

        asymmetry = (jet1->pt - jet2->pt)/(jet1->pt + jet2->pt);
        rel_r = jet1->pt / jet2->pt;

        pt.Set(jet2->pt * cos(jet2->phi),jet2->pt * sin(jet2->phi));
        mpf_r = 1 + (met.Px()*pt.Px() + met.Py()*pt.Py())/(pt.Px()*pt.Px() + pt.Py()*pt.Py());
      }

    }
    float pt_ave = (jet1->pt + jet2->pt)/2;

    event.pt_ave = pt_ave;

    event.probejet_eta    = probejet_eta;
    event.probejet_phi    = probejet_phi;
    event.probejet_pt     = probejet_pt;
    event.probejet_ptRaw  = probejet_ptRaw;

    event.barreljet_eta   = barrel_eta;
    event.barreljet_phi   = barrel_phi;
    event.barreljet_pt    = barrel_pt;
    event.barreljet_ptRaw = barrel_ptRaw;

    event.asymmetry = asymmetry;
    event.rel_r = rel_r;
    event.mpf_r = mpf_r;

    float alpha = 0.;
    if (njets > 2) {
      baconhep::TJet* jet3 = (baconhep::TJet*)js[2];
      alpha = (2*(jet3->pt))/(jet1->pt + jet2->pt);
    }
    event.alpha = alpha;

    if(!sel.DiJet()) return false;

    h_nocuts->fill(event);


    if(!sel.DiJetAdvanced()) return false;

    h_dijet->fill(event);


    h_match->fill(event);

    if(!is_mc){
      if(!sel.Trigger()) return false;
    }

    // fill histos after dijet event selection
    //h_sel->fill(event);

    //   if(!sel.goodPVertex()) return false;


    //! for alpha < 0.2
    //   for( unsigned int k=0; k < alpha_range.size()-1; ++k ){
    //     if ((alpha>=alpha_range[k])&&(alpha<alpha_range[k+1])) {
    //         for( unsigned int j=0; j < eta_range.size()-1; ++j ){
    //             if ((fabs(probejet_eta)>=eta_range[j])&&(fabs(probejet_eta)<eta_range[j+1])) {
    for( unsigned int i=0; i < pt_range.size()-1; ++i ){
      if ((event.pt_ave >= pt_range[i])&&(event.pt_ave < pt_range[i+1])) {
	h_pt_bins[i].fill(event, ran);//j*pt_range.size()+i
	//                         h_pt_bins[j*(pt_range.size()-1)+i].fill(event, ran);//j*pt_range.size()+i

	//                         h_pt_bins[k*(eta_range.size()-1)*(pt_range.size()-1)+j*(pt_range.size()-1)+i].fill(event, ran);//j*pt_range.size()+i
	//cout <<"eta range = "<< eta_range[j]<<" - "<< eta_range[j+1]<< "pt range = "<< pt_range[i]<<" - "<< pt_range[i+1]<<endl;
	//cout <<"eta value = "<< fabs(probejet_eta) << " pt value = "<< probejet_pt <<endl;
      }
    }
    //             }
    //         }
    //     }
    //   }

    if (event.alpha < 0.4) {
      for( unsigned int i=0; i < eta_range.size()-1; ++i ){
	if ((fabs(event.probejet_eta)>=eta_range[i])&&(fabs(event.probejet_eta)<eta_range[i+1])) h_eta_bins_a04[i].fill(event, ran);
      }
    }

    if (event.alpha < 0.3) {
      for( unsigned int i=0; i < eta_range.size()-1; ++i ){
	if ((fabs(event.probejet_eta)>=eta_range[i])&&(fabs(event.probejet_eta)<eta_range[i+1])) h_eta_bins_a03[i].fill(event, ran);
      }
    }

    if (event.alpha < 0.2) {
      for( unsigned int i=0; i < eta_range.size()-1; ++i ){
	if ((fabs(event.probejet_eta)>=eta_range[i])&&(fabs(event.probejet_eta)<eta_range[i+1])) h_eta_bins[i].fill(event, ran);
      }
      h_sel->fill(event);
    }

    if (event.alpha < 0.2) {
      for( unsigned int j=0; j < eta_range.size()-1; ++j ){
	if ((fabs(event.probejet_eta)>=eta_range[j])&&(fabs(event.probejet_eta)<eta_range[j+1])) {
	  for( unsigned int i=0; i < pt_range.size()-1; ++i ){
	    if ((event.pt_ave>=pt_range[i])&&(event.pt_ave<pt_range[i+1])) {
	      h_noalpha_bins[j*(pt_range.size()-1)+i].fill(event, ran);//j*pt_range.size()+i
	    }
	  }
	}
      }

      for( unsigned int j=0; j < eta_range_mikko_barrel.size()-1; ++j ){
	if ((fabs(event.probejet_eta)>=eta_range_mikko_barrel[j])&&(fabs(event.probejet_eta)<eta_range_mikko_barrel[j+1])) {
	  for( unsigned int i=0; i < pt_range.size()-1; ++i ){
	    if ((event.pt_ave>=pt_range[i])&&(event.pt_ave<pt_range[i+1])) {
	      h_mikko_barrel[j*(pt_range.size()-1)+i].fill(event, ran);//j*pt_range.size()+i
	    }
	  }
	}
      }

    }


    //! for alpha < 0.1

    if (event.alpha < 0.1) {
      for( unsigned int i=0; i < pt_range.size()-1; ++i ){
	if ((event.pt_ave >= pt_range[i])&&(event.pt_ave < pt_range[i+1])) h_pt_bins_a01[i].fill(event, ran);//j*pt_range.size()+i
      }
      for( unsigned int i=0; i < eta_range.size()-1; ++i ){
	if ((fabs(event.probejet_eta)>=eta_range[i])&&(fabs(event.probejet_eta)<eta_range[i+1])) h_eta_bins_a01[i].fill(event, ran);
      }
    }

    //     if ((eventInfo->runNum == 208307) || (eventInfo->runNum == 208339)|| (eventInfo->runNum == 208341)|| (eventInfo->runNum == 208351)|| (eventInfo->runNum == 208353)){
    //     if ((fabs(jet1->eta)>0.)&&(fabs(jet1->eta)<1.3)) {
    //                 cout << "1)eta <1.3: jet # 1 : pt53: "<<jet1->pt<<" pt_raw53: "<<jet1->ptRaw<<" f53=pt53/pt_raw53: "<<jet1->pt/jet1->ptRaw<<'\n';
    //                 cout << "2)eta <1.3: jet # 1 : pt53: "<<jet1->pt<<" pt_raw53: "<<jet1->ptRaw<<" f53=pt_raw53/pt53: "<<jet1->ptRaw/jet1->pt<<'\n';    }
    //cSetTree.fillTree(eventInfo->evtNum, eventInfo->runNum, pt_ave, probejet_eta, barrel_eta, probejet_phi, barrel_phi, probejet_pt, barrel_pt, jet1->pt, jet2->pt, jet1->phi, jet2->phi, jet1->eta, jet2->eta, jet1->ptRaw, jet2->ptRaw, probejet_ptRaw, barrel_ptRaw);

    //     }

    return true;
  }

  // as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
  // make sure the ExampleModule is found by class name. This is ensured by this macro:
  UHH2_REGISTER_ANALYSIS_MODULE(TestModule)

}
