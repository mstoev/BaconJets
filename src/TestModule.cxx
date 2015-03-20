#include <iostream>
#include <memory>
#include <stdlib.h>

#include "UHH2/core/include/AnalysisModule.h"
#include "UHH2/core/include/Event.h"
#include "../include/JECAnalysisHists.h"

#include "UHH2/BaconJets/include/selection.h"

#include "UHH2/bacondataformats/interface/TJet.hh"
#include "UHH2/bacondataformats/interface/TEventInfo.hh"
#include "UHH2/bacondataformats/interface/BaconAnaDefs.hh"

#include "TClonesArray.h"
#include "TString.h"

using namespace std;
using namespace uhh2;

namespace uhh2BaconJets {

class TestModule: public AnalysisModule {
public:

    explicit TestModule(Context & ctx);
    virtual bool process(Event & event) override;

private:

  Event::Handle<TClonesArray> h_jets;
  Event::Handle<baconhep::TEventInfo> h_eventInfo;

  std::unique_ptr<Hists> h_nocuts, h_sel, h_dijet, h_test;
  std::vector<double> eta_range;
  std::vector<JECAnalysisHists> h_eta_bins_alpha05, h_eta_bins;

  Selection sel;

};


TestModule::TestModule(Context & ctx) :
    sel(ctx)
{
  h_jets = ctx.declare_event_input<TClonesArray>("Jet04");
  h_eventInfo = ctx.declare_event_input<baconhep::TEventInfo>("Info");

  h_nocuts.reset(new JECAnalysisHists(ctx,"noCuts"));
  h_dijet.reset(new JECAnalysisHists(ctx,"diJet"));
  h_sel.reset(new JECAnalysisHists(ctx,"Selection"));

  eta_range = {0, 0.261, 0.522, 0.763, 0.957, 1.131, 1.305, 1.479, 1.93, 2.322, 2.411, 2.5, 2.853, 2.964, 3.139, 3.489, 5.191};
  // int size = sizeof(eta_range)/sizeof(double); // to get size of string object
  // eta_range.size() // to get size of vector

  std::vector<std::string> range_name;
  for( unsigned int i=0; i < eta_range.size(); ++i ){
    char buffer [50];
    sprintf (buffer, "%5.3f", eta_range[i]);
    range_name.push_back(buffer);
  }
  for( unsigned int i=0; i < eta_range.size()-1; ++i ){
    h_eta_bins_alpha05.push_back(JECAnalysisHists(ctx,(std::string)("alpha05_eta_"+range_name[i]+"_"+range_name[i+1])));
  }
  for( unsigned int i=0; i < eta_range.size()-1; ++i ){
    h_eta_bins.push_back(JECAnalysisHists(ctx,(std::string)("eta_"+range_name[i]+"_"+range_name[i+1])));
  }
  //h_eta_bin1.reset(new JECAnalysisHists(ctx,"eta_bin1"));
  //h_eta_bin2.reset(new JECAnalysisHists(ctx,"eta_bin2"));
}


bool TestModule::process(Event & event) {

  sel.SetEvent(event);

  const TClonesArray & js = event.get(h_jets);

  baconhep::TJet* jet1 = (baconhep::TJet*)js[0];
  baconhep::TJet* jet2 = (baconhep::TJet*)js[1];

  double probejet_eta = -99.;

  int ran = rand();
  int numb = ran % 2 + 1;
  if ((fabs(jet1->eta)<1.3)&&(fabs(jet2->eta)<1.3)) {
    if(numb==1){
        probejet_eta = jet2->eta;
    }
    if(numb==2){
        probejet_eta = jet1->eta;
    }
  } else if ((fabs(jet1->eta)<1.3)||(fabs(jet2->eta)<1.3)){
    if(fabs(jet1->eta)<1.3){
        probejet_eta = jet2->eta;
    }
    else{
        probejet_eta = jet1->eta;
    }
  }

  if(!sel.DiJet()) return false;

  h_nocuts->fill(event);

  if(!sel.DiJetAdvancedAlpha05()) return false;

  for( unsigned int i=0; i < eta_range.size()-1; ++i ){
    if ((fabs(probejet_eta)>=eta_range[i])&&(fabs(probejet_eta)<eta_range[i+1])) h_eta_bins_alpha05[i].fill(event, ran);

  }

  if(!sel.DiJetAdvanced()) return false;

  h_dijet->fill(event);

  if(!sel.Trigger()) return false;

  // fill histos after dijet event selection
  h_sel->fill(event);



  for( unsigned int i=0; i < eta_range.size()-1; ++i ){
    if ((fabs(probejet_eta)>=eta_range[i])&&(fabs(probejet_eta)<eta_range[i+1])) h_eta_bins[i].fill(event, ran);

  }


  return true;
}

// as we want to run the ExampleCycleNew directly with AnalysisModuleRunner,
// make sure the ExampleModule is found by class name. This is ensured by this macro:
UHH2_REGISTER_ANALYSIS_MODULE(TestModule)

}
