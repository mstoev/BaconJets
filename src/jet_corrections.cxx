#include "UHH2/bacon/include/jet_corrections.h"

#include "UHH2/bacondataformats/interface/TJet.hh"
using namespace std;
namespace uhh2bacon {

JetCorrections::JetCorrections(uhh2::Context & ctx) :
    context(ctx),
    event(0)
{
  h_jets = context.declare_event_input<TClonesArray>("Jet05");
  h_eventInfo = context.declare_event_input<baconhep::TEventInfo>("Info");
 // h_jetsout = context.declare_event_output<TClonesArray>("Jet05");
}

void JetCorrections::SetEvent(uhh2::Event& evt)
{
   event = &evt;
   assert(event);
}

bool JetCorrections::JetMatching()
{
  assert(event);

  //const TClonesArray & jsout = event->set(h_jetsout);
  const TClonesArray & js = event->get(h_jets);
  const baconhep::TEventInfo & info = event->get(h_eventInfo);


  baconhep::TEventInfo* eventInfo= new baconhep::TEventInfo(info);
  assert(eventInfo);

  int njets = js.GetEntries();

  // njets >= 2
  if (njets<2) return false;

  baconhep::TJet* jet1 = (baconhep::TJet*)js[0];
  baconhep::TJet* jet2 = (baconhep::TJet*)js[1];

  double delra_R_jet1 = pow(pow(jet1->genphi - jet1->phi,2) + pow(jet1->geneta - jet1->eta,2),0.5);
  double delra_R_jet2 = pow(pow(jet2->genphi - jet2->phi,2) + pow(jet2->geneta - jet2->eta,2),0.5);

  if ((delra_R_jet1 > 0.3) || (delra_R_jet2 > 0.3)) return false;
  //std::cout << " runNum: "<< eventInfo->runNum<< " evtNum: "<< eventInfo->evtNum <<" delra_R_jet1 = "<< delra_R_jet1 << std::endl;

 return true;
}

bool JetCorrections::JetResolutionSmearer()
{
  assert(event);

  const TClonesArray & js = event->get(h_jets);
  const baconhep::TEventInfo & info = event->get(h_eventInfo);

  baconhep::TEventInfo* eventInfo= new baconhep::TEventInfo(info);
  assert(eventInfo);

  //numbers taken from https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
  // from 8TeV JER measurement.
  constexpr const size_t n = 7;
  static float eta_hi[n] = {0.5, 1.1, 1.7, 2.3, 2.8, 3.2, 5.0};
  static float c_nominal[n] = {1.079, 1.099, 1.121, 1.208, 1.254, 1.395, 1.056};
//   static float c_up[n] = {1.105, 1.127, 1.150, 1.254, 1.316, 1.458, 1.247};
//   static float c_down[n] = {1.053, 1.071, 1.092, 1.162, 1.192, 1.332, 0.865};

  Int_t njets = js.GetEntries();
  baconhep::TJet* jet1 = (baconhep::TJet*)js[0];

  for(int i=0; i < njets; ++i) {
    baconhep::TJet* jet = (baconhep::TJet*)js[i];
    float genpt = jet->genpt;

    //ignore unmatched jets (which have zero vector) or jets with very low pt:
    if(genpt < 15.0) continue;

    float recopt = jet->pt;
    float abseta = fabs(jet->eta);
    size_t ieta = 0;
    while(ieta < n && eta_hi[ieta] < abseta) ++ieta;
    if(ieta == n) ieta = n-1;
    float c;
    c = c_nominal[ieta];
    float new_pt = std::max(0.0f, genpt + c * (recopt - genpt));

    if (i == 0)std::cout <<"jet corr: jet1->pt " << jet1->pt<< std::endl;
    if (i == 0)std::cout <<"jet corr: jet->pt " << new_pt<< std::endl;

    jet->pt = new_pt;

    }
    return true;
}

bool JetCorrections::FullJetCorrections()
{
    return JetMatching()&&JetResolutionSmearer();

}

JetCorrections::~JetCorrections()
{
}

}