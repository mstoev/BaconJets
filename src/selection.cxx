#include "UHH2/BaconJets/include/selection.h"

#include "UHH2/bacondataformats/interface/TJet.hh"

namespace uhh2BaconJets {

Selection::Selection(uhh2::Context & ctx) :
    context(ctx),
    event(0)
{
  h_jets = context.declare_event_input<TClonesArray>("Jet04");
  h_eventInfo = context.declare_event_input<baconhep::TEventInfo>("Info");
}

void Selection::SetEvent(uhh2::Event& evt)
{
   event = &evt;
   assert(event);
}

bool Selection::Trigger()
{
  assert(event);

  const TClonesArray & js = event->get(h_jets);
  const baconhep::TEventInfo & info = event->get(h_eventInfo);

  const baconhep::TJet * jet = dynamic_cast<const baconhep::TJet*>(js[0]);
  assert(jet);

  baconhep::TEventInfo* eventInfo= new baconhep::TEventInfo(info);
  assert(eventInfo);

  baconhep::TJet* jet1 = (baconhep::TJet*)js[0];
  baconhep::TJet* jet2 = (baconhep::TJet*)js[1];

 double avePt = (jet1->pt + jet2->pt)/2;

 bool trigger40fired = false;
 bool trigger80fired = false;
 bool trigger140fired = false;
 bool trigger200fired = false;
 bool trigger260fired = false;
 bool trigger320fired = false;
 bool trigger400fired = false;


 if(eventInfo->triggerBits[0]==1)  trigger40fired = true;
 if(eventInfo->triggerBits[1]==1)  trigger80fired = true;
 if(eventInfo->triggerBits[2]==1)  trigger140fired = true;
 if(eventInfo->triggerBits[3]==1)  trigger200fired = true;
 if(eventInfo->triggerBits[4]==1)  trigger260fired = true;
 if(eventInfo->triggerBits[5]==1)  trigger320fired = true;
 if(eventInfo->triggerBits[6]==1)  trigger400fired = true;

 if (avePt < 49) return false;

 if (avePt >= 49  && avePt < 92 && trigger40fired) return true;
 if (avePt >= 92  && avePt < 157 && trigger80fired) return true;
 if (avePt >= 157 && avePt < 218 && trigger140fired) return true;
 if (avePt >= 218 && avePt < 280 && trigger200fired) return true;
 if (avePt >= 280 && avePt < 346 && trigger260fired) return true;
 if (avePt >= 346 && avePt < 432 && trigger320fired) return true;
 if (avePt >= 432 && trigger400fired) return true;

//   double triggerThreshold[6] = {92,157,218,280,346,432};
//   bool trigger_fired = false;
//   for (int j = 1; j < 7; j++) {
//     if(eventInfo->triggerBits[j]==1) {
//         if(avePt < triggerThreshold[j-1]){
//             trigger_fired = true;
//             break;
//         }
//     }
//   }
//   if(eventInfo->triggerBits.none() || !trigger_fired) return false;

 return false;
}

bool Selection::DiJet()
{
  assert(event);

  const TClonesArray & js = event->get(h_jets);

  const baconhep::TJet * jet = dynamic_cast<const baconhep::TJet*>(js[0]);
  assert(jet);

  int njets = js.GetEntries();

  // njets >= 2
  if (njets>=2) return true;

  return false;
}
bool Selection::DiJetAdvancedAlpha05()
{
  assert(event);

  const TClonesArray & js = event->get(h_jets);

  const baconhep::TJet * jet = dynamic_cast<const baconhep::TJet*>(js[0]);
  assert(jet);

  int njets = js.GetEntries();

  // njets >= 2
  if (njets<2) return false;

  baconhep::TJet* jet1 = (baconhep::TJet*)js[0];
  baconhep::TJet* jet2 = (baconhep::TJet*)js[1];

  // at least one barrel jet
  if((fabs(jet1->eta)>=1.3) && (fabs(jet2->eta)>=1.3)) return false; 

  // delta phi < 2.7
  double deltaPhi = std::min(std::abs(double(jet1->phi) - double(jet2->phi)),2*M_PI-std::abs(jet2->phi - jet1->phi));
  if (deltaPhi < 2.7) return false;

  // |asymm| < 0.7
  if ((fabs(jet2->pt - jet1->pt) / (jet2->pt + jet1->pt)) > 0.7) return false;

  // p_t,rel < 0.5
  if (njets>2){
    baconhep::TJet* jet3 = (baconhep::TJet*)js[2];
    if ((2*(jet3->pt))/(jet1->pt + jet2->pt) > 0.5) return false;
  }

 return true;
}

bool Selection::DiJetAdvanced()
{
  assert(event);

  const TClonesArray & js = event->get(h_jets);

  const baconhep::TJet * jet = dynamic_cast<const baconhep::TJet*>(js[0]);
  assert(jet);

  int njets = js.GetEntries();

  // njets >= 2
  if (njets<2) return false;
  baconhep::TJet* jet1 = (baconhep::TJet*)js[0];
  baconhep::TJet* jet2 = (baconhep::TJet*)js[1];

  // p_t,rel < 0.2
  if (njets>2){
    baconhep::TJet* jet3 = (baconhep::TJet*)js[2];
    if ((2*(jet3->pt))/(jet1->pt + jet2->pt) > 0.2) return false;
  }

 return true;
}

bool Selection::FullSelection()
{
    return Trigger()&&DiJet()&&DiJetAdvanced();
}

Selection::~Selection()
{
}

}
