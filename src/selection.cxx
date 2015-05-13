#include "UHH2/BaconJets/include/selection.h"

#include "UHH2/bacondataformats/interface/TJet.hh"
#include "UHH2/bacondataformats/interface/TVertex.hh"
namespace uhh2bacon {

Selection::Selection(uhh2::Context & ctx) :
    context(ctx),
    event(0)
{
  h_jets = context.declare_event_input<TClonesArray>("Jet05");
  h_eventInfo = context.declare_event_input<baconhep::TEventInfo>("Info");
  h_pv = context.declare_event_input<TClonesArray>("PV");

//   h_jets = context.declare_event_autput<TClonesArray>("Jet05");
//   h_eventInfo = context.declare_event_autput<baconhep::TEventInfo>("Info");
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
//   const baconhep::TJet * jet = dynamic_cast<const baconhep::TJet*>(js[0]);
//   assert(jet);

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

 if (avePt < 66) return false;

 if (avePt >= 66  && avePt < 107 && trigger40fired) return true;
 if (avePt >= 107 && avePt < 191 && trigger80fired) return true;
 if (avePt >= 191 && avePt < 240 && trigger140fired) return true;
 if (avePt >= 240 && avePt < 306 && trigger200fired) return true;
 if (avePt >= 306 && avePt < 379 && trigger260fired) return true;
 if (avePt >= 379 && avePt < 468 && trigger320fired) return true;
 if (avePt >= 468 && trigger400fired) return true;

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
  const baconhep::TEventInfo & info = event->get(h_eventInfo);
  baconhep::TEventInfo* eventInfo= new baconhep::TEventInfo(info);
  assert(eventInfo);
//   const baconhep::TJet * jet = dynamic_cast<const baconhep::TJet*>(js[0]);
//   assert(jet);

  int njets = js.GetEntries();

  // njets >= 2
  if (njets>=2) return true;
  //std::cout << eventInfo->triggerBits << std::endl;
  return false;
}
bool Selection::DiJetAdvanced()
{
  assert(event);

  const TClonesArray & js = event->get(h_jets);

//   const baconhep::TJet * jet = dynamic_cast<const baconhep::TJet*>(js[0]);
//   assert(jet);

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
  if (fabs((jet2->pt - jet1->pt) / (jet2->pt + jet1->pt)) > 0.7) return false;

  // p_t,rel < 0.2
  if (njets>2){
    baconhep::TJet* jet3 = (baconhep::TJet*)js[2];
    //if ((2*(jet3->pt))/(jet1->pt + jet2->pt) > 0.2) return false;
  }

 return true;
}

bool Selection::goodPVertex()
{
  assert(event);
  const baconhep::TEventInfo & info = event->get(h_eventInfo);
  baconhep::TEventInfo* eventInfo= new baconhep::TEventInfo(info);
  const TClonesArray & pvs = event->get(h_pv);
    //        std::cout << "good prim vertex  = "<< eventInfo->pvz<<std::endl;
  Int_t nvertices = pvs.GetEntries();
  //std::cout << "1 number of  vertex  = "<< nvertices<<std::endl;
  // require in the event that there is at least one reconstructed vertex
  if(nvertices<=0) return false;
  float nPrVer = 0;
  // pick the first (i.e. highest sum pt) verte
  //const reco::Vertex* theVertex=&(vertices_h->front());

  for (int i=0; i<nvertices; i++){
        baconhep::TVertex* vertices = (baconhep::TVertex*)pvs[i];
        // require that the vertex meets certain criteria
        
        if((vertices->ndof>5) && (fabs(vertices->z)<24.0)  && (vertices->chi2>0.05)){
            nPrVer++;
         //   std::cout << "prim vertex  = "<< vertices->z<<std::endl;
        }
  }
         //std::cout << "2 number of  vertex  = "<< nPrVer<<std::endl;

 return true;
}


bool Selection::FullSelection()
{
    return Trigger()&&DiJet()&&DiJetAdvanced()&&goodPVertex();

}

Selection::~Selection()
{
}

}
