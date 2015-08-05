#include "UHH2/BaconJets/include/selection.h"

#include "UHH2/bacondataformats/interface/TJet.hh"
#include "UHH2/bacondataformats/interface/TVertex.hh"
#include "UHH2/BaconJets/include/constants.h"

#include "TVector2.h"

namespace uhh2bacon {

Selection::Selection(uhh2::Context & ctx) :
    context(ctx),
    event(0)
{
  h_jets = context.declare_event_input<TClonesArray>("AK4PFCHS");
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


bool Selection::Trigger(uhh2::Event& evt)
{
    assert(event);

    const TClonesArray & js = event->get(h_jets);
    const baconhep::TEventInfo & info = event->get(h_eventInfo);
//   const baconhep::TJet * jet = dynamic_cast<const baconhep::TJet*>(js[0]);
//   assert(jet);

    baconhep::TEventInfo* eventInfo= new baconhep::TEventInfo(info);
    assert(eventInfo);



    bool trigger40fired = false;
    bool trigger60fired = false;
    bool trigger80fired = false;
    bool trigger140fired = false;
    bool trigger200fired = false;
    bool trigger260fired = false;
    bool trigger320fired = false;
    bool trigger400fired = false;
    bool trigger500fired = false;
    ///triggerNames  briggerBits
//     HLT_DiPFJetAve140 = triggerBits[1]
//     HLT_DiPFJetAve200 = triggerBits[3]
//     HLT_DiPFJetAve260 = triggerBits[5]
//     HLT_DiPFJetAve320 = triggerBits[7]
//     HLT_DiPFJetAve40  = triggerBits[8]
//     HLT_DiPFJetAve400 = triggerBits[9]
//     HLT_DiPFJetAve500 = triggerBits[10]
//     HLT_DiPFJetAve60  = triggerBits[11]
//     HLT_DiPFJetAve80  = triggerBits[13]


    if(eventInfo->triggerBits[8]==1)  trigger40fired = true;
    if(eventInfo->triggerBits[11]==1) trigger60fired = true;
    if(eventInfo->triggerBits[13]==1) trigger80fired = true;
    if(eventInfo->triggerBits[1]==1)  trigger140fired = true;
    if(eventInfo->triggerBits[3]==1)  trigger200fired = true;
    if(eventInfo->triggerBits[5]==1)  trigger260fired = true;
    if(eventInfo->triggerBits[7]==1)  trigger320fired = true;
    if(eventInfo->triggerBits[9]==1)  trigger400fired = true;
    if(eventInfo->triggerBits[10]==1) trigger500fired = true;


    if (evt.pt_ave < s_Pt_Ave40_cut) return false;
    if (evt.pt_ave >= s_Pt_Ave40_cut  && evt.pt_ave < s_Pt_Ave60_cut  && trigger40fired) return true;
    if (evt.pt_ave >= s_Pt_Ave60_cut  && evt.pt_ave < s_Pt_Ave80_cut  && trigger60fired) return true;
    if (evt.pt_ave >= s_Pt_Ave80_cut  && evt.pt_ave < s_Pt_Ave140_cut && trigger80fired) return true;
    if (evt.pt_ave >= s_Pt_Ave140_cut && evt.pt_ave < s_Pt_Ave200_cut /*&& trigger140fired*/) return true;
    if (evt.pt_ave >= s_Pt_Ave200_cut && evt.pt_ave < s_Pt_Ave260_cut && trigger200fired) return true;
    if (evt.pt_ave >= s_Pt_Ave260_cut && evt.pt_ave < s_Pt_Ave320_cut && trigger260fired) return true;
    if (evt.pt_ave >= s_Pt_Ave320_cut && evt.pt_ave < s_Pt_Ave400_cut && trigger320fired) return true;
    if (evt.pt_ave >= s_Pt_Ave400_cut && evt.pt_ave < s_Pt_Ave500_cut && trigger400fired) return true;
    if (evt.pt_ave >= s_Pt_Ave500_cut && trigger500fired) return true;


 return false;
}

bool Selection::DiJet()
{
    assert(event);

    const TClonesArray & js = event->get(h_jets);
    const baconhep::TEventInfo & info = event->get(h_eventInfo);
    baconhep::TEventInfo* eventInfo= new baconhep::TEventInfo(info);
    assert(eventInfo);

    int njets = js.GetEntries();

    // njets >= 2
    if (njets>=2) return true;

 return false;
}
bool Selection::DiJetAdvanced(uhh2::Event& evt)
{
  assert(event);

  const TClonesArray & js = event->get(h_jets);
  const baconhep::TEventInfo & info = event->get(h_eventInfo);
  //   const baconhep::TJet * jet = dynamic_cast<const baconhep::TJet*>(js[0]);
  //   assert(jet);

  baconhep::TEventInfo* eventInfo= new baconhep::TEventInfo(info);
  assert(eventInfo);
  //   const baconhep::TJet * jet = dynamic_cast<const baconhep::TJet*>(js[0]);
  //   assert(jet);

  int njets = js.GetEntries();

  // njets >= 2
  if (njets < 2) return false;

  baconhep::TJet* jet1 = (baconhep::TJet*)js[0];
  baconhep::TJet* jet2 = (baconhep::TJet*)js[1];
  
  // at least one barrel jet
  if((fabs(jet1->eta) >= s_eta_barr) && (fabs(jet2->eta) >= s_eta_barr)) return false; 
  
  // delta phi < 2.7
  double deltaPhi = std::abs(TVector2::Phi_mpi_pi(jet1->phi - jet2->phi));
  if (deltaPhi < s_delta_phi) return false;
  
  // |asymm| < 0.7
  if (fabs((evt.jet2_pt - evt.jet1_pt) / (evt.jet2_pt + evt.jet1_pt)) > s_asymm) return false;

  // p_t,rel < 0.2
  if (njets>2){
    baconhep::TJet* jet3 = (baconhep::TJet*)js[2];
    if ((2*(evt.jet3_pt))/(evt.jet1_pt + evt.jet2_pt) > s_pt_rel) return false;
  }
  
  return true;
}

bool Selection::goodPVertex()
{
    assert(event);

    const TClonesArray & pvs = event->get(h_pv);

    Int_t nvertices = pvs.GetEntries();
    // require in the event that there is at least one reconstructed vertex
    if(nvertices<=0) return false;
    float nPrVer = 0;
    // pick the first (i.e. highest sum pt) verte
    for (int i=0; i<nvertices; i++){
        baconhep::TVertex* vertices = (baconhep::TVertex*)pvs[i];
        // require that the vertex meets certain criteria

        if((vertices->nTracksFit > s_n_PvTracks) && (fabs(vertices->z) < s_n_Pv_z) && (fabs(vertices->y) < s_n_Pv_xy) && (fabs(vertices->x) < s_n_Pv_xy)){
            nPrVer++;
        }
    }

 return true;
}





bool Selection::FullSelection()
{
    return DiJet()&&goodPVertex();
}

Selection::~Selection()
{
}

}
