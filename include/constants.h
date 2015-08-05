#ifndef  CONSTANTS_H
#define  CONSTANTS_H

/** \brief Binning **/
// static std::vector<double>   eta_range  = {0, 0.261, 0.522, 0.763, 0.957, 1.131, 1.305, 1.479, 1.93, 2.322, 2.411, 2.5, 2.853, 2.964, 3.139, 3.489, 5.191};
// static std::vector<double>   pt_range   = {66, 107, 191, 240, 306, 379, 468, 900};
// static std::vector<double>   alpha_range= {0., 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25};


static std::vector<double>   eta_range  = {0, 0.25, 0.5, 0.75, 1.0, 1.3, 1.6, 1.9, 2.2, 2.5, 2.7, 3.0, 3.2, 5.0};//{0, 1.3, 2.5, 3.0, 5.0};
static std::vector<double>   eta_range_mikko_barrel = {0, 0.8, 1.3};
//static std::vector<double>   pt_range   = {55, 76, 93, 172, 232, 300, 366, 453, 558};
static std::vector<double>   pt_range   = {43, 80, 88, 135, 223, 290, 365, 448, 561};
// static std::vector<double>   pt_range   = {57, 80, 100, 166, 231, 290, 365, 478, 556};

static std::vector<double>   alpha_range= {0., 0.1, 0.125, 0.15, 0.175, 0.2, 0.225, 0.25};

/** \brief Dijet event selection **/
// barrel region (|eta| < 1.3)
static float s_eta_barr = 1.3;
// two back-to-back leading jets (delta_phi(j1,j2) = min(|phi1 - phi2|, 2PI - |phi2 - phi1|) > 2.7)
static float s_delta_phi = 2.7;
// cut on the asymmetry for events with two jets  |(j2->pt - j1->pt /(j2->pt + j1->pt)| < 0.70
static float s_asymm = 0.7;
// relative third jet fraction pt_rel = 2*j3_pt/(j1_pt + j2_pt) < 0.2
static float s_pt_rel = 0.4;

/** \brief good Primary Vertex reconstruction **/
// more than four tracks
static float s_n_PvTracks = 4;
// PV is located within 24cm in z vertex
static float s_n_Pv_z = 24.0;
// PV is located within 2cm in xy direction from the nominal interaction point
 static float s_n_Pv_xy = 2.0;


/** \brief The trigger thresholds of pt_ave **/
/// Used for the 2012 data
// static float s_Pt_Ave40_cut   = 66; //Hennings: 62
// static float s_Pt_Ave80_cut   = 107;//Hennings: 107
// static float s_Pt_Ave140_cut  = 191;//Hennings: 175
// static float s_Pt_Ave200_cut  = 240;//Hennings: 242
// static float s_Pt_Ave260_cut  = 306;//Hennings: 310
// static float s_Pt_Ave320_cut  = 379;//Hennings: 379
// static float s_Pt_Ave400_cut  = 468;//Hennings: 467
/// Used for the 2015 data
static float s_Pt_Ave40_cut   = 43; //Hennings: 62
static float s_Pt_Ave60_cut   = 80; //Hennings: 62
static float s_Pt_Ave80_cut   = 88;//Hennings: 107
static float s_Pt_Ave140_cut  = 135;//Hennings: 175
static float s_Pt_Ave200_cut  = 223;//Hennings: 242
static float s_Pt_Ave260_cut  = 290;//Hennings: 310
static float s_Pt_Ave320_cut  = 365;//Hennings: 379
static float s_Pt_Ave400_cut  = 448;//Hennings: 467
static float s_Pt_Ave500_cut  = 561;//Hennings: 467

/** \brief Jet Resolution Smearering **/
// doing the matching from GEN to RECO
static float s_delta_R   = 0.3; 
//constant numberstaken from https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
// from 8TeV JER measurement.
constexpr const size_t n = 7;
static float eta_hi[n]    = {0.5, 1.1, 1.7, 2.3, 2.8, 3.2, 5.0};
static float c_nominal[n] = {1.079, 1.099, 1.121, 1.208, 1.254, 1.395, 1.056};
static float c_up[n]      = {1.105, 1.127, 1.150, 1.254, 1.316, 1.458, 1.247};
static float c_down[n]    = {1.053, 1.071, 1.092, 1.162, 1.192, 1.332, 0.865};


/** \brief Dijet event weighting **/
// in each pt-region defined by the trigger thresholds # of events in MC and in data should equal
static float scale_factor1[] = {2.14563, 4.27422,7.57198,  10.7035,17.899, 15.6095, 4.64629};
static float scale_factor2[] = {242.617, 3514.77, 49549.8, 180540, 693177, 1.83449e+06, 1.63241e+07};



#endif
