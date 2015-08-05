// system headers
#include<iostream>
#include <string>
#include<iostream>
#include<fstream>
#include <sys/stat.h>
#include <map>
#include <algorithm>
// ROOT headers
#include<TCollection.h>
#include<TMath.h>
#include<TFile.h>
#include<TRandom3.h>
#include<TSystem.h>
#include<TVector3.h>
#include <vector>

#include "UHH2/BaconJets/include/data_corrections.h"
#include "UHH2/BaconJets/include/constants.h"
#include "UHH2/bacondataformats/interface/TJet.hh"
using namespace std;

namespace uhh2bacon {
    int NPOINTS;
    map<float, map<float, float> > m_corr_l1;
//     map<float, float> m_corr_l1;

    float eta1_jet[82];
    float eta2_jet[82];
    float L1corr[82];

DataCorr::DataCorr(uhh2::Context & ctx) :
    context(ctx),
    event(0)
{
    // open the file
    TString fullname = "/nfs/dust/cms/user/mstoev/CMSSW_7_2_1_patch4/src/UHH2/BaconJets/conf/Summer15_50nsV2_DATA_L1FastJet_AK4PFchs.txt";
    ifstream file(fullname);
    if (file.good()) {
        cout << "INFO: opened " << fullname << endl;
    } else {
        cout << "ERROR: could not open " << fullname << endl;
        abort();
    }

    // loop over the file lines
    string line;
    while ( file.good() ) {

        // read each line and convert it to TString
        getline (file,line);
        TString line_str = line;

        // tokenize it
        TObjArray * tokens = line_str.Tokenize(" ");
        if (tokens -> IsEmpty()) continue;
        TString token_eta1 = ((TObjString*) tokens->At(0))  -> GetString();
        TString token_eta2 = ((TObjString*) tokens->At(1))  -> GetString();
        TString token1     = ((TObjString*) tokens->At(10)) -> GetString();
//             cout << "token_eta1: " << token_eta1 << " token1 " << token1 << endl;

        eta1_jet[NPOINTS] = token_eta1.Atof();
        eta2_jet[NPOINTS] = token_eta2.Atof();
        L1corr[NPOINTS]   = token1.Atof();
        cout << "eta1_jet: " << eta1_jet[NPOINTS]  << "eta2_jet: " << eta2_jet[NPOINTS]<< " L1corr " << L1corr[NPOINTS] << endl;

//         m_corr_l1[eta1_jet[NPOINTS]] = L1corr[NPOINTS];
        m_corr_l1[eta1_jet[NPOINTS]][eta2_jet[NPOINTS]] = L1corr[NPOINTS];

        // increment the point counter
        NPOINTS++;
    }

    cout << "INFO: " << NPOINTS << " entries found in " << fullname << endl;

    file.close();

}

void DataCorr::SetEvent(uhh2::Event& evt)
{
    event = &evt;
    assert(event);

}
float  DataCorr::jet12getL1correction(float Eta1, float Eta2, float & j2L1corr) {
    assert(event);

    for (int i=0; i<NPOINTS; ++i){
        if ((Eta2 > eta1_jet[i]) && (Eta2 < eta2_jet[i])) {
            j2L1corr = m_corr_l1[eta1_jet[i]][eta2_jet[i]];
        }
        if ((Eta1 > eta1_jet[i]) && (Eta1 < eta2_jet[i])) {
            return m_corr_l1[eta1_jet[i]][eta2_jet[i]];
        }
    }

}

float  DataCorr::jet3getL1correction(float Eta) {
    assert(event);

    for (int i=0; i<NPOINTS; ++i){

        if ((Eta > eta1_jet[i]) && (Eta < eta2_jet[i])) {
//           cout << "JET3: eta from ev: "<<Eta<<" eta from txt: "<<eta_jet[i] <<" l1 corr = "<< L1corr[i]<< endl;
//           cout <<" JET3: CORR : "<<m_corr_l1[eta_jet[i]]<<endl;
         return m_corr_l1[eta1_jet[i]][eta2_jet[i]];        }
    }

}

DataCorr::~DataCorr()
{
}

}
