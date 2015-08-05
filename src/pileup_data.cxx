// system headers
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

#include "UHH2/BaconJets/include/pileup_data.h"
#include "UHH2/BaconJets/include/constants.h"
#include "UHH2/bacondataformats/interface/TJet.hh"
using namespace std;

namespace uhh2bacon {

map<int, map<int, double> > m_PU;

      double PU ;

      int run ;
      int ls;
PileupData::PileupData(uhh2::Context & ctx) :
    context(ctx),
    event(0)
{

  string line;
  ifstream file("pileup_new_json_low_pu.txt");

  if (file.is_open()){
//     cout << "ok" << endl;

    //read first line
    getline(file, line);

    //loop over lines in file
    while ( getline(file,line) ){

      string str, run_str, ls_str;
      int delim_pos;

      //loop over strings in line
      for (int string_num=0; (delim_pos = line.find(",")) != -1; string_num++){

        str = line.substr(0, delim_pos);
        line.erase(0, delim_pos + 1);

        if (string_num == 0)  //first string holds run number
          run_str = str.substr(0, str.find(":"));

        else if (string_num == 1) //second string has ls
          ls_str = str.substr(0, str.find(":"));
      }
      //last part of line  has PU info
       PU = stod( line );

       run = stoi( run_str );
       ls = stoi( ls_str );
//         if ((run == Run) && (ls == Ls)){
      m_PU[run][ls] = PU;
//     cout << "in class: run : "<< run<<" ls: "<<ls<<" PU : "<< PU<<endl;

//         }
    }
    file.close();
  }
  else
    cout << "Unable to open file" << endl; 

}

void PileupData::SetEvent(uhh2::Event& evt)
{
    event = &evt;
    assert(event);

}

float  PileupData::getDataPU(int Run, int Ls) {

    assert(event);

   // if (ev_weighting_factor!=0) return      ev_weighting_factor;
      //  if ((run == Run) && (ls == Ls)){
     return m_PU[Run][Ls];
     //   }
}

PileupData::~PileupData()
{
}

}
