#include<stdlib.h> 
#include<stdio.h>
#include<math.h>
#include<iostream>
#include<time.h>
#include<vector>
#include "TMath.h"
#include "TVector3.h" 
#include "TLorentzVector.h" 
#include "TTree.h" 
#include "TFile.h" 
#include "TRandom3.h"
#include "TComplex.h"
#include "LinkDef.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TMinuit.h"
#include "TStyle.h"
#include "TROOT.h"
#include <TVector3.h>
#include <TLorentzVector.h>
using namespace std;
using namespace ROOT::Math;



void fastalpha(){
  int nevent;
  int first = 0;
  int n; //nstrip and so on
  char namecore[200]="";
  char filename[200];
  const char* core;



  //constant----------------------------
  const double pi = TMath::Pi();
 

  //Plate  
  double error;
  //  double max[10] = {0};
  //  double min[10] = {0};
  int test = 0;
  int accept = 0;
  int accepted = 0;
  double drcut = atoi(getenv("DRCUT"));
  int dlocalx = 0;

  //cout <<"Please Enter the Filename before \".root\""<<endl;
  //cin >> namecore;
  strcpy(namecore,getenv("NAMECORE"));
  core = namecore;
  char *ret;
  char namefordata[200] = "";
  ret = strstr(Form("%s",namecore),"drcut");
  if(ret!=NULL){
    int drcutlength = strlen(Form("%s",ret));
    int namecorelength = strlen(Form("%s",namecore));
    strncat(namefordata,Form("%s",namecore),namecorelength-drcutlength);
    dlocalx = 1;
  }else{
    strcat(namefordata,Form("%s",namecore));
    //drcut = 1000000;
  }

  strcpy(filename,core);
  strcat(filename,"tag.root");
  TFile *file = new TFile(filename);
  TTree *tree = (TTree*)file->Get("tag");
  //Treesetting

  int ndata = 0;
  double x = 0;
  double y = 0;
  double z = 0;
  double mux = 0;
  double muy = 0;
  double muz = 0;
  double mut = 0;
  double mux0 = 0;
  double muy0 = 0;
  double muz0 = 0;
  double mut0 = 0;
  double lx = 0;
  double ly = 0;
  int i = 0;
  int j = 0;
  int k = 0;
  double t = 0;
  int pos = 0;
  float charge = 0;
  int layer = 0;
  int dial = 0;
  int sign = 0;
  double eta = 0;
  double mueta = 0;
  double muentryeta = 0;
  double dr0 = 0;
  int sector=0;
  int band = 0;
 
  tree->SetBranchAddress("ndata",&ndata);
  tree->SetBranchAddress("x",&x);
  tree->SetBranchAddress("y",&y);
  tree->SetBranchAddress("z",&z);
  tree->SetBranchAddress("mux",&mux);
  tree->SetBranchAddress("muy",&muy);
  tree->SetBranchAddress("muz",&muz);
  tree->SetBranchAddress("mut",&mut);
  tree->SetBranchAddress("mux0",&mux0);
  tree->SetBranchAddress("muy0",&muy0);
  tree->SetBranchAddress("muz0",&muz0);
  tree->SetBranchAddress("mut0",&mut0);
  tree->SetBranchAddress("lx",&lx);
  tree->SetBranchAddress("ly",&ly);
  tree->SetBranchAddress("t",&t);
  tree->SetBranchAddress("i",&i);
  tree->SetBranchAddress("j",&j);
  tree->SetBranchAddress("k",&k);
  tree->SetBranchAddress("pos",&pos);
  tree->SetBranchAddress("charge",&charge);
  tree->SetBranchAddress("layer",&layer);
  tree->SetBranchAddress("dial",&dial);
  tree->SetBranchAddress("sign",&sign);
  tree->SetBranchAddress("eta",&eta);
  tree->SetBranchAddress("mueta",&mueta);
  tree->SetBranchAddress("muentryeta",&muentryeta);
  tree->SetBranchAddress("dr0",&dr0);
  tree->SetBranchAddress("sector",&sector);
  tree->SetBranchAddress("band",&band);

  nevent = tree->GetEntries();
  cout<<nevent<<" events"<<endl; 



  n = 0;
  first = 0;


  core = namecore;
  strcpy(filename,core);
  strcat(filename,"fast.root");
  //tree creation
  TFile *filefast = new TFile(filename,"recreate");
  TTree *fast = new TTree("fast","fast");

  //Plate
  int platendata = 0;                                                                       
  double platemueta = 0;
  double platemuentryeta = 0;
  int platelayer[33] = {0};
  double platetime[33] = {0};
  int platedial = 0;
  int platei = 0;
  int platenhit = 0;
  int plateneta = 0;
  int platenstereo = 0;
  double platex[33] = {0};
  double platey[33] = {0};
  double platez[33] = {0};
  double platemux[33] = {0};
  double platemuy[33] = {0};
  double platemuz[33] = {0};
  double platemut[33] = {0};
  double platemux0[33] = {0};
  double platemuy0[33] = {0};
  double platemuz0[33] = {0};
  double platemut0[33] = {0};
  double platelx[33] = {0};
  double plately[33] = {0};
  double platedr0[33] = {0};
  int platesign[33] = {0};
  int plateband = 0;
  int platesector = 0;
  fast->Branch("ndata",&platendata);
  fast->Branch("mueta",&platemueta);
  fast->Branch("muentryeta",&platemuentryeta);
  fast->Branch("x",platex,"x[33]/D");
  fast->Branch("y",platey,"y[33]/D");
  fast->Branch("z",platez,"z[33]/D");
  fast->Branch("mux",platemux,"mux[33]/D");
  fast->Branch("muy",platemuy,"muy[33]/D");
  fast->Branch("muz",platemuz,"muz[33]/D");
  fast->Branch("mut",platemut,"mut[33]/D");
  fast->Branch("mux0",platemux0,"mux0[33]/D");
  fast->Branch("muy0",platemuy0,"muy0[33]/D");
  fast->Branch("muz0",platemuz0,"muz0[33]/D");
  fast->Branch("mut0",platemut0,"mut0[33]/D");
  fast->Branch("lx",platelx,"lx[33]/D");
  fast->Branch("ly",plately,"ly[33]/D");
  fast->Branch("dr0",platedr0,"dr0[33]/D");
  fast->Branch("dial",&platedial);
  fast->Branch("layer",platelayer,"layer[33]/I");
  fast->Branch("time",platetime,"time[33]/D");
  fast->Branch("sign",platesign,"sign[33]/I");
  fast->Branch("i",&platei);
  fast->Branch("nhit",&platenhit);
  fast->Branch("neta",&plateneta);
  fast->Branch("nstereo",&platenstereo);
  fast->Branch("sector",&platesector);
  fast->Branch("band",&plateband);

  int m = 0;
  int f = 0;
  int ok = 0;
  int dialcheck[16]={0};
  int check = 0;
  int bandmemory = 0;

  tree->GetEntry(0);
  while(n < nevent){

    if((n+1)%100000 == 0) cout<<(n+1)/1.0/nevent<<endl;
    ok = 0;
    
    first = 0;
    platedial = 10;//fatal error

    for(int d = 1;d<17;d++){
      dialcheck[d-1]=0;
    }

    while(m == i){

      ok = 1;
      platei = i;
 
      //initiallize
      
      for(int d = 0; d < 33; d++){
	platex[d] = 0;
	platey[d] = 0;
	platez[d] = 0;
	platemux[d] = 0;
	platemuy[d] = 0;
	platemuz[d] = 0;
	platemut[d] = 0;
	platemux0[d] = 0;
	platemuy0[d] = 0;
	platemuz0[d] = 0;
	platemut0[d] = 0;
	platelx[d] = 0;
	plately[d] = 0;
	platedr0[d] = -7000;
	platesign[d] = 0;
	platelayer[d] = 0;
	platetime[d] = 0;
      }
     
      if(first == 0){ 
	platedial = dial;
	first = 1;
	platendata = ndata;
	platemueta = mueta;
	platemuentryeta = muentryeta;
      }

      if(first == 1){
	for(int d = 1;d < 17;d++){
	  if(dial == d && dialcheck[d-1]!=1){
	    dialcheck[d-1] = 1;
	  }
	}
      }

      bandmemory = band;
      plateband = band;
      platesector = sector;
      if(dlocalx==1){
	platesector = 3;
	bandmemory = 0;
	plateband = 0;
      }

      while(m==i && band == bandmemory){
	for(int d = 0; d < 33; d++){
	  
	  if(platelayer[d]==0){
	    
	    if(layer+16 == d){
	      platex[d] = x;
	      platey[d] = y;
	      platez[d] = z;
	      platemux[d] = mux;
	      platemuy[d] = muy;
	      platemuz[d] = muz;
	      platemut[d] = mut;
	      platemux0[d] = mux0;
	      platemuy0[d] = muy0;
	      platemuz0[d] = muz0;
	      platemut0[d] = mut0;
	      platelx[d] = lx;
	      plately[d] = ly;
	      platedr0[d] = dr0;
	      platesign[d] = sign;
	      platelayer[d] = 1;
	      platetime[d] = t;
	    }
	  }	
		  
	}
      
	n = n+1;
	if(n < nevent){
	  tree->GetEntry(n);
	}
	if(n == nevent){
	  m = 0;
	  band = -1;
	}
      }

      if(ok==1){
	platenhit = 0;
	plateneta = 0;
	platenstereo = 0;
	
	for(int d=0;d<33;d++){
	  if(platelayer[d]==1){
	    platenhit += platelayer[d];
	    if(platesign[d]!=0){
	      platenstereo += platelayer[d];
	    }else{
	      plateneta += platelayer[d];
	    }
	  }
	}
	
	test = 0;
	check = 0;
	for(int d=1;d<17;d++){
	  test += dialcheck[d-1];
	  if(dialcheck[d-1]!=0){check += d;}
	}     
	if(test == 1){
	}else{
	  if(test == 2){
	    if(dialcheck[0]==1 && dialcheck[15]==1){
	      platedial = 17;
	    }else{
	      platedial = (check-1)/2+16;
	      if(dialcheck[(check+1)/2-1]==0 || dialcheck[(check-1)/2-1]==0){
		platedial = 0;
	      }
	    }
	  }else{
	    platedial = 0;
	  }
	}

	fast->Fill();
	f = f+1;
	accept = 0;
	for(int d=0;d<33;d++){
	  if(platedr0[d]!=-7000&&fabs(platedr0[d])>drcut) accept = -1;
	}
	if(accept==0) accepted += 1;
	
      }
      
    }//while mi last


    m++;    
  }
  
  
  
  fast->Write();
  filefast->Close();
  
  cout <<"nfast"<< f <<endl;
  cout <<"valid events"<< accepted <<endl;
  
}
