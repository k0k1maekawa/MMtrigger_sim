//first program, so this is only for the deployment of vector<vector< > >
#include<stdlib.h> 
#include<stdio.h>
#include<math.h>
#include<iostream>
#include<time.h>
#include<vector>
#include<algorithm>
#include<functional>
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

double intersecX(double a1, double b1, double c1, double a2, double b2, double c2, double z) {
  return a1+(z-c1)/(c2-c1)*(a2-a1);
}

double intersecY(double a1, double b1, double c1, double a2, double b2, double c2, double z) {
  return b1+(z-c1)/(c2-c1)*(b2-b1);
}

double distance1(double a, double b, double c, double d1, double d2, double d3, double x, double y, double z) {
  return ((a-x)*(a-x)+(b-y)*(b-y)+(c-z)*(c-z))-((a-x)*d1+(b-y)*d2+(c-z)*d3)*((a-x)*d1+(b-y)*d2+(c-z)*d3);
}

int sign(double x){
  int sgn=0;
  if(x>0) sgn=1;
  if(x<0) sgn=-1;
  return sgn;
}

void readalpha(){
  string MML="MML";
  int nevent;
  int first; //first or not
  int n; //nstrip and so on
  char namecore[300]="";
  const double lowPtedge = 80000;//cut off low pt muon (signal muon remain) The unit is MeV.
  const int Ndata = atoi(getenv("NDATA"));
  const int shortcut = 0;//for debug

  //constant----------------------------
  const double pi = TMath::Pi();
  /*  const double alpha = 0.007297;
  const double E_cms = 5.29; 
  const double M_tau = 1.777; 
  const double M_rho = 0.7754;
  const double M_pi0 = 0.13497;
  const double M_pi = 0.13957;
  const double M_mu = 0.10565836668;
  const double M_e = 0.000510999;
  const double G = 0.0000116637;
  const double e = 1.60217657*pow(10,-19);
  const double b = -0.091;
  const double cosc = 0.97;*/

  //Plate  
  double error;
  double max[10] = {0};
  double min[10] = {0};
  double test = 0;
  double layerz[33] = {0};
  int layersign[33] = {0};
  int intlayer=0;
  int intphi=0;
  int inteta=0;
  int intposition=0;
  int intvmm=0;
  //layer:phi:module(eta):id

  //cout <<"Please Enter the Filename before \".root\""<<endl;
  //cin >> namecore;
  strcpy(namecore,getenv("NAMECORE"));
  cout<<namecore<<".root"<<endl;
  char *ret;
  char namefordata[300] = "";
  ret = strstr(Form("%s",namecore),"drcut");
  if(ret!=NULL){
    int drcutlength = strlen(Form("%s",ret));
    int namecorelength = strlen(Form("%s",namecore));
    strncat(namefordata,Form("%s",namecore),namecorelength-drcutlength);
  }else{
    strcat(namefordata,Form("%s",namecore));
  }
  TFile *file[Ndata];
  file[0] = new TFile(Form("%s.root",namefordata));
  for(int d=1;d<Ndata;d++){
    file[d] = new TFile(Form("%s_%d.root",namefordata,d));
  }
  TTree *tree = (TTree*)file[0]->Get("NSWHitsTree");
  //Treesetting
  vector<double> *vertex_X=0;
  vector<double> *vertex_Y=0;
  vector<double> *vertex_Z=0;
  vector<double> *hitspos_X=0;
  vector<double> *hitspos_Y=0;
  vector<double> *hitspos_Z=0;
  vector<double> *hitstime=0;
  vector<double> *local_hitspos_X=0;
  vector<double> *local_hitspos_Y=0;
  vector<double> *local_hitspos_Z=0;
  vector<double> *truth_Eta=0;
  vector<double> *truth_Phi=0;
  vector<double> *truth_Pt=0;
  vector<int>    *truth_Pdg =0;
  vector<vector<double> > *Digits_Gstrip_X=0;
  vector<vector<double> > *Digits_Gstrip_Y=0;
  vector<vector<double> > *Digits_Gstrip_Z=0;
  vector<vector<double> > *Digits_Lstrip_X=0;
  vector<vector<double> > *Digits_Lstrip_Y=0;


  //  vector<vector<double> > *Digits_stripPosition=0;                             
  vector<vector<double> > *Digits_time=0;
  vector<vector<float> > *Digits_charge=0;
  vector<double> *Digits_stripForTrigger=0;
  vector<double> *Digits_stripTimeForTrigger=0;
  vector<double> *truth_MuEntry_Eta=0;
  vector<double> *truth_MuEntry_Phi=0;
  vector<int> *channel=0;
  //new variable
  vector<vector<double> >  *Digits_stripResponse_time=0;
  vector<vector<int> >  *Digits_stripResponse_stripPosition=0;
  vector<vector<double> >  *Digits_stripresponse_stripGposX=0;
  vector<vector<double> >  *Digits_stripResponse_stripGposY=0;
  vector<vector<double> >  *Digits_stripResponse_stripGposZ=0;
  vector<int> *hitspdg=0;
  vector<string> *Digits_Sector=0;
  vector<int> *Digits_stationEta=0;
  vector<int> *Digits_multiplet=0;
  vector<int> *Digits_gasgap=0;
  vector<int> *Digits_stationPhi=0;
  vector<string> *Hits_Sector=0;
  vector<int> *Hits_side=0;
  vector<int> *Hits_multilayer=0;
  vector<int> *Hits_layer=0;
  vector<int> *Hits_stationPhi=0;
  vector<double> *Hits_kineticEnergy=0;
  vector<double> *Hits_DirectionZ=0;
  vector<vector<int> > *Digits_position=0;
  vector<vector<int> > *Trigger_position=0;
  vector<vector<float> > *Trigger_time=0;
  vector<vector<float> > *Trigger_charge=0;
  vector<vector<int> > *Trigger_VMM=0;
  vector<vector<int> > *Trigger_MMFE=0;
 



  tree->SetBranchAddress("TruthVertex_X",&vertex_X);
  tree->SetBranchAddress("TruthVertex_Y",&vertex_Y);
  tree->SetBranchAddress("TruthVertex_Z",&vertex_Z);
  tree->SetBranchAddress("Hits_MM_hitGlobalPositionX",&hitspos_X);
  tree->SetBranchAddress("Hits_MM_hitGlobalPositionY",&hitspos_Y);
  tree->SetBranchAddress("Hits_MM_hitGlobalPositionZ",&hitspos_Z);
  tree->SetBranchAddress("Hits_MM_globalTime",&hitstime);
  tree->SetBranchAddress("Hits_MM_hitLocalPositionX",&local_hitspos_X);
  tree->SetBranchAddress("Hits_MM_hitLocalPositionY",&local_hitspos_Y);
  tree->SetBranchAddress("Hits_MM_hitLocalPositionZ",&local_hitspos_Z);
  tree->SetBranchAddress("TruthParticle_Eta",&truth_Eta);
  tree->SetBranchAddress("TruthParticle_Pt" ,&truth_Pt);
  tree->SetBranchAddress("TruthParticle_Phi",&truth_Phi);
  tree->SetBranchAddress("TruthParticle_Pdg",&truth_Pdg);
  tree->SetBranchAddress("Hits_MM_particleEncoding" ,&hitspdg); //hitspdg           
  tree->SetBranchAddress("Digits_MM_stripGposX",&Digits_Gstrip_X);
  tree->SetBranchAddress("Digits_MM_stripGposY",&Digits_Gstrip_Y);
  tree->SetBranchAddress("Digits_MM_stripGposZ",&Digits_Gstrip_Z);
  tree->SetBranchAddress("Digits_MM_stripLposX",&Digits_Lstrip_X);
  tree->SetBranchAddress("Digits_MM_stripLposY",&Digits_Lstrip_Y);

  tree->SetBranchAddress("Digits_MM_time",&Digits_time);
  tree->SetBranchAddress("Digits_MM_charge",&Digits_charge);

  //new
  
  tree->SetBranchAddress("Digits_MM_stripResponse_time",&Digits_stripResponse_time);
  //  tree->SetBranchAddress("Digits_MM_stripResponse_charge",&Digits_stripResponse_charge);
  tree->SetBranchAddress("Digits_MM_stripResponse_stripPosition",&Digits_stripResponse_stripPosition);
  tree->SetBranchAddress("Digits_MM_stripresponse_stripGposX",&Digits_stripresponse_stripGposX);
  tree->SetBranchAddress("Digits_MM_stripResponse_stripGposY",&Digits_stripResponse_stripGposY);
  tree->SetBranchAddress("Digits_MM_stripResponse_stripGposZ",&Digits_stripResponse_stripGposZ);
  //  tree->SetBranchAddress("Digits_MM_stripForTrigger",&Digits_stripForTrigger);
  //  tree->SetBranchAddress("Digits_MM_stripTimeForTrigger",&Digits_stripTimeForTrigger);
  tree->SetBranchAddress("MuEntry_Particle_Eta",&truth_MuEntry_Eta);
  tree->SetBranchAddress("MuEntry_Particle_Phi",&truth_MuEntry_Phi);
  tree->SetBranchAddress("Digits_MM_channel",&channel);
  tree->SetBranchAddress("Digits_MM_stationName",&Digits_Sector);
  tree->SetBranchAddress("Digits_MM_stationEta",&Digits_stationEta);
  tree->SetBranchAddress("Digits_MM_multiplet",&Digits_multiplet);
  tree->SetBranchAddress("Digits_MM_gas_gap",&Digits_gasgap);
  tree->SetBranchAddress("Digits_MM_stationPhi",&Digits_stationPhi);
  tree->SetBranchAddress("Hits_MM_sim_stationName",&Hits_Sector);
  tree->SetBranchAddress("Hits_MM_sim_side",&Hits_side);
  tree->SetBranchAddress("Hits_MM_sim_multilayer",&Hits_multilayer);
  tree->SetBranchAddress("Hits_MM_sim_layer",&Hits_layer);
  tree->SetBranchAddress("Hits_MM_sim_stationPhi",&Hits_stationPhi);
  tree->SetBranchAddress("Hits_MM_kineticEnergy",&Hits_kineticEnergy);
  tree->SetBranchAddress("Hits_MM_hitGlobalDirectionZ",&Hits_DirectionZ);
  tree->SetBranchAddress("Digits_MM_stripPosition",&Digits_position);
  tree->SetBranchAddress("Digits_MM_position_trigger",&Trigger_position);
  tree->SetBranchAddress("Digits_MM_time_trigger",&Trigger_time);
  tree->SetBranchAddress("Digits_MM_charge_trigger",&Trigger_charge);
  tree->SetBranchAddress("Digits_MM_VMM_id_trigger",&Trigger_VMM);
  tree->SetBranchAddress("Digits_MM_MMFE_VMM_id_trigger",&Trigger_MMFE);

  nevent = tree->GetEntries();
    
  first = 0;
  n = 0;

  if(namecore[0]=='M'){
    for(int i = 0; i < nevent; i++){
      if((i+1)%100000 == 0) cout<<i+1<<endl;
      tree->GetEntry(i);
      
      for(int j = 0; j < Digits_Gstrip_Z->size(); j++){
	for(int k=0; k<Digits_Gstrip_Z->at(j).size(); k++){
	  if(Digits_Gstrip_Z->at(j).at(k)!=0){
	    
	    intlayer = (2*(Digits_stationEta->at(j)>0)-1)*(8*(Digits_Sector->at(j)==MML)+4*(Digits_multiplet->at(j)-1)+Digits_gasgap->at(j))+16;
	    
	    layerz[intlayer] = Digits_Gstrip_Z->at(j).at(k);
	    layersign[intlayer] = sign(Digits_Lstrip_Y->at(j).at(k));
	    
	  }
	}
      }
      if(layerz[17]!=0&&layerz[18]!=0&&layerz[19]!=0&&layerz[20]!=0&&layerz[21]!=0&&layerz[22]!=0&&layerz[23]!=0&&layerz[24]!=0&&layerz[25]!=0&&layerz[26]!=0&&layerz[27]!=0&&layerz[28]!=0&&layerz[29]!=0&&layerz[30]!=0&&layerz[31]!=0&&layerz[32]!=0)break;
    }
  }else{
    for(int i = 0; i < nevent; i++){
      if((i+1)%100000 == 0) cout<<i+1<<endl;
      tree->GetEntry(i);
      
      for(int j = 0; j < Digits_Gstrip_Z->size(); j++){
	for(int k=0; k<Digits_Gstrip_Z->at(j).size(); k++){
	  if(Digits_Gstrip_Z->at(j).at(k)!=0){
	    
	    intlayer = (2*(Digits_stationEta->at(j)>0)-1)*(8*(Digits_Sector->at(j)==MML)+4*(Digits_multiplet->at(j)-1)+(5-Digits_gasgap->at(j)))+16;
	    
	    layerz[intlayer] = Digits_Gstrip_Z->at(j).at(k);
	    layersign[intlayer] = -sign(intlayer-16)*sign(Digits_Lstrip_Y->at(j).at(k));
	    
	  }
	}
      }
    }
  }
  
  for(int d=0;d<33;d++){
    cout<<"layer"<<d-16<<"\t"<<layerz[d]<<"\t"<<layersign[d]<<endl;
  }
  
  
  cout<<nevent<<" events"<<endl;
  first = 0; //reset
  n = 0;
  TFile *filexyz = new TFile(Form("%sxyz.root",namefordata),"recreate");
  TTree *xyz = new TTree("xyz","xyz");
  int ndata = 0;
  int nmu[33]={0};
  double mux[33]={0};
  double muy[33]={0};
  double muz[33]={0};
  double mut[33]={0};
  double mur[33]={0};
  double mufirst[33]={0};
  int mudial[33]={0};
  double mux0[33]={0};
  double muy0[33]={0};
  double muz0[33]={0};
  double mut0[33]={0};
  double mur0[33]={0};
  int mudial0[33]={0};
  double phifraction = 0;
  double hit_phimemory = 0;
  double hit_thetamemory = 0;
  double hit_xmemory = 0;
  double hit_ymemory = 0;
  double hit_zmemory = 0;
  int nmuhit = 0;
  int nmuhits = 0;
  vector<double> hitmux(100);
  vector<double> hitmuy(100);
  vector<double> hitmuz(100);
  vector<double> hitmut(100);

  double lphi;
  double l1;
  double l2;
  double l3;

  int bbb=0;

  double strip_phifraction=0;
  double hit_phifraction=0;
  double hit0_phifraction=0;

  ///////j
  double plate_mu_eta;
  double plate_muentry_eta;

  ///////k
  double platex;
  double platey;
  double platez;
  double plater;
  double platephi;
  double platemux;
  double platemuy;
  double platemuz;
  double platemut;
  //double platemupos;
  double platemur;
  int platemudial;
  double platemux0;
  double platemuy0;
  double platemuz0;
  double platemut0;
  double platemur0;
  int platemudial0;
  double platelx;
  double plately;
  //  double platerx;
  //  double platery;
  //  double platerz;
  int f1;//first integer
  int f2;//second integer
  int f3;//third integer
  double platet;
  //double platert;
  int platepos;
  int platevmm;
  int platemmfe;
  //  int platerpos;
  float platecharge;
  int platelayer;
  int plateupdown;
  int platedial;
  int platesign;
  double plateeta;
  double platedr;
  double platedr0;
  double platedly;
  double platepror;
  //width
  double zwidth = 6;
  double phiwidth = 0.18;
  int ons = 0;
  int on[3] = {0}; 
  int non[3] = {0};
  xyz->Branch("ndata",&ndata);
  xyz->Branch("mueta",&plate_mu_eta);
  xyz->Branch("muentryeta",&plate_muentry_eta);
  xyz->Branch("eta",&plateeta);
  xyz->Branch("x",&platex);
  xyz->Branch("y",&platey);
  xyz->Branch("z",&platez);
  xyz->Branch("mux",&platemux);
  xyz->Branch("muy",&platemuy);
  xyz->Branch("muz",&platemuz);
  xyz->Branch("mut",&platemut);
  xyz->Branch("mudial",&platemudial);
  xyz->Branch("mur",&platemur);
  xyz->Branch("mux0",&platemux0);
  xyz->Branch("muy0",&platemuy0);
  xyz->Branch("muz0",&platemuz0);
  xyz->Branch("mut0",&platemut0);
  xyz->Branch("mur0",&platemur0);
  xyz->Branch("mudial0",&platemudial0);
  xyz->Branch("lx",&platelx);
  xyz->Branch("ly",&plately);
  xyz->Branch("i",&f1);
  xyz->Branch("j",&f2);
  xyz->Branch("k",&f3);
  xyz->Branch("t",&platet);
  xyz->Branch("pos",&platepos);
  xyz->Branch("vmm",&platevmm);
  xyz->Branch("charge",&platecharge);
  xyz->Branch("layer",&platelayer);
  xyz->Branch("dial",&platedial);
  xyz->Branch("sign",&platesign);
  xyz->Branch("dr",&platedr);
  xyz->Branch("pror",&platepror);
  xyz->Branch("dr0",&platedr0);
  xyz->Branch("dly",&platedly);
  xyz->Branch("hit_phifraction",&hit_phifraction);
  xyz->Branch("hit0_phifraction",&hit0_phifraction);
  xyz->Branch("strip_phifraction",&strip_phifraction);
  xyz->Branch("vmm",&platevmm);
  xyz->Branch("mmfe",&platemmfe);
  xyz->Branch("updown",&plateupdown);

  f1 = 0;
  for(int nd=0;nd<Ndata;nd++){
    ndata = nd;
    tree = (TTree*)file[ndata]->Get("NSWHitsTree");
    nevent = tree->GetEntries();    
    cout<<"DATA"<<ndata<<" "<<nevent<<"event"<<endl;
    tree->SetBranchAddress("TruthVertex_X",&vertex_X);
    tree->SetBranchAddress("TruthVertex_Y",&vertex_Y);
    tree->SetBranchAddress("TruthVertex_Z",&vertex_Z);
    tree->SetBranchAddress("Hits_MM_hitGlobalPositionX",&hitspos_X);
    tree->SetBranchAddress("Hits_MM_hitGlobalPositionY",&hitspos_Y);
    tree->SetBranchAddress("Hits_MM_hitGlobalPositionZ",&hitspos_Z);
    tree->SetBranchAddress("Hits_MM_globalTime",&hitstime);
    tree->SetBranchAddress("Hits_MM_hitLocalPositionX",&local_hitspos_X);
    tree->SetBranchAddress("Hits_MM_hitLocalPositionY",&local_hitspos_Y);
    tree->SetBranchAddress("Hits_MM_hitLocalPositionZ",&local_hitspos_Z);
    tree->SetBranchAddress("TruthParticle_Eta",&truth_Eta);
    tree->SetBranchAddress("TruthParticle_Pt" ,&truth_Pt);
    tree->SetBranchAddress("TruthParticle_Phi",&truth_Phi);
    tree->SetBranchAddress("TruthParticle_Pdg",&truth_Pdg);
    tree->SetBranchAddress("Hits_MM_particleEncoding" ,&hitspdg); //hitspdg                         
    tree->SetBranchAddress("Digits_MM_stripGposX",&Digits_Gstrip_X);
    tree->SetBranchAddress("Digits_MM_stripGposY",&Digits_Gstrip_Y);
    tree->SetBranchAddress("Digits_MM_stripGposZ",&Digits_Gstrip_Z);
    tree->SetBranchAddress("Digits_MM_stripLposX",&Digits_Lstrip_X);
    tree->SetBranchAddress("Digits_MM_stripLposY",&Digits_Lstrip_Y);

    tree->SetBranchAddress("Digits_MM_time",&Digits_time);
    tree->SetBranchAddress("Digits_MM_charge",&Digits_charge);

    //new                                                                                                 

    tree->SetBranchAddress("Digits_MM_stripResponse_time",&Digits_stripResponse_time);
    //  tree->SetBranchAddress("Digits_MM_stripResponse_charge",&Digits_stripResponse_charge);            
    tree->SetBranchAddress("Digits_MM_stripResponse_stripPosition",&Digits_stripResponse_stripPosition);
    tree->SetBranchAddress("Digits_MM_stripresponse_stripGposX",&Digits_stripresponse_stripGposX);
    tree->SetBranchAddress("Digits_MM_stripResponse_stripGposY",&Digits_stripResponse_stripGposY);
    tree->SetBranchAddress("Digits_MM_stripResponse_stripGposZ",&Digits_stripResponse_stripGposZ);
    //  tree->SetBranchAddress("Digits_MM_stripForTrigger",&Digits_stripForTrigger);                      
    //  tree->SetBranchAddress("Digits_MM_stripTimeForTrigger",&Digits_stripTimeForTrigger);              
    tree->SetBranchAddress("MuEntry_Particle_Eta",&truth_MuEntry_Eta);
    tree->SetBranchAddress("MuEntry_Particle_Phi",&truth_MuEntry_Phi);
    tree->SetBranchAddress("Digits_MM_channel",&channel);
    tree->SetBranchAddress("Digits_MM_stationName",&Digits_Sector);
    tree->SetBranchAddress("Digits_MM_stationEta",&Digits_stationEta);
    tree->SetBranchAddress("Digits_MM_multiplet",&Digits_multiplet);
    tree->SetBranchAddress("Digits_MM_gas_gap",&Digits_gasgap);
    tree->SetBranchAddress("Digits_MM_stationPhi",&Digits_stationPhi);
    tree->SetBranchAddress("Hits_MM_sim_stationName",&Hits_Sector);
    tree->SetBranchAddress("Hits_MM_sim_side",&Hits_side);
    tree->SetBranchAddress("Hits_MM_sim_multilayer",&Hits_multilayer);
    tree->SetBranchAddress("Hits_MM_sim_layer",&Hits_layer);
    tree->SetBranchAddress("Hits_MM_sim_stationPhi",&Hits_stationPhi);
    tree->SetBranchAddress("Hits_MM_kineticEnergy",&Hits_kineticEnergy);
    tree->SetBranchAddress("Hits_MM_hitGlobalDirectionZ",&Hits_DirectionZ);
    tree->SetBranchAddress("Digits_MM_stripPosition",&Digits_position);
    tree->SetBranchAddress("Digits_MM_position_trigger",&Trigger_position);
    tree->SetBranchAddress("Digits_MM_time_trigger",&Trigger_time);
    tree->SetBranchAddress("Digits_MM_charge_trigger",&Trigger_charge);
    tree->SetBranchAddress("Digits_MM_VMM_id_trigger",&Trigger_VMM);
    tree->SetBranchAddress("Digits_MM_MMFE_VMM_id_trigger",&Trigger_MMFE);
    if(shortcut==1)cout<<"set End"<<endl;

    for(int i = 0; i < nevent; i++){
      if(shortcut == 1 && i==0) i = nevent - 1;
      if((i+1)%100000 == 0) cout<<i+1<<endl;
    tree->GetEntry(i);
    //cout<<"Entry"<<i<<endl;

    ons = 0;
    for(int c=0;c<3;c++){
      on[c]=0;
    }
    
    if(truth_Eta->size() > 0){
      if(truth_MuEntry_Eta->size() > 0){
	plate_mu_eta = truth_Eta->at(0);
	plate_muentry_eta = truth_MuEntry_Eta->at(0);
      }else{
	plate_muentry_eta = 0;
      }
    }else{
      plate_mu_eta = 0;
    }

    for(int d=0;d<33;d++){
      nmu[d] = 0;
      mux[d] = 0;
      muy[d] = 0;
      muz[d] = 0;
      mut[d] = 0;
      mur[d] = 0;
      mufirst[d] = 0;
      mudial[d] = 0;
      mur0[d] = 0;
      mudial0[d] = 0;
    }

    hit_phimemory = 0;
    hit_thetamemory = 0;
    hit_xmemory = 0;
    hit_ymemory = 0;
    hit_zmemory = 0;
    
    hitmux.erase(hitmux.begin(),hitmux.end());
    hitmuy.erase(hitmuy.begin(),hitmuy.end());
    hitmuz.erase(hitmuz.begin(),hitmuz.end());
    hitmut.erase(hitmut.begin(),hitmut.end());


    for(int c=0;c<hitspos_X->size();c++){
      if(hitspdg->at(c)==13&&hitstime->at(c)>22&&hitstime->at(c)<32&&sign(plate_mu_eta)==sign(hitspos_Z->at(c))){
        for(int d=0;d<33;d++){
          if(fabs(hitspos_Z->at(c)-layerz[d])<zwidth && lowPtedge < Hits_kineticEnergy->at(c)*sqrt(1-Hits_DirectionZ->at(c)*Hits_DirectionZ->at(c))){
	    hitmux.push_back(hitspos_X->at(c));
	    hitmuy.push_back(hitspos_Y->at(c));
	    hitmuz.push_back(hitspos_Z->at(c));
	    hitmut.push_back(hitstime->at(c));
	    hit_xmemory += hitspos_X->at(c);
	    hit_ymemory += hitspos_Y->at(c);
	    hit_zmemory += hitspos_Z->at(c);

	    mux[d] += hitspos_X->at(c);
	    muy[d] += hitspos_Y->at(c);
	    muz[d] += hitspos_Z->at(c);
	    if(mufirst[d]==0){
	      mut[d] = hitstime->at(c);
	      mufirst[d] = 1;
	    }
	    if(hitstime->at(c)<mut[d]){
	      mut[d] = hitstime->at(c);
	    }
	    nmu[d] += 1;
	    nmuhit += 1;
	    nmuhits += 1;
	  }
	}
      }
    }

    
    for(int d=0;d<33;d++){
      if(nmu[d]!=0){
	mux[d] /= nmu[d];
	muy[d] /= nmu[d];
	muz[d] /= nmu[d];
	lphi = atan2(muy[d],mux[d]);
	if(lphi<0) lphi = lphi+2*pi;
	if(abs(d-16)>0&&abs(d-16)<9){
	  mudial[d] = round(lphi/pi*4.0+0.5)*2;
	}
	if(abs(d-16)>8&&abs(d-16)<17){
	  mudial[d] = round(lphi/pi*4.0)*2+1; 
	}
	if(mudial[d]==17) mudial[d]=1;
	mur[d] = hypot(mux[d],muy[d])*cos(atan2(muy[d],mux[d])-(mudial[d]-1)*pi/8.0);
      }
    }
    
    if(nmuhit!=0){
      hit_xmemory /= nmuhit;
      hit_ymemory /= nmuhit;
      hit_zmemory /= nmuhit;
    }
    
    hit_phimemory = atan2(hit_ymemory,hit_xmemory);
    hit_thetamemory = atan2(sqrt(hit_xmemory*hit_xmemory+hit_ymemory*hit_ymemory),hit_zmemory);


    for(int d=0;d<33;d++){
      if(hitmux.size()>1){
	mux0[d] = intersecX(hitmux.at(0),hitmuy.at(0),hitmuz.at(0),hitmux.at(hitmux.size()-1),hitmuy.at(hitmuy.size()-1),hitmuz.at(hitmuz.size()-1),layerz[d]);
	muy0[d] = intersecY(hitmux.at(0),hitmuy.at(0),hitmuz.at(0),hitmux.at(hitmux.size()-1),hitmuy.at(hitmuy.size()-1),hitmuz.at(hitmuz.size()-1),layerz[d]);
	muz0[d] = layerz[d];
	mut0[d] = intersecX(hitmut.at(0),hitmuy.at(0),hitmuz.at(0),hitmut.at(hitmux.size()-1),hitmuy.at(hitmuy.size()-1),hitmuz.at(hitmuz.size()-1),layerz[d]);
	lphi = atan2(muy0[d],mux0[d]);
        if(lphi<0) lphi = lphi+2*pi;
        if(abs(d-16)>0&&abs(d-16)<9){
          mudial0[d] = round(lphi/pi*4.0+0.5)*2;
	}
        if(abs(d-16)>8&&abs(d-16)<17){
          mudial0[d] = round(lphi/pi*4.0)*2+1;
        }
        if(mudial0[d]==17) mudial0[d]=1;
        mur0[d] = hypot(mux0[d],muy0[d])*cos(atan2(muy0[d],mux0[d])-(mudial0[d]-1)*pi/8.0);
      }else{
        mux0[d] = vertex_X->at(0)+(layerz[d]-vertex_Z->at(0))*tan(2.0*atan(exp(-truth_Eta->at(0))))*cos(truth_Phi->at(0));
        muy0[d] = vertex_Y->at(0)+(layerz[d]-vertex_Z->at(0))*tan(2.0*atan(exp(-truth_Eta->at(0))))*sin(truth_Phi->at(0));
        muz0[d] = layerz[d];
	mut0[d] = (layerz[d]-vertex_Z->at(0))/299.8/cos(2.0*atan(exp(-truth_Eta->at(0))));
	lphi = atan2(muy0[d],mux0[d]);
        if(lphi<0) lphi = lphi+2*pi;
        if(abs(d-16)>0&&abs(d-16)<9){
          mudial0[d] = round(lphi/pi*4.0+0.5)*2;
        }
        if(abs(d-16)>8&&abs(d-16)<17){
          mudial0[d] = round(lphi/pi*4.0)*2+1;
        }
        if(mudial0[d]==17) mudial0[d]=1;
        mur0[d] = hypot(mux0[d],muy0[d])*cos(atan2(muy0[d],mux0[d])-(mudial0[d]-1)*pi/8.0);
	/*mux0[d] = 0;
	muy0[d] = 0;
	muz0[d] = layerz[d];
	mut0[d] = 0;*/
      }
    }


    
    for(int j = 0; j < Trigger_position->size(); j++){
      for(int k = 0; k < Trigger_position->at(j).size(); k++){
	//if(0<=Trigger_time->at(j).at(k)&&Trigger_time->at(j).at(k)<200){
	ons = 1;

	///////k	
	//f1 = i;
	f2 = j;
	f3 = k;
	intlayer = (2*(Digits_stationEta->at(j)>0)-1)*(8*(Digits_Sector->at(j)==MML)+4*(Digits_multiplet->at(j)-1)+Digits_gasgap->at(j));
	intphi = Digits_stationPhi->at(j);
	inteta = Digits_stationEta->at(j);
	intposition = Trigger_position->at(j).at(k);
	//intvmm = Trigger_VMM->at(j).at(k);
	
	bbb = 0;
	for(int r=0;r<Digits_Gstrip_X->at(j).size();r++){
	  if(intposition==Digits_position->at(j).at(r)){
	    bbb = 1;	
	    platex = Digits_Gstrip_X->at(j).at(r);
	    platey = Digits_Gstrip_Y->at(j).at(r);
	    platez = Digits_Gstrip_Z->at(j).at(r);
	    platelx = Digits_Lstrip_X->at(j).at(r);
	    plately = Digits_Lstrip_Y->at(j).at(r);
	  }
	  if(bbb==1)break;
	}
	
	if(bbb==0){
	  //cout<<"selecting error of art"<<endl;
	  //cout<<"i "<<i<<"j "<<j<<"k "<<k<<"triggerposition "<<intposition<<endl;
	}	
	
	plater = sqrt(platex*platex+platey*platey);
	platephi = atan2(platey,platex);
	plateeta = -log( tan( atan2( sqrt(platex*platex+platey*platey),platez )/2.0) );	

	//if(Trigger_time->size() > j){
	//if(Trigger_time->at(j).size() > 0){
	platet = Trigger_time->at(j).at(k);
	    //}else{
	    //platet = 0;
	    //}
	    //}else{platet = 0;}	
	



	//platert = Digits_stripResponse_time->at(j).at(k);   
	platepos = Trigger_position->at(j).at(k);
	platevmm = Trigger_VMM->at(j).at(k);
	platemmfe = Trigger_MMFE->at(j).at(k);                          
	//	platerpos = Digits_stripResponse_stripPosition->at(j).at(k);
	//	platerx = Digits_stripresponse_stripGposX->at(j).at(k);               
	//	platery = Digits_stripResponse_stripGposY->at(j).at(k);                
	//	platerz = Digits_stripResponse_stripGposZ->at(j).at(k);
	
	platecharge = Trigger_charge->at(j).at(k);

	//layer::the number of its given layer
	platelayer = 0;

	for(int d = 0; d < 33; d++){
	  if( fabs(platez-layerz[d]) < zwidth){
	    platelayer = d-16;//not from 0
	  }
	}

	plateupdown = Digits_stationEta->at(j);
        platemux = mux[platelayer+16];
        platemuy = muy[platelayer+16];
        platemuz = muz[platelayer+16];
	platemut = mut[platelayer+16];
	platemur = mur[platelayer+16];
	platemudial = mudial[platelayer+16];
        platemux0 = mux0[platelayer+16];
        platemuy0 = muy0[platelayer+16];
        platemuz0 = muz0[platelayer+16];
	platemut0 = mut0[platelayer+16];
	platemur0 = mur0[platelayer+16];
        platemudial0 = mudial0[platelayer+16];

	if(platemut!=0){platet-=platemut;}else{platet-=platemut0;}

	//dial::the location of phi
	platedial = 0;//mean error
      	if( platephi>=0 ){
	  if(abs(platelayer)>0&&abs(platelayer)<9){
	    platedial = round(platephi/pi*4.0+0.5)*2;
	  }
          if(abs(platelayer)>8&&abs(platelayer)<17){
            platedial = round(platephi/pi*4.0)*2+1;
          }
	}else{
	  if(abs(platelayer)>0&&abs(platelayer)<9){
            platedial = round((platephi+2*pi)/pi*4.0+0.5)*2;
          }
          if(abs(platelayer)>8&&abs(platelayer)<17){
            platedial = round((platephi+2*pi)/pi*4.0)*2+1;
	  }
	}
	if(platedial == 17) platedial = 1;
	
	strip_phifraction = atan2(platey,platex);
	if(abs(platelayer)>8){
	  while( abs(strip_phifraction) > pi/8.0 ){
	    if(strip_phifraction>0){ strip_phifraction -= pi/4.0;
	    }else{ strip_phifraction += pi/4.0;}
	  }
	}else{
	  strip_phifraction += pi/8.0;
	  while( abs(strip_phifraction) > pi/8.0 ){
	    if(strip_phifraction>0){ strip_phifraction -= pi/4.0;
	    }else{ strip_phifraction += pi/4.0; }
	  }
	}

	hit_phifraction = 0;
	if(platemuz!=0){
	  hit_phifraction = atan2(platemuy,platemux);
	  if(abs(platelayer)>8){
	    while( abs(hit_phifraction) > pi/8.0 ){
	      if(hit_phifraction>0){ hit_phifraction -= pi/4.0;
	      }else{ hit_phifraction += pi/4.0;}
	    }
	  }else{
	    hit_phifraction += pi/8.0;
	    while( abs(hit_phifraction) > pi/8.0 ){
	      if(hit_phifraction>0){ hit_phifraction -= pi/4.0;
            }else{ hit_phifraction += pi/4.0; }
	    }
	  }
	}

        
	hit0_phifraction = atan2(platemuy,platemux);
	if(abs(platelayer)>8){
	  while( abs(hit0_phifraction) > pi/8.0 ){
	    if(hit0_phifraction>0){ hit0_phifraction -= pi/4.0;
	    }else{ hit0_phifraction += pi/4.0;}
	  }
	}else{
	  hit0_phifraction += pi/8.0;
	  while( abs(hit0_phifraction) > pi/8.0 ){
	    if(hit0_phifraction>0){ hit0_phifraction -= pi/4.0;
	    }else{ hit0_phifraction += pi/4.0; }
	  }
	}
	
	
	
	//sign::eta(==0) and stereo +1,-1
	platepror = hypot(platex,platey)*(cos(strip_phifraction)+layersign[platelayer+16]*sin(strip_phifraction)*tan(1.5*pi/180));
	platesign = layersign[platelayer+16];
	platedr = platepror-sqrt(platemux*platemux+platemuy*platemuy)*(cos(hit_phifraction)+layersign[platelayer+16]*sin(hit_phifraction)*tan(1.5*pi/180));
	platedr0 = platepror-sqrt(platemux0*platemux0+platemuy0*platemuy0)*(cos(hit0_phifraction)+layersign[platelayer+16]*sin(hit0_phifraction)*tan(1.5*pi/180));
	
	if(platemudial!=platedial){
	  if(plate_mu_eta*platez>0){//require same A/C side
	    if(hitmux.size()>0){//require hit exists
	      if(platemuz!=0){//require hit exists in same layer
		platedr = hypot(platex-platemux,platey-platemuy);
		platedr0 = hypot(platex-platemux0,platey-platemuy0);
	      }else{//another Sector
                platedr = hypot(platex-platemux0,platey-platemuy0);//virtual
                platedr0 = hypot(platex-platemux0,platey-platemuy0);
	      }
	    }else{//there's no hits
	      platedr = hypot(platex-platemux0,platey-platemuy0);//virtual
	      platedr0 = hypot(platex-platemux0,platey-platemuy0);
	    }
	  }else{//opposite side
	    platedr = 7000;
	    platedr0 = 7000;
	  }
	  
	}

	platedly = platedr0*layersign[platelayer+16]/tan(1.5/180*pi);

	if(platedr>1000) on[0]=1;
	if(platedr>500&&platedr<1000) on[1]=1;
	if(platedr<500) on[2]=1;
	
	if(truth_Pdg->at(0)==13&&platez>0){
	  if(hitmux.size()<=1){
	    platemux0 = 0;
	    platemuy0 = 0;
	    platemut0 = 0;
	    platemur0 = 0;
	    platemudial0 = 0;
	  }//neds hits more than 1	
	  xyz->Fill(); 
	  n = n+1;
	  //cout<<n<<endl;
	}	
	//}
      }
    }
    

      if(on[0]==1) non[0] += 1;
      if(on[0]==0&&on[1]==1) non[1] += 1;
      if(on[0]==0&&on[1]==0&&on[2]==1) non[2] += 1;

      f1+=1;
  }
    file[ndata]->Close();  
  }
  
  xyz->Write();
  filexyz->Close();
  
  
  for(int c=0;c<3;c++){
    cout<<"non"<<c<<non[c]<<endl;
  }

  cout <<"totalevents"<< f1 <<endl;
  cout <<"nstrip"<< n <<endl;
  cout <<"nmuhits"<< nmuhits <<endl;
  
}
