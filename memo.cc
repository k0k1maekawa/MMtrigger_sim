//The original root file is very large and this programs aims to memo
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
#include "TFitter.h"
#include "TMinuit.h"
using namespace std;
using namespace ROOT::Math;

double interseclX(double a1, double b1, double c1, double a2, double b2, double c2, double z) {
  return a1+(z-c1)/(c2-c1)*(a2-a1);
}

double interseclY(double a1, double b1, double c1, double a2, double b2, double c2, double z) {
  return b1+(z-c1)/(c2-c1)*(b2-b1);
}

double intersecX(double a1, double b1, double c1, double theta, double phi, double z) {
  return a1+(z-c1)*tan(theta)*cos(phi);
}

double intersecY(double a1, double b1, double c1, double theta, double phi, double z) {
  return b1+(z-c1)*tan(theta)*sin(phi);
}

double distance1(double a, double b, double c, double d1, double d2, double d3, double x, double y, double z) {
  return ((a-x)*(a-x)+(b-y)*(b-y)+(c-z)*(c-z))-((a-x)*d1+(b-y)*d2+(c-z)*d3)*((a-x)*d1+(b-y)*d2+(c-z)*d3);
}


double distance2(double a, double b, double c, double d1, double d2, double d3, double x, double y, double z, double l1, double l2, double l3) {

  double a1d1_a2d1 = 0;
  double a1d2_a2d2 = 0;
  double cosgamma = 0;

  a1d1_a2d1 = (a*d1+b*d2+c*d3)-(x*d1+y*d2+z*d3);
  a1d2_a2d2 = (a*l1+b*l2+c*l3)-(x*l1+y*l2+z*l3);
  cosgamma = d1*l1+d2*l2+d3*l3;
  return ((a-x)*(a-x)+(b-y)*(b-y)+(c-z)*(c-z))-(a1d1_a2d1*a1d1_a2d1-2.0*cosgamma*a1d1_a2d1*a1d2_a2d2+a1d2_a2d2*a1d2_a2d2)/(1-cosgamma*cosgamma);

}


void trackFunction(int& nDim, double* gout, double& result, double par[], int flg) {
  double sum = 0;
  double a = par[32];
  double b = par[33];
  double theta = par[30];
  double phi = par[31];
  double d1 = sin(theta)*cos(phi);
  double d2 = sin(theta)*sin(phi);
  double d3 = cos(theta);
  double x=0;
  double y=0;
  double z=0;
  double ini=0;

  for(int p=0;p<34;p++){
    ini = par[p];
  }

  for(int i=0;i<10;i++){
    if(fabs(par[3*i])>0.1||fabs(par[3*i+1])>0.1||fabs(par[3*i+2])>0.1){

      x = par[3*i];
      y = par[3*i+1];
      z = par[3*i+2];
      sum += distance1(a,b,0,d1,d2,d3,x,y,z);

    }
  }

  result = sum;
}


void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg) {
  double sum = 0;
  double a = par[62];
  double b = par[63];
  double c = par[64];
  double theta = par[60];
  double phi = par[61];
  double d1 = sin(theta)*cos(phi);
  double d2 = sin(theta)*sin(phi);
  double d3 = cos(theta);
  double x=0;
  double y=0;
  double z=0;
  double l1=0;
  double l2=0;
  double l3=0;
  double ini=0;



  for(int p=0;p<65;p++){
    ini = par[p];
  }

  for(int i=0;i<10;i++){
    if(fabs(par[6*i])>0.1||fabs(par[6*i+1])>0.1||fabs(par[6*i+2])>0.1){

      x = par[6*i];
      y = par[6*i+1];
      z = par[6*i+2];
      l1 = par[6*i+3];
      l2 = par[6*i+4];
      l3 = par[6*i+5];

      sum += distance2(a, b, c, d1, d2, d3, x, y, z, l1, l2, l3);


    }
  }

  result = sum;
}

int signfunc(double x){
  int sgn=0;
  if(x>0) sgn=1;
  if(x<0) sgn=-1;
  return sgn;
}

void memo(){
  string MML="MML";
  int nevent;
  int first; //first or not
  int n; //nstrip and so on
  char namecore[300]="";
  const double lowPtedge = 80000;
  const int nvertexcrit = 100000;//cut off events in which overlaying BG has failed
  const int Ndata = atoi(getenv("NDATA"));
  const int Mu = atoi(getenv("MU"));
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

  //width                                                                                                                                          
  double zwidth = 6;
  double phiwidth = 0.18;

  TH1D *BChist[4];
  for(int pp=0;pp<4;pp++){
    BChist[pp] = new TH1D(Form("BChist%d",pp),Form("BChist%d",pp),Mu+40,-0.5,Mu+39.5);
  }
  cout<<"from "<<-0.5<<" to "<<Mu+39.5<<endl;

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
  unsigned int vertex_n=0;
  vector<double> *vertex_X=0;
  vector<double> *vertex_Y=0;
  vector<double> *vertex_Z=0;
  vector<int> *vertex_Id=0;

  unsigned int *truth_n=0;
  vector<double> *truth_Eta=0;
  vector<double> *truth_Phi=0;
  vector<double> *truth_Pt=0;
  vector<int> *truth_Pdg=0;
  vector<int> *truth_Status=0;
  vector<int> *truth_Production_vertex_id=0;
  vector<int> *truth_End_vertex_id=0;

  vector<double> *hitspos_X=0;
  vector<double> *hitspos_Y=0;
  vector<double> *hitspos_Z=0;
  vector<double> *hitspos_R=0;
  vector<double> *hitstime=0;
  vector<double> *local_hitspos_X=0;
  vector<double> *local_hitspos_Y=0;
  vector<double> *local_hitspos_Z=0;
  vector<int> *hitspdg=0;
  vector<string> *Hits_Sector=0;
  vector<int> *Hits_side=0;
  vector<int> *Hits_multilayer=0;
  vector<int> *Hits_layer=0;
  vector<int> *Hits_stationPhi=0;
  vector<double> *Hits_kineticEnergy=0;
  vector<double> *Hits_DirectionZ=0;

  vector<vector<double> > *Digits_Gstrip_X=0;
  vector<vector<double> > *Digits_Gstrip_Y=0;
  vector<vector<double> > *Digits_Gstrip_Z=0;
  vector<vector<double> > *Digits_Lstrip_X=0;
  vector<vector<double> > *Digits_Lstrip_Y=0;
  vector<vector<double> > *Digits_time=0;
  vector<vector<float> > *Digits_charge=0;
  vector<double> *Digits_stripForTrigger=0;
  vector<double> *Digits_stripTimeForTrigger=0;
  vector<double> *truth_MuEntry_Eta=0;
  vector<double> *truth_MuEntry_Phi=0;
  vector<int> *channel=0;
  vector<vector<double> >  *Digits_stripResponse_time=0;
  vector<vector<int> >  *Digits_stripResponse_stripPosition=0;
  vector<vector<double> >  *Digits_stripresponse_stripGposX=0;
  vector<vector<double> >  *Digits_stripResponse_stripGposY=0;
  vector<vector<double> >  *Digits_stripResponse_stripGposZ=0;
  vector<string> *Digits_Sector=0;
  vector<int> *Digits_stationEta=0;
  vector<int> *Digits_multiplet=0;
  vector<int> *Digits_gasgap=0;
  vector<int> *Digits_stationPhi=0;





  tree->SetBranchAddress("TruthVertex_n",&vertex_n);
  tree->SetBranchAddress("TruthVertex_X",&vertex_X);
  tree->SetBranchAddress("TruthVertex_Y",&vertex_Y);
  tree->SetBranchAddress("TruthVertex_Z",&vertex_Z);
  tree->SetBranchAddress("TruthVertex_Id",&vertex_Id);
  tree->SetBranchAddress("TruthParticle_n",&truth_n);
  tree->SetBranchAddress("TruthParticle_Eta",&truth_Eta);
  tree->SetBranchAddress("TruthParticle_Phi" ,&truth_Phi);
  tree->SetBranchAddress("TruthParticle_Pt",&truth_Pt);
  tree->SetBranchAddress("TruthParticle_Pdg",&truth_Pdg);
  tree->SetBranchAddress("TruthParticle_Status",&truth_Status);
  tree->SetBranchAddress("TruthParticle_Production_vertex_id",&truth_Production_vertex_id);
  tree->SetBranchAddress("TruthParticle_End_vertex_id",&truth_End_vertex_id);

  tree->SetBranchAddress("Hits_MM_hitGlobalPositionX",&hitspos_X);
  tree->SetBranchAddress("Hits_MM_hitGlobalPositionY",&hitspos_Y);
  tree->SetBranchAddress("Hits_MM_hitGlobalPositionZ",&hitspos_Z);
  tree->SetBranchAddress("Hits_MM_hitGlobalPositionR",&hitspos_R);
  tree->SetBranchAddress("Hits_MM_globalTime",&hitstime);
  tree->SetBranchAddress("Hits_MM_hitLocalPositionX",&local_hitspos_X);
  tree->SetBranchAddress("Hits_MM_hitLocalPositionY",&local_hitspos_Y);
  tree->SetBranchAddress("Hits_MM_hitLocalPositionZ",&local_hitspos_Z);
  tree->SetBranchAddress("Hits_MM_particleEncoding" ,&hitspdg); //hitspdg      
  tree->SetBranchAddress("Hits_MM_sim_stationName",&Hits_Sector);
  tree->SetBranchAddress("Hits_MM_sim_side",&Hits_side);
  tree->SetBranchAddress("Hits_MM_sim_multilayer",&Hits_multilayer);
  tree->SetBranchAddress("Hits_MM_sim_layer",&Hits_layer);
  tree->SetBranchAddress("Hits_MM_sim_stationPhi",&Hits_stationPhi);
  tree->SetBranchAddress("Hits_MM_kineticEnergy",&Hits_kineticEnergy);
  tree->SetBranchAddress("Hits_MM_hitGlobalDirectionZ",&Hits_DirectionZ);

  tree->SetBranchAddress("Digits_MM_stripGposX",&Digits_Gstrip_X);
  tree->SetBranchAddress("Digits_MM_stripGposY",&Digits_Gstrip_Y);
  tree->SetBranchAddress("Digits_MM_stripGposZ",&Digits_Gstrip_Z);
  tree->SetBranchAddress("Digits_MM_stripLposX",&Digits_Lstrip_X);
  tree->SetBranchAddress("Digits_MM_stripLposY",&Digits_Lstrip_Y);
  tree->SetBranchAddress("Digits_MM_time",&Digits_time);
  tree->SetBranchAddress("Digits_MM_charge",&Digits_charge);
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



  TCanvas *c1 = new TCanvas("vertex_n","vertex_n");
  tree->Draw("TruthVertex_n");
  c1->Print("vertex_n_check.pdf");
  TCanvas *c2 = new TCanvas("particle_n","particle_n");
  tree->Draw("TruthParticle_n");
  c2->Print("particle_n_check.pdf");

  TH1F *h;
  double nbg;
  double nmiddle;
  double nsingle;

  tree->Draw("TruthVertex_n>>h",Form("TruthVertex_n>%d",nvertexcrit),"goff");
  h = (TH1F*)gDirectory->Get("h");
  nbg = h->GetEntries();
  cout<<"nbg"<<nbg<<endl;

  tree->Draw("TruthVertex_n>>h",Form("TruthVertex_n>1&&TruthVertex_n<=%d",nvertexcrit),"goff");
  h = (TH1F*)gDirectory->Get("h");
  nmiddle = h->GetEntries();
  cout<<"nmiddle"<<nmiddle<<endl;

  tree->Draw("TruthVertex_n>>h","TruthVertex_n==1","goff");
  h = (TH1F*)gDirectory->Get("h");
  nsingle = h->GetEntries();
  cout<<"nsingle"<<nsingle<<endl;  

  int nBC;
  double layerz[33] = {0};
  int layersign[33] = {0};
  int intlayer=0;

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
            layersign[intlayer] = signfunc(Digits_Lstrip_Y->at(j).at(k));

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
            layersign[intlayer] = -signfunc(intlayer-16)*signfunc(Digits_Lstrip_Y->at(j).at(k));

          }
        }
      }
    }
  }

  for(int d=0;d<33;d++){
    cout<<"layer"<<d-16<<"\t"<<layerz[d]<<"\t"<<layersign[d]<<endl;
  }



  first = 0; //reset
  n = 0;
  TFile *filememo = new TFile(Form("%smemo.root",namefordata),"recreate");
  TTree *memo = new TTree("memo","memo");

  int ndata = 0;
  int platenvertex = 0;

  int nmu[33]={0};
  double mux[33]={0};
  double muy[33]={0};
  double muz[33]={0};
  double mut[33]={0};
  double hiteta[33] = {0};
  double hitphi[33] = {0};
  double mufirst[33]={0};
  double mudial[33]={0};
  double mux0[33]={0};
  double muy0[33]={0};
  double muz0[33]={0};
  double mut0[33]={0};
  double hiteta0[33] ={0};
  double hitphi0[33] ={0};
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

  double platex0=0;
  double platey0=0;
  double platez0=0;
  double plateeta0=0;
  double platephi0=0;
  double platept0=0;
  int platepdg0=0;
  double plate_muentry_eta=0;
  //for fit
  double hit_ltheta[2]={0};
  double hit_gtheta[2]={0};
  double hit_ltheta3d=0;
  double hit_lphi3d=0;
  double hit_gtheta3d[2]={0};
  double hit_gphi3d[2]={0};
  double hit_eta[2]={0};
  double hit_phi[2]={0};
  double hit_dtheta[2]={0};

  double detecZ[2] = {0};
  detecZ[0] = (layerz[20]+layerz[21])/2.0;
  detecZ[1] = (layerz[28]+layerz[29])/2.0;
  double detecX = 0;
  double detecY = 0;


  double p1 = -1;
  TFitter* track = new TFitter(34);
  track->ExecuteCommand("SET PRINTOUT",&p1,1);
  track->SetFCN(trackFunction);
  TF1 *linear = new TF1("lin","tan([0])*(x-[1])");
  TF1 *linearl = new TF1("linl","[0]*(x-[1])");

  //vector def                                                                                         
  vector<double> hitxofe(10);
  vector<double> hityofe(10);
  vector<double> hitxyofe(10);
  vector<double> hitxofla(10);
  vector<double> hityofla(10);
  vector<double> hitxyofla(10);
  vector<double> hitxofsm(10);
  vector<double> hityofsm(10);
  vector<double> hitxyofsm(10);
  vector<double> hitexofe(10);
  vector<double> hiteyofe(10);
  vector<double> hit3dxofe(10);
  vector<double> hit3dyofe(10);
  vector<double> hit3dzofe(10);
  for(int c=0;c<10;c++){
    hitexofe.at(c) = 0;
    hiteyofe.at(c) = 0;
  }

  memo->Branch("ndata",&ndata);
  memo->Branch("nvertex",&platenvertex);//int
  memo->Branch("mux",mux,"mux[33]/D");
  memo->Branch("muy",muy,"muy[33]/D");
  memo->Branch("muz",muz,"muz[33]/D");
  memo->Branch("mut",mut,"mut[33]/D");
  memo->Branch("hiteta",hiteta,"hiteta[33]/D");
  memo->Branch("hitphi",hitphi,"hitphi[33]/D");
  memo->Branch("mux0",mux0,"mux0[33]/D");
  memo->Branch("muy0",muy0,"muy0[33]/D");
  memo->Branch("muz0",muz0,"muz0[33]/D");
  memo->Branch("mut0",mut0,"mut0[33]/D");
  memo->Branch("hiteta0",hiteta0,"hiteta0[33]/D");
  memo->Branch("hitphi0",hitphi0,"hitphi0[33]/D");
  memo->Branch("x0",&platex0);
  memo->Branch("y0",&platey0);
  memo->Branch("z0",&platez0);
  memo->Branch("vertex_n",&vertex_n);//unsigned int
  memo->Branch("pdg0",&platepdg0);
  memo->Branch("eta0",&plateeta0);
  memo->Branch("pt0" ,&platept0);
  memo->Branch("phi0",&platephi0);
  memo->Branch("hit_ltheta",hit_ltheta,"hit_ltheta[2]/D");
  memo->Branch("hit_gtheta",hit_gtheta,"hit_gtheta[2]/D");
  memo->Branch("hit_ltheta3d",&hit_ltheta3d);
  memo->Branch("hit_lphi3d",&hit_lphi3d);
  memo->Branch("hit_gtheta3d",hit_gtheta3d,"hit_gtheta3d[2]/D");
  memo->Branch("hit_gphi3d",hit_gphi3d,"hit_gphi3d[2]/D");
  memo->Branch("hit_eta",hit_eta,"hit_eta[2]/D");
  memo->Branch("hit_phi",hit_phi,"hit_phi[2]/D");
  memo->Branch("hit_dtheta",hit_dtheta,"hit_dtheta[2]/D");

  for(int nd=0;nd<Ndata;nd++){
    ndata = nd;
    tree = (TTree*)file[ndata]->Get("NSWHitsTree");
    nevent = tree->GetEntries();
    cout<<"DATA"<<ndata<<" "<<nevent<<"event"<<endl;
    tree->SetBranchAddress("TruthVertex_n",&vertex_n);
    tree->SetBranchAddress("TruthVertex_X",&vertex_X);
    tree->SetBranchAddress("TruthVertex_Y",&vertex_Y);
    tree->SetBranchAddress("TruthVertex_Z",&vertex_Z);
    tree->SetBranchAddress("TruthVertex_Id",&vertex_Id);
    tree->SetBranchAddress("TruthParticle_n",&truth_n);
    tree->SetBranchAddress("TruthParticle_Eta",&truth_Eta);
    tree->SetBranchAddress("TruthParticle_Phi" ,&truth_Phi);
    tree->SetBranchAddress("TruthParticle_Pt",&truth_Pt);
    tree->SetBranchAddress("TruthParticle_Pdg",&truth_Pdg);
    tree->SetBranchAddress("TruthParticle_Status",&truth_Status);
    tree->SetBranchAddress("TruthParticle_Production_vertex_id",&truth_Production_vertex_id);
    tree->SetBranchAddress("TruthParticle_End_vertex_id",&truth_End_vertex_id); 

    tree->SetBranchAddress("Hits_MM_hitGlobalPositionX",&hitspos_X);
    tree->SetBranchAddress("Hits_MM_hitGlobalPositionY",&hitspos_Y);
    tree->SetBranchAddress("Hits_MM_hitGlobalPositionZ",&hitspos_Z);
    tree->SetBranchAddress("Hits_MM_hitGlobalPositionR",&hitspos_R);
    tree->SetBranchAddress("Hits_MM_globalTime",&hitstime);
    tree->SetBranchAddress("Hits_MM_hitLocalPositionX",&local_hitspos_X);
    tree->SetBranchAddress("Hits_MM_hitLocalPositionY",&local_hitspos_Y);
    tree->SetBranchAddress("Hits_MM_hitLocalPositionZ",&local_hitspos_Z);
    tree->SetBranchAddress("Hits_MM_particleEncoding" ,&hitspdg); //hitspdg      
    tree->SetBranchAddress("Hits_MM_sim_stationName",&Hits_Sector);
    tree->SetBranchAddress("Hits_MM_sim_side",&Hits_side);
    tree->SetBranchAddress("Hits_MM_sim_multilayer",&Hits_multilayer);
    tree->SetBranchAddress("Hits_MM_sim_layer",&Hits_layer);
    tree->SetBranchAddress("Hits_MM_sim_stationPhi",&Hits_stationPhi);
    tree->SetBranchAddress("Hits_MM_kineticEnergy",&Hits_kineticEnergy);
    tree->SetBranchAddress("Hits_MM_hitGlobalDirectionZ",&Hits_DirectionZ);

    tree->SetBranchAddress("Digits_MM_stripGposX",&Digits_Gstrip_X);
    tree->SetBranchAddress("Digits_MM_stripGposY",&Digits_Gstrip_Y);
    tree->SetBranchAddress("Digits_MM_stripGposZ",&Digits_Gstrip_Z);
    tree->SetBranchAddress("Digits_MM_stripLposX",&Digits_Lstrip_X);
    tree->SetBranchAddress("Digits_MM_stripLposY",&Digits_Lstrip_Y);
    tree->SetBranchAddress("Digits_MM_time",&Digits_time);
    tree->SetBranchAddress("Digits_MM_charge",&Digits_charge);
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



    for(int i = 0; i < nevent; i++){
      if((i+1)%100000 == 0) cout<<i+1<<endl;
      tree->GetEntry(i);
      //cout<<"Entry"<<i<<endl;
      platenvertex = vertex_n;
      
      nBC = 0;
      for(int c=0;c<truth_Status->size();c++){
	if(truth_Status->at(c)==4){
	  nBC += 1;
	}
      }
      if(nBC%2!=0)cout<<"odd number proton"<<endl;
      nBC = nBC/2;
      //cout<<nBC<<endl;
      BChist[3]->Fill(nBC);
     
      if(platenvertex>nvertexcrit)BChist[0]->Fill(nBC);
      if(platenvertex>1&&platenvertex<=nvertexcrit)BChist[1]->Fill(nBC);
      if(platenvertex==1)BChist[2]->Fill(nBC);

      if(vertex_X->size() > 0){
	platex0 = vertex_X->at(0);
	platey0 = vertex_Y->at(0);
	platez0 = vertex_Z->at(0);
      }else{
	cout<<"Vertex Error"<<endl;	
      }
      if(truth_Eta->size() > 0){
	plateeta0 = truth_Eta->at(0);
	platephi0 = truth_Phi->at(0);
	platept0 = truth_Pt->at(0);
	platepdg0 = truth_Pdg->at(0);
      }else{
        plateeta0 = 0;
	platephi0 = 2*pi;
	platept0 = 0;
	platepdg0 = 0;
      }
      if(truth_MuEntry_Eta->size() > 0){
	plate_muentry_eta = truth_MuEntry_Eta->at(0);
      }else{
	plate_muentry_eta = 0;
      }
      
      
      //hits
      for(int d=0;d<33;d++){
	nmu[d] = 0;
	mux[d] = 0;
	muy[d] = 0;
	muz[d] = 0;
	mut[d] = 0;
	mufirst[d] = 0;
	mudial[d] = 0;
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
      hitxofla.erase(hitxofla.begin(),hitxofla.end());
      hityofla.erase(hityofla.begin(),hityofla.end());
      hitxyofla.erase(hitxyofla.begin(),hitxyofla.end());
      hitxofsm.erase(hitxofsm.begin(),hitxofsm.end());
      hityofsm.erase(hityofsm.begin(),hityofsm.end());
      hitxyofsm.erase(hitxyofsm.begin(),hitxyofsm.end());
      
      for(int c=0;c<hitspos_X->size();c++){
	if(hitspdg->at(c)==truth_Pdg->at(0)&&hitstime->at(c)>22&&hitstime->at(c)<32&&signfunc(plateeta0)==signfunc(hitspos_Z->at(c))){
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
	if(signfunc(plateeta0)==signfunc(d-16)){
	  if(hitmux.size()>1){
	    mux0[d] = interseclX(hitmux.at(0),hitmuy.at(0),hitmuz.at(0),hitmux.at(hitmux.size()-1),hitmuy.at(hitmuy.size()-1),hitmuz.at(hitmuz.size()-1),layerz[d]);
	    muy0[d] = interseclY(hitmux.at(0),hitmuy.at(0),hitmuz.at(0),hitmux.at(hitmux.size()-1),hitmuy.at(hitmuy.size()-1),hitmuz.at(hitmuz.size()-1),layerz[d]);
	    muz0[d] = layerz[d];
	    mut0[d] = interseclX(hitmut.at(0),hitmuy.at(0),hitmuz.at(0),hitmut.at(hitmux.size()-1),hitmuy.at(hitmuy.size()-1),hitmuz.at(hitmuz.size()-1),layerz[d]);

	    phifraction = atan2(muy0[d],mux0[d]);
	    if(abs(d-16)>=9){
	      while( fabs(phifraction) > pi/8.0 ){
		if(phifraction>0){ phifraction -= pi/4.0;
		}else{ phifraction += pi/4.0; }
	      }
	    }else{
	      phifraction += pi/8.0;
	      while( fabs(phifraction) > pi/8.0 ){
		if(phifraction>0){ phifraction -= pi/4.0;
		}else{ phifraction += pi/4.0; }
	      }
	    }
	    if(layersign[d]==0){
	      /*in 2D, comparison for Xstrip is only needed, so do not reckon 1.5deg*/
	      if(abs(d-16)>=9){
		hitxofla.push_back(muz0[d]);
		hityofla.push_back(hypot(mux0[d],muy0[d])*cos(phifraction));
		hitxyofla.push_back(muz0[d]*hypot(mux0[d],muy0[d])*cos(phifraction));
	      }else{
		hitxofsm.push_back(muz0[d]);
                hityofsm.push_back(hypot(mux0[d],muy0[d])*cos(phifraction));
                hitxyofsm.push_back(muz0[d]*hypot(mux0[d],muy0[d])*cos(phifraction));
	      }
	    }
	 
	  }else{
	    /*mux0[d] = vertex_X->at(0)+(layerz[d]-vertex_Z->at(0))*tan(2.0*atan(exp(-truth_Eta->at(0))))*cos(truth_Phi->at(0));
	      muy0[d] = vertex_Y->at(0)+(layerz[d]-vertex_Z->at(0))*tan(2.0*atan(exp(-truth_Eta->at(0))))*sin(truth_Phi->at(0));
	      muz0[d] = layerz[d];
	      mut0[d] = (layerz[d]-vertex_Z->at(0))/299.8/cos(2.0*atan(exp(-truth_Eta->at(0))));*/
	    mux0[d] = 0;        
	    muy0[d] = 0;                                                     
	    muz0[d] = layerz[d];                                                                                     
	    mut0[d] = 0;
	  }
	  hiteta[d] = -log(tan(atan2(hypot(mux[d],muy[d]),muz[d])/2.0));
	  hitphi[d] = atan2(muy[d],mux[d]);
	  if(muz[d]==0){
	    hiteta[d]=0;
	    hitphi[d]=2*pi;
	  }
	  
	  hiteta0[d] = -log(tan(atan2(hypot(mux0[d],muy0[d]),muz0[d])/2.0));
	  hitphi0[d] = atan2(muy0[d],mux0[d]);
	  if(mux0[d]==0&&muy0[d]==0){
	    hiteta0[d]=0;
	    hitphi0[d]=2*pi;
	  }
	}else{
	  hiteta[d] = 0;
	  hitphi[d] = 0;
	  hiteta0[d] = 0;
	  hitphi0[d] = 0;
	  mux0[d] = 0;
	  muy0[d] = 0;
	  muz0[d] = 0;
	  mut0[d] = 0;
	}
      }
      
      //truths
      hit_ltheta3d=0;
      hit_lphi3d=0;
      for(int d=0;d<2;d++){
	hit_ltheta[d]=0;
	hit_gtheta[d]=0;
	hit_gtheta3d[d]=0;
	hit_gphi3d[d]=0;
	hit_eta[d]=0;
	hit_phi[d]=0;
	hit_dtheta[d]=0;
      }
      if(hitmux.size()>1){
	for(int d=0;d<2;d++){
	  if(d==1){
	    hit_ltheta[d] = atan2(hityofla.at(hityofla.size()-1)-hityofla.at(0),hitxofla.at(hitxofla.size()-1)-hitxofla.at(0));
	    detecY = interseclY(hitxofla.at(0),hityofla.at(0),hitxofla.at(0),hitxofla.at(hitxofla.size()-1),hityofla.at(hityofla.size()-1),hitxofla.at(hitxofla.size()-1),detecZ[d]);
	    hit_gtheta[d] = atan2(detecY,detecZ[d]);
	    hit_eta[d] = -log(tan(hit_gtheta[d]/2.0));
	  }
	  if(d==0){
	    hit_ltheta[d] = atan2(hityofsm.at(hityofsm.size()-1)-hityofsm.at(0),hitxofsm.at(hitxofsm.size()-1)-hitxofsm.at(0));
	    detecY = interseclY(hitxofsm.at(0),hityofsm.at(0),hitxofsm.at(0),hitxofsm.at(hitxofsm.size()-1),hityofsm.at(hityofsm.size()-1),hitxofsm.at(hitxofsm.size()-1),detecZ[d]);
	    hit_gtheta[d] = atan2(detecY,detecZ[d]);
	    hit_eta[d] = -log(tan(hit_gtheta[d]/2.0));
	  }
	}
	hit_ltheta3d = atan2(hypot(hitmux.at(hitmux.size()-1)-hitmux.at(0),hitmuy.at(hitmuy.size()-1)-hitmuy.at(0)),hitmuz.at(hitmuz.size()-1)-hitmuz.at(0));
	hit_lphi3d = atan2(hitmuy.at(hitmuy.size()-1)-hitmuy.at(0),hitmux.at(hitmux.size()-1)-hitmux.at(0));
	for(int d=0;d<2;d++){
	  detecX = interseclX(hitmux.at(0),hitmuy.at(0),hitmuz.at(0),hitmux.at(hitmux.size()-1),hitmuy.at(hitmuy.size()-1),hitmuz.at(hitmuz.size()-1),detecZ[d]);
	  detecY = interseclY(hitmux.at(0),hitmuy.at(0),hitmuz.at(0),hitmux.at(hitmux.size()-1),hitmuy.at(hitmuy.size()-1),hitmuz.at(hitmuz.size()-1),detecZ[d]);
	  hit_gtheta3d[d] = atan2(hypot(detecX,detecY),detecZ[d]);
	  hit_gphi3d[d] = atan2(detecY,detecX);
	  hit_phi[d] = hit_gphi3d[d];
	  hit_dtheta[d] = hit_ltheta3d - hit_gtheta3d[d];
	}
      }

      memo->Fill();
      n+= 1;  
      
    }
    file[ndata]->Close();
  }
 
  memo->Write();
  filememo->Close();
   
  cout <<"nmemo"<< n <<endl;
  TCanvas *BCcan[4];  
  for(int ff=0;ff<4;ff++){
    BCcan[ff] = new TCanvas(Form("BCcan%d",ff),Form("BCcan%d",ff));
    BChist[ff]->Draw();
    if(ff==0)BCcan[ff]->Print(Form("BChist_bg_%s_check.pdf",namefordata));
    if(ff==1)BCcan[ff]->Print(Form("BChist_middle_%s_check.pdf",namefordata));
    if(ff==2)BCcan[ff]->Print(Form("BChist_single_%s_check.pdf",namefordata));
    if(ff==3)BCcan[ff]->Print(Form("BChist_all_%s_check.pdf",namefordata));
  }
  
}
