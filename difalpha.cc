#include "/home/maekawa/atlasstyle-00-03-05/AtlasStyle.C"
#include "/home/maekawa/atlasstyle-00-03-05/AtlasLabels.C"
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
#include <iomanip>
#include "TFitter.h"
#include "TGraphErrors.h"
using namespace std;
using namespace ROOT::Math;

int SectorAllpat(int n,int *vec){
  int Allpat = 0;
  int layeron[8]={0};
  int d = 0;
  for(int i=0;i<n;i++){
    d = *vec;
    layeron[(d-1)%8] = 1;
    ++vec;
  }
  int dlack=0;
  for(int d=7;d>-1;d--){
    if(layeron[d]==0){
      Allpat += (d+1)*pow(10,dlack);
      dlack += 1;
    }
  }
  return Allpat;
}


int SectorUVpat(int n,int *vec){
  int UVpat = 0;
  int layeron[8]={0};
  int d = 0;
  for(int i=0;i<n;i++){
    d = *vec;
    layeron[(d-1)%8] = 1;
    ++vec;
  }
  int dlack=0;
  for(int d=5;d>1;d--){
    if(layeron[d]==0){
      UVpat += (d-1)*pow(10,dlack);
      dlack += 1;
    }
  }
  return UVpat;
}

int SectorXpat(int n,int *vec){
  int Xpat = 0;
  int layeron[8]={0};
  int d = 0;
  for(int i=0;i<n;i++){
    d = *vec;
    layeron[(d-1)%8] = 1;
    ++vec;
  }
  int dlack = 0;
  for(int d=7;d>5;d--){
    if(layeron[d]==0){
      Xpat += (d-3)*pow(10,dlack);
      dlack += 1;
    }
  }
  for(int d=1;d>-1;d--){
    if(layeron[d]==0){
      Xpat += (d+1)*pow(10,dlack);
      dlack += 1;
    }
  }
  return Xpat;
}


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

void XFitFunction(int& nDim, double* gout, double& result, double par[], int flg) {

  double sum = 0;
  double a = par[8];
  double b = par[9];
  double r=0;
  double z=0;
  double ini=0;


  for(int p=0;p<10;p++){
    ini = par[p];
  }

  for(int i=0;i<4;i++){
    if(fabs(par[2*i])>0.1||fabs(par[2*i+1])>0.1){

      z = par[2*i];
      r = par[2*i+1];
      sum += sqrt((r-a*z-b)*(r-a*z-b));
  
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


void difalpha(){
  SetAtlasStyle();
  int nevent;
  int drcut=0;
  int first = 0;
  int n; //nstrip and so on
  double eta1 = 1.2;
  double eta2 = 2.9;
  char namecore[200]="";
  char filename[200];
  const int shortcut = 0;
  const char* core;
  const int checknum=10;
  const double zwidth=6;
  const int fixon = 1;//fix acceptance definition
  const int checkexp = 1;//express errors
  const int bg = atoi(getenv("BG"));
  int nvertexcrit = 0;
  if(bg==1) nvertexcrit = 100000;
  const int uvfit = 1;//"0" meancs no cut for UV, "1" means 2U(2V) is valid but 1U(1V) is not valid, "2" means 2U(2V) and these two < etaslopewidth is valid.  
  const double slopewidth = atof(getenv("ETASLOPEWIDTH"));
  const int robust = 0;//robust algorithm
  const int majority = 0;//majority algorithm
  const int analytic = 1;//analytic 2d(Xonly) fit  
  const int l3d = 1;//calculate local_3d
  const double drcutin = 1;//for tight cut of the signal from secondary particles  The unit is mm.
  const int Ndata = atoi(getenv("NDATA"));
  const int truthtracktype = 0; //the way you treat the truth track
  //0:tie mu start-point with mu end-point 1:fit for all mupoints 2:treat as curve 
  const int UVtruth = 1; //0:easy virtual plain (in R-Z 2d plane)  1:virtual plain (calculated in 3d)
  const int Xtruth = 1; //0:chi^2 fit 1:virtual plain 
  //constant----------------------------
  const double pi = TMath::Pi();
  const double lowPtedge = 80000;

  if(bg==0) cout<<"Single Mu"<<endl;
  if(bg==1) cout<<"BG"<<endl;
  if(uvfit==0) cout<<"xuvfit"<<endl;
  if(uvfit==1) cout<<"uvfit"<<endl;
  if(uvfit==2) cout<<"exceeding uvfit"<<endl;

  //Plate  
  double error;
  //  double max[10] = {0};
  //  double min[10] = {0};
  double test = 0;
  int alert[checknum] = {0};

  double aofz0=0;
  double bofz0=0;
  double zintercept3d=0;
  double gtheta3d=0;
  double ltheta3d=0;
  double gphi3d=0;
  double lphi3d=0;
  double hit_aofz0=0;
  double hit_bofz0=0;
  double truth_aofz0=0;
  double truth_bofz0=0;
  double hit_3dtheta=0;
  double hit_3dphi=0;
  double hit_gtheta3d=0;
  double hit_gphi3d=0;
  TFitter* minuit = new TFitter(65);
  double p1 = -1;
  minuit->ExecuteCommand("SET PRINTOUT",&p1,1);
  minuit->SetFCN(minuitFunction);
  TFitter* track = new TFitter(34);
  track->ExecuteCommand("SET PRINTOUT",&p1,1);
  track->SetFCN(trackFunction);
  TFitter* xfit = new TFitter(10);
  xfit->ExecuteCommand("SET PRINTOUT",&p1,1);
  xfit->SetFCN(XFitFunction);

  double delta[5] = {0};
  TF1 *linear = new TF1("lin","tan([0])*(x-[1])");
  TF1 *linearl = new TF1("linl","[0]*(x-[1])");
  //vector def
  vector<double> xofe(4);
  vector<double> yofe(4);
  vector<double> xyofe(4);
  vector<int> dofe(4);
  vector<int> dofs(4);
  vector<double> rofs(4);
  vector<int> signofs(4);
  vector<double> exofe(4);
  vector<double> eyofe(4);
  vector<double> x3d(10);
  vector<double> y3d(10);
  vector<double> z3d(10);
  vector<int> d3d(10);
  vector<double> l1ofe(10);
  vector<double> l2ofe(10);
  vector<double> l3ofe(10);
  vector<double> hitxofe(10);
  vector<double> hityofe(10);
  vector<double> hitxyofe(10);
  vector<double> hitexofe(10);
  vector<double> hiteyofe(10);
  vector<double> hit3dxofe(10);
  vector<double> hit3dyofe(10);
  vector<double> hit3dzofe(10);
  for(int c=0;c<4;c++){
    exofe.at(c) = 5.;
    eyofe.at(c) = 0.435;
  }
  for(int c=0;c<10;c++){
    hitexofe.at(c) = 0;
    hiteyofe.at(c) = 0;
  }

  //for diffuse
  int nmu[33];
  double mur[33];
  double mut[33];
  double flac;
  int nother;
  int nmuhit;

  double layerz[33] = {0};
  double layerzcent[33] = {0};
  int layersign[33] = {0};

  double ltheta=0;
  double b=0;
  double errltheta=0;
  double errb=0;
  double hit_gtheta=0;
  double hit_ltheta=0;
  double gtheta=0;
  double errgtheta=0;
  int dd=0;
  double resid=0;
  double zintercept=0;
  double rintercept=0;
  double hit_zintercept=0;
  double hit_rintercept=0;
  string MML = "MML";

  TH1D *residhist[6];
  //TH2D *residhist[2];
  residhist[0] = new TH1D("residual_error","residual error for hit_ltheta",500,-0.002,0.002);
  //residhist[0] = new TH2D("residual_error","residual error for hit_ltheta",500,-0.002,0.002,10000,-8000,8000);
  residhist[0]->GetXaxis()->SetTitle("residual_error/mm");
  residhist[1] = new TH1D("residual_error","residual error for hit_gtheta",500,-5.,5.);
  //residhist[1] = new TH2D("residual_error","residual error for hit_gtheta",500,-5.,5.,10000,-8000,8000);
  residhist[1]->GetXaxis()->SetTitle("residual_error/mm");
  residhist[2] = new TH1D("residual_error","residual error for 3dhit_ltheta",500,-0.01,0.04);
  residhist[2]->GetXaxis()->SetTitle("residual_error/mm");
  residhist[3] = new TH1D("vetex_to_track","distance of vertex to track",400,-1,19);
  residhist[3]->GetXaxis()->SetTitle("distance of vertex to track/mm");
  residhist[4] = new TH1D("residual_error","residual error for ltheta",500,-5,5);
  residhist[4]->GetXaxis()->SetTitle("residual_error/mm");
  residhist[5] = new TH1D("residual_error","residual error for gtheta",500,-5,5);
  residhist[5]->GetXaxis()->SetTitle("residual_error/mm");


  TH2D *intercepthist[4];//intercepthist
  intercepthist[0] = new TH2D("hit_xyatz0-truth","hit_xyofz0-truth",500,-500.,500.,500,-500.,500.);
  intercepthist[0]->GetXaxis()->SetTitle("hit_xatz0-truth/mm");
  intercepthist[0]->GetYaxis()->SetTitle("hit_yatz0-truth/mm");
  intercepthist[1] = new TH2D("hitintercept_and_hit_dtheta","hitintercept and hit_dtheta",500,-2500.,2500.,500,-0.1,0.1);
  intercepthist[1]->GetXaxis()->SetTitle("zintercept/mm");
  intercepthist[1]->GetYaxis()->SetTitle("dtheta/rad");
  intercepthist[2] = new TH2D("xyatz0-truth","xyatz0-truth",500,-500.,500.,500,-500.,500.);
  intercepthist[2]->GetXaxis()->SetTitle("xatz0-truth/mm");
  intercepthist[2]->GetYaxis()->SetTitle("yatz0-truth/mm");
  intercepthist[3] = new TH2D("hitdtheta_and_dtheta","hit_dtheta and dtheta",500,-0.1,0.1,500,-1.,1.);
  intercepthist[3]->GetXaxis()->SetTitle("hit_dtheta/rad");
  intercepthist[3]->GetYaxis()->SetTitle("dtheta/rad");




  cout <<"Please Enter the Filename before \".root\""<<endl;
  strcpy(namecore,getenv("NAMECORE"));
  drcut = atof(getenv("DRCUT"));
  TFile *filefast = new TFile(Form("%sfast.root",namecore));
  TTree *tree = (TTree*)filefast->Get("fast");

  char *ret;
  char namefordata[200] = "";
  ret = strstr(Form("%s",namecore),"drcut");
  if(ret!=NULL){
    int drcutlength = strlen(Form("%s",ret));
    int namecorelength = strlen(Form("%s",namecore));
    strncat(namefordata,Form("%s",namecore),namecorelength-drcutlength);
  }else{
    strcat(namefordata,Form("%s",namecore));
    //drcut = 1000000;
  }
  TFile *filememo = new TFile(Form("%smemo.root",namefordata));
  TTree *memo = (TTree*)filememo->Get("memo");
  TFile *file2[Ndata];
  file2[0] = new TFile(Form("%s.root",namefordata));
  for(int d=1;d<Ndata;d++){
    file2[d] = new TFile(Form("%s_%d.root",namefordata,d));
  }
  int ninc[Ndata]={0};
  int Ninc = 0;
  for(int d=0;d<Ndata;d++){
    ninc[d] = Ninc;
    cout<<"for"<<d<<"th data "<<ninc[d]<<"event"<<endl;
    TTree *data = (TTree*)file2[d]->Get("NSWHitsTree");
    Ninc += data->GetEntries();
  }
  TTree *data = (TTree*)file2[0]->Get("NSWHitsTree");
  if(memo->GetEntries()!=Ninc)cout<<"Sum up Error!!"<<endl;
  TCanvas *c4 = new TCanvas("");
  TCanvas *c5 = new TCanvas("");



  //Treesetting
  double hitphi0[33] = {0};
  double hiteta0[33] = {0};
  double hitz[33] = {0};
  double hitx0[33] = {0};
  double hity0[33] = {0};
  double hitz0[33] = {0};
  double hitt0[33] = {0};

  memo->SetBranchAddress("hiteta0",hiteta0);
  memo->SetBranchAddress("hitphi0",hitphi0);
  memo->SetBranchAddress("muz",hitz);
  memo->SetBranchAddress("mux0",hitx0);
  memo->SetBranchAddress("muy0",hity0);
  memo->SetBranchAddress("muz0",hitz0);
  memo->SetBranchAddress("mut0",hitt0);

  int ndata = 0;
  int i = 0;
  double time[33] = {0};
  int layer[33] = {0};
  double x[33] = {0};
  double y[33] = {0};
  double z[33] = {0};
  //int nmu[33]={0};
  double mufirst[33]={0};
  /*  double mux[33] = {0};
  double muy[33] = {0};
  double muz[33] = {0};
  double mut[33] = {0};
  double mux0[33] = {0};
  double muy0[33] = {0};
  double muz0[33] = {0};
  double mut0[33] = {0};*/
  double lx[33] = {0};
  double ly[33] = {0};
  double dr0[33] = {0};
  int sign[33] = {0};
  int dial = 0;
  double mueta = 0;
  double muentryeta = 0;
  int nhit = 0;
  int neta = 0;
  int nstereo = 0;
  int sector=0;
  int band = 0;

  tree->SetBranchAddress("ndata",&ndata);
  tree->SetBranchAddress("i",&i);
  tree->SetBranchAddress("layer",layer);
  tree->SetBranchAddress("time",time);
  tree->SetBranchAddress("dial",&dial);
  tree->SetBranchAddress("mueta",&mueta);
  tree->SetBranchAddress("muentryeta",&muentryeta);
  tree->SetBranchAddress("nhit",&nhit);
  tree->SetBranchAddress("neta",&neta);
  tree->SetBranchAddress("nstereo",&nstereo);
  tree->SetBranchAddress("x",x);
  tree->SetBranchAddress("y",y);
  tree->SetBranchAddress("z",z);
  /*tree->SetBranchAddress("mux",mux);
  tree->SetBranchAddress("muy",muy);
  tree->SetBranchAddress("muz",muz);
  tree->SetBranchAddress("mut",mut);
  tree->SetBranchAddress("mux0",mux0);
  tree->SetBranchAddress("muy0",muy0);
  tree->SetBranchAddress("muz0",muz0);
  tree->SetBranchAddress("mut0",mut0);*/
  tree->SetBranchAddress("lx",lx);
  tree->SetBranchAddress("ly",ly);
  tree->SetBranchAddress("dr0",dr0);
  tree->SetBranchAddress("sign",sign);
  tree->SetBranchAddress("sector",&sector);
  tree->SetBranchAddress("band",&band);

  //int numbers

  int efef[33]={0};
  //layer
  int ef[33][3]={{0}};
  //layer:time

  int en[5][3]={{0}};
  //neta:time (8 layer hit 4 stereo layer hit) 
  int esn[5][5][3]={{{0}}};
  //neta:nstereo:time (8 layer hit 4 stereo layer hit) 
  double esdelta[5][5][3]={{{0}}};
  double edelta[5][3]={{0}};
  double hit_esdelta[5][5][3]={{{0}}};
  double hit_edelta[5][3]={{0}};
  double sigma_sltheta[5][5][3]={{{0}}};
  double sigma_ltheta[5][3]={{0}};
  double sigma_sgtheta[5][5][3]={{{0}}};
  double sigma_gtheta[5][3]={{0}};

  int def[33]={0};//denomi of each
  //layer
  int dnf[17]={0};//denomi of nthrough
  //nthrough
  int dsnf[17][9]={{0}};//denomi of nthrough and nstereo 
  //nthrough:nstereothrough 
  double x0=0;
  double y0=0;
  double z0=0;
  double eta0=0;
  double phi0=0;
  double pt0=0;
  int pdg0=0;
  unsigned int vertex_n=0;

  memo->SetBranchAddress("x0",&x0);
  memo->SetBranchAddress("y0",&y0);
  memo->SetBranchAddress("z0",&z0);
  memo->SetBranchAddress("vertex_n",&vertex_n);
  memo->SetBranchAddress("pdg0",&pdg0);
  memo->SetBranchAddress("eta0",&eta0);
  memo->SetBranchAddress("pt0" ,&pt0);
  memo->SetBranchAddress("phi0",&phi0);

  //for fitting

  int nfast = tree->GetEntries();
  nevent = data->GetEntries();
  cout<<nfast<<"nfast"<<endl;
  cout<<Ninc<<"all events"<<endl; 

  TH1D *phihist[33];
  TH1D *etahist[33];
  char histname1[32], histname2[32];
  for(int i=0; i<33; i++){
    sprintf(histname1, "phi%d", i);
    sprintf(histname2, "h%d phi", i);
    phihist[i] = new TH1D(histname1, histname2, round(5.608*sqrt(nevent)/100)*100, -3.2, 3.2);
  }

  for(int i=0; i<33; i++){
    sprintf(histname1, "eta%d", i);
    sprintf(histname2, "h%d eta", i);
    etahist[i] = new TH1D(histname1, histname2, round(5.608*sqrt(nevent)/100)*100, -3.0, 3.0);
  }




  n = 0;
  first = 0;



  int f = 0;
  TGraphErrors *mint[Ninc];
  TGraphErrors *hit[Ninc];

  //get layer info
  int nlayerz = 0;
  for(int i = 0; i < nfast; i++){
    cout<<i<<endl<<endl;
    cout<<nlayerz<<endl;
    if((i+1)%100000 == 0) cout<<i+1<<endl;
    tree->GetEntry(i);
    for(int d=0;d<33;d++){
      if(layer[d]==1){
	if(layerz[d]==0) nlayerz+=1;
        layerz[d] = z[d];
        layersign[d] = sign[d];
      }
    }
    if(layerz[17]!=0&&layerz[18]!=0&&layerz[19]!=0&&layerz[20]!=0&&layerz[21]!=0&&layerz[22]!=0&&layerz[23]!=0&&layerz[24]!=0&&layerz[25]!=0&&layerz[26]!=0&&layerz[27]!=0&&layerz[28]!=0&&layerz[29]!=0&&layerz[30]!=0&&layerz[31]!=0&&layerz[32]!=0)break;
  }
  nlayerz = 0;
  for(int i = 0; i < nevent; i++){
    if((i+1)%100000 == 0) cout<<i+1<<endl;
    memo->GetEntry(i);
    for(int d=0;d<33;d++){
      if(fabs(hitz[d]-layerz[d])<zwidth&&fabs(hitz[d]-layerz[d])>2.5){
	if(layerzcent[d]==0)nlayerz+=1;
	layerzcent[d] = layerz[d]+signfunc(hitz[d]-layerz[d])*2.5;
      }
    }
    if(nlayerz>=16) break;
  }

  for(int d=0;d<33;d++){
    cout<<"layer"<<d-16<<"\t"<<layerz[d]<<"\t"<<layerzcent[d]<<"\t"<<layersign[d]<<endl;
  }

  double per = TMath::ErfcInverse(0.01);//95%                                                                                                     
  double a=0;
  double errora=0;
  double v[8]={0};
  double w[8]={0};
  double phibound[16][4]={{0}};
  double etabound[2][33]={{0}};
  double generaletabound[2][4]={{0}};
  double generalphiwidth[4]={0};
  double e[16]={0};
  int first_eta[33] = {0};

  if(fixon==0){
    vector<int> *PDGID =0;
    vector<double> *Hits_kineticEnergy=0;
    vector<double> *Hits_DirectionZ=0;
    vector<double> *hitspos_X=0;
    vector<double> *hitspos_Y=0;
    vector<double> *hitspos_Z=0;
    vector<double> *hitspos_R=0;
    vector<double> *hitstime=0;
    data->SetBranchAddress("Hits_MM_hitGlobalPositionX",&hitspos_X);
    data->SetBranchAddress("Hits_MM_hitGlobalPositionY",&hitspos_Y);
    data->SetBranchAddress("Hits_MM_hitGlobalPositionZ",&hitspos_Z);
    data->SetBranchAddress("Hits_MM_hitGlobalPositionR",&hitspos_R);
    data->SetBranchAddress("Hits_MM_globalTime",&hitstime);
    data->SetBranchAddress("Hits_MM_particleEncoding" ,&PDGID); //PDGID                                                                    
    data->SetBranchAddress("Hits_MM_kineticEnergy",&Hits_kineticEnergy);
    data->SetBranchAddress("Hits_MM_hitGlobalDirectionZ",&Hits_DirectionZ);
  for(int nn=0;nn<data->GetEntries();nn++){
    data->GetEntry(nn);
    
    if((nn+1)%100000 == 0) cout<<(nn+1)/1.0/data->GetEntries()<<endl;
    
    for(int c=0;c<hitspos_X->size();c++){
      if(PDGID->at(c)==13){    
	for(int d=0;d<33;d++){
	  if(fabs(hitspos_Z->at(c)-layerz[d])<zwidth){
	    phihist[d]->Fill( atan2(hitspos_Y->at(c),hitspos_X->at(c)) );
	    etahist[d]->Fill( -log(tan(atan2(hitspos_R->at(c),hitspos_Z->at(c))/2.0)) );
	  }
	}   
      }
    }
    f = f+1;
  }      
  

  TF1 *func = new TF1("f","[0]*TMath::Erfc([1]*(x-[2]))");

  for(int l=-8;l<8;l++){
    for(int d=0;d<8;d++){
      func->SetParameters(phihist[d]->GetBinContent(phihist[d]->GetMaximumBin())/2.0,30*pow(-1,l),(2*l+1)/16.0*pi);
      func->SetParLimits(0,phihist[d]->GetBinContent(phihist[d]->GetMaximumBin())/4.0,phihist[d]->GetBinContent(phihist[d]->GetMaximumBin())/2.0);
      func->SetParLimits(2,2*l/16.0*pi,2*(l+1)/16.0*pi);
      phihist[d]->Fit("f","Q","goff",2*l/16.0*pi,2*(l+1)/16.0*pi);
      cout<<"laC_side"<<d-16<<endl;
      a = func->GetParameter(1);
      errora = func->GetParError(1);
      v[d] = func->GetParameter(2)-per/a;
      error = func->GetParError(2);
      error = sqrt(error*error+per*per*errora*errora/a/a/a/a);
      w[d] = error*error;
      cout<<func->GetChisquare()/func->GetNDF()<<endl;
      phibound[l+8][0]=TMath::Mean(8,&v[0]);
    }
  }

  for(int l=-8;l<8;l++){
    for(int d=8;d<16;d++){
      func->SetParameters(phihist[d]->GetBinContent(phihist[d]->GetMaximumBin())/2.0,-30*pow(-1,l),(2*l+1)/16.0*pi);
      func->SetParLimits(0,phihist[d]->GetBinContent(phihist[d]->GetMaximumBin())/4.0,phihist[d]->GetBinContent(phihist[d]->GetMaximumBin())/2.0);
      func->SetParLimits(2,2*l/16.0*pi,2*(l+1)/16.0*pi);
      phihist[d]->Fit("f","Q","goff",2*l/16.0*pi,2*(l+1)/16.0*pi);
      cout<<"smC_side"<<d-16<<endl;
      a = func->GetParameter(1);
      errora = func->GetParError(1);
      v[d-8] = func->GetParameter(2)-per/a;
      error =func->GetParError(2);
      error =sqrt(error*error+per*per*errora*errora/a/a/a/a);
      w[d-8] = error*error;

      cout<<func->GetChisquare()/func->GetNDF()<<endl;
      phibound[l+8][1]=TMath::Mean(8,&v[0]);
    }
  }

  for(int l=-8;l<8;l++){
    for(int d=17;d<25;d++){
      func->SetParameters(phihist[d]->GetBinContent(phihist[d]->GetMaximumBin())/2.0,-30*pow(-1,l),(2*l+1)/16.0*pi);
      func->SetParLimits(0,phihist[d]->GetBinContent(phihist[d]->GetMaximumBin())/4.0,phihist[d]->GetBinContent(phihist[d]->GetMaximumBin())/2.0);
      func->SetParLimits(2,2*l/16.0*pi,2*(l+1)/16.0*pi);
      phihist[d]->Fit("f","Q","goff",2*l/16.0*pi,2*(l+1)/16.0*pi);
      cout<<"smA_side"<<d-16<<endl;
      a = func->GetParameter(1);
      errora = func->GetParError(1);
      v[d-17] = func->GetParameter(2)-per/a;
      error =func->GetParError(2);
      error =sqrt(error*error+per*per*errora*errora/a/a/a/a);
      w[d-17] = error*error;
      cout<<func->GetChisquare()/func->GetNDF()<<endl;
      phibound[l+8][2]=TMath::Mean(8,&v[0]);
    }
  }


  for(int l=-8;l<8;l++){
    for(int d=25;d<33;d++){
      func->SetParameters(phihist[d]->GetBinContent(phihist[d]->GetMaximumBin())/2.0,30*pow(-1,l),(2*l+1)/16.0*pi);
      func->SetParLimits(0,phihist[d]->GetBinContent(phihist[d]->GetMaximumBin())/4.0,phihist[d]->GetBinContent(phihist[d]->GetMaximumBin())/2.0);
      func->SetParLimits(2,2*l/16.0*pi,2*(l+1)/16.0*pi);
      phihist[d]->Fit("f","Q","goff",2*l/16.0*pi,2*(l+1)/16.0*pi);
      cout<<"laA_side"<<d-16<<endl;
      a = func->GetParameter(1);
      errora = func->GetParError(1);
      v[d-25] = func->GetParameter(2)-per/a;
      error =func->GetParError(2);
      error =sqrt(error*error+per*per*errora*errora/a/a/a/a);
      w[d-25] = error*error;
      cout<<func->GetChisquare()/func->GetNDF()<<endl;
      phibound[l+8][3]=TMath::Mean(8,&v[0]);
    }
  }
  
 
  for(int d=0;d<16;d++){
    func->SetParameters(etahist[d]->GetBinContent(etahist[d]->GetMaximumBin())/2.0,-10,-2.8);
    func->SetParLimits(0,etahist[d]->GetBinContent(etahist[d]->GetMaximumBin())/4.0,etahist[d]->GetBinContent(etahist[d]->GetMaximumBin())/2.0);
    func->SetParLimits(2,-3.0,-2.5);
    etahist[d]->Fit("f","Q","",-3.0,-2.5);
    cout<<"C_side"<<d-16<<endl;
    a = func->GetParameter(1);
    etabound[0][d] = func->GetParameter(2)-per/a;
    cout<<func->GetChisquare()/func->GetNDF()<<endl;
   }
  for(int d=0;d<16;d++){
    func->SetParameters(etahist[d]->GetBinContent(etahist[d]->GetMaximumBin())/2.0,10,-1.24);
    func->SetParLimits(0,etahist[d]->GetBinContent(etahist[d]->GetMaximumBin())/4.0,etahist[d]->GetBinContent(etahist[d]->GetMaximumBin())/2.0);
    func->SetParLimits(2,-1.5,-1.0);
    etahist[d]->Fit("f","Q","",-1.5,-1.0);
    cout<<"C_side"<<d-16<<endl;
    a = func->GetParameter(1);
    etabound[1][d] = func->GetParameter(2)-per/a;
    cout<<func->GetChisquare()/func->GetNDF()<<endl;
  }
  for(int d=17;d<33;d++){
    func->SetParameters(etahist[d]->GetBinContent(etahist[d]->GetMaximumBin())/2.0,-10,1.24);
    func->SetParLimits(0,etahist[d]->GetBinContent(etahist[d]->GetMaximumBin())/4.0,etahist[d]->GetBinContent(etahist[d]->GetMaximumBin())/2.0);    
    func->SetParLimits(2,1.0,1.5);
    etahist[d]->Fit("f","Q","",1.0,1.5);
    cout<<"A_side"<<d-16<<endl;
    a = func->GetParameter(1);
    etabound[0][d] = func->GetParameter(2)-per/a;
    cout<<func->GetChisquare()/func->GetNDF()<<endl;
  }
  for(int d=17;d<33;d++){
    func->SetParameters(etahist[d]->GetBinContent(etahist[d]->GetMaximumBin())/2.0,10,2.8);
    func->SetParLimits(0,etahist[d]->GetBinContent(etahist[d]->GetMaximumBin())/4.0,etahist[d]->GetBinContent(etahist[d]->GetMaximumBin())/2.0);
    func->SetParLimits(2,2.5,3.0);
    etahist[d]->Fit("f","Q","",2.5,3.0);
    cout<<"A_side"<<d-16<<endl;
    a = func->GetParameter(1);
    etabound[1][d] = func->GetParameter(2)-per/a;
    cout<<func->GetChisquare()/func->GetNDF()<<endl;
  }
  
  
  cout<<"phi"<<endl<<endl;

  for(int c=0;c<4;c++){
    cout<<c<<endl;
    for(int d=0;d<16;d++){
      cout<<d<<"\t"<<phibound[d][c]<<endl;
    }
  }

  cout<<"eta"<<endl<<endl;  

  for(int d=0;d<2;d++){
    for(int c=0;c<33;c++){
      cout<<c<<"\t"<<etabound[d][c]<<endl;
    }
  }
  int generalbound = 0;
  for(int d=0;d<33;d++){
    if(d==0||d==8||d==17||d==25){
      generaletabound[0][generalbound] = etabound[0][d];
      generaletabound[1][generalbound] = etabound[1][d];
      generalbound += 1;
    }
    if(generaletabound[0][generalbound]>etabound[0][d]) generaletabound[0][generalbound] = etabound[0][d];
    if(generaletabound[1][generalbound]<etabound[1][d]) generaletabound[1][generalbound] = etabound[1][d];
  }

  for(int c=0;c<4;c++){
    for(int d=0;d<16;d++){
      if(fabs(c-1.5)>1){//large                                                                                                                        
        if(d==0)generalphiwidth[c]=fabs(phibound[d][c]-pi/4.0*floor((d-7)/2.0));
        if(fabs(phibound[d][c]-pi/4.0*floor((d-7)/2.0))<generalphiwidth[c]){
          generalphiwidth[c]=fabs(phibound[d][c]-pi/4.0*floor((d-7)/2.0));
        }
      }else{
        if(d==0)generalphiwidth[c]=fabs(phibound[d][c]-pi/4.0*floor((d-8)/2.0)-pi/8.0);
        if(fabs(phibound[d][c]-pi/4.0*floor((d-8)/2.0)-pi/8.0)<generalphiwidth[c]){
          generalphiwidth[c]=fabs(phibound[d][c]-pi/4.0*floor((d-8)/2.0)-pi/8.0);
        }
      }
    }
  }
  }

  if(fixon==1){
    generaletabound[0][0] = 0;
    generaletabound[0][1] = 0;
    //generaletabound[0][2] = 1.32748;                                                                                                                 
    //generaletabound[0][3] = 1.33956;                                                                                                                 
    generaletabound[0][2] = 1.35;
    generaletabound[0][3] = 1.35;
    generaletabound[1][0] = 0;
    generaletabound[1][1] = 0;
    //generaletabound[1][2] = 2.66647;                                                                                                                 
    //generaletabound[1][3] = 2.69886;                                                                                                                 
    generaletabound[1][2] = 2.65;
    generaletabound[1][3] = 2.65;

    generalphiwidth[0] = 0;
    generalphiwidth[1] = 0;
    generalphiwidth[2] = 0.13788;
    generalphiwidth[3] = 0.245034;
  }
  cout<<"eta"<<endl<<endl;
  for(int d=0;d<2;d++){
    for(int c=0;c<4;c++){
      cout<<c<<"\t"<<generaletabound[d][c]<<endl;
    }
  }

  cout<<"phi"<<endl<<endl;
  for(int c=0;c<4;c++){
    cout<<c<<"\t"<<generalphiwidth[c]<<endl;
  }


  
  cout <<"neff"<< f <<endl;



  //tree creation                                                               


  TFile *filed3 = new TFile(Form("%sd3.root",namecore),"recreate");
  TTree *d3 = new TTree("d3","d3");


  n = 0;
  first = 0;
  nevent = memo->GetEntries();


  //Plate                                                                      


  f=0;
  int platendata = 0;
  vector<double> hitmux(300);
  vector<double> hitmuy(300);
  vector<double> hitmuz(300);
  int nn = 0;
  int plateef;
  int platenf=0;
  int platesnf=0;
  int plateneta=0;
  int Uon = 0;
  int Von = 0;
  int Ulegit = 0;
  int Vlegit = 0;
  int Uvalid = 0;
  int Vvalid = 0;
  double URZmemory = 0;
  double VRZmemory = 0;
  int plateth=0;//through
  int platesth=0;//stereo through
  int platenetath=0;//eta through
  int platei=0;
  int platenUVin=0;
  int platenXin=0;
  int plateXpat=0;
  int plateUVpat=0;
  int platePhipat=0;
  int platePhiXpat=0;
  double platex0=0;
  double platey0=0;
  double platez0=0;
  double plater0=0;
  double platetheta0=0;
  double plateeta0=0;
  double plateentrytheta=0;
  double platephi0=0;
  double plateentryphi=0;
  double dr2d=0;
  double piphi=0;
  vector<double> plateresid2dl(4);
  vector<double> plateresid2dg(4);
  int on = 0;
  int cut;
  int hitbutcut=0;
  double phifraction=0;
  double strip_phifraction=0;
  double phimemory=0;
  double lphi=0;
  int dall = 0;
  int dstereo = 0;
  double hit_phimemory=0;
  double hit_thetamemory=0;
  double hit_xmemory=0;
  double hit_ymemory=0;
  double hit_zmemory=0;
  double hit_gx=0;
  double hit_gy=0;
  double hit_gz=0;
  double strip_phimemory=0;
  double strip_xmemory=0;
  double strip_ymemory=0;
  double strip_zmemory=0;
  double hit_dphi=0;
  double dt2d=0;
  double platetimewindow = 0;
  double estimated_phifraction=0;
  double rhat=0;
  int platesector=0;
  int platesectoron=0;
  double platechi2_ndf=0;
  double platedr0=0;
  double plateetadr0=0;
  double platetmin=0;
  int fakeornot = 0;
  double rplus = 0;
  double rminus = 0;
  int dplus = 0;
  int dminus = 0;
  double dialphi = 0;

  vector<double> vec_dr2d(300);
  vector<double> vec_zintercept(300);
  vector<double> vec_rintercept(300);
  vector<double> vec_ltheta(300);
  vector<double> vec_gtheta(300);
  vector<double> vec_ltheta3d(300);
  vector<double> vec_lphi3d(300);
  vector<double> vec_gtheta3d(300);
  vector<double> vec_gphi3d(300);
  vector<double> vec_zintercept3d(300);
  vector<double> vec_strip_phimemory(300);
  vector<double> vec_timewindow(300);
  vector<int> vec_nlayer(300);
  vector<int> vec_nUV(300);
  vector<int> vec_nX(300);
  vector<int> vec_nUVin(300);
  vector<int> vec_nXin(300);
  vector<int> vec_Xpat(300);
  vector<int> vec_UVpat(300);
  vector<int> vec_Phipat(300);
  vector<int> vec_PhiXpat(300);
  vector<int> vec_Uvalid(300);
  vector<int> vec_Vvalid(300);
  vector<double> vec_estimated_phifraction(300);
  vector<int> vec_sector(300);
  vector<int> vec_sectoron(300);
  vector<double> vec_chi2_ndf(300);
  vector<double> vec_dr0(300);
  vector<double> vec_etadr0(300);
  vector<double> vec_tmin(300);
  vector<int> vec_fakeornot(300);
  vector<double> vec_hit_zintercept3d(300);
  vector<double> vec_hit_zintercept(300);
  vector<double> vec_hit_rintercept(300);
  vector<double> vec_hit_ltheta(300);
  vector<double> vec_hit_gtheta(300);
  vector<double> vec_hit_3dtheta(300);
  vector<double> vec_hit_3dphi(300);
  vector<double> vec_hit_gtheta3d(300);
  vector<double> vec_hit_gphi3d(300);
  vector<double> vec_hit_phimemory(300);
  vector<double> vec_hit_thetamemory(300);
  vec_dr2d.erase(vec_dr2d.begin(),vec_dr2d.end());
  vec_zintercept3d.erase(vec_zintercept3d.begin(),vec_zintercept3d.end());
  vec_zintercept.erase(vec_zintercept.begin(),vec_zintercept.end());
  vec_rintercept.erase(vec_rintercept.begin(),vec_rintercept.end());
  vec_ltheta.erase(vec_ltheta.begin(),vec_ltheta.end());
  vec_gtheta.erase(vec_gtheta.begin(),vec_gtheta.end());
  vec_ltheta3d.erase(vec_ltheta3d.begin(),vec_ltheta3d.end());
  vec_lphi3d.erase(vec_lphi3d.begin(),vec_lphi3d.end());
  vec_gtheta3d.erase(vec_gtheta3d.begin(),vec_gtheta3d.end());
  vec_gphi3d.erase(vec_gphi3d.begin(),vec_gphi3d.end());
  vec_strip_phimemory.erase(vec_strip_phimemory.begin(),vec_strip_phimemory.end());
  vec_timewindow.erase(vec_timewindow.begin(),vec_timewindow.end());
  vec_nlayer.erase(vec_nlayer.begin(),vec_nlayer.end());
  vec_nUV.erase(vec_nUV.begin(),vec_nUV.end());
  vec_nX.erase(vec_nX.begin(),vec_nX.end());
  vec_nUVin.erase(vec_nUVin.begin(),vec_nUVin.end());
  vec_nXin.erase(vec_nXin.begin(),vec_nXin.end());
  vec_Xpat.erase(vec_Xpat.begin(),vec_Xpat.end());
  vec_UVpat.erase(vec_UVpat.begin(),vec_UVpat.end());
  vec_Phipat.erase(vec_Phipat.begin(),vec_Phipat.end());
  vec_PhiXpat.erase(vec_PhiXpat.begin(),vec_PhiXpat.end());
  vec_Uvalid.erase(vec_Uvalid.begin(),vec_Uvalid.end());
  vec_Vvalid.erase(vec_Vvalid.begin(),vec_Vvalid.end());
  vec_estimated_phifraction.erase(vec_estimated_phifraction.begin(),vec_estimated_phifraction.end());
  vec_sector.erase(vec_sector.begin(),vec_sector.end());
  vec_sectoron.erase(vec_sectoron.begin(),vec_sectoron.end());
  vec_chi2_ndf.erase(vec_chi2_ndf.begin(),vec_chi2_ndf.end());
  vec_dr0.erase(vec_dr0.begin(),vec_dr0.end());
  vec_etadr0.erase(vec_etadr0.begin(),vec_etadr0.end());
  vec_tmin.erase(vec_tmin.begin(),vec_tmin.end());
  vec_fakeornot.erase(vec_fakeornot.begin(),vec_fakeornot.end());
  vec_hit_zintercept.erase(vec_hit_zintercept.begin(),vec_hit_zintercept.end());
  vec_hit_rintercept.erase(vec_hit_rintercept.begin(),vec_hit_rintercept.end());
  vec_hit_ltheta.erase(vec_hit_ltheta.begin(),vec_hit_ltheta.end());
  vec_hit_gtheta.erase(vec_hit_gtheta.begin(),vec_hit_gtheta.end());
  vec_hit_3dtheta.erase(vec_hit_3dtheta.begin(),vec_hit_3dtheta.end());
  vec_hit_3dphi.erase(vec_hit_3dphi.begin(),vec_hit_3dphi.end());
  vec_hit_gtheta3d.erase(vec_hit_gtheta3d.begin(),vec_hit_gtheta3d.end());
  vec_hit_gphi3d.erase(vec_hit_gphi3d.begin(),vec_hit_gphi3d.end());
  vec_hit_phimemory.erase(vec_hit_phimemory.begin(),vec_hit_phimemory.end());
  vec_hit_thetamemory.erase(vec_hit_thetamemory.begin(),vec_hit_thetamemory.end());

  d3->Branch("ndata",&platendata);
  d3->Branch("i",&platei);
  d3->Branch("x0",&platex0);
  d3->Branch("y0",&platey0);
  d3->Branch("z0",&platez0);
  d3->Branch("r0",&plater0);
  d3->Branch("theta0",&platetheta0);
  d3->Branch("entrytheta",&plateentrytheta);
  d3->Branch("phi0",&platephi0);
  d3->Branch("entryphi",&plateentryphi);
  d3->Branch("dr2d",&vec_dr2d);
  d3->Branch("zcept",&vec_zintercept);
  d3->Branch("rcept",&vec_rintercept);
  d3->Branch("hit_zcept",&hit_zintercept);
  d3->Branch("hit_rcept",&hit_rintercept);
  d3->Branch("ltheta",&vec_ltheta);
  d3->Branch("gtheta",&vec_gtheta);
  d3->Branch("hit_ltheta",&hit_ltheta);
  d3->Branch("hit_gtheta",&hit_gtheta);
  d3->Branch("ltheta3d",&vec_ltheta3d);
  d3->Branch("lphi3d",&vec_lphi3d);
  d3->Branch("hit_ltheta3d",&hit_3dtheta);
  d3->Branch("hit_lphi3d",&hit_3dphi);
  d3->Branch("gtheta3d",&vec_gtheta3d);
  d3->Branch("gphi3d",&vec_gphi3d);
  d3->Branch("hit_gtheta3d",&hit_gtheta3d);
  d3->Branch("hit_gphi3d",&hit_gphi3d);
  d3->Branch("zcept3d",&vec_zintercept3d);
  d3->Branch("hit_posphi",&hit_phimemory);
  d3->Branch("hit_postheta",&hit_thetamemory);
  d3->Branch("strip_posphi",&vec_strip_phimemory);
  d3->Branch("timewindow",&vec_timewindow);
  d3->Branch("nth",&plateth);
  d3->Branch("nUVth",&platesth);
  d3->Branch("nXth",&platenetath);
  d3->Branch("nlayer",&vec_nlayer);
  d3->Branch("nUV",&vec_nUV);
  d3->Branch("nX",&vec_nX);
  d3->Branch("nUVin",&vec_nUVin);
  d3->Branch("nXin",&vec_nXin);
  d3->Branch("Xpat",&vec_Xpat);
  d3->Branch("UVpat",&vec_UVpat);
  d3->Branch("Phipat",&vec_Phipat);
  d3->Branch("PhiXpat",&vec_PhiXpat);
  d3->Branch("Uvalid",&vec_Uvalid);
  d3->Branch("Vvalid",&vec_Vvalid);
  d3->Branch("estimated_phifraction",&vec_estimated_phifraction);
  d3->Branch("sector",&vec_sector);
  d3->Branch("sectoron",&vec_sectoron);
  d3->Branch("chi2_ndf",&vec_chi2_ndf);
  d3->Branch("dr0",&vec_dr0);
  d3->Branch("etadr0",&vec_etadr0);
  d3->Branch("tmin",&vec_tmin);
  d3->Branch("fakeornot",&vec_fakeornot);
  cut = Ninc;

  int firsttime = 0;
  int truthcheck = 0;
  int la_sm = 1;
  int legitimacy = 0;
  int acceptance=0;
  int validity=0;
  int layeron[33] = {0};
  int allhitsectoron[2] = {0};
  int bandmemory = 0;
  int Xfirst = 0;
  int UVfirst = 0;
  int Xsecond = 0;
  int Xfirston = 0;
  int UVfirston = 0;
  int Xsecondon = 0; 

  n = 0;
  int m = 0;
  tree->GetEntry(0);
  if(shortcut==1){
    n = nevent-100;
    m = 13100;
    tree->GetEntry(13100);
  }
  while(m < nfast){

    if(checkexp==1)cout<<"n"<<n<<"i"<<i<<endl;//check1
    truthcheck = 0;//for truth

    if((n+1)%100000 == 0) cout<<(n+1)/1.0/Ninc<<endl;

    on = 0;    

    vec_dr2d.erase(vec_dr2d.begin(),vec_dr2d.end());
    vec_zintercept.erase(vec_zintercept.begin(),vec_zintercept.end());
    vec_rintercept.erase(vec_rintercept.begin(),vec_rintercept.end());
    vec_ltheta.erase(vec_ltheta.begin(),vec_ltheta.end());
    vec_gtheta.erase(vec_gtheta.begin(),vec_gtheta.end());
    vec_ltheta3d.erase(vec_ltheta3d.begin(),vec_ltheta3d.end());
    vec_lphi3d.erase(vec_lphi3d.begin(),vec_lphi3d.end());
    vec_gtheta3d.erase(vec_gtheta3d.begin(),vec_gtheta3d.end());
    vec_gphi3d.erase(vec_gphi3d.begin(),vec_gphi3d.end());
    vec_zintercept3d.erase(vec_zintercept3d.begin(),vec_zintercept3d.end());
    vec_strip_phimemory.erase(vec_strip_phimemory.begin(),vec_strip_phimemory.end());
    vec_timewindow.erase(vec_timewindow.begin(),vec_timewindow.end());
    vec_nlayer.erase(vec_nlayer.begin(),vec_nlayer.end());
    vec_nUV.erase(vec_nUV.begin(),vec_nUV.end());
    vec_nX.erase(vec_nX.begin(),vec_nX.end());
    vec_nUVin.erase(vec_nUVin.begin(),vec_nUVin.end());
    vec_nXin.erase(vec_nXin.begin(),vec_nXin.end());
    vec_Xpat.erase(vec_Xpat.begin(),vec_Xpat.end());
    vec_UVpat.erase(vec_UVpat.begin(),vec_UVpat.end());
    vec_Phipat.erase(vec_Phipat.begin(),vec_Phipat.end());
    vec_PhiXpat.erase(vec_PhiXpat.begin(),vec_PhiXpat.end());
    vec_Uvalid.erase(vec_Uvalid.begin(),vec_Uvalid.end());
    vec_Vvalid.erase(vec_Vvalid.begin(),vec_Vvalid.end());
    vec_estimated_phifraction.erase(vec_estimated_phifraction.begin(),vec_estimated_phifraction.end());
    vec_sector.erase(vec_sector.begin(),vec_sector.end());
    vec_sectoron.erase(vec_sectoron.begin(),vec_sectoron.end());
    vec_chi2_ndf.erase(vec_chi2_ndf.begin(),vec_chi2_ndf.end());
    vec_dr0.erase(vec_dr0.begin(),vec_dr0.end());
    vec_etadr0.erase(vec_etadr0.begin(),vec_etadr0.end());
    vec_tmin.erase(vec_tmin.begin(),vec_tmin.end());
    vec_fakeornot.erase(vec_fakeornot.begin(),vec_fakeornot.end());
    vec_hit_zintercept.erase(vec_hit_zintercept.begin(),vec_hit_zintercept.end());
    vec_hit_rintercept.erase(vec_hit_rintercept.begin(),vec_hit_rintercept.end());
    vec_hit_ltheta.erase(vec_hit_ltheta.begin(),vec_hit_ltheta.end());
    vec_hit_gtheta.erase(vec_hit_gtheta.begin(),vec_hit_gtheta.end());
    vec_hit_3dtheta.erase(vec_hit_3dtheta.begin(),vec_hit_3dtheta.end());
    vec_hit_3dphi.erase(vec_hit_3dphi.begin(),vec_hit_3dphi.end());
    vec_hit_gtheta3d.erase(vec_hit_gtheta3d.begin(),vec_hit_gtheta3d.end());
    vec_hit_gphi3d.erase(vec_hit_gphi3d.begin(),vec_hit_gphi3d.end());
    vec_hit_phimemory.erase(vec_hit_phimemory.begin(),vec_hit_phimemory.end());
    vec_hit_thetamemory.erase(vec_hit_thetamemory.begin(),vec_hit_thetamemory.end());

    while(n == i){//while ni start
      cout<<"i"<<i<<endl;
      platendata = ndata;
      cout<<ndata<<endl;
      firsttime = 0;//for hit tracking
      if(truthcheck == 0){
	if(checkexp==1)cout<<"truthcheck"<<endl;//check2
      	truthcheck = 1;
	memo->GetEntry(n);
	if(checkexp==1)cout<<"mueta"<<mueta<<" : "<<eta0<<endl;//datacheck
	allhitsectoron[0] = 0;
	allhitsectoron[1] = 0;
	for(int d=0;d<33;d++){layeron[d] = 0;} 
	plateth = 0;
	platesth = 0;
	acceptance = 0;
	legitimacy = 0;
	validity = 0;
	if(checkexp==1)cout<<"Vertex_x"<<vertex_n<<endl;
	if(checkexp==1)cout<<"MU +/-"<<pdg0<<endl;
	if(mueta>eta1 && mueta<eta2 && pdg0==13 && vertex_n>nvertexcrit && acceptance==0){
	  if(checkexp==1)cout<<"ACCEPT"<<endl;
	  platetheta0 = 2.0*atan(exp(-mueta));
	  plateeta0 = mueta;
	  if(mueta<0){
	    for(int d=0;d<8;d++){
	      if(fabs(hitphi0[d]-round(hitphi0[d]/pi*4.0)/4.0*pi)<generalphiwidth[0]){
		if((hiteta0[d]<generaletabound[1][0] && hiteta0[d]>1.7 ) || ( 1.5>hiteta0[d] && hiteta0[d]>generaletabound[0][0])){
		  on = 1;
		  layeron[d] = 1;
		  def[d] += 1;
		  plateth += 1;   
		  if(layersign[d]!=0){
		    platesth += 1;
		  }
		}
	      }
	    }//large C_side
	    
	    for(int d=8;d<16;d++){
	      if(fabs(hitphi0[d]-round((hitphi0[d]+pi/8.0)/pi*4.0)/4.0*pi+pi/8.0)<generalphiwidth[1]){
		if((hiteta0[d]<generaletabound[1][1] && hiteta0[d]>1.7 ) || ( 1.5>hiteta0[d] && hiteta0[d]>generaletabound[0][1])){
		  on = 1;
		  layeron[d] = 1;
		  def[d] += 1;
		  plateth += 1;
		  if(layersign[d]!=0){
		    platesth += 1;
		  }
		}
	      }
	    }//small C_side
	  }
	  
	  if(mueta>0){
	    for(int d=17;d<25;d++){
	      if(fabs(hitphi0[d]-round((hitphi0[d]+pi/8.0)/pi*4.0)/4.0*pi+pi/8.0)<generalphiwidth[2]){
		if((hiteta0[d]<generaletabound[1][2] && hiteta0[d]>1.7 ) || ( 1.5>hiteta0[d] && hiteta0[d]>generaletabound[0][2])){
		  on = 1;
		  layeron[d] = 1;
		  def[d] += 1;
		  plateth += 1;
		  if(layersign[d]!=0){
		    platesth += 1;
		  }
		}
	      }
	    }//small A_side         
	    
	    for(int d=25;d<33;d++){
	      if(fabs(hitphi0[d]-round(hitphi0[d]/pi*4.0)/4.0*pi)<generalphiwidth[3]){
		if((hiteta0[d]<generaletabound[1][3] && hiteta0[d]>1.7 ) || ( 1.5>hiteta0[d] && hiteta0[d]>generaletabound[0][3])){
		  on = 1;
		  layeron[d] = 1;
		  def[d] += 1;
		  plateth += 1;
		  if(layersign[d]!=0){
		    platesth += 1;
		  }
		}
	      }
	    }//large A_side    
	  }
	}
	
	platenetath = plateth-platesth;    
	
	for(int b=0;b<17;b++){
	  if(plateth==b){
	    dnf[b] += 1;
	    for(int h=0;h<9;h++){
	      if(platesth==h){
		dsnf[b][h] += 1;
	      }
	    }
	  }
	}      
	
	if(layeron[17]+layeron[18]+layeron[23]+layeron[24]+layeron[19]+layeron[20]+layeron[21]+layeron[22]==8){
	  allhitsectoron[0] = 1;
	}
	if(layeron[25]+layeron[26]+layeron[31]+layeron[32]+layeron[27]+layeron[28]+layeron[29]+layeron[30]==8){
	  allhitsectoron[1] = 1;
	}
      }//truthcheck    
      
      if(plateth==8){
	if(platesth==4){
	  if(checkexp==1)cout<<"allhit"<<endl;//check3
	  //data->GetEntry(n);
	  if(sector==3) sector = allhitsectoron[1];
	  platesector = sector;
	  platesectoron = allhitsectoron[sector];
	  if(mueta != eta0) alert[4]+=1;
	  phimemory = phi0;
	  platephi0 = phi0;
	  plateentryphi = 0;
	  plateeta0 = mueta;
	  platetheta0 = 2.0*atan(exp(-mueta));
	  plateentrytheta = 2.0*atan(exp(-muentryeta));
	  platex0 = x0;
	  platey0 = y0;
	  platez0 = z0;
	  plater0 = sqrt(x0*x0+y0*y0+z0*z0);
	  platei = n;
	  truth_aofz0 = platex0+(0-platez0)/cos(platetheta0)*sin(platetheta0)*cos(platephi0);
	  truth_bofz0 = platey0+(0-platez0)/cos(platetheta0)*sin(platetheta0)*sin(platephi0);
	  
	  /*	for(int c=0;c<hitspos_X->size();c++){
		if(PDGID->at(c)!=13){
		mux.push_back();
		muy.push_back();
		muz.push_back();
		}
		}*/
	  
	  for(int s=0;s<3;s++){
	    platetimewindow = 25*(s+1);
	    if(checkexp==1){
	      if(allhitsectoron[0]==1)cout<<"sm"<<endl;
	      if(allhitsectoron[1]==1)cout<<"la"<<endl;
	      cout<<"band"<<band;
	      if(sector==0)cout<<" sector sm";
	      if(sector==1)cout<<" sector la";
	      cout<<" "<<neta<<"X"<<nstereo<<"UV"<<endl;//check5   
	      cout<<"window"<<platetimewindow<<endl;//check5
	    }
	    platenf = 0;
	    platesnf = 0;
	    Uon = 0;
	    Ulegit = 0;
	    Uvalid = 0;
	    URZmemory = 0;
	    Von = 0;
	    Vlegit = 0;
	    Vvalid = 0;
	    VRZmemory = 0;
	    for(int d=0;d<33;d++){
	      if(time[d]!=0&&layeron[d]==1){
		if(time[d]!=0 && time[d] < (s+1)*25){
		  ef[d][s] += 1;
		  platenf += 1;
		  if(layersign[d]!=0){
		    platesnf += 1;
		    if(layersign[d]==1){Uon = 1; Ulegit = 1;}else{Von = 1; Vlegit = 1;} 
		  }
		}
	      }
	    }      
	    
	    if(platenf>8)alert[0]+=1;      
	    
	    plateneta = platenf - platesnf;
	    
	    //initialize 
	    estimated_phifraction = 0;
	    strip_phimemory = 0;
	    strip_xmemory = 0;
	    strip_ymemory = 0;
	    strip_zmemory = 0;
	    
	    dd = 0;dall = 0;dstereo = 0;
	    xofe.erase(xofe.begin(),xofe.end());
	    yofe.erase(yofe.begin(),yofe.end());
	    xyofe.erase(xyofe.begin(),xyofe.end());
	    dofe.erase(dofe.begin(),dofe.end());
	    x3d.erase(x3d.begin(),x3d.end());
	    y3d.erase(y3d.begin(),y3d.end());
	    z3d.erase(z3d.begin(),z3d.end());
	    d3d.erase(d3d.begin(),d3d.end());
	    dofs.erase(dofs.begin(),dofs.end());
	    rofs.erase(rofs.begin(),rofs.end());
	    signofs.erase(signofs.begin(),signofs.end());
	    platenXin = 0;
	    platenUVin = 0;
	    plateXpat = 0;
	    plateUVpat = 0;
	    platePhipat = 0;
	    platePhiXpat = 0;
	    platechi2_ndf = 0;	
	    platedr0 = 0;
	    plateetadr0 = 0;
	    platetmin = 0;
	    fakeornot = 1;	    

	    for(int c=0;c<10;c++){
	      l1ofe.at(c)=0;
	      l2ofe.at(c)=0;      
	      l3ofe.at(c)=0;
	    }
	    
	    for(int l=2;l<5;l++){
	      if(platenf-platesnf==l){
		en[l][s] += 1;
		for(int k=1;k<5;k++){
		  if(platesnf==k){
		    esn[l][k][s] += 1;

		    if(uvfit==2){
		      Ulegit = 0;
		      Vlegit = 0;
		      Uvalid = 0;
		      Vvalid = 0;
		      URZmemory = 0;
		      VRZmemory = 0;
		      for(int d=0;d<33;d++){
			if(time[d]!=0 && time[d] < 25*(s+1) && layeron[d]==1){
			  //cout<<d<<endl;
			  if(layersign[d]!=0){
			    phifraction = atan2(y[d],x[d]);
                            if(abs(d-16)>8){
                              while( abs(phifraction) > pi/8.0 ){
                                if(phifraction>0){ phifraction -= pi/4.0;
                                }else{ phifraction += pi/4.0;}
                              }
                            }else{
                              phifraction += pi/8.0;
                              while( abs(phifraction) > pi/8.0 ){
                                if(phifraction>0){ phifraction -= pi/4.0;
                                }else{ phifraction += pi/4.0; }
                              }
                            }
			    if(layersign[d]==1){
			      if(URZmemory != 0){
				cout<<URZmemory-hypot(x[d],y[d])*(cos(phifraction)+layersign[d]*sin(strip_phifraction)*tan(1.5*pi/180))/layerz[d]<<" : "<<slopewidth/2.0<<endl;
				if(fabs(hypot(x[d],y[d])*(cos(phifraction)+layersign[d]*sin(strip_phifraction)*tan(1.5*pi/180))/layerz[d]-URZmemory)<=slopewidth/2.0){
				  Ulegit = 1;
				  Uvalid = 1;
				} 
			      }
			      URZmemory = hypot(x[d],y[d])*(cos(phifraction)+layersign[d]*sin(strip_phifraction)*tan(1.5*pi/180))/layerz[d];
			    }
			    if(layersign[d]==-1){
			      if(VRZmemory != 0){
				if(fabs(hypot(x[d],y[d])*(cos(phifraction)+layersign[d]*sin(strip_phifraction)*tan(1.5*pi/180))/layerz[d]-VRZmemory)<=slopewidth/2.0){
				  Vlegit = 1;
				  Vvalid = 1;
				}
			      }
                              VRZmemory = hypot(x[d],y[d])*(cos(phifraction)+layersign[d]*sin(strip_phifraction)*tan(1.5*pi/180))/layerz[d];
                            }
			  }
			}
		      }
		      if(Ulegit==0&&Vlegit==0){
			Ulegit = 1;
			Vlegit = 1;
		      }
		    }

		    for(int d=0;d<33;d++){
		      //eta
		      if(time[d]!=0 && time[d] < 25*(s+1) && layeron[d]==1){
			if(layersign[d]==0){
			  if(dd>3){cout<<"e"<<dd<<"layer"<<d-16<<endl;}else{
			    xofe.push_back(layerz[d]);
			    dofe.push_back(d);
			    if(abs(d-16)>8){
			      la_sm = 2;
			    }else{
			      la_sm = 1;
			    }
			    platedr0 += dr0[d]*dr0[d];
			    plateetadr0 += dr0[d]*dr0[d];
			    platetmin += time[d];
			    if(fabs(dr0[d])<drcut){
			      fakeornot = 0;
			    }
			    if(fabs(dr0[d])<drcutin){
			      platenXin += 1;
			    }
			    yofe.push_back(sqrt(x[d]*x[d]+y[d]*y[d]));//as r
			    xyofe.push_back(layerz[d]*sqrt(x[d]*x[d]+y[d]*y[d]));
			    //strip_phimemory += atan2(y[d],x[d]);
			    strip_xmemory += x[d];
			    strip_ymemory += y[d];
			    strip_zmemory += z[d];
			    //exofe.at(dd) = 5.;
			    //eyofe.at(dd) = 0.435;
			    dialphi = atan2(y[d],x[d]);
			    dialphi = round(dialphi/pi*8.0)*pi/8.0;
			    if(uvfit==0||Uon==0||Von==0||Ulegit==0||Vlegit==0){
			      x3d.push_back(x[d]);
			      y3d.push_back(y[d]);
			      z3d.push_back(layerz[d]);
			      d3d.push_back(d);
			      lphi = atan2(y[d],x[d]);
			      lphi = round(lphi/pi*8.0)*pi/8.0+pi/2.0;
			      //cout<<"layer"<<d<<"\tlphi"<<(lphi-pi/2.0)/pi*8.0<<endl;
			      l1ofe.at(dall) = cos(lphi);
			      l2ofe.at(dall) = sin(lphi);
			      l3ofe.at(dall) = 0;
			      dall += 1;//for xuv/uv fit		  
			    }
			    dd += 1;//for eta
			  }
			}else{
			  if((Ulegit!=0&&layersign[d]==1)||(Vlegit!=0&&layersign[d]==-1)){
			    platedr0 += dr0[d]*dr0[d];
			    platetmin += time[d];
			    if(fabs(dr0[d])<drcutin){
			      platenUVin += 1;
			    }
			    x3d.push_back(x[d]);
			    y3d.push_back(y[d]);
			    z3d.push_back(layerz[d]);
			    d3d.push_back(d);
			    dofs.push_back(d);//stereo
			    phifraction = atan2(y[d],x[d]);
			    if(la_sm==2){
			      while( abs(phifraction) > pi/8.0 ){
				if(phifraction>0){ phifraction -= pi/4.0;
				}else{ phifraction += pi/4.0;}
			      }
			    }else{
			      phifraction += pi/8.0;
			      while( abs(phifraction) > pi/8.0 ){
				if(phifraction>0){ phifraction -= pi/4.0;
				}else{ phifraction += pi/4.0; }
			      }
			    }
			    rofs.push_back(hypot(x[d],y[d])*(cos(phifraction)+layersign[d]*sin(strip_phifraction)*tan(1.5*pi/180)));//stereo//degpoint
			    signofs.push_back(layersign[d]);
			    lphi = atan2(y[d],x[d]);
			    lphi = round(lphi/pi*8.0)*pi/8.0+pi/2.0+layersign[d]*1.5/180.0*pi;//degpoint
			    l1ofe.at(dall) = cos(lphi);
			    l2ofe.at(dall) = sin(lphi);
			    l3ofe.at(dall) = 0;
			    dstereo += 1;
			    dall += 1;
			  }
			}
		      }
		    }
		    if(dd!=l){alert[2]+=1;}
		    if(((uvfit==2&&Ulegit!=0&&Vlegit!=0&&Uon!=0&&Von!=0)&&(dall!=k))||((uvfit==1&&Uon!=0&&Von!=0)&&(dall!=k))||((uvfit==0||Uon==0||Von==0)&&(dall!=l+k))){alert[3]+=1;}
		    if(dd==l){
		      plateetadr0 /= 1.0*l;
		      platedr0 /= 1.0*(l+dstereo);
		      platetmin /= 1.0*(l+dstereo);
		      if(analytic==0){
			if(robust==0){
			  mint[n] = new TGraphErrors(l,&xofe[0],&yofe[0],&exofe[0],&eyofe[0]);
			  linear->SetParameters(2*atan(exp(-(eta1+eta2)/2.0)),0);
			  if(mueta>0){		      
			    linear->SetParLimits(0,0.,pi/2.0);
			  }else{
			    linear->SetParLimits(0,pi/2.0,pi);
			  }
			  linear->FixParameter(1,0);
			  mint[n]->Fit("lin","Q","",-8000,8000);
			  if(fabs(linear->GetParameter(0)/pi*2.0-round(linear->GetParameter(0)/pi*2.0))<0.00001)alert[7]+=1;
			  gtheta = linear->GetParameter(0);
			  linear->ReleaseParameter(1);
			  mint[n]->Fit("lin","Q","",-8000,8000);
			  if(fabs(linear->GetParameter(0)/pi*2.0-round(linear->GetParameter(0)/pi*2.0))<0.00001)alert[7]+=1;
			  ltheta = linear->GetParameter(0);
			  if(l>2){platechi2_ndf = linear->GetChisquare()/1.0/(l-2);}else{platechi2_ndf = linear->GetChisquare()/0.1;}
			  esdelta[l][k][s] += fabs(gtheta-ltheta);
			  edelta[l][s] += fabs(gtheta-ltheta);
			  zintercept = linear->GetParameter(1);
			  rintercept = -zintercept*tan(ltheta);
			}
			if(robust==1){
			  for(int mr=0;mr<l;mr++){
			    xfit->ReleaseParameter(mr*2);
			    xfit->ReleaseParameter(mr*2+1);
			    xfit->SetParameter(mr*2,Form("x[%d]",mr),xofe.at(mr),1,-8000,8000);
			    xfit->SetParameter(mr*2+1,Form("y[%d]",mr),yofe.at(mr),1,0,6000);
			    xfit->FixParameter(mr*2);
			    xfit->FixParameter(mr*2+1);
			  }
			  for(int mr=l;mr<4;mr++){
			    xfit->ReleaseParameter(mr*2);
			    xfit->ReleaseParameter(mr*2+1);
			    xfit->SetParameter(mr*2,Form("z[%d]",mr),0,1,-8000,8000);
			    xfit->SetParameter(mr*2+1,Form("r[%d]",mr),0,1,0,6000);
			    xfit->FixParameter(mr*2);
			    xfit->FixParameter(mr*2+1);
			  }
			  if(mueta>0){
			    xfit->SetParameter(8,"a",tan(2*atan(exp(-(eta1+eta2)/2.0))),0.5,0.0,1.0);
			  }else{
			    xfit->SetParameter(8,"a",tan(2*atan(exp(-(eta1+eta2)/2.0))),0.5,-1.0,0.0);
			  }
			  xfit->SetParameter(9,"b",0,300,-2000,2000);
			  xfit->FixParameter(9);
			  xfit->ExecuteCommand("MIGRAD",0,0);
			  if(fabs(atan(xfit->GetParameter(8))/pi*2.0-round(atan(xfit->GetParameter(8))/pi*2.0))<0.00001)alert[7]+=1;
                          //gtheta = atan(xfit->GetParameter(8));

			  xfit->ReleaseParameter(9);
			  xfit->SetParameter(9,"b",0,300,-2000,2000);
			  xfit->ExecuteCommand("MIGRAD",0,0);
			  /*while(fabs(fabs(minuit->GetParameter(61))-pi)<0.00000001){
			    piphi = minuit->GetParameter(61);
			    cout<<piphi<<endl;
			    minuit->SetParameter(61,"phi",-piphi*0.99,0.1,-pi,pi);
			    minuit->ExecuteCommand("MIGRAD",0,0);
			    //cout<<platephi0<<" : "<<minuit->GetParameter(61)<<endl;//for phi
			    }
			  minuit->ExecuteCommand("MIGRAD",0,0);
			  while(fabs(fabs(minuit->GetParameter(61))-pi)<0.00000001){
			    piphi = minuit->GetParameter(61);
			    cout<<piphi<<endl;
			    minuit->SetParameter(61,"phi",-piphi*0.99,0.1,-pi,pi);
			    minuit->ExecuteCommand("MIGRAD",0,0);
			    //cout<<platephi0<<" : "<<minuit->GetParameter(61)<<endl;//for phi*/                                         
			  if(fabs(atan(xfit->GetParameter(8))/pi*2.0-round(atan(xfit->GetParameter(8))/pi*2.0))<0.00001)alert[7]+=1;
                          ltheta = atan(xfit->GetParameter(8));
                          //if(l>2){platechi2_ndf = xfit->GetChisquare()/1.0/(l-2);}else{platechi2_ndf = xfit->GetChisquare()/0.1;}
                          esdelta[l][k][s] += fabs(gtheta-ltheta);
                          edelta[l][s] += fabs(gtheta-ltheta);
			  rintercept = xfit->GetParameter(9);
                          zintercept = -rintercept/tan(ltheta);

			  xfit->ReleaseParameter(8);
			  xfit->ReleaseParameter(9);
			  gtheta = atan(TMath::Mean(l,&xyofe[0])/(TMath::RMS(l,&xofe[0])*TMath::RMS(l,&xofe[0])*(l-1)/l+TMath::Mean(l,&xofe[0])*TMath::Mean(l,&xofe[0])));

			}	   
		      }else{
			if(checkexp==1){
			  for(int ppp=0;ppp<l;ppp++){
			    cout<<xofe[ppp]<<"\t"<<yofe[ppp]<<endl;
			  }
			}
			//cout<<TMath::RMS(l,&xofe[0])<<endl;
			if(analytic==1){
			  gtheta = atan(TMath::Mean(l,&xyofe[0])/(TMath::RMS(l,&xofe[0])*TMath::RMS(l,&xofe[0])*(l-1)/l+TMath::Mean(l,&xofe[0])*TMath::Mean(l,&xofe[0])));
			  ltheta = atan((TMath::Mean(l,&xyofe[0])-TMath::Mean(l,&xofe[0])*TMath::Mean(l,&yofe[0]))/(TMath::RMS(l,&xofe[0])*TMath::RMS(l,&xofe[0])*(l-1)/l));
			  if(checkexp==1){cout<<"gthetacheck"<<gtheta<<endl;
			    cout<<"lthetacheck"<<ltheta<<endl;}
			  if(l>2){
			    platechi2_ndf = TMath::RMS(l,&yofe[0])*TMath::RMS(l,&yofe[0])*(l-1)/l-(TMath::Mean(l,&xyofe[0])-TMath::Mean(l,&xofe[0])*TMath::Mean(l,&yofe[0]))*(TMath::Mean(l,&xyofe[0])-TMath::Mean(l,&xofe[0])*TMath::Mean(l,&yofe[0]))/TMath::RMS(l,&xofe[0])/TMath::RMS(l,&xofe[0])/(l-1)*l/(l-2);
			  }else{
			    platechi2_ndf = TMath::RMS(l,&yofe[0])*TMath::RMS(l,&yofe[0])*(l-1)/l-(TMath::Mean(l,&xyofe[0])-TMath::Mean(l,&xofe[0])*TMath::Mean(l,&yofe[0]))*(TMath::Mean(l,&xyofe[0])-TMath::Mean(l,&xofe[0])*TMath::Mean(l,&yofe[0]))/TMath::RMS(l,&xofe[0])/TMath::RMS(l,&xofe[0])/(l-1)*l/0.1;
			  }
			  esdelta[l][k][s] += fabs(gtheta-ltheta);
			  edelta[l][s] += fabs(gtheta-ltheta);
			  rintercept = TMath::Mean(l,&yofe[0])-TMath::Mean(l,&xofe[0])*tan(ltheta);
			  zintercept = -rintercept/tan(ltheta);
			}
			if(analytic==2){
			  mint[n] = new TGraphErrors(l,&xofe[0],&yofe[0],&exofe[0],&eyofe[0]);
			  linearl->SetParameters(tan(2*atan(exp(-(eta1+eta2)/2.0))),0);
			  if(mueta>0){
			    linearl->SetParLimits(0,0.,1000.);
			  }else{
			    linearl->SetParLimits(0,-1000.,0.);
			  }
			  linearl->FixParameter(1,0);
			  mint[n]->Fit("linl","QC","",-8000,8000);
			  //if(fabs(linearl->GetParameter(0)/pi*2.0-round(linearl->GetParameter(0)/pi*2.0))<0.00001)alert[7]+=1;
			  gtheta = atan(linearl->GetParameter(0));
			  linearl->ReleaseParameter(1);
			  mint[n]->Fit("linl","QC","",-8000,8000);
			  //if(fabs(linearl->GetParameter(0)/pi*2.0-round(linearl->GetParameter(0)/pi*2.0))<0.00001)alert[7]+=1;
			  ltheta = atan(linearl->GetParameter(0));
			  if(l>2){platechi2_ndf = linearl->GetChisquare()/1.0/(l-2);}else{platechi2_ndf = linearl->GetChisquare()/0.1;}
			  esdelta[l][k][s] += fabs(gtheta-ltheta);
			  edelta[l][s] += fabs(gtheta-ltheta);
			  zintercept = linearl->GetParameter(1);
			  rintercept = -zintercept*tan(ltheta);
			}
		      }

		      plateresid2dl.erase(plateresid2dl.begin(),plateresid2dl.end());
		      plateresid2dg.erase(plateresid2dg.begin(),plateresid2dg.end());
		      for(int mr=0;mr<dd;mr++){
			resid = yofe.at(mr)-tan(ltheta)*xofe.at(mr)+tan(ltheta)*zintercept;
			residhist[4]->Fill(/*xofe.at(mr),*/resid);
			plateresid2dl.push_back(resid);
			resid = yofe.at(mr)-tan(gtheta)*xofe.at(mr);
			residhist[5]->Fill(/*xofe.at(mr),*/resid);
			plateresid2dg.push_back(resid);
		      }
		      if(uvfit==0||Uon==0||Von==0||Ulegit==0||Vlegit==0){
			for(int pp=0;pp<dofs.size();pp++){
			  rhat = tan(ltheta)*layerz[dofs.at(pp)]-tan(ltheta)*zintercept;
			  resid = rofs.at(pp)-rhat;
			  estimated_phifraction += signofs.at(pp)*atan(resid/rhat/tan(1.5/180*pi));
			}
			if(dofs.size()!=0)estimated_phifraction /= dofs.size();
		      }else{
			rplus = 0;
			rminus = 0;
			dplus = 0; dminus = 0;
			for(int pp=0;pp<dofs.size();pp++){
			  if(signofs.at(pp)==1){
			    dplus += 1;
			    rplus += rofs.at(pp); 
			  }
			  if(signofs.at(pp)==-1){
                            dminus += 1;
                            rminus += rofs.at(pp);
                          }
			}
			rplus /= 1.0*dplus;
			rminus /= 1.0*dminus;
			rhat = (rplus+rminus)/2.0;
			resid = rplus-rhat;
			estimated_phifraction = atan(resid/rhat/tan(1.5/180*pi));
		      }

		      if(platesnf!=dofs.size())alert[5] += 1;
		      strip_xmemory /= l;
		      strip_ymemory /= l;
		      strip_zmemory /= l;
		      strip_phimemory = atan2(strip_ymemory,strip_xmemory);

		    }
		    
		    if(checkexp==1)cout<<"minuit"<<endl;//check6
		    plateXpat = SectorXpat(dofe.size(),&dofe[0]);
		    platePhipat = SectorAllpat(d3d.size(),&d3d[0]);
		    platePhiXpat = SectorXpat(d3d.size(),&d3d[0]);
		    plateUVpat = SectorUVpat(d3d.size(),&d3d[0]);
		    cout<<plateXpat<<"\t"<<plateUVpat<<"\t"<<platePhipat<<"\t"<<platePhiXpat<<endl;
		    
		    for(int mr=0;mr<dall;mr++){
		      minuit->ReleaseParameter(mr*6);
		      minuit->ReleaseParameter(mr*6+1);
		      minuit->ReleaseParameter(mr*6+2);
		      minuit->ReleaseParameter(mr*6+3);
		      minuit->ReleaseParameter(mr*6+4);
		      minuit->ReleaseParameter(mr*6+5);
		      minuit->SetParameter(mr*6,Form("x[%d]",mr),x3d.at(mr),1,-4500,4500);
		      minuit->SetParameter(mr*6+1,Form("y[%d]",mr),y3d.at(mr),1,-4500,4500);
		      minuit->SetParameter(mr*6+2,Form("z[%d]",mr),z3d.at(mr),1,-8000,8000);
		      minuit->SetParameter(mr*6+3,Form("l1[%d]",mr),l1ofe.at(mr),1,-100.,100.);
		      minuit->SetParameter(mr*6+4,Form("l2[%d]",mr),l2ofe.at(mr),1,-100.,100.);
		      minuit->SetParameter(mr*6+5,Form("l3[%d]",mr),l3ofe.at(mr),1,-100.,100.);
		      minuit->FixParameter(mr*6);
		      minuit->FixParameter(mr*6+1);
		      minuit->FixParameter(mr*6+2);
		      minuit->FixParameter(mr*6+3);
		      minuit->FixParameter(mr*6+4);
		      minuit->FixParameter(mr*6+5);
		    }
		    
		    for(int mr=dall;mr<10;mr++){
		      minuit->ReleaseParameter(mr*6);
		      minuit->ReleaseParameter(mr*6+1);
		      minuit->ReleaseParameter(mr*6+2);
		      minuit->ReleaseParameter(mr*6+3);
		      minuit->ReleaseParameter(mr*6+4);
		      minuit->ReleaseParameter(mr*6+5);
		      minuit->SetParameter(mr*6,Form("x[%d]",mr),0,1,-4500,4500);
		      minuit->SetParameter(mr*6+1,Form("y[%d]",mr),0,1,-4500,4500);
		      minuit->SetParameter(mr*6+2,Form("z[%d]",mr),0,1,-8000,8000);
		      minuit->SetParameter(mr*6+3,Form("l1[%d]",mr),l1ofe.at(mr),1,-100.,100.);
		      minuit->SetParameter(mr*6+4,Form("l2[%d]",mr),l2ofe.at(mr),1,-100.,100.);
		      minuit->SetParameter(mr*6+5,Form("l3[%d]",mr),l3ofe.at(mr),1,-100.,100.);
		      minuit->FixParameter(mr*6);
		      minuit->FixParameter(mr*6+1);
		      minuit->FixParameter(mr*6+2);
		      minuit->FixParameter(mr*6+3);
		      minuit->FixParameter(mr*6+4);
		      minuit->FixParameter(mr*6+5);
		    }
		    minuit->SetParameter(62,"a",0,300,-2000,2000);
		    minuit->SetParameter(63,"b",0,300,-2000,2000);
		    minuit->FixParameter(62);
		    minuit->FixParameter(63);

		    if(l3d==1){
		      minuit->SetParameter(64,"c",zintercept,300,-2000,2000);
		      if(z3d.at(0)>0){
			if(uvfit==0||Uon==0||Von==0||Ulegit==0||Vlegit==0){
			  minuit->SetParameter(60,"theta",atan(tan(ltheta)/cos(estimated_phifraction)),0.1,0.,pi/2.0);
			}else{
			  minuit->SetParameter(60,"theta",atan(signfunc(mueta)*rhat/(layerz[la_sm*8+9]+layerz[la_sm*8+16])*2.0)/cos(estimated_phifraction),0.1,0.,pi/2.0);
			}
		      }else{
			if(uvfit==0||Uon==0||Von==0||Ulegit==0||Vlegit==0){
			  minuit->SetParameter(60,"theta",atan(tan(ltheta)/cos(estimated_phifraction)),0.1,pi/2.0,pi);
			}else{
			  minuit->SetParameter(60,"theta",atan(signfunc(mueta)*rhat/(layerz[la_sm*8+9]+layerz[la_sm*8+16])*2.0)/cos(estimated_phifraction),0.1,pi/2.0,pi);
			}
		      }
		      minuit->SetParameter(61,"phi",dialphi+estimated_phifraction,0.1,-pi,pi);
		      //cout<<platephi0<<" : "<<dialphi+estimated_phifraction<<endl;//for phi
		      minuit->FixParameter(60);
		      minuit->FixParameter(61);
		      minuit->ExecuteCommand("MIGRAD",0,0);
		      minuit->ReleaseParameter(60);
		      minuit->ReleaseParameter(61);
		      minuit->ExecuteCommand("MIGRAD",0,0);
		      //cout<<platephi0<<" : "<<minuit->GetParameter(61)<<endl;//for phi    
		      while(fabs(fabs(minuit->GetParameter(61))-pi)<0.00000001){
			piphi = minuit->GetParameter(61);
			cout<<piphi<<endl;
			minuit->SetParameter(61,"phi",-piphi*0.99,0.1,-pi,pi);
			minuit->ExecuteCommand("MIGRAD",0,0);
			//cout<<platephi0<<" : "<<minuit->GetParameter(61)<<endl;//for phi  
		      }
		      minuit->ExecuteCommand("MIGRAD",0,0);
		      while(fabs(fabs(minuit->GetParameter(61))-pi)<0.00000001){
			piphi = minuit->GetParameter(61);
			cout<<piphi<<endl;
			minuit->SetParameter(61,"phi",-piphi*0.99,0.1,-pi,pi);
			minuit->ExecuteCommand("MIGRAD",0,0);
			//cout<<platephi0<<" : "<<minuit->GetParameter(61)<<endl;//for phi      
		      }
		      minuit->ExecuteCommand("MIGRAD",0,0);
		      while(fabs(fabs(minuit->GetParameter(61))-pi)<0.00000001){
			piphi = minuit->GetParameter(61);
			cout<<piphi<<endl;
			minuit->SetParameter(61,"phi",-piphi*0.99,0.1,-pi,pi);
			minuit->ExecuteCommand("MIGRAD",0,0);
			//cout<<platephi0<<" : "<<minuit->GetParameter(61)<<endl;//for phi \
			
		      }
		      
		      
		      ltheta3d = minuit->GetParameter(60);
		      lphi3d = minuit->GetParameter(61);
		      //cout<<platephi0<<" : "<<lphi3d<<endl;//for phi    
		      zintercept3d = minuit->GetParameter(64);
		      aofz0 = intersecX(0,0,zintercept3d,gtheta3d,gphi3d,0);
		      bofz0 = intersecY(0,0,zintercept3d,gtheta3d,gphi3d,0);
		    }

		    minuit->SetParameter(64,"c",0,300,-2000,2000);
		    minuit->FixParameter(64);
		    if(z3d.at(0)>0){
		      if(uvfit==0||Uon==0||Von==0||Ulegit==0||Vlegit==0){
			minuit->SetParameter(60,"theta",atan(tan(gtheta)/cos(estimated_phifraction)),0.1,0.,pi/2.0);
		      }else{
			minuit->SetParameter(60,"theta",atan(signfunc(mueta)*rhat/(layerz[la_sm*8+9]+layerz[la_sm*8+16])*2.0)/cos(estimated_phifraction),0.1,0.,pi/2.0);
		      }
		    }else{
		      if(uvfit==0||Uon==0||Von==0||Ulegit==0||Vlegit==0){
			minuit->SetParameter(60,"theta",atan(tan(gtheta)/cos(estimated_phifraction)),0.1,pi/2.0,pi);
		      }else{
			minuit->SetParameter(60,"theta",atan(signfunc(mueta)*rhat/(layerz[la_sm*8+9]+layerz[la_sm*8+16])*2.0)/cos(estimated_phifraction),0.1,pi/2.0,pi);
		      }
		    }
		    minuit->SetParameter(61,"phi",dialphi+estimated_phifraction,0.1,-pi,pi);
		    //cout<<platephi0<<" : "<<dialphi+estimated_phifraction<<endl;//for phi
   
		    minuit->ExecuteCommand("MIGRAD",0,0);
		    while(fabs(fabs(minuit->GetParameter(61))-pi)<0.00000001){
                      piphi = minuit->GetParameter(61);
                      cout<<piphi<<endl;
                      minuit->SetParameter(61,"phi",-piphi*0.99,0.1,-pi,pi);
                      minuit->ExecuteCommand("MIGRAD",0,0);
                      //cout<<platephi0<<" : "<<minuit->GetParameter(61)<<endl;//for phi                         
		    }
                    minuit->ExecuteCommand("MIGRAD",0,0);
		    while(fabs(fabs(minuit->GetParameter(61))-pi)<0.00000001){
                      piphi = minuit->GetParameter(61);
                      cout<<piphi<<endl;
                      minuit->SetParameter(61,"phi",-piphi*0.99,0.1,-pi,pi);
                      minuit->ExecuteCommand("MIGRAD",0,0);
                      //cout<<platephi0<<" : "<<minuit->GetParameter(61)<<endl;//for phi 		      			   
		    }
		    minuit->ExecuteCommand("MIGRAD",0,0);
		    while(fabs(fabs(minuit->GetParameter(61))-pi)<0.00000001){
                      piphi = minuit->GetParameter(61);
                      cout<<piphi<<endl;
                      minuit->SetParameter(61,"phi",-piphi*0.99,0.1,-pi,pi);
                      minuit->ExecuteCommand("MIGRAD",0,0);
                      //cout<<platephi0<<" : "<<minuit->GetParameter(61)<<endl;//for phi     \
                                                                                              
                    }
		    minuit->ExecuteCommand("MIGRAD",0,0);
		    while(fabs(fabs(minuit->GetParameter(61))-pi)<0.00000001){
                      piphi = minuit->GetParameter(61);
                      cout<<piphi<<endl;
                      minuit->SetParameter(61,"phi",-piphi*0.99,0.1,-pi,pi);
                      minuit->ExecuteCommand("MIGRAD",0,0);
                      //cout<<platephi0<<" : "<<minuit->GetParameter(61)<<endl;//for phi     \
                                                                                              
                    }
		    minuit->ReleaseParameter(62);
		    minuit->ReleaseParameter(63);                  
		    minuit->ReleaseParameter(64);
		    
		    gtheta3d = minuit->GetParameter(60);
		    gphi3d = minuit->GetParameter(61);
		    if(firsttime==0){//start of hit
		      if(checkexp==1)cout<<"hitstart"<<endl;//check4
		      firsttime = 1;
		      //data->GetEntry(n);		    
		      //initialize
		      hit_phimemory = 0;
		      hit_thetamemory = 0;
		      hit_xmemory = 0;
		      hit_ymemory = 0;
		      hit_zmemory = 0;
		      hit_dphi = 0;
		    
		      for(int d=0;d<33;d++){
			mur[d] = 0;
			mut[d] = 0;
		      }
				      		      
		      //initialize                                                
		      int ddd = 0;
		      
		      hitxofe.erase(hitxofe.begin(),hitxofe.end());
		      hityofe.erase(hityofe.begin(),hityofe.end());
		      hitxyofe.erase(hitxyofe.begin(),hitxyofe.end());
		      hit3dxofe.erase(hit3dxofe.begin(),hit3dxofe.end());
		      hit3dyofe.erase(hit3dyofe.begin(),hit3dyofe.end());
		      hit3dzofe.erase(hit3dzofe.begin(),hit3dzofe.end());
		      
		      //for hitpoints
		      for(int d=0;d<33;d++){

			if(signfunc(d-16)==signfunc(mueta)&&(abs(d-16)>8)+1==la_sm){
			  if(hitt0[d]==0)alert[6]+=1;
			  phifraction = atan2(hity0[d],hitx0[d]);
			  if(la_sm==2){	
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
			    hitxofe.push_back(hitz0[d]);
			    hityofe.push_back(hypot(hitx0[d],hity0[d])*cos(phifraction));
			    hitxyofe.push_back(hitz0[d]*hypot(hitx0[d],hity0[d])*cos(phifraction));
			  }
			  hit_xmemory += hitx0[d];
			  hit_ymemory += hity0[d];
			  hit_zmemory += hitz0[d];
			  mur[d] = hypot(hitx0[d],hity0[d])*(cos(phifraction)+layersign[d]*sin(phifraction)*tan(1.5*pi/180));
			  mut[d] = hitt0[d];

			  hit3dxofe.push_back(hitx0[d]);
			  hit3dyofe.push_back(hity0[d]);
			  hit3dzofe.push_back(hitz0[d]);
			  ddd += 1;
			} 
			
		      }

		      if(hitxofe.size()!=4)alert[6]+=1;
		      if(ddd!=8)alert[6]+=1;
		      if(ddd!=0){
                        hit_xmemory /= ddd;
                        hit_ymemory /= ddd;
                        hit_zmemory /= ddd;
                        hit_phimemory = atan2(hit_ymemory,hit_xmemory);
                        hit_thetamemory = atan2(hypot(hit_xmemory,hit_ymemory),hit_zmemory);
                      }		      

		      //for 2dtrack
		      if(hitxofe.size()>1){
			if(analytic==0){
			  hit[n] = new TGraphErrors(hitxofe.size(),&hitxofe[0],&hityofe[0],&hitexofe[0],&hiteyofe[0]);
			  linear->SetParameters(2*atan(exp(-(eta1+eta2)/2.0)),0);
			  linear->SetParLimits(0,0.,pi);
			  linear->FixParameter(1,0);
			  hit[n]->Fit("lin","Q","",-8000,8000);
			  hit_gtheta = linear->GetParameter(0);
			  linear->ReleaseParameter(1);
			  hit[n]->Fit("lin","Q","",-8000,8000);
			  hit_ltheta = linear->GetParameter(0);
			  hit_zintercept = linear->GetParameter(1);
			  hit_rintercept = -hit_zintercept*tan(hit_ltheta);
			  if(Xtruth==1){
                            if(la_sm==1)hit_gtheta = atan(tan(hit_ltheta)+hit_rintercept/(layerz[20]+layerz[21])*2.0);
                            if(la_sm==2)hit_gtheta = atan(tan(hit_ltheta)+hit_rintercept/(layerz[28]+layerz[29])*2.0);
                          }
			  intercepthist[1]->Fill(hit_zintercept,hit_ltheta-hit_gtheta);
			  intercepthist[3]->Fill(hit_ltheta-hit_gtheta,ltheta-gtheta);
			}else{
			  hit_gtheta = atan(TMath::Mean(4,&hitxyofe[0])/(TMath::RMS(4,&hitxofe[0])*TMath::RMS(4,&hitxofe[0])*(4-1)/4+TMath::Mean(4,&hitxofe[0])*TMath::Mean(4,&hitxofe[0])));
			  hit_ltheta = atan((TMath::Mean(4,&hitxyofe[0])-TMath::Mean(4,&hitxofe[0])*TMath::Mean(4,&hityofe[0]))/(TMath::RMS(4,&hitxofe[0])*TMath::RMS(4,&hitxofe[0])*(4-1)/4));
			  if(checkexp==1){cout<<"hit_gthetacheck"<<hit_gtheta<<endl;
			    cout<<"hit_lthetacheck"<<hit_ltheta<<endl;}
			  
			  //platechi2_ndf = TMath::RMS(4,&hityofe[0])*TMath::RMS(4,&hityofe[0])*(4-1)/4-(TMath::Mean(4,&hitxyofe[0])-TMath::Mean(4,&hitxofe[0])*TMath::Mean(4,&hityofe[0]))*(TMath::Mean(4,&hitxyofe[0])-TMath::Mean(4,&hitxofe[0])*TMath::Mean(4,&hityofe[0]))/TMath::RMS(4,&hitxofe[0])/TMath::RMS(4,&hitxofe[0])/(4-1)*4/(4-2);
			  
			  hit_rintercept = TMath::Mean(4,&hityofe[0])-TMath::Mean(4,&hitxofe[0])*tan(hit_ltheta);
			  hit_zintercept = -hit_rintercept/tan(hit_ltheta);
			  if(Xtruth==1){
                            if(la_sm==1)hit_gtheta = atan(tan(hit_ltheta)+hit_rintercept/(layerz[20]+layerz[21])*2.0);
                            if(la_sm==2)hit_gtheta = atan(tan(hit_ltheta)+hit_rintercept/(layerz[28]+layerz[29])*2.0);
                          }
			  intercepthist[1]->Fill(hit_zintercept,hit_ltheta-hit_gtheta);
			  intercepthist[3]->Fill(hit_ltheta-hit_gtheta,ltheta-gtheta);
			}
			//for 3dtrack
			for(int mr=0;mr<ddd;mr++){
			  track->ReleaseParameter(mr*3);
			  track->ReleaseParameter(mr*3+1);
			  track->ReleaseParameter(mr*3+2);
			  track->SetParameter(mr*3,Form("x[%d]",mr),hit3dxofe.at(mr),1,-4500,4500);
			  track->SetParameter(mr*3+1,Form("y[%d]",mr),hit3dyofe.at(mr),1,-4500,4500);
			  track->SetParameter(mr*3+2,Form("z[%d]",mr),hit3dzofe.at(mr),1,-8000,8000);
			  track->FixParameter(mr*3);
			  track->FixParameter(mr*3+1);
			  track->FixParameter(mr*3+2);
			}
			for(int mr=ddd;mr<10;mr++){
			  track->ReleaseParameter(mr*3);
			  track->ReleaseParameter(mr*3+1);
			  track->ReleaseParameter(mr*3+2);
			  track->SetParameter(mr*3,Form("x[%d]",mr),0,1,-4500,4500);
			  track->SetParameter(mr*3+1,Form("y[%d]",mr),0,1,-4500,4500);
			  track->SetParameter(mr*3+2,Form("z[%d]",mr),0,1,-8000,8000);
			  track->FixParameter(mr*3);
			  track->FixParameter(mr*3+1);
			  track->FixParameter(mr*3+2);
			}
			
			track->SetParameter(32,"a",truth_aofz0,100,-2000,2000);
			track->SetParameter(33,"b",truth_bofz0,100,-2000,2000);
			if(hit3dzofe.at(0)>0){track->SetParameter(30,"theta",platetheta0,0.1,0.,pi/2.0);}else{
			  track->SetParameter(30,"theta",platetheta0,0.1,pi/2.0,pi);}
			track->SetParameter(31,"phi",platephi0,0.1,-pi,pi);
			track->FixParameter(30);
			track->FixParameter(31);
			track->ExecuteCommand("MIGRAD",0,0);
			track->ReleaseParameter(30);
			track->ReleaseParameter(31);
			track->ExecuteCommand("MIGRAD",0,0);
			track->ExecuteCommand("MIGRAD",0,0);
			track->ExecuteCommand("MIGRAD",0,0);
			track->ExecuteCommand("MIGRAD",0,0);
			
			hit_aofz0 = track->GetParameter(32);
			hit_bofz0 = track->GetParameter(33);
			hit_3dtheta = track->GetParameter(30);
			hit_3dphi = track->GetParameter(31);
			//cout<<platephi0<<" : "<<dialphi+estimated_phifraction<<" : "<<hit_3dphi<<endl;//for phi   
			intercepthist[0]->Fill(hit_aofz0-truth_aofz0,hit_bofz0-truth_bofz0);
			/*cout<<"hit_aofz0"<<hit_aofz0<<endl;
			  cout<<"hit_bofz0"<<hit_bofz0<<endl;
			  cout<<hit_3dtheta<<"\t"<<hit_ltheta<<endl;
			  cout<<hit_3dphi<<"\t"<<phimemory<<endl;*/
			
			//for 2dtrack resid
			for(int mr=0;mr<hitxofe.size();mr++){
			  resid = hityofe.at(mr)-tan(hit_ltheta)*hitxofe.at(mr)+tan(hit_ltheta)*hit_zintercept;
			  residhist[0]->Fill(/*hitxofe.at(mr),*/resid);
			  resid = hityofe.at(mr)-tan(hit_gtheta)*hitxofe.at(mr);
			  residhist[1]->Fill(/*hitxofe.at(mr),*/resid);
			}
			
			//for 3dtrack resid(distance)
			for(int mr=0;mr<ddd;mr++){
			  resid = sqrt(distance1(hit_aofz0, hit_bofz0, 0, sin(hit_3dtheta)*cos(hit_3dphi), sin(hit_3dtheta)*sin(hit_3dphi), cos(hit_3dtheta), hit3dxofe.at(mr), hit3dyofe.at(mr), hit3dzofe.at(mr)));	  
			  residhist[2]->Fill(resid);
			}
			
			resid = sqrt(distance1(hit_aofz0, hit_bofz0, 0, sin(hit_3dtheta)*cos(hit_3dphi), sin(hit_3dtheta)*sin(hit_3dphi), cos(hit_3dtheta), x0, y0, z0));
			
			residhist[3]->Fill(resid);
			
		      }
		    }//last of hit
		    if(checkexp==1){
		      cout<<platendata<<"th DATA Entry"<<platei-ninc[platendata]<<endl;//datacheck
		      cout<<"mueta"<<plateeta0<<" : "<<eta0<<endl;
		    }
		    if(hit3dxofe.size()!=0){
		      hit_gz = TMath::Mean(z3d.size(),&z3d[0]);
		      if(UVtruth==1){
			if(la_sm==1)hit_gz = (layerz[20]+layerz[21])/2.0;
			if(la_sm==2)hit_gz = (layerz[28]+layerz[29])/2.0;
		      }
		      hit_gx = interseclX(hit3dxofe.at(0),hit3dyofe.at(0),hit3dzofe.at(0),hit3dxofe.at(hit3dxofe.size()-1),hit3dyofe.at(hit3dyofe.size()-1),hit3dzofe.at(hit3dzofe.size()-1),hit_gz);
		      hit_gy = interseclY(hit3dxofe.at(0),hit3dyofe.at(0),hit3dzofe.at(0),hit3dxofe.at(hit3dxofe.size()-1),hit3dyofe.at(hit3dyofe.size()-1),hit3dzofe.at(hit3dzofe.size()-1),hit_gz);
		      hit_gphi3d = atan2(hit_gy,hit_gx);
		      hit_gtheta3d = atan2(hypot(hit_gx,hit_gy),hit_gz);
		    }else{
		      hit_gphi3d = 0;hit_gtheta3d=0;alert[8]+=1;
		      cout<<"SIZE0error "<<platendata<<"th DATA Entry"<<platei-ninc[platendata]<<endl;
		      cout<<"mueta"<<plateeta0<<" : "<<eta0<<endl;
		    }
		    
		    cout<<platei<<endl;
                    intercepthist[3]->Fill(aofz0-hit_aofz0,bofz0-hit_bofz0);
                    cout<<"aofz0-hit_aofz0"<<fabs(aofz0-hit_aofz0)<<endl;
                    cout<<"bofz0-hit_bofz0"<<fabs(bofz0-hit_bofz0)<<endl;
                    cout<<fabs(gtheta3d-hit_gtheta3d)<<endl;
                    cout<<fabs(gphi3d-hit_gphi3d)<<endl;		  
		    cout<<"lthetacheck"<<fabs(ltheta-hit_ltheta)<<endl;
                    cout<<"gthetacheck"<<fabs(gtheta-hit_gtheta)<<endl;
		    delta[0] += ltheta - hit_ltheta;
		    //cout<<"hit_gtheta-eta"<<hit_gtheta-atan(tan(2*atan(exp(-eta)))*cos(phi))<<endl;
		    
		    delta[1] += gtheta - hit_gtheta;
		    nn += 1;
		    esdelta[l][k][s] += fabs(gtheta-ltheta);
		    hit_esdelta[l][k][s] += fabs(hit_gtheta-hit_ltheta);
		    sigma_sltheta[l][k][s] += (ltheta-hit_ltheta)*(ltheta-hit_ltheta);
		    sigma_sgtheta[l][k][s] += (gtheta-hit_gtheta)*(gtheta-hit_gtheta);
		    edelta[l][s] += fabs(gtheta-ltheta);
		    hit_edelta[l][s] += fabs(hit_gtheta-hit_ltheta);
		    sigma_ltheta[l][s] += (ltheta-hit_ltheta)*(ltheta-hit_ltheta);
		    sigma_gtheta[l][s] += (gtheta-hit_gtheta)*(gtheta-hit_gtheta);
		    
		    //for dr2d dt2d
		    dr2d = 0;
		    dt2d = 0;
		    
		    for(int mr=0;mr<dd;mr++){
		      for(int d=0;d<33;d++){
			if(fabs(xofe.at(mr)-layerz[d])<zwidth){
			  dr2d += (yofe.at(mr)-tan(hit_ltheta)*(layerz[d]-hit_zintercept))*(yofe.at(mr)-tan(hit_ltheta)*(layerz[d]-hit_zintercept));
			  dt2d += (time[d]-mut[d])*(time[d]-mut[d]);
			  if(d==16)alert[1]+=1;
			}
		      }
		    }
		    if(dd!=0){
		      dr2d = sqrt(dr2d/dd);
		      dt2d = sqrt(dt2d/dd); 
		    }
		    
		    
		    //d3->Fill();
		    on = 1;
		    f = f+1;
		    //cout<<platephi0<<" : "<<dialphi+estimated_phifraction<<" : "<<hit_3dphi<<endl;//for phi
		    vec_dr2d.push_back(dr2d);
		    vec_zintercept.push_back(zintercept);
		    vec_rintercept.push_back(rintercept);
		    vec_hit_zintercept.push_back(hit_zintercept);
		    vec_hit_rintercept.push_back(hit_rintercept);
		    vec_ltheta.push_back(ltheta);
		    vec_gtheta.push_back(gtheta);
		    vec_hit_ltheta.push_back(hit_ltheta);
		    vec_hit_gtheta.push_back(hit_gtheta);
		    vec_ltheta3d.push_back(ltheta3d);
		    vec_lphi3d.push_back(lphi3d);
		    vec_hit_3dtheta.push_back(hit_3dtheta);
		    vec_hit_3dphi.push_back(hit_3dphi);
		    vec_gtheta3d.push_back(gtheta3d);
		    vec_gphi3d.push_back(gphi3d);
		    vec_hit_gtheta3d.push_back(hit_gtheta3d);
		    vec_hit_gphi3d.push_back(hit_gphi3d);
		    vec_zintercept3d.push_back(zintercept3d);
		    vec_hit_phimemory.push_back(hit_phimemory);
		    vec_hit_thetamemory.push_back(hit_thetamemory);
		    vec_strip_phimemory.push_back(strip_phimemory);
		    vec_timewindow.push_back(platetimewindow);
		    vec_nlayer.push_back(platenf);
		    vec_nUV.push_back(platesnf);
		    vec_nX.push_back(plateneta);
		    vec_nUVin.push_back(platenUVin);
                    vec_nXin.push_back(platenXin);
		    vec_Xpat.push_back(plateXpat);
		    vec_UVpat.push_back(plateUVpat);
		    vec_Phipat.push_back(platePhipat);
		    vec_PhiXpat.push_back(platePhiXpat);
		    vec_Uvalid.push_back(Uvalid);
		    vec_Vvalid.push_back(Vvalid);
		    vec_estimated_phifraction.push_back(estimated_phifraction);
		    vec_sector.push_back(platesector);
		    vec_sectoron.push_back(platesectoron);
		    vec_chi2_ndf.push_back(platechi2_ndf);
		    vec_dr0.push_back(platedr0);
		    vec_etadr0.push_back(plateetadr0);
		    vec_tmin.push_back(platetmin);
		    vec_fakeornot.push_back(fakeornot);
		    
		  }
		}
	      }
	    }
	  }
	}
      }//last of s
      
      
      if(on==1){
	cut -= 1;
      }
      if(nhit!=0 && on ==0){
	hitbutcut += 1; 
      }
      
      m += 1;
      if(m < nfast){
	tree->GetEntry(m);
	cout<<"entry"<<n<<"nfast"<<m<<" : "<<ndata<<" : "<<platendata<<endl;
	if(ndata!=platendata){
	  cout<<"Change!!"<<endl;
	}
      }
      if(m == nfast){
	n = 0;
      }
    }//while ni last
    if(vec_fakeornot.size()>0){
      d3->Fill();
    }
    
    n += 1;
    //f = f+1;
  }      
  d3->Write();
  filed3->Close();  

  //cout<<"ltheta-hit_ltheta"<<delta[0]/nn;
  //cout<<"hit_gtheta-entryeta"<<delta[1]/nn; 
  const int tractnum = 5;
  int tract[tractnum][2]={{0}};
  tract[0][0]=2;tract[0][1]=1;
  tract[1][0]=2;tract[1][1]=2;
  tract[2][0]=3;tract[2][1]=2;
  tract[3][0]=3;tract[3][1]=3;
  tract[4][0]=4;tract[4][1]=4;

  TH1D *xuvlhist[3];
  TH1D *xuvghist[3];

  for(int s=0; s<3; s++){
    xuvlhist[s] = new TH1D(Form("xuvldist%d",s+1),Form("%dns",(s+1)*25),5,-0.5,4.5);
    xuvghist[s] = new TH1D(Form("xuvgdist%d",s+1),Form("%dns",(s+1)*25),5,-0.5,4.5);
  }

  TCanvas *c1 = new TCanvas("c1","c1");
  c1->Divide(3,1);
  TCanvas *c2 = new TCanvas("c2","c2");
  c2->Divide(3,1);

  char binname[200];
  const char *binnamead;
  double p=0;  

  cout <<"neff"<< f <<endl;
  cout <<"cut event"<< cut/1.0/Ninc <<endl;
  cout <<"hit but cut event"<< hitbutcut/1.0/Ninc <<endl;
  for(int s=0;s<3;s++){

    cout <<"time < "<<(s+1)*25<<"ns"<<endl;
    for(int l=1;l<5;l++){
      if(en[l][s]!=0){
        cout <<"eta"<<l<<"層"<<endl;
        cout <<edelta[l][s]/en[l][s]<<endl;
        cout <<hit_edelta[l][s]/en[l][s]<<endl;
        cout <<sqrt(sigma_ltheta[l][s])/en[l][s]<<endl;
        cout <<sqrt(sigma_gtheta[l][s])/en[l][s]<<endl;
      }
      for(int k=1;k<5;k++){
        if(esn[l][k][s]!=0){
          cout <<"stereo"<<k<<"層"<<endl;
          cout <<esdelta[l][k][s]/esn[l][k][s]<<endl;
          cout <<hit_esdelta[l][k][s]/esn[l][k][s]<<endl;
          cout <<sqrt(sigma_sltheta[l][k][s])/esn[l][k][s]<<endl;
          cout <<sqrt(sigma_sgtheta[l][k][s])/esn[l][k][s]<<endl;
        }
      }
    }


    for(int j=0;j<5;j++){
      p = sqrt(sigma_sltheta[tract[j][0]][tract[j][1]][s])/1.0/esn[tract[j][0]][tract[j][1]][s];
      if(esn[tract[j][0]][tract[j][1]][s]!=0){xuvlhist[s]->SetBinContent(j+1,p);}else{xuvlhist[s]->SetBinContent(j+1,0);}
      //xuvlhist[s]->SetBinError(j+1,sqrt(p*(1-p)/dsnf[8][4]));
      sprintf(binname, "%dX,%dUV",tract[j][0],tract[j][1]);
      binnamead = binname;
      xuvlhist[s]->GetXaxis()->SetBinLabel(j+1,binnamead);
    }
  
    for(int j=0;j<5;j++){
      p = sqrt(sigma_sgtheta[tract[j][0]][tract[j][1]][s])/1.0/esn[tract[j][0]][tract[j][1]][s];
      if(esn[tract[j][0]][tract[j][1]][s]!=0){xuvghist[s]->SetBinContent(j+1,p);}else{xuvghist[s]->SetBinContent(j+1,0);}
      //xuvlhist[s]->SetBinError(j+1,sqrt(p*(1-p)/dsnf[8][4]));                                                                 
      sprintf(binname, "%dX,%dUV",tract[j][0],tract[j][1]);
      binnamead = binname;
      xuvghist[s]->GetXaxis()->SetBinLabel(j+1,binnamead);
    }


    c1->cd(s+1);
    gStyle->SetOptStat(000000000);
    xuvlhist[s]->GetXaxis()->SetTitle("CT");    
    xuvlhist[s]->GetXaxis()->SetTitleOffset(1);
    xuvlhist[s]->GetXaxis()->SetLabelSize(0.06);
    xuvlhist[s]->GetXaxis()->SetTitleSize(0.04);
    xuvlhist[s]->GetYaxis()->SetTitle("sigma_ltheta");
    xuvlhist[s]->GetYaxis()->SetTitleOffset(1.2);
    xuvlhist[s]->GetYaxis()->SetTitleSize(0.04);
    xuvlhist[s]->Draw("");
                                                                               
    c2->cd(s+1);
    gStyle->SetOptStat(000000000);
    xuvghist[s]->GetXaxis()->SetTitle("CT");
    xuvghist[s]->GetXaxis()->SetTitleOffset(1);
    xuvghist[s]->GetXaxis()->SetLabelSize(0.06);
    xuvghist[s]->GetXaxis()->SetTitleSize(0.04);
    xuvghist[s]->GetYaxis()->SetTitle("sigma_gtheta");
    xuvghist[s]->GetYaxis()->SetTitleOffset(1.2);
    xuvghist[s]->GetYaxis()->SetTitleSize(0.04);                                                      
    xuvghist[s]->Draw("");                                                                   
  

   
  }

  TCanvas *c[10];
  c[0] = new TCanvas("resid_hit_ltheta");
  residhist[0]->Draw();
  c[1] = new TCanvas("resid_hit_gtheta");
  residhist[1]->Draw();
  c[2] = new TCanvas("hitintercept_vertex");
  intercepthist[0]->Draw();
  c[3] = new TCanvas("hitintercept_dtheta");
  intercepthist[1]->Draw();
  c[4] = new TCanvas("resid_3dltheta");
  residhist[2]->Draw();
  c[5] = new TCanvas("vertex_to_track");
  residhist[3]->Draw();
  c[6] = new TCanvas("resid_ltheta");
  residhist[4]->Draw();
  c[7] = new TCanvas("resid_gtheta");
  residhist[5]->Draw();
  c[8] = new TCanvas("intercept_vertex");
  intercepthist[2]->Draw();
  c[9] = new TCanvas("hitdtheta_dtheta");
  intercepthist[3]->Draw();

  for(int mm=0;mm<10;mm++){
    c[mm]->Print(Form("resid_ana%d.pdf",mm+1)); 
  }

  cout<<"error\t";
  for(int c=0;c<checknum;c++){
    cout<<alert[c]<<"\t";
  }
  cout<<endl;

}

