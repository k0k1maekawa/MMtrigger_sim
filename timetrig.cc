#include "/home/maekawa/atlasstyle-00-03-05/AtlasStyle.C"
#include "/home/maekawa/atlasstyle-00-03-05/AtlasLabels.C"
#include<stdlib.h> 
#include<stdio.h>
#include<math.h>
#include<iostream>
#include<time.h>
#include<vector>
#include<algorithm>
#include<functional>
#include<iomanip>
#include<string>
#include<fstream>
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
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TFitter.h"
#include <TVector3.h>
#include <TLorentzVector.h>
#include "TLegend.h"
#include "Math/Math.h"
#include "Math/ProbFuncMathCore.h"//CDF            
#include "Math/QuantFuncMathCore.h"//Quantile    
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

int signfunc(double x){
  int sgn=0;
  if(x>0) sgn=1;
  if(x<0) sgn=-1;
  return sgn;
}

int Dialturn(int d,int mu){
  int ans = 0;
  int Z = 16;
  if(d%2==1){
    if(d>=mu){
      ans = d-mu+1;
    }else{
      ans = d-mu+1+Z;
    }
  }
  if(d%2==0){
    if(d>=mu){
      ans = d-mu+2;
    }else{
      ans = d-mu+2+Z;
    }
  }
  if(d%2!=mu%2)ans = 0;
  if(d==0)ans = 0;
  return ans;
}


void timetrig(){
  SetAtlasStyle();
  int nevent;
  int nclean = 18300;

  int n; //nstrip and so on
  char namecore[300]="";



  //constant----------------------------
  const double pi = TMath::Pi();
  const double alpha = 2.0*ROOT::Math::normal_cdf(-1);

  //Overall
  const int option = 1;//0:whatever you like 
  //1:trigger 2:test in 3X3UV&50ns 3:eff for each timewindow 4:fakerate 5:etadependence 6:etaslopewidth optimization 7:eventdebug

  //etadep
  double eta1 = 1.35;
  double eta2 = 2.65;
  double etapitch = 0.2;
  //selection
  const int biason = 1;
  const int Ron = 0;
  const int timeon = 0;
  const int randomon = 0;
  //loop
  int etadep = 0;
  //const int BCloop = 0;
  int widthloop = 0;
  int xuvloop = 0;
  int timeloop = 0;
  //eraser
  const int mmfeon = 1;//8 hits from 4MMFE8 selection is implemented here
  //1 means from high R, -1 means from low R, -2 means select by time (time selection is not actual but for debugging)
  const int duplicityerase = 1;//erasing of duplicity between track segments (eliminate all the same track segment)
  int dthetacut = 0;//doing dtheta cut or not
  const int combinationerase = 1;//for recognition , not real. 1 is default.
  //other
  int threeX3UVdefault = 0;
  int fourX4UVdefault = 0;
  const int truthcheck = 1;
  const int bg = atoi(getenv("BG"));
  const int MU = atoi(getenv("MU"));
  int nvertexcrit = 0;
  if(bg==1) nvertexcrit = 100000;
  const int fixframe=1;
  const int plan=1;
  //BCID
  int BCstart = atoi(getenv("BCstart"));
  //BCstart = 0;    
  const int backcheck = 0;//view until nth BC after
  const int BCIDmethod = 0;//0:1stBC,1:2ndBC,2:last BC
  //view
  int viewSL[2]={0};
  viewSL[0] = 1;
  viewSL[1] = 1;
  int viewACC = 0;//1:see 1sector, 2:see both
  int viewMuDial = 0;
  //debug
  int dontroot = 0;
  int eventdebug = 0;
  int fortest = 0;
  int etadeplog = 0;

  //1:trigger 2:test in 3X3UV&50ns 3:eff for each timewindow 4:fakerate 5:etadependence 6:etaslopewidth optimization 7:eventdebug    
  if(option==1){etadep=0;widthloop=0;xuvloop=0;timeloop=0;threeX3UVdefault=0;viewSL[0]=1;viewSL[1]=1;viewACC=0;viewMuDial=0;}
  if(option==2){etadep=0;widthloop=0;xuvloop=0;timeloop=0;threeX3UVdefault=1;viewSL[0]=1;dontroot=1;}
  if(option==3){etadep=0;widthloop=0;xuvloop=1;timeloop=1;viewSL[0]=1;viewSL[1]=1;viewACC=2;viewMuDial=1;}
  if(option==4){etadep=0;widthloop=0;xuvloop=1;timeloop=0;viewSL[0]=0;viewSL[1]=1;viewACC=2;viewMuDial=1;}
  if(option==5){etadep=1;widthloop=0;xuvloop=1;timeloop=0;viewSL[0]=1;viewSL[1]=1;viewACC=2;viewMuDial=1;}
  if(option==6){etadep=0;widthloop=1;xuvloop=1;timeloop=0;viewSL[0]=1;viewSL[1]=1;viewACC=2;viewMuDial=1;}
  if(option==7){etadep=0;widthloop=0;xuvloop=0;timeloop=0;fourX4UVdefault=1;viewSL[0]=1;viewSL[1]=1;viewACC=0;viewMuDial=0;dontroot=1;eventdebug=1;}
  //declare property
  if(etadep==1){
    cout<<"eta "<<eta1<<"-"<<eta2<<endl;
    cout<<"etadependence"<<endl;
  }else{
    cout<<"eta 1.35-2.65"<<endl;
    cout<<"overall"<<endl;
  }
  if(biason==1)cout<<"bias select"<<endl;
  if(biason==-1)cout<<"antibias select"<<endl;
  if(Ron==1)cout<<"large R select"<<endl;
  if(Ron==-1)cout<<"small R select"<<endl;
  if(timeon==1)cout<<"timing select"<<endl;
  if(randomon==1)cout<<"random select"<<endl;
  if(widthloop==1)cout<<"co-width loop"<<endl;
  if(xuvloop==1)cout<<"xuv loop"<<endl;
  if(timeloop==1)cout<<"timewindow loop"<<endl;
  if(mmfeon==1)cout<<"mmfe large R selection"<<endl;
  if(mmfeon==-1)cout<<"mmfe small R selection"<<endl;
  if(mmfeon==-2)cout<<"mmfe time selection"<<endl;
  if(duplicityerase==1)cout<<"duplicity erasing"<<endl;
  if(dthetacut!=0)cout<<"dtheta cut ";
  if(dthetacut==1)cout<<"analytically"<<endl;
  if(dthetacut==2)cout<<"taking average"<<endl;
  if(combinationerase==1)cout<<"combination erasing"<<endl;
  if(threeX3UVdefault==1)cout<<"using default"<<endl;
  if(fourX4UVdefault==1)cout<<"using 4X4UVdefault"<<endl;
  if(truthcheck==1)cout<<"checking acceptance"<<endl;
  cout<<BCstart<<"BCstart"<<endl;
  if(backcheck!=0)cout<<"checking succeeding tracks for "<<backcheck<<"BC"<<endl;
  if(fixframe==1)cout<<"fixedframe"<<endl;
  if(plan==1){cout<<"Plan A"<<endl;}else{cout<<"Plan B"<<endl;}
  if(bg==0){cout<<"Single Muon"<<endl;}else{cout<<"with BG"<<endl;}
  if(viewSL[0]=1){cout<<"Seeing Small Sector"<<endl;}else{cout<<"Ignoring Small Sector"<<endl;}
  if(viewSL[1]=1){cout<<"Seeing Large Sector"<<endl;}else{cout<<"Ignoring Large Sector"<<endl;}
  if(viewACC==0)cout<<"vewing even Not Accept region tracks"<<endl;
  if(viewACC==1)cout<<"vewing Small and Large for Accept region tracks"<<endl;
  if(viewACC==1)cout<<"vewing Small or Lagre for Accept region tracks"<<endl;
  if(viewMuDial==1){cout<<"view Only Muon hit Phi-Wedge"<<endl;}else{cout<<"view all Phi-Wedge"<<endl;}
  int alertnum = 10;
  int alert[alertnum] = {0};
  int defbb = 0;
  int deftt = 0;
  int defaa = 0;
  int defcc = 0;
  double error;
  const int express=0;
  //const double slopewidth=0.00075;
  double slopewidth = atof(getenv("ETASLOPEWIDTH"));
  //cout<<"slopewidth"<<slopewidth<<endl;
  const double stereowidth=0.004;
  double window=50;
  int nX=3;
  int nUV=3;
  double phiwidth[2];
  phiwidth[0] = 0.13788;
  phiwidth[1] = 0.245034;//no nboth
  //  double max[10] = {0};
  //  double min[10] = {0};
  double test = 0;
  double zero[1000] = {0};

  vector<double> etavec(300);
  etavec.erase(etavec.begin(),etavec.end());

  if(etadep==1){
    for(double etadd=eta1;etadd<eta2;etadd=etadd+etapitch){
      if(etadd+etapitch<=1.5||etadd>=1.7){
	etavec.push_back(etadd+etapitch/2.0);
      }
    }
  }else{
    etavec.push_back(2.0);
    etapitch = 1.3;
  }
  cout<<"etapitch "<<etapitch<<endl;

  vector<double> loopvec(300);
  loopvec.erase(loopvec.begin(),loopvec.end());

  if(widthloop==1){
    defbb = 2;
    if(fixframe==1)defbb = 4;
    loopvec.push_back(0.00025);
    loopvec.push_back(0.0005);
    loopvec.push_back(0.00075);
    loopvec.push_back(0.001);
    loopvec.push_back(0.0015);
    loopvec.push_back(0.003);
    loopvec.push_back(0.006);
    loopvec.push_back(0.012);
  }else{
    //loopvec.push_back(atof(getenv("ETASLOPEWIDTH")));
    if(fixframe==0)loopvec.push_back(0.00075);
    if(fixframe==1&&plan==2)loopvec.push_back(0.000375);
    if(fixframe==1&&plan==1&&MU==0)loopvec.push_back(0.0015);
    if(fixframe==1&&plan==1&&MU==80)loopvec.push_back(0.0015);
    //    if(fixframe==1&&plan==1&&MU==160)loopvec.push_back(0.00075);
    if(fixframe==1&&plan==1&&MU==160)loopvec.push_back(0.0015); 
  }

  vector<double> timevec(300);
  timevec.erase(timevec.begin(),timevec.end());
  if(timeloop==1){
    deftt = 1;
    timevec.push_back(25);
    timevec.push_back(50);
    timevec.push_back(75);
  }else{
    if(widthloop+xuvloop+etadep==0){
      if(timeon==1&&threeX3UVdefault==0){
	timevec.push_back(75);
      }else{
	timevec.push_back(50);
      }
    }else{
      timevec.push_back(50);
    }
  }
 
  double drcut = atoi(getenv("DRCUT"));
  //drcut = 100;
  cout<<"drcut"<<drcut<<endl;
  //cout <<"Please Enter the Filename before \".root\""<<endl;
  //cin >> namecore;
  strcpy(namecore,getenv("NAMECORE"));
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

  TFile *file = new TFile(Form("%sxyz.root",namefordata));
  TTree *tree = (TTree*)file->Get("xyz");
  TFile *filememo = new TFile(Form("%smemo.root",namefordata));
  TTree *memo = (TTree*)filememo->Get("memo");

  //Hist
  int tractnum = 6;
  int tract[tractnum][2]={{0}};
  tract[0][0]=2;tract[0][1]=1;
  tract[1][0]=2;tract[1][1]=2;
  tract[2][0]=3;tract[2][1]=2;
  tract[3][0]=3;tract[3][1]=3;
  tract[4][0]=4;tract[4][1]=3;
  tract[5][0]=4;tract[5][1]=4;
  if(xuvloop == 0){
    tractnum = 1;
    if(etadep+widthloop+timeloop==0&&threeX3UVdefault==0&&fourX4UVdefault==0){
      tract[0][0]=2;tract[0][1]=1;
    }else{
      tract[0][0]=3;tract[0][1]=3;
      if(fourX4UVdefault==1){
	tract[0][0]=4;tract[0][1]=4;
      }
    }
  }
  if(fortest == 1){
    tractnum = 2;
    tract[0][0]=2;tract[0][1]=2;
    tract[1][0]=3;tract[1][1]=3;
    tract[2][0]=4;tract[2][1]=4;
  }
  //graph declare           
  TH1D *mudial4MMFE8hist[2][4];   
  TH1D *bandhist[tractnum];
  TH1D *drdist[tractnum];  
  TH1D *drdistresol[tractnum];
  TGraph *xuv[tractnum];
  TGraph *xuvresol[tractnum];
  TGraph *xuvall;
  TGraph *xuvallresol;
  TGraphErrors *nseff[tractnum];
  TGraphErrors *nsresoleff[tractnum];
  TGraph *nseffall;
  TGraphErrors *nseffetadep[tractnum];
  TGraphErrors *nsresoleffetadep[tractnum];
  TGraph *nseffetadepall;
  TGraph *nsresoleffetadepall;

  //Memory
  double truthrate[etavec.size()][loopvec.size()][tractnum][timevec.size()]={{{{0}}}};
  double trackrate[etavec.size()][loopvec.size()][tractnum][timevec.size()]={{{{0}}}};
  double fakerate[etavec.size()][loopvec.size()][tractnum][timevec.size()]={{{{0}}}};
  double eff[etavec.size()][loopvec.size()][tractnum][timevec.size()]={{{{0}}}};
  double resoleff[etavec.size()][loopvec.size()][tractnum][timevec.size()]={{{{0}}}};
  double effcc[etavec.size()][loopvec.size()][tractnum][timevec.size()]={{{{0}}}};
  double effError[etavec.size()][loopvec.size()][tractnum][timevec.size()]={{{{0}}}};
  double resoleffError[etavec.size()][loopvec.size()][tractnum][timevec.size()]={{{{0}}}};
  double effccError[etavec.size()][loopvec.size()][tractnum][timevec.size()]={{{{0}}}};
  int Ntruth[etavec.size()][loopvec.size()][tractnum][timevec.size()]={{{{0}}}};

  //Treesetting
  int ndata = 0;
  double x = 0;
  double y = 0;
  double z = 0;
  double pror = 0;
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
  double dr = 0;
  double dr0 = 0;
  int mudial = 0;
  int xcombi = 0;
  int uvcombi = 0;
  int vmm = 0;
  int mmfe = 0;
  int updown = 0;
  int c[8] = {0};
  int checkmmfe = 0;
  int band = 0;
  int bands = 0;
  int banddetect = 0;
  int nbands = 0;
  int dialmemory = 0;
  int mmfecount = 0;
  int mmfeRank = 0;
  double etastart = 0;
  double etastop = 0;
  int sector=0;
  tree->SetBranchAddress("ndata",&ndata);
  tree->SetBranchAddress("x",&x);
  tree->SetBranchAddress("y",&y);
  tree->SetBranchAddress("z",&z);
  tree->SetBranchAddress("pror",&pror);
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
  tree->SetBranchAddress("dr",&dr);
  tree->SetBranchAddress("mudial",&mudial);
  tree->SetBranchAddress("dr0",&dr0);
  tree->SetBranchAddress("vmm",&vmm);
  tree->SetBranchAddress("mmfe",&mmfe);
  tree->SetBranchAddress("updown",&updown);

  double hiteta0[33] = {0};
  double hitphi0[33] = {0};
  int vertex_n=100;

  memo->SetBranchAddress("hiteta0",hiteta0);
  memo->SetBranchAddress("hitphi0",hitphi0);
  memo->SetBranchAddress("nvertex",&vertex_n);

  nevent = tree->GetEntries();

  cout<<nevent<<" events"<<endl; 
  n = 0;




  //tree creation
  if(widthloop+xuvloop+etadep+timeloop+dontroot!=0)strcpy(namecore,"loop");
  TFile *filetag = new TFile(Form("%stag.root",namecore),"recreate");
  TTree *xyztag = new TTree("tag","tag");
  TFile *filesel = new TFile(Form("%ssel.root",namecore),"recreate");
  TTree *sel = new TTree("sel","sel");

  //Plate     
  int platendata = 0;            
  double platemueta = 0;
  double platemuentryeta = 0;
  double platex = 0;
  double platey = 0;
  double platez = 0;
  double platepror = 0;
  double platemux = 0;
  double platemuy = 0;
  double platemuz = 0;
  double platemut = 0;
  double platemux0 = 0;
  double platemuy0 = 0;
  double platemuz0 = 0;
  double platemut0 = 0;
  double platelx = 0;
  double plately = 0;
  int platei = 0;
  int platej = 0;
  int platek = 0;
  double platet = 0;
  int platepos = 0;
  float platecharge = 0;
  int platelayer = 0;
  int platedial = 0;
  int platesign = 0;
  double plateeta = 0;
  double platedr = 0;
  double platedr0 = 0;
  double platemudial = 0;
  int ffirst[33] = {0};
  int hiton[33] = {0};
  int hitsectoron[2] = {0};
  int allhitsectoron[2] = {0};
  double hiteta = 0;
  double hitphi = 0;
  int layeron[33] = {0};
  int sectoron[2] = {0};
  vector<double> selmueta(300);
  vector<double> selmuentryeta(300);
  vector<double> seltime(300);
  vector<double> selx(300);
  vector<double> sely(300);
  vector<double> selz(300);
  vector<double> selpror(300);
  vector<double> selslope(300);//most important
  vector<double> selmux(300);
  vector<double> selmuy(300);
  vector<double> selmuz(300);
  vector<double> selmut(300);
  vector<double> selmux0(300);
  vector<double> selmuy0(300);
  vector<double> selmuz0(300);
  vector<double> selmut0(300);
  vector<double> sellx(300);
  vector<double> selly(300);
  vector<int> seli(300);
  vector<int> selj(300);
  vector<int> selk(300);
  vector<double> selpos(300);
  vector<float> selcharge(300);
  vector<int> seldial(300);
  vector<int> selsign(300);
  vector<double> seleta(300);
  vector<double> seldr(300);
  vector<double> seldr0(300);
  vector<double> selmudial(300);
  vector<int> sellayer(300);
  vector<int> selsector(300);
  vector<int> selmmfe(300);
  vector<int> selupdown(300);
  vector<double> slope(300);
  vector<int> number(300);
  int layerselect[8]={0};
  double selectX[8]={0};
  vector<double> xofe(4);
  vector<double> yofe(4);
  vector<double> xyofe(4);
  double plateeta0 = 0;
  double platephi0 = 0;
  double platetheta0 = 0;
  double platetheta = 0;
  double platephi = 0;
  double plater = 0;
  double strip_phifraction = 0;
  double estimated_phifraction = 0;
  double minestimated_phi[2] = {0};
  double piphi = 0;
  double platemin = 0;
  double platedrsum = 0;
  double gtheta = 0;
  double ltheta = 0;
  double dtheta = 0;
  int l = 0;
  double theta3d = 0;
  double phi3d = 0;
  double mintheta3d[2] = {0};
  double minphi3d[2] = {0};
  int plateone = 0;
 
  int hitdial[2] = {0};//1-16, 0 means no hit
  int etacandidate=0;
  int stereocandidate=0;
  int etalayer=0;
  int stereolayer=0;
  int status = 0;
  int drstatus = 0;
  int finalstatus = 0;
  int total = 0;
  int totalreject[10] = {0};
  int track[2] = {0};
  int tracks = 0;
  int faketracks = 0;
  int truthtracks = 0;
  int sumoftrack = 0;
  int denomi = 0;
  int fakeornot = 0;
  int detection = 0;
  int acceptdetection = 0;
  int accepttracks = 0;
  int acceptfaketracks = 0;
  int accepttruthtracks = 0;
  int acceptcombis = 0;
  int acceptfakecombis = 0;
  int accepttruthcombis = 0;
  int mudialdetection = 0;
  int mudialtracks = 0;
  int mudialfaketracks = 0;
  int mudialtruthtracks = 0;
  int mudialcombis = 0;
  int mudialfakecombis = 0;
  int mudialtruthcombis = 0;
  double plateslope=0;
  double plateedge=0;
  int platestart=0;
  int etastopkk = 0;
  int stereostop = 0;
  int etadialcounter[16] = {0};
  int stereodialcounter[16] = {0};
  int n4X4UV = 0;
  int neffcc = 0;
  int ntruth = 0;
  int nsignal = 0;
  int nusing = 0;
  int truenX = 0;
  int truenUV = 0;
  int naXbUV = 0;
  int nvalid = 0;
  int nshortage = 0;
  int nboth = 0;
  int noverlap = 0;
  int nBC[10] = {0};
  int naccept[2] = {0};
  int nmudial4MMFE8[2][4] = {{0}};
  int mudial4MMFE8State[2][4][8] = {{{0}}};
  int mudial4MMFE8Number[2][4][8] = {{{0}}};
  int nmudial8HitsOver[2][4] = {{0}};
  int acceptsumof4MMFE8Number[16][4]={{0}};
  int acceptsumofSectorTrack[16]={0};
  int ntrue = 0;
  int nfake = 0;
  int ndead = 0;
  int nreject = 0;
  int ndetect = 0;
  int ncl = 0;
  int ntruedetect = 0;
  int combi = 0;
  int combicheck = 0;
  int combidenomi = 0;
  int combixdenomi = 0;
  int combiuvdenomi = 0;
  int sumofcombi = 0;
  int sumofxcombi = 0;
  int sumofuvcombi = 0;
  int sumoftruth = 0;
  int sumoffake = 0;
  int acceptsumoftrack = 0;
  int acceptsumoffake = 0;
  int acceptsumoftruth = 0;
  int acceptsumofcombi = 0;
  int acceptsumoffakecombi = 0;
  int acceptsumoftruthcombi = 0;
  int mudialsumoftrack = 0;
  int mudialsumoffake = 0;
  int mudialsumoftruth = 0;
  int mudialsumofcombi = 0;
  int mudialsumoffakecombi = 0;
  int mudialsumoftruthcombi = 0;
  int higher[33] = {0};
  int lower[33] = {0};
  int chosen[33] = {0};
  int platetracknum = {0};
  int selection[2][8] = {{0}};
  vector<int> selected(8);
  vector<int> smselectedmemory(8);
  vector<int> laselectedmemory(8);
  int allhit = 0;
  int truth = 0;
  int truedetection = 0;
  int cl = 0;
  int aXbUV = 0;
  int deadtime = 1;
  double chi2[2]={0};
  double amin=0;
  double legitimacy[2]={0};
  double drsum=0;
  double edm=0;
  double errdef=0;
  int nvpar=0;
  int nparx=0;
  int maxnX[2]={0};
  int platesector=0;
  int duplicity = 0;
  double minetadr = 0;
  int kEnd = 0;
  int mframe = 0;
  double frameedge = 0;
  int dual = 0;
  int BCID = 0;
  int BCIDrecog = 0;
  int N_Succ = 0;
  xyztag->Branch("ndata",&platendata);
  xyztag->Branch("mueta",&platemueta);
  xyztag->Branch("muentryeta",&platemuentryeta);
  xyztag->Branch("x",&platex);
  xyztag->Branch("y",&platey);
  xyztag->Branch("z",&platez);
  xyztag->Branch("pror",&platepror);
  xyztag->Branch("strip_phifrac",&strip_phifraction);
  xyztag->Branch("mux",&platemux);
  xyztag->Branch("muy",&platemuy);
  xyztag->Branch("muz",&platemuz);
  xyztag->Branch("mut",&platemut);
  xyztag->Branch("mux0",&platemux0);
  xyztag->Branch("muy0",&platemuy0);
  xyztag->Branch("muz0",&platemuz0);
  xyztag->Branch("mut0",&platemut0);
  xyztag->Branch("lx",&platelx);
  xyztag->Branch("ly",&plately);
  xyztag->Branch("i",&platei);
  xyztag->Branch("j",&platej);
  xyztag->Branch("k",&platek);
  xyztag->Branch("t",&platet);
  xyztag->Branch("pos",&platepos);
  xyztag->Branch("charge",&platecharge);
  xyztag->Branch("layer",&platelayer);
  xyztag->Branch("dial",&platedial);
  xyztag->Branch("sign",&platesign);
  xyztag->Branch("eta",&plateeta);
  xyztag->Branch("dr",&platedr);
  xyztag->Branch("dr0",&platedr0);
  xyztag->Branch("mudial",&platemudial);
  xyztag->Branch("slope",&plateslope);
  xyztag->Branch("edge",&plateedge);
  xyztag->Branch("start",&platestart);
  xyztag->Branch("tracknum",&platetracknum);
  xyztag->Branch("combi",&combi);
  xyztag->Branch("xcombi",&xcombi);
  xyztag->Branch("uvcombi",&uvcombi);
  xyztag->Branch("sector",&platesector);
  xyztag->Branch("band",&bands);
  xyztag->Branch("mframe",&mframe);
  xyztag->Branch("BCID",&BCID);


  int m = 0;
  double selplatedr0[33]={0};
  double etadrsum=0;
  int plateetadr=0;

  sel->Branch("ndata",&platendata);
  sel->Branch("i",&m);
  sel->Branch("ave_time",&platemin);
  sel->Branch("all_dr2",&platedrsum);
  sel->Branch("eta_dr2",&plateetadr);
  sel->Branch("mueta",&plateeta0);
  sel->Branch("mutheta",&platetheta0);
  sel->Branch("muphi",&platephi0);
  sel->Branch("bandeta",&plateeta);
  sel->Branch("bandtheta",&platetheta);
  sel->Branch("dialphi",&platephi);
  sel->Branch("phifrac",&estimated_phifraction);
  sel->Branch("theta3d",&theta3d);
  sel->Branch("phi3d",&phi3d);
  sel->Branch("nX",&etalayer);
  sel->Branch("nUV",&stereolayer);
  sel->Branch("one",&plateone);
  sel->Branch("fakeornot",&fakeornot);
  sel->Branch("band",&bands);
  sel->Branch("dr0",selplatedr0,"dr0[33]/D");
  sel->Branch("sector",&platesector);
  sel->Branch("mframe",&mframe);
  sel->Branch("BCID",&BCID);
  
  //int m = 0;
  int f = 0;
  int ok = 0;
  int viewon = 1;

  int layersign[33]={0};
  for(int rr=0;rr<1000;rr++){
    tree->GetEntry(rr);
    layersign[layer+16] = sign;
  }

  for(int cc=0;cc<etavec.size();cc++){//loop of eta
    etastart = etavec.at(cc)-etapitch/2.0;
    if(etastart<1.35) etastart = 1.35;
    if(etastart>1.5&&etastart<1.7) etastart = 1.7;
    etastop = etavec.at(cc)+etapitch/2.0;
    if(etastop>1.5&&etastop<1.7) etastop = 1.5;
    if(etastop>2.65) etastop = 2.65;
    etavec.at(cc) = (etastart+etastop)/2.0;
    cout<<etastart<<" "<<etastop<<endl;
    cout<<etavec.at(cc)<<endl;

  for(int tt=0;tt<timevec.size();tt++){//loop of window
 
    if(tt==0){
      for(int secint=0;secint<2;secint++){
	for(int d=0;d<4;d++){
	  mudial4MMFE8hist[secint][d] = new TH1D(Form("Sector%dEta%d_hithist",secint,d),Form("Sector%dEta%d_hithist",secint,d),30,-0.5,29.5);
	  mudial4MMFE8hist[secint][d]->GetXaxis()->SetTitle("hits/4MMFE8");
	}
      }
    }

    window = timevec.at(tt);
    cout<<Form("%.0lfns",window)<<endl;

  for(int aa=0;aa<tractnum;aa++){//loop of xuv

    nX = tract[aa][0];
    nUV = tract[aa][1];
    cout<<nX<<"X"<<nUV<<"UV"<<endl;
    if(tt==0){
      bandhist[aa] = new TH1D(Form("%dX%dUV_bandhist",tract[aa][0],tract[aa][1]),Form("%dX%dUV_bandhist",tract[aa][0],tract[aa][1]),30,-0.5,29.5);
      bandhist[aa]->GetXaxis()->SetTitle("faketracks");
      drdist[aa] = new TH1D(Form("%dX%dUV_etadr0dist",tract[aa][0],tract[aa][1]),Form("%dX%dUV_etadr0dist",tract[aa][0],tract[aa][1]),100,0.,100.);
      drdist[aa]->GetXaxis()->SetTitle("etadr0[mm]");
      drdistresol[aa] = new TH1D(Form("%dX%dUV_etadr0dist",tract[aa][0],tract[aa][1]),Form("%dX%dUV_etadr0dist",tract[aa][0],tract[aa][1]),100,0.,100.);
      drdistresol[aa]->GetXaxis()->SetTitle("etadr0[mm]");
      //cout<<bandhist[aa]->GetEntries()<<endl;
    }

  for(int bb=0;bb<loopvec.size();bb++){//loop of co-width start
    slopewidth = loopvec.at(bb);
    cout<<"slopewidth"<<slopewidth<<endl;

    //initialize    
    m=0;
    f=0;
    n=0;
    
    total = 0;
    for(int d=0;d<10;d++){
      totalreject[d] = 0;
      nBC[d] = 0;
    }
    naccept[0] = 0;
    naccept[1] = 0;
    for(int d=0;d<4;d++){
      nmudial4MMFE8[0][d] = 0;
      nmudial4MMFE8[1][d] = 0;
      nmudial8HitsOver[0][d] = 0;
      nmudial8HitsOver[1][d] = 0;
      for(int pp=0;pp<16;pp++){
	acceptsumof4MMFE8Number[pp][d] = 0;
	if(d==0)acceptsumofSectorTrack[pp] = 0;
      }
    }//means see only 1sector
    sumoftrack = 0;
    denomi = 0;
    n4X4UV = 0;
    neffcc = 0;
    ntruth = 0;
    nsignal = 0;
    nusing = 0;
    naXbUV = 0;
    nvalid = 0;
    nshortage = 0;
    nboth = 0;
    noverlap = 0;
    ntrue = 0;
    nfake = 0;
    ndead = 0;
    nreject = 0;
    ndetect = 0;
    ncl = 0;
    ntruedetect = 0;
    combidenomi = 0;
    combixdenomi = 0;
    combiuvdenomi = 0;
    sumofcombi = 0;
    sumofxcombi = 0;
    sumofuvcombi = 0;
    sumoftruth = 0;
    sumoffake = 0;
    acceptsumoftrack = 0;
    acceptsumoffake = 0;
    acceptsumoftruth = 0;
    acceptsumofcombi = 0;
    acceptsumoffakecombi = 0;
    acceptsumoftruthcombi = 0;
    mudialsumoftrack = 0;
    mudialsumoffake = 0;
    mudialsumoftruth = 0;
    mudialsumofcombi = 0;
    mudialsumoffakecombi = 0;
    mudialsumoftruthcombi = 0;
    nbands = 0;    

    tree->GetEntry(0);
    while(n < nevent){
      
      ok = 0;

      //initiallize
      drstatus = 1;//mean rejected               
      truth = 0;
      finalstatus = 0;//"ok"
      allhit = 0;//mean muon allhit event   
      for(int d = 0; d < 33; d++){
	ffirst[d] = 0;
      }
      for(int d = 0; d < 33; d++){
	layeron[d] = 0;
	hiton[d] = 0;
      }
      for(int d=0;d<4;d++){
	for(int pp=0;pp<8;pp++){
	  mudial4MMFE8State[0][d][pp] = 0;
	  mudial4MMFE8State[1][d][pp] = 0;
	  mudial4MMFE8Number[0][d][pp] = 0;
	  mudial4MMFE8Number[1][d][pp] = 0;
	}
      }
      hitdial[0] = 0;
      hitdial[1] = 0;
      sectoron[0] = 0;//small
      sectoron[1] = 0;//large
      hitsectoron[0] = 0;
      hitsectoron[1] = 0;
      allhitsectoron[0] = 0;
      allhitsectoron[1] = 0;
      selmueta.erase(selmueta.begin(),selmueta.end());
      selmuentryeta.erase(selmuentryeta.begin(),selmuentryeta.end());
      seltime.erase(seltime.begin(),seltime.end());
      sellayer.erase(sellayer.begin(),sellayer.end());
      selx.erase(selx.begin(),selx.end());    
      sely.erase(sely.begin(),sely.end());
      selz.erase(selz.begin(),selz.end());
      selpror.erase(selpror.begin(),selpror.end());    
      selslope.erase(selslope.begin(),selslope.end());
      selmux.erase(selmux.begin(),selmux.end());
      selmuy.erase(selmuy.begin(),selmuy.end());
      selmuz.erase(selmuz.begin(),selmuz.end());
      selmut.erase(selmut.begin(),selmut.end());
      selmux0.erase(selmux0.begin(),selmux0.end());
      selmuy0.erase(selmuy0.begin(),selmuy0.end());
      selmuz0.erase(selmuz0.begin(),selmuz0.end());
      selmut0.erase(selmut0.begin(),selmut0.end());
      sellx.erase(sellx.begin(),sellx.end());
      selly.erase(selly.begin(),selly.end());
      seli.erase(seli.begin(),seli.end());
      selj.erase(selj.begin(),selj.end());
      selk.erase(selk.begin(),selk.end());
      selpos.erase(selpos.begin(),selpos.end());
      selcharge.erase(selcharge.begin(),selcharge.end());
      seldial.erase(seldial.begin(),seldial.end());
      selsign.erase(selsign.begin(),selsign.end());
      seleta.erase(seleta.begin(),seleta.end());
      seldr.erase(seldr.begin(),seldr.end());
      seldr0.erase(seldr0.begin(),seldr0.end());
      selmudial.erase(selmudial.begin(),selmudial.end());
      selsector.erase(selsector.begin(),selsector.end());
      selmmfe.erase(selmmfe.begin(),selmmfe.end());
      selupdown.erase(selupdown.begin(),selupdown.end());
      
      slope.erase(slope.begin(),slope.end());
      number.erase(number.begin(),number.end());
      
      while(m == i){//start of mi     

	if((n+1)%1000000 == 0) cout<<(n+1)/1.0/nevent<<endl;
	if(ok==0){
	  memo->GetEntry(m);
	  platendata = ndata;
	  
	  if(truthcheck==1){      	
	    for(int d=0;d<33;d++){
	      //if(d-16>0)cout<<d-16<<" "<<hiteta0[d]<<endl;
	      if((hiteta0[d]<1.5||1.7<hiteta0[d])&&(etastart<hiteta0[d]&&hiteta0[d]<etastop)){
		//cout<<hiteta0[d]<<endl;
		//cout<<d-16<<" "<<fabs(hitphi0[d]-(round((hitphi0[d]+((32-d)/8)*pi/8.0)*4.0/pi)*pi/4.0-((32-d)/8)*pi/8.0))<<endl;
		if(fabs(hitphi0[d]-(round((hitphi0[d]+((32-d)/8)*pi/8.0)*4.0/pi)*pi/4.0-((32-d)/8)*pi/8.0))<phiwidth[(d-17)/8]){
		hiton[d] = 1;
		//cout<<hitphi<<" : "<<layer<<" : "<<fabs(hitphi-(round((hitphi+((16-layer)/8)*pi/8.0)*4.0/pi)*pi/4.0-((16-layer)/8)*pi/8.0))<<endl;        
		}
	      }
	    }
	    //if((etastart<mueta&&mueta<1.5)||(1.7<mueta&&mueta<etastop)) ok = 1;
	    if(hiton[17]+hiton[18]+hiton[23]+hiton[24]>=nX){
	      if(hiton[19]+hiton[20]+hiton[21]+hiton[22]>=nUV){
		hitsectoron[0] = 1;
	      }
	    }
	    if(hiton[17]+hiton[18]+hiton[23]+hiton[24]+hiton[19]+hiton[20]+hiton[21]+hiton[22]==8){
	      truth = 1;
	      allhit = 1;
	      allhitsectoron[0] = 1;
	      naccept[0] += 1;
	    }
	    if(allhitsectoron[0]==1){
	      if(hitphi0[20]<0){
		hitdial[0] = round((hitphi0[20]+2*pi)/pi*4.0+0.5)*2;
	      }else{
		hitdial[0] = round(hitphi0[20]/pi*4.0+0.5)*2;
	      }
	      if(hitdial[0]==17) hitdial[0]=1;
	    }
	    if(hiton[25]+hiton[26]+hiton[31]+hiton[32]>=nX){
	      if(hiton[27]+hiton[28]+hiton[29]+hiton[30]>=nUV){
		hitsectoron[1] = 1;
	      }
	    }
	    if(hiton[25]+hiton[26]+hiton[31]+hiton[32]+hiton[27]+hiton[28]+hiton[29]+hiton[30]==8){
	      truth = 1;
	      allhit = 1;
	      allhitsectoron[1] = 1;
	      naccept[1] += 1;
	    }
	    if(allhitsectoron[1]==1){
	      if(hitphi0[28]<0){
		hitdial[1] = round((hitphi0[28]+2*pi)/pi*4.0)*2+1;
	      }else{
		hitdial[1] = round(hitphi0[28]/pi*4.0)*2+1;
	      }
	      if(hitdial[1]==17) hitdial[1]=1;
	    }

	  }
	}
	ok = 1;	

	checkmmfe = -1;//mean no selection for mmfe
	//cout<<checkmmfe<<endl;

	if(abs(layer)>8){sector=1;}else{sector=0;}

	viewon = 1;//mean view  
	if(viewSL[sector]==0)viewon = 0;
	if(viewACC==1&&allhitsectoron[sector]==0)viewon = 0;
	if(viewACC==2&&truth==0)viewon = 0;
	if(viewMuDial==1){
	  if(hitdial[sector]==0){
	    viewon = 0;
	  }else{
	    if(dial!=hitdial[sector]){
	      viewon = 0;
	    }
	  }
	}
	if(layer!=0&&viewon==1){
	  if(viewMuDial==1&&hitdial[sector]==0) alert[1]+=1;
	  if(allhitsectoron[sector]==1&&dial==hitdial[sector]&&BCstart*25<=t&&t<25*(BCstart+1)){
	    nmudial4MMFE8[sector][(mmfe+80*(updown-1))/32] += 1;
	    //accept && mudial && in 25ns
	    mudial4MMFE8Number[sector][(mmfe+80*(updown-1))/32][layer-sector*8-1] += 1;
	  }
	  if(allhitsectoron[sector]==1&&BCstart*25<=t&&t<25*(BCstart+1)){
	    //cout<<Dialturn(dial,hitdial[sector])<<endl;
	    acceptsumof4MMFE8Number[Dialturn(dial,hitdial[sector])-1][(mmfe+80*(updown-1))/32] += 1;
	  }
	  if(mmfeon!=0){
	    mmfecount = 0;
	    if(sellayer.end()!=find(sellayer.begin(),sellayer.end(),layer)){
	      for(int qq=0;qq<sellayer.size();qq++){
		if(sellayer.at(qq)==layer&&floor(seltime.at(qq)/25)==floor(t/25)&&(selmmfe.at(qq)+80*(selupdown.at(qq)-1))/32==(mmfe+80*(updown-1))/32&&selupdown.at(qq)==updown&&seldial.at(qq)==dial){
		  mmfecount += 1;
		}
	      }
	    }
	    if(mmfecount==8){
	      if(express==1){cout<<(mmfe+80*(updown-1))/32+1<<"thFiber"<<floor(t/25)<<"BC"<<endl;}
	      if(allhitsectoron[sector]==1&&dial==hitdial[sector]&&BCstart*25<=t&&t<25*(BCstart+1)){
		if(mudial4MMFE8State[sector][(mmfe+80*(updown-1))/32][layer-sector*8-1]==0){
		  mudial4MMFE8State[sector][(mmfe+80*(updown-1))/32][layer-sector*8-1] = 1;
		  nmudial8HitsOver[sector][(mmfe+80*(updown-1))/32] += 1;
		}
	      }
	      checkmmfe = -2;//mean reject
	      if(mmfeon==1||mmfeon==-1)mmfeRank = pror;
	      if(mmfeon==-2)mmfeRank = t;
	      for(int qq=0;qq<sellayer.size();qq++){
       if(sellayer.at(qq)==layer&&floor(seltime.at(qq)/25)==floor(t/25)&&(selmmfe.at(qq)+80*(selupdown.at(qq)-1))/32==(mmfe+80*(updown-1))/32&&selupdown.at(qq)==updown){
		  if(mmfeon==1){
		    if(selpror.at(qq)<mmfeRank){
		      checkmmfe = qq;
		      mmfeRank = selpror.at(qq);
		    }
		  }
		  if(mmfeon==-1){
  if(selpror.at(qq)>mmfeRank){
    checkmmfe = qq;
    mmfeRank = selpror.at(qq);
  }
         }
		  if(mmfeon==-2){
           if(seltime.at(qq)>mmfeRank){
                      checkmmfe = qq;
                      mmfeRank = seltime.at(qq);
                    }
                  }
		}
              }
	      //if(seltime.at(qq)>t){ checkmmfe = qq; }else{checkmmfe = -2;/*mean reject*/}       
	    }
	  }else{checkmmfe = -1;}
	  if(checkmmfe==-1){
	    seltime.push_back(t);
	    selslope.push_back(pror/z);
	    slope.push_back(pror/z);
	    selx.push_back(x);
	    sely.push_back(y);
	    selz.push_back(z);
	    selpror.push_back(pror);
	    selmux.push_back(mux);
	    selmuy.push_back(muy);
	    selmuz.push_back(muz);
	    selmut.push_back(mut);
	    selmux0.push_back(mux0);
	    selmuy0.push_back(muy0);
	    selmuz0.push_back(muz0);
	    selmut0.push_back(mut0);
	    sellx.push_back(lx);
	    selly.push_back(ly);
	    seli.push_back(i);
	    selj.push_back(j);
	    selk.push_back(k);
	    selpos.push_back(pos);
	    selcharge.push_back(charge);
	    seldial.push_back(dial);
	    selsign.push_back(sign);
	    seleta.push_back(eta);
	    selmueta.push_back(mueta);
	    selmuentryeta.push_back(muentryeta);
	    seldr.push_back(dr);
	    seldr0.push_back(dr0);
	    selmudial.push_back(mudial);
	    sellayer.push_back(layer);
	    selmmfe.push_back(mmfe);
	    selupdown.push_back(updown);
	    if(abs(layer)>8){selsector.push_back(1);}else{selsector.push_back(0);}
	  }else{
	    //cout<<checkmmfe<<endl;
	    if(checkmmfe!=-2){
	      seltime.at(checkmmfe) = t;
	      selslope.at(checkmmfe) = pror/z;
	      slope.at(checkmmfe) = pror/z;
	      selx.at(checkmmfe) = x;
	      sely.at(checkmmfe) = y;
	      selz.at(checkmmfe) = z;
	      selpror.at(checkmmfe) = pror;
	      selmux.at(checkmmfe) = mux;
	      selmuy.at(checkmmfe) = muy;
	      selmuz.at(checkmmfe) = muz;
	      selmut.at(checkmmfe) = mut;
	      selmux0.at(checkmmfe) = mux0;
	      selmuy0.at(checkmmfe) = muy0;
	      selmuz0.at(checkmmfe) = muz0;
	      selmut0.at(checkmmfe) = mut0;
	      sellx.at(checkmmfe) = lx;
	      selly.at(checkmmfe) = ly;
	      seli.at(checkmmfe) = i;
	      selj.at(checkmmfe) = j;
	      selk.at(checkmmfe) = k;
	      selpos.at(checkmmfe) = pos;
	      selcharge.at(checkmmfe) = charge;
	      seldial.at(checkmmfe) = dial;
	      selsign.at(checkmmfe) = sign;
	      seleta.at(checkmmfe) = eta;
	      selmueta.at(checkmmfe) = mueta;
	      selmuentryeta.at(checkmmfe) = muentryeta;
	      seldr.at(checkmmfe) = dr;
	      seldr0.at(checkmmfe) = dr0;
	      selmudial.at(checkmmfe) = mudial;
	      sellayer.at(checkmmfe) = layer;
	      selmmfe.at(checkmmfe) = mmfe;
	      selupdown.at(checkmmfe) = updown;
	      if(layer>8){selsector.at(checkmmfe) = 1;}else{selsector.at(checkmmfe) = 0;}
	    }
	  }
	}      
	
	//cout<<n<<endl;
	n = n+1;
	if(n < nevent){
	  tree->GetEntry(n);
	}
	if(n == nevent){
	  m = 0;
	}
      }//last of mi while
      
      //warning m!=i please using m
      
      if(ok==1&&vertex_n>nvertexcrit){
	if(truth==1&&tt==0&&aa==defaa&&bb==defbb){
	  for(int secint=0;secint<2;secint++){
	    for(int d=0;d<4;d++){
	      for(int pp=0;pp<8;pp++){
		if(allhitsectoron[secint]==1)mudial4MMFE8hist[secint][d]->Fill(mudial4MMFE8Number[secint][d][pp]);
		//cout<<mudial4MMFE8Number[secint][d]<<endl;     
	      }  
	    }
	  }
	}
	nsignal += 1;
	if(seltime.size()!=0)nusing += 1;
	drstatus = 1;//mean rejected
	test = 0;
	BCIDrecog = BCstart;
	N_Succ = 0;	

	for(int qq=0;qq<sellayer.size();qq++){
	  if(fabs(seldr0.at(qq))<drcut&&BCstart*25<=seltime.at(qq)&&seltime.at(qq)<BCstart*25+window) layeron[sellayer.at(qq)+16] = 1; 	
	}
	

	if(layeron[17]+layeron[18]+layeron[23]+layeron[24]>=nX){
	  if(layeron[19]+layeron[20]+layeron[21]+layeron[22]>=nUV){
	    //if(hitsectoron[0]==1) naXbUV += 1;
	    nvalid += 1;
	    test = 1;
	    drstatus = 0;
	    if(allhitsectoron[0]==1) neffcc += 1;
	  }
	}

	if(layeron[25]+layeron[26]+layeron[31]+layeron[32]>=nX){
	  if(layeron[27]+layeron[28]+layeron[29]+layeron[30]>=nUV){
	    if(test==0) nvalid += 1;
	    if(test==1) nboth += 1;
	    test = 1;
	    drstatus = 0;
	    if(allhitsectoron[1]==1) neffcc += 1;
	  }
	}
	
	if(sellayer.size()!=0){
	  if(test==0){
	    nshortage += 1;
	    drstatus = 1;
	  }
	  if(truth==1) ntruth += 1;
	  if(allhit==1) n4X4UV += 1;
	}
	//if(allhit==1&&test==1) neffcc += 1;      
	
	//sort
	sort(slope.begin(),slope.end());
	for(int k=0;k<slope.size();k++){
	  for(int a=0;a<selslope.size();a++){
	    if(slope.at(k)==selslope.at(a)){
	      if(number.end()==find(number.begin(),number.end(),a)){
		number.push_back(a);
	      }
	    }
	  }
	}

	detection = 0;             
        acceptdetection = 0;
	accepttracks = 0;
	acceptfaketracks = 0;
	accepttruthtracks = 0;
	acceptcombis = 0;
	acceptfakecombis = 0;
	accepttruthcombis = 0;
	mudialdetection = 0;
        mudialtracks = 0;
        mudialfaketracks = 0;
        mudialtruthtracks = 0;
        mudialcombis = 0;
        mudialfakecombis = 0;
        mudialtruthcombis = 0;

	//cout<<"coincidence"<<endl;
	for(int BCtime=BCstart;BCtime<=backcheck+BCstart;BCtime++){	
	  BCID = BCtime;
	  //cout<<"BCID"<<BCID<<endl;

	//coincidence
	//initialize
	for(int d=0;d<16;d++){
	  etadialcounter[d]=0;
	  stereodialcounter[d]=0;
	}
	
	track[0] = -1; track[1] = -1;//mean reject
	tracks = 0;
	band = -1;
	mframe = 0;
	bands = 0;
	chi2[0] = 0; chi2[1] = 0;
	legitimacy[0] = 0; legitimacy[1] = 0;
	maxnX[0] = 0; maxnX[1] = 0;
	for(int d=0;d<8;d++){
	  selection[0][d] = -1;//mean reject
          selection[1][d] = -1;
	}
	faketracks = 0;
	truthtracks = 0;
	cl = 0;
	aXbUV = 0;
	sectoron[0] = 0;//small
	sectoron[1] = 0;//large      
	for(int k = 0; k < slope.size(); k++){//start of loop of k
	  frameedge = 0;
	  kEnd = 0;
	  dual = 0;
	  while(kEnd==0){
	  status = 0;//track is "ok"
	  etacandidate = 0;
	  stereocandidate = 0;
	  etalayer = 0;
	  stereolayer = 0;
	  if(fixframe==0){
	    frameedge = slope.at(k);
	    kEnd=1;
	  }
	  if(fixframe==1){
            if(dual==0)mframe = floor(slope.at(k)*2.0/slopewidth-0.5)-1;
	    //cout<<mframe<<endl;
            frameedge = (mframe/2.0+0.25)*slopewidth;
	    dual += 1;
            if(dual==2)kEnd=1;
	  }

	  for(int kk=k;kk<slope.size();kk++){
	    if(BCtime*25<=seltime.at(number.at(k))&&seltime.at(number.at(k))<BCtime*25+window&&selsign.at(number.at(k))==0&&selsign.at(number.at(kk))==0&&BCtime*25<=seltime.at(number.at(kk))&&seltime.at(number.at(kk))<BCtime*25+window&&selsector.at(number.at(k))==selsector.at(number.at(kk))&&seldial.at(number.at(k))==seldial.at(number.at(kk))&&frameedge<=slope.at(kk)&&slope.at(kk)<frameedge+slopewidth){
	      etacandidate += 1;
	    }
	  }	  
	  
	  if(etacandidate >= nX){
	    //reset
	    for(int d = 0; d < 33; d++){
	      layeron[d] = 0;
	      higher[d] = 0;
	      lower[d] = 0;
	      chosen[d] = 0;
	    }
	    for(int d=0;d<8;d++){
	      layerselect[d] = -1;
	      selectX[d] = 0;
	    }
	    
	    for(int kk=k;kk<slope.size();kk++){
	      if(BCtime*25<=seltime.at(number.at(kk))&&seltime.at(number.at(kk))<BCtime*25+window&&selsign.at(number.at(kk))==0&&selsector.at(number.at(k))==selsector.at(number.at(kk))&&seldial.at(number.at(k))==seldial.at(number.at(kk))&&frameedge<=slope.at(kk)&&slope.at(kk)<frameedge+slopewidth){

		if(layerselect[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1]==-1){
		  layerselect[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1] = kk;
		  if(timeon==1)selectX[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1] = seltime.at(number.at(kk));
		  if(biason==1)selectX[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1] = (1-sellayer.at(number.at(kk))%2*2)*kk;
		  if(biason==-1)selectX[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1] = -(1-sellayer.at(number.at(kk))%2*2)*kk;
		}
		if(timeon==1){
		  if(seltime.at(number.at(kk)) < selectX[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1]){
		    layerselect[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1] = kk;
		    selectX[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1] = seltime.at(number.at(kk));
		  }
		}
		if(biason==1){
                  if((1-sellayer.at(number.at(kk))%2*2)*kk < selectX[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1]){
                    layerselect[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1] = kk;
                    selectX[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1] = (1-sellayer.at(number.at(kk))%2*2)*kk;
                  }
		}
		if(biason==-1){
                  if( -(1-sellayer.at(number.at(kk))%2*2)*kk < selectX[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1]){
                    layerselect[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1] = kk;
                    selectX[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1] = -(1-sellayer.at(number.at(kk))%2*2)*kk;
                  }
                }
		if(Ron==1){
		  layerselect[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1] = kk;
		}
		if(Ron==-1){
		  if(layerselect[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1]==-1){
		    layerselect[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1] = kk;
		  }
		}
		if(randomon==1){

		}
		if(layeron[sellayer.at(number.at(kk))+16]==0){
		  etalayer += 1;
		  layeron[sellayer.at(number.at(kk))+16] += 1;
		  //lower[sellayer.at(number.at(kk))+16] = kk;
		}else{
		  status = 1;//need to select error
		  layeron[sellayer.at(number.at(kk))+16] += 1;
		//higher[sellayer.at(number.at(kk))+16] = kk;
		} 
		etastopkk = kk;
	      }
	    }
	    
	    if(etalayer >= nX){
	      for(int kk=0;kk<slope.size();kk++){
		if(BCtime*25<=seltime.at(number.at(kk))&&seltime.at(number.at(kk))<BCtime*25+window&&selsign.at(number.at(kk))!=0&&selsector.at(number.at(k))==selsector.at(number.at(kk))&&seldial.at(number.at(k))==seldial.at(number.at(kk))&&slope.at(kk)<frameedge+slopewidth+stereowidth&&slope.at(kk)>=frameedge-stereowidth){
		  if(layerselect[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1]==-1){
		    layerselect[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1] = kk;
		    if(timeon==1)selectX[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1] = seltime.at(number.at(kk));
		    if(biason==1)selectX[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1] = (1-sellayer.at(number.at(kk))%2*2)*kk;
		    if(biason==-1)selectX[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1] = -(1-sellayer.at(number.at(kk))%2*2)*kk;
		  }
		  if(timeon==1){
		    if(seltime.at(number.at(kk)) < selectX[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1]){
		      layerselect[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1] = kk;
		      selectX[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1] = seltime.at(number.at(kk));
		    }
		  }
		  if(biason==1){
                    if((1-sellayer.at(number.at(kk))%2*2)*kk < selectX[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1]){
                      layerselect[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1] = kk;
                      selectX[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1] = (1-sellayer.at(number.at(kk))%2*2)*kk;
                    }
                  }
		  if(biason==-1){
                    if(-(1-sellayer.at(number.at(kk))%2*2)*kk < selectX[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1]){
                      layerselect[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1] = kk;
                      selectX[sellayer.at(number.at(kk))-8*selsector.at(number.at(k))-1] = -(1-sellayer.at(number.at(kk))%2*2)*kk;
                    }
                  }

		  if(layeron[sellayer.at(number.at(kk))+16]==0){
		    stereolayer += 1;
		    layeron[sellayer.at(number.at(kk))+16] += 1;
		  }else{
		    status = 1;//need to select error
		    layeron[sellayer.at(number.at(kk))+16] += 1;
		  }
		  stereostop = kk;
		}
	      }
	      
	      if(stereolayer >= nUV){
		//initialize
		banddetect = 0;
		combi = 1;
		xcombi = 1;
		uvcombi = 1;
		for(int d=0;d<33;d++){
		  selplatedr0[d] = -7000;
		}
		//cout<<m<<" dial "<<seldial.at(number.at(k))<<" "<<etastopkk<<" "<<etadialcounter[seldial.at(number.at(k))-1]<<endl;
		if(tracks>0){
		  if(etastopkk<=etadialcounter[seldial.at(number.at(k))-1]&&stereostop<=stereodialcounter[seldial.at(number.at(k))-1]){
		    //cout<<m<<" dial "<<seldial.at(number.at(k))<<" "<<etastopkk<<" "<<etadialcounter[seldial.at(number.at(k))-1]<<endl;
		    if(tracks == 1)finalstatus = 4;
		    status = 2;//one track include another
		  }else{
		    if(etastopkk>etadialcounter[seldial.at(number.at(k))-1]) etadialcounter[seldial.at(number.at(k))-1]=etastopkk;
		    if(stereostop>stereodialcounter[seldial.at(number.at(k))-1]) stereodialcounter[seldial.at(number.at(k))-1]=stereostop;
		    if(abs(sellayer.at(number.at(k)))<9){
		      sectoron[0] = 1;
		    }else{
		      sectoron[1] = 1;
		    }
		    if(sectoron[0]==1&&sectoron[1]==1&&status!=1&&tracks==1&&finalstatus!=5){
		      //tracks -= 1;
		      finalstatus = 5;
		    }
		    if(tracks>0) finalstatus = 3;
		    if(status==1){ //finalstatus = 2;
		      combi = 1;
		      xcombi = 1;
		      uvcombi = 1;
		      for(int d=1;d<17;d++){
			if(selsector.at(number.at(k))==(d-1)/8&&layeron[d+16]!=0){
			  combi *= layeron[d+16];
			  if(layersign[d+16]==0){xcombi *= layeron[d+16];}else{uvcombi *= layeron[d+16];}
			}
		      }
		    }
		  }
		}else{
		  etadialcounter[seldial.at(number.at(k))-1]=etastopkk;
		  stereodialcounter[seldial.at(number.at(k))-1]=stereostop;
		  if(abs(sellayer.at(number.at(k)))<9){
		    sectoron[0] = 1;
		  }else{
		    sectoron[1] = 1;
		  }
		  if(status==1){ //finalstatus = 2;
		    combi = 1;
		    xcombi = 1;
		    uvcombi= 1;
		    for(int d=1;d<17;d++){
		      if(selsector.at(number.at(k))==(d-1)/8&&layeron[d+16]!=0){
			combi *= layeron[d+16];
			if(layersign[d+16]==0){xcombi *= layeron[d+16];}else{uvcombi *= layeron[d+16];}
		      }
		    }
		  }
		}
		
		if(status == 1){
		  for(int d=1;d<17;d++){
		    if((d-1)/8==selsector.at(number.at(k))){
		      if(layersign[d+16]==0&&layeron[d+16]>2){
			deadtime = 0;
			//cout<<m<<" : "<<d<<" : "<<layeron[d+16]<<endl;
		      }
		      if(layeron[d+16]>2&&status!=4){
			status = 3;
			//cout<<d<<" : "<<layeron[d+16]<<endl;
			finalstatus = 6;
		      }
		      if(layersign[d+16]!=0){
			if(layeron[d+16]>1){
			  status = 4;//etacluster
			  finalstatus = 7;
			}
		      }
		    }
		  }
		}
		
		if(status==1)finalstatus = 2;
		if(status==3)finalstatus = 7;
		if(status==4)finalstatus = 8;		
		selected.erase(selected.begin(),selected.end());
		if(dthetacut!=0){		
		  xofe.erase(xofe.begin(),xofe.end());
		  yofe.erase(yofe.begin(),yofe.end());
		  xyofe.erase(xyofe.begin(),xyofe.end());
		}
		for(int d=0;d<8;d++){
		  if(layerselect[d]!=-1){
		    selected.push_back(layerselect[d]);
		    if(dthetacut!=0){
		      if(d==0||d==1||d==6||d==7){
			xofe.push_back(selz.at(number.at(layerselect[d])));
			yofe.push_back(selpror.at(number.at(layerselect[d])));
			xyofe.push_back(selz.at(number.at(layerselect[d]))*selpror.at(number.at(layerselect[d])));
		      }
		    }
		  }
		}

		if(duplicityerase != 0){
		  duplicity = 0;
		  if(tracks==0){
		    smselectedmemory.erase(smselectedmemory.begin(),smselectedmemory.end());
		    laselectedmemory.erase(laselectedmemory.begin(),laselectedmemory.end());
		  }

		  if(selsector.at(number.at(k))==0){
		    if(selected.size()==smselectedmemory.size()){
		      for(int dd=0;dd<selected.size();dd++){
			if(selected.at(dd)==smselectedmemory.at(dd)) duplicity += 1;
		      }
		      if(selected.size()==duplicity){
			duplicity = 1;
		      }else{
			duplicity = 0;
		      }
		    }
		  }

		  if(selsector.at(number.at(k))==1){
                    if(selected.size()==laselectedmemory.size()){
                      for(int dd=0;dd<selected.size();dd++){
                        if(selected.at(dd)==laselectedmemory.at(dd)) duplicity += 1;
                      }
                      if(selected.size()==duplicity){
                        duplicity = 1;
                      }else{
                        duplicity = 0;
                      }
                    }
                  }
		
		  //if(duplicity==1)cout<<m<<" : duplicity"<<duplicity<<endl;
		  if(selsector.at(number.at(k))==0){
		    smselectedmemory.erase(smselectedmemory.begin(),smselectedmemory.end());
		    for(int d=0;d<8;d++){
		      if(layerselect[d]!=-1)smselectedmemory.push_back(layerselect[d]);
		    }
		  }
                  if(selsector.at(number.at(k))==1){
                    laselectedmemory.erase(laselectedmemory.begin(),laselectedmemory.end());
                    for(int d=0;d<8;d++){
                      if(layerselect[d]!=-1)laselectedmemory.push_back(layerselect[d]);
                    }
                  }

		  if(duplicity == 1) status = 2;//mean ignore
		}

		if(duplicity==0&&dthetacut!=0){
		  if(xofe.size()!=etalayer)cout<<"xofeError"<<endl;
		  l = etalayer;
		  if(dthetacut==1){
		    ltheta = atan((TMath::Mean(l,&xyofe[0])-TMath::Mean(l,&xofe[0])*TMath::Mean(l,&yofe[0]))/(TMath::RMS(l,&xofe[0])*TMath::RMS(l,&xofe[0])*(l-1)/l));
		    gtheta = atan(TMath::Mean(l,&xyofe[0])/(TMath::RMS(l,&xofe[0])*TMath::RMS(l,&xofe[0])*(l-1)/l+TMath::Mean(l,&xofe[0])*TMath::Mean(l,&xofe[0])));
		    dtheta = ltheta-gtheta;
		    //cout<<"dtheta"<<dtheta<<endl;
		    if(fabs(dtheta)>0.015) status = 2;
		  }
		  if(dthetacut==2){
                    ltheta = atan((TMath::Mean(l,&xyofe[0])-TMath::Mean(l,&xofe[0])*TMath::Mean(l,&yofe[0]))/(TMath::RMS(l,&xofe[0])*TMath::RMS(l,&xofe[0])*(l-1)/l));
		    gtheta = 0;
		    for(int d=0;d<l;d++){
		      gtheta += yofe.at(d)/xofe.at(d);
		    }
                    gtheta = atan(gtheta/l);
                    dtheta = ltheta-gtheta;
                    //cout<<"dtheta"<<dtheta<<endl; 
		    if(fabs(dtheta)>0.015) status = 2;
                  }

		}

		if(status!=2){
		  if(selected.size()==etalayer+stereolayer){
		    truedetection = 1;
		    truenX = 0;
		    truenUV = 0;
		    fakeornot = 1;
		    amin = 0;
		    //initialized
		    estimated_phifraction = 0;
		    drsum = 0;
		    etadrsum = 0;				    
		    
		    platetheta0 = atan2(hypot(selmux.at(number.at(k)),selmuy.at(number.at(k))),selmuz.at(number.at(k)));
		    plateeta0 = -log(tan(platetheta0/2.0));
		    platephi0 = atan2(selmuy.at(number.at(k)),selmux.at(number.at(k)));
		    plateeta = seleta.at(number.at(k));
		    platetheta = 2.0*atan(exp(-plateeta));
		    platephi = (seldial.at(number.at(k))-1)*pi/8.0;
		    for(int cc=0;cc<selected.size();cc++){
		      amin += seltime.at(number.at(selected.at(cc)));
		      drsum += seldr0.at(number.at(selected.at(cc)))*seldr0.at(number.at(selected.at(cc)));
		      if(selsign.at(number.at(selected.at(cc)))==0) etadrsum += seldr0.at(number.at(selected.at(cc)))*seldr0.at(number.at(selected.at(cc)));
                      platet = seltime.at(number.at(selected.at(cc)));
		      platex = selx.at(number.at(selected.at(cc)));
		      platey = sely.at(number.at(selected.at(cc)));
		      platez = selz.at(number.at(selected.at(cc)));
		      platepror = selpror.at(number.at(selected.at(cc)));
		      platemux = selmux.at(number.at(selected.at(cc)));
                      platemuy = selmuy.at(number.at(selected.at(cc)));
                      platemuz = selmuz.at(number.at(selected.at(cc)));
                      platemut = selmut.at(number.at(selected.at(cc)));
                      platemux0 = selmux0.at(number.at(selected.at(cc)));
                      platemuy0 = selmuy0.at(number.at(selected.at(cc)));
                      platemuz0 = selmuz0.at(number.at(selected.at(cc)));
                      platemut0 = selmut0.at(number.at(selected.at(cc)));
		      platelx = sellx.at(number.at(selected.at(cc)));
                      plately = selly.at(number.at(selected.at(cc)));
                      platei = seli.at(number.at(selected.at(cc)));
                      platej = selj.at(number.at(selected.at(cc)));
                      platek = selk.at(number.at(selected.at(cc)));
                      platepos = selpos.at(number.at(selected.at(cc)));
                      platecharge = selcharge.at(number.at(selected.at(cc)));
		      platedial = seldial.at(number.at(selected.at(cc)));
		      if(platedial!=seldial.at(number.at(k))){
			cout<<"dialerror"<<platedial<<" : "<<seldial.at(number.at(k))<<endl;
			cout<<selsector.at(number.at(selected.at(cc)))<<" : "<<selsector.at(number.at(k))<<endl;
			alert[0] += 1;
		      }
		      platesector = selsector.at(number.at(selected.at(cc)));
		      platelayer = sellayer.at(number.at(selected.at(cc)));
		      plateeta = seleta.at(number.at(selected.at(cc)));
		      platemueta = selmueta.at(number.at(selected.at(cc)));
		      platemuentryeta = selmuentryeta.at(number.at(selected.at(cc)));
		      platedr = seldr.at(number.at(selected.at(cc)));
		      platesign = selsign.at(number.at(selected.at(cc)));
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
		      plater = hypot(platex,platey)*(cos(strip_phifraction)+platesign*sin(strip_phifraction)*tan(1.5*pi/180));
		      //platepror = plater;
		      estimated_phifraction += platesign*atan((plater-platez*tan(platetheta))/platez/tan(platetheta)/tan(1.5/180*pi));
		      platedr0 = seldr0.at(number.at(selected.at(cc)));
		      selplatedr0[platelayer-16] = platedr0;
		      if(fabs(platedr0)<drcut){
			//banddetect = 1;
			if(platesign==0){
			  banddetect = 1;
			  fakeornot = 0;
			  detection = 1;
			  if(allhitsectoron[platesector]==1){
			    truenX += 1;
			    acceptdetection = 1;
			    mudialdetection = 1;
			  }
			}
			if(platesign!=0&&allhitsectoron[platesector]==1)truenUV += 1;
		      }else{ if(allhitsectoron[platesector]==1){truedetection = 0;}}
		      platemudial = selmudial.at(number.at(selected.at(cc)));
		      plateslope = selslope.at(number.at(selected.at(cc)));
		      plateedge = frameedge;
		      platestart = sellayer.at(number.at(k));
		      platetracknum = tracks;
		      if(platet<BCtime*25 || platet>=BCtime*25+window || platelayer == 0){
			cout<<"error0"<<endl;
		      }
		      //cout<<"m"<<m<<"tracks"<<tracks<<"f"<<f<<endl;                                                                   
		      if(timeloop+widthloop+xuvloop+etadep+dontroot==0)xyztag->Fill();
		      f = f+1;

		    }
		    
		    estimated_phifraction /= stereolayer;
		    
		    platemin = amin/selected.size();
		    platedrsum = sqrt(drsum/selected.size());
		    plateetadr = sqrt(etadrsum/etalayer);
		    if(express==1){
		      cout<<"chi2/ndf"<<amin/selected.size()<<endl;
		    }

		    /*if( track[selsector.at(number.at(k))]==-1 || (amin/selected.size()<chi2[selsector.at(number.at(k))]&&track[selsector.at(number.at(k))]==k) || (((etalayer>maxnX[selsector.at(number.at(k))])||(etalayer==maxnX[selsector.at(number.at(k))]&&drsum/selected.size()<legitimacy[selsector.at(number.at(k))]))&&track[selsector.at(number.at(k))]!=k) ){
		      track[selsector.at(number.at(k))] = k;
		      maxnX[selsector.at(number.at(k))] = etalayer;
		      chi2[selsector.at(number.at(k))] = amin/selected.size();
		      legitimacy[selsector.at(number.at(k))] = drsum/selected.size();
		      mintheta3d[selsector.at(number.at(k))] = theta3d;
		      minphi3d[selsector.at(number.at(k))] = phi3d;
		      minestimated_phi[selsector.at(number.at(k))] = estimated_phifraction;
		      for(int d=0;d<8;d++){
			if(layerselect[d]!=-1){ selection[selsector.at(number.at(k))][d] = layerselect[d]; }else{selection[selsector.at(number.at(k))][d] = -1;}
		      }
		      }*/
		    plateeta = seleta.at(number.at(k));
		    platetheta = 2.0*atan(exp(-plateeta));
		    platephi = (seldial.at(number.at(k))-1)*pi/8.0;	    
		    plateone = 0;
		    if(timeloop+widthloop+xuvloop+etadep+dontroot==0)sel->Fill();
		    if((band!=k&&fixframe==0)||((band!=mframe||seldial.at(number.at(k)!=dialmemory))&&fixframe==1)){
                      if(fixframe==0)band = k;
		      if(fixframe==1)band = mframe;
		      bands += 1;
		      dialmemory = seldial.at(number.at(k));
                      nbands += 1;
                    }
		    tracks += 1;
		    if(tracks != bands)cout<<"errortracks dial";
		    //cout<<seldial.at(number.at(k))<<"tracks"<<tracks<<"bands"<<bands<<"band"<<band<<"mframe"<<mframe<<endl;
		    if(allhitsectoron[selsector.at(number.at(k))]==1){
		      if(accepttracks==0)minetadr = plateetadr;
		      if(minetadr>plateetadr)minetadr = plateetadr;
		      accepttracks += 1;
		      acceptsumofSectorTrack[Dialturn(seldial.at(number.at(k)),hitdial[selsector.at(number.at(k))])-1] += 1;
		      if(eventdebug==1){
			if(seldial.at(number.at(k))!=hitdial[selsector.at(number.at(k))]){
			  cout<<"n"<<m<<" dial"<<seldial.at(number.at(k))<<" ";
			  if(selsector.at(number.at(k))==1){cout<<"Large"<<endl;}else{cout<<"Small"<<endl;}
			  cout<<" pror"<<selpror.at(number.at(k))<<" X"<<selx.at(number.at(k))<<":"<<sely.at(number.at(k))<<":"<<selz.at(number.at(k))<<" time"<<seltime.at(number.at(k))<<endl;
			}
		      }
		      if(tt==deftt&&bb==defbb&&truth==1){
			if(combinationerase==0){
			  drdist[aa]->Fill(plateetadr,combi);
			}else{
			  drdist[aa]->Fill(plateetadr);
			}
		      }
		    }

		  }else{cout<<"error selected.size()!=etalayer+stereolayer"<<endl;}//selected.size()==etalayer+stereolayer
		  if(tracks>0&&truenX>=nX&&truenUV>=nUV&&truth==1) aXbUV=1; 
		  if(combi != xcombi*uvcombi) cout<<"error combi != xcombi*uvcombi"<<endl;
		  if(combi > 1) cl = 1;
		  sumofcombi += combi;
		  if(allhitsectoron[selsector.at(number.at(k))]==1){
		    acceptcombis += combi;
		    if(tracks>0&&banddetect==0){
		      acceptfakecombis += combi;
		    }
		    if(tracks>0&&banddetect==1){
                      accepttruthcombis += combi;
                    }
		  }
		  combidenomi += 1;
		  if(tracks>0&&banddetect==0){
		    faketracks += 1;
		    if(allhitsectoron[selsector.at(number.at(k))]==1){
		      acceptfaketracks += 1;
		    }
		  }
		  if(tracks>0&&banddetect==1){
		    truthtracks += 1;
		    if(allhitsectoron[selsector.at(number.at(k))]==1){
		      accepttruthtracks += 1;
		      if(N_Succ==0){
			N_Succ = 1; 
			BCIDrecog = BCID;
		      }
		      if(BCIDmethod==1){
			if(BCID==BCIDrecog+1){
			  N_Succ += 1;
			  if(N_Succ==2){
			    BCIDrecog = BCID;
			  }
			}
		      }
		      if(BCIDmethod==2){
                        if(BCID==BCIDrecog+1){
			  BCIDrecog = BCID;
			}
                      }

                    }
                  }
     
		  if(xcombi>1){
		    sumofxcombi += xcombi;
		    combixdenomi += 1;
		  }
		  if(uvcombi>1){
		    sumofuvcombi += uvcombi;
		    combiuvdenomi += 1;
		  }
		}
	      }
	    }
	  }//etacandidate>=nX
	  if(dual==1)mframe += 1;
	  }//last of kEnd==0 while
	  
	}//last of loop of k

	if(BCtime==BCstart){	
	  if(cl==1) ncl += 1;
	  if(tracks>0){
	    sumoftrack += tracks;
	    denomi += 1;
	    if(truth==1){
	      acceptsumoftrack += accepttracks;
	      acceptsumoffake += acceptfaketracks;
	      acceptsumoftruth += accepttruthtracks;
	      acceptsumofcombi += acceptcombis;
	      acceptsumoffakecombi += acceptfakecombis;
	      acceptsumoftruthcombi += accepttruthcombis;
	    }
	  }
	  if(sectoron[0]==1&&sectoron[1]==1){noverlap += 1;}
	  if(tracks>0&&acceptdetection==1&&truth==1) ndetect += 1;
	  //if(truth==1&&aa==0&&acceptdetection==0) cout<<"0nf"<<m<<endl;
	  if(tracks>0&&truedetection==1&&truth==1) ntruedetect += 1;
	  if(aXbUV==1) naXbUV += 1;
	  sumoffake += faketracks;
	  sumoftruth += truthtracks;
	  if(tracks == 0) finalstatus = 1;//no track
	  if(tracks > 0) total += 1;
	  if(finalstatus!=0) totalreject[finalstatus-1] += 1;
	  
	  //if(tracks>0&&drstatus==0) ntrue += 1;
	  //if(tracks>0&&drstatus==1) nfake += 1;
	  if(tracks>0&&detection==1) ntrue += 1;
	  if(tracks>0&&detection==0) nfake += 1;               
	  if(tracks==0&&drstatus==0) ndead += 1;
	  if(tracks==0&&drstatus==1) nreject += 1;
	  if(tt==deftt&&bb==defbb&&truth==1){
	    if(combinationerase==0){
	      bandhist[aa]->Fill(acceptfakecombis);
	    }else{
	      bandhist[aa]->Fill(acceptfaketracks);
	    }
	    if(accepttracks>0){
	      drdistresol[aa]->Fill(minetadr);
	    }
	  }
	}
        }//last of back check
	if(acceptdetection==1&&truth==1){
	  nBC[BCIDrecog-BCstart] += 1;
	}
      }//last of ok
      m++;    
    }

    if(etadep+widthloop+xuvloop+timeloop==0){
      cout<<"4MMFE8 "<<endl;
      cout<<"Small\t"<<"Large"<<endl;
      for(int d=3;d>-1;d--){
        cout<<mudial4MMFE8hist[0][d]->GetMean()*mudial4MMFE8hist[0][d]->GetEntries()/8.0/naccept[0]<<"\t"<<mudial4MMFE8hist[1][d]->GetMean()*mudial4MMFE8hist[1][d]->GetEntries()/8.0/naccept[1]<<endl;
      }
      cout<<endl;

      cout<<"8HitsOver "<<endl;
      cout<<"Small\t"<<"Large"<<endl;
      for(int d=3;d>-1;d--){
        cout<<mudial4MMFE8hist[0][d]->Integral(10,30)/8.0/naccept[0]*100.<<"%\t"<<mudial4MMFE8hist[1][d]->Integral(10,30)/8.0/naccept[1]*100.<<"%"<<endl;
      }
      cout<<endl;
    
      cout<<"8<="<<bandhist[aa]->Integral(9,30)*100.0/ntruth<<"%"<<endl;

      TCanvas *cMMFE[2][4];
      for(int secint=0;secint<2;secint++){
	for(int d=0;d<4;d++){
	  cMMFE[secint][d] = new TCanvas(Form("Sector%dEta%d",secint,d),Form("Sector%dEta%d",secint,d));
          mudial4MMFE8hist[secint][d]->Draw();
	  cMMFE[secint][d]->Print(Form("MMFEhits_Sector%dEta%d_%.2lf-%.2lf_%dBC%s.pdf",secint,d,etastart,etastop,BCstart,namecore));
	}
      }
      if(dontroot==0){
	filetag->cd();
	xyztag->Write();
	filesel->cd();
	sel->Write();
	filetag->Close();
	filesel->Close();
      }
    }
   
    cout<<BCstart<<"BCstart"<<endl;
    cout<<"eta "<<etastart<<"-"<<etastop<<endl;
    cout<<"slopewidth"<<slopewidth<<endl;
    cout<<nX<<"X"<<nUV<<"UV"<<endl;
    cout<<window<<"ns"<<endl<<endl;

    if(mmfeon==1){
      cout<<"4MMFE8 "<<endl;
      cout<<"Small\t"<<"Large"<<endl;
      for(int d=3;d>-1;d--){
	cout<<nmudial4MMFE8[0][d]/8.0/naccept[0]<<"\t"<<nmudial4MMFE8[1][d]/8.0/naccept[1]<<endl;
      }
      cout<<endl;
      
      cout<<"8HitsOver "<<endl;
      cout<<"Small\t"<<"Large"<<endl;
      for(int d=3;d>-1;d--){
        cout<<nmudial8HitsOver[0][d]/8.0/naccept[0]*100.<<"%\t"<<nmudial8HitsOver[1][d]/8.0/naccept[1]*100.<<"%"<<endl;
      }
      cout<<endl;

      cout<<"Small "<<naccept[0]<<" accepted"<<endl;
      cout<<"Large "<<naccept[1]<<" accepted"<<endl;

      cout<<"Dial Hit"<<endl;
      for(int pp=0;pp<16;pp++){
	cout<<"dial"<<pp+1<<" "<<acceptsumof4MMFE8Number[pp][0]/8.0/naccept[(pp+1)%2]<<"\t"<<endl;
      }
      cout<<endl;

      cout<<"Not Muon hit "<<endl;
      cout<<"Small\t"<<"Large"<<endl;
      for(int d=3;d>-1;d--){
        cout<<acceptsumof4MMFE8Number[9][d]/8.0/naccept[0]<<"\t"<<acceptsumof4MMFE8Number[8][d]/8.0/naccept[1]<<endl;
      }
      cout<<endl;

      cout<<"Dial Track"<<endl;
      for(int pp=0;pp<16;pp++){
        cout<<"dial"<<pp+1<<" "<<acceptsumofSectorTrack[pp]/1.0/naccept[(pp+1)%2]<<"\t"<<endl;
      }
      cout<<endl;

    }
    cout<<"nxyz "<<f<<endl<<endl; 
    cout<<"ntruth "<<ntruth<<endl;
    cout<<"naXbUV "<<naXbUV<<endl;
    cout<<"ndetect "<<ndetect<<endl;
    //cout<<"trueefficiency "<<ntruedetect/1.0/ntruth<<endl;
    cout<<"efficiency "<<ndetect/1.0/ntruth<<endl;
    cout<<"eff.cc "<<neffcc/1.0/ntruth<<endl;//for comparizon to eff.cc     
    cout<<"resolefficiency"<<naXbUV/1.0/ntruth<<endl<<endl;
    cout<<"ave-tracks "<<acceptsumoftrack/1.0/ntruth<<endl;
    cout<<"faketracks "<<acceptsumoffake/1.0/ntruth<<endl;
    if(etadep+widthloop+xuvloop+timeloop+dontroot!=0&&tt==deftt&&bb==defbb)cout<<"8<="<<bandhist[aa]->Integral(9,30)*100.0/ntruth<<"%"<<endl;
    cout<<"truthtrack "<<acceptsumoftruth/1.0/ntruth<<endl;
    cout<<"alltruthtrack "<<sumoftruth/1.0/ntruth<<endl;
    cout<<"ave-combis "<<acceptsumofcombi/1.0/ntruth<<endl;
    cout<<"fakecombis "<<acceptsumoffakecombi/1.0/ntruth<<endl;
    cout<<"truthcombis "<<acceptsumoftruthcombi/1.0/ntruth<<endl<<endl;
    if(backcheck!=0){
      for(int d=0;d<=backcheck;d++){
	cout<<BCstart+d<<"\t"<<nBC[d]/1.0/ntruth<<endl;
      }
    }
    cout<<"nsignal "<<nsignal<<endl;
    cout<<"nuse "<<nusing<<endl;
    cout<<"nvalid "<<nvalid<<endl;
    cout<<"accept "<<total<<endl;
    cout<<"ntrue "<<ntrue<<endl;
    cout<<"nfake "<<nfake<<endl;
    //cout<<"needs maximizing logic "<<totalreject[3]<<endl;
    //cout<<"nboth "<<nboth<<endl;
    //cout<<"both large and small "<<noverlap<<endl;
    //cout<<"nshortage "<<nshortage<<endl;
    cout<<"no candidate "<<totalreject[0]<<endl;
    cout<<"ndead "<<ndead<<endl;
    cout<<"ntrue "<<nreject<<endl;
    cout<<"needs selecting logic from VMMs "<<ncl<<endl;
    cout<<"multiple tracks (cluster) "<<totalreject[2]<<endl;
    //cout<<"fake track"<<totalreject[5]<<endl;
    //cout<<"fake track "<<sumoffake<<endl;    
    cout<<"3 candidates "<<totalreject[6]<<endl;
    cout<<"etacluster "<<totalreject[7]<<endl;
    //cout<<"average tracks "<<sumoftrack/1.0/total<<endl;
    //cout<<"average combis "<<sumofcombi/1.0/combidenomi<<endl;
    cout<<"etacombis "<<combixdenomi<<endl;
    cout<<"stereocombis "<<combiuvdenomi<<endl;
    if(deadtime==0)cout<<"no deadtime"<<endl;
    if(deadtime!=0)cout<<"deadtime on"<<endl;
    if(nbands!=combidenomi)cout<<"error"<<endl;    
    //cout<<"nbands"<<nbands<<endl;
    //cout<<"combidenomi"<<combidenomi<<endl;

    //Memory
    Ntruth[cc][bb][aa][tt] = ntruth;
    effcc[cc][bb][aa][tt] = neffcc/1.0/ntruth*100.;
    eff[cc][bb][aa][tt] = ndetect/1.0/ntruth*100.;
    resoleff[cc][bb][aa][tt] = naXbUV/1.0/ntruth*100.;
    if(effcc[cc][bb][aa][tt]!=0&&effcc[cc][bb][aa][tt]!=100){effccError[cc][bb][aa][tt] = sqrt(effcc[cc][bb][aa][tt]/ntruth*(100-effcc[cc][bb][aa][tt]));}else{effccError[cc][bb][aa][tt] = 0;}
    if(eff[cc][bb][aa][tt]!=0&&eff[cc][bb][aa][tt]!=100){effError[cc][bb][aa][tt] = sqrt(eff[cc][bb][aa][tt]/ntruth*(100-eff[cc][bb][aa][tt]));}else{effError[cc][bb][aa][tt] = 0;}
    if(resoleff[cc][bb][aa][tt]!=0&&resoleff[cc][bb][aa][tt]!=100){resoleffError[cc][bb][aa][tt] = sqrt(resoleff[cc][bb][aa][tt]/ntruth*(100-resoleff[cc][bb][aa][tt]));}else{resoleffError[cc][bb][aa][tt] = 0;}
   
    if(combinationerase==0){
      fakerate[cc][bb][aa][tt] = acceptsumoffakecombi*1.0/ntruth;
      trackrate[cc][bb][aa][tt] = acceptsumofcombi*1.0/ntruth;
      truthrate[cc][bb][aa][tt] = acceptsumoftruthcombi*1.0/ntruth;
    }else{
      fakerate[cc][bb][aa][tt] = acceptsumoffake*1.0/ntruth;
      trackrate[cc][bb][aa][tt] = acceptsumoftrack*1.0/ntruth;
      truthrate[cc][bb][aa][tt] = acceptsumoftruth*1.0/ntruth;
    }
    /*cout<<"nsignal"<<nsignal/1.0/nvalid<<endl;
      cout<<"nvalid"<<nvalid/1.0/nvalid<<endl;
      cout<<"nshortage"<<nshortage/1.0/nvalid<<endl;
      cout<<"nboth"<<nboth/1.0/nvalid<<endl;
      cout<<"nxyz "<<f/1.0/nvalid<<endl;
      cout<<"accept "<<total/1.0/nvalid<<endl;
      cout<<"ntrue"<<ntrue/1.0/nvalid<<endl;
      cout<<"nfake"<<nfake/1.0/nvalid<<endl;
      cout<<"needs maximizing logic "<<totalreject[3]/1.0/nvalid<<endl;
      cout<<"both large and small "<<totalreject[4]/1.0/nvalid<<endl;
      cout<<"no candidate "<<totalreject[0]/1.0/nvalid<<endl;
      cout<<"ndead"<<ndead/1.0/nvalid<<endl;
      cout<<"ntrue"<<nreject/1.0/nvalid<<endl;
      cout<<"needs selecting logic from VMM "<<totalreject[1]/1.0/nvalid<<endl;
      cout<<"multiple tracks (cluster) "<<totalreject[2]/1.0/nvalid<<endl;
      cout<<"fake track"<<totalreject[5]/1.0/nvalid<<endl;
      cout<<"3 candidates"<<totalreject[6]/1.0/nvalid<<endl;
      cout<<"etacluster"<<totalreject[7]/1.0/nvalid<<endl;
      cout<<"average tracks "<<sumoftrack/1.0/denomi<<endl;
      cout<<"average combis "<<sumofcombi/1.0/(totalreject[1]+totalreject[6]+totalreject[7])<<endl;*/

  }//loop of co-width last
  }//loop of xuv last
  }//loop of window last
  // }//loop of eta last

  if(xuvloop==1){
    if(widthloop==1){
      for(int aa=0;aa<tractnum;aa++){
	xuv[aa] = new TGraph(loopvec.size());
        xuvresol[aa] = new TGraph(loopvec.size());
	if(aa==0) xuvall = new TGraph(loopvec.size()*tractnum);
	if(aa==0) xuvallresol = new TGraph(loopvec.size()*tractnum);
	for(int bb=0;bb<loopvec.size();bb++){
	  xuv[aa]->SetPoint(bb,fakerate[cc][bb][aa][deftt],eff[cc][bb][aa][deftt]);
	  xuvresol[aa]->SetPoint(bb,fakerate[cc][bb][aa][deftt],resoleff[cc][bb][aa][deftt]);
	  xuvall->SetPoint(aa*loopvec.size()+bb,fakerate[cc][bb][aa][deftt],eff[cc][bb][aa][deftt]);
	  xuvallresol->SetPoint(aa*loopvec.size()+bb,fakerate[cc][bb][aa][deftt],resoleff[cc][bb][aa][deftt]);
	}
      }
      TCanvas *c1 = new TCanvas("eff-fake","eff-fake");
      xuvall->SetMarkerColor(0);
      xuvall->GetXaxis()->SetTitle("faketracks");
      xuvall->GetYaxis()->SetTitle("efficiency[%]");
      xuvall->Draw("AP");
      TLegend *xuvlegend = new TLegend(0.6,0.2,0.9,0.5);
      for(int aa=0;aa<tractnum;aa++){
	xuv[aa]->SetLineColor(aa+1);
	xuv[aa]->SetLineWidth(2);
	xuvlegend->AddEntry(xuv[aa],Form("%dX%dUV",tract[aa][0],tract[aa][1]),"lep");
	xuv[aa]->Draw("same");
      }
      xuvlegend->Draw();
      c1->Print(Form("eff-fake_eta%.2lf-%.2lf_%.0lfns%s.pdf",etastart,etastop,timevec.at(deftt),namecore));

      xuvallresol->SetMarkerColor(0);
      xuvallresol->GetXaxis()->SetTitle("faketracks");
      xuvallresol->GetYaxis()->SetTitle("efficiency[%]");
      xuvallresol->Draw("AP");
      for(int aa=0;aa<tractnum;aa++){
        xuvresol[aa]->SetLineColor(aa+1);
        xuvresol[aa]->SetLineWidth(2);
	xuvresol[aa]->Draw("same");
      }
      xuvlegend->Draw();
      c1->Print(Form("eff-fakeresol_eta%.2lf-%.2lf_%.0lfns%s.pdf",etastart,etastop,timevec.at(deftt),namecore));

    }
    
    
    TCanvas *c2 = new TCanvas("faketracks","faketracks");
    TLegend *drlegend = new TLegend(0.9,0.9,1.0,1.0);
    for(int aa=0;aa<tractnum;aa++){
      drdist[aa]->Scale(1.0/ntruth);
      drdist[aa]->GetYaxis()->SetRangeUser(0.00005,2);
      drdist[aa]->SetLineColor(aa+1);
      drdist[aa]->SetLineWidth(2);
      drdist[aa]->SetLineWidth(2);
      drlegend->AddEntry(drdist[aa],Form("%dX%dUV",tract[aa][0],tract[aa][1]),"lep");
      gPad->SetLogy();                                                          
      if(aa==0){drdist[aa]->Draw("");}else{drdist[aa]->Draw("same");}
    }
    drlegend->Draw();
    c2->Print(Form("drdist_eta%.2lf-%.2lf_slope%.5lf_%.0lfns%s.pdf",etastart,etastop,loopvec.at(defbb),timevec.at(deftt),namecore));

    for(int aa=0;aa<tractnum;aa++){
      drdistresol[aa]->Scale(1.0/ntruth);
      drdistresol[aa]->GetYaxis()->SetRangeUser(0.00005,2);
      drdistresol[aa]->SetLineColor(aa+1);
      drdistresol[aa]->SetLineWidth(2);
      drdistresol[aa]->SetLineWidth(2);
      gPad->SetLogy();              
      if(aa==0){drdistresol[aa]->Draw("");}else{drdistresol[aa]->Draw("same");}
    }
    drlegend->Draw();
    c2->Print(Form("drdistresol_eta%.2lf-%.2lf_slope%.5lf_%.0lfns%s.pdf",etastart,etastop,loopvec.at(defbb),timevec.at(deftt),namecore));

    TLegend *bandlegend = new TLegend(0.9,0.9,1.0,1.0);
    for(int aa=0;aa<tractnum;aa++){
      bandhist[aa]->SetLineColor(aa+1);
      bandhist[aa]->SetLineWidth(2);
      bandhist[aa]->SetLineWidth(2);
      bandlegend->AddEntry(bandhist[aa],Form("%dX%dUV",tract[aa][0],tract[aa][1]),"lep");
      gPad->SetLogy();
      if(aa==0){bandhist[aa]->Draw("");}else{bandhist[aa]->Draw("same");}
    }
    bandlegend->Draw();
    c2->Print(Form("faketrack_eta%.2lf-%.2lf_slope%.5lf_%.0lfns%s.pdf",etastart,etastop,loopvec.at(defbb),timevec.at(deftt),namecore));


    cout<<"Efficiency"<<endl;
    cout<<"\t";
    for(int j=0;j<tractnum;j++){
      cout<<Form("%dX%dUV \t",tract[j][0],tract[j][1]);
    }
    cout<<endl;
    cout<<"Efficiency";
    for(int j=0;j<tractnum;j++){
      cout<<"\t"<<eff[defcc][defbb][j][deftt];
    }
    cout<<endl;
    cout<<endl;

    cout<<"8<= band"<<endl;
    cout<<"\t";
    for(int j=0;j<tractnum;j++){
      cout<<Form("%dX%dUV \t",tract[j][0],tract[j][1]);
    }
    cout<<endl;
    cout<<"8<=";
    for(int j=0;j<tractnum;j++){
      cout<<"\t"<<bandhist[j]->Integral(9,30)*1.0/ntruth;
    }
    cout<<endl;
    cout<<endl;

    cout<<"Track Rate"<<endl;
    cout<<"\t";
    for(int j=0;j<tractnum;j++){
      cout<<Form("%dX%dUV \t",tract[j][0],tract[j][1]);
    }
    cout<<endl;
    cout<<"truth";
    for(int j=0;j<tractnum;j++){
      cout<<"\t"<<truthrate[cc][defbb][j][deftt];
    }
    cout<<endl;
    cout<<"fake";
    for(int j=0;j<tractnum;j++){
      cout<<"\t"<<fakerate[cc][defbb][j][deftt];
    }
    cout<<endl;
    cout<<"sum";
    for(int j=0;j<tractnum;j++){
      cout<<"\t"<<trackrate[cc][defbb][j][deftt];
    }
    cout<<endl;
    cout<<endl;

    cout<<"\t";
    for(int j=0;j<tractnum;j++){
      cout<<Form("%dX%dUV \t",tract[j][0],tract[j][1]);
    }
    cout<<endl;
    cout<<"<1mm";
    for(int j=0;j<tractnum;j++){
      cout<<"\t"<<drdist[j]->Integral(1,1);
    }
    cout<<endl;
    cout<<"1mm< <10mm";
    for(int j=0;j<tractnum;j++){
      cout<<"\t"<<drdist[j]->Integral(2,10);
    }
    cout<<endl;   
    cout<<"10mm< <100mm";
    for(int j=0;j<tractnum;j++){
      cout<<"\t"<<drdist[j]->Integral(11,100);
    }
    cout<<endl;
    cout<<"100mm<";
    for(int j=0;j<tractnum;j++){
      cout<<"\t"<<drdist[j]->Integral(101,101);
    }   
    cout<<endl;
    cout<<endl;

    cout<<"\t";
    for(int j=0;j<tractnum;j++){
      cout<<Form("%dX%dUV \t",tract[j][0],tract[j][1]);
    }
    cout<<endl;
    cout<<"<1mm";
    for(int j=0;j<tractnum;j++){
      cout<<"\t"<<drdistresol[j]->Integral(1,1);
    }
    cout<<endl;
    cout<<"1mm< <10mm";
    for(int j=0;j<tractnum;j++){
      cout<<"\t"<<drdistresol[j]->Integral(2,10);
    }
    cout<<endl;
    cout<<"10mm< <100mm";
    for(int j=0;j<tractnum;j++){
      cout<<"\t"<<drdistresol[j]->Integral(11,100);
    }
    cout<<endl;
    cout<<"100mm<";
    for(int j=0;j<tractnum;j++){
      cout<<"\t"<<drdistresol[j]->Integral(101,101);
    }
    cout<<endl;



    if(timeloop==1){
      for(int aa=0;aa<tractnum;aa++){
	nseff[aa] = new TGraphErrors(timevec.size(),&timevec[0],&eff[cc][defbb][aa][0],&zero[0],&effError[cc][defbb][aa][0]);
	nseff[aa]->GetYaxis()->SetRangeUser(0,110);
	
	if(aa==0)nseffall = new TGraph(tractnum*timevec.size());
	for(int tt=0;tt<timevec.size();tt++){
	  nseffall->SetPoint(aa*timevec.size()+tt,timevec.at(tt),eff[cc][defbb][aa][tt]);
	}
      }
      TCanvas *c3 = new TCanvas("nseff","nseff");
      nseffall->GetYaxis()->SetRangeUser(0,110);
      nseffall->SetMarkerColor(0);
      nseffall->GetXaxis()->SetTitle("timewindow[ns]");
      nseffall->GetYaxis()->SetTitle("efficiency[%]");
      nseffall->Draw("AP");
      TLegend *nsefflegend = new TLegend(0.6,0.2,0.9,0.5);
      for(int aa=0;aa<tractnum;aa++){
	nseff[aa]->SetLineColor(aa+1);
	nseff[aa]->SetLineWidth(2);
	nsefflegend->AddEntry(nseff[aa],Form("%dX%dUV",tract[aa][0],tract[aa][1]),"lep");
	nseff[aa]->Draw("same");
      }
      
      nsefflegend->Draw();
      c3->Print(Form("nseffetaallxuv_eta%.2lf-%.2lf_slope%.5lf_%s.pdf",etastart,etastop,loopvec.at(defbb),namecore));
      
      double p = 0;
      cout<<endl;
      cout<<"dLoclaX Efficiency"<<endl;
      cout<<"\t"<<"25ns \t"<<"50ns \t"<<"75ns \t"<<endl;
      for(int j=0;j<tractnum;j++){
	cout<<Form("%dX%dUV",tract[j][0],tract[j][1])<<"\t";
	for(int s=0;s<timevec.size();s++){
	  p = effcc[cc][defbb][j][s]/100.;
	  if(ntruth==0)p=0;
	  if(p==0){
	    cout<<fixed<<setprecision(0)<<p*100.0<<"(+-"<<sqrt(p*(1-p)/ntruth)*100.0<<")%\t";
	  }else{
	    cout<<fixed<<setprecision(-floor(log10(sqrt(p*(1-p)/ntruth)*100.0)))<<p*100.0<<"(+-"<<sqrt(p*(1-p)/ntruth)*100.0<<")%\t";
	  }
	}
	cout<<endl;
      }
      cout<<endl;
      cout<<"Slope Coincidence Efficiency"<<endl;
      cout<<"\t"<<"25ns \t"<<"50ns \t"<<"75ns \t"<<endl;
      for(int j=0;j<tractnum;j++){
	cout<<Form("%dX%dUV",tract[j][0],tract[j][1])<<"\t";
	for(int s=0;s<timevec.size();s++){
	  p = eff[cc][defbb][j][s]/100.;
	  if(ntruth==0)p=0;
	  if(p==0){
	    cout<<fixed<<setprecision(0)<<p*100.0<<"(+-"<<sqrt(p*(1-p)/ntruth)*100.0<<")%\t";
	  }else{
	    cout<<fixed<<setprecision(-floor(log10(sqrt(p*(1-p)/ntruth)*100.0)))<<p*100.0<<"(+-"<<sqrt(p*(1-p)/ntruth)*100.0<<")%\t";
	  }
	}
	cout<<endl;
      }
      cout<<endl;
      cout<<"For Resolution Efficiency"<<endl;
      cout<<"\t"<<"25ns \t"<<"50ns \t"<<"75ns \t"<<endl;
      for(int j=0;j<tractnum;j++){
	cout<<Form("%dX%dUV",tract[j][0],tract[j][1])<<"\t";
	for(int s=0;s<timevec.size();s++){
	  p = resoleff[cc][defbb][j][s]/100.;
	  if(ntruth==0)p=0;
	  if(p==0){
	    cout<<fixed<<setprecision(0)<<p*100.0<<"(+-"<<sqrt(p*(1-p)/ntruth)*100.0<<")%\t";
	  }else{
	    cout<<fixed<<setprecision(-floor(log10(sqrt(p*(1-p)/ntruth)*100.0)))<<p*100.0<<"(+-"<<sqrt(p*(1-p)/ntruth)*100.0<<")%\t";
	  }
	}
	cout<<endl;
      }
      
      
      
    }
  }
  }//loop of eta last
 
  if(etadep==1){
    for(int aa=0;aa<tractnum;aa++){
      nseffetadep[aa] = new TGraphErrors(etavec.size());
      nseffetadep[aa]->GetYaxis()->SetRangeUser(80,100);
      if(aa==0)nseffetadepall = new TGraph(tractnum*etavec.size());
      for(int cc=0;cc<etavec.size();cc++){
	nseffetadep[aa]->SetPoint(cc,etavec.at(cc),eff[cc][defbb][aa][deftt]);
	nseffetadep[aa]->SetPointError(cc,0.,effError[cc][defbb][aa][deftt]);
	nseffetadepall->SetPoint(aa*etavec.size()+cc,etavec.at(cc),eff[cc][defbb][aa][deftt]);
      }
    }
    TCanvas *c10 = new TCanvas("nseffetadep","nseffetadep");
    nseffetadepall->GetYaxis()->SetRangeUser(80,100);
    nseffetadepall->SetMarkerColor(0);
    nseffetadepall->GetXaxis()->SetTitle("eta");
    nseffetadepall->GetYaxis()->SetTitle("efficiency[%]");
    nseffetadepall->Draw("AP");
    TLegend *nseffetadeplegend = new TLegend(0.6,0.2,0.9,0.4);
    for(int aa=0;aa<tractnum;aa++){
      nseffetadep[aa]->SetLineColor(aa+1);
      nseffetadep[aa]->SetLineWidth(2);
      nseffetadeplegend->AddEntry(nseffetadep[aa],Form("%dX%dUV",tract[aa][0],tract[aa][1]),"lep");
      nseffetadep[aa]->Draw("same");
    }
    nseffetadeplegend->Draw();
    c10->Print(Form("nseffetadepxuv_eta%.2lf-%.2lf_slope%.5lf_%.0lf_%s.pdf",eta1,eta2,loopvec.at(defbb),timevec.at(deftt),namecore));

    for(int aa=0;aa<tractnum;aa++){
      nsresoleffetadep[aa] = new TGraphErrors(etavec.size());
      nsresoleffetadep[aa]->GetYaxis()->SetRangeUser(80,100);
      if(aa==0)nsresoleffetadepall = new TGraph(tractnum*etavec.size());
      for(int cc=0;cc<etavec.size();cc++){
	nsresoleffetadep[aa]->SetPoint(cc,etavec.at(cc),resoleff[cc][defbb][aa][deftt]);
	nsresoleffetadep[aa]->SetPointError(cc,0.,resoleffError[cc][defbb][aa][deftt]);
	nsresoleffetadepall->SetPoint(aa*etavec.size()+cc,etavec.at(cc),resoleff[cc][defbb][aa][deftt]);
      }
    }
    TCanvas *c11 = new TCanvas("nsresoleffetadep","nsresoleffetadep");
    nsresoleffetadepall->GetYaxis()->SetRangeUser(80,100);
    nsresoleffetadepall->SetMarkerColor(0);
    nsresoleffetadepall->GetXaxis()->SetTitle("timewindow[ns]");
    nsresoleffetadepall->GetYaxis()->SetTitle("efficiency[%]");
    nsresoleffetadepall->Draw("AP");
    TLegend *nsresoleffetadeplegend = new TLegend(0.6,0.2,0.9,0.4);
    for(int aa=0;aa<tractnum;aa++){
      nsresoleffetadep[aa]->SetLineColor(aa+1);
      nsresoleffetadep[aa]->SetLineWidth(2);
      nsresoleffetadeplegend->AddEntry(nsresoleffetadep[aa],Form("%dX%dUV",tract[aa][0],tract[aa][1]),"lep");
      nsresoleffetadep[aa]->Draw("same");
    }
    nsresoleffetadeplegend->Draw();
    c11->Print(Form("nsresoleffetadepxuv_eta%.2lf-%.2lf_slope%.5lf_%.0lf_%s.pdf",eta1,eta2,loopvec.at(defbb),timevec.at(deftt),namecore));

  }
  if(etadeplog==1){

    string filename;
    string asymmfilename;
    for(int aa=0;aa<tractnum;aa++){
      nX = tract[aa][0];
      nUV = tract[aa][1];
      cout<<nX<<"X"<<nUV<<"UV"<<endl;
      for(int tt=0;tt<timevec.size();tt++){
	window = timevec.at(tt);
	cout<<Form("%.0lfns",window)<<endl;
	filename = Form("%.0lfns%dX%dUVetadepMu%d.txt",window,nX,nUV,MU);
	asymmfilename = Form("%.0lfns%dX%dUVetadepMu%dAsymm.txt",window,nX,nUV,MU);
	std::ofstream outputfile("eee.txt");
	std::ofstream asymmputfile("ooo.txt");
	//ofstream outputfile(Form("%.0lfns%dX%dUVetadepMu%d.txt",window,nX,nUV,MU));
	//ofstream asymmputfile(Form("%.0lfns%dX%dUVetadepMu%dAsymm.txt",window,nX,nUV,MU));
	for(int cc=0;cc<etavec.size();cc++){
	  cout<<"ntruth"<<Ntruth[cc][defbb][aa][tt]<<endl;
	  // etastart = etavec.at(cc)-etapitch/2.0;
	  //etastop = etavec.at(cc)+etapitch/2.0;
	  //cout<<etastart<<" "<<etastop<<endl;
	  cout<<fixed<<setprecision(6)<<etavec.at(cc)<<" "<<eff[cc][defbb][aa][tt]<<" "<<0<<" "<<effError[cc][defbb][aa][tt]<<endl;
	  outputfile<<fixed<<setprecision(6)<<etavec.at(cc)<<" "<<eff[cc][defbb][aa][tt]<<" "<<0<<" "<<effError[cc][defbb][aa][tt]<<endl;
	  cout<<fixed<<setprecision(6)<<etavec.at(cc)<<" "<<eff[cc][defbb][aa][tt]<<" "<<0<<" "<<0<<" "<<-ROOT::Math::beta_quantile(alpha/2.0,eff[cc][defbb][aa][tt]*Ntruth[cc][defbb][aa][tt]/100.0+1,Ntruth[cc][defbb][aa][tt]-eff[cc][defbb][aa][tt]*Ntruth[cc][defbb][aa][tt]/100.0+1)*100.0+eff[cc][defbb][aa][tt]<<" "<<ROOT::Math::beta_quantile_c(alpha/2.0,eff[cc][defbb][aa][tt]*Ntruth[cc][defbb][aa][tt]/100.0+1,Ntruth[cc][defbb][aa][tt]-eff[cc][defbb][aa][tt]*Ntruth[cc][defbb][aa][tt]/100.0+1)*100.0-eff[cc][defbb][aa][tt]<<endl;
	  asymmputfile<<fixed<<setprecision(6)<<etavec.at(cc)<<" "<<eff[cc][defbb][aa][tt]<<" "<<0<<" "<<0<<" "<<-ROOT::Math::beta_quantile(alpha/2.0,eff[cc][defbb][aa][tt]*Ntruth[cc][defbb][aa][tt]/100.0+1,Ntruth[cc][defbb][aa][tt]-eff[cc][defbb][aa][tt]*Ntruth[cc][defbb][aa][tt]/100.0+1)*100.0+eff[cc][defbb][aa][tt]<<" "<<ROOT::Math::beta_quantile_c(alpha/2.0,eff[cc][defbb][aa][tt]*Ntruth[cc][defbb][aa][tt]/100.0+1,Ntruth[cc][defbb][aa][tt]-eff[cc][defbb][aa][tt]*Ntruth[cc][defbb][aa][tt]/100.0+1)*100.0-eff[cc][defbb][aa][tt]<<endl;
	  /*truthrate[etavec.size()][loopvec.size()][tractnum][timevec.size()]={{{{0}}}};
	  trackrate[etavec.size()][loopvec.size()][tractnum][timevec.size()]={{{{0}}}};
	  fakerate[etavec.size()][loopvec.size()][tractnum][timevec.size()]={{{{0}}}};
	  eff[etavec.size()][loopvec.size()][tractnum][timevec.size()]={{{{0}}}};
	  resoleff[etavec.size()][loopvec.size()][tractnum][timevec.size()]={{{{0}}}};
	  effcc[etavec.size()][loopvec.size()][tractnum][timevec.size()]={{{{0}}}};
	  effError[etavec.size()][loopvec.size()][tractnum][timevec.size()]={{{{0}}}};
	  resoleffError[etavec.size()][loopvec.size()][tractnum][timevec.size()]={{{{0}}}};
	  effccError[etavec.size()][loopvec.size()][tractnum][timevec.size()]={{{{0}}}};*/
	}
	outputfile.close();
	asymmputfile.close();
      }
    }
  }
  cout<<endl;
  cout<<"Alert"<<endl;
  for(int d=0;d<alertnum;d++){
    cout<<alert[d]<<"\t";
  }    

}
