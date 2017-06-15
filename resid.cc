#include "/home/maekawa/atlasstyle-00-03-05/AtlasStyle.C"
#include "/home/maekawa/atlasstyle-00-03-05/AtlasLabels.C"
#include <iomanip>
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
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraph2D.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLegend.h"
#include "TPaveText.h"
using namespace std;

double Dphi(int sector,double phi){
  const double pi = TMath::Pi();
  double dphi = phi; 
  if(sector==0) dphi += pi/8.0;
  while(dphi>pi/8.0||dphi<=-pi/8.0){
    if(dphi>pi/8.0) dphi -= pi/4.0;
    if(dphi<=-pi/8.0) dphi += pi/4.0;
  }
  return dphi;
}


void resid(){
  SetAtlasStyle();
  cout<<Dphi(0,1.5)<<endl;
  const double pi = TMath::Pi();
  double bound1 = -0.05;
  double bound2 = 0.05;
  int nbin = 200;
  double dx = (bound2-bound1)/nbin;
  string varname = "";
  string filevarname = "";
  const char inex[100] = ">="; //inclusive/exclusive  
  const int resolon = 1; //for resolution or just seeing all tracks
  //in using "resolon = 1", the only one best track is selected for residual in each event
  vector<double> bound1s(30);
  vector<double> bound2s(30);
  vector<double> require1s(30);
  vector<double> require2s(30);
  vector<string> varnames(30);
  vector<string> filevarnames(30);
  vector<string> axisvarnames(30);
  vector<int> fitornot(30);
  vector<int> requireon(30);
  vector<string> units(30);
  const int deftt=1;//default time window
  const int defX=3;//default X
  const int defUV=3;//defaul UV
  const int omit = 1;//omit other than globaleta globalphi dtheta
  string geta = "2dgeta";
  string gphi = "3dgphi";
  varnames.erase(varnames.begin(),varnames.end());
  filevarnames.erase(filevarnames.begin(),filevarnames.end());
  axisvarnames.erase(axisvarnames.begin(),axisvarnames.end());
  bound1s.erase(bound1s.begin(),bound1s.end());
  bound2s.erase(bound2s.begin(),bound2s.end());
  require1s.erase(require1s.begin(),require1s.end());
  require2s.erase(require2s.begin(),require2s.end());
  fitornot.erase(fitornot.begin(),fitornot.end());
  requireon.erase(requireon.begin(),requireon.end());
  units.erase(units.begin(),units.end());


  varnames.push_back("ltheta-gtheta-hit_ltheta3d+hit_gtheta3d"); filevarnames.push_back("dtheta");        
  bound1s.push_back(-0.02);   bound2s.push_back(0.02);                                                   
  require1s.push_back(-0.001);   require2s.push_back(0.001);  
  axisvarnames.push_back("Residual of #Delta#theta"); units.push_back("[rad]");                                 
  fitornot.push_back(1);
  requireon.push_back(1);

  /*varnames.push_back("ltheta3d-gtheta3d-hit_ltheta3d+hit_gtheta3d"); filevarnames.push_back("dtheta3d");
  bound1s.push_back(-0.02);   bound2s.push_back(0.02);
  require1s.push_back(-0.02);   require2s.push_back(0.02);
  axisvarnames.push_back("#sigma(#Delta#theta)"); units.push_back("[rad]");
  fitornot.push_back(1);
  requireon.push_back(0);
  */

  /*varnames.push_back("gtheta-hit_gtheta"); filevarnames.push_back("2dg");
  bound1s.push_back(-0.0002);   bound2s.push_back(0.0002);
  require1s.push_back(-0.0002);   require2s.push_back(0.0002);
  axisvarnames.push_back("#sigma(#theta_{global})"); units.push_back("[rad]");
  fitornot.push_back(1);
  requireon.push_back(0);*/

  varnames.push_back("ltheta-hit_ltheta"); filevarnames.push_back("2dl");
  bound1s.push_back(-0.02);   bound2s.push_back(0.02);
  require1s.push_back(-0.02);   require2s.push_back(0.02);
  axisvarnames.push_back("#sigma(#theta_{local})"); units.push_back("[rad]");
  fitornot.push_back(1);
  requireon.push_back(0);

  /*varnames.push_back("ltheta3d-hit_ltheta3d"); filevarnames.push_back("3dltheta");
  bound1s.push_back(-0.05);   bound2s.push_back(0.05);
  require1s.push_back(-0.05);   require2s.push_back(0.05);
  axisvarnames.push_back("#sigma(#theta_{local})"); units.push_back("[rad]");
  fitornot.push_back(1);
  requireon.push_back(0);*/

  /*varnames.push_back("asin(sin(lphi3d)*cos(hit_lphi3d)-cos(lphi3d)*sin(hit_lphi3d))"); filevarnames.push_back("3dlphi");
  bound1s.push_back(-0.05);   bound2s.push_back(0.05);
  require1s.push_back(-0.05);   require2s.push_back(0.05);
  axisvarnames.push_back("#sigma(#phi_{local})"); units.push_back("[rad]");
  fitornot.push_back(1);
  requireon.push_back(0);*/

  varnames.push_back("gtheta3d-hit_gtheta3d"); filevarnames.push_back("3dgtheta");
  bound1s.push_back(-0.0005);   bound2s.push_back(0.0005);
  require1s.push_back(-0.0005);   require2s.push_back(0.0005);
  axisvarnames.push_back("#sigma(#theta_{global})"); units.push_back("[rad]");
  fitornot.push_back(0);
  requireon.push_back(0);

  varnames.push_back("-log(tan(gtheta/2.0))+log(tan(hit_gtheta/2.0))"); filevarnames.push_back("2dgeta");
  filevarnames.push_back("2dgeta");
  bound1s.push_back(-0.0005);   bound2s.push_back(0.0005);
  require1s.push_back(-0.005);   require2s.push_back(0.005);
  axisvarnames.push_back("Residual of #eta"); units.push_back("");
  fitornot.push_back(1);
  requireon.push_back(1);

  /*varnames.push_back("-log(tan(atan(tan(gtheta)/cos(Dphi(sector,hit_gphi3d)))/2.0))+log(tan(atan(tan(hit_gtheta)/cos(Dphi(sector,gphi3d)))/2.0))");
  filevarnames.push_back("2dgeta");
  bound1s.push_back(-0.005);   bound2s.push_back(0.005);
  require1s.push_back(-0.005);   require2s.push_back(0.005);
  axisvarnames.push_back("Residual of #eta"); units.push_back("");
  fitornot.push_back(1);
  requireon.push_back(1);*/

  /*varnames.push_back("-log(tan(gtheta3d/2.0))+log(tan(hit_gtheta3d/2.0))"); filevarnames.push_back("3dgeta");
  bound1s.push_back(-0.002);   bound2s.push_back(0.002);
  require1s.push_back(-0.002);   require2s.push_back(0.002);
  axisvarnames.push_back("#sigma(#eta_{global})"); units.push_back("");
  fitornot.push_back(0);
  requireon.push_back(0);*/

  varnames.push_back("asin(sin(gphi3d)*cos(hit_gphi3d)-cos(gphi3d)*sin(hit_gphi3d))"); filevarnames.push_back("3dgphi");
  bound1s.push_back(-0.05);   bound2s.push_back(0.05);
  require1s.push_back(-0.02);   require2s.push_back(0.02);
  axisvarnames.push_back("Residual of #phi"); units.push_back("[rad]");
  fitornot.push_back(1);
  requireon.push_back(1);
 
  varnames.push_back("dr2d"); filevarnames.push_back("dr2d");
  bound1s.push_back(0);   bound2s.push_back(2);
  require1s.push_back(0);   require2s.push_back(1);
  axisvarnames.push_back("#sigma(LocalX)"); units.push_back("[mm]");
  fitornot.push_back(0);
  requireon.push_back(0);


  int tractnum = 5;
  int tract[tractnum][2]={{0}};
  tract[0][0]=2;tract[0][1]=1;
  tract[1][0]=2;tract[1][1]=2;
  tract[2][0]=3;tract[2][1]=2;
  tract[3][0]=3;tract[3][1]=3;
  tract[4][0]=4;tract[4][1]=4;

  double sigma[5][5][3]={{{0}}};
  double tractsigma[tractnum][3]={{0}};
  double tracttail[tractnum][3]={{0}};

  char namecore[120]="";
  strcpy(namecore,getenv("NAMECORE"));

  TF1 *func = new TF1("gaussian","[2]*exp(-(x-[0])*(x-[0])/2.0/[1]/[1])",-1,1);
  func->SetLineColorAlpha(2,0.7);
  TFile *file1 = new TFile(Form("%sd3.root",namecore));
  TTree *d3 = (TTree*)file1->Get("d3");
  char *ret;
  char namefordata[120] = "";
  ret = strstr(Form("%s",namecore),"drcut");
  if(ret!=NULL){
    int drcutlength = strlen(Form("%s",ret));
    int namecorelength = strlen(Form("%s",namecore));
    strncat(namefordata,Form("%s",namecore),namecorelength-drcutlength);
  }else{
    strcat(namefordata,Form("%s",namecore));
  }

  TFile *file2 = new TFile(Form("%s.root",namefordata));
  TTree *data = (TTree*)file2->Get("NSWHitsTree");

  TCanvas *c1 = new TCanvas("","");
  TCanvas *c[5][5][3];
  TH1F *h;

  for(int pp=0;pp<varnames.size();pp++){
    varname = varnames.at(pp);
    filevarname = filevarnames.at(pp);
    bound1 = bound1s.at(pp);
    bound2 = bound2s.at(pp);
    //initialize
    for(int j=0;j<tractnum;j++){
      for(int s=0;s<3;s++){
	tractsigma[j][s]=0;
	tracttail[j][s]=0;
      }
    }
    
    //for each
    if(requireon.at(pp)==1){
    TLine *l1[filevarnames.size()];
    TLine *l2[filevarnames.size()];

   for(int s=deftt;s<deftt+1;s++){
     for(int nX=defX;nX<defX+1;nX++){
      for(int nUV=defUV;nUV<defUV+1;nUV++){
	c1->cd();
	if(resolon==1){
          d3->Draw(Form("%s>>h(%d,%lf,%lf)",varname.c_str(),nbin,bound1,bound2),Form("etadr0==MinIf$(etadr0,%s)&&timewindow==%d&&fakeornot==0",Form("nX%s%d&&nUV%s%d&&timewindow==%d&&fakeornot==0",inex,nX,inex,nUV,25*(s+1)),25*(s+1)),"goff");
	}else{
	  d3->Draw(Form("%s>>h(%d,%lf,%lf)",varname.c_str(),nbin,bound1,bound2),Form("nX%s%d&&nUV%s%d&&timewindow==%d",inex,nX,inex,nUV,25*(s+1)),"goff");
	}
	TLegend *fit = new TLegend(0.6,0.6,0.85,0.7,"","nbNDC");
	fit->SetBorderSize(0);
	fit->SetFillColorAlpha(0,0);
	h = (TH1F*)gDirectory->Get("h");
	if(h->GetEntries()!=0){
	  h->Scale(1.0/h->GetEntries());
	  c[nX][nUV][s] = new TCanvas(Form("%dX%dUV-%dns",nX,nUV,25*(s+1)),Form("%dX%dUV-%dns",nX,nUV,25*(s+1)));
	  h->SetTitle(Form("%dX%dUV-%dns %s",nX,nUV,25*(s+1),varname.c_str()));
	  sigma[nX][nUV][s] = sqrt(h->GetRMS()*h->GetRMS()+h->GetMean()*h->GetMean());
	  if(h->GetEntries()>500 && fitornot.at(pp)==1){
	    func->SetParameters((bound1+bound2)/2.0,(bound2-bound1)/20.0,h->GetBinContent(h->GetMaximumBin()),pi);
	    func->SetParLimits(0,-fabs(bound2-bound1)/10.0,fabs(bound2-bound1)/10.0);
	    func->SetParLimits(1,0,fabs(bound2-bound1)/5.0);
	    func->SetParLimits(2,h->GetBinContent(h->GetMaximumBin())*0.9,h->GetBinContent(h->GetMaximumBin())*1.1);
	    func->FixParameter(3,pi);
	    h->Fit("gaussian","","",bound1,bound2);
	    if(Form("%s",filevarnames.at(pp).c_str())=="dtheta")h->Fit("gaussian","","",-0.002,0.002);
	    func->ReleaseParameter(3);
	    cout<<func->GetChisquare()/func->GetNDF()<<endl;
	    //if(func->GetChisquare()/func->GetNDF()<300&&func->GetChisquare()/func->GetNDF()>0.5)sigma[nX][nUV][s] = func->GetParameter(1);
	    sigma[nX][nUV][s] = func->GetParameter(1);  
	  }
	  c[nX][nUV][s]->cd();
	  gStyle->SetOptStat("");
	  h->GetXaxis()->SetTitle(Form("%s %s",axisvarnames.at(pp).c_str(),units.at(pp).c_str()));
	  if(filevarnames.at(pp).c_str()==geta){
	    h->GetYaxis()->SetTitle(Form("entries / %.1lf#times10^{-5}",(bound2-bound1)*100000.0/nbin));
	  }else{
	    h->GetYaxis()->SetTitle(Form("entries / %.1lf mrad",(bound2-bound1)*1000.0/nbin));
	  }
	  //h->SetLineWidth(3);
	  //func->SetLineColorAlpha(2,0.7);
	  //func->SetLineStyle(2);
	  //func->SetLineWidth(3);
	  if(nX==3&&nUV==3&&s==1){
	    cout<<"AAAAA"<<filevarnames.at(pp)<<endl;}
	  cout<<sigma[nX][nUV][s]<<endl;
	  if(h->GetEntries()>500 && fitornot.at(pp)==1){ 
	    cout<<"tail"<<(1-(h->Integral(h->FindBin(func->GetParameter(0)-3*sigma[nX][nUV][s]),h->FindBin(func->GetParameter(0)+3*sigma[nX][nUV][s])))/1.0)*100<<endl;
	  }else{
	    cout<<"tail"<<(1-(h->Integral(h->FindBin(h->GetMean()-3*sigma[nX][nUV][s]),h->FindBin(h->GetMean()+3*sigma[nX][nUV][s])))/1.0)*100<<endl;
	  }
	  cout<<"mean"<<func->GetParameter(0)<<endl;
	  gStyle->SetOptStat("");
          //if(pp==9)gStyle->SetOptStat("nemr");
	  gPad->SetBottomMargin(0.2);
	  gPad->SetRightMargin(0.15);
	  //h->GetXaxis()->SetLabelSize(0.07);
	  //h->GetXaxis()->SetTitleSize(0.07);
	  h->GetXaxis()->SetNdivisions(509);
	  h->Draw();
	  if(requireon.at(pp)==1){
	    l1[pp] = new TLine(func->GetParameter(0)+require1s.at(pp),0,func->GetParameter(0)+require1s.at(pp),1.01*h->GetBinContent(h->GetMaximumBin()));
	    l2[pp] = new TLine(func->GetParameter(0)+require2s.at(pp),0,func->GetParameter(0)+require2s.at(pp),1.01*h->GetBinContent(h->GetMaximumBin()));
	    l1[pp]->SetLineColor(2);
	    l2[pp]->SetLineColor(2);
	    func->SetLineColorAlpha(2,0.7);
	    //l1[pp]->Draw("same");
	    //l2[pp]->Draw("same");
	    cout<<"out of requirements"<<(1-(h->Integral(h->FindBin(func->GetParameter(0)+require1s.at(pp)),h->FindBin(func->GetParameter(0)+require2s.at(pp))))/1.0)*100<<endl;
	  }
	  func->SetMarkerSize(0);

	  TLine *l = new TLine(0,1,0,1);
	  l->SetLineColor(2); 
	  fit->AddEntry(func,Form("Gaussian fit"),"l");
	  double below = 0.47;
	  //if(filevarnames.at(pp).c_str()==gphi){ below = 0.39; }
	  TPaveText *sig = new TPaveText(0.6,below,0.8,0.55,"NBNDC");
	  if(filevarnames.at(pp).c_str()==geta){
	    sig->AddText(Form("#sigma = %.1lf#times10^{-5}",sigma[nX][nUV][s]*100000.0));
	  }else{
	    sig->AddText(Form("#sigma = %.1lf mrad",sigma[nX][nUV][s]*1000.0));
	    /*if(filevarnames.at(pp).c_str()==gphi){
	      sig->AddText(Form("#mu = %.1lf mrad",func->GetParameter(0)*1000.0));
	    }*/
	  }
	  sig->SetFillColor(0);
	  fit->Draw();
	  sig->Draw("same");
	  c[nX][nUV][s]->Print(Form("resid%s%dX%dUV-%dns.pdf",filevarname.c_str(),nX,nUV,25*(s+1)));
	}
      }
    }
  }
    
   /*
  TH1D *xuvlhist[3];
  for(int s=0; s<3; s++){
    xuvlhist[s] = new TH1D(Form("xuvldist%d",s+1),Form("%dns",(s+1)*25),5,-0.5,4.5);
  }
  TCanvas *c2 = new TCanvas("c2","c2");
  c2->Divide(3,1);

  gStyle->SetOptStat(000000000);     
  for(int s=0;s<3;s++){
    for(int j=0;j<tractnum;j++){
      xuvlhist[s]->GetXaxis()->SetBinLabel(j+1,Form("%dX,%dUV",tract[j][0],tract[j][1]));
       if(sigma[tract[j][0]][tract[j][1]][s]!=0){
	xuvlhist[s]->SetBinContent(j+1,sigma[tract[j][0]][tract[j][1]][s]);
	xuvlhist[s]->SetBinError(j+1,0.00000000001);
      }else{
	xuvlhist[s]->SetBinContent(j+1,0);
	xuvlhist[s]->SetBinError(j+1,0);
      }
    }
    c2->cd(s+1);
    xuvlhist[s]->GetXaxis()->SetTitle("CT");
    xuvlhist[s]->GetXaxis()->SetTitleOffset(1);
    xuvlhist[s]->GetXaxis()->SetLabelSize(0.06);
    xuvlhist[s]->GetXaxis()->SetTitleSize(0.04);
    xuvlhist[s]->GetYaxis()->SetTitle("sigma_ltheta");
    xuvlhist[s]->GetYaxis()->SetTitleOffset(1.2);
    xuvlhist[s]->GetYaxis()->SetTitleSize(0.04);
    //gStyle->SetOptStat(000000000);
    xuvlhist[s]->SetTitleSize(0.20);
    xuvlhist[s]->Draw("E1");
  }

  if(resolon==1){
    c2->Print(Form("resid%s%sresol.pdf",filevarname.c_str(),namecore));
  }else{
    c2->Print(Form("resid%s%s.pdf",filevarname.c_str(),namecore));
  }

  for(int s=0;s<3;s++){
    cout<<"s="<<25*(s+1)<<"ns"<<endl;
    cout<<"\tnX2\tnX3\tnX4"<<endl;
    for(int nUV=1;nUV<5;nUV++){
      cout<<"nUV"<<nUV;
      for(int nX=2;nX<5;nX++){
	cout<<"\t"<<sigma[nX][nUV][s];
      }
      cout<<endl;
    }
  }
   */
    }
  //for each
    if(omit==0){
  TH1F *h1[tractnum][3];
  TCanvas *c4;
  c4 = new TCanvas("c4","c4",2000,900);
  c4->Divide(tractnum,3);
  for(int s=0;s<3;s++){
    for(int j=0;j<tractnum;j++){
      //h1[j][s] = new TH1F(Form("%dX%dUV%dns",tract[j][0],tract[j][1],25*(s+1)),Form("%dX%dUV%dns",tract[j][0],tract[j][1],25*(s+1)),nbin,bound1,bound2);
      c1->cd();
      if(resolon==1){
	d3->Draw(Form("%s>>h1[%d][%d](%d,%lf,%lf)",varname.c_str(),j,s,nbin,bound1,bound2),Form("etadr0==MinIf$(etadr0,%s)&&timewindow==%d&&-log(tan(theta0/2.0))>1.2&&-log(tan(theta0/2.0))<2.8",Form("nX%s%d&&nUV%s%d&&timewindow==%d",inex,tract[j][0],inex,tract[j][1],25*(s+1)),25*(s+1)),"goff");
      }else{
	d3->Draw(Form("%s>>h1[%d][%d](%d,%lf,%lf)",varname.c_str(),j,s,nbin,bound1,bound2),Form("nX%s%d&&nUV%s%d&&timewindow==%d",inex,tract[j][0],inex,tract[j][1],25*(s+1)),"goff");
      }
      h1[j][s] = (TH1F*)gDirectory->Get(Form("h1[%d][%d]",j,s));
      if(h1[j][s]->GetEntries()!=0){
	h1[j][s]->SetTitle(Form("%dX%dUV-%dns",tract[j][0],tract[j][1],25*(s+1)));
	h1[j][s]->SetTitleSize(0.20);
	if(h1[j][s]->GetEntries()>300 && fitornot.at(pp)==1){
	  func->SetParameters((bound1+bound2)/2.0,(bound2-bound1)/20.0,h1[j][s]->GetBinContent(h1[j][s]->GetMaximumBin()),pi);
	  func->SetParLimits(0,-fabs(bound2-bound1)/10.0,fabs(bound2-bound1)/10.0);
	  func->SetParLimits(1,0,fabs(bound2-bound1)/5.0);
	  func->SetParLimits(2,h1[j][s]->GetBinContent(h1[j][s]->GetMaximumBin())*0.9,h1[j][s]->GetBinContent(h1[j][s]->GetMaximumBin())*1.1);
	  func->FixParameter(3,pi);
	  h1[j][s]->Fit("gaussian","","",bound1,bound2);
	  //if(pp==0)h1[j][s]->Fit("gaussian","","",-0.002,0.002);
	  func->ReleaseParameter(3);
	  cout<<func->GetChisquare()/func->GetNDF()<<endl;
	  //if(func->GetChisquare()/func->GetNDF()<300&&func->GetChisquare()/func->GetNDF()>0.5)
	  tractsigma[j][s] = func->GetParameter(1);
	  cout<<tractsigma[j][s]<<endl;
	  tracttail[j][s] = (1-(h1[j][s]->Integral(h1[j][s]->FindBin(func->GetParameter(0)-3*tractsigma[j][s]),h1[j][s]->FindBin(func->GetParameter(0)+3*tractsigma[j][s])))/1.0/h1[j][s]->GetEntries())*100;
	  cout<<25*(s+1)<<"ns"<<tract[j][0]<<"X"<<tract[j][1]<<"UV"<<endl;
	  cout<<h1[j][s]->Integral(h1[j][s]->FindBin(func->GetParameter(0)-3*tractsigma[j][s]),h1[j][s]->FindBin(func->GetParameter(0)+3*tractsigma[j][s]))<<endl;
	  cout<<h1[j][s]->Integral(h1[j][s]->FindBin(func->GetParameter(0)-4*tractsigma[j][s]),h1[j][s]->FindBin(func->GetParameter(0)+4*tractsigma[j][s]))<<endl;
	  cout<<h1[j][s]->Integral(h1[j][s]->FindBin(func->GetParameter(0)-5*tractsigma[j][s]),h1[j][s]->FindBin(func->GetParameter(0)+5*tractsigma[j][s]))<<endl;
	  cout<<tracttail[j][s]<<endl;
	}else{
	  tractsigma[j][s] = h1[j][s]->GetRMS();
          cout<<tractsigma[j][s]<<endl;
	  tracttail[j][s] = (1-(h1[j][s]->Integral(h1[j][s]->FindBin(h1[j][s]->GetMean()-3*tractsigma[j][s]),h1[j][s]->FindBin(h1[j][s]->GetMean()+3*tractsigma[j][s])))/1.0/h1[j][s]->GetEntries())*100;
          cout<<tracttail[j][s]<<endl;
	}
	c4->cd(j+1+s*tractnum);
	//gStyle->SetOptStat("nemiou");
	func->SetLineColorAlpha(2,0.5);
	h1[j][s]->Draw();
      }
    }
    if(resolon==1){
      c4->Print(Form("resid%sdist%sresol.pdf",filevarname.c_str(),namecore));
    }else{
      c4->Print(Form("resid%sdist%s.pdf",filevarname.c_str(),namecore));
    }
  }

  cout<<filevarname.c_str()<<endl;
  cout<<"\t";
  for(int j=0;j<tractnum;j++){
    cout<<Form("%dX%dUV \t",tract[j][0],tract[j][1]);
  }
  cout<<endl;
  for(int s=0;s<3;s++){
    cout<<Form("%dns",25*(s+1));
    for(int j=0;j<tractnum;j++){
      cout<<"\t"<<fixed<<setprecision(-floor(log10(tractsigma[j][s]))-2)<<1000.0*tractsigma[j][s]<<fixed<<setprecision(-floor(log10(tracttail[j][s]))+1)<<"(tail"<<tracttail[j][s]<<"%)";
    }
    cout<<endl;
  }
    }

  }
}
