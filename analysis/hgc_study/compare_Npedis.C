#include <string>
#include <iterator>
#include <algorithm>
#include <iostream> 
#include <fstream>
#include <cmath> 
#include <stdio.h>
#include <vector>
#include <math.h> 
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TPaveText.h>
#include <TText.h>
#include <TSystem.h>
#include <TArc.h>

using namespace std;

void compare_Npedis()
{
gROOT->Reset();
// gStyle->SetPalette(57);
gStyle->SetOptStat(0);
gStyle->SetPadColor(0);
  
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin(0.1);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.15);  
  
  gStyle->SetLabelSize(0.05,"xyz"); // size of axis values
  gStyle->SetTitleSize(0.06,"xyz");  
  gStyle->SetTitleOffset(0.7,"y");
  gStyle->SetTitleOffset(1,"x");    
  gStyle->SetTitleSize(0.07,"t");   

// char* label[4]={
// "alone,#pi-,z=0cm,#theta=8deg",
// "alone,#pi-,z=0cm,#theta=14.8deg old",    
// "alone,#pi-,z=0cm,#theta=8deg",    
// "alone,#pi-,z=0cm,#theta=14.8deg old",    
// };

// const int n=4;

// char* input_filename[n]={
// "output/hgc_moved/background3sectordouble/solid_SIDIS_He3_hgc_moved_pim_z350_p2.5_theta8.0_1e5_output.root",
// "output/hgc_moved/background3sectordouble/solid_SIDIS_He3_hgc_moved_km_z350_p2.5_theta8.0_1e5_output.root",
// "output/hgc_moved/background3sectordouble/solid_SIDIS_He3_hgc_moved_pim_z350_p7.5_theta14.5_1e5_output.root",
// "output/hgc_moved/background3sectordouble/solid_SIDIS_He3_hgc_moved_km_z350_p7.5_theta14.5_1e5_output.root",
// };
// char* label[2]={
// "#pi^{-}",
// "K^{-}",    
// };

// // char* dir1="/volatile/halla/solid/sim/solid_gemc/SIDIS_He3_JLAB_VERSION_2.5/pass3/";
// char* dir1="solid_SIDIS_He3_full_lgc/";
// 
// // char* dir2="backgroundno_nodecaytrack/";
// // char* dir2="background3sector_nodecaytrack/";
// char* dir2="background3sectordouble_nodecaytrack/";
// 
// const int n=6;
// char* input_filename[n]={
// "solid_SIDIS_He3_moved_full_ele_z350_p2.5_theta8.0_1e6_output.root",
// "solid_SIDIS_He3_moved_full_pim_z350_p2.5_theta8.0_1e6_output.root",
// "solid_SIDIS_He3_moved_full_ele_z350_p3.0_theta12.0_1e6_output.root",
// "solid_SIDIS_He3_moved_full_pim_z350_p3.0_theta12.0_1e6_output.root",
// "solid_SIDIS_He3_moved_full_ele_z350_p4.5_theta14.5_1e6_output.root",
// "solid_SIDIS_He3_moved_full_pim_z350_p4.5_theta14.5_1e6_output.root",
// };
// char* label[2]={
// "e^{-}",
// "#pi^{-}",    
// };

char* dir1="./";

char* dir2="";

const int n=5;
char* input_filename[n]={
"solid_PVDIS_LD2_moved_full_ele_z10_p2.5_theta22_1e6_output_allsector.root",
"solid_PVDIS_LD2_moved_full_ele_z10_p2.5_theta22_1e6_output_3sector.root",
"solid_PVDIS_LD2_moved_full_ele_z10_p2.5_theta22_1e6_output_2sectorsmaller.root",
"solid_PVDIS_LD2_moved_full_ele_z10_p2.5_theta22_1e6_output_2sectorlarger.root",
"solid_PVDIS_LD2_moved_full_ele_z10_p2.5_theta22_1e6_output_1sector.root",
};
char* label[n]={
"all","3","3","3","3",
};


TCanvas *c = new TCanvas("compare","compare",1600,800);
TFile *input[n];
TH1F *h[n];
TLegend* leg= new TLegend(0.85, 0.75, 0.95, 0.9);
for(int j=0;j<n;j++){
  c->cd(j+1); 
  gPad->SetLogy(1);
  
  input[j]=new TFile(Form("%s%s%s",dir1,dir2,input_filename[j]));
  if (input[j]->IsZombie()) {
    cout << "Error opening file " << input_filename[j] << endl;
    exit(-1);
  }
  else cout << "open file " << input_filename[j] << endl;
  
//   char hstname[100];
//   sprintf(hstname,"%s_%i_%i",hst[i],hit_id[i],pid[i]);    
//   cout << hstname << endl;
  h[j]=(TH1F*) input[j]->Get("hnpe");

  h[j]->SetAxisRange(1,1e6,"Y");    
//   h[2*j+i]->SetAxisRange(0,100,"X");  
  h[j]->SetTitle(";Npe;counts");  
  h[j]->SetLineColor(j+1);    
  
  if (j==0) h[j]->Draw("E");
  else h[j]->Draw("same E");
  
  leg->AddEntry(h[j],label[j],"l");  
}
leg->Draw();

c->SaveAs("compare.png");

}