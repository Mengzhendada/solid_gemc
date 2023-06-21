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
#include <TROOT.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TPaveText.h>
#include <TText.h>
#include <TSystem.h>
#include <TArc.h>

using namespace std;

void compare()
{
  gROOT->Reset();
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
// const int m=4;
// char* input_filename[m]={
// "output_theta8_z0_blockoff_fieldoff_1e4_output.root",
// "output_theta8_z0_blockon_fieldoff_1e4_output.root",
// "output_theta8_z0_blockoff_fieldon_1e4_output.root",
// "output_theta8_z0_blockon_fieldon_1e4_output.root",
// };
// 
// char* label[m]={
// "blockoff and fieldoff",
// "blockon and fieldoff",
// "blockoff and fieldon",
// "blockon and fieldon",
// };
  
    const int m=2;
char* input_filename[m]={
"output_theta14.8_z0_blockoff_fieldoff_1e4_output.root",
"output_theta8_z0_blockoff_fieldoff_1e4_output.root",
};

char* label[m]={
"theta14.8_z0_blockoff_fieldoff",
"theta8_z0_blockoff_fieldoff",
};
int color[m]={4,1};

//   const int m=2;
// char* input_filename[m]={
// "output_theta8_z0_blockoff_fieldoff_1e4_output.root",
// "output_theta8_z0_blockoff_fieldon_1e4_output.root",
// };
// 
// char* label[m]={
// "theta8_z0_blockoff_fieldoff",
// "theta8_z0_blockoff_fieldon",
// };
// int color[m]={1,2};
  
//     const int m=2;
// char* input_filename[m]={
// "output_theta8_z0_blockoff_fieldon_1e4_output.root",
// "output_theta8_z0_blockon_fieldon_1e4_output.root",
// };
// 
// char* label[m]={
// "theta8_z0_blockoff_fieldon",
// "theta8_z0_blockon_fieldon",
// };
// int color[m]={2,3};

TCanvas *c = new TCanvas("compare","compare",1000,800);

TH1F *hcount_total_p=new TH1F("hcount_total_p","photoelectron;P (GeV);count",600,2,8);
hcount_total_p->SetAxisRange(2,8,"X");
hcount_total_p->SetAxisRange(0,50,"Y");   
hcount_total_p->Draw();

TFile *input[m];
TGraph *g[m];  
for(int i=0;i<m;i++){
//     if(i==3 || i>6) continue;
  input[i]=new TFile(input_filename[i]);
  if (input[i]->IsZombie()) {
    cout << "Error opening ratefile " << input_filename[i] << endl;
    exit(-1);
  }
  else cout << "open file " << input_filename[i] << endl;
  
//   c->cd(i+1);
  
//   char hstname[100];
//   sprintf(hstname,"%s_%i_%i",hst[i],hit_id[i],pid[i]);    
//   cout << hstname << endl;
  g[i]=(TGraph*) input[i]->Get("gcount_total_p");
//   h[i]->SetAxisRange(ymin,ymax,"Y");    
//   h[i]->SetAxisRange(xmin,xmax,"X");  
//   h[i]->SetTitle(Form("%s %s",title,label[i]));  
  g[i]->SetMarkerStyle(8);
  g[i]->SetMarkerSize(0.15*(m-i));    
  g[i]->SetMarkerColor(color[i]);
  g[i]->SetLineColor(color[i]);
  g[i]->Draw("same P");
  
//   cout << h[i]->Integral() << endl;
//     input.Close();   
}
TLegend* leg = new TLegend(0.7, 0.9-0.05*m, 0.95, 0.95);
for(int i=0;i<m;i++){
//   leg->AddEntry(h[i], Form("%s   %02f",label[i],h[i]->Integral()),"l");  
  leg->AddEntry(g[i], label[i],"l");    
}
leg->Draw();

}