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

void compare_1D()
{
  gROOT->Reset();
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
// const int m=4;
// char* input_filename[m]={
// "output_z0_theta8_blockoff_fieldoff_1e4_output.root",
// "output_z0_theta8_blockon_fieldoff_1e4_output.root",
// "output_z0_theta8_blockoff_fieldon_1e4_output.root",
// "output_z0_theta8_blockon_fieldon_1e4_output.root",
// };
// 
// char* label[m]={
// "blockoff and fieldoff",
// "blockon and fieldoff",
// "blockoff and fieldon",
// "blockon and fieldon",
// };

const int m=11;  
char* var[m]={"2.5","3.0","3.5","4.0","4.5","5.0","5.5","6.0","6.5","7.0","7.5"};    

const int n=1;
string input_filename[n][m];
for(int i=0;i<m;i++) {
  char name[200]; 
// pmtmove10cm_cone10cmtilt65deg
  
  sprintf(name,Form("data/JLAB_VERSION_1.3/NH3_shift20cmdown/output_pim_z0_p%s_blockoff_1e6.root",var[i]));  
  input_filename[0][i]=name;   
}
char* label[n]={
"(new) alone,#pi-,block off,z=0cm",  
};
int MarkerStyle[n]={4,4,4,26,26,26};
int color[n]={1,2,3,1,2,3};

const int N_bin_theta=50;
const int N_bin_phi=180;

TH1F *hcount_total_p[n];
for(int j=0;j<n;j++){
// hcount_total_p[j]=new TH1F(Form("hcount_total_p_%i",j),"photoelectron;P (GeV);count",m,1.75,8.25);
hcount_total_p[j]=new TH2F(Form("hcount_total_pe_%i",j),";P (GeV);photoelectron count (sim*0.5)",m,2.25,7.75);  
// hcount_total_p->SetAxisRange(2,8,"X");
hcount_total_p[j]->SetAxisRange(0,50,"Y");   
hcount_total_p[j]->SetMarkerStyle(MarkerStyle[j]);
hcount_total_p[j]->SetMarkerSize(2);
hcount_total_p[j]->SetMarkerColor(color[j]);
hcount_total_p[j]->SetLineColor(color[j]);
}

TCanvas *c = new TCanvas("compare","compare",1600,800);
c->Divide(n,1);
TFile *input[n][m];

TH1F *h[n][m][N_bin_theta][N_bin_phi];
TLegend* leg[n];
for(int j=0;j<n;j++){
  c->cd(j+1);  
  leg[j]= new TLegend(0.7, 0.95-0.05*m, 0.95, 0.95);  
for(int i=0;i<m;i++){
//   cout << j << " " << i << endl;  
//     if(i==3 || i>6) continue;
  input[j][i]=new TFile(input_filename[j][i].c_str());
//   cout << " " << input_filename[j][i] << endl;
  if (input[j][i]->IsZombie()) {
    cout << "Error opening ratefile " << input_filename[j][i] << endl;
    exit(-1);
  }
  else cout << "ok open file " << input_filename[j][i] << endl;
  
for(int bin_theta=0;bin_theta<50;bin_theta++){
  for(int bin_phi=0;bin_phi<180;bin_phi++){  
//   char hstname[100];
//   sprintf(hstname,"%s_%i_%i",hst[i],hit_id[i],pid[i]);    
//   cout << hstname << endl;
  h[j][i][bin_theta][bin_phi]=(TH1F*) input[j][i]->Get("hcount");
//   h[i]->SetAxisRange(ymin,ymax,"Y");    
//   h[i]->SetAxisRange(xmin,xmax,"X");  
//   h[i]->SetTitle(Form("%s %s",title,label[i]));
  h[j][i]->SetLineColor(i);  
  if (i==0) h[j][i]->Draw();  
  else h[j][i]->Draw("same");
//   h[i]->SetMarkerStyle(8);
//   g[i]->SetMarkerSize(0.15*(m-i));    
//   g[i]->SetMarkerColor(color[i]);
//   g[i]->SetLineColor(color[i]);
//   g[i]->Draw("same P");
  hcount_total_p[j]->SetBinContent(i+1,h[j][i]->GetMean());
  hcount_total_p[j]->SetBinError(i+1,h[j][i]->GetRMS());  
  
//   cout << h[i]->Integral() << endl;
//     input.Close();
//   leg->AddEntry(h[i], Form("%s   %02f",input_filename[i],h[i]->GetMean()),"l");  
  leg[j]->AddEntry(h[j][i], Form("%02f %02f",h[j][i]->GetMean(),h[j][i]->GetRMS()),"l");    
//   leg->AddEntry(g[i], label[i],"l");    
  
}
leg[j]->Draw();
}


TLegend* legend;
  legend= new TLegend(0.5, 0.98-0.05*n, 0.95, 0.98);  
TCanvas *c1 = new TCanvas("c","c",1000,800);
for(int j=0;j<n;j++){
if (j==0) hcount_total_p[j]->Draw("E1");
else hcount_total_p[j]->Draw("E same");
  legend->AddEntry(hcount_total_p[j], Form("%s",label[j]),"lep");    
}
legend->Draw();


}