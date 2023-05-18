#include <iomanip>
#include <cstdlib>
#include <iostream> 
#include <fstream> 
#include <string>

#include "TTree.h"
#include "TStyle.h"
#include "TFile.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPad.h"

using namespace std;

void compare_e_pibggen(string type,string filename_e,string filename_h,string hadron="pim"){
   gStyle->SetPalette(kBird);
   gStyle->SetOptStat(0);
  gStyle->SetLabelSize(0.04,"xyz"); // size of axis values
  gStyle->SetTitleSize(0.04,"xyz");   
  gStyle->SetTitleSize(0.07,"t");    
   gStyle->SetPaintTextFormat("4.1f");   

   char expname[200];
   double min_theta,max_theta,min_p,max_p;
   double min_Q2=0,max_Q2=14,min_x=0,max_x=1;   
    double rej_FA_low,rej_FA_high,rej_LA;    
   if (type=="solid_JPsi_LH2") {
     sprintf(expname, "JPsi");     
     min_theta=8;max_theta=28;min_p=0.5;max_p=7;
     rej_FA_low=1e-4;
     rej_FA_high=5e-3;     
     rej_LA=5e-3;
   }
   else if (type=="solid_SIDIS_He3") {
     sprintf(expname, "3He");
     min_theta=8;max_theta=24;min_p=1;max_p=7;     
     rej_FA_low=1e-4;
     rej_FA_high=5e-3;     
     rej_LA=5e-3;     
   }
   else if (type=="solid_PVDIS_LD2") {
     sprintf(expname, "LD");
     min_theta=22;max_theta=35;min_p=1.5;max_p=5;   
     rej_FA_low=1e-4;
     rej_FA_high=5e-3;     
   }
   else {
     cout << "unknown type" << endl;
     return;
  }

//   cout << min_theta << " " << max_theta << " " << min_p << " " << max_p << " " << endl;
  
  string hadron_title,e_title,hadron_filename;
  if (hadron=="pim") {hadron_title="#pi^{-}"; e_title="e^{-}"; hadron_filename="pion_m";}
  else if (hadron=="pip") {hadron_title="#pi^{+}"; e_title="e^{+}"; hadron_filename="pion_p";}
  else if (hadron=="pi0") {hadron_title="#pi^{0}"; e_title="e^{-}"; hadron_filename="pion_0";}  
  else if (hadron=="p") {
    hadron_title="p";  e_title="e^{+}"; hadron_filename="proton";
    rej_FA_high=rej_FA_low;	//proton has much higher Cherenkov threshold
  }
  else {
     cout << "unknown hadron" << endl;
     return;
  }  
  

   //TFile *file_e=new TFile("/work/halla/solid/zwzhao/solid/solid_svn/solid/evgen/eicRate_20101102/output/rate_solid_PVDIS_LD2_eDIS_1e6.root");
   // TFile *file_e=new TFile("/work/halla/solid/zwzhao/solid/solid_svn/solid/evgen/eicRate_20101102/output/rate_solid_PVDIS_LD2_eDIS_50GeV_1e6.root");
   TFile *file_e=new TFile(filename_e.c_str());
   TTree *T_e = (TTree*) file_e->Get("T");

   TCanvas *cPTheta_e = new TCanvas("cPTheta_e","cPTheta_e",1800,1000);
   T_e->Draw("sqrt(px*px+py*py+pz*pz):theta>>hPTheta_e(20,5,45,20,0,10)","rate/1e3","text colz");       
   TH2F *hPTheta_e = (TH2F*)gPad->GetPrimitive("hPTheta_e");
   hPTheta_e->SetTitle("e^{-} rate (kHz);#theta (deg);P (GeV/c)");      
//    hPTheta_e->GetXaxis()->SetRangeUser(min_theta,max_theta);
//    hPTheta_e->GetYaxis()->SetRangeUser(min_p,max_p);
   cPTheta_e->SaveAs(Form("%s_rate_PTheta_e.pdf",type.c_str()));     
   {
   int bin,binmin_theta,binmin_p,binzmin,binmax_theta,binmax_p,binzmax;  
   bin=hPTheta_e->FindBin(min_theta,min_p);  
   hPTheta_e->GetBinXYZ(bin,binmin_theta,binmin_p,binzmin);   
   bin=hPTheta_e->FindBin(max_theta,max_p);  
   hPTheta_e->GetBinXYZ(bin,binmax_theta,binmax_p,binzmax);   
   cout << "e- rate in kHz "<< hPTheta_e->Integral(binmin_theta,binmax_theta,binmin_p,binmax_p) << endl;   
   }   

   TCanvas *cQ2x_e = new TCanvas("cQ2x_e","cQ2x_e",1800,1000);      
   T_e->Draw("Q2:x>>hQ2x_e(10,0,1,14,0,14)","rate/1e3","text colz");        
   TH2F *hQ2x_e = (TH2F*)gPad->GetPrimitive("hQ2x_e");
   hQ2x_e->SetTitle("e^{-} rate (kHz);x;Q^{2} (GeV/c)^{2}");      
//    hPTheta_e->GetXaxis()->SetRangeUser(min_theta,max_theta);
//    hPTheta_e->GetYaxis()->SetRangeUser(min_p,max_p);
   cQ2x_e->SaveAs(Form("%s_rate_Q2x_e.pdf",type.c_str()));     
   
   //get pi
   TH2F *hPTheta_pi=new TH2F("hPTheta_pi","hPTheta_pi",20,5,45,20,0,10);
   TH2F *hQ2x_pi=new TH2F("hQ2x_pi","hQ2x_pi",10,0,1,14,0,14);   
   TH2F *hPTheta_pi_rej=new TH2F("hPTheta_pi_rej","hPTheta_pi_rej",20,5,45,20,0,10);
   TH2F *hQ2x_pi_rej=new TH2F("hQ2x_pi_rej","hQ2x_pi_rej",10,0,1,14,0,14);   
   
   int filecount;
    if (type=="solid_JPsi_LH2" || type=="solid_SIDIS_He3" ) filecount=1000;   
    else if (type=="solid_PVDIS_LD2") filecount=1800;   
//     filecount=10;    
    
   for(int i=0; i<filecount;i++){
      string filename_hPTheta_name=Form("%s/LUND/hallD_%s_%s_1k_%i.lund",filename_h.c_str(),hadron_filename.c_str(),expname,i+1);
      ifstream input(filename_hPTheta_name.c_str());
      int eventcont=1000;
      if (input.good()) {} //cout << "open file " << filename_hPTheta_name << " OK" << "\r";
      else {cout << "can't open the file" << endl; return;}

      // reads and discard the first 13 lines of text
//       char textline[200];
//       for (int k = 0; k<13; k++) input.getline(textline,200);

      double  header[10];
      double  particles[14];

      for (int j = 0; j<eventcont; j++) {
//defined in https://github.com/JeffersonLab/evgen_bggen/blob/master/scripts/HallD_gen_lund.cc
// myfile[pid[111]]<< "1" << " \t " << Wprate  << " \t " << Wmrate  << " \t " << "0"  << " \t " << "0" << " \t "  << fEvXbj << " \t " << y  << " \t " << fEvW2  << " \t " << fEvQ2  << " \t " << targ_A[targ_index]*totalRate[pid[111]]/nLundEntries << endl;
//   myfile[pid[111]]<< " \t " << "1" << " \t " << pidcharge[111] << " \t " << "1" << " \t " << 111 << " \t " << "0" << " \t " << "0" << " \t " << pout[j][0] << " \t " << pout[j][1] << " \t " << pout[j][2] << " \t " << Ef << " \t " << pidmass[111] << " \t " << fEvV[0]*100  << " \t " << fEvV[1]*100 << " \t " << vz+targ_offset[targ_index]*100 << endl;	
// 1        0       0       0       0       0       0       0       0       231093
//          1       -1      1       -211    0       0       -0.0675592      0.440683        0.421785        0.629403        0.13957         0.0881298       -0.18	963      -311.584	
	      input >> header[0] >> header[1] >> header[2] >> header[3] >> header[4] >> header[5] >> header[6] >> header[7] >> header[8] >> header[9];
	      input >> particles[0] >> particles[1] >> particles[2] >> particles[3] >> particles[4] >> particles[5] >> particles[6] >> particles[7] >> particles[8] >> particles[9] >> particles[10] >> particles[11] >> particles[12] >> particles[13];
	      
// 	      cout << header[9] << " " << particles[6] << " " << particles[7] << " " << particles[8] << endl;
	      double rate=header[9];
	      double px=particles[6],py=particles[7],pz=particles[8];	      
	      double p=sqrt(px*px+py*py+pz*pz);
	      double theta=acos(pz/p)/3.1416*180;
	      double Q2=4*11*p*sin(theta/2/180*3.1416)*sin(theta/2/180*3.1416);
	      double x=Q2/2/0.938/(11-p);

	      double rate_rej=rate;
      	      if (type=="solid_JPsi_LH2" || type=="solid_SIDIS_He3" ) {		
		if (theta <15) {//FA
		  if (p<4) rate_rej=rate*rej_FA_low;
		  else  rate_rej=rate*rej_FA_high;
		}
		else rate_rej=rate*rej_LA;
	      }
	      else if (type=="solid_PVDIS_LD2") {      
		  if (p<4) rate_rej=rate*rej_FA_low;
		  else  rate_rej=rate*rej_FA_high;
	      }
 
	      hPTheta_pi->Fill(theta,p,rate);
	      hQ2x_pi->Fill(x,Q2,rate);	      
      	      hPTheta_pi_rej->Fill(theta,p,rate_rej);
	      hQ2x_pi_rej->Fill(x,Q2,rate_rej);	      
    }
    
   }
   
   cout << "pi bggen files " << filecount << endl;

    //in kHz and normalized by filecount   
   hPTheta_pi->Scale(1e-3/filecount);
   hQ2x_pi->Scale(1e-3/filecount);
   hPTheta_pi_rej->Scale(1e-3/filecount);  
   hQ2x_pi_rej->Scale(1e-3/filecount);     
   
   TCanvas *cPTheta_pi = new TCanvas("cPTheta_pi","cPTheta_pi",1800,1000);
   hPTheta_pi->SetTitle(Form("%s rate (kHz);#theta (deg);P (GeV/c)",hadron_title.c_str()));       
//    hPTheta_pi->GetXaxis()->SetRangeUser(min_theta,max_theta);   
//    hPTheta_pi->GetYaxis()->SetRangeUser(min_p,max_p);     
   hPTheta_pi->Draw("text colz");
//    hPTheta_pi->SetMaximum(1e13);
//    hPTheta_pi->SetMinimum(1e-3);
   cPTheta_pi->SaveAs(Form("%s_rate_PTheta_%s.pdf",type.c_str(),hadron.c_str()));         
   {
   int bin,binmin_theta,binmin_p,binzmin,binmax_theta,binmax_p,binzmax;  
   bin=hPTheta_pi->FindBin(min_theta,min_p);  
   hPTheta_pi->GetBinXYZ(bin,binmin_theta,binmin_p,binzmin);   
   bin=hPTheta_pi->FindBin(max_theta,max_p);  
   hPTheta_pi->GetBinXYZ(bin,binmax_theta,binmax_p,binzmax);   
   cout << "pi rate in kHz "<< hPTheta_pi->Integral(binmin_theta,binmax_theta,binmin_p,binmax_p) << endl;   
   }   

   TCanvas *cQ2x_pi = new TCanvas("cQ2x_pi","cQ2x_pi",1800,1000);
   hQ2x_pi->SetTitle(Form("%s rate (kHz);x;Q^{2} (GeV/c)^{2}",hadron_title.c_str()));       
   hQ2x_pi->Draw("text colz");
//    hQ2x_pi->SetMaximum(1e13);
//    hQ2x_pi->SetMinimum(1e-3);
   cQ2x_pi->SaveAs(Form("%s_rate_Q2x_%s.pdf",type.c_str(),hadron.c_str()));   
   
   // make pi over e plot
      TCanvas *cPTheta_pie = new TCanvas("cPTheta_pie","cPTheta_pie",1800,1000);
      TH2F *hPTheta_pie=(TH2F*) hPTheta_pi->Clone();      
      hPTheta_pie->Divide(hPTheta_e);
      hPTheta_pie->SetTitle(Form("%s/%s ratio;#theta (deg);P (GeV/c)",hadron_title.c_str(),e_title.c_str()));             
//       hPTheta_pie->GetXaxis()->SetRangeUser(min_theta,max_theta);
//       hPTheta_pie->GetYaxis()->SetRangeUser(min_p,max_p);
      hPTheta_pie->Draw("text colz");
      cPTheta_pie->SaveAs(Form("%s_ratio_PTheta_%se.pdf",type.c_str(),hadron.c_str()));         
      
      TCanvas *cQ2x_pie = new TCanvas("cQ2x_pie","cQ2x_pie",1800,1000);
      TH2F *hQ2x_pie=(TH2F*) hQ2x_pi->Clone();      
      hQ2x_pie->Divide(hQ2x_e);
      hQ2x_pie->SetTitle(Form("%s/%s ratio;x;Q^{2} (GeV/c)^{2}",hadron_title.c_str(),e_title.c_str()));       
//       hQ2x_pie->GetXaxis()->SetRangeUser(min_theta,max_theta);   
      hQ2x_pie->Draw("text colz");
      cQ2x_pie->SaveAs(Form("%s_ratio_Q2x_%se.pdf",type.c_str(),hadron.c_str()));         
   
      
   // make pi over e plot with rejection    
   gStyle->SetPaintTextFormat("2.2f");    
   
      TCanvas *cPTheta_pie_rej = new TCanvas("cPTheta_pie_rej","cPTheta_pie_rej",1800,1000);   
      TH2F *hPTheta_pie_rej=(TH2F*) hPTheta_pi_rej->Clone();      
      hPTheta_pie_rej->Divide(hPTheta_e);           
      hPTheta_pie_rej->Scale(1e2); //percent

      //change text size
      if (type=="solid_JPsi_LH2" || type=="solid_SIDIS_He3" ) hPTheta_pie_rej->SetMarkerSize(1);
      else if (type=="solid_PVDIS_LD2") hPTheta_pie_rej->SetMarkerSize(2);      

// old way doing rejection
//       TH2F *hPTheta_pie_rej=(TH2F*) hPTheta_pie->Clone();             
// 	  for (int i = 1; i<=20; i++) {	
// 	    for (int j = 1; j<=20; j++) {
// 	      double x=hPTheta_pie_rej->GetXaxis()->GetBinCenter(i);
// 	      double y=hPTheta_pie_rej->GetYaxis()->GetBinCenter(j);	
// 
// 	      if (type=="solid_JPsi_LH2" || type=="solid_SIDIS_He3" ) {		
// 		if (x <15) {//FA
// 		  if (y <4) hPTheta_pie_rej->SetBinContent(i,j,hPTheta_pie_rej->GetBinContent(i,j)*rej_FA_low);
// 		  else  hPTheta_pie_rej->SetBinContent(i,j,hPTheta_pie_rej->GetBinContent(i,j)*rej_FA_high);
// 		}
// 		else hPTheta_pie_rej->SetBinContent(i,j,hPTheta_pie_rej->GetBinContent(i,j)*rej_LA); //LA	      
// 	      }
// 	      else if (type=="solid_PVDIS_LD2") {      
// 		  if (y <4) hPTheta_pie_rej->SetBinContent(i,j,hPTheta_pie_rej->GetBinContent(i,j)*rej_FA_low);
// 		  else  hPTheta_pie_rej->SetBinContent(i,j,hPTheta_pie_rej->GetBinContent(i,j)*rej_FA_high); 
// 	      }
// 	      
// 	    }
// 	  }      
      
      //clean up for PTheta
      	  for (int i = 1; i<=20; i++) {	
	    for (int j = 1; j<=20; j++) {
	      double x=hPTheta_pie_rej->GetXaxis()->GetBinCenter(i);
	      double y=hPTheta_pie_rej->GetYaxis()->GetBinCenter(j);	
	      
	      if (type=="solid_SIDIS_He3" ) {  	      
		if (x==16 && y>6) hPTheta_pie_rej->SetBinContent(i,j,0); // cleanup
		if (x==18 && y>5) hPTheta_pie_rej->SetBinContent(i,j,0); // cleanup	      
	      }
	      else if (type=="solid_PVDIS_LD2") {
		if (x==24 && y>4.5) hPTheta_pie_rej->SetBinContent(i,j,0); // cleanup	      		
	      }
	    }
	  }
	  
      hPTheta_pie_rej->SetTitle(Form("%s/%s ratio (\%) after %s rejection;#theta (deg);P (GeV/c)",hadron_title.c_str(),e_title.c_str(),hadron_title.c_str()));       
      hPTheta_pie_rej->GetXaxis()->SetRangeUser(min_theta,max_theta);   
      hPTheta_pie_rej->GetYaxis()->SetRangeUser(min_p,max_p);                       
      hPTheta_pie_rej->Draw("text colz");

// lines in P vs theta at 11GeV beam
  TF1 *W_2=new TF1("W_2","(0.938*0.938-4+2*0.938*11)/(2*0.938+4*11*sin(x/180.*3.1415926/2)*sin(x/180.*3.1415926/2))",0,180);
   W_2->SetLineColor(kOrange);
  TF1 *W_3=new TF1("W_3","(0.938*0.938-9+2*0.938*11)/(2*0.938+4*11*sin(x/180.*3.1415926/2)*sin(x/180.*3.1415926/2))",0,180);
   W_3->SetLineColor(kOrange);
   W_3->SetLineStyle(2);      
   TF1 *Q2_1=new TF1("Q2_1","1./(4*11*sin(x/180.*3.1415926/2)*sin(x/180.*3.1415926/2))",0,180); 
   Q2_1->SetLineColor(kRed);
   TF1 *Q2_2=new TF1("Q2_2","2./(4*11*sin(x/180.*3.1415926/2)*sin(x/180.*3.1415926/2))",0,180); 
   Q2_2->SetLineColor(kRed);
   Q2_2->SetLineStyle(2);         
   TF1 *x_35=new TF1("x_35","(0.35*2*0.938*11)/(0.35*2*0.938+4*11*sin(x/180.*3.1415926/2)*sin(x/180.*3.1415926/2))",0,180); 
   x_35->SetLineColor(kRed);
   TF1 *x_55=new TF1("x_55","(0.55*2*0.938*11)/(0.55*2*0.938+4*11*sin(x/180.*3.1415926/2)*sin(x/180.*3.1415926/2))",0,180); 
   x_55->SetLineColor(kRed);
   x_55->SetLineStyle(2);
   
  {
   TLegend* leg = new TLegend(0.8, 0.7, 0.9, 0.9);      
   if (type=="solid_SIDIS_He3" ) {
    W_2->Draw("same");
    leg->AddEntry(W_2, "W=2","l");    
    W_3->Draw("same");
    leg->AddEntry(W_3, "W=3","l");    
    Q2_1->Draw("same");
    leg->AddEntry(Q2_1, "Q^{2}=1","l");        
    Q2_2->Draw("same");
    leg->AddEntry(Q2_2, "Q^{2}=2","l");        
    leg->Draw();    
    
    TLine *line_p = new TLine(15,3.5,25,3.5);
    line_p->SetLineWidth(6);
    line_p->SetLineColor(kGreen);        
    line_p->Draw("same");
    TLine *line_theta = new TLine(15,1,15,3.5);    
    line_theta->SetLineWidth(6);
    line_theta->SetLineColor(kGreen);        
    line_theta->Draw("same");    
   }
   else if (type=="solid_PVDIS_LD2") {      
    W_2->Draw("same");
    leg->AddEntry(W_2, "W=2","l");    
    W_3->Draw("same");
    leg->AddEntry(W_3, "W=3","l");    
    x_35->Draw("same");
    leg->AddEntry(x_35, "x=0.35","l");        
    x_55->Draw("same");
    leg->AddEntry(x_55, "x=0.55","l");        
    leg->Draw();    
   }
   }
   
    cPTheta_pie_rej->SaveAs(Form("%s_ratio_PTheta_%se_rej.pdf",type.c_str(),hadron.c_str()));         
   
      TCanvas *cQ2x_pie_rej = new TCanvas("cQ2x_pie_rej","cQ2x_pie_rej",1800,1000);	  
      TH2F *hQ2x_pie_rej=(TH2F*) hQ2x_pi_rej->Clone();      
      hQ2x_pie_rej->Divide(hQ2x_e);           
      hQ2x_pie_rej->Scale(1e2); //percent   
      	  
      hQ2x_pie_rej->SetTitle(Form("%s/%s ratio (\%) after %s rejection;x;Q^{2} (GeV/c)^{2}",hadron_title.c_str(),e_title.c_str(),hadron_title.c_str()));       
//       hQ2x_pie_rej->GetXaxis()->SetRangeUser(min_theta,max_theta);   
//       hQ2x_pie_rej->GetYaxis()->SetRangeUser(min_p,max_p);                       
      hQ2x_pie_rej->Draw("text colz");

    cQ2x_pie_rej->SaveAs(Form("%s_ratio_Q2x_%se_rej.pdf",type.c_str(),hadron.c_str()));               
    
    //physics cut   
//    TCanvas *cPTheta_pie_rej_cut = new TCanvas("pie_rej_cut","pie_rej_cut",1800,1000);  
//     TH2F *hPTheta_pie_rej_cut=(TH2F*) hPTheta_pie_rej->Clone();              
// 
//     for (int i = 1; i<=20; i++) {	
// 	    for (int j = 1; j<=20; j++) {
// 	      double x=hPTheta_pie_rej_cut->GetXaxis()->GetBinCenter(i);
// 	      double y=hPTheta_pie_rej_cut->GetYaxis()->GetBinCenter(j);	
// 	      
// 	      if (type=="solid_SIDIS_He3" ) {  	      
// 		if (x>15 && y<3.5) hPTheta_pie_rej_cut->SetBinContent(i,j,0); // cut LA at 3.5GeV
// 		if (x<9 && y<3.5) hPTheta_pie_rej_cut->SetBinContent(i,j,0); // cut FA Q2=1
// 		if (x<11 && y<2.5) hPTheta_pie_rej_cut->SetBinContent(i,j,0); // cut FA Q2=1
// 		if (x<13 && y<1.5) hPTheta_pie_rej_cut->SetBinContent(i,j,0); // cut FA Q2=1 
// 	      }
// // 	      if (type=="solid_JPsi_LH2" ) {  	      
// // 		if (y<2.5) hPTheta_pie_rej_cut->SetBinContent(i,j,0); // for decay lepton 
// // 	      }
// // 	      if (type=="solid_PVDIS_LD2" ) {  	      
// // 		if (y<1.5) hPTheta_pie_rej_cut->SetBinContent(i,j,0); 
// // 	      }
// 	      
// 	    }
//     }
//    hPTheta_pie_rej_cut->Draw("text colz");
//   
//     cPTheta_pie_rej_cut->SaveAs(Form("%s_rate_%se_rej_cut.pdf",type.c_str(),hadron.c_str()));         
    
}
