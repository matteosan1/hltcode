#include "TFile.h"
#include "TTree.h"
#include "TProfile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TEfficiency.h"

#include <iostream>

float deltaPhi(Float_t p1, Float_t p2) {

  float dp = p1-p2; 
  if (dp > float(TMath::Pi())) 
    dp-=(2*float(TMath::Pi()));  
  return dp;
}

float deltaR(Float_t e1, Float_t p1, Float_t e2, Float_t p2) {

  float dp = fabs(p1-p2); 
  if (dp > float(TMath::Pi())) 
    dp-=(2*float(TMath::Pi()));  
  return (e1-e2)*(e1-e2) + dp*dp;
}

void plotMacro(const char* filename, bool barrel=true) {
  TFile* file = new TFile(filename);
  TTree* t = (TTree*)file->Get("tree");

  Int_t n, el_n;
  Float_t el_eta[10], el_phi[10], el_e[10];
  Float_t eta[10], e[10], phi[10];
  Float_t sieie[10], el_sieie[10];
  Float_t hoe[10], el_hoe[10];
  Float_t hcal[10], el_hcal[10];
  Float_t ecal[10], el_ecal[10];
  Float_t tkiso[10], el_tkiso[10];
  Float_t dphi[10], el_dphi[10];
  Float_t deta[10], el_deta[10];

  t->SetBranchAddress("n", &n);
  t->SetBranchAddress("eta", &eta);
  t->SetBranchAddress("e", &e);
  t->SetBranchAddress("phi", &phi);
  t->SetBranchAddress("sieie", &sieie);
  t->SetBranchAddress("hoe", &hoe);
  t->SetBranchAddress("hcal", &hcal);
  t->SetBranchAddress("ecal", &ecal);
  t->SetBranchAddress("tkiso", &tkiso);
  t->SetBranchAddress("deta", &deta);
  t->SetBranchAddress("dphi", &dphi);

  t->SetBranchAddress("el_n",     &el_n);
  t->SetBranchAddress("el_eta",   &el_eta);
  t->SetBranchAddress("el_e",    &el_e);
  t->SetBranchAddress("el_phi",   &el_phi);
  t->SetBranchAddress("el_sieie", &el_sieie);
  t->SetBranchAddress("el_hcal",  &el_hcal);
  t->SetBranchAddress("el_hoe",   &el_hoe);
  t->SetBranchAddress("el_ecal",   &el_ecal);
  t->SetBranchAddress("el_tkiso",   &el_tkiso);
  t->SetBranchAddress("el_deta", &el_deta);
  t->SetBranchAddress("el_dphi", &el_dphi);

  Int_t entries = t->GetEntries();

  TH1F* heR = new TH1F("heR", "heR", 80, 20., 100.);
  TH1F* heO = new TH1F("heO", "heO", 80, 20., 100.);
  TH1F* heD = new TH1F("heD", "heD", 100, -.1, .1);
  TH1F* hsieieR = new TH1F("hsieieR", "hsieieR", 100, 0.0, 0.035);
  TH1F* hsieieO = new TH1F("hsieieO", "hsieieO", 100, 0.0, 0.035);
  TH1F* hsieieD = new TH1F("hsieieD", "hsieieD", 100, -0.1, 0.1);
  TH1F* hhoeR = new TH1F("hhoeR", "hhoeR", 100, 0., 10.);
  TH1F* hhoeO = new TH1F("hhoeO", "hhoeO", 100, 0., 10.);
  TH1F* hhoeD = new TH1F("hhoeD", "hhoeD", 100, -0.1, 0.1);
  TH1F* hhcalR = new TH1F("hhcalR", "hhcalR", 100, 0., 5.);
  TH1F* hhcalO = new TH1F("hhcalO", "hhcalO", 100, 0., 5.);
  TH1F* hhcalD = new TH1F("hhcalD", "hhcalD", 100, -0.1, 0.1);
  TH1F* hecalR = new TH1F("hecalR", "hecalR", 100, 0., 10.);
  TH1F* hecalO = new TH1F("hecalO", "hecalO", 100, 0., 10.);
  TH1F* hecalD = new TH1F("hecalD", "hecalD", 100, -.5, .5);
  TH1F* htkisoR = new TH1F("htkisoR", "htkisoR", 100, 0., 5.);
  TH1F* htkisoO = new TH1F("htkisoO", "htkisoO", 100, 0., 5.);
  TH1F* htkisoD = new TH1F("htkisoD", "htkisoD", 100, -.1, .1);
  TH1F* hdetaR = new TH1F("hdetaR", "hdetaR", 100, 0., 0.01);
  TH1F* hdetaO = new TH1F("hdetaO", "hdetaO", 100, 0., 0.01);
  TH1F* hdetaD = new TH1F("hdetaD", "hdetaD", 100, -0.005, 0.005);
  TH1F* hdphiR = new TH1F("hdphiR", "hdphiR", 100, 0., .05);
  TH1F* hdphiO = new TH1F("hdphiO", "hdphiO", 100, 0., .05);
  TH1F* hdphiD = new TH1F("hdphiD", "hdphiD", 100, -0.1, 0.1);
  
  for (int z=0; z<entries; ++z) {
    t->GetEntry(z);

    for (int j=0; j<n; j++) {
      if (e[j]*2*atan(exp(-eta[j])) < 20.)
	continue;

      if ((fabs(eta[j]) > 1.479 && barrel) || (fabs(eta[j])<1.479 && !barrel))
      	  continue;
      Int_t index = -1;
      Float_t minDr = 0.1;
      for (int i=0; i<el_n; ++i) {
	Float_t dR = deltaR(eta[j], phi[j], el_eta[i], el_phi[i]);
	if (dR < minDr) {
	  minDr = dR;
	  index = i;
	}
      }
      
      if (index != -1) {
	float thetaR = 2*atan(exp(-el_eta[index]));
	float thetaO = 2*atan(exp(-eta[j]));
	float etR = el_e[index]*sin(thetaR);
	float etO = e[j]*sin(thetaO);
	heR->Fill(etR);
	heO->Fill(etO);
	heD->Fill((etO-etR)/etR);
	//heR->Fill(el_e[index]);
	//heO->Fill(e[j]);
	//heD->Fill(e[j]-el_e[index]);
	

	if (sieie[j] < 1000.) { 
	  hsieieR->Fill(el_sieie[index]);
	  hsieieO->Fill(sieie[j]);
	  hsieieD->Fill((sieie[j]-el_sieie[index])/el_sieie[index]);
	}
	if (hoe[j] < 1000.) {
	  //hhoeR->Fill(el_hoe[index]/el_e[index]);
	  //hhoeO->Fill(hoe[j]/e[j]);
	  hhoeR->Fill(el_hoe[index]);
	  hhoeO->Fill(hoe[j]);
	  if (el_hoe[index] == 0 && hoe[j] == 0) 
	    hhoeD->Fill(0);
	  else if (el_hoe[index] == 0 && hoe[j] != 0) 
	    hhoeD->Fill(1);
	  else 
	    //hhoeD->Fill((hoe[j]/e[j]-el_hoe[index]/el_e[index])/el_hoe[index]*el_e[index]);
	    hhoeD->Fill((hoe[j]-el_hoe[index])/el_hoe[index]);
	}
	if (hcal[j] < 1000.) {
	  //hhcalR->Fill(el_hcal[index]/etR);
	  //hhcalO->Fill(hcal[j]/etO);
	  hhcalR->Fill(el_hcal[index]);
	  hhcalO->Fill(hcal[j]);
	  if (el_hcal[index] == 0 && hcal[j] == 0)
	    hhcalD->Fill(0);
	  else if(el_hcal[index] == 0 && hcal[j] != 0)
	    hhcalD->Fill(1);
	  else
	    hhcalD->Fill((hcal[j]-el_hcal[index])/el_hcal[index]);
	}
	if (ecal[j] < 1000.) {
	  //hecalR->Fill(el_ecal[index]/etR);
	  //hecalO->Fill(ecal[j]/etO);
	  hecalR->Fill(el_ecal[index]);
	  hecalO->Fill(ecal[j]);
	  
	  if (el_ecal[index] == 0 && ecal[j] == 0)
	    hecalD->Fill(0);
	  else if (el_ecal[index] == 0 && ecal[j] != 0)
	    hecalD->Fill(1);
	  else
	    hecalD->Fill((ecal[j]-el_ecal[index])/el_ecal[index]);
	  //hecalD->Fill((ecal[j]-el_ecal[index]));
	}
	if (tkiso[j] < 1000.) {
	  //htkisoR->Fill(el_tkiso[index]/etR);
	  //htkisoO->Fill(tkiso[j]/etO);
	  htkisoR->Fill(el_tkiso[index]);
	  htkisoO->Fill(tkiso[j]);
	  if (el_tkiso[index] == 0 && tkiso[j] == 0)
	    htkisoD->Fill(0);
	  else if (el_tkiso[index] == 0 && tkiso[j] != 0)
	    htkisoD->Fill(1);
	  else
	    htkisoD->Fill((tkiso[j]-el_tkiso[index])/el_tkiso[index]);
	}

	if (deta[j] < 1000.) {
	  //htkisoR->Fill(el_tkiso[index]/etR);
	  //htkisoO->Fill(tkiso[j]/etO);
	  hdetaR->Fill(fabs(el_deta[index]));
	  hdetaO->Fill(deta[j]);
	  //if (el_deta[index] == 0 && deta[j] == 0)
	  //  hdetaD->Fill(0);
	  //else if (el_deta[index] == 0 && deta[j] != 0)
	  //  hdetaD->Fill(1);
	  //else
	  //hdetaD->Fill((deta[j]/etO-el_deta[index]/etR)/el_deta[index]*etR);
	  hdetaD->Fill((deta[j]-fabs(el_deta[index])));///fabs(el_deta[index]));
	}

	if (dphi[j] < 1000.) {
	  //hdphiR->Fill(el_dphi[index]/etR);
	  //hdphiO->Fill(dphi[j]/etO);
	  hdphiR->Fill(fabs(el_dphi[index]));
	  hdphiO->Fill(dphi[j]);
	  //if (el_dphi[index] == 0 && dphi[j] == 0)
	  //  hdphiD->Fill(0);
	  //else if (el_dphi[index] == 0 && dphi[j] != 0)
	  //  hdphiD->Fill(1);
	  //else
	  //hdphiD->Fill((dphi[j]/etO-el_dphi[index]/etR)/el_dphi[index]*etR);
	  hdphiD->Fill((dphi[j]-el_dphi[index]));///el_dphi[index]);
	}
      }
    }
  }

  TCanvas* c0 = new TCanvas("c0", "c0");
  c0->Divide(2,1);
  c0->cd(1);
  heR->SetMarkerStyle(20);
  heO->SetFillColor(kRed);
  heR->Draw("PE");
  heO->Draw("SAME");
  heR->Draw("PESAME");
  c0->cd(2);
  heD->SetBinContent(heD->GetNbinsX(), heD->GetBinContent(heD->GetNbinsX())+ heD->GetBinContent(heD->GetNbinsX()+1));
  heD->SetBinContent(1, heD->GetBinContent(0)+ heD->GetBinContent(1));
  heD->Draw();
  heD->GetXaxis()->SetNdivisions(505);

  TCanvas* c1 = new TCanvas("c1", "c1");
  c1->Divide(2,1);
  c1->cd(1);
  hsieieR->SetMarkerStyle(20);
  hsieieO->SetFillColor(kRed);
  hsieieR->Draw("PE");
  hsieieO->Draw("SAME");
  hsieieR->Draw("PESAME");
  hsieieR->GetXaxis()->SetNdivisions(505);
  c1->cd(2);
  hsieieD->SetBinContent(hsieieD->GetNbinsX(), hsieieD->GetBinContent(hsieieD->GetNbinsX())+ hsieieD->GetBinContent(hsieieD->GetNbinsX()+1));
  hsieieD->SetBinContent(1, hsieieD->GetBinContent(0)+ hsieieD->GetBinContent(1));
  hsieieD->Draw();
  hsieieD->GetXaxis()->SetNdivisions(505);

  TCanvas* c2 = new TCanvas("c2", "c2");
  c2->Divide(2,1);
  c2->cd(1);
  hhoeR->SetMarkerStyle(20);
  hhoeO->SetFillColor(kRed);
  hhoeR->Draw("PE");
  hhoeO->Draw("SAME");
  hhoeR->Draw("PESAME");
  c2->GetPad(1)->SetLogy(1);
  c2->cd(2);
  hhoeD->SetBinContent(hhoeD->GetNbinsX(), hhoeD->GetBinContent(hhoeD->GetNbinsX())+ hhoeD->GetBinContent(hhoeD->GetNbinsX()+1));
  hhoeD->SetBinContent(1, hhoeD->GetBinContent(0)+ hhoeD->GetBinContent(1));
  hhoeD->Draw();
  c2->GetPad(2)->SetLogy(1);
  hhoeD->GetXaxis()->SetNdivisions(505);

  TCanvas* c3 = new TCanvas("c3", "c3");
  c3->Divide(2,1);
  c3->cd(1);
  hhcalR->SetMarkerStyle(20);
  hhcalO->SetFillColor(kRed);
  hhcalR->Draw("PE");
  hhcalO->Draw("SAME");
  hhcalR->Draw("PESAME");
  c3->GetPad(1)->SetLogy(1);
  c3->cd(2);
  hhcalD->SetBinContent(hhcalD->GetNbinsX(), hhcalD->GetBinContent(hhcalD->GetNbinsX())+ hhcalD->GetBinContent(hhcalD->GetNbinsX()+1));
  hhcalD->SetBinContent(1, hhcalD->GetBinContent(0)+ hhcalD->GetBinContent(1));
  hhcalD->Draw();
  c3->GetPad(2)->SetLogy(1);
  hhcalD->GetXaxis()->SetNdivisions(505);

  TCanvas* c4 = new TCanvas("c4", "c4");
  c4->Divide(2,1);
  c4->cd(1);
  hecalR->SetMarkerStyle(20);
  hecalO->SetFillColor(kRed);
  hecalR->Draw("PE");
  hecalO->Draw("SAME");
  hecalR->Draw("PESAME");
  c4->cd(2);
  hecalD->SetBinContent(hecalD->GetNbinsX(), hecalD->GetBinContent(hecalD->GetNbinsX())+ hecalD->GetBinContent(hecalD->GetNbinsX()+1));
  hecalD->SetBinContent(1, hecalD->GetBinContent(0)+ hecalD->GetBinContent(1));
  hecalD->Draw();

  TCanvas* c5 = new TCanvas("c5", "c5");
  c5->Divide(2,1);
  c5->cd(1);
  htkisoR->SetMarkerStyle(20);
  htkisoO->SetFillColor(kRed);
  htkisoR->Draw("PE");
  htkisoO->Draw("SAME");
  htkisoR->Draw("PESAME");
  c5->GetPad(1)->SetLogy(1);
  c5->cd(2);
  htkisoD->SetBinContent(htkisoD->GetNbinsX(), htkisoD->GetBinContent(htkisoD->GetNbinsX())+ htkisoD->GetBinContent(htkisoD->GetNbinsX()+1));
  htkisoD->SetBinContent(1, htkisoD->GetBinContent(0)+ htkisoD->GetBinContent(1));
  htkisoD->Draw();
  c5->GetPad(2)->SetLogy(1);
  htkisoD->GetXaxis()->SetNdivisions(505);

  TCanvas* c6 = new TCanvas("c6", "c6");
  c6->Divide(2,1);
  c6->cd(1);
  hdetaR->SetMarkerStyle(20);
  hdetaO->SetFillColor(kRed);
  hdetaR->Draw("PE");
  hdetaO->Draw("SAME");
  hdetaR->Draw("PESAME");
  hdetaR->GetXaxis()->SetNdivisions(505);
  c6->cd(2);
  hdetaD->SetBinContent(hdetaD->GetNbinsX(), hdetaD->GetBinContent(hdetaD->GetNbinsX())+ hdetaD->GetBinContent(hdetaD->GetNbinsX()+1));
  hdetaD->SetBinContent(1, hdetaD->GetBinContent(0)+ hdetaD->GetBinContent(1));
  hdetaD->Draw();
  hdetaD->GetXaxis()->SetNdivisions(505);

  TCanvas* c7 = new TCanvas("c7", "c7");
  c7->Divide(2,1);
  c7->cd(1);
  hdphiR->SetMarkerStyle(20);
  hdphiO->SetFillColor(kRed);
  hdphiR->Draw("PE");
  hdphiO->Draw("SAME");  
  hdphiR->Draw("PESAME");
  hdphiR->GetXaxis()->SetNdivisions(505);
  c7->cd(2);
  hdphiD->SetBinContent(hdphiD->GetNbinsX(), hdphiD->GetBinContent(hdphiD->GetNbinsX())+ hdphiD->GetBinContent(hdphiD->GetNbinsX()+1));
  hdphiD->SetBinContent(1, hdphiD->GetBinContent(0)+ hdphiD->GetBinContent(1));
  hdphiD->Draw();
  hdphiD->GetXaxis()->SetNdivisions(505);
}
