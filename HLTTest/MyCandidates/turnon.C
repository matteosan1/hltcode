void turnon() {
  
  TChain* t = new TChain("tree");

  t->Add("turnon_D.root");
  //t->Add("root://eoscms//eos/cms/store/caf/user/sani/turnon_RunD/hltD_367_1_aS1.root");
  //t->Add("root://eoscms//eos/cms/store/caf/user/sani/turnon_RunD/hltD_368_1_APx.root");
  //t->Add("root://eoscms//eos/cms/store/caf/user/sani/turnon_RunD/hltD_371_1_teg.root");
  //t->Add("root://eoscms//eos/cms/store/caf/user/sani/turnon_RunD/hltD_372_1_TLp.root");

//TTree* t = (TTree*)f->Get("tree");

  Float_t pt1, pt2, nvtx1, nvtx2, eta1, eta2;
  Float_t ptden, nvtxden, etaden;
  double xbin[29] = {0,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40,45,50,55,60,70,80,100,1000};
  TH1F* hpt1   = new TH1F("hpt1",   "hpt1", 29, xbin);
  TH1F* hpt2   = new TH1F("hpt2",   "hpt2", 29, xbin);
  TH1F* hptden = new TH1F("hptden", "hptden", 29, xbin);
  TH1F* heta1   = new TH1F("eta1",   "eta1",   25, -2.5, 2.5);
  TH1F* heta2   = new TH1F("eta2",   "eta2",   25, -2.5, 2.5);
  TH1F* hetaden = new TH1F("etaden", "etaden", 25, -2.5, 2.5);
  TH1F* hnvtx1   = new TH1F("hnvtx1",   "hnvtx1",   25, 0, 50);
  TH1F* hnvtx2   = new TH1F("hnvtx2",   "hnvtx2",   25, 0, 50);
  TH1F* hnvtxden = new TH1F("hnvtxden", "hnvtxden", 25, 0, 50);

  t->SetBranchAddress("pt1",   &pt1);
  t->SetBranchAddress("pt2",   &pt2);
  t->SetBranchAddress("ptden", &ptden);
  t->SetBranchAddress("eta1",   &eta1);
  t->SetBranchAddress("eta2",   &eta2);
  t->SetBranchAddress("etaden", &etaden);
  t->SetBranchAddress("nvtx1",   &nvtx1);
  t->SetBranchAddress("nvtx2",   &nvtx2);
  t->SetBranchAddress("nvtxden", &nvtxden);
  
  Int_t entries = t->GetEntries();
  
  for(int z=0; z<entries; z++) {
    t->GetEntry(z);
 
    if (fabs(eta1)>1.479)
      hnvtxden->Fill(nvtxden);
    if (fabs(eta1)>1.479 || fabs(eta2)>1.479)
      hptden->Fill(ptden);
    hetaden->Fill(etaden);

    if (pt1 != 999. && fabs(eta1)>1.479)
      hpt1->Fill(pt1);
    if (pt2 != 999.  && fabs(eta2)>1.479)
      hpt2->Fill(pt2);

    if (eta1 != 999.)
      heta1->Fill(eta1);
    if (eta2 != 999.)
      heta2->Fill(eta2);

    if (nvtx1 != 999.  && fabs(eta1)<1.479)
      hnvtx1->Fill(nvtx1);
    if (nvtx2 != 999.  && fabs(eta2)<1.479)
      hnvtx2->Fill(nvtx2);
    
  }

  hpt1->Sumw2();
  hpt2->Sumw2();
  hptden->Sumw2();
  hnvtx1->Sumw2();
  hnvtx2->Sumw2();
  hnvtxden->Sumw2();
  heta1->Sumw2();
  heta2->Sumw2();
  hetaden->Sumw2();
  
  TCanvas* c0 = new TCanvas("c0", "c0");
  c0->Divide(2, 1);
  c0->cd(1);
  hpt1->Divide(hptden);
  hpt1->Draw("PE");
  c0->cd(2);
  hpt2->Divide(hptden);
  hpt2->Draw("PE");


  TCanvas* c1 = new TCanvas("c1", "c1");
  c1->Divide(2, 1);
  c1->cd(1);
  heta1->Divide(hetaden);
  heta1->Draw("PE");
  c1->cd(2);
  heta2->Divide(hetaden);
  heta2->Draw("PE");


  TCanvas* c2 = new TCanvas("c2", "c2");
  c2->Divide(2, 1);
  c2->cd(1);
  hnvtx1->Divide(hnvtxden);
  hnvtx1->Draw("PE");
  c2->cd(2);
  hnvtx2->Divide(hnvtxden);
  hnvtx2->Draw("PE");
}  
  
