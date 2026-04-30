void CircularDalitz_oPs_angle() {

    TChain *hitsTree = new TChain("Hits");
    hitsTree->Add("out_*.root");

    if (hitsTree->GetEntries() == 0) {
        std::cerr << "No entries found!" << std::endl;
        return;
    }

    Int_t EventID;
    Int_t GammaIndex;
    Float_t HitX, HitY, HitZ;

    hitsTree->SetBranchAddress("EventID", &EventID);
    hitsTree->SetBranchAddress("GammaIndex", &GammaIndex);
    hitsTree->SetBranchAddress("HitX", &HitX);
    hitsTree->SetBranchAddress("HitY", &HitY);
    hitsTree->SetBranchAddress("HitZ", &HitZ);

    Long64_t nentries = hitsTree->GetEntries();

    std::map<int, std::map<int, TVector3>> events;

    // getting direction vectors of interaction
    for (Long64_t i = 0; i < nentries; ++i) {

        hitsTree->GetEntry(i);

        if (GammaIndex >= 0 && GammaIndex < 3) {
            //std::cout << HitX << " " << HitY << " " << HitZ << std::endl;
            // only storing first interaction per gamma
            if (events[EventID].count(GammaIndex) == 0) {

                TVector3 dir(HitX, HitY, HitZ);

                events[EventID][GammaIndex] = dir.Unit();
                
            }
        }
    }

    TH2F *hDalitz = new TH2F(
        "hDalitz",
        "Dalitz plot for o-Ps #rightarrow 3#gamma;#theta_{12} (deg);#theta_{23} (deg)",
        300, 0, 180,
        300, 0, 180);

    int filledEvents = 0;

    for (auto &evt : events) {

        if (evt.second.size() == 3) {

            TVector3 d0 = evt.second[0];
            TVector3 d1 = evt.second[1];
            TVector3 d2 = evt.second[2];
            
            // converting to degrees
            double theta12 = d0.Angle(d1) * 180.0 / TMath::Pi();
            double theta23 = d1.Angle(d2) * 180.0 / TMath::Pi();

            hDalitz->Fill(theta12, theta23);
            filledEvents++;
        }
    }
    //for (auto &evt : events) {
      //  if (evt.second.size() == 3) {
        //    TVector3 d0 = evt.second[0];
          //  TVector3 d1 = evt.second[1];
            //TVector3 d2 = evt.second[2];
            //double theta12 = d0.Angle(d1) * 180.0 / TMath::Pi();
            //double theta23 = d1.Angle(d2) * 180.0 / TMath::Pi();
            //std::cout << "theta12 = " << theta12 << ", theta23 = " << theta23 << std::endl;
            //break; // just first event
        // }
    // }
    //std::cout << "Filled events: " << filledEvents << std::endl;
    //std::cout << "Histogram entries: " << hDalitz->GetEntries() << std::endl;

    TCanvas *c = new TCanvas("c", "Angular Dalitz", 800, 800);
    
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.12);
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.12);

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kRainBow);
    hDalitz->Draw("COLZ");


    c->SaveAs("dalitz_oPs_angular.pdf");
}
