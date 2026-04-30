void R_histogram_fit()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    TFile *file = TFile::Open("out_QE1.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open out.root\n";
        return;
    }

    TTree *Ntuple = (TTree*)file->Get("Ntuple");
    if (!Ntuple) { std::cerr << "Error: Ntuple not found!\n"; return; }

    Double_t Gamma1_Theta, Gamma2_Theta, dPhi;
    Ntuple->SetBranchAddress("Gamma1_Theta", &Gamma1_Theta);
    Ntuple->SetBranchAddress("Gamma2_Theta", &Gamma2_Theta);
    Ntuple->SetBranchAddress("dPhi", &dPhi);
    
    TH1F *h_80_84 = new TH1F("dPhi", "dPhi for 80 < Gamma1_Theta, Gamma2_Theta < 84", 50, -200, 200);
    
    const double step = 4.0;
    const int Nbins = 180 / step;
    
    TH2F *h_R = new TH2F("h_R",
                           "R(Theta1,Theta2);Theta1 [deg];Theta2 [deg];R",
                           Nbins, 0, 180,
                           Nbins, 0, 180);

    std::vector<std::vector<TH1F*>> H(Nbins,
         std::vector<TH1F*>(Nbins, nullptr));

    for (int i = 0; i < Nbins; i++) {
        for (int j = 0; j < Nbins; j++) {

            TString name = Form("h_%d_%d", i, j);
            TString title = Form("dPhi for %d-%d deg, %d-%d deg",
                                  i*4, (i+1)*4, j*4, (j+1)*4);
            H[i][j] = new TH1F(name, title, 50, -200, 200);
        }
    }

    Long64_t nEntries = Ntuple->GetEntries();
    for (Long64_t k = 0; k < nEntries; k++) {
        Ntuple->GetEntry(k);

        if (Gamma1_Theta < 0 || Gamma1_Theta > 180 || Gamma2_Theta < 0 || Gamma2_Theta > 180)
            continue;

        int i = int(Gamma1_Theta / step);
        int j = int(Gamma2_Theta / step);

        H[i][j]->Fill(dPhi);
        
        if (Gamma1_Theta > 80 && Gamma2_Theta < 84 &&
            Gamma2_Theta > 80 && Gamma1_Theta < 84) {
            h_80_84->Fill(dPhi);
        }
    }

    TF1 *fit80 = new TF1("fit80", "[0]*cos(2*TMath::DegToRad()*x) + [1]", -180, 180);
    fit80->SetLineColor(kRed);
    fit80->SetLineWidth(2);
    h_80_84->Fit(fit80, "R");
    TCanvas *c80 = new TCanvas("c80","DeltaPhi fit: 80–84 deg", 800, 600);
    h_80_84->SetTitle("80–84 deg, 80–84 deg;#Delta#phi [deg];Counts");
    h_80_84->Draw("E");
    fit80->Draw("SAME");
    double A = fit80->GetParameter(0);
    double B = fit80->GetParameter(1);
    double R = (B - A) / (B + A);
    TLegend *leg1 = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg1->AddEntry(h_80_84, Form("A = %.3f, B = %.3f, R = %.3f", A, B, R), "l");
    leg1->Draw();
    std::cout << "Enhancement R (80–84 deg) = " << R << std::endl;
    //c80->SaveAs("dPhi_80_84_test.pdf");
    

    TF1 *fitFcn = new TF1("fitFcn", "[0]*cos(2*TMath::DegToRad()*x) + [1]", -180, 180);

    for (int i = 0; i < Nbins; i++) {
        for (int j = 0; j < Nbins; j++) {

            TH1F *h = H[i][j];

            // if (h->GetEntries() < 100)
                // continue;

            h->Fit(fitFcn, "Q");

            double A = fitFcn->GetParameter(0);
            double B = fitFcn->GetParameter(1);

            double R = (B - A) / (B + A);
            std::cout << "A, B, R = "
                      << A << ", "
                      << B << ", "
                      << R << std::endl;


            double theta1_center = (i + 0.5)*step;
            double theta2_center = (j + 0.5)*step;
            std::cout << "theta1 centre = " << theta1_center << std::endl;
            std::cout << "theta2 centre = " << theta2_center << std::endl;
            h_R->SetBinContent(i+1, j+1, R);
        }
    }

 
    TCanvas *c = new TCanvas("c", "R(Theta1,Theta2)", 1000, 800);
    c->SetRightMargin(0.15);

    h_R->GetZaxis()->SetTitle("R");
    h_R->GetZaxis()->SetTitleOffset(1.3);
    h_R->GetXaxis()->SetTitleOffset(1.8);
    h_R->GetYaxis()->SetTitleOffset(1.8);
    h_R->Draw("LEGO2");   // or "LEGO2", "SURF3"

    c->Update();
    c->SaveAs("R_3D_surface.pdf");

    file->Close();
}

