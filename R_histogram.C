void R_histogram() {

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // Open file
    TFile *file = TFile::Open("out_QE.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open out.root\n";
        return;
    }

    TTree *Ntuple = (TTree*)file->Get("Ntuple");
    if (!Ntuple) {
        std::cerr << "Error: Ntuple not found\n";
        return;
    }

    Double_t t1, t2, dphi;
    Ntuple->SetBranchAddress("Gamma1_Theta", &t1);
    Ntuple->SetBranchAddress("Gamma2_Theta", &t2);
    Ntuple->SetBranchAddress("dPhi", &dphi);

    // 5 deg binning
    const int Nbins = 180/5;
    TH3F *h_Rmap = new TH3F("h_Rmap",
                            "R(Theta1,Theta2);Theta1 [deg];Theta2 [deg];R",
                            Nbins, 0, 180,
                            Nbins, 0, 180,
                            200, -3.0, 3.0);

    auto computeKa = [&](double th1, double th2) {
        double c1 = cos(th1 * TMath::DegToRad());
        double c2 = cos(th2 * TMath::DegToRad());
        double n1 = pow(1 - c1, 3) + 2;
        double n2 = pow(1 - c2, 3) + 2;
        double d1 = pow(2 - c1, 3);
        double d2 = pow(2 - c2, 3);
        return (n1/d1) * (n2/d2);
    };

    auto computeKb = [&](double th1, double th2) {
        double s1 = sin(th1 * TMath::DegToRad());
        double s2 = sin(th2 * TMath::DegToRad());
        double c1 = cos(th1 * TMath::DegToRad());
        double c2 = cos(th2 * TMath::DegToRad());
        return (s1*s1)/pow(2-c1,2) * (s2*s2)/pow(2-c2,2);
    };

    Long64_t N = Ntuple->GetEntries();
    const double r0 = 2.8179403262e-15;
    const double factor = pow(r0,4) / 16.0;

    for (Long64_t i = 0; i < N; i++) {
        Ntuple->GetEntry(i);

        if (t1 < 0 || t1 > 180 || t2 < 0 || t2 > 180)
            continue;

        double Ka = computeKa(t1, t2);
        double Kb = computeKb(t1, t2);
        double A = Ka * factor;
        double B = Kb * factor;
        double P_perp  = -A + B;
        double P_para  =  A + B;
        double R = P_perp / P_para;

        h_Rmap->Fill(t1, t2, R);
    }

    // Draw the 3D histogram
    TCanvas *c = new TCanvas("c", "R map 3D", 1000, 800);
    c->cd();
    h_Rmap->SetContour(80);
    h_Rmap->Draw("LEGO2");
    c->Update();
    c->SaveAs("theta1_theta2_3d.pdf");
    file->Close();
}
