void some_theta() {
    gStyle->SetOptStat(0);  // hide the mean/std box
    gStyle->SetOptFit(0);
    // Open the ROOT file change name if need be
    TFile *file = TFile::Open("out_QE1.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Could not open file out.root!" << std::endl;
        return;
    }
    file->ls();

    // Access the tree named qepet
    TTree *Ntuple = (TTree*)file->Get("Ntuple");
    if (!Ntuple) {
        std::cerr << "Error: Ntuple does not exist in the file!" << std::endl;
        return;
    }

    // Define branches for Gamma1_Theta, Gamma2_Theta, and dPhi
    Double_t Gamma1_Theta, Gamma2_Theta, dPhi;
    Ntuple->SetBranchAddress("Gamma1_Theta", &Gamma1_Theta);
    Ntuple->SetBranchAddress("Gamma2_Theta", &Gamma2_Theta);
    Ntuple->SetBranchAddress("dPhi", &dPhi);

    TH1F *hist_1 = new TH1F("dPhi_1", "dPhi for 0 < Gamma1_Theta, Gamma2_Theta < 180", 150, -200, 200);
    TH1F *hist_2 = new TH1F("dPhi_2", "dPhi for 60 < Gamma1_Theta, Gamma2_Theta < 100", 70, -200, 200);
    TH1F *hist_3 = new TH1F("dPhi_3", "dPhi for 80 < Gamma1_Theta, Gamma2_Theta < 84", 50, -200, 200);


    Double_t Theta0_min = 0.0; 
    Double_t Theta0_max = 180.0;
    Double_t Theta1_min = 60;
    Double_t Theta1_max = 100.0;
    Double_t Theta2_min = 80;
    Double_t Theta2_max = 84; 

    
    Long64_t nEntries = Ntuple->GetEntries();
    for (Long64_t j = 0; j < nEntries; j++) {
        Ntuple->GetEntry(j);
    
   //  std::cout << "Entry " << j << ": Gamma1_Theta = " << Gamma1_Theta << ", Gamma2_Theta = " << Gamma2_Theta << std::endl;
    
        if (Gamma1_Theta > Theta0_min && Gamma1_Theta < Theta0_max &&
            Gamma2_Theta > Theta0_min && Gamma2_Theta < Theta0_max) {
            hist_1->Fill(dPhi);
        }
        
        if (Gamma1_Theta > Theta1_min && Gamma1_Theta < Theta1_max &&
            Gamma2_Theta > Theta1_min && Gamma2_Theta < Theta1_max) {
            hist_2->Fill(dPhi);                                       
        }

        if (Gamma1_Theta > Theta2_min && Gamma1_Theta < Theta2_max &&
            Gamma2_Theta > Theta2_min && Gamma2_Theta < Theta2_max) {
            hist_3->Fill(dPhi);
        }
    }

    TF1 *fitFunc1 = new TF1("fitFunc1", "[0]*cos(2*TMath::DegToRad()*x) + [1]", -180, 180);
    TF1 *fitFunc2 = new TF1("fitFunc2", "[0]*cos(2*TMath::DegToRad()*x) + [1]", -180, 180);
    TF1 *fitFunc3 = new TF1("fitFunc3", "[0]*cos(2*TMath::DegToRad()*x) + [1]", -180, 180);

    fitFunc1->SetLineColor(kRed);
    fitFunc2->SetLineColor(kBlue);
    fitFunc3->SetLineColor(kGreen+2);
    
    TCanvas *c1 = new TCanvas("c1", "0 < Theta < 180", 1200, 800);
    hist_1->SetLineColor(kRed);
    hist_1->SetMinimum(74000);
    hist_1->Draw("E1 HIST");
    hist_1->Fit(fitFunc1, "R");
    fitFunc1->SetLineColor(kRed);
    fitFunc1->SetLineWidth(3);
    fitFunc1->Draw("same");
    double A1 = fitFunc1->GetParameter(0);
    double B1 = fitFunc1->GetParameter(1);
    double p01 = fitFunc1->GetParError(0);
    double p11 = fitFunc1->GetParError(1);
    double R1 = (-A1 + B1) / (B1 + A1); 
    double dR1 = sqrt(pow((-2*B1) / pow(B1 + A1, 2) * p01, 2) + pow((2*A1) / pow(B1 + A1, 2) * p11, 2));
    TLegend *leg1 = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg1->AddEntry(hist_1, Form("#splitline{A=%.3f,  B=%.3f,}{ R=%.3f #pm %.3f}", A1, B1, R1, dR1), "l");
    leg1->Draw();
    c1->SaveAs("fit_dphi_0_180_QE2.pdf");

    TCanvas *c2 = new TCanvas("c2", "60 < Theta < 100", 1200, 800);
    hist_2->SetLineColor(kBlue);
    hist_2->SetMinimum(4400);
    hist_2->Draw("E1 HIST");
    hist_2->Fit(fitFunc2, "R");
    fitFunc2->SetLineColor(kBlue);
    fitFunc2->SetLineWidth(2);
    fitFunc2->Draw("same");
    double A2 = fitFunc2->GetParameter(0);
    double B2 = fitFunc2->GetParameter(1);
    double p02 = fitFunc2->GetParError(0);               
    double p12 = fitFunc2->GetParError(1);               
    double dR2 = sqrt(pow((-2*B2) / pow(B2 + A2, 2) * p02, 2) + pow((2*A2) / pow(B2 + A2, 2) * p12, 2));
    double R2 = (-A2 + B2)/(B2 + A2);
    TLegend *leg2 = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg2->AddEntry(hist_2, Form("#splitline{A=%.3f,  B=%.3f,}{ R=%.3f #pm %.3f}", A2, B2, R2, dR2), "l");
    leg2->Draw();
    c2->SaveAs("fit_dphi_60_100_70.pdf");
    
    TCanvas *c3 = new TCanvas("c3", "80 < Theta < 84", 1200, 800);
    hist_3->SetLineColor(kGreen+2);
    //hist_3->SetMaximum(50);
    hist_3->Draw("E1 HIST");
    hist_3->Fit(fitFunc3, "R");
    fitFunc3->SetLineColor(kGreen +2);
    fitFunc3->SetLineWidth(2);
    fitFunc3->Draw("same");
    double A3 = fitFunc3->GetParameter(0);
    double B3 = fitFunc3->GetParameter(1);
    double p03 = fitFunc3->GetParError(0);
    double p13 = fitFunc3->GetParError(1);
    double R3 = (- A3 + B3) / (B3 + A3);
    double dR3 = sqrt(pow((-2*B3) / pow(B3 + A3, 2) * p03, 2) + pow((2*A3) / pow(B3 + A3, 2) * p13, 2));
    TLegend *leg3 = new TLegend(0.65, 0.75, 0.88, 0.88);
    leg3->AddEntry(hist_3, Form("#splitline{A=%.3f,  B=%.3f,}{ R=%.3f #pm %.3f}", A3, B3, R3, dR3), "l");
    leg3->Draw();
    c3->SaveAs("fit_dphi_80_84_QE2.pdf");

    file->Close();
  
}
