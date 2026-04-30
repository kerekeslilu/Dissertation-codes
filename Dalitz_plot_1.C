void CircularDalitz_oPs() {

    TChain *hitsTree = new TChain("Hits");
    hitsTree->Add("out_*.root");

    if (hitsTree->GetEntries() == 0) {
        std::cerr << "No entries found!" << std::endl;
        return;
    }

    Int_t EventID;
    Int_t GammaIndex;
    Float_t VertexE;

    hitsTree->SetBranchAddress("EventID", &EventID);
    hitsTree->SetBranchAddress("GammaIndex", &GammaIndex);
    hitsTree->SetBranchAddress("VertexE", &VertexE);

    std::set<int> countedGammas;

    TH2F *hDalitz = new TH2F(
        "hDalitz",
        "Energy Correlations in o-Ps #rightarrow 3#gamma (Circular Dalitz Plot);\
        X = (E_{2}/E_{tot} - E_{1}/E_{tot}) / #sqrt{2};\
        Y = (2E_{3}/E_{tot} - E_{1}/E_{tot} - E_{2}/E_{tot}) / #sqrt{6}",
        300, -0.5, 0.5,
        300, -0.5, 0.5
    );

    hDalitz->GetZaxis()->SetTitle("Counts");

    Int_t currentEvent = -1;
    double E[3] = {0.0, 0.0, 0.0};

    Long64_t nentries = hitsTree->GetEntries();

    for (Long64_t i = 0; i < nentries; ++i) {

        hitsTree->GetEntry(i);

        if (EventID != currentEvent && currentEvent != -1) {

            if (countedGammas.size() == 3) {

                double Etot = E[0] + E[1] + E[2];

                if (Etot > 0) {

                    //std::sort(E, E + 3);

                    double eps1 = E[0] / Etot;
                    double eps2 = E[1] / Etot;
                    double eps3 = E[2] / Etot;

                    double X = (eps2 - eps1) / std::sqrt(2.0);
                    double Y = (2.0 * eps3 - eps1 - eps2) / std::sqrt(6.0);

                    hDalitz->Fill(X, Y);
                }
            }

            // Reset for next event
            E[0] = E[1] = E[2] = 0.0;
            countedGammas.clear();
        }

        currentEvent = EventID;

        if (GammaIndex >= 0 && GammaIndex < 3) {
            if (countedGammas.find(GammaIndex) == countedGammas.end()) {
                E[GammaIndex] = VertexE;
                countedGammas.insert(GammaIndex);
            }
        }
    }

    // Process final event
    if (countedGammas.size() == 3) {

        double Etot = E[0] + E[1] + E[2];

        if (Etot > 0) {

            std::sort(E, E + 3);

            double eps1 = E[0] / Etot;
            double eps2 = E[1] / Etot;
            double eps3 = E[2] / Etot;

            double X = (eps2 - eps1) / std::sqrt(2.0);
            double Y = (2.0 * eps3 - eps1 - eps2) / std::sqrt(6.0);

            hDalitz->Fill(X, Y);
        }
    }

    /*TCanvas *c = new TCanvas("c", "Circular Dalitz Plot", 1600, 1600);
    c->SetLeftMargin(0.12);
    c->SetRightMargin(0.15);
    c->SetBottomMargin(0.12);
    c->SetTopMargin(0.08);

    gStyle->SetOptStat(0);
    gStyle->SetPalette(kViridis);
    c->SetFixedAspectRatio(1);
    hDalitz->Draw("COLZ");
    c->Update();
    c->SaveAs("dalitz_oPs_circular_final_1.png");*/
    TCanvas *c = new TCanvas("c", "Dalitz Plot for 3-Gamma Energy distribution", 1600, 1600);
    hDalitz->Draw("COLZ");
    c->SetFixedAspectRatio();
    c->Update();
    c->SaveAs("dalitz_oPs_final.png");
}
