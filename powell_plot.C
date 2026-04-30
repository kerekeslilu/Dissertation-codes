void powell_plot() {

    const double Emax = 0.5110;   // keV = m_e c^2
    const int Nbins = 100;

    TChain *hitsTree = new TChain("Hits");
    hitsTree->Add("out_*.root");

    if (hitsTree->GetEntries() == 0) {
        std::cerr << "No entries found!" << std::endl;
        return;
    }

    Int_t EventID;
    Int_t GammaIndex;
    Float_t VertexE;
    Char_t CreationProcess[50];

    hitsTree->SetBranchAddress("CreationProcess", CreationProcess);
    hitsTree->SetBranchAddress("EventID", &EventID);
    hitsTree->SetBranchAddress("GammaIndex", &GammaIndex);
    hitsTree->SetBranchAddress("VertexE", &VertexE);

    TH1F *hData = new TH1F(
        "hData",
        "Ore-Powell Energy Spectrum;E_{#gamma}/m_{e}c^{2};(1/N) dN/dx",
        Nbins, 0.0, 1.0);

    Long64_t nentries = hitsTree->GetEntries();

    // Avoid double counting photons
    std::set<std::pair<int,int>> counted;

    for (Long64_t i = 0; i < nentries; ++i) {

        hitsTree->GetEntry(i);

        if (std::string(CreationProcess) != "annihil")
            continue;

        std::pair<int,int> key = {EventID, GammaIndex};

        if (counted.find(key) != counted.end())
            continue;

        counted.insert(key);

        double x = VertexE / Emax;

        if (x > 0.0 && x < 1.0)
            hData->Fill(x);
    }

    if (hData->Integral("width") > 0)
        hData->Scale(1.0 / hData->Integral("width"));


    TF1 *fTheory = new TF1("fTheory",
        "[0]*((2-x)/x + (1-x)*x/((2-x)*(2-x))"
        "- 2*(1-x)*(1-x)*log(1-x)/pow(2-x,3)"
        "+ 2*(1-x)*log(1-x)/(x*x))",
        0.001, 0.999);   // avoid singularities

    fTheory->SetParameter(0,1.0);

    // Normalize theory numerically
    double norm = fTheory->Integral(0.001,0.999);
    fTheory->SetParameter(0,1.0/norm);

    fTheory->SetLineColorAlpha(kBlue, 0.5);
    fTheory->SetLineWidth(3);
    fTheory->SetNpx(2000);

    TCanvas *c = new TCanvas("c","Ore-Powell Comparison",800,600);

    gStyle->SetOptStat(0);

    hData->SetLineWidth(2);
    hData->SetMarkerStyle(20);
    hData->Draw("E");

    fTheory->Draw("same");

    TLegend *leg = new TLegend(0.60,0.75,0.88,0.88);
    leg->AddEntry(hData,"Simulation","lep");
    leg->AddEntry(fTheory,"Ore-Powell theory","l");
    leg->Draw();

    c->SaveAs("powell_plot.pdf[300]");

    std::cout << "Unique photons counted: "
              << counted.size() << std::endl;
}
