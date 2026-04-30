void filter() {
    #include <map>
    #include <set>

    TChain *t = new TChain("Hits");
    t->Add("out_*.root");
    //TTree *t = (TTree*)f->Get("Hits");
    std::map<int, std::set<int>> eventGammas;

    int EventID, GammaIndex;
    t->SetBranchAddress("EventID", &EventID);
    t->SetBranchAddress("GammaIndex", &GammaIndex);

    Long64_t nentries = t->GetEntries();

    for (Long64_t i = 0; i < nentries; i++) {
        t->GetEntry(i);
        eventGammas[EventID].insert(GammaIndex);
    }

    TFile *out = new TFile("filtered.root", "RECREATE");
    TTree *newTree = t->CloneTree(0);

    for (Long64_t i = 0; i < nentries; i++) {
        t->GetEntry(i);
        if (eventGammas[EventID].size() == 3) {
            newTree->Fill();
        }
    }

    newTree->Write();
    out->Close();
}
