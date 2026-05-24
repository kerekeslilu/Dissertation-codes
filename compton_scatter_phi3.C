#include <iostream>
#include <map>
#include <cmath>
#include <cstring>
#include <vector>

#include "TChain.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLegend.h"
#include "TLatex.h"

void analyze_decay_angles()
{
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);

    
    TChain *hitsTree = new TChain("Hits");
    hitsTree->Add("out_*.root");

    if (hitsTree->GetEntries() == 0)
    {
        std::cerr << "Error: No entries found in out_*.root" << std::endl;
        return;
    }

    
    Int_t   EventID;
    Int_t   GammaIndex;
    Char_t  HitProcess[20];
    Float_t VertexX, VertexY, VertexZ;
    Float_t HitX, HitY, HitZ;
    Float_t Theta; // This is the Compton scattering angle

    hitsTree->SetBranchAddress("EventID",    &EventID);
    hitsTree->SetBranchAddress("GammaIndex", &GammaIndex);
    hitsTree->SetBranchAddress("HitProcess", HitProcess);
    hitsTree->SetBranchAddress("VertexX", &VertexX);
    hitsTree->SetBranchAddress("VertexY", &VertexY);
    hitsTree->SetBranchAddress("VertexZ", &VertexZ);
    hitsTree->SetBranchAddress("HitX", &HitX);
    hitsTree->SetBranchAddress("HitY", &HitY);
    hitsTree->SetBranchAddress("HitZ", &HitZ);
    hitsTree->SetBranchAddress("Theta", &Theta);

    struct HitData
    {
        TVector3 vertex;
        TVector3 hit1;
        TVector3 hit2;
        double theta;
        bool hasFirst = false;
        bool hasSecond = false;
    };

    std::map<int, std::map<int, HitData>> events;


    Long64_t nentries = hitsTree->GetEntries();
    std::cout << "Processing " << nentries << " entries..." << std::endl;

    for (Long64_t i = 0; i < nentries; ++i)
    {
        hitsTree->GetEntry(i);

        // Only process Compton events
        if (strcmp(HitProcess, "compt") != 0) continue;

        TVector3 vertex(VertexX, VertexY, VertexZ);
        TVector3 hit(HitX, HitY, HitZ);

        auto &g = events[EventID][GammaIndex];

        if (!g.hasFirst)
        {
            g.vertex = vertex;
            g.hit1   = hit;
            g.theta  = Theta;
            g.hasFirst = true;
        }
        else if (!g.hasSecond)
        {
            g.hit2 = hit;
            g.hasSecond = true;
        }
    }

    // 
    // Bins from 0 to 180 degrees
    TH1F *hAng12 = new TH1F("hAng12", "Opening Angle #gamma_{1}-#gamma_{2};Angle (deg);Events", 180, 0, 180);
    TH1F *hAng23 = new TH1F("hAng23", "Opening Angle #gamma_{2}-#gamma_{3};Angle (deg);Events", 180, 0, 180);
    TH1F *hAng13 = new TH1F("hAng13", "Opening Angle #gamma_{1}-#gamma_{3};Angle (deg);Events", 180, 0, 180);

    TH1F *hAng12_sel = new TH1F("hAng12_sel", "Opening Angle #gamma_{1}-#gamma_{2} (Sel 82-85^{#circ});Angle (deg);Events", 180, 0, 180);
    TH1F *hAng23_sel = new TH1F("hAng23_sel", "Opening Angle #gamma_{2}-#gamma_{3} (Sel 82-85^{#circ});Angle (deg);Events", 180, 0, 180);
    TH1F *hAng13_sel = new TH1F("hAng13_sel", "Opening Angle #gamma_{1}-#gamma_{3} (Sel 82-85^{#circ});Angle (deg);Events", 180, 0, 180);

    int totalEvents = 0;
    int selectedEvents = 0;

   
    for (auto &evt : events)
    {
        // Require exactly 3 gammas
        if (evt.second.size() != 3) continue;
        if (evt.second.find(0) == evt.second.end() ||
            evt.second.find(1) == evt.second.end() ||
            evt.second.find(2) == evt.second.end()) continue;

        auto &g0 = evt.second[0];
        auto &g1 = evt.second[1];
        auto &g2 = evt.second[2];

        
        if (!g0.hasSecond || !g1.hasSecond || !g2.hasSecond) continue;


        TVector3 kin0 = (g0.hit1 - g0.vertex).Unit();
        TVector3 kin1 = (g1.hit1 - g1.vertex).Unit();
        TVector3 kin2 = (g2.hit1 - g2.vertex).Unit();

     
        double ang12 = kin0.Angle(kin1) * TMath::RadToDeg();
        double ang23 = kin1.Angle(kin2) * TMath::RadToDeg();
        double ang13 = kin0.Angle(kin2) * TMath::RadToDeg();

   
        hAng12->Fill(ang12);
        hAng23->Fill(ang23);
        hAng13->Fill(ang13);

        
        bool passSel = (g0.theta >= 82 && g0.theta <= 85) &&
                       (g1.theta >= 82 && g1.theta <= 85) &&
                       (g2.theta >= 82 && g2.theta <= 85);

        if (passSel)
        {
            hAng12_sel->Fill(ang12);
            hAng23_sel->Fill(ang23);
            hAng13_sel->Fill(ang13);
            selectedEvents++;
        }
        totalEvents++;
    }

    std::cout << "Total valid 3-gamma events: " << totalEvents << std::endl;
    std::cout << "Events passing Compton selection (82-85): " << selectedEvents << std::endl;


    auto printAvg = [](TH1F* h, const char* name) {
        if (h->GetEntries() > 0) {
            std::cout << name << " Mean: " << h->GetMean() << " +/- " << h->GetRMS() << " deg" << std::endl;
        } else {
            std::cout << name << " No entries." << std::endl;
        }
    };

    std::cout << "\n--- Averages (No Selection) ---" << std::endl;
    printAvg(hAng12, "Angle 1-2");
    printAvg(hAng23, "Angle 2-3");
    printAvg(hAng13, "Angle 1-3");

    std::cout << "\n--- Averages (Compton Selection 82-85) ---" << std::endl;
    printAvg(hAng12_sel, "Angle 1-2");
    printAvg(hAng23_sel, "Angle 2-3");
    printAvg(hAng13_sel, "Angle 1-3");


    TCanvas *c = new TCanvas("c", "Decay Opening Angles", 1600, 1200);
    c->Divide(2, 2);


    c->cd(1);
    hAng12->SetLineColor(kBlue);
    hAng23->SetLineColor(kRed);
    hAng13->SetLineColor(kGreen+2);
    
    hAng12->Draw("HIST");
    hAng23->Draw("HIST SAME");
    hAng13->Draw("HIST SAME");

    TLegend *leg1 = new TLegend(0.6, 0.6, 0.88, 0.75);
    leg1->AddEntry(hAng12, "#gamma_{1} - #gamma_{2}", "l");
    leg1->AddEntry(hAng23, "#gamma_{2} - #gamma_{3}", "l");
    leg1->AddEntry(hAng13, "#gamma_{1} - #gamma_{3}", "l");
    leg1->SetBorderSize(0);
    leg1->Draw();

    TLatex *tex1 = new TLatex(0.15, 0.85, "No Selection");
    tex1->SetTextSize(0.05);
    tex1->Draw();


    c->cd(2);
    if (selectedEvents > 0) {
        hAng12_sel->SetLineColor(kBlue);
        hAng23_sel->SetLineColor(kRed);
        hAng13_sel->SetLineColor(kGreen+2);

        hAng12_sel->Draw("HIST");
        hAng23_sel->Draw("HIST SAME");
        hAng13_sel->Draw("HIST SAME");

        TLegend *leg2 = new TLegend(0.6, 0.6, 0.88, 0.75);
        leg2->AddEntry(hAng12_sel, "#gamma_{1} - #gamma_{2}", "l");
        leg2->AddEntry(hAng23_sel, "#gamma_{2} - #gamma_{3}", "l");
        leg2->AddEntry(hAng13_sel, "#gamma_{1} - #gamma_{3}", "l");
        leg2->SetBorderSize(0);
        leg2->Draw();

        TLatex *tex2 = new TLatex(0.15, 0.85, "Compton #theta #in [82, 85]^{#circ}");
        tex2->SetTextSize(0.05);
        tex2->Draw();
    } else {
        TLatex *texNo = new TLatex(0.5, 0.5, "No events passed selection");
        texNo->SetTextAlign(22);
        texNo->Draw();
    }

    c->cd(3);
    TH1F *hDiff = (TH1F*)hAng12->Clone("hDiff");
    hDiff->Add(hAng23, -1);
    hDiff->SetTitle("Difference (Angle 12 - Angle 23);Angle Diff (deg);Counts");
    hDiff->SetLineColor(kBlack);
    hDiff->Draw("HIST");
    gPad->SetGridx();

    c->cd(4);
    if (selectedEvents > 0) {
        TH1F *hDiffSel = (TH1F*)hAng12_sel->Clone("hDiffSel");
        hDiffSel->Add(hAng23_sel, -1);
        hDiffSel->SetTitle("Difference (Angle 12 - Angle 23) Sel;Angle Diff (deg);Counts");
        hDiffSel->SetLineColor(kMagenta);
        hDiffSel->Draw("HIST");
        gPad->SetGridx();
    }

    c->SaveAs("decay_opening_angles.png");
    std::cout << "\nPlot saved as decay_opening_angles.png" << std::endl;
}
