#include <iostream>
#include <map>
#include <cmath>
#include <cstring>

#include "TChain.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TStyle.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLegend.h"

void compton_scatter_phi()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    TChain *hitsTree = new TChain("Hits");
    hitsTree->Add("out_*.root");

    if (hitsTree->GetEntries() == 0)
    {
        std::cerr << "No entries found!" << std::endl;
        return;
    }

    Int_t   EventID;
    Int_t   GammaIndex;
    Char_t  HitProcess[20];
    Float_t VertexX, VertexY, VertexZ;
    Float_t HitX, HitY, HitZ;
    Float_t Theta;

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

    Long64_t nentries = hitsTree->GetEntries();

    struct HitData
    {
        TVector3 vertex;
        TVector3 hit1;
        TVector3 hit2;
        double theta;
        bool hasFirst=false;
        bool hasSecond=false;
    };

    std::map<int, std::map<int, HitData>> events;

    // --- Step 1: Build Event Map ---
    for (Long64_t i=0;i<nentries;++i)
    {
        hitsTree->GetEntry(i);
        TVector3 vertex(VertexX,VertexY,VertexZ);
        TVector3 hit(HitX,HitY,HitZ);
        auto &g = events[EventID][GammaIndex];

        if(strcmp(HitProcess,"compt")==0 && !g.hasFirst)
        {
            g.vertex = vertex;
            g.hit1   = hit;
            g.theta  = Theta;
            g.hasFirst = true;
            continue;
        }
        if(g.hasFirst && !g.hasSecond)
        {
            g.hit2 = hit;
            g.hasSecond = true;
        }
    }

    // --- Step 2: Define Histograms ---
    // Added histograms for Symmetric Selection
    TH1F *hDphi12_sym = new TH1F("hDphi12_sym","#Delta#phi_{12} (Symmetric Selection);#Delta#phi (deg);Counts",100,-180,180);
    TH1F *hDphi23_sym = new TH1F("hDphi23_sym","#Delta#phi_{23} (Symmetric Selection);#Delta#phi (deg);Counts",100,-180,180);
    TH1F *hDphi13_sym = new TH1F("hDphi13_sym","#Delta#phi_{13} (Symmetric Selection);#Delta#phi (deg);Counts",100,-180,180);

    TH1F *hDphi12_sel = new TH1F("hDphi12_sel","#Delta#phi_{12} (80^{#circ} < #theta < 90^{#circ});#Delta#phi (deg);Counts",50,-180,180);
    TH1F *hDphi23_sel = new TH1F("hDphi23_sel","#Delta#phi_{23} (80^{#circ} < #theta < 90^{#circ});#Delta#phi (deg);Counts",50,-180,180);
    TH1F *hDphi13_sel = new TH1F("hDphi13_sel","#Delta#phi_{13} (80^{#circ} < #theta < 90^{#circ});#Delta#phi (deg);Counts",50,-180,180);

    auto deltaPhi = [](double a,double b)
    {
        double d=a-b;
        d=std::remainder(d,360.0);
        return d;
    };

    int filledEvents=0;
    int symEvents=0;

    // --- Step 3: Analysis Loop ---
    for(auto &evt:events)
    {
        if(evt.second.size()!=3) continue;

        auto &g0 = evt.second[0];
        auto &g1 = evt.second[1];
        auto &g2 = evt.second[2];

        if(!g0.hasSecond || !g1.hasSecond || !g2.hasSecond) continue;

        // Define Kinematics
        TVector3 kin0  = (g0.hit1 - g0.vertex).Unit();
        TVector3 kout0 = (g0.hit2 - g0.hit1).Unit();
        TVector3 kin1  = (g1.hit1 - g1.vertex).Unit();
        TVector3 kout1 = (g1.hit2 - g1.hit1).Unit();
        TVector3 kin2  = (g2.hit1 - g2.vertex).Unit();
        TVector3 kout2 = (g2.hit2 - g2.hit1).Unit();

        // Calculate Angles Between Photons (for Symmetric Cut)
        double angle12 = kin0.Angle(kin1);
        double angle23 = kin1.Angle(kin2);
        double angle13 = kin0.Angle(kin2);

        // Calculate Decay Plane Normal
        TVector3 nDecay = (kin0-kin1).Cross(kin0-kin2).Unit();
        if (nDecay.Mag() == 0) continue; // Skip collinear events

        // Compute Phi for each photon
        auto computePhi = [&](TVector3 kin,TVector3 kout)
        {
            TVector3 x = nDecay.Cross(kin).Unit();
            TVector3 y = kin.Cross(x);
            double phi = atan2(kout.Dot(y), kout.Dot(x));
            return phi*TMath::RadToDeg();
        };

        double phi0 = computePhi(kin0,kout0);
        double phi1 = computePhi(kin1,kout1);
        double phi2 = computePhi(kin2,kout2);

        double dphi12 = deltaPhi(phi0,phi1);
        double dphi23 = deltaPhi(phi1,phi2);
        double dphi13 = deltaPhi(phi0,phi2);

        // --- CUT 1: Symmetric Kinematics (Dalitz Center) ---
        // Approx 120 degrees = 2.094 rad. Window 115-125 deg (2.0 - 2.2 rad)
        bool isSymmetric = (angle12 > 2.0 && angle12 < 2.2 &&
                            angle23 > 2.0 && angle23 < 2.2 &&
                            angle13 > 2.0 && angle13 < 2.2);
        
        if (isSymmetric) {
            hDphi12_sym->Fill(dphi12);
            hDphi23_sym->Fill(dphi23);
            hDphi13_sym->Fill(dphi13);
            symEvents++;
        }

        // --- CUT 2: Polar Angle Selection ---
        if (g0.theta >= 80 && g0.theta <= 90 &&
            g1.theta >= 80 && g1.theta <= 90 &&
            g2.theta >= 80 && g2.theta <= 90)
        {
            hDphi12_sel->Fill(dphi12);
            hDphi23_sel->Fill(dphi23);
            hDphi13_sel->Fill(dphi13);
        }
        
        filledEvents++;
    }

    std::cout<<"Total Valid Events: "<<filledEvents<<std::endl;
    std::cout<<"Symmetric Events: "<<symEvents<<std::endl;

    // --- Step 4: Fitting ---
    TF1 *fitFcn = new TF1("fitFcn","[0]*cos(2*TMath::DegToRad()*x)+[1]",-180,180);
    
    // Fit Symmetric Histograms
    TH1F* hSymList[] = {hDphi12_sym, hDphi23_sym, hDphi13_sym};
    const char* namesSym[] = {"12", "23", "13"};
    
    TCanvas *c_sym = new TCanvas("c_sym","Symmetric Selection",2400,1000);
    c_sym->Divide(3,1);
    
    for(int i=0; i<3; i++) {
        c_sym->cd(i+1);
        hSymList[i]->Draw("E1 HIST");
        TF1* fit = (TF1*)fitFcn->Clone(Form("fit_sym_%s", namesSym[i]));
        hSymList[i]->Fit(fit, "R");
        
        double A = fit->GetParameter(0);
        double B = fit->GetParameter(1);
        double R = (B+A)/(B-A);
        
        TLegend *leg = new TLegend(0.65, 0.78, 0.88, 0.88);
        leg->AddEntry("", Form("R = %.3f", R), "");
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->Draw();
    }
    c_sym->SaveAs("dphi_symmetric_selection.pdf");

    // Fit Theta Selection Histograms (Existing logic)
    TH1F* hSelList[] = {hDphi12_sel, hDphi23_sel, hDphi13_sel};
    TCanvas *c_sel = new TCanvas("c_sel","Theta Selection",2400,1000);
    c_sel->Divide(3,1);

    for(int i=0; i<3; i++) {
        c_sel->cd(i+1);
        hSelList[i]->Draw("E1 HIST");
        TF1* fit = (TF1*)fitFcn->Clone(Form("fit_sel_%s", namesSym[i]));
        hSelList[i]->Fit(fit, "R");
        
        double A = fit->GetParameter(0);
        double B = fit->GetParameter(1);
        double R = (B+A)/(B-A);

        TLegend *leg = new TLegend(0.65, 0.78, 0.88, 0.88);
        leg->AddEntry("", Form("R = %.3f", R), "");
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->Draw();
    }
    c_sel->SaveAs("dphi_theta_selection.pdf");
}
