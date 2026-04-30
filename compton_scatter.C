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

void compton_scatter() {

    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    TChain *hitsTree = new TChain("Hits");
    hitsTree->Add("out_*.root");

    if (hitsTree->GetEntries() == 0) {
        std::cerr << "No entries found!" << std::endl;
        return;
    }

    Int_t   EventID;
    Int_t   GammaIndex;
    Char_t  HitProcess[20];
    Float_t Theta;
    Float_t Phi;

    hitsTree->SetBranchAddress("EventID",    &EventID);
    hitsTree->SetBranchAddress("GammaIndex", &GammaIndex);
    hitsTree->SetBranchAddress("HitProcess", HitProcess);
    hitsTree->SetBranchAddress("Theta",      &Theta);
    hitsTree->SetBranchAddress("Phi",        &Phi);

    Long64_t nentries = hitsTree->GetEntries();

    // Store first Compton interaction (Theta, Phi)
    std::map<int, std::map<int, std::pair<double,double>>> events;

    for (Long64_t i = 0; i < nentries; ++i) {

        hitsTree->GetEntry(i);

        if (strcmp(HitProcess,"compt") != 0)
            continue;

        if (events[EventID].count(GammaIndex) == 0) {
            events[EventID][GammaIndex] = std::make_pair(Theta, Phi);
        }
    }

    TH1F *hDphi12 = new TH1F("hDphi12",
        "#Delta#phi_{12};#Delta#phi (deg);Counts",
        180, -180, 180);

    TH1F *hDphi23 = new TH1F("hDphi23",
        "#Delta#phi_{23};#Delta#phi (deg);Counts",
        180, -180, 180);

    auto deltaPhi = [](double a, double b) {
        double d = a - b;
        d = std::remainder(d, 360.0);
        return d;
    };

    int filledEvents = 0;

    for (auto &evt : events) {

        if (evt.second.size() == 3) {

            double theta0 = evt.second[0].first;
            double phi0   = evt.second[0].second;

            double theta1 = evt.second[1].first;
            double phi1   = evt.second[1].second;

            double theta2 = evt.second[2].first;
            double phi2   = evt.second[2].second;

            // Apply theta cuts (160°–180°)
            if (theta0 < 160 || theta0 > 180) continue;
            if (theta1 < 160 || theta1 > 180) continue;
            if (theta2 < 160 || theta2 > 180) continue;

            double dphi12 = deltaPhi(phi0, phi1);
            double dphi23 = deltaPhi(phi1, phi2);

            hDphi12->Fill(dphi12);
            hDphi23->Fill(dphi23);

            filledEvents++;
        }
    }

    std::cout << "Filled events: " << filledEvents << std::endl;

    TF1 *fitFcn12 = new TF1(
        "fitFcn12",
        "[0]*cos(2*TMath::DegToRad()*x) + [1]",
        -180, 180);

    fitFcn12->SetParameters(
        0.1 * hDphi12->GetMaximum(),
        hDphi12->GetMaximum());

    hDphi12->Fit(fitFcn12, "R");

    double A12  = fitFcn12->GetParameter(0);
    double B12  = fitFcn12->GetParameter(1);
    double eA12 = fitFcn12->GetParError(0);
    double eB12 = fitFcn12->GetParError(1);

    double R12  = (B12 - A12) / (B12 + A12);

    double eR12 = 2.0 / pow(B12 + A12,2) *
        sqrt( pow(B12*eA12,2) + pow(A12*eB12,2) );

    std::cout << "\nPair 1–2" << std::endl;
    std::cout << "A12 = " << A12 << " ± " << eA12 << std::endl;
    std::cout << "B12 = " << B12 << " ± " << eB12 << std::endl;
    std::cout << "Enhancement factor R12 = "
              << R12 << " ± " << eR12 << std::endl;

    TF1 *fitFcn23 = new TF1(
        "fitFcn23",
        "[0]*cos(2*TMath::DegToRad()*x) + [1]",
        -180, 180);

    fitFcn23->SetParameters(
        0.1 * hDphi23->GetMaximum(),
        hDphi23->GetMaximum());

    hDphi23->Fit(fitFcn23, "R");

    double A23  = fitFcn23->GetParameter(0);
    double B23  = fitFcn23->GetParameter(1);
    double eA23 = fitFcn23->GetParError(0);
    double eB23 = fitFcn23->GetParError(1);

    double R23  = (B23 - A23) / (B23 + A23);

    double eR23 = 2.0 / pow(B23 + A23,2) *
        sqrt( pow(B23*eA23,2) + pow(A23*eB23,2) );

    std::cout << "\nPair 2–3" << std::endl;
    std::cout << "A23 = " << A23 << " ± " << eA23 << std::endl;
    std::cout << "B23 = " << B23 << " ± " << eB23 << std::endl;
    std::cout << "Enhancement factor R23 = "
              << R23 << " ± " << eR23 << std::endl;

    TCanvas *c = new TCanvas("c", "Delta Phi", 1200, 500);
    c->Divide(2,1);

    c->cd(1);
    hDphi12->Draw();

    c->cd(2);
    hDphi23->Draw();

    c->SaveAs("dphi_compton_R.png");
}
