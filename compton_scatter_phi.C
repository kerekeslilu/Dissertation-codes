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

/*
    does not used stored phi for delta phi relations.
    calculates a normal to the decay plane BEFORE the gammas have scattered, and
    relates phi to that normal vector AFTER they have compton scattered
*/


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
    Float_t Theta;    //compton polar scattering angle

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
        TVector3 vertex;   // gamma creation position
        TVector3 hit;      // first interaction position
        double theta;      // compton polar angle
    };

    std::map<int, std::map<int, HitData>> events;

    //record the first compton scatter
    for (Long64_t i = 0; i < nentries; ++i)
    {
        hitsTree->GetEntry(i);

        // only keep compton scatters
        if (strcmp(HitProcess,"compt") != 0)
            continue;

        // if this gamma has not yet been stored for this event
        if (events[EventID].count(GammaIndex) == 0)
        {
            HitData data;

            data.vertex = TVector3(VertexX,VertexY,VertexZ);
            data.hit    = TVector3(HitX,HitY,HitZ);
            data.theta  = Theta;
            
            events[EventID][GammaIndex] = data;
        }
    }

    TH1F *hDphi12 = new TH1F(
        "hDphi12",
        "#Delta#phi_{12};#Delta#phi (deg);Counts",
        180, -180, 180);

    TH1F *hDphi23 = new TH1F(
        "hDphi23",
        "#Delta#phi_{23};#Delta#phi (deg);Counts",
        180, -180, 180);
    
    TH1F *hDphi13 = new TH1F(
        "hDphi13",
        "#Delta#phi_{13};#Delta#phi (deg);Counts",
        180, -180, 180);

    auto deltaPhi = [](double a, double b)
    {
        double d = a - b;
        d = std::remainder(d,360.0);
        return d;
    };
    int filledEvents = 0;

    for (auto &evt : events)
    {

        // only accept 3 gamma events
        if (evt.second.size() != 3)
            continue;

        auto g0 = evt.second[0];
        auto g1 = evt.second[1];
        auto g2 = evt.second[2];

        // limit the theta differences which we accept 160-180 to emulate backscatter
        //if (g0.theta < 160 || g0.theta > 180) continue;
        //if (g1.theta < 160 || g1.theta > 180) continue;
        //if (g2.theta < 160 || g2.theta > 180) continue;
        
        // incoming gamma directions
        TVector3 k0 = (g0.hit - g0.vertex).Unit();
        TVector3 k1 = (g1.hit - g1.vertex).Unit();
        TVector3 k2 = (g2.hit - g2.vertex).Unit();

        // reconstruct decay plane using the three gamma directions
        TVector3 nDecay = (k0 - k1).Cross(k0 - k2).Unit();

        // phi relative to decay plane
        auto computePhi = [&](TVector3 kin, double thetaDeg)
        {
            TVector3 perp = kin.Orthogonal().Unit();

            // create scattered direction vector
            TVector3 kout = kin;
            kout.Rotate(thetaDeg*TMath::DegToRad(), perp);
            // local coordinate system vectors
            TVector3 x = nDecay.Cross(kin).Unit();
            TVector3 y = kin.Cross(x);


            // azimuthal angle

            double phi = atan2(kout.Dot(y),kout.Dot(x));

            return phi * TMath::RadToDeg();
        };

        // phi from the normal and the azimuthal scatter angle
        double phi0 = computePhi(k0,g0.theta);
        double phi1 = computePhi(k1,g1.theta);
        double phi2 = computePhi(k2,g2.theta);

        // delta phi correlations
        double dphi12 = deltaPhi(phi0,phi1);
        double dphi23 = deltaPhi(phi1,phi2);
        double dphi13 = deltaPhi(phi0,phi2);
    

        hDphi12->Fill(dphi12);
        hDphi23->Fill(dphi23);
        hDphi13->Fill(dphi13);

        filledEvents++;
    }

    std::cout << "Filled events: " << filledEvents << std::endl;

    // fit distributions with cos(2dphi)
    /*
    TF1 *fitFcn12 = new TF1(
        "fitFcn12",
        "[0]*cos(2*TMath::DegToRad()*x) + [1]",
        -180,180);

    fitFcn12->SetParameters(
        0.1*hDphi12->GetMaximum(),
        hDphi12->GetMaximum());

    hDphi12->Fit(fitFcn12,"R");

    TF1 *fitFcn23 = new TF1(
        "fitFcn23",
        "[0]*cos(2*TMath::DegToRad()*x) + [1]",
        -180,180);

    fitFcn23->SetParameters(
        0.1*hDphi23->GetMaximum(),
        hDphi23->GetMaximum());

    hDphi23->Fit(fitFcn23,"R");

    // calculate enhancement factor
    double A12 = fitFcn12->GetParameter(0);
    double B12 = fitFcn12->GetParameter(1);
    double R12 = (B12 - A12)/(B12 + A12);
    std::cout << "\nPair 12 enhancement factor R = " << R12 << std::endl;

    double A23 = fitFcn23->GetParameter(0);
    double B23 = fitFcn23->GetParameter(1);
    double R23 = (B23 - A23)/(B23 + A23);
    std::cout << "Pair 23 enhancement factor R = " << R23 << std::endl;
     */
    

    TCanvas *c = new TCanvas("c","Delta Phi",2400,1000);
    c->Divide(3,1);
    
    c->cd(1);
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.14);
    hDphi12->Draw();

    c->cd(2);
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.14);
    hDphi23->Draw();

    c->cd(3);
    gPad->SetLeftMargin(0.14);
    gPad->SetBottomMargin(0.14);
    hDphi13->Draw();
    
    c->SaveAs("dphi_compton_reconstructed_phi.png");
}
