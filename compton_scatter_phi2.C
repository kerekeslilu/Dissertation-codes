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
    //Float_t Phi;

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
    //hitsTree->SetBranchAddress("Phi", &Phi);

    Long64_t nentries = hitsTree->GetEntries();

    struct HitData
    {
        TVector3 vertex;
        TVector3 hit1;
        TVector3 hit2;

        double theta;
        //double phi;

        bool hasFirst=false;
        bool hasSecond=false;
    };

    std::map<int, std::map<int, HitData>> events;

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
            //g.phi    = Phi;

            g.hasFirst = true;
            continue;
        }

        if(g.hasFirst && !g.hasSecond)
        {
            g.hit2 = hit;
            g.hasSecond = true;
        }
    }

    // NO SELECTION
    TH1F *hDphi12 = new TH1F("hDphi12","#Delta#phi_{12} (No Selection, No QE);#Delta#phi (deg);Counts",160,-180,180);
    TH1F *hDphi23 = new TH1F("hDphi23","#Delta#phi_{23} (No Selection, No QE);#Delta#phi (deg);Counts",160,-180,180);
    TH1F *hDphi13 = new TH1F("hDphi13","#Delta#phi_{13} (No Selection, No QE);#Delta#phi (deg);Counts",160,-180,180);

    // WITH SELECTION
    TH1F *hDphi12_sel = new TH1F("hDphi12_sel","#Delta#phi_{12} (80^{#circ} < #theta < 90^{#circ}, No QE);#Delta#phi (deg);Counts",40,-180,180);
    TH1F *hDphi23_sel = new TH1F("hDphi23_sel","#Delta#phi_{23} (80^{#circ} < #theta < 90^{#circ}, No QE);#Delta#phi (deg);Counts",40,-180,180);
    TH1F *hDphi13_sel = new TH1F("hDphi13_sel","#Delta#phi_{13} (80^{#circ} < #theta < 90^{#circ}, No QE);#Delta#phi (deg);Counts",40,-180,180);

    //TH1F *hPol_mm1 = new TH1F("hPol_mm1","Polarisation m_{s}=-1;#phi (deg);Counts",180,-180,180);
    //TH1F *hPol_m0  = new TH1F("hPol_m0","Polarisation m_{s}=0;#phi (deg);Counts",180,-180,180);
    //TH1F *hPol_mp1 = new TH1F("hPol_mp1","Polarisation m_{s}=+1;#phi (deg);Counts",180,-180,180);

    auto deltaPhi = [](double a,double b)
    {
        double d=a-b;
        d=std::remainder(d,360.0);
        return d;
    };

    int filledEvents=0;

    for(auto &evt:events)
    {
        if(evt.second.size()!=3)
            continue;

        auto &g0 = evt.second[0];
        auto &g1 = evt.second[1];
        auto &g2 = evt.second[2];

        if(!g0.hasSecond || !g1.hasSecond || !g2.hasSecond)
            continue;

        // limit the theta differences which we accept 60-120 to emulate backscatter
        /*if (g0.theta < 60 || g0.theta > 120) continue;
        if (g1.theta < 60 || g1.theta > 120) continue;
        if (g2.theta < 60 || g2.theta > 120) continue;*/

        TVector3 kin0  = (g0.hit1 - g0.vertex).Unit();
        TVector3 kout0 = (g0.hit2 - g0.hit1).Unit();

        TVector3 kin1  = (g1.hit1 - g1.vertex).Unit();
        TVector3 kout1 = (g1.hit2 - g1.hit1).Unit();

        TVector3 kin2  = (g2.hit1 - g2.vertex).Unit();
        TVector3 kout2 = (g2.hit2 - g2.hit1).Unit();

        TVector3 nDecay = (kin0-kin1).Cross(kin0-kin2).Unit();

        auto computePhi = [&](TVector3 kin,TVector3 kout)
        {
            TVector3 x = nDecay.Cross(kin).Unit();
            TVector3 y = kin.Cross(x);

            double phi = atan2(kout.Dot(y),kout.Dot(x));

            return phi*TMath::RadToDeg();
        };

        double phi0 = computePhi(kin0,kout0);
        double phi1 = computePhi(kin1,kout1);
        double phi2 = computePhi(kin2,kout2);

        double dphi12 = deltaPhi(phi0,phi1);
        double dphi23 = deltaPhi(phi1,phi2);
        double dphi13 = deltaPhi(phi0,phi2);

        // NO selection
        hDphi12->Fill(dphi12);
        hDphi23->Fill(dphi23);
        hDphi13->Fill(dphi13);

        // WITH selection
        if (g0.theta >= 80 && g0.theta <= 90 &&
            g1.theta >= 80 && g1.theta <= 90 &&
            g2.theta >= 80 && g2.theta <= 90)
        {
            hDphi12_sel->Fill(dphi12);
            hDphi23_sel->Fill(dphi23);
            hDphi13_sel->Fill(dphi13);
        }

        // fill polarisation from simulation phi
        //hPol_m0->Fill(g0.phi);
        //hPol_mp1->Fill(g1.phi);
        //hPol_mm1->Fill(g2.phi);

        filledEvents++;
    }

    std::cout<<"Filled events: "<<filledEvents<<std::endl;

    TF1 *fitFcn = new TF1("fitFcn","[0]*cos(2*TMath::DegToRad()*x)+[1]",-180,180);

    TF1 *fit12=(TF1*)fitFcn->Clone("fit12");
    TF1 *fit23=(TF1*)fitFcn->Clone("fit23");
    TF1 *fit13=(TF1*)fitFcn->Clone("fit13");

    TF1 *fit12_sel=(TF1*)fitFcn->Clone("fit12_sel");
    TF1 *fit23_sel=(TF1*)fitFcn->Clone("fit23_sel");
    TF1 *fit13_sel=(TF1*)fitFcn->Clone("fit13_sel");

    /*TF1 *fit_m0 =(TF1*)fitFcn->Clone("fit_m0");
    TF1 *fit_mp1=(TF1*)fitFcn->Clone("fit_mp1");
    TF1 *fit_mm1=(TF1*)fitFcn->Clone("fit_mm1");*/

    hDphi12->Fit(fit12,"R");
    hDphi23->Fit(fit23,"R");
    hDphi13->Fit(fit13,"R");

    hDphi12_sel->Fit(fit12_sel,"R");
    hDphi23_sel->Fit(fit23_sel,"R");
    hDphi13_sel->Fit(fit13_sel,"R");

   // hPol_m0->Fit(fit_m0,"R");
    //hPol_mp1->Fit(fit_mp1,"R");
    //hPol_mm1->Fit(fit_mm1,"R");

    auto getR = [](TF1* f)
    {
        double A=f->GetParameter(0);
        double B=f->GetParameter(1);

        double eA=f->GetParError(0);
        double eB=f->GetParError(1);

        double R=(B+A)/(B-A);
        double eR = R*sqrt(pow((-2*B) / pow(B + A, 2) * eA, 2) + pow((2*A) / pow(B + A, 2) * eB, 2));
        std::cout << "A = " << A << ", B = " << B << ", R = " << R << std::endl;

        return std::pair<double,double>(R,eR);
    };

    auto R12 = getR(fit12);
    auto R23 = getR(fit23);
    auto R13 = getR(fit13);

    auto R12_sel = getR(fit12_sel);
    auto R23_sel = getR(fit23_sel);
    auto R13_sel = getR(fit13_sel);

    TCanvas *c = new TCanvas("c","Delta Phi (No Selection)",2400,1000);
    c->Divide(3,1);

    c->cd(1);
    hDphi12->Draw("E1 HIST");
    fit12->Draw("same");
    {
        /*TLegend *leg = new TLegend(0.65, 0.78, 0.88, 0.88);
        leg->AddEntry((TObject*)0,Form("Events = %.0f", hDphi12->GetEntries()),"");
        leg->AddEntry((TObject*)0,Form("R = %.3f #pm %.3f",R12.first,R12.second),"");
        leg->SetTextSize(0.045);
        leg->SetMargin(0.1);
        leg->SetEntrySeparation(0.01);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->Draw();*/
    }

    c->cd(2);
    hDphi23->Draw("E1 HIST");
    fit23->Draw("same");
    {
        /*TLegend *leg = new TLegend(0.65, 0.78, 0.88, 0.88);
        leg->AddEntry((TObject*)0,Form("Events = %.0f", hDphi23->GetEntries()),"");
        leg->AddEntry((TObject*)0,Form("R = %.3f #pm %.3f",R23.first,R23.second),"");
        leg->SetTextSize(0.045);
        leg->SetMargin(0.1);
        leg->SetEntrySeparation(0.01);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->Draw();*/
    }

    c->cd(3);
    hDphi13->Draw("E1 HIST");
    fit13->Draw("same");
    {
        /*TLegend *leg = new TLegend(0.65, 0.78, 0.88, 0.88);
        leg->AddEntry((TObject*)0,Form("Events = %.0f", hDphi13->GetEntries()),"");
        leg->AddEntry((TObject*)0,Form("R = %.3f #pm %.3f",R13.first,R13.second),"");
        leg->SetTextSize(0.045);
        leg->SetMargin(0.1);
        leg->SetEntrySeparation(0.01);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->Draw();*/
    }

    c->SaveAs("dphi_compton_noq_0_180_50.pdf");

    TCanvas *c_sel = new TCanvas("c_sel","Delta Phi (80-90 Selection)",2400,1000);
    c_sel->Divide(3,1);

    c_sel->cd(1);
    hDphi12_sel->Draw("E1 HIST");
    fit12_sel->Draw("same");
    {
        /*TLegend *leg = new TLegend(0.65, 0.78, 0.88, 0.88);
        leg->AddEntry((TObject*)0, Form("Events = %.0f", hDphi12_sel->GetEntries()), "");
        leg->AddEntry((TObject*)0, Form("R = %.3f #pm %.3f", R12_sel.first, R12_sel.second), "");
        leg->SetTextSize(0.045);
        leg->SetMargin(0.1);
        leg->SetEntrySeparation(0.01);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->Draw();*/
    }

    c_sel->cd(2);
    hDphi23_sel->Draw("E1 HIST");
    fit23_sel->Draw("same");
    {
        /*TLegend *leg = new TLegend(0.65, 0.78, 0.88, 0.88);
        leg->AddEntry((TObject*)0, Form("Events = %.0f", hDphi23_sel->GetEntries()), "");
        leg->AddEntry((TObject*)0, Form("R = %.3f #pm %.3f", R23_sel.first, R23_sel.second), "");
        leg->SetTextSize(0.045);
        leg->SetMargin(0.1);
        leg->SetEntrySeparation(0.01);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->Draw();*/
    }

    c_sel->cd(3);
    hDphi13_sel->Draw("E1 HIST");
    fit13_sel->Draw("same");
    {
        /*TLegend *leg = new TLegend(0.65, 0.78, 0.88, 0.88);
        leg->AddEntry((TObject*)0, Form("Events = %.0f", hDphi13_sel->GetEntries()), "");
        leg->AddEntry((TObject*)0, Form("R = %.3f #pm %.3f", R13_sel.first, R13_sel.second), "");
        leg->SetTextSize(0.045);
        leg->SetMargin(0.1);
        leg->SetEntrySeparation(0.01);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->Draw();*/
    }

    c_sel->SaveAs("dphi_compton_noq_80_90_50.pdf");
}
