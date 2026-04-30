#include "TCanvas.h"
#include "TPad.h"
#include "TView.h"
#include "TLine.h"
#include "TText.h"
#include "TArc.h"
#include "TMath.h"
#include "TVector3.h"
#include "TStyle.h"
#include "TH3F.h"
#include "TH2F.h"
#include "TPolyLine3D.h"

void draw_phi_geometry() {
    gStyle->SetOptStat(0);
    gStyle->SetCanvasBorderMode(0);

    TCanvas *c = new TCanvas("c", "Phi Calculation Geometry", 1600, 800);
    c->Divide(2, 1);

    // ============================================================
    // PAD 1: 3D Perspective View
    // ============================================================
    c->cd(1);
    
    // 1. Create dummy 3D histogram to establish the 3D frame
    TH3F *hFrame = new TH3F("hFrame", "3D Decay Geometry;X;Y;Z", 1, -2.5, 2.5, 1, -2.5, 2.5, 1, -2.5, 2.5);
    hFrame->SetStats(0);
    hFrame->Draw(" ");
    
    // 2. Set View Angles
    gPad->SetTheta(30);
    gPad->SetPhi(45);
    gPad->Update();

    // --- Define Vectors ---
    TVector3 k_in(0, 0, 1.8);
    TVector3 n_decay(0, 1.2, 0.6);
    n_decay = n_decay.Unit() * 1.8;

    TVector3 x_prime = n_decay.Cross(k_in);
    if (x_prime.Mag() == 0) x_prime = TVector3(1, 0, 0);
    x_prime = x_prime.Unit() * 1.8;

    TVector3 y_prime = k_in.Cross(x_prime);
    y_prime = y_prime.Unit() * 1.8;

    double phi_demo = 45.0 * TMath::DegToRad();
    double theta_scatt = 60.0 * TMath::DegToRad();
    
    double k_out_x = sin(theta_scatt) * cos(phi_demo);
    double k_out_y = sin(theta_scatt) * sin(phi_demo);
    double k_out_z = cos(theta_scatt);
    
    TVector3 x_hat = x_prime.Unit();
    TVector3 y_hat = y_prime.Unit();
    TVector3 z_hat = k_in.Unit();

    TVector3 k_out = (x_hat * k_out_x) + (y_hat * k_out_y) + (z_hat * k_out_z);
    k_out = k_out.Unit() * 1.8;

    // --- Draw Vectors using TPolyLine3D ---
    auto drawVec3D = [](Double_t x2, Double_t y2, Double_t z2, Color_t col, Int_t width = 2, Int_t style = 1) {
        Double_t pts[6] = {0, 0, 0, x2, y2, z2};
        TPolyLine3D *line = new TPolyLine3D(2, pts);
        line->SetLineColor(col);
        line->SetLineWidth(width);
        line->SetLineStyle(style);
        line->Draw();
    };

    drawVec3D(k_in.X(), k_in.Y(), k_in.Z(), kBlack, 3);
    drawVec3D(n_decay.X(), n_decay.Y(), n_decay.Z(), kRed, 2, 2);
    drawVec3D(x_prime.X(), x_prime.Y(), x_prime.Z(), kBlue, 2);
    drawVec3D(y_prime.X(), y_prime.Y(), y_prime.Z(), kGreen+2, 2);
    drawVec3D(k_out.X(), k_out.Y(), k_out.Z(), kMagenta, 3);

    // --- Draw Labels (2D Text positioned to align with 3D view) ---
    TText *tKin = new TText(0.1, 1.6, "k_{in}"); tKin->SetTextSize(0.04); tKin->Draw();
    TText *tN = new TText(0.6, 1.0, "n_{Decay}"); tN->SetTextColor(kRed); tN->Draw();
    TText *tX = new TText(-0.8, 0.5, "x'"); tX->SetTextColor(kBlue); tX->Draw();
    TText *tY = new TText(0.2, -0.5, "y'"); tY->SetTextColor(kGreen+2); tY->Draw();
    TText *tOut = new TText(-0.5, -0.8, "k_{out}"); tOut->SetTextColor(kMagenta); tOut->Draw();

    TText *tTitle1 = new TText(0, 2.2, "3D Geometry (Rotatable)");
    tTitle1->SetTextSize(0.05); tTitle1->SetTextAlign(22); tTitle1->Draw();

    // ============================================================
    // PAD 2: 2D Projection
    // ============================================================
    c->cd(2);
    
    TH2F *hFrame2D = new TH2F("hFrame2D", "Projection Perpendicular to k_{in};x';y'", 100, -2, 2, 100, -2, 2);
    hFrame2D->SetStats(0);
    hFrame2D->Draw(" ");

    gPad->SetFixedAspectRatio(kTRUE);
    
    // Axes
    TLine *lX = new TLine(-1.8, 0, 1.8, 0);
    lX->SetLineColor(kBlue); lX->SetLineWidth(2); lX->Draw();
    
    TLine *lY = new TLine(0, -1.8, 0, 1.8);
    lY->SetLineColor(kGreen+2); lY->SetLineWidth(2); lY->Draw();

    // FIX: Renamed variables to tX2 and tY2 to avoid redefinition errors
    TText *tX2 = new TText(1.9, 0, "x'");
    tX2->SetTextColor(kBlue); tX2->SetTextAlign(12); tX2->Draw();
    
    TText *tY2 = new TText(0.1, 1.9, "y'");
    tY2->SetTextColor(kGreen+2); tY2->Draw();

    // Vector
    double projLen = 1.2;
    double endX = projLen * cos(phi_demo);
    double endY = projLen * sin(phi_demo);

    TLine *lProj = new TLine(0, 0, endX, endY);
    lProj->SetLineColor(kMagenta); lProj->SetLineWidth(3); lProj->Draw();
    
    // Arrowhead
    double headLen = 0.15;
    double angle = atan2(endY, endX);
    TLine *h1 = new TLine(endX, endY, endX - headLen * cos(angle - 0.5), endY - headLen * sin(angle - 0.5));
    TLine *h2 = new TLine(endX, endY, endX - headLen * cos(angle + 0.5), endY - headLen * sin(angle + 0.5));
    h1->SetLineColor(kMagenta); h1->SetLineWidth(3); h1->Draw();
    h2->SetLineColor(kMagenta); h2->SetLineWidth(3); h2->Draw();

    TText *tProj = new TText(endX*1.1, endY*1.1, "k_{out}");
    tProj->SetTextColor(kMagenta); tProj->Draw();

    // Arc
    TArc *arc = new TArc(0, 0, 0.6, 0, TMath::RadToDeg()*phi_demo);
    arc->SetLineColor(kBlack); arc->SetLineWidth(2); arc->SetFillStyle(0); arc->Draw();
    
    TText *tPhi = new TText(0.7 * cos(phi_demo/2.0), 0.7 * sin(phi_demo/2.0), "#phi");
    tPhi->SetTextSize(0.06); tPhi->Draw();

    TText *tDesc = new TText(0, -2.2, "#phi = atan2(k_{out}#cdot;y', k_{out}#cdot;x')");
    tDesc->SetTextSize(0.04); tDesc->SetTextAlign(22); tDesc->Draw();

    c->Update();
    c->SaveAs("phi_geometry_schematic.pdf");
}
