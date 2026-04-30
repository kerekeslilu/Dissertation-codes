#include <TChain.h>
#include <map>
#include <vector>
#include <iostream>
#include <cmath>
#include <string>

struct Vec3 {
    double x, y, z;
};

Vec3 make_dir(double vx, double vy, double vz,
              double hx, double hy, double hz) {
    return {hx - vx, hy - vy, hz - vz};
}

double angle(const Vec3& a, const Vec3& b) {
    double dot = a.x*b.x + a.y*b.y + a.z*b.z;
    double magA = std::sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
    double magB = std::sqrt(b.x*b.x + b.y*b.y + b.z*b.z);
    return std::acos(dot / (magA * magB)); // radians
}

void count() {
    TChain c("Hits");
    c.Add("out_*.root");

    int EventID, GammaIndex;
    char CreationProcess[50];

    float VertexX, VertexY, VertexZ;
    float HitX, HitY, HitZ;

    c.SetBranchAddress("EventID", &EventID);
    c.SetBranchAddress("GammaIndex", &GammaIndex);
    c.SetBranchAddress("CreationProcess", CreationProcess);

    c.SetBranchAddress("VertexX", &VertexX);
    c.SetBranchAddress("VertexY", &VertexY);
    c.SetBranchAddress("VertexZ", &VertexZ);

    c.SetBranchAddress("HitX", &HitX);
    c.SetBranchAddress("HitY", &HitY);
    c.SetBranchAddress("HitZ", &HitZ);

    // Store first hit direction per (event, gamma)
    std::map<int, std::map<int, Vec3>> eventGammas;

    Long64_t n = c.GetEntries();

    for (Long64_t i = 0; i < n; i++) {
        c.GetEntry(i);

        if (std::string(CreationProcess) != "annihil") continue;

        // Only store the FIRST hit per gamma
        if (eventGammas[EventID].count(GammaIndex) == 0) {
            eventGammas[EventID][GammaIndex] =
                make_dir(VertexX, VertexY, VertexZ,
                         HitX, HitY, HitZ);
        }
    }

    double sum12 = 0, sum23 = 0, sum13 = 0;
    int n12 = 0, n23 = 0, n13 = 0;

    for (auto& [evt, gammas] : eventGammas) {
        // We expect GammaIndex = 0,1,2 for 3-photon events
        if (gammas.size() < 3) continue;

        auto& g = gammas;

        if (g.count(0) && g.count(1)) {
            sum12 += angle(g[0], g[1]);
            n12++;
        }
        if (g.count(1) && g.count(2)) {
            sum23 += angle(g[1], g[2]);
            n23++;
        }
        if (g.count(0) && g.count(2)) {
            sum13 += angle(g[0], g[2]);
            n13++;
        }
    }

    std::cout << "Average opening angles (radians):\n";
    if (n12) std::cout << "1-2: " << sum12 / n12 << "\n";
    if (n23) std::cout << "2-3: " << sum23 / n23 << "\n";
    if (n13) std::cout << "1-3: " << sum13 / n13 << "\n";
}
