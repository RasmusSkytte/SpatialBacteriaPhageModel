#include "code.hpp"
#include "numThreads.hpp"
#include <omp.h>

using namespace std;

int main(int argc, char** argv){

    // Default parameters
    vector<double> T_i = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};

    double P_0 = 6*pow(10,-1); // PFU / µm^2

    double T = 20; 	// Simulation length

    // Default parameters
    int nGrid    = 100;
    int beta     = 400;
    double delta = 0.003;
    double eta   = 1.32*1e4;
    double tau   = 1.0;
    double g_max = 60.0/31.0;
    double L     = 1e4;
    double H     = 4e2;

    int boundaryModifier = 10;
    omp_set_dynamic(0);
    #pragma omp parallel for collapse(4) schedule(dynamic) num_threads(NUM_THREADS)
    for (int j = -1; j < 75; j++) {
        for (int k = -2; k < 3; k++) {
            for (int i = 0; i < T_i.size(); i++) {
                for (int l = 1; l < 3; l++) {

                    double D_P;
                    if (l == 1) {
                        D_P   = 3e3;
                    }  else if (l == 2) {
                        D_P   = 1e4;
                    }

                    double P = P_0 / H * 1e12;

                    Colonies3D s(1e12 / (L * L * H), P );

                    s.SetRngSeed(j);

                    s.SetLength(L);
                    s.SetHeight(H);
                    s.SetGridSize(nGrid);
                    s.PhageBurstSize(beta);
                    s.PhageDecayRate(delta);
                    s.PhageAdsorptionRate(eta);
                    s.PhageLatencyTime(tau);
                    s.PhageDiffusionConstant(D_P);
                    s.SetAlpha(0.5);

                    s.PhageInvasionStartTime(T_i[i]);
                    s.CellGrowthRate(g_max);
                    s.SetSamples(100);

                    s.SimulateExperimentalConditions();
                    s.ReducedBoundary(boundaryModifier);


                    char buffer[80];                                  // Create a buffer to store the strings

                    string pathName = "";
                    pathName += "Experiment_Conditions/beta_";
                    pathName += to_string(beta);

                    pathName += "_delta_";
                    snprintf(buffer, sizeof(buffer), "%.3f", delta);
                    pathName += string(buffer);

                    pathName += "_eta_";
                    snprintf(buffer, sizeof(buffer), "%.0f", eta);
                    pathName += string(buffer);

                    pathName += "_tau_";
                    snprintf(buffer, sizeof(buffer), "%.2f", tau);
                    pathName += string(buffer);

                    pathName += "_D_P_";
                    snprintf(buffer, sizeof(buffer), "%.0f", D_P);
                    pathName += string(buffer);

                    pathName += "_nGrid_";
                    pathName += to_string(nGrid);

                    pathName += "_boundaryModifier_";
                    pathName += to_string(boundaryModifier);

                    if (k == -2) {
                        pathName += "_Model_3";
                        s.DisableShielding();
                        s.DisablesClustering();
                        s.SetAlpha(0.0);
                    } else if (k == -1) {
                        pathName += "_full_noSheilding";
                        s.DisableShielding();
                    } else if (k == 0) {
                        pathName += "_full_zeta_0.1";
                        s.SurfacePermeability(0.1);
                    } else if (k == 1) {
                        pathName += "_full_zeta_0.2";
                        s.SurfacePermeability(0.2);
                    } else if (k == 2) {
                        pathName += "_full_zeta_0.4";
                        s.SurfacePermeability(0.4);
                    }

                    pathName += "/T_i_";
                    snprintf(buffer, sizeof(buffer), "%.1f", T_i[i]);
                    pathName += string(buffer);

                    pathName += "/";
                    pathName += to_string(j);

                    // Check if data is already made
                    string path_s = "data/"; // Data folder name
                    path_s += pathName;
                    path_s += "/Completed.txt";

                    // Check if run exists and is completed
                    struct stat info;
                    if (stat(path_s.c_str(), &info) == 0 && S_ISREG(info.st_mode)) {
                        continue;
                    }

                    s.SetPath(pathName);

                    if (j >= 0) {
                        s.FastExit();
                    } else {
                        s.SetSamples(2);
                        s.ExportAll();
                    }

                    // Start the simulation
                    s.Run(T+T_i[i]);
                }
            }
        }
    }
    return 0;
}
