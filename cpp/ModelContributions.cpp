#include "code.hpp"
#include "numThreads.hpp"
#include <omp.h>

using namespace std;

int main(int argc, char** argv){

    // Set the type of output
    int model = 4;           // Lotka Volterra with time delay, space and colony level protection

    // Default parameters
    int nGrid 	 = 50;
    int beta  	 = 100;
    double delta = 0.1;
    double eta   = 1e4;
    double tau   = 0.5;

    for (int i = 0; i < argc; i++) {
        if (i == 1) {
            model = atoi(argv[i]);
        }
    }

    // Default parameters

    // Powers of 10 and half-powers of 10
    vector<double> B_0 = {pow(10, 0.00), pow(10, 0.50), pow(10, 1.00), pow(10, 1.50), pow(10, 2.00), pow(10, 2.50), pow(10, 3.00), pow(10, 3.50), pow(10, 4.00), pow(10, 4.50), pow(10, 5.00), pow(10, 5.50), pow(10, 6.00), pow(10, 6.50), pow(10, 7.00), pow(10, 7.50), pow(10, 8.00), pow(10, 8.50), pow(10, 9.00)};
    vector<double> P_0 = {pow(10, 0.00), pow(10, 0.50), pow(10, 1.00), pow(10, 1.50), pow(10, 2.00), pow(10, 2.50), pow(10, 3.00), pow(10, 3.50), pow(10, 4.00), pow(10, 4.50), pow(10, 5.00), pow(10, 5.50), pow(10, 6.00), pow(10, 6.50), pow(10, 7.00), pow(10, 7.50), pow(10, 8.00), pow(10, 8.50), pow(10, 9.00)};

    double T = 15; 	// Simulation length

    omp_set_dynamic(0);
    #pragma omp parallel for collapse(3) schedule(dynamic) num_threads(NUM_THREADS)
    for (int k = 0; k < 6; k++) {
        for (int j = P_0.size()-1; j >= 0; j--) {
            for (int i = 0; i < B_0.size(); i++) {

                Colonies3D s(B_0[i],P_0[j]);

                char buffer[80];                                  // Create a buffer to store the strings

                string pathName = "";
                pathName += "Model_";
                pathName += to_string(model);

                if (k == 0) {
                    pathName += "_onlyClustering";
                    s.SetAlpha(0);
                    s.DisableShielding();
                } else if (k == 1) {
                    pathName += "_onlyShielding";
                    s.SetAlpha(0);
                    s.DisablesClustering();
                } else if (k == 2) {
                    pathName += "_onlyReadsorption";
                    s.DisablesClustering();
                    s.DisableShielding();
                } else if (k == 3) {
                    pathName += "_noClustering";
                    s.DisablesClustering();
                } else if (k == 4) {
                    pathName += "_noShielding";
                    s.DisableShielding();
                } else if (k == 5) {
                    pathName += "_noReadsorption";
                    s.SetAlpha(0);
                } else if (k == 6) {
                    pathName += "_ReducedBeta";
                    s.ReducedBurstSize();
                }

                pathName += "/beta_";
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

                pathName += "_nGrid_";
                pathName += to_string(nGrid);

                pathName += "/CFU_1e";
                snprintf(buffer, sizeof(buffer), "%.2f", log10(B_0[i]));
                pathName += string(buffer);

                pathName += "_PFU_1e";
                snprintf(buffer, sizeof(buffer), "%.2f", log10(P_0[j]));
                pathName += string(buffer);

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

                s.SetRngSeed(5);

                // Set the static parameters
                s.SetGridSize(nGrid);
                s.PhageBurstSize(beta);
                s.PhageDecayRate(delta);
                s.PhageAdsorptionRate(eta);
                s.PhageLatencyTime(tau);

                s.SetSamples(100);

                s.FastExit();

                // Start the simulation
                s.Run(T);
            }
        }
    }

    return 0;
}
