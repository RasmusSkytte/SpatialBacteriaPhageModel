#include "code.hpp"
#include "numThreads.hpp"
#include <omp.h>

using namespace std;

int main(int argc, char** argv){

    // Set the type of output
    int model = 1;

    // Default parameters
    int nGrid 	 = 50;
    int beta  	 = 100;
    double delta = 1.0/10.0;
    double eta   = 1e4;
    double tau   = 0.5;

    for (int i = 0; i < argc; i++) {
        if (i == 1) {
            model = atoi(argv[i]);
        }
    }

    // Powers of 10 and half-powers of 10
    vector<double> B_0 = {pow(10, 0.00), pow(10, 0.50), pow(10, 1.00), pow(10, 1.50), pow(10, 2.00), pow(10, 2.50), pow(10, 3.00), pow(10, 3.50), pow(10, 4.00), pow(10, 4.50), pow(10, 5.00), pow(10, 5.50), pow(10, 6.00), pow(10, 6.50), pow(10, 7.00), pow(10, 7.50), pow(10, 8.00), pow(10, 8.50), pow(10, 9.00)};
    vector<double> P_0 = {pow(10, 0.00), pow(10, 0.50), pow(10, 1.00), pow(10, 1.50), pow(10, 2.00), pow(10, 2.50), pow(10, 3.00), pow(10, 3.50), pow(10, 4.00), pow(10, 4.50), pow(10, 5.00), pow(10, 5.50), pow(10, 6.00), pow(10, 6.50), pow(10, 7.00), pow(10, 7.50), pow(10, 8.00), pow(10, 8.50), pow(10, 9.00)};

    double T = 20; 	// Simulation length


    // Set tau to zero for some models
    if ((model == 1) or (model == 8) or (model == 9))  tau = 0.0;


    #pragma omp parallel for collapse(2) schedule(dynamic) num_threads(NUM_THREADS)
    for (int j = P_0.size()-1; j >= 0; j--) {
        for (int i = B_0.size()-1; i >= 0; i--) {

            Colonies3D s(B_0[i],P_0[j]);

            char buffer[80];                                  // Create a buffer to store the strings
            string pathName = "";
            pathName += "Model_";
            pathName += to_string(model);

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

            if (model > 2) {
                pathName += "_nGrid_";
                pathName += to_string(nGrid);
            }

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
            s.PhageBurstSize(beta);
            s.PhageDecayRate(delta);
            s.PhageAdsorptionRate(eta);
            s.FastExit();
            s.SetSamples(100);

            // Set the latency
            s.PhageLatencyTime(tau);

            // Model 1
            if (model == 1) {               // Simple model

                // Set the size of the space
                s.SetGridSize(1);

                // Set the time step
                s.SetTimeStep(2e-3);

                // Disable colony level mechanisms
                s.SetAlpha(0.0);
                s.DisablesClustering();
                s.DisableShielding();

            } else if (model == 2) {        // Model with time delay

                // Set the size of the space
                s.SetGridSize(1);

                // Set the time step
                s.SetTimeStep(2e-3);

                // Disable colony level mechanisms
                s.SetAlpha(0.0);
                s.DisablesClustering();
                s.DisableShielding();

            } else if (model == 3) {        // Model with time delay and space

                // Set the size of the space
                s.SetGridSize(nGrid);

                // Disable colony level mechanisms
                s.SetAlpha(0.0);
                s.DisablesClustering();
                s.DisableShielding();

            } else if (model == 4) {        // Model with time delay, space and colony level protection

                // Set the size of the space
                s.SetGridSize(nGrid);

            } else if (model == 5) {        // Model with time delay, space and colony level protection (weak shielding)

                // Set the size of the space
                s.SetGridSize(nGrid);

                // Set the zeta parameter
                s.SurfacePermeability(0.2);


            } else if (model == 6) {        // Model with time delay, space and colony level protection (strong readsoprtion)

                // Set the size of the space
                s.SetGridSize(nGrid);

                // Set the alpha value
                s.SetAlpha(0.95);

            } else if (model == 7) {        // Model with time delay, space and colony level protection (complete readsorption)

                // Set the size of the space
                s.SetGridSize(nGrid);

                // Set the alpha value
                s.SetAlpha(1.0);

            } else if (model == 8) {        // Model with space

                // Set the size of the space
                s.SetGridSize(nGrid);

                // Disable colony level mechanisms
                s.SetAlpha(0.0);
                s.DisablesClustering();
                s.DisableShielding();

            } else if (model == 9) {        // Model with space and colony level protection

                // Set the size of the space
                s.SetGridSize(nGrid);

            }

            // Start the simulation
            s.Run(T);

        }
    }

    return 0;
}
