#include "code.hpp"
#include "numThreads.hpp"
#include <omp.h>

using namespace std;

int main(int argc, char** argv){

    // Set the type of output
    // int model = 1;           // Simple Lotka Volterra
    // int model = 2;           // Lotka Volterra with time delay
    // int model = 3;           // Lotka Volterra with time delay and space
    int model = 4;           // Lotka Volterra with time delay, space and colony level protection

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

    double B_0 = pow(10, 4.0);
    double P_0 = pow(10, 5.0);

    double T = 20; 	// Simulation length

    Colonies3D s(B_0,P_0);

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
    snprintf(buffer, sizeof(buffer), "%.2f", log10(B_0));
    pathName += string(buffer);

    pathName += "_PFU_1e";
    snprintf(buffer, sizeof(buffer), "%.2f", log10(P_0));
    pathName += string(buffer);

    // Check if data is already made
    string path_s = "data/"; // Data folder name
    path_s += pathName;
    path_s += "/Completed.txt";

    // Check if run exists and is completed
    struct stat info;
    if (stat(path_s.c_str(), &info) == 0 && S_ISREG(info.st_mode)) {
        return 0;
    }

    s.SetPath(pathName);

    s.SetRngSeed(5);

    // Set the static parameters
    s.PhageBurstSize(beta);
    s.PhageDecayRate(delta);
    s.PhageAdsorptionRate(eta);
    s.SetSamples(100);

    // Model 1
    if (model == 1) {               // Simple Lotka Volterra

        // Set the size of the space
        s.SetGridSize(1);

        // Set the latency
        s.PhageLatencyTime(0.0);

        // Set the time step
        s.SetTimeStep(2e-3);

        // Disable colony level mechanisms
        s.SetAlpha(0.0);
        s.DisablesClustering();
        s.DisableShielding();

    } else if (model == 2) {        // Lotka Volterra with time delay

        // Set the size of the space
        s.SetGridSize(1);

        // Set the latency
        s.PhageLatencyTime(tau);

        // Set the time step
        s.SetTimeStep(2e-3);

        // Disable colony level mechanisms
        s.SetAlpha(0.0);
        s.DisablesClustering();
        s.DisableShielding();

    } else if (model == 3) {        // Lotka Volterra with time delay and space

        // Set the size of the space
        s.SetGridSize(nGrid);

        // Set the latency
        s.PhageLatencyTime(tau);

        // Disable colony level mechanisms
        s.SetAlpha(0.0);
        s.DisablesClustering();
        s.DisableShielding();

    } else if (model == 4) {        // Lotka Volterra with time delay, space and colony level protection

        // Set the size of the space
        s.SetGridSize(nGrid);

        // Set the latency
        s.PhageLatencyTime(tau);

    }

    // Start the simulation
    s.Run(T);

    return 0;
}
