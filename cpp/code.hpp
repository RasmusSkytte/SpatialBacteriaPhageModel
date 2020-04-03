#ifndef COLONIES3DDEF
#define COLONIES3DDEF

#define ARMA_NO_DEBUG

#include <iostream>         // Input and output
#include <iomanip>          // Input and output formatting
#include <fstream>          // File streams

#include <armadillo>        // Matrix library

#include <random>           // Random numbers
#include <math.h>           // Mathmatical constants
#include <algorithm>        // Mathmatical constants

#include <vector>           // C++ standard vector
#include <string.h>         // Strings

#include <cassert>          // Assertions

#include <sys/types.h>      // Packages for the directory
#include <sys/stat.h>       //    information,
#include <dirent.h>         //    handling and etc.
#include <unistd.h>
#include <ctime>            // Time functions

/* Class to contain simulation parameters, and which drives the simulation forward */
class Colonies3D {
 private:

    double B_0;                     // [CFU/mL] Initial density of bacteria
    double P_0;                     // [PFU/mL] Initial density of phages

    double K;                       //          Carrying capacity
    double n_0;                     // [1/ml]   Initial nutrient level (number of bacteria per ml)

    double L;                       // [µm]     Side-length of simulation array
    double H;                       // [µm]     Height of the simulation array
    int    nGridXY;                 //          Number of gridpoints
    int    nGridZ;                  //          Number of gridpoints

    double nSamp;                   //          Number of samples to save per simulation hour

    double g;                       // [1/hour] Growth rate for the cells

    double alpha;                   //          Percentage of phages which escape the colony upon lysis
    int    beta;                    //          Multiplication factor phage
    double eta;                     //          Adsorption coefficient
    double delta;                   // [1/hour] Rate of phage decay
    double r;                       //          Constant used in the time-delay mechanism
    double zeta;                    //          permeability of colony surface

    double D_B;                     // [µm^2/hour] Diffusion constant for the cells
    double D_P;                     // [µm^2/hour] Diffusion constant for the phage
    double D_n;                     // [µm^2/hour] Diffusion constant for the nutrient

    double lambdaB;                 // Probability of cell to jump to neighbour point
    double lambdaP;                 // Probability of phage to jump to neighbour point

    double T;                       // [hours]  Current time
    double dT;                      // [hours]  Time-step size
    double T_end;                   // [hours]  End time of simulation
    double T_i;                     // [hours]  Time when the phage infections begins (less than 0 disables phage infection)

    double initialOccupancy;        // Number of gridpoints occupied initially;

    bool   exit;                    // Boolean to control early exit

    bool   Warn_g;                  //
    bool   Warn_r;                  //
    bool   Warn_eta;                // Booleans to keep track of warnings
    bool   Warn_delta;              //
    bool   Warn_density;            //
    bool   Warn_fastGrowth;         //

    bool   experimentalConditions;  // Booleans to control simulation type

    bool   clustering;              // When false, the ((B+I)/nC)^(1/3) factor is removed.
    bool   shielding;               // When true the simulation uses the shielding function (full model)
    bool   reducedBeta;             // When true the simulation modifies the burst size by the growthfactor

    bool   reducedBoundary;         // When true, bacteria are spawned at X = 0 and Y = 0. And phages are only spawned within nGrid/s boxes from (0,0,z).
    int    s;

    int    timeScaleSeperation;     // Indicates the difference in time scale between diffusion of nutrient

    bool   fastExit;                // Stop simulation when all cells are dead

    bool   exportAll;               // Boolean to export everything, not just populationsize

    double rngSeed;                 // The seed for the random number generator
    std::mt19937 rng;               // Mersenne twister, random number generator
    std::uniform_real_distribution  <double> rand;
    std::normal_distribution        <double> randn;

    std::ofstream f_B;              // Filestream to save configuration of sucebtible cells
    std::ofstream f_I;              // Filestream to save configuration of infected cells
    std::ofstream f_P;              // Filestream to save configuration of phages
    std::ofstream f_n;              // Filestream to save configuration of nutrient
    std::ofstream f_N;              // Filestream to save number agents
    std::ofstream f_log;            // Filestream to save log.txt

    std::string path;               // Sets the path to store in

    // Coordinates of agents in the simulation
    arma::Cube<double> B;           // Sensitive bacteria
    arma::Cube<double> P;           // Phages
    arma::Cube<double> I0;          // Infected bacteria
    arma::Cube<double> I1;          // Infected bacteria
    arma::Cube<double> I2;          // Infected bacteria
    arma::Cube<double> I3;          // Infected bacteria
    arma::Cube<double> I4;          // Infected bacteria
    arma::Cube<double> I5;          // Infected bacteria
    arma::Cube<double> I6;          // Infected bacteria
    arma::Cube<double> I7;          // Infected bacteria
    arma::Cube<double> I8;          // Infected bacteria
    arma::Cube<double> I9;          // Infected bacteria
    arma::Cube<double> nC;          // Number of colonies in gridpoint

    arma::Cube<double> B_new;       // Sensitive bacteria
    arma::Cube<double> P_new;       // Phages
    arma::Cube<double> I0_new;      // Infected bacteria
    arma::Cube<double> I1_new;      // Infected bacteria
    arma::Cube<double> I2_new;      // Infected bacteria
    arma::Cube<double> I3_new;      // Infected bacteria
    arma::Cube<double> I4_new;      // Infected bacteria
    arma::Cube<double> I5_new;      // Infected bacteria
    arma::Cube<double> I6_new;      // Infected bacteria
    arma::Cube<double> I7_new;      // Infected bacteria
    arma::Cube<double> I8_new;      // Infected bacteria
    arma::Cube<double> I9_new;      // Infected bacteria

    // Nutrient matrix
    arma::cube   nutrient;
    arma::cube   nutrient_new;

    // Occupancy of grid
    arma::Cube<double> Occ;

    // Laplacian for nutrient diffusion
    double alphaXY;
    double alphaZ;

 public:
    // Constructers
    explicit    Colonies3D(double B_0, double P_0);                           // Direct constructer

    // Driver
    int         Run(double T_end);                                              // Controls the evaluation of the simulation

 private:
    void        Initialize();                                                   // Initialize the simulation
    void        spawnBacteria();                                                // Spawns the bacteria
    void        spawnPhages();                                                  // Spawns the phages
    void        ComputeTimeStep();                                              // Computes the size of the time-step needed
    double      ComputeEvents(double n, double p, int flag);                    // Returns the number of events ocurring for given n and p
    void        ComputeDiffusion(double n, double lambda,                       // Computes how many particles has moved to neighbouing points
                    double* n_0, double* n_u, double* n_d, double* n_l, double* n_r, double* n_f, double* n_b, int flag);

 public:
    void        SetLength(double L);                                            // Set the side-length of the simulation
    void        SetHeight(double H);                                            // Set the height of the simulation
    void        SetGridSize(double nGrid);                                      // Set the number of gridpoints
    void        SetTimeStep(double dT);                                         // Set the time step size
    void        SetSamples(int nSamp);                                          // Set the number of output samples

    void        PhageInvasionStartTime(double T_i);                             // Sets the time when the phages should start infecting

    void        CellGrowthRate(double g);                                       // Sets the maximum growthrate
    void        CellCarryingCapacity(double K);                                 // Sets the carrying capacity
    void        CellDiffusionConstant(double D_B);                              // Sets the diffusion constant of the phages

    void        PhageBurstSize(int beta);                                       // Sets the size of the bursts
    void        PhageAdsorptionRate(double eta);                                // sets the adsorption parameter gamma
    void        PhageDecayRate(double delta);                                   // Sets the decay rate of the phages
    void        PhageInfectionRate(double r);                                   // Sets rate of the infection increaasing in stage
    void        PhageDiffusionConstant(double D_P);                             // Sets the diffusion constant of the phages
    void        PhageLatencyTime(double tau);                                   // Sets latency time of the phage (r and tau are related by r = 10 / tau)

    void        SurfacePermeability(double zeta);                               // Sets the permeability of the surface

    void        InitialNutrient(double n_0);                                    // Sets the amount of initial nutrient
    void        NutrientDiffusionConstant(double D_n);                          // Sets the nutrient diffusion rate

    void        SimulateExperimentalConditions();                               // Sets the simulation to spawn phages at top layer and only have x-y periodic boundaries

    void        DisableShielding();                                             // Sets shielding bool to false
    void        DisablesClustering();                                           // Sets clustering bool to false
    void        ReducedBurstSize();                                             // Sets the simulation to limit beta as n -> 0

    void        ReducedBoundary(int s);                                         // Sets the reduced boundary bool to true and the value of s

    void        SetAlpha(double alpha);                                         // Sets the value of alpha

 private:
    // Helping functions
    int         RandI(int n);                                                   // Returns random integer less than n
    double      Rand(double n);                                                 // Returns random double less than n
    double      RandN(double m, double s);                                      // Returns random normal dist. number with mean m and variance s^2
    double      RandP(double l);                                                // Returns poisson dist. number with mean l

 public:
    void        SetRngSeed(int n);                                              // Sets the seed of the random number generator

 private:
    void        WriteLog();                                                     // Write a log.txt file

 public:
    void        FastExit();                                                     // Stop simulation when all cells are dead
    void        ExportAll();                                                    // Sets the simulation to export everything

 private:
    void        ExportData(double t);                                           // Master function to export the data

    // Data handling
    void        OpenFileStream(std::ofstream& stream,                           // Open filstream if not allready opened
                    std::string& fileName);
    std::string GeneratePath();                                                 // Generates a save path for datafiles

 public:
    void        SetFolderNumber(int number);                                    // Sets the folder number (useful when running parralel code)
    void        SetPath(std::string& path);                                     // Sets the folder path (useful when running parralel code)

    // Get properties
    std::string GetPath();                                                      // Returns the save path
    int         GetTime();                                                      // Returns the time
    int         GetDeltaT();                                                    // Returns the time-step dT


    // Clean up
    void        DeleteFolder();                                                 // Delete the data folder
 private:
    void        DeleteFolderTree(const char* directory_name);                   // Delete folders recursively

 public:
    ~Colonies3D();                                                            // Destructor

};

#endif

