#include "code.hpp"

using namespace arma;
using namespace std;

// Constructers /////////////////////////////////////////////////////////////////////////
// Direct constructer
Colonies3D::Colonies3D(double B_0, double P_0){

    // Store the initial densities
    this->B_0 = B_0;
    this->P_0 = P_0;

    // Set some default parameters (initlize some default objects)
    K                       = 1.0 / 5.0;//          Half-Speed constant
    n_0                     = 1e9;      // [1/ml]   Initial nutrient level (Carrying capacity per ml)

    L                       = 1e4;      // [µm]     Side-length of simulation array
    H                       = L;        // [µm]     Height of the simulation array
    nGridXY                 = 50;       //          Number of gridpoints
    nGridZ                  = nGridXY;  //          Number of gridpoints

    nSamp                   = 10;       //          Number of samples to save per simulation hour

    g                       = 2;        // [1/h]    Doubling rate for the cells

    alpha                   = 0.5;      //          Percentage of phages which reinfect the colony upon lysis
    beta                    = 100;      //          Multiplication factor phage
    eta                     = 1e4;      // [µm^3/h] Adsorption coefficient
    delta                   = 1.0/10.0; // [1/h]    Rate of phage decay
    r                       = 10.0/0.5; //          Constant used in the time-delay mechanism
    zeta                    = 1.0;      //          permeability of colony surface

    D_P                     = 1e4;      // [µm^2/h] Diffusion constant for the phage
    D_B                     = 0;//D_P/20;   // [µm^2/h] Diffusion constant for the cells
    D_n                     = 25e5;     // [µm^2/h] Diffusion constant for the nutrient

    T                       = 0;        // [h]      Current time
    dT                      = -1;       // [h]      Time-step size (-1 to compute based on fastest diffusion rate)
    T_end                   = 0;        // [h]      End time of simulation
    T_i                     = -1;       // [h]      Time when the phage infections begins (less than 0 disables phage infection)

    initialOccupancy        = 0;        // Number of gridpoints occupied initially;

    exit                    = false;    // Boolean to control early exit

    Warn_g                  = false;    //
    Warn_r                  = false;    //
    Warn_eta                = false;    // Booleans to keep track of warnings
    Warn_delta              = false;    //
    Warn_density            = false;    //
    Warn_fastGrowth         = false;    //

    experimentalConditions  = false;    // Booleans to control simulation type

    clustering              = true;     // When false, the ((B+I)/nC)^(1/3) factor is removed.
    shielding               = true;     // When true the simulation uses the shielding function (full model)
    reducedBeta             = false;    // When true the simulation modifies the burst size by the growthfactor

    singleInternalState     = false;     // Boolean to toggle how many internal infected states are used

    reducedBoundary         = false;    // When true, bacteria are spawned at X = 0 and Y = 0. And phages are only spawned within nGrid boxes from (0,0,z).
    s                       = 1;

    fastExit                = false;     // Stop simulation when all cells are dead

    exportAll               = false;    // Boolean to export everything, not just populationsize

    rngSeed                 = -1;       // Random number seed  ( set to -1 if unused )

};


// Controls the evaluation of the simulation
int Colonies3D::Run(double T_end) {

    this->T_end = T_end;

    // Get start time
    time_t  tic;
    time(&tic);

    // Generate a path
    path = GeneratePath();

    // Initialize the simulation matrices
    Initialize();

    // Export data
    ExportData(T);

    // Determine the number of samples to take
    int nSamplings = nSamp*T_end;

    // Loop over samplings
    for (int n = 0; n < nSamplings; n++) {
        if (exit) break;

        // Determine the number of timesteps between sampings
        int nStepsPerSample = static_cast<int>(round(1 / (nSamp *  dT)));

        for (int t = 0; t < nStepsPerSample; t++) {
            if (exit) break;

            // Increase time
            T += dT;

            // Spawn phages
            if ((T_i >= 0) and (abs(T - T_i) < dT / 2)) {
                spawnPhages();
                T_i = -1;
            }

            // Reset density counter
            double maxOccupancy = 0.0;

	        //#pragma omp parallel for collapse(3) schedule(dynamic)
            for (uword k = 0; k < nGridZ; k++ ) {
                for (uword j = 0; j < nGridXY; j++ ) {
                    for (uword i = 0; i < nGridXY; i++) {

                        // Copy values from memory only once
                        double tmpOcc       = Occ(i, j, k);
                        double tmpNutrient  = nutrient(i, j, k);
                        double tmpB         = B(i, j, k);
                        double tmpI0        = I0(i, j, k);
                        double tmpI1        = I1(i, j, k);
                        double tmpI2        = I2(i, j, k);
                        double tmpI3        = I3(i, j, k);
                        double tmpI4        = I4(i, j, k);
                        double tmpI5        = I5(i, j, k);
                        double tmpI6        = I6(i, j, k);
                        double tmpI7        = I7(i, j, k);
                        double tmpI8        = I8(i, j, k);
                        double tmpI9        = I9(i, j, k);
                        double tmpP         = P(i, j, k);
                        double tmpNC        = nC(i, j, k);

                        // Skip empty sites
                        if ((tmpOcc < 1) and (tmpP < 1)) continue;

                        // Record the maximum observed density
                        if (tmpOcc > maxOccupancy) maxOccupancy = tmpOcc;

                        // Compute the growth modifier
                        double growthModifier = tmpNutrient / (tmpNutrient + K);

                        // Compute beta
                        double Beta = beta;
                        if (reducedBeta) {
                            Beta *= growthModifier;
                        }

                        double p = 0;
                        double N = 0;
                        double M = 0;

                        // Birth //////////////////////////////////////////////////////////////////////
                        p = g*growthModifier*dT;
                        if (tmpNutrient < 1) p = 0;
                        N = ComputeEvents(tmpB, p, 1);


                        // Ensure there is enough nutrient
                        if ( N > tmpNutrient ) {
                            N = round( tmpNutrient );
                        }

                        // Update count
                        B_new(i, j, k) += N;
                        tmpNutrient = max(0.0, tmpNutrient - N);

                        // Increase Infections ////////////////////////////////////////////////////////
                        if (r > 0.0) {

                            p = r*growthModifier*dT;
                            if (tmpNutrient < 1) p = 0;

                            N = ComputeEvents(tmpI9, p, 2);  // Bursting events

                            if (N > tmpI9) N = tmpI9;        // If more bacteria than present are set to burst, round down

                            // Update count
                            tmpI9    = max(0.0, tmpI9  - N);
                            tmpOcc   = max(0.0, tmpOcc - N);
                            P_new(i, j, k) += round( (1 - alpha) * Beta * N);   // Phages which escape the colony
                            M = round(alpha * Beta * N);                        // Phages which reinfect the colony

                            // Non-bursting events
                            if (not singleInternalState) {
                                N = ComputeEvents(tmpI8, p, 2); tmpI8 = max(0.0, tmpI8 - N); tmpI9 += N;
                                N = ComputeEvents(tmpI7, p, 2); tmpI7 = max(0.0, tmpI7 - N); tmpI8 += N;
                                N = ComputeEvents(tmpI6, p, 2); tmpI6 = max(0.0, tmpI6 - N); tmpI7 += N;
                                N = ComputeEvents(tmpI5, p, 2); tmpI5 = max(0.0, tmpI5 - N); tmpI6 += N;
                                N = ComputeEvents(tmpI4, p, 2); tmpI4 = max(0.0, tmpI4 - N); tmpI5 += N;
                                N = ComputeEvents(tmpI3, p, 2); tmpI3 = max(0.0, tmpI3 - N); tmpI4 += N;
                                N = ComputeEvents(tmpI2, p, 2); tmpI2 = max(0.0, tmpI2 - N); tmpI3 += N;
                                N = ComputeEvents(tmpI1, p, 2); tmpI1 = max(0.0, tmpI1 - N); tmpI2 += N;
                                N = ComputeEvents(tmpI0, p, 2); tmpI0 = max(0.0, tmpI0 - N); tmpI1 += N;
                            }
                        }

                        // Infections /////////////////////////////////////////////////////////////////
                        if ((tmpOcc >= 1) and (tmpP >= 1)) {
                            double s;   // The factor which modifies the adsorption rate
                            double n;   // The number of targets the phage has

                            if (clustering) {   // Check if clustering is enabled
                                s = pow(tmpOcc / tmpNC, 1.0 / 3.0);
                                n = tmpNC;
                            } else {            // Else use mean field computation
                                s = 1.0;
                                n = tmpOcc;
                            }

                            // Compute the number of hits
                            if (eta * s * dT >= 1) { // In the diffusion limited case every phage hits a target
                                N = tmpP;
                            } else {
                                p = 1 - pow(1 - eta * s * dT, n);  // Probability hitting any target
                                N = ComputeEvents(tmpP, p, 4);     // Number of targets hit
                            }

                            // If bacteria were hit, update events
                            if ((N + M) >= 1) {

                                tmpP    = max(0.0, tmpP - N);     // Update count

                                double S;
                                if (shielding) {
                                    // Absorbing medium model
                                    double d = pow(tmpOcc / tmpNC, 1.0 / 3.0) - pow(tmpB / tmpNC, 1.0 / 3.0);
                                    S = exp(-zeta*d); // Probability of hitting susceptible target

                                } else {
                                    // Well mixed model
                                    S = tmpB / tmpOcc;
                                }

                                p = max(0.0, min(tmpB / tmpOcc, S)); // Probability of hitting susceptible target
                                N = ComputeEvents(N + M, p, 4);                  // Number of targets hit

                                if (N > tmpB) N = tmpB;              // If more bacteria than present are set to be infected, round down

                                // Update the counts
                                tmpB = max(0.0, tmpB - N);

                                if (r > 0.0) {
                                    if (not singleInternalState) {
                                        I0_new(i, j, k) += N;
                                    } else {
                                        I9_new(i, j, k) += N;
                                    }
                                } else {
                                    P_new(i, j, k) += N * (1 - alpha) * Beta;
                                }
                            }
                        }

                        // Phage Decay ////////////////////////////////////////////////////////////////
                        p = delta*dT;
                        N = ComputeEvents(tmpP, p, 5);
                        if (N > tmpP) N = tmpP;
                        tmpP = max(0.0, tmpP - N);

                        // Movement ///////////////////////////////////////////////////////////////////
                        if (nGridXY > 1) {

                            // Update positions
                            uword ip, jp, kp, im, jm, km;

                            if (i + 1 >= nGridXY/s) ip = i + 1 - nGridXY/s;
                            else ip = i + 1;

                            if (i == 0) im = nGridXY/s - 1;
                            else im = i - 1;

                            if (j + 1 >= nGridXY/s) jp = j + 1 - nGridXY/s;
                            else jp = j + 1;

                            if (j == 0) jm = nGridXY/s - 1;
                            else jm = j - 1;

                            if (not experimentalConditions) {   // Periodic boundaries in Z direction

                                if (k + 1 >= nGridZ) kp = k + 1 - nGridZ;
                                else kp = k + 1;

                                if (k == 0) km = nGridZ - 1;
                                else km = k - 1;

                            } else {    // Reflective boundaries in Z direction

                                if (k + 1 >= nGridZ) kp = k;
                                else kp = k + 1;

                                if (k == 0) km = k;
                                else km = k - 1;

                            }

                            // CELLS
                            B_new(i, j, k) += tmpB;
                            if (r > 0.0) {
                                I0_new(i, j, k) += tmpI0;
                                I1_new(i, j, k) += tmpI1;
                                I2_new(i, j, k) += tmpI2;
                                I3_new(i, j, k) += tmpI3;
                                I4_new(i, j, k) += tmpI4;
                                I5_new(i, j, k) += tmpI5;
                                I6_new(i, j, k) += tmpI6;
                                I7_new(i, j, k) += tmpI7;
                                I8_new(i, j, k) += tmpI8;
                                I9_new(i, j, k) += tmpI9;
                            }

                            // Update counts
                            double n_0; // No movement
                            double n_u; // Up
                            double n_d; // Down
                            double n_l; // Left
                            double n_r; // Right
                            double n_f; // Front
                            double n_b; // Back

                            // PHAGES
                            ComputeDiffusion(tmpP, lambdaP, &n_0, &n_u, &n_d, &n_l, &n_r, &n_f, &n_b, 3);
                            P_new(i, j, k)  += n_0;
                            P_new(ip, j, k) += n_u;
                            P_new(im, j, k) += n_d;
                            P_new(i, jp, k) += n_r;
                            P_new(i, jm, k) += n_l;
                            P_new(i, j, kp) += n_f;
                            P_new(i, j, km) += n_b;

			                // NUTRIENT
                            nutrient(i, j, k) = tmpNutrient;

                        } else {

                            // CELLS
                            B_new += tmpB;

                            if (r > 0.0) {
                                I0_new += tmpI0;
                                I1_new += tmpI1;
                                I2_new += tmpI2;
                                I3_new += tmpI3;
                                I4_new += tmpI4;
                                I5_new += tmpI5;
                                I6_new += tmpI6;
                                I7_new += tmpI7;
                                I8_new += tmpI8;
                                I9_new += tmpI9;
                            }

                            // PHAGES
                            P_new += tmpP;

			                // NUTRIENT
                            nutrient = tmpNutrient;

                        }
                    }
                }
            }

            //#pragma omp parallel for collapse(3) schedule(dynamic)
	        for (uword k = 0; k < nGridZ; k++ ) {
                for (uword j = 0; j < nGridXY; j++ ) {
                    for (uword i = 0; i < nGridXY; i++) {

                        if (nGridXY > 1) {

                            double tmpNutrient = nutrient(i, j, k);
                            double inflowNutrient = 0.0;
                            double outflowNutrient = 0.0;

                            // Update positions
                            uword ip, jp, kp, im, jm, km;

                            if (i + 1 >= nGridXY) ip = i + 1 - nGridXY;
                            else ip = i + 1;

                            if (i == 0) im = nGridXY - 1;
                            else im = i - 1;

                            if (j + 1 >= nGridXY) jp = j + 1 - nGridXY;
                            else jp = j + 1;

                            if (j == 0) jm = nGridXY - 1;
                            else jm = j - 1;

                            // XY Diffusion is always with periodic boundaries:
                            inflowNutrient  += alphaXY * nutrient(ip, j, k);
                            inflowNutrient  += alphaXY * nutrient(im, j, k);
                            inflowNutrient  += alphaXY * nutrient(i, jp, k);
                            inflowNutrient  += alphaXY * nutrient(i, jm, k);
                            outflowNutrient += 4 * alphaXY * tmpNutrient;

                            // Z direction needs special attention
                            if (not experimentalConditions) {   // Periodic boundaries in Z direction

                                if (k + 1 >= nGridZ) kp = k + 1 - nGridZ;
 				                else kp = k + 1;

                                if (k == 0) km = nGridZ - 1;
                                else km = k - 1;

                                inflowNutrient  += alphaZ  * nutrient(i, j, kp);
                                inflowNutrient  += alphaZ  * nutrient(i, j, km);
				                outflowNutrient += 2 * alphaZ * tmpNutrient;

                            } else {    // Reflective boundaries in Z direction

                                if (k + 1 < nGridZ) {
                                    // (i,j,k) is not an (upper) edge point
                                    kp = k + 1;
                                    inflowNutrient += alphaZ  * nutrient(i, j, kp);
                                    outflowNutrient += alphaZ  * tmpNutrient;
                                }

                                if (k > 0) {
                                    // (i,j,k) is not an (lower) edge point
                                    km = k - 1;
                                    inflowNutrient += alphaZ  * nutrient(i, j, km);
                                    outflowNutrient += alphaZ  * tmpNutrient;
                                }
                            }

			                // Inflow minus outflow
                            nutrient_new(i, j, k) = tmpNutrient + inflowNutrient - outflowNutrient;


                        } else {
                            // Nutrient
                            nutrient_new = nutrient;
                        }
                    }
                }
            }

            // Update arrays
            B.swap(B_new);                  B_new.zeros();
            I0.swap(I0_new);                I0_new.zeros();
            I1.swap(I1_new);                I1_new.zeros();
            I2.swap(I2_new);                I2_new.zeros();
            I3.swap(I3_new);                I3_new.zeros();
            I4.swap(I4_new);                I4_new.zeros();
            I5.swap(I5_new);                I5_new.zeros();
            I6.swap(I6_new);                I6_new.zeros();
            I7.swap(I7_new);                I7_new.zeros();
            I8.swap(I8_new);                I8_new.zeros();
            I9.swap(I9_new);                I9_new.zeros();
            P.swap(P_new);                  P_new.zeros();
            nutrient.swap(nutrient_new);    nutrient_new.zeros();

            // Update occupancy
            Occ = B + I0 + I1 + I2 + I3 + I4 + I5 + I6 + I7 + I8 + I9;

            // Update nC
            uvec I = find(Occ < nC);
            nC(I) = Occ(I);

            if ((maxOccupancy > L * L * H / (nGridXY * nGridXY * nGridZ)) and (!Warn_density)) {
                cout << "\tWarning: Maximum Density Large!" << "\n";
                f_log  << "Warning: Maximum Density Large!" << "\n";
                Warn_density = true;
            }
        }

        // 2) There are no more alive cells
        // -> Stop simulation

        if ((fastExit) and (accu(Occ) < 1)) {
            exit = true;
        }

        // 3) The food is on average less than one per gridpoint
        // and the maximal nutrient at any point in space is less than 1

        if (fastExit) {
            if  ((accu(nutrient) < nGridZ*pow(nGridXY,2)) and (nutrient.max() < 0.5)) {
                exit = true;
            }
        }

        // Store the state
        ExportData(T);

        // Check for nutrient stability
        assert(accu(nutrient) >= 0);
        assert(accu(nutrient) <= n_0 * L * L * H);
    }

    // Get stop time
    time_t  toc;
    time(&toc);

    // Calculate time difference
    float seconds = difftime(toc, tic);
    float hours   = floor(seconds/3600);
    float minutes = floor(seconds/60);
    minutes -= hours*60;
    seconds -= minutes*60 + hours*3600;

    cout << "\n";
    cout << "\tSimulation complete after ";
    if (hours > 0.0)   cout << hours   << " hours and ";
    if (minutes > 0.0) cout << minutes << " minutes and ";
    cout  << seconds << " seconds." << "\n";

    std::ofstream f_out;
    f_out.open(GetPath() + "/Completed.txt",fstream::trunc);
    f_out << "\tSimulation complete after ";
    if (hours > 0.0)   f_out << hours   << " hours and ";
    if (minutes > 0.0) f_out << minutes << " minutes and ";
    f_out  << seconds << " seconds." << "\n";
    f_out.flush();
    f_out.close();

    // Write sucess to log
    if (exit) {
        f_log << ">>Simulation completed with exit flag<<" << "\n";
    } else {
        f_log << ">>Simulation completed without exit flag<<" << "\n";
    }

    if (true) {
    std::ofstream f_timing;
    // cout << "\n";
    f_timing << "\t"       << setw(3) << difftime(toc, tic) << " s of total time" << "\n";

    f_timing.flush();
    f_timing.close();
    // cout << "\t----------------------------------------------------"<< "\n" << "\n" << "\n";
    }

    if (exit) {
        return 1;
    } else {
        return 0;
    }
}

// Initialize the simulation
void Colonies3D::Initialize() {

    // Set the random number generator seed
    if (rngSeed >= 0.0) {
        rng.seed( rngSeed );
    } else {
        static std::random_device rd;
        rng.seed(rd());
    }

    // Compute nGridZ
    if (L != H) {
        nGridZ = round(H / L * nGridXY);
        H = nGridZ * L / nGridXY;
    } else {
        nGridZ = nGridXY;
    }

    // Allocate the arrays
    B.zeros(nGridXY, nGridXY, nGridZ);
    I0.zeros(nGridXY, nGridXY, nGridZ);
    I1.zeros(nGridXY, nGridXY, nGridZ);
    I2.zeros(nGridXY, nGridXY, nGridZ);
    I3.zeros(nGridXY, nGridXY, nGridZ);
    I4.zeros(nGridXY, nGridXY, nGridZ);
    I5.zeros(nGridXY, nGridXY, nGridZ);
    I6.zeros(nGridXY, nGridXY, nGridZ);
    I7.zeros(nGridXY, nGridXY, nGridZ);
    I8.zeros(nGridXY, nGridXY, nGridZ);
    I9.zeros(nGridXY, nGridXY, nGridZ);
    P.zeros(nGridXY, nGridXY, nGridZ);
    nC.zeros(nGridXY, nGridXY, nGridZ);

    B_new.zeros(nGridXY, nGridXY, nGridZ);
    I0_new.zeros(nGridXY, nGridXY, nGridZ);
    I1_new.zeros(nGridXY, nGridXY, nGridZ);
    I2_new.zeros(nGridXY, nGridXY, nGridZ);
    I3_new.zeros(nGridXY, nGridXY, nGridZ);
    I4_new.zeros(nGridXY, nGridXY, nGridZ);
    I5_new.zeros(nGridXY, nGridXY, nGridZ);
    I6_new.zeros(nGridXY, nGridXY, nGridZ);
    I7_new.zeros(nGridXY, nGridXY, nGridZ);
    I8_new.zeros(nGridXY, nGridXY, nGridZ);
    I9_new.zeros(nGridXY, nGridXY, nGridZ);
    P_new.zeros(nGridXY, nGridXY, nGridZ);

    nutrient.zeros(nGridXY, nGridXY, nGridZ);
    nutrient_new.zeros(nGridXY, nGridXY, nGridZ);

    Occ.zeros(nGridXY, nGridXY, nGridZ);

    // Compute the step size
    double dXY = L / (double)nGridXY;
    double dZ  = H / (double)nGridZ;
    double dV  = dXY * dXY * dZ;

    // Initialize nutrient
    nutrient.fill( n_0 / 1e12 * dV ); // (n_0 / 1e12) is the nutrient per µm^Z

    // Compute the size of the time step
    ComputeTimeStep();

    // Compute the diffusion alphas
    alphaXY = D_n * dT / pow(dXY, 2);
    alphaZ  = D_n * dT / pow(dZ, 2);

    assert(2 * alphaXY <= 1);
    assert(2 * alphaZ  <= 1);

    // Store the parameters
    WriteLog();

    // Convert parameters to match gridpoint ///////
    // Adjust eta to match volume
    // eta/V is the number of collisions per hour for a single target
    eta = eta / dV;    // Number of collisions per gridpoint per hour

    // Adjust carrying capacity
    K = K * n_0 / 1e12 * dV;   // Determine the Monod growth factor used in n(i,j)/(n(i,j)+K)

    // Initialize the bacteria and phage populations
    spawnBacteria();
    if (T_i <= dT) {
        spawnPhages();
        T_i = -1;
    }
}

// Spawns the bacteria
void Colonies3D::spawnBacteria() {

    // Determine the number of cells to spawn
    double nBacteria = round(L * L * H * B_0 / 1e12);

    // Average bacteria per gridpoint
    double avgBacteria = nBacteria / (nGridXY * nGridXY * nGridZ);

    // Keep track of the number of cells spawned
    double numB = 0;

    // Initialize cell and phage populations
    if (nBacteria > (nGridXY * nGridXY * nGridZ)) {
        for (uword k = 0; k < nGridZ; k++) {
            for (uword j = 0; j < nGridXY; j++) {
                for (uword i = 0; i < nGridXY; i++) {

                    // Compute the number of bacteria to land in this gridpoint
                    double BB = RandP(avgBacteria);
                    if (BB < 1) continue;

                    // Store the number of clusters in this gridpoint
                    nC(i, j, k) = BB;

                    // Add the bacteria
                    B(i, j, k) = BB;
                    numB += BB;
                }
            }
        }
    }

    // Correct for underspawning
    while (numB < nBacteria) {

        // Choose random point in space
        uword i = RandI(nGridXY - 1);
        uword j = RandI(nGridXY - 1);
        uword k = RandI(nGridZ  - 1);

        if (reducedBoundary) {
            i = 0;
            j = 0;
        }

        // Add the bacteria
        B(i, j, k)++;
        numB++;
        nC(i, j, k)++;


    }

    // Correct for overspawning
    while (numB > nBacteria) {
        uword i = RandI(nGridXY - 1);
        uword j = RandI(nGridXY - 1);
        uword k = RandI(nGridZ  - 1);

        if (B(i, j, k) < 1) continue;

        B(i, j, k)--;
        numB--;
        nC(i, j, k)--;
    }

    // Count the initial occupancy
    uvec nz = find(B);
    initialOccupancy = nz.n_elem;

    // Determine the occupancy
    for (uword k = 0; k < nGridZ; k++ ) {
        for (uword j = 0; j < nGridXY; j++ ) {
            for (uword i = 0; i < nGridXY; i++) {
                Occ(i, j, k) = B(i, j, k) + I0(i, j, k) + I1(i, j, k) + I2(i, j, k) + I3(i, j, k) + I4(i, j, k) + I5(i, j, k) + I6(i, j, k) + I7(i, j, k) + I8(i, j, k) + I9(i, j, k);
            }
        }
    }

    // Write output
    if (nBacteria == 1) {

        uvec h = find(B);

        std::ofstream f_out;
        f_out.open(GetPath() + "/Depth.txt",fstream::trunc);
        f_out << "Bacteria initialized at Z = " << H / nGridZ * (0.5 + h(0) / (int)pow(nGridXY,2)) << endl;
        f_out.flush();
        f_out.close();
    }
}

// Spawns the phages
void Colonies3D::spawnPhages() {

    // Store the values of L and nGrid
    double L = this->L;
    int nGridXY = this->nGridXY;

    if (reducedBoundary) {
        L = round(L / s);
        nGridXY = round(static_cast<double>(nGridXY) / s);
    }

     // Determine the number of phages to spawn
    double nPhages = (double)round(L * L * H * P_0 / 1e12);

    // Apply generic spawning
    if (not experimentalConditions) {

        double numP = 0;
        if (nPhages <= nGridXY * nGridXY * nGridZ ) {
            for (double n = 0; n < nPhages; n++) {
                P(RandI(nGridXY - 1), RandI(nGridXY - 1), RandI(nGridZ - 1))++;
                numP++;
            }
        } else {
            for (uword k = 0; k < nGridZ; k++ ) {
                for (uword j = 0; j < nGridXY; j++ ) {
                    for (uword i = 0; i < nGridXY; i++) {
                        P(i, j, k) = RandP(nPhages / (nGridXY * nGridXY * nGridZ));
                        numP += P(i, j, k);
                    }
                }
            }
            // Correct for overspawning
            while (numP > nPhages) {
                uword i = RandI(nGridXY - 1);
                uword j = RandI(nGridXY - 1);
                uword k = RandI(nGridZ - 1);

                if (P(i, j, k) > 0) {
                    P(i, j, k)--;
                    numP--;
                }
            }
            // Correct for underspawning
            while (numP < nPhages) {
                uword i = RandI(nGridXY - 1);
                uword j = RandI(nGridXY - 1);
                uword k = RandI(nGridZ - 1);

                P(i, j, k)++;
                numP++;
            }
        }

    } else { // Apply scenario specific settings

        // Determine the number of phages to spawn
        double nPhages = (double)round(L * L * H * P_0 / 1e12);
        double nGridXY = this->nGridXY;

        double numP = 0;
        if (nPhages <= nGridXY * nGridXY) {
            for (double n = 0; n < nPhages; n++) {
                P(RandI(nGridXY - 1), RandI(nGridXY - 1), nGridZ - 1)++;
                numP++;
            }
        } else {
            for (uword j = 0; j < nGridXY; j++ ) {
                for (uword i = 0; i < nGridXY; i++ ) {
                        P(i, j, nGridZ - 1) = RandP(nPhages / (nGridXY * nGridXY * nGridZ));
                        numP += P(i, j, nGridZ - 1);
                }
            }
            // Correct for overspawning
            while (numP > nPhages) {
                uword i = RandI(nGridXY - 1);
                uword j = RandI(nGridXY - 1);

                if (P(i, j, nGridZ - 1) > 0) {
                    P(i, j, nGridZ - 1)--;
                    numP--;
                }
            }
            // Correct for underspawning
            while (numP < nPhages) {
                uword i = RandI(nGridXY - 1);
                uword j = RandI(nGridXY - 1);

                P(i, j, nGridZ - 1)++;
                numP++;
            }
        }
    }
}

// Computes the size of the time-step needed
void Colonies3D::ComputeTimeStep() {

    if (this->dT > 0) return;

    // Compute the step size
    double dXY = L / (double)nGridXY;
    double dZ  = H / (double)nGridZ;
    assert(dXY == dZ);
    double dx  = dXY;

    // Compute the time-step size
    int limiter = 0;

    double dT = min(pow(10,-2), 1 / nSamp);
    double dt;

    // Compute time-step limit set by D_P (LambdaP < 0.1)
    if (D_P > 0) {
        dt = pow(dx, 2) * 0.1 / (2 * D_P);
        if (dt < dT) {
            dT = dt;
            limiter = 1;
        }
    }

    // Compute time-step limit set by D_B (LambdaP < 0.1)
    if (D_B > 0) {
        dt = pow(dx, 2) * 0.1 / (2 * D_B);
        if (dt < dT) {
            dT = dt;
            limiter = 2;
        }
    }

    // Compute time-step limit set by D_n (D_n *dT/pow(dx,2) < 1/8)
    dt = pow(dx, 2) / (8 * D_n);
    if (dt < dT) {

        dT = dt;
        limiter = 3;
    }

    // Compute time-step limit set by r (r*dT < 0.25)
    if (r > 0.0) {
        dt = 0.25 / r;
        if (dt < dT) {
            dT = dt;
            limiter = 4;
        }
    }

    // Compute time-step limit set by g (g*dT < 0.1)
    dt = 0.1 / g;
    if (dt < dT) {
        dT = dt;
        limiter = 5;
    }


    // Compute time-step limit set by delta (delta*dT < 0.1)
    dt = 0.1 / delta;
    if (dt < dT) {

        dT = dt;
        limiter = 6;
    }

    // Get the order of magnitude of the timestep
    double m = floor(log10(dT));

    // Round remainder to 1, 2 or 5
    double r = round(dT * pow(10, -m));
    if (r >= 5)      dT = 5*pow(10, m);
    else if (r >= 2) dT = 2*pow(10, m);
    else             dT =   pow(10, m);

    if ( not (dT == this->dT)) {
        this->dT = dT;

        if (limiter == 1)      cout << "\tdT is Limited by D_P" << "\n";
        else if (limiter == 2) cout << "\tdT is Limited by D_B" << "\n";
        else if (limiter == 3) cout << "\tdT is Limited by D_n" << "\n";
        else if (limiter == 4) cout << "\tdT is Limited by r" << "\n";
        else if (limiter == 5) cout << "\tdT is Limited by g" << "\n";
        else if (limiter == 6) cout << "\tdT is Limited by delta" << "\n";
    }

    // Compute the jumping probabilities
    lambdaB = 2 * D_B * dT / pow(dx, 2);
    if (lambdaB > 0.1) {
        cout << "lambdaB = " << lambdaB << "\n";
        assert(lambdaB <= 0.1);
    }

    lambdaP = 2 * D_P * dT / pow(dx, 2);
    if (lambdaP > 0.1) {
        cout << "lambdaP = " << lambdaP << "\n";
        assert(lambdaP <= 0.1);
    }

}

// Returns the number of events ocurring for given n and p
double Colonies3D::ComputeEvents(double n, double p, int flag) {

    // Trivial cases
    if (p >= 1) return n;
    if (p == 0) return 0.0;
    if (n < 1)  return 0.0;

    double N = RandP(n*p);

    return round(N);
}

// Computes how many particles has moved to neighbouring points
void Colonies3D::ComputeDiffusion(double n, double lambda, double* n_0, double* n_u, double* n_d, double* n_l, double* n_r, double* n_f, double* n_b, int flag) {

    // Reset positions
    *n_0 = 0.0;
    *n_u = 0.0;
    *n_d = 0.0;
    *n_l = 0.0;
    *n_r = 0.0;
    *n_f = 0.0;
    *n_b = 0.0;

    // Trivial case
    if (n < 1) return;

    // Check if diffusion should occur
    if ((lambda == 0) or (nGridXY == 1)) {
        *n_0 = n;
        return;
    }

    if (lambda*n < 5) {   // Compute all movement individually

        for (int k = 0; k < round(n); k++) {

            double r = rand(rng);
            if       (r <    lambda)                     (*n_u)++;  // Up movement
            else if ((r >=   lambda) and (r < 2*lambda)) (*n_d)++;  // Down movement
            else if ((r >= 2*lambda) and (r < 3*lambda)) (*n_l)++;  // Left movement
            else if ((r >= 3*lambda) and (r < 4*lambda)) (*n_r)++;  // Right movement
            else if ((r >= 4*lambda) and (r < 5*lambda)) (*n_f)++;  // Forward movement
            else if ((r >= 5*lambda) and (r < 6*lambda)) (*n_b)++;  // Backward movement
            else                                         (*n_0)++;  // No movement

        }


    } else {

        // Compute the number of agents which move
        //double N = RandP(3*lambda*n); // Factor of 3 comes from 3D

        *n_u = RandP(0.5*lambda*n);
        *n_d = RandP(0.5*lambda*n);
        *n_l = RandP(0.5*lambda*n);
        *n_r = RandP(0.5*lambda*n);
        *n_f = RandP(0.5*lambda*n);
        *n_b = RandP(0.5*lambda*n);
        *n_0 = n - (*n_u + *n_d + *n_l + *n_r + *n_f + *n_b);
    }

    *n_u = round(*n_u);
    *n_d = round(*n_d);
    *n_l = round(*n_l);
    *n_r = round(*n_r);
    *n_f = round(*n_f);
    *n_b = round(*n_b);
    *n_0 = n - (*n_u + *n_d + *n_l + *n_r + *n_f + *n_b);

    assert(*n_0 >= 0);
    assert(*n_u >= 0);
    assert(*n_d >= 0);
    assert(*n_l >= 0);
    assert(*n_r >= 0);
    assert(*n_f >= 0);
    assert(*n_b >= 0);
    assert(fabs(n - (*n_0 + *n_u + *n_d + *n_l + *n_r + *n_f + *n_b)) < 1);

}


// Settings /////////////////////////////////////////////////////////////////////////////
void Colonies3D::SetLength(double L){this->L=L;}                                 // Set the side-length of the simulation
void Colonies3D::SetHeight(double H) {this->H=H;}                                // Set the height of the simulation}
void Colonies3D::SetGridSize(double nGrid){this->nGridXY=nGrid;}                 // Set the number of gridpoints
void Colonies3D::SetTimeStep(double dT){this->dT=dT;}                            // Set the time step size
void Colonies3D::SetSamples(int nSamp){this->nSamp=nSamp;}                       // Set the number of output samples

void Colonies3D::PhageInvasionStartTime(double T_i){this->T_i=T_i;}              // Sets the time when the phages should start infecting

void Colonies3D::CellGrowthRate(double g){this->g=g;}                            // Sets the maximum growthrate
void Colonies3D::CellCarryingCapacity(double K){this->K=K;}                      // Sets the carrying capacity
void Colonies3D::CellDiffusionConstant(double D_B){this->D_B=D_B;}               // Sets the diffusion constant of the phages

void Colonies3D::PhageBurstSize(int beta){this->beta=beta;}                      // Sets the size of the bursts
void Colonies3D::PhageAdsorptionRate(double eta){this->eta=eta;}                 // sets the adsorption parameter eta
void Colonies3D::PhageDecayRate(double delta){this->delta=delta;}                // Sets the decay rate of the phages
void Colonies3D::PhageInfectionRate(double r){this->r=r;}                        // Sets rate of the infection increaasing in stage
void Colonies3D::PhageDiffusionConstant(double D_P){this->D_P=D_P;}              // Sets the diffusion constant of the phages

// Sets latency time of the phage (r and tau are related by r = 10 / tau)
void Colonies3D::PhageLatencyTime(double tau) {
    if (tau > 0.0) r = 10 / tau;
    else r = 0.0;
}

void Colonies3D::SurfacePermeability(double zeta){this->zeta=zeta;}             // Sets the permeability of the surface

void Colonies3D::InitialNutrient(double n_0){this->n_0=n_0;}                    // Sets the amount of initial nutrient
void Colonies3D::NutrientDiffusionConstant(double D_n){this->D_n=D_n;}          // Sets the nutrient diffusion rate

void Colonies3D::SimulateExperimentalConditions(){experimentalConditions=true;} // Sets the simulation to spawn phages at top layer and only have x-y periodic boundaries

void Colonies3D::DisableShielding(){shielding=false;}                           // Sets shielding bool to false
void Colonies3D::DisablesClustering(){clustering=false;}                        // Sets clustering bool to false
void Colonies3D::ReducedBurstSize(){reducedBeta=true;}                          // Sets the simulation to limit beta as n -> 0
void Colonies3D::SingleInternalState(){singleInternalState=true;}               // Enables the use of a single internal state

// Sets the reduced boundary bool to true and the value of s
void Colonies3D::ReducedBoundary(int s) {
    this->s = s;
    reducedBoundary = true;
}

void Colonies3D::SetAlpha(double alpha){this->alpha=alpha;}                     // Sets the value of alpha

// Helping functions ////////////////////////////////////////////////////////////////////
// Returns random integter between 0 and n
int Colonies3D::RandI(int n) {

    // Set limit on distribution
    uniform_int_distribution <int> distr(0, n);

    return distr(rng);
}

// Returns random double between 0 and n
double Colonies3D::Rand(double n) {

    // Set limit on distribution
    uniform_real_distribution <double> distr(0, n);

    return distr(rng);
}

// Returns random normal dist. number with mean m and variance s^2
double Colonies3D::RandN(double m, double s) {

    // Set limit on distribution
    normal_distribution <double> distr(m, s);

    return distr(rng);
}

// Returns poisson dist. number with mean l
double Colonies3D::RandP(double l) {

    // Set limit on distribution
    poisson_distribution <long long> distr(l);

    return distr(rng);
}

// Sets the seed of the random number generator
void Colonies3D::SetRngSeed(int n) {
    rngSeed = n;
}

// Write a log.txt file
void Colonies3D::WriteLog() {
    if ((not f_log.is_open()) and (not exit)) {

        // Open the file stream and write the command
        f_log.open(path + "/log.txt", fstream::trunc);

        // Store the initial densities
        f_log << "B_0 = " << fixed << setw(12)  << B_0      << "\n";    // Initial density of bacteria
        f_log << "P_0 = " << fixed << setw(12)  << P_0      << "\n";    // Initial density of phages
        f_log << "n_0 = " << fixed << setw(12)  << n_0      << "\n";    // Initial density of nutrient
        f_log << "K = "   << fixed << setw(12)  << K        << "\n";    // Carrying capacity
        f_log << "L = "                         << L        << "\n";    // Side-length of simulation array
        f_log << "H = "                         << H        << "\n";    // height of simulation array
        f_log << "nGridXY = "                   << nGridXY  << "\n";    // Number of gridpoints
        f_log << "nGridZ = "                    << nGridZ   << "\n";    // Number of gridpoints
        f_log << "nSamp = "                     << nSamp    << "\n";    // Number of samples to save per simulation hour
        f_log << "g = "                         << g        << "\n";    // Growth rate for the cells
        f_log << "alpha = "                     << alpha    << "\n";    // Reinfection Percentage
        f_log << "beta = "                      << beta     << "\n";    // Multiplication factor phage
        f_log << "eta = "                       << eta      << "\n";    // Adsorption coefficient
        f_log << "delta = "                     << delta    << "\n";    // Rate of phage decay
        f_log << "r = "                         << r        << "\n";    // Constant used in the time-delay mechanism
        f_log << "zeta = "                      << zeta     << "\n";    // Permeability of surface
        f_log << "D_B = "                       << D_B      << "\n";    // Diffusion constant for the cells
        f_log << "D_P = "                       << D_P      << "\n";    // Diffusion constant for the phage
        f_log << "D_n = "                       << D_n      << "\n";    // Diffusion constant for the nutrient
        f_log << "dT = "                        << dT       << "\n";    // Time-step size
        f_log << "T_end = "                     << T_end    << "\n";    // Time when the simulation stops

        f_log << "rngSeed = "                   << rngSeed  << "\n";    // Random number seed  ( set to -1 if unused )

        f_log << "s = "                         << s        << "\n";    // The reduction of the phage boundary                       = 1;

        f_log << "experimentalConditions = "    << experimentalConditions   << "\n";
        f_log << "clustering = "                << clustering               << "\n";
        f_log << "shielding = "                 << shielding                << "\n";
        f_log << "reducedBeta = "               << reducedBeta              << "\n";
        f_log << "reducedBoundary = "           << reducedBoundary          << endl;

    }
}

// File outputs /////////////////////////////////////////////////////////////////////////

// Stop simulation when all cells are dead
void Colonies3D::FastExit(){fastExit=true;}

// Sets the simulation to export everything
void Colonies3D::ExportAll(){exportAll=true;}

// Master function to export the data
void Colonies3D::ExportData(double t){

    // Verify the file stream is open
    string fileName = "PopulationSize";
    OpenFileStream(f_N, fileName);

    // Writes the time, number of cells, number of infected cells, number of phages
    f_N << fixed    << setprecision(2);
    f_N << setw(6)  << t       << "\t";
    f_N << setw(12) << round(accu(B))    << "\t";
    f_N << setw(12) << round(accu(I0) + accu(I1) + accu(I2) + accu(I3) + accu(I4) + accu(I5) + accu(I6) + accu(I7) + accu(I8) + accu(I9))    << "\t";
    f_N << setw(12) << round(accu(P))    << "\t";

    uvec nz = find(B);
    f_N << setw(12) << static_cast<double>(nz.n_elem) / initialOccupancy << "\t";
    f_N << setw(12) << n_0 / 1e12 * pow(L, 2) * H - accu(nutrient) << "\t";
    f_N << setw(12) << round(accu(nC)) << endl;

    if (exportAll) {
        // Save the position data
        // Verify the file stream is open
        fileName = "CellDensity";
        OpenFileStream(f_B, fileName);

        fileName = "InfectedDensity";
        OpenFileStream(f_I, fileName);

        fileName = "PhageDensity";
        OpenFileStream(f_P, fileName);

        fileName = "NutrientDensity";
        OpenFileStream(f_n, fileName);

        // Write file as MATLAB would a 3D matrix!
        // row 1 is x_vector, for y_1 and z_1
        // row 2 is x_vector, for y_2 and z_1
        // row 3 is x_vector, for y_3 and z_1
        // ...
        // When y_vector for x_n has been printed, it goes:
        // row n+1 is x_vector, for y_1 and z_2
        // row n+2 is x_vector, for y_2 and z_2
        // row n+3 is x_vector, for y_3 and z_2
        // ... and so on

        // Loop over z
        for (int z = 0; z < nGridZ; z++) {

            // Loop over x
            for (int x = 0; x < nGridXY; x++) {

                // Loop over y
                for (int y = 0; y < nGridXY - 1; y++) {

                    f_B << setw(6) << round(B(x,y,z)) << "\t";
                    f_P << setw(6) << round(P(x,y,z)) << "\t";
                    double nI = round(I0(x,y,z) + I1(x,y,z) + I2(x,y,z) + I3(x,y,z) + I4(x,y,z) + I5(x,y,z) + I6(x,y,z) + I7(x,y,z) + I8(x,y,z) + I9(x,y,z));
                    f_I << setw(6) << nI       << "\t";
                    f_n << setw(6) << nutrient(x,y,z) << "\t";
                }

                // Write last line ("\n" instead of tab)
                f_B << setw(6) << round(B(x,nGridXY - 1,z)) << "\n";
                f_P << setw(6) << round(P(x,nGridXY - 1,z)) << "\n";
                double nI = round(I0(x,nGridXY - 1,z) + I1(x,nGridXY - 1,z) + I2(x,nGridXY - 1,z) + I3(x,nGridXY - 1,z) + I4(x,nGridXY - 1,z) + I5(x,nGridXY - 1,z) + I6(x,nGridXY - 1,z) + I7(x,nGridXY - 1,z) + I8(x,nGridXY - 1,z) + I9(x,nGridXY - 1,z));
                f_I << setw(6) << nI                        << "\n";
                f_n << setw(6) << round(nutrient(x,nGridXY - 1,z)) << "\n";
            }
        }
    }

}

// Open filstream if not allready opened
void Colonies3D::OpenFileStream(ofstream& stream, string& fileName) {

    // Check that if file stream is open.
    if ((not stream.is_open()) and (not exit)) {

        // Debug info
        cout << "\tSaving data to file: " << path << "/" << fileName << ".txt" << "\n";

        // Check if the output file exists
        struct stat info;
        string streamPath;
        streamPath = path + "/" + fileName + ".txt";

        // Open the file stream
        stream.open(streamPath, fstream::trunc);

        // Check stream is open
        if ((not exit) and (not stream.is_open())) {
            cerr << "\t>>Could not open filestream \"" << streamPath << "\"! Exiting..<<" << "\n";
            f_log <<  ">>Could not open filestream \"" << streamPath << "\"! Exiting..<<" << "\n";
            exit = true;
        };

        // Write meta data to the data file
        stream << "Datatype: "  << fileName << "\n";
    }
}

// Generates a save path for datafiles
string Colonies3D::GeneratePath() {

    // Generate a directory path
    string prefix = "data";    // Data folder name

    // Create the path variable
    string path_s = prefix;

    // Check if user has specified numbered folder
    if (path.empty()) {

        // Check if path exists
        struct stat info;
        if (not(stat(path_s.c_str(), &info) == 0 && S_ISDIR(info.st_mode))) {
            // Create path if it does not exist
            mkdir(path_s.c_str(), 0700);
        }

        // Loop over folders in date folder, to find current number
        int currentNumerateFolder = 1;
        DIR *dir;
        if ((dir = opendir (path_s.c_str())) != NULL) {
            struct dirent *ent;
            while ((ent = readdir (dir)) != NULL) {
                if (ent->d_type == DT_DIR) {
                    // Skip . or ..
                    if (ent->d_name[0] == '.') {continue;}
                    currentNumerateFolder++;        // Increment folder number
                }
            }
            closedir (dir);
        }

        // Append numerate folder
        path_s += "/";
        path_s += to_string(currentNumerateFolder);

        // Check if path exists
        if (not(stat(path_s.c_str(), &info) == 0 && S_ISDIR(info.st_mode))) {
            // Create path if it does not exist
            mkdir(path_s.c_str(), 0700);
        }

    } else {    // User has specified a path

        // This path maybe more than one layer deep, so attempt to make it recursively
        int len = path.length();

        // Boolean to see name of first folder
        bool firstFolder = true;

        string folder = "";
        for (int i = 0; i < len; i++) {
            folder += path[i]; // Append char to folder name

            // If seperator is found or if end of path is reached, construct folder
            if ((path[i] == '/') or (i == len - 1)) {

                // If seperator is found, remove it:
                if (path[i] == '/') folder.pop_back();

                // Check if this is the first subfolder
                if (firstFolder) {
                    firstFolder = false;

                    // Check if first folder contains date format
                    if (not ((folder.length() == 10) and(folder[4] == '-') and (folder[7] == '-'))) {

                        // Check if path exists
                        struct stat info;
                        if (not(stat(path_s.c_str(), &info) == 0 && S_ISDIR(info.st_mode))) {
                            // Create path if it does not exist
                            mkdir(path_s.c_str(), 0700);
                        }
                    }
                }

                // Append folder to path
                path_s += "/";
                path_s += folder;

                // Make folder
                struct stat info;
                if (not(stat(path_s.c_str(), &info) == 0 && S_ISDIR(info.st_mode)))
                { // Create path if it does not exist
                    mkdir(path_s.c_str(), 0700);
                }

                folder = ""; // Reset folder
            }
        }
    }

    return path_s;
}

// Sets the folder number (useful when running parralel code)
void Colonies3D::SetFolderNumber(int number) {path = to_string(number);}

// Sets the folder path (useful when running parralel code)
void Colonies3D::SetPath(std::string& path) {this->path = path;}

// Returns the save path
std::string Colonies3D::GetPath() {
    return path;
}


// Clean up /////////////////////////////////////////////////////////////////////////////

// Delete the data folder
void Colonies3D::DeleteFolder() {
    DeleteFolderTree(path.c_str());
}

// Delete folders recursively
void Colonies3D::DeleteFolderTree(const char* directory_name) {

    DIR*            dp;
    struct dirent*  ep;
    char            p_buf[512] = {0};


    dp = opendir(directory_name);

    while ((ep = readdir(dp)) != NULL) {
        // Skip self dir "."
        if (strcmp(ep->d_name, ".") == 0 || strcmp(ep->d_name, "..") == 0) continue;

        sprintf(p_buf, "%s/%s", directory_name, ep->d_name);

        // Is the path a folder?
        struct stat s_buf;
        int IsDirectory = -1;
        if (stat(p_buf, &s_buf)){
            IsDirectory = 0;
        } else {
            IsDirectory = S_ISDIR(s_buf.st_mode);
        }

        // If it is a folder, go recursively into
        if (IsDirectory) {
            DeleteFolderTree(p_buf);
        } else {    // Else delete the file
            unlink(p_buf);
        }
    }

    closedir(dp);
    rmdir(directory_name);
}

// Destructor
Colonies3D::~Colonies3D() {

    // Close filestreams
    if (f_B.is_open()) {
        f_B.flush();
        f_B.close();
    }
    if (f_I.is_open()) {
        f_I.flush();
        f_I.close();
    }
    if (f_P.is_open()) {
        f_P.flush();
        f_P.close();
    }
    if (f_N.is_open()) {
        f_N.flush();
        f_N.close();
    }
    if (f_log.is_open()) {
        f_log.flush();
        f_log.close();
    }
}
