#include <stdlib.h>
#include <chrono>
#include <stdio.h>
#include <iostream>
#include <random>
#include <math.h>
#include <time.h>
#include <vector>
#include <string>
//#include <mcheck.h>
#include <assert.h>

#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define POW5(x) ((x)*(x)*(x)*(x)*(x))

//
// Parameters
//
// Force field constants
const double eps_CH4 = 119.84;  // LJ epsilon, TraPPE, K methane = 148.0 Ar = 119.84
const double sig_CH4 = 3.405;  // LJ sigma, TraPPE, A methane = 3.73 Ar = 3.405
const double r_c = 12;  // A, LJ cutoff radius methane  = 12.8
const double r_c2 = r_c * r_c;  // LJ cutoff radius, squared
const double eps_C = 28.0; // LJ epsilon of graphene C
const double sig_C = 3.4; // LJ sigma of graphene C

// Box size
const double L_x = 40;
const double L_y = 40;
const double L_z = 18.5;  // length of box, make twice the cutoff radius 40
const double volume = L_x*L_y*L_z; // volume of sys
// const double volume = L*L*(L+graphiteThickness); // volume of graphite sys

//Graphene structural parameters
const double Delta_C = 3.35; // Interlayer spacing in A
const double graphiteArea = L_x*L_y; // A^2
const double graphiteThickness = 3.35*2; //A
const double rho = 0.382; // Number density per graphite wall 576.0/graphiteArea



// Thermodynamic parameters
// const double T = 273.15;  // temperature, K 298.0

// Monte Carlo simulation parameters
const double max_move_distance = 1;  // translation step size, A
const int sample_every = 10; // sample every this many MC moves

const bool verbose = false;  // verbose mode
const bool write_to_file = true;  // write sim results to file?
const bool write_to_position_file = true; // write position file
const bool write_to_numberDensity_file = true; // write number density file

struct Particle {
    // particle structure for storing methane positions
    double x, y, z; 
};

double GuestGuestEnergy(int which_ch4, std::vector<Particle> methanes) {
    // Compute guest-guest energy of particle id which_ch4
    double E_gg = 0.0; //guest-guest energy
    for (int i = 0; i < methanes.size(); i++) {
        // do not count interactions with itself!
        if (i == which_ch4)
            continue;
        
        // distances in each coordinate
        double x_d = methanes[i].x - methanes[which_ch4].x;
        double y_d = methanes[i].y - methanes[which_ch4].y;
        double z_d = methanes[i].z - methanes[which_ch4].z;
        
        // apply the nearest image convention for PBC
        x_d = x_d - L_x * round (x_d / L_x);
        y_d = y_d - L_y * round (y_d / L_y);
        z_d = z_d - L_z * round (z_d / L_z);
        
        // distance squared
        double r2 = x_d*x_d + y_d*y_d + z_d*z_d;
        
        if (r2 <= r_c2) {
            // if within cutoff, add energy according to LJ potential
            double sigovrr6 = CUBE(sig_CH4 * sig_CH4 / r2);
            E_gg += 4.0 * eps_CH4 * sigovrr6 * (sigovrr6 - 1.0);
        }
    }
    return E_gg;
}

double HostGuestEnergy(int which_ch4, std::vector<Particle> methanes) {
    double E_hg = 0.0; //guest-guest energy
    double eps_C_guest = pow((eps_C*eps_CH4),0.5);
    double sig_C_guest = (sig_C+sig_CH4)/2;
    // formula 2*M_PI*sig_C_guest*sig_C_guest*eps_C_guest*rho*Delta_C*(2/5*pow(sig_C_guest/z, 10)-pow(sig_C_guest/z, 4)-pow(sig_C_guest, 4)/(3*Delta_C*pow((z+0.61*Delta_C),3))
    
    double z_d = sig_C_guest;

    if (methanes[which_ch4].z <= (L_z/2)){
        z_d = methanes[which_ch4].z;
    }else{
        z_d = L_z-methanes[which_ch4].z;
    }

    // distance squared
    double r2 = z_d*z_d;

    if (r2 <= r_c2) {
        // if within cutoff, add energy according to LJ potential
        double term1 = 2*M_PI*sig_C_guest*sig_C_guest*eps_C_guest*rho;
        double term2 = 2*POW5(sig_C_guest/z_d)*POW5(sig_C_guest/z_d)/5; //pow(sig_C_guest/z_d, 10)/5
        double term3 = SQR(sig_C_guest/z_d)*SQR(sig_C_guest/z_d); //pow(sig_C_guest/z_d, 4)
        double term4 = SQR(sig_C_guest)*SQR(sig_C_guest)/(3*Delta_C*(CUBE(z_d+0.61*Delta_C))); //pow(sig_C_guest, 4), pow((z_d+0.61*Delta_C),3)
        E_hg = term1*(term2-term3-term4);
        // if (z_d <= 1.0){
        //     printf("%.3f\t%.3f\t%.3f\t%.3f\n",term1,term2,term3,term4);
        // }
    }
    
    return E_hg;
}

void TotalNumberOfMolecules(std::vector<Particle> molecules, std::vector<double> zPositions,std::vector<int> &totalNumOfMethanesAtZ){    
    for(int i=0;i<zPositions.size()-1;i++){
        int totNum = 0;
        for(int j=0;j<molecules.size();j++){
            if((molecules[j].z >= zPositions[i]) & (molecules[j].z < zPositions[i+1])){
                totNum++;
            }
        }
        totalNumOfMethanesAtZ[i] += totNum;
    }
}

double TotalEnergyOfSystem(std::vector<Particle> molecules){
    double E_gg_total = 0.0;
    double E_hg_total = 0.0;

    for (int i=0;i<molecules.size();i++){
        for (int j=i+1;j<molecules.size()-1;j++){
            // distances in each coordinate
            double x_d = molecules[i].x - molecules[j].x;
            double y_d = molecules[i].y - molecules[j].y;
            double z_d = molecules[i].z - molecules[j].z;
            
            // apply the nearest image convention for PBC
            x_d = x_d - L_x * round (x_d / L_x);
            y_d = y_d - L_y * round (y_d / L_y);
            z_d = z_d - L_z * round (z_d / L_z);
            
            // distance squared
            double r2 = x_d*x_d + y_d*y_d + z_d*z_d;
            
            if (r2 <= r_c2) {
                // if within cutoff, add energy according to LJ potential
                double sigovrr6 = CUBE(sig_CH4 * sig_CH4 / r2);
                E_gg_total += 4.0 * eps_CH4 * sigovrr6 * (sigovrr6 - 1.0);
            }
        }
        // std::cout << "E_hg_total = " << E_hg_total << "\n";
        E_hg_total += HostGuestEnergy(i,molecules);
    }
    double E_total = E_gg_total+E_hg_total;
    return E_total;
}

int main(int argc, char *argv[]) {
    //
    // Check for correct usage, then take input arguments
    //
    if (argc != 6) {
        printf("Run as:\n");
        printf("./this_code x temperature y z adsorbate \n");
        printf("\tx = Number of MC cycles \n");
        printf("\ty = fugacity Pa\n");
        printf("\tz = pressure Pa\n");
        exit(EXIT_FAILURE);
    }
    // energy of background field (K)
    // double U_0 = atof(argv[1]); 
    // Number of Monte Carlo cycles
    int N_cycles = atoi(argv[1]);
    assert(N_cycles > 0);
    double origPressure = atof(argv[4]);
    double T = atof(argv[2]);
    std::vector<double> pressures(1); // (Pa)
    pressures[0] = atof(argv[3]);
    char* adsorbate = argv[5];
    // printf("%f",3*(graphiteArea/(2.46*4.26)*4));
    // printf("%d and %f \n", N_cycles, U_0);
    //
    // Set up random number generators
    //
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::mt19937 generator (seed);
    // uniformly distributed real no in [0,1]
    std::uniform_real_distribution<double> uniform01(0.0, 1.0); 
    // uniformly distributed int; initialize with bogus range, will change later
    std::uniform_int_distribution<int> uniformint(0, 10); 
    // For picking a move: insertion, deletion, translation
    std::uniform_int_distribution<int> movepicker(0, 2); 
    
    // Number of trials to run for equilibration
    int N_equil_trials = (int)(0.5 * N_cycles); 
    printf("Simulating L_z = %.3f", L_z);
    printf("\t%d Monte Carlo cycles, %d for equilibration\n", N_cycles,
        N_equil_trials);
    
    // pressures of isotherm (actually, fugacities, computed by PREOS)
    // PREOS: https://github.com/CorySimon/PREOS
    // pressures[0] = 0.998 * 100000; 
    // pressures[1] = 5.726 * 100000; 
    // pressures[2] = 32.452 * 100000; 
    // pressures[3] = 56.7 * 100000;
    printf("Number of pressures = %lu\n", pressures.size());
    
    // store N_avg methane at pressure P here (the isotherm!)
    std::vector<double> N_of_P(pressures.size(), 0.0);

    // Compute isotherm, N(P)
    for (int i_P = 0; i_P < pressures.size(); i_P++) {
        // vector for storing methane molecule positions
        std::vector<Particle> methanes;

        // This is for calculating number density across z-directions
        std::vector<double> zPositions;
        int zPosNum = 200;
        double zPosSpacing = L_z/(double) zPosNum;
        double z0 = 0.0;
        for(int zz = 0; zz<zPosNum; zz++){
            zPositions.push_back(z0);
            z0 += zPosSpacing;
        }
        printf("Length of zPositions is %d\n",zPositions.size());
        std::vector<int> totalNumOfMethanesAtZ(zPositions.size()-1,0);
        std::vector<double> avgNumOfMethanesAtZ(zPositions.size()-1,0);

        double P = pressures[i_P]; // get pressure (Pa)
        printf("GCMC simulation at P = %f Pa\n", P);

        // Initialize statistics
        double N_avg = 0.0;
        int N_samples = 0;
        int N_accept = 0;
        int N_ch4 = 0;

        int MC_counter = 0;

        // set up output file
        char outputname[1024];
        sprintf(outputname, "./Graphene%sResults/SimulationResults/SimulationResults_%.2f_%.2f.txt",adsorbate,origPressure,T);
        FILE *outfile;
        if (write_to_file) { 
            outfile = fopen(outputname, "w");
            if (outfile == NULL) {
                printf("ERROR: could not open output file %s\n", outputname);
                exit(EXIT_FAILURE);
            }
            fprintf(outfile, "L_z = %f A\n", L_z);
            fprintf(outfile, "N_cycles = %d , N_equil = %d\n", 
                N_cycles, N_equil_trials);
            fprintf(outfile, "Pressure(bar),Loading(vSTP/v),Loading(N_avg),Qst(J/mol)\n");
        }
  
        // for calculating qst
        
        double U_total_avg = 0.0;
        double U_X_N = 0.0;
        double N_sqr = 0.0;

        for (int cycle=0; cycle < N_cycles; cycle++) {
            int N_inner_cycles = (N_ch4 > 20) ? N_ch4 : 20; //Ensuring the minimum number of molecules in the system is 20
            // printf("Cycle Number is : %d and Number of inner cycle is %d\n",cycle,N_inner_cycles);
            for (int i = 0; i < N_inner_cycles; i++) { // The MC moves will be applied to the number of existing molecules in the system at most
                MC_counter += 1;
                
                //
                // take samples
                //
                if ((cycle >= N_equil_trials) & (MC_counter % sample_every == 0)) {
                    N_avg += N_ch4; //divide by N_samples later
                    N_samples++;
                    TotalNumberOfMolecules(methanes,zPositions,totalNumOfMethanesAtZ);
                    // For qst calculation
                    double U_total = TotalEnergyOfSystem(methanes);
                    U_total_avg += U_total;
                    U_X_N += U_total*N_ch4;
                    N_sqr += pow(N_ch4,2);
                }
          
                // Choose insertion, deletion, or translation
                int which_move = movepicker(generator);

                //
                // Insertion
                //
                if (which_move == 0) {
                    // Create new methane molecule at random position in box
                    Particle new_methane;
                    new_methane.x = L_x * uniform01(generator);
                    new_methane.y = L_y * uniform01(generator);
                    new_methane.z = L_z * uniform01(generator);

                    methanes.push_back(new_methane);  // add methane
                    
                    // get energy of new methane
                    double E_gg = GuestGuestEnergy(N_ch4, methanes);
                    double E_hg = HostGuestEnergy(N_ch4, methanes);
                    // if (cycle == N_cycles-1){
                    //     printf("E_gg=%.3f , E_hg=%.3f\n",E_gg,E_hg);
                    // }

                    // dE of this move. energy with other guests + field energy
                    double dE = E_gg + E_hg;  

                    double acceptance_insertion = P * volume / 
                        ((N_ch4 + 1) * 1.3806488e7 * T) * exp(-dE / T);
                    if (uniform01(generator) < acceptance_insertion) { 
                        // accept insertion
                        N_ch4++;
                        N_accept++;
                    } 
                    else {
                        // remove methane from vector
                        methanes.pop_back();
                    }
                }  // end insertion
          
                //
                //Deletion
                //
                if (which_move == 1 && N_ch4 > 0) {
                    // choose which methane molecule to attempt to delete
                    decltype(uniformint.param()) new_range(0, N_ch4 - 1);
                    uniformint.param(new_range);
                    int which_ch4 = uniformint(generator);

                    if (verbose) 
                        printf("DEBUG: trying deletion of particle %d\n",
                            which_ch4);
                    
                    // calculate energy of this methane
                    double E_gg = GuestGuestEnergy(which_ch4, methanes);
                    double E_hg = HostGuestEnergy(which_ch4, methanes);

                    // printf("E_gg=%.3f , E_hg=%.3f\n",E_gg,E_hg);
                    // calculate dE of propose deletion
                    double dE = E_gg + E_hg;

                    double acceptance_del = (N_ch4 * 1.3806488e7 * T) / 
                        (P * volume) * exp(dE / T);
                    if (uniform01(generator) < acceptance_del) { 
                        // accept deletion
                        // erase this guest from the adsorbates vector
                        methanes.erase(methanes.begin() + which_ch4);
                        
                        N_ch4--;
                        N_accept++;

                        if (verbose) 
                            printf("DEBUG: deleted! now there are %d particles\n", N_ch4);
                    }
                }
          
                //
                //Translation
                //
                if (which_move == 2 && N_ch4 > 0) {
                    // choose which methane to attempt to translate
                    decltype(uniformint.param()) new_range(0, N_ch4 - 1);
                    uniformint.param(new_range);
                    int which_ch4 = uniformint(generator);

                    if (verbose) 
                        printf("DEBUG: trying translation of particle ID %d\n", which_ch4);

                    // Calculate energy at old configuration
                    double E_old = GuestGuestEnergy(which_ch4, methanes)+HostGuestEnergy(which_ch4, methanes);

                    // store old position
                    Particle old_position = methanes[which_ch4];
            
                    // Perturb coordinates by:
                    double dx = max_move_distance * (uniform01(generator) - 0.5);
                    double dy = max_move_distance * (uniform01(generator) - 0.5);
                    double dz = max_move_distance * (uniform01(generator) - 0.5);
                    
                    // Move this methane molecule
                    // If move outside of box, use periodic BCs
                    methanes[which_ch4].x = fmod(methanes[which_ch4].x + dx,  L_x);
                    methanes[which_ch4].y = fmod(methanes[which_ch4].y + dy,  L_y);
                    methanes[which_ch4].z = fmod(methanes[which_ch4].z + dz,  L_z);
           
                    if (methanes[which_ch4].x < 0.0)
                        methanes[which_ch4].x += L_x;
                    if (methanes[which_ch4].y < 0.0)
                        methanes[which_ch4].y += L_y;
                    if (methanes[which_ch4].z < 0.0)
                        methanes[which_ch4].z += L_z;

                    // Calculate energy at new configuration
                    double E_new = GuestGuestEnergy(which_ch4, methanes)+HostGuestEnergy(which_ch4, methanes);

                    // here U_0 cancels
                    double dE = E_new - E_old; 
                    double acceptance_move = exp(-dE / T);
                    if (uniform01(generator) < acceptance_move) {
                        // keep new position
                        N_accept++;
                    }
                    else {
                        // restore old position
                        methanes[which_ch4] = old_position;
                    }
                } // end translation
            } // end inner cycle loop
            // set up output file for storing molecule positions
            if ((cycle >= N_equil_trials) & (cycle%1000==0)){
                char outputNamePositions[1024];
                sprintf(outputNamePositions,"./Graphene%sResults/MoleculePositions/position_%d_%.2f_%.2f.csv",adsorbate,cycle,origPressure,T);
                FILE *positionFile;
                if (write_to_position_file) { 
                    positionFile = fopen(outputNamePositions, "w");
                    if (positionFile == NULL) {
                        printf("ERROR: could not open output file %s\n", outputNamePositions);
                        exit(EXIT_FAILURE);
                    }
                    fprintf(positionFile, "x, y, z\n");
                    for (int mol=0; mol <N_ch4;mol++){
                        fprintf(positionFile,"%.3f, %.3f, %.3f\n",methanes[mol].x,methanes[mol].y,methanes[mol].z);
                    }
                    fclose(positionFile);
                }
            }
        } // end outer cycle loop
        printf("\tDone. %d/%d MC moves accepted.\n", N_accept, MC_counter);
        printf("\t%d samples taken.\n", N_samples);

        N_avg = N_avg / ((double) N_samples);
        N_of_P[i_P] = N_avg;
        printf("\tAverage %s molecules at pressure %.2f: %.3f\n", adsorbate, P, N_avg);
        printf("\tIG predicted: N_avg = %f\n" , P * volume / T / 8.314 *6.022e-7);
        double vSTPv = N_avg / volume / 6.022 * 1e7 * 22.4 / 1000.0;

        // for qst calculation 
        U_total_avg = U_total_avg / ((double) N_samples);
        U_X_N += U_X_N / ((double) N_samples);
        N_sqr += N_sqr / ((double) N_samples);

        double qst = 8.314*T - (U_X_N - U_total_avg*N_avg)/(N_sqr - pow(N_avg,2))*1.38*6.023;

        if (write_to_file)
            fprintf(outfile, "%f,%f,%f,%f\n", pressures[i_P], vSTPv, N_avg, qst);
        if (write_to_file)
            fclose(outfile);

        // For z direction molecular distribution
        for(int n=0;n<totalNumOfMethanesAtZ.size();n++){
            avgNumOfMethanesAtZ[n] = (double) totalNumOfMethanesAtZ[n]/N_samples/(L_x*L_y*zPosSpacing);
        }

        char outputNameNumberDensity[1024];
        sprintf(outputNameNumberDensity,"./Graphene%sResults/NumberDensity/NumberDensity_%.2f_%.2f.csv",adsorbate,origPressure,T);
        FILE *numberDensityFile;
        if (write_to_numberDensity_file) { 
            numberDensityFile = fopen(outputNameNumberDensity, "w");
            if (numberDensityFile == NULL) {
                printf("ERROR: could not open output file %s\n", outputNameNumberDensity);
                exit(EXIT_FAILURE);
            }
            fprintf(numberDensityFile, "zPositions, avgNumOfMethanesAtZ\n");
            for (int n=0; n < avgNumOfMethanesAtZ.size();n++){
                fprintf(numberDensityFile,"%.3f, %.3f\n",zPositions[n],avgNumOfMethanesAtZ[n]);
            }
            fclose(numberDensityFile);
        }
    }  // end pressure loop
    

    printf("\nProgram complete\n\n");
}
