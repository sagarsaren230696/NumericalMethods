#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <string>
#include <iostream>
// #include <random>
// #include <assert.h>
// #include <chrono>

//
// Parameters
//
// Force field constants
const double eps_CH4 = 148.0;  // LJ epsilon, TraPPE, K
const double sig_CH4 = 3.73;  // LJ sigma, TraPPE, A
const double r_c = 12.8;  // A, LJ cutoff radius
const double r_c2 = r_c * r_c;  // LJ cutoff radius, squared

// Box size
const double L = r_c * 2;  // length of box, make twice the cutoff radius
const double volume = pow(L, 3); // volume of sys

// Thermodynamic parameters
const double T = 298.0;  // temperature, K

// Monte Carlo simulation parameters
const double max_move_distance = 0.1;  // translation step size, A
const int sample_every = 10; // sample every this many MC moves

const bool verbose = false;  // verbose mode
const bool write_to_file = true;  // write sim results to file?

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
        x_d = x_d - L * round (x_d / L);
        y_d = y_d - L * round (y_d / L);
        z_d = z_d - L * round (z_d / L);
        
        // distance squared
        double r2 = x_d*x_d + y_d*y_d + z_d*z_d;
        
        if (r2 <= r_c2) {
            // if within cutoff, add energy according to LJ potential
            double sigovrr6 = pow(sig_CH4 * sig_CH4 / r2, 3);
            E_gg += 4.0 * eps_CH4 * sigovrr6 * (sigovrr6 - 1.0);
        }
    }
    return E_gg;
}

int main(int argc, char *argv[]) {
    if (argc != 3) {
        printf("Run as:\n");
        printf("./this_code x y\n");
        printf("\tx = energy of background field (K)\n");
        printf("\ty = number of Monte Carlo cycles\n");
        exit(EXIT_FAILURE);
    }

    // // energy of background field (K)
    // double U_0 = atof(argv[1]); 
    // // Number of Monte Carlo cycles
    // int N_cycles = atoi(argv[2]);
    // int N_equil_trials = (int)(0.5 * N_cycles); 
    // printf("%d and %f \n", N_cycles, U_0);
    // // assert(N_cycles > 0);
    

    //
    // Set up random number generators
    //
    // unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    // std::mt19937 generator (seed);
    // // uniformly distributed real no in [0,1]
    // std::uniform_real_distribution<double> uniform01(0.0, 1.0); 
    // // uniformly distributed int; initialize with bogus range, will change later
    // std::uniform_int_distribution<int> uniformint(0, 10); 
    // // For picking a move: insertion, deletion, translation
    // std::uniform_int_distribution<int> movepicker(0, 2); 
    
    // Number of trials to run for equilibration
    // printf("Simulating L = %.3f", L);
    // printf("\t%d Monte Carlo cycles, %d for equilibration\n", N_cycles,
    //     N_equil_trials);
    
    // pressures of isotherm (actually, fugacities, computed by PREOS)
    // PREOS: https://github.com/CorySimon/PREOS
    

    try{
        std::vector<double> pressures(4); // (Pa)
        pressures[0] = 0.998 * 100000; 
        pressures[1] = 5.726 * 100000; 
        pressures[2] = 32.452 * 100000; 
        pressures[3] = 56.7 * 100000;
        // printf("%f",pressures[0]);
        std::cout << pressures[1];
    }catch(const std::bad_array_new_length &e) {
        std::cout << e.what();
    }
    // printf("Number of pressures = %lu\n", pressures.size());
    
    // store N_avg methane at pressure P here (the isotherm!)
    // std::vector<double> N_of_P(pressures.size(), 0.0);
    // char outputname[1024];
    // sprintf(outputname, "./results_empty_box/SIM_U%.2f.txt", U_0);
    // FILE *outfile;
    // if (true) { 
    //     outfile = fopen(outputname, "w");
    //     if (outfile == NULL) {
    //         printf("ERROR: could not open output file %s\n", outputname);
    //         exit(EXIT_FAILURE);
    //     }
    //     fprintf(outfile, "L = %f A\n", 12.8);
    //     fprintf(outfile, "N_cycles = %d , N_equil = %d\n", 
    //         N_cycles, N_equil_trials);
    //     fprintf(outfile, "Pressure(bar),Loading(vSTP/v),Loading(N_avg)\n");
    // }
    
    return 0;
}