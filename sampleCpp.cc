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

//
// Parameters
//
// Force field constants
const double eps_CH4 = 148.0;  // LJ epsilon, TraPPE, K
const double sig_CH4 = 3.73;  // LJ sigma, TraPPE, A
const double r_c = 12.8;  // A, LJ cutoff radius
const double r_c2 = r_c * r_c;  // LJ cutoff radius, squared
const double eps_C = 28.0; // LJ epsilon of graphene C
const double sig_C = 3.4; // LJ sigma of graphene C

//Graphene structural parameters
const double Delta_C = 3.35; // Interlayer spacing in A
const double graphiteArea = 40*40; // A^2
const double graphiteThickness = 3.35*2; //A
const double rho = 1728/(graphiteArea*graphiteThickness); // Number density per graphite wall


// Box size
const double L = 40;  // length of box, make twice the cutoff radius
const double volume = pow(L, 3); // volume of sys

// Thermodynamic parameters
const double T = 298.0;  // temperature, K

// Monte Carlo simulation parameters
const double max_move_distance = 0.1;  // translation step size, A
const int sample_every = 10; // sample every this many MC moves

const bool verbose = false;  // verbose mode
const bool write_to_file = true;  // write sim results to file?
const bool write_to_position_file = true; // write position file

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

double HostGuestEnergyTest(double z_d) {
    double E_hg = 0.0; //guest-guest energy
    double eps_C_guest = pow((eps_C*eps_CH4),0.5);
    double sig_C_guest = (sig_C+sig_CH4)/2;
    // formula 2*M_PI*sig_C_guest*sig_C_guest*eps_C_guest*rho*Delta_C*(2/5*pow(sig_C_guest/z, 10)-pow(sig_C_guest/z, 4)-pow(sig_C_guest, 4)/(3*Delta_C*pow((z+0.61*Delta_C),3))

    // distance squared
    double r2 = z_d*z_d;

    if (r2 <= r_c2) {
        // if within cutoff, add energy according to LJ potential
        double term1 = 2*M_PI*sig_C_guest*sig_C_guest*eps_C_guest*rho*Delta_C;
        double term2 = 2*pow(sig_C_guest/z_d, 10)/5;
        double term3 = pow(sig_C_guest/z_d, 4);
        double term4 = pow(sig_C_guest, 4)/(3*Delta_C*pow((z_d+0.61*Delta_C),3));
        E_hg = term1*(term2-term3-term4);
        if (z_d <= 1.0){
            printf("%.3f\t%.3f\t%.3f\t%.3f\t%.4f\n",term1,term2,term3,term4,2*pow(sig_C_guest/z_d, 10)/5);
        }
    }
    
    return E_hg;
}

int main(int argc, char *argv[]){
    double E_hg = HostGuestEnergyTest(12);
    printf("%f",E_hg);
}