#include <vector>
#include <complex>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <time.h>
#include <random>
#include <cmath>

#include "class_bfield.h"

using namespace std;

int Bfield::generate_turbulence_at_point(double x, double y, double z){
    vector<double> e_k(n_k);
    
    vector< complex<double> > delta_B(3);
    vector< complex<double> > F(n_k), epsilon_x(n_k), epsilon_y(n_k), epsilon_z(n_k), B_x_k(n_k), B_y_k(n_k), B_z_k(n_k);

    // This function takes a position r, and generate a "random" turbulence to the B-field

    // Follow the articles calculation of delta_B (or delta_Omega, where delta_B = delta_Omega * mc/q)

    for(int i = 0; i < n_k; i++){

        // e_k is the exponent of delta_B
        e_k[i] = k[i] * (st[i]*cp[i]*x + st[i]*sp[i]*y + ct[i]*z) + b[i];

        // F = B(k)*e^i(kz+beta), B_k = B(k)
        F[i].real(B_k[i] * cos(e_k[i]));  F[i].imag(B_k[i] * sin(e_k[i]));

        //Epsilon is a vector, so it is split into epsilon_(x/y/z), where (x/y/z) is the direction in the non-rotated coord.sys.
        epsilon_x[i].real( ca[i] * ct[i] * cp[i]);          epsilon_x[i].imag(-s[i] * sa[i] * sp[i]);
        epsilon_y[i].real( ca[i] * ct[i] * sp[i]);          epsilon_y[i].imag( s[i] * sa[i] * cp[i]);
        epsilon_z[i].real(-ca[i] * st[i]);                  epsilon_z[i].imag(0.0);

        //Calculate B(k) in the x-, y- and z-directions
        B_x_k[i] = F[i] * epsilon_x[i];                      // Complex number multiplication!
        B_y_k[i] = F[i] * epsilon_y[i];
        B_z_k[i] = F[i] * epsilon_z[i];
    }

    for(int i = 0; i < delta_B.size(); i++){                 // Need to set delta_B to 0 before we start summing
        delta_B[i].imag(0.0); delta_B[i].real(0.0);
    }

        //Calculates the final delta_B in the x-, y- and z-directions.
    for(int i = 0; i < n_k; i++){
        delta_B[0] += B_x_k[i];
        delta_B[1] += B_y_k[i];
        delta_B[2] += B_z_k[i];
    }

// As we are concerned about the real value of the B-field we need to calculate the real part of delta_B

    for(int i = 0; i < delta_B.size(); i++){
        turbAtPoint[i] = real(delta_B[i]);
    }

// B = m*c*Omega/q, Omega = Omega_0 + dOmega. Check units!
    return 0;
}

/**** END GENERATE TURBULENCE ***/


/********** INITIALIZE **********/

void Bfield::initialize_turbulence(){
    if (!turbulence_is_initialized){
        initialize_phases();
        initialize_normalization();

        turbulence_is_initialized = true;
    }
}

void Bfield::reinitialize_turbulence(){

    initialize_phases();
    initialize_normalization();

    turbulence_is_initialized = true;
    
}


/********************************/
/****** INITIALIZE PHASES *******/

void Bfield::initialize_phases(){
    double r;

    for(int i = 0; i < n_k; i++){                           // From the article:
        a[i] = two_pi * ran.doub();                         // a = alpha
        b[i] = two_pi * ran.doub();                         // b = beta
        p[i] = two_pi * ran.doub();                         // p = phi
        t[i] = M_PI * ran.doub();
        //t[i] = std::acos(1 - 2*ran.doub());                 // t = theta with p(t) = sin(t)/2, 0<t<pi

        s[i] = 1.0;                                         // s = sign in epsilon (+-)
        r = ran.doub();
        if(r > 0.5) s[i] = -1.0;

        cp[i] = cos(p[i]);                                  //sines and cosines of theta, phi and alpha
        sp[i] = sin(p[i]);
        ct[i] = cos(t[i]);
        st[i] = sin(t[i]);
        ca[i] = cos(a[i]);
        sa[i] = cos(a[i]);
    }

}
/***** END INITIALIZE PHASES ****/
/********************************/

/********************************/
/*** INITIALIZE NORMALIZATION ***/

void Bfield::initialize_normalization(){

    double logk_diff, dlog_k;

    logk_diff = log10(k_max) - log10(k_min);
    dlog_k = logk_diff/(n_k - 1);


    //The for loop spreads the values of k equally along the logarithmic scale.
    for(int i = 0; i < n_k; i++){
        k[i] = pow(10, log10(k_min) + i * dlog_k);
    }

/*
    Normalization: The energy density S should stay the same for any number of wavevectors (n_k) used.
    Energy density is given by:
    S = sum( B(k)^2 / 8pi )

    The total energy density can be normalized by setting S = (B_0)^2 / 8pi, which gives

    (B_0)^2 / 8pi = B(k_min)^2 / 8pi * sum[ (k / k_min)^-gamma ]; B(k_min) = dB_min

    B(k_min) = sqrt{ (B_0)^2 / sum[ (k/k_min)^-gamma ]}

    Remember to check units against giacalone & jokipii 1994, where they use Omega instead of B.
*/
  
  double sum = 0;
  for(int i = 0; i < n_k; i++){
      sum += pow(k[i]/k_min, -gamma);
  }
  dB_min = B_rms_turb * sqrt(2 / sum);
  for(int i = 0; i < n_k; i++){
      B_k[i] = dB_min * pow(k[i] / k_min, -gamma/2);
  }

  //Normalization is now done after the value given for B_rms_turb (default 1.0).  
}
/****** END INITIALIZE NORM *****/
/********************************/


Bfield::Bfield() : ran(15321){
    
    // See class declaration for vector descriptions.

    turbAtPoint.resize(3);

    im.imag(1.0); im.real(0.0);

    a.resize(n_k); b.resize(n_k); p.resize(n_k); t.resize(n_k); s.resize(n_k);              
    ca.resize(n_k); sa.resize(n_k); cp.resize(n_k); sp.resize(n_k); ct.resize(n_k); 
    st.resize(n_k);                                                                         

    B_k.resize(n_k); k.resize(n_k);
}



Bfield::Bfield( Trajectory_initializer &init) : B_0(init.B_0), B_rms_turb(init.B_rms_turb), ran(15321 + init.seed){
    
    // See class declaration for vector descriptions.
    n_k = init.n_k;
  
    turbAtPoint.resize(3);

    lambda_max = init.lambda_max;
    lambda_min = init.lambda_min;

    k_min = two_pi / lambda_max;                     //Smallest wavenumber
    k_max = two_pi / lambda_min;                     //Largest wavenumber

    im.imag(1.0); im.real(0.0);

    a.resize(n_k); b.resize(n_k); p.resize(n_k); t.resize(n_k); s.resize(n_k);              
    ca.resize(n_k); sa.resize(n_k); cp.resize(n_k); sp.resize(n_k); ct.resize(n_k); 
    st.resize(n_k);                                                                         

    B_k.resize(n_k); k.resize(n_k);
    initialize_turbulence();

}

Bfield::~Bfield(){
    //do nothing
}
