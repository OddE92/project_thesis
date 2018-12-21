/**********************************************

 * This program is a translation of the structure.f90 code written by Kristian Joten Andersen
 * for his thesis "Charged Particle Trajectories in the Local Superbubble".

 * The program was translated to C++ by Odd-Einar C. Nervik for his project thesis, under the guidance of
 * Michael Kachlerieß.
 
 * Some of the comments in the code are taken directly from the structure.f90-program.
 * The program follows the method shown in:

 *      Charged-particle motion in multidimensional magnetic-field turbulence
 *      Giacalone, J. and Jokipii, J. R.
 *      Department of Planetary Sciences, university of Arizona, Tucson AZ
 *      May 1994

 * The final calculation of the turbulent field follows the equation in:
 
 *      The transport of cosmic rays across a turbulent magnetic field
 *      Giacalone, J. and Jokipii, J. R.
 *      Department of Planetary Sciences, university of Arizona, Tucson AZ
 *      February 1999

 * Which has a different rotation and translation direction than the 1994 article.

 * The output of the program are two files. The first three numbers are the spatial coordinates x, y and z.
 * The last three numbers are Bx, By and Bz respectively.
 * *** The part where a box volume with a magnetic field has been deprecated. Now the field is generated at the provided point instead.
 
 * The second file is scaled for improved readability when plotted.

**********************************************/

#ifndef CLASS_BFIELD
#define CLASS_BFIELD

#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <complex>
#include <vector>

#include "ran.h"
#include "cpp_calculate_trajectory/class_trajectory_initializer.h"

using std::vector; using std::complex;

/********************************/
/********* CLASS BFIELD *********/

class Bfield{
public:
    vector<double> turbAtPoint;


    //Constructors and destructor
    Bfield();
    Bfield(Trajectory_initializer &init);
    ~Bfield();

    //Functions for generating the B-field
    template <class Bfield_func_x, class Bfield_func_y, class Bfield_func_z>
    void generate_bfield_at_point(Bfield_func_x &Bx, Bfield_func_y &By, Bfield_func_z &Bz, double t, double &bx, double &by, double &bz,
                                  double x, double y, double z, bool do_gen_turb_at_point);

    //General functions
    void generate_turbulence(int x);
    int generate_turbulence_at_point(double x, double y, double z);
    void initialize_turbulence();
    void reinitialize_turbulence();



    Ran ran;

private:

    //functions to initialize the turbulence. Only run these once!
    void initialize_phases();
    void initialize_normalization();

    //Bools to check if functions need to be run or not.
    bool turbulence_is_initialized = false;
    

/********** CONSTANTS **********/

    //pi = M_PI
    const double two_pi = 2 * M_PI;
    const double mtopc = 3.24078e-17;

    complex<double> im;                                     //For working with complex numbers

    //Magnetic field
    int n_k = 50;                                           //Nr. of modes used to generate the turbulence
    double B_0 = 1.0;
    double B_rms_turb = 1.0;                                //Normalized magnetic field,
    const double gamma = 5.0 / 3.0;                         //Power law for the fluctuation spectrum
    double lambda_min = 0.2;                                //Smallest wavelength, in parsecs
    double lambda_max = 10.0;                               //Largest wavelength, in parsecs
    double k_min = two_pi / lambda_max;                     //Smallest wavenumber
    double k_max = two_pi / lambda_min;                     //Largest wavenumber

/********* END CONSTANTS ********/

/********* RANDOM PHASES ********/
    
    vector<double> a, b, p, t, s;                           //a = alpha, b = beta, p = phi, t = theta, s = sign
    vector<double> ca, sa, cp, sp, ct, st;                  //c = cos, s = sin; ca = cos(a(n_k));

/******* END RANDOM PHASES ******/

/***** NORMALIZED PARAMETERS ****/

    vector<double> B_k, k;
    double dB_min;

/*** END NORMALIZED PARAMETERS **/

};

/********** GENERATE FIELD AT POINT **********/
template <class Bfield_func_x, class Bfield_func_y, class Bfield_func_z>
void Bfield::generate_bfield_at_point(Bfield_func_x &Bx, Bfield_func_y &By, Bfield_func_z &Bz, double t, double &bx, double &by, double &bz,
                                      double x, double y, double z, bool do_gen_turb_at_point){

/* 
    This functions takes the position and time, generates the turbulence (without reinitializing it) 
    and calculates the regular field at the position and time. B_0 is used to scale the regular field.
*/
    
    if(!turbulence_is_initialized){ 
        initialize_turbulence();
    }
    if(do_gen_turb_at_point) {
        generate_turbulence_at_point(x, y, z);
    }

    bx = B_0 * Bx(t, x, y, z) + turbAtPoint[0];
    by = B_0 * By(t, x, y, z) + turbAtPoint[1];
    bz = B_0 * Bz(t, x, y, z) + turbAtPoint[2];

}
/********** END GENERATE FIELD AT POINT **********/

#endif