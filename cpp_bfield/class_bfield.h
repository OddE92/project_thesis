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
    vector< vector<double> > B, turbulencevec, B0vec;
    vector<double> Bvec_scal, turbAtPoint;


    //Constructors and destructor
    Bfield();
    Bfield(double xyz_max, double xyz_min, int numGridPoints, bool Bgenerate_turbulence);
    Bfield(Trajectory_initializer &init);
    ~Bfield();

    //Functions for generating the B-field
    template <class Bfield_func_x, class Bfield_func_y, class Bfield_func_z>
    void generate_bfield(double t, Bfield_func_x &Bx, Bfield_func_y &By, Bfield_func_z &Bz);
    template <class Bfield_func_x, class Bfield_func_y, class Bfield_func_z>
    void generate_bfield_at_point(Bfield_func_x &Bx, Bfield_func_y &By, Bfield_func_z &Bz, double t, double &bx, double &by, double &bz,
                                  double x, double y, double z, bool do_gen_turb_at_point);

    //General functions
    void generate_turbulence(int x);
    int generate_turbulence_at_point(double x, double y, double z);
    void initialize_turbulence();

    int get_bfield(double& bx, double& by, double& bz, const double x, const double y, const double z);

    void write_turbulence_to_files();
    void write_B_to_file();

    Ran ran;

private:

    //functions to initialize the turbulence. Only run these once!
    void initialize_phases();
    void initialize_normalization();
    void initialize_coords();

    //Bools to check if functions need to be run or not.
    bool generate_bfield_is_run = false;
    bool turbulence_is_initialized = false;
    

/********** CONSTANTS **********/

    //pi = M_PI
    const double two_pi = 2 * M_PI;

    complex<double> im;                                     //For working with complex numbers

    //Magnetic field
    const int n_k = 50;                                     //Nr. of modes used to generate the turbulence
    double B_0 = 1.0;
    double B_rms_turb = 1.0;                                //Normalized magnetic field,
    const double gamma = 5.0 / 3.0;                         //Power law for the fluctuation spectrum
    const double lambda_min = 0.2;                          //Smallest wavelength
    const double lambda_max = 10.0;                         //Largest wavelength
    const double k_min = two_pi / lambda_max;               //Smallest wavenumber
    const double k_max = two_pi / lambda_min;               //Largest wavenumber

/********* END CONSTANTS ********/

/********* RANDOM PHASES ********/
    
    vector<double> a, b, p, t, s;                           //a = alpha, b = beta, p = phi, t = theta, s = sign
    vector<double> ca, sa, cp, sp, ct, st;                  //c = cos, s = sin; ca = cos(a(n_k));

/******* END RANDOM PHASES ******/

/***** NORMALIZED PARAMETERS ****/

    vector<double> B_k, k;
    double dB_min;

/*** END NORMALIZED PARAMETERS **/

/********* COORDINATES **********/

    int c_size = 10, c2, c3;
    vector< vector<double> > coord;                         //Vector of length c3, with 3 coordinates in each row
    double scal = 0.3, cscal = 1.0;                         //Scaling factors for the plot.
    double xyz_max, xyz_min;                                //max and min value for the coordinates (square box volume)

/******* END COORDINATES ********/
/********************************/
};

/********************************/
/******* GENERATE B-FIELD *******/
template <class Bfield_func_x, class Bfield_func_y, class Bfield_func_z>
void Bfield::generate_bfield(double t, Bfield_func_x &Bx, Bfield_func_y &By, Bfield_func_z &Bz){
    
    if(!turbulence_is_initialized)
    initialize_turbulence();                                                //Setup for generating turbulence

    for(int i = 0; i < c3; i++){                                            //For each point:
        B0vec[i][0] = Bx(t, coord[i][0], coord[i][1], coord[i][2]);         //Generate field in x-direction
        B0vec[i][1] = By(t, coord[i][0], coord[i][1], coord[i][2]);         //Generate field in y-direction
        B0vec[i][2] = Bz(t, coord[i][0], coord[i][1], coord[i][2]);         //Generate field in z-direction
        generate_turbulence(i);                                             //Generate turbulence in all directions

        B[i][0] = B0vec[i][0] + turbulencevec[i][0];                        //The total field is the sum of the turbulence
        B[i][1] = B0vec[i][1] + turbulencevec[i][1];                        //and background field
        B[i][2] = B0vec[i][2] + turbulencevec[i][2];
    }
    generate_bfield_is_run = true;
}
/***** END GENERATE B-FIELD *****/

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
