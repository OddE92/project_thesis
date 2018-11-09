#ifndef CLASS_TRAJECTORY
#define CLASS_TRAJECTORY

#include "cpp_bfield/class_bfield.h"
#include "cpp_calculate_trajectory/class_trajectory_initializer.h"

/********************************/
/*********** FUNCTORS ***********/
/*
    Set the functions in each directions in accordance with the desired magnetic field.
*/
struct Bfield_func_x{
    double operator()(double t, double x, double y, double z){
        return 0;
    }
};
struct Bfield_func_y{
    double operator()(double t, double x, double y, double z){
        return 0;
    }
};
struct Bfield_func_z{
    double operator()(double t, double x, double y, double z){
        return 1;
    }
};
/********* END FUNCTORS *********/

class Trajectory: public Bfield {
public:
    vector<double> x, y, z, vx, vy, vz;
    vector<double> D_ij;

    bool do_gen_turb_at_point = false;

    double t, q = 1.0, m = 1.0;
    double min_err, max_err;
    double vx0, vy0, vz0;
    double R_Larmor;                                //Larmor-radius (gyroradius)
    double gamma_l;                                 //Lorentz-factor
    
    double E;
    double unit_coeff;
    double v_tot_init;
    const double mtopc = 3.24078e-17;
    const double c = 2.99792458e8;
    
    int RK4t();
    int RK4t_step_size_control();
    int RK4t_step_size_control_nw();

    Trajectory();
    Trajectory(Trajectory_initializer &init);
    ~Trajectory();

    Bfield_func_x Bx;
    Bfield_func_y By;
    Bfield_func_z Bz;

    double t_start = 0.0, t_end = 1000.0, dt = 0.1;
    double qm;
};



#endif