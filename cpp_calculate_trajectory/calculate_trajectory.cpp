#include "class_trajectory.h"
#include "cpp_bfield/class_bfield.h" 
#include "cpp_runge-kutta/RK4.h"

#include <iostream>
#include <time.h>
#include <algorithm>
#include <cstdlib>
#include <ctime>

int main(void){
    clock_t begin = clock();

    Trajectory_initializer init;

    init.procID = 0;

    std::srand(time(0));
    init.seed = 0; //std::rand();                                                   //Sets the seed for the RNG. Set to 0 
                                                                                    //to generate equal results

    init.n_k = 50;                                                                 //#modes used to generate TMF
    
    init.t_end_y = 1e2;                                                             //Set time in years
    init.t_start = 0.0; 
    init.t_end = 31557600.0 * init.t_end_y; 
    init.dt = 0.1;                                                                  //Set time start and end here, and the timestep

/*
    Set B_0 as the regular magnetic field strength, given in microgauss.
    Set B_rms_turb as the turbulent magnetic field RMS-strength, given in microgauss.
*/
    init.B_0 = 1.0; init.B_rms_turb = 0.1;                             //B_0 is the regular field, B_rms_turb is the RMS of the turb field
    init.lambda_max = 10.0; init.lambda_min = 0.2;                     //set these in parsecs
/*
    Set initial conditions. E is the total energy given in eV.
    x0, y0 and z0 are the initial coordinates, given in parsec.
    vx0, vy0 and vz0 are the initial velocities, given as a percentage (must add to 1).
    The speed in each direction is calculated on the basis of the total energy and the
    given percentages.
*/
    init.E = 1.0 * 1e12;                                                
    init.x0 = 0.0; init.y0 = 0.0; init.z0 = 0.0;                        //Initial conditions
    init.vx0 = 0.9; init.vy0 = 0.0; init.vz0 = 0.1;

/* 
    Initialize charge Q, where Q = q*e, q being an integer and e the elementary charge
    Mass is given in MeV/c^2
*/
    init.q = -1.0; init.m = 0.5109989;                                  //Set the charge and mass
                                                                        
    init.max_err = 1.0e-03;                                             //Set min and max error
    init.min_err = 1.0e-08;

    init.generate_turbulence = true;                                    //This does decide if you generate a turbulence or not

    //Trajectory trajectory1(init);
    Trajectory trajectory2(init);

    //trajectory.generate_bfield(trajectory.t, trajectory.Bx, trajectory.By, trajectory.Bz);
    //trajectory.write_turbulence_to_files();
    //trajectory.write_B_to_file();

    //trajectory1.RK4t();
    trajectory2.RK4t_step_size_control();

    clock_t end = clock();
    double elapsed_secs = double(end - begin)/CLOCKS_PER_SEC;

    std::cout << "Program ended in " << elapsed_secs << " seconds \n";

}