#include "class_trajectory.h"
#include "cpp_bfield/class_bfield.h" 
#include "cpp_runge-kutta/RK4.h"
#include "cpp_generate_samples/calculate_eigenvalues.h"

#include <iostream>
#include <time.h>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <random>

int main(void){
    clock_t begin = clock();

    Trajectory_initializer init;
    
    std::srand(time(0));  //rand() +                                    //Add this line to seed for "more" randomness
    int seed = rand() + 1000;

    Ran rng(seed);

    double ran1, ran2, ran3, ran4, ran5, ran6, ranTot;
    //Pick three random numbers to generate direction of particle
    ran1 = rng.doub(); ran2 = rng.doub(); ran3 = rng.doub();
    ranTot = ran1 + ran2 + ran3;


    //Pick three random numbers to set the signs
    ran4 = rng.doub(); ran5 = rng.doub(); ran6 = rng.doub();
    if(ran4 > 0.5){ ran4 = -1.0; }else{ ran4 = 1.0; }
    if(ran5 > 0.5){ ran5 = -1.0; }else{ ran5 = 1.0; }
    if(ran6 > 0.5){ ran6 = -1.0; }else{ ran6 = 1.0; }


    init.procID = 0;

    std::srand(time(0));
    init.seed = 1;//rand();                                                             //Sets the seed for the RNG. Set to 0 
                                                                                    //to generate equal results

    init.n_k = 500;                                                                 //#modes used to generate TMF
    
    init.t_end_y = 1e5;                                                             //Set time in years
    init.D_ij_size = log10(init.t_end_y) * 9 + 1;
    init.t_start = 0.0; 
    init.t_end = 31557600.0 * init.t_end_y; 
    //init.t_end = 1000.0;
    init.dt = 0.1;                                                                  //Set time start and end here, and the timestep

/*
    Set B_0 as the regular magnetic field strength, given in microgauss.
    Set B_rms_turb as the turbulent magnetic field RMS-strength, given in microgauss.
*/
    init.B_0 = 0.0; init.B_rms_turb = 4.0;                              //B_0 is the regular field, B_rms_turb is the RMS of the turb field
    init.lambda_max = 10; init.lambda_min = 0.0027;                     //set these in parsecs
/*
    Set initial conditions. E is the total energy given in eV.
    x0, y0 and z0 are the initial coordinates, given in parsec.
    vx0, vy0 and vz0 are the initial velocities, given as a percentage (must add to 1).
    The speed in each direction is calculated on the basis of the total energy and the
    given percentages.
*/
    init.E = 1.0 * 1e17;                                                
    init.x0 = 0.0; init.y0 = 0.0; init.z0 = 0.0;                        //Initial conditions
    init.vx0 = 0.9; init.vy0 = 0.0; init.vz0 = 0.1;
    init.N = 1;
/* 
    Initialize charge Q, where Q = q*e, q being an integer and e the elementary charge
    Mass is given in MeV/c^2
*/
    init.q = 1.0; init.m = 938.2720813;                                 //Set the charge and mass
                                                                        
    init.max_err = 1.0e-05;                                             //Set min and max error
    init.min_err = 1.0e-08;

    init.generate_turbulence = true;                                    //This does decide if you generate a turbulence or not
    
    //Trajectory trajectory1(init);
    Trajectory trajectory2(init);
    
    trajectory2.vx[0] = ran4 * (ran1/ranTot) * trajectory2.v_tot_init;
    trajectory2.vy[0] = ran5 * (ran2/ranTot) * trajectory2.v_tot_init;
    trajectory2.vz[0] = ran6 * (ran3/ranTot) * trajectory2.v_tot_init;

    trajectory2.vx[0] = 0.9*trajectory2.v_tot_init;
    trajectory2.vy[0] = 0.0*trajectory2.v_tot_init;
    trajectory2.vz[0] = 0.1*trajectory2.v_tot_init;

    //trajectory.generate_bfield(trajectory.t, trajectory.Bx, trajectory.By, trajectory.Bz);
    //trajectory.write_turbulence_to_files();
    //trajectory.write_B_to_file();

    std::cout << "vx0: " << trajectory2.vx[0] << " vy0: " << trajectory2.vy[0] << " vz0: " << trajectory2.vz[0] << std::endl;
    trajectory2.RK4t_step_size_control();
    
 /*       
    trajectory2.RK4t_step_size_control_nw();

    for(int i = 0; i < trajectory2.D_ij.size(); i += 6){
      for(int j = 0; j < 6; j++){

        //i will represent the time in logarithmic years (i.e t_y = 10^(i/6))
        trajectory2.D_ij[i + j] = trajectory2.D_ij[i + j] / (2.0 * init.N * trajectory2.D_ij_time[i/6]);
        //D_ij has now been calculated for one instance of the magnetic field.

      }
    }//End for i

    std::vector<double> eigenvalues_current(3,0);

    for(int i = 0; i < trajectory2.D_ij.size() - 1; i += 6){
           
        calculate_eigenvalues_3x3_sym(trajectory2.D_ij, i, eigenvalues_current);
        cout << eigenvalues_current[0] << ' ' << eigenvalues_current[1] << ' ' << eigenvalues_current[2] << '\n';
    
    }//End for i 
        
 */   
/*
    ofstream test_rng;
    test_rng.open("data/test_rng.dat");

    if(!test_rng) std::cout << "couldn't open test_rng.dat. \n";

    for(int i = 0; i < 10000; i++){
        test_rng << rng.doub() << ' ' << rng.doub() << '\n';
    }

    test_rng.close();
*/

    ofstream fileMF;
    fileMF.open("data/test_mf.dat");

    if(!fileMF) cout<<"couldn't open test_mf.dat\n";

    for(double x = 0; x < 500; x += 0.1){
        trajectory2.generate_turbulence_at_point(x,0,0);
        fileMF << trajectory2.turbAtPoint[0] << ' ' << trajectory2.turbAtPoint[1] << ' ' << trajectory2.turbAtPoint[2] << '\n';
    }

    init.lambda_max = 150;
    Trajectory trajectory3(init);
    
    for(double x = 0; x < 500; x += 0.1){
        trajectory3.generate_turbulence_at_point(x,0,0);
        fileMF << trajectory3.turbAtPoint[0] << ' ' << trajectory3.turbAtPoint[1] << ' ' << trajectory3.turbAtPoint[2] << '\n';
    }

    fileMF.close();

    clock_t end = clock();
    double elapsed_secs = double(end - begin)/CLOCKS_PER_SEC;

    std::cout << "Program ended in " << elapsed_secs << " seconds \n";

}