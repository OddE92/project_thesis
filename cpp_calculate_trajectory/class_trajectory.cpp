#include "class_trajectory.h"
#include "cpp_bfield/class_bfield.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <stdio.h>

using std::cout; using std::endl; using std::ofstream; using std::string;

double frexp10(double arg, int * exp)
{
   *exp = (arg == 0) ? 0 : 1 + (int)std::floor(std::log10(std::fabs(arg) ) );
   return arg * std::pow(10 , -(*exp));    
}

/******************************************/
/************ RK4 STEP SIZE CTRL **********/
int Trajectory::RK4t_step_size_control(){

/*
    The step size controls uses two RK4-solvers. The first one goes from t -> t+dt in one step.
    The second goes from t -> t + dt/2, then from t + dt/2 -> t + dt. This should in theory give 
    a more accurate approximation. The error is calculated as difference in the solution between the 
    half and whole-step approximations. In this equation the total speed is used for the truncation error.
    v(dt2) - v(dt1) = truncation error estimate, where
    v = sqrt( vx^2 + vy^2 + vz^2 )    
*/

    t = t_start;
    double bx, by, bz;                                      //To hold field at the given point
    double dt2 = dt/2;                                      //Step size for size controlled part
    double truncation_error;                                //To hold truncation error
    double total_velocity_dt1, total_velocity_dt2;          //To calculate truncation error
    int nbad = 0, ntoolow = 0;                              //Counters to check number of steps where dt changes

    vector<double> kvxdt1(4), kvxdt2(4), lvydt1(4), lvydt2(4), hvzdt1(4), hvzdt2(4);    //Runge-kutta coefficients
    double vxdt1, vydt1, vzdt1, vxdt2, vydt2, vzdt2;        //to hold speed and check error

    ofstream file;
    ofstream fileTurb;
/*
    double mantissa; int exp;                               //Generate a unique filename.
    mantissa = frexp10(E, &exp);
    
    string fileName = "data/Traj_E" + std::to_string(static_cast<int>(10*mantissa)) + "e" + std::to_string(exp - 1);
    fileName += "_vx0=" + std::to_string(static_cast<int>(100*vx0))
              + "_vy0=" + std::to_string(static_cast<int>(100*vy0))
              + "_vz0=" + std::to_string(static_cast<int>(100*vz0));
    fileName += "_pID=" + std::to_string(procID) + ".dat"; 
    
    file.open(fileName);
*/
    string fileName = "cpp_runge-kutta/RK4_approx_stepsizectrl.dat";       //For running test

    file.open(fileName);
    fileTurb.open("cpp_runge-kutta/turbulence.dat");
    //file << setprecision(20);

    if(!file){
        cout << "Couldn't open " << fileName << endl;
    }

    file << x[0] << ' ' << y[0] << ' ' << z[0] << '\n';             //write initial position to file

    //Generate the bfield at the starting point
    generate_bfield_at_point(Bx, By, Bz, t, bx, by, bz, x[0], y[0], z[0], do_gen_turb_at_point);
    fileTurb << turbAtPoint[0] << ' ' << turbAtPoint[1] << ' ' << turbAtPoint[2] << '\n';
    int i;
    for (i = 0; t <= t_end; i++){
        dt2 = dt / 2;

        //Standard implementation of a RK4-method for 3 coupled ODEs.
        kvxdt1[0] = dt  * unit_coeff * (vy[0] * bz - vz[0] * by);
        lvydt1[0] = dt  * unit_coeff * (vz[0] * bx - vx[0] * bz);
        hvzdt1[0] = dt  * unit_coeff * (vx[0] * by - vy[0] * bx);
        kvxdt2[0] = dt2 * unit_coeff * (vy[0] * bz - vz[0] * by);
        lvydt2[0] = dt2 * unit_coeff * (vz[0] * bx - vx[0] * bz);
        hvzdt2[0] = dt2 * unit_coeff * (vx[0] * by - vy[0] * bx);

        kvxdt1[1] = dt  * unit_coeff * ((vy[0] + lvydt1[0] / 2) * bz - (vz[0] + hvzdt1[0] / 2) * by);
        lvydt1[1] = dt  * unit_coeff * ((vz[0] + hvzdt1[0] / 2) * bx - (vx[0] + kvxdt1[0] / 2) * bz);
        hvzdt1[1] = dt  * unit_coeff * ((vx[0] + kvxdt1[0] / 2) * by - (vy[0] + lvydt1[0] / 2) * bx);
        kvxdt2[1] = dt2 * unit_coeff * ((vy[0] + lvydt2[0] / 2) * bz - (vz[0] + hvzdt2[0] / 2) * by);
        lvydt2[1] = dt2 * unit_coeff * ((vz[0] + hvzdt2[0] / 2) * bx - (vx[0] + kvxdt2[0] / 2) * bz);
        hvzdt2[1] = dt2 * unit_coeff * ((vx[0] + kvxdt2[0] / 2) * by - (vy[0] + lvydt2[0] / 2) * bx);

        kvxdt1[2] = dt  * unit_coeff * ((vy[0] + lvydt1[1] / 2) * bz - (vz[0] + hvzdt1[1] / 2) * by);
        lvydt1[2] = dt  * unit_coeff * ((vz[0] + hvzdt1[1] / 2) * bx - (vx[0] + kvxdt1[1] / 2) * bz);
        hvzdt1[2] = dt  * unit_coeff * ((vx[0] + kvxdt1[1] / 2) * by - (vy[0] + lvydt1[1] / 2) * bx);
        kvxdt2[2] = dt2 * unit_coeff * ((vy[0] + lvydt2[1] / 2) * bz - (vz[0] + hvzdt2[1] / 2) * by);
        lvydt2[2] = dt2 * unit_coeff * ((vz[0] + hvzdt2[1] / 2) * bx - (vx[0] + kvxdt2[1] / 2) * bz);
        hvzdt2[2] = dt2 * unit_coeff * ((vx[0] + kvxdt2[1] / 2) * by - (vy[0] + lvydt2[1] / 2) * bx);

        kvxdt1[3] = dt  * unit_coeff * ((vy[0] + lvydt1[2]) * bz - (vz[0] + hvzdt1[2]) * by);
        lvydt1[3] = dt  * unit_coeff * ((vz[0] + hvzdt1[2]) * bx - (vx[0] + kvxdt1[2]) * bz);
        hvzdt1[3] = dt  * unit_coeff * ((vx[0] + kvxdt1[2]) * by - (vy[0] + lvydt1[2]) * bx);
        kvxdt2[3] = dt2 * unit_coeff * ((vy[0] + lvydt2[2]) * bz - (vz[0] + hvzdt2[2]) * by);
        lvydt2[3] = dt2 * unit_coeff * ((vz[0] + hvzdt2[2]) * bx - (vx[0] + kvxdt2[2]) * bz);
        hvzdt2[3] = dt2 * unit_coeff * ((vx[0] + kvxdt2[2]) * by - (vy[0] + lvydt2[2]) * bx);

/***** velocity calculation *****/
        vxdt1 = vx[0] + (1.0 / 6) * (kvxdt1[0] + 2 * (kvxdt1[1] + kvxdt1[2]) + kvxdt1[3]);
        vydt1 = vy[0] + (1.0 / 6) * (lvydt1[0] + 2 * (lvydt1[1] + lvydt1[2]) + lvydt1[3]);
        vzdt1 = vz[0] + (1.0 / 6) * (hvzdt1[0] + 2 * (hvzdt1[1] + hvzdt1[2]) + hvzdt1[3]);
        vxdt2 = vx[0] + (1.0 / 6) * (kvxdt2[0] + 2 * (kvxdt2[1] + kvxdt2[2]) + kvxdt2[3]);
        vydt2 = vy[0] + (1.0 / 6) * (lvydt2[0] + 2 * (lvydt2[1] + lvydt2[2]) + lvydt2[3]);
        vzdt2 = vz[0] + (1.0 / 6) * (hvzdt2[0] + 2 * (hvzdt2[1] + hvzdt2[2]) + hvzdt2[3]);
/********************************/
//Continue another halfstep for dt2

        kvxdt2[0] = dt2 * unit_coeff * (vydt2 * bz - vzdt2 * by);
        lvydt2[0] = dt2 * unit_coeff * (vzdt2 * bx - vxdt2 * bz);
        hvzdt2[0] = dt2 * unit_coeff * (vxdt2 * by - vydt2 * bx);

        kvxdt2[1] = dt2 * unit_coeff * ((vydt2 + lvydt2[0] / 2) * bz - (vzdt2 + hvzdt2[0] / 2) * by);
        lvydt2[1] = dt2 * unit_coeff * ((vzdt2 + hvzdt2[0] / 2) * bx - (vxdt2 + kvxdt2[0] / 2) * bz);
        hvzdt2[1] = dt2 * unit_coeff * ((vxdt2 + kvxdt2[0] / 2) * by - (vydt2 + lvydt2[0] / 2) * bx);

        kvxdt2[2] = dt2 * unit_coeff * ((vydt2 + lvydt2[1] / 2) * bz - (vzdt2 + hvzdt2[1] / 2) * by);
        lvydt2[2] = dt2 * unit_coeff * ((vzdt2 + hvzdt2[1] / 2) * bx - (vxdt2 + kvxdt2[1] / 2) * bz);
        hvzdt2[2] = dt2 * unit_coeff * ((vxdt2 + kvxdt2[1] / 2) * by - (vydt2 + lvydt2[1] / 2) * bx);

        kvxdt2[3] = dt2 * unit_coeff * ((vydt2 + lvydt2[2]) * bz - (vzdt2 + hvzdt2[2]) * by);
        lvydt2[3] = dt2 * unit_coeff * ((vzdt2 + hvzdt2[2]) * bx - (vxdt2 + kvxdt2[2]) * bz);
        hvzdt2[3] = dt2 * unit_coeff * ((vxdt2 + kvxdt2[2]) * by - (vydt2 + lvydt2[2]) * bx);

        vxdt2 = vxdt2 + (1.0 / 6) * (kvxdt2[0] + 2 * (kvxdt2[1] + kvxdt2[2]) + kvxdt2[3]);
        vydt2 = vydt2 + (1.0 / 6) * (lvydt2[0] + 2 * (lvydt2[1] + lvydt2[2]) + lvydt2[3]);
        vzdt2 = vzdt2 + (1.0 / 6) * (hvzdt2[0] + 2 * (hvzdt2[1] + hvzdt2[2]) + hvzdt2[3]);

/***********************************************************/
/********** CALCULATE AND CHECK TRUNCATION ERROR ***********/
        total_velocity_dt1 = sqrt(pow(vxdt1, 2) + pow(vydt1, 2) + pow(vzdt1, 2));
        total_velocity_dt2 = sqrt(pow(vxdt2, 2) + pow(vydt2, 2) + pow(vzdt2, 2));
        truncation_error = std::abs(total_velocity_dt2 - total_velocity_dt1);

        if (truncation_error <= max_err && truncation_error >= min_err){    //if error is within desired interval
            //Set next velocity
            vx[1] = vxdt2;
            vy[1] = vydt2;
            vz[1] = vzdt2;

            //Calculate the next position
            x[1] = x[0] + mtopc * vx[0] * dt;
            y[1] = y[0] + mtopc * vy[0] * dt;
            z[1] = z[0] + mtopc * vz[0] * dt;

            file << x[1] << ' ' << y[1] << ' ' << z[1] << '\n';             //write position to file

            t += dt;       

            //Generate the field at the new point.
            generate_bfield_at_point(Bx, By, Bz, t, bx, by, bz, x[1], y[1], z[1], do_gen_turb_at_point);
            fileTurb << turbAtPoint[0] << ' ' << turbAtPoint[1] << ' ' << turbAtPoint[2] << '\n';

            vx[0] = vx[1]; vy[0] = vy[1]; vz[0] = vz[1];                    //Setup for next iteration
            x[0] = x[1]; y[0] = y[1]; z[0] = z[1];

        }else if(truncation_error > max_err){           //If the truncation error is too large, we need to decrease the stepsize
            dt = dt2;
            i--;                                        //decrease i to redo the step
            nbad++;

        }else if(truncation_error < min_err){           //If the truncation error is too small, take a larger step next time.
            //Set next velocity
            vx[1] = vxdt2;
            vy[1] = vydt2;
            vz[1] = vzdt2;

            //Calculate the next position
            x[1] = x[0] + mtopc * vx[0] * dt;
            y[1] = y[0] + mtopc * vy[0] * dt;
            z[1] = z[0] + mtopc * vz[0] * dt;

            file << x[1] << ' ' << y[1] << ' ' << z[1] << '\n';             //write position to file

            t += dt; 
            dt = 2*dt;                                                      //double the step size
            ntoolow++;

            //Generate the field at the new point.
            generate_bfield_at_point(Bx, By, Bz, t, bx, by, bz, x[1], y[1], z[1], do_gen_turb_at_point);
            fileTurb << turbAtPoint[0] << ' ' << turbAtPoint[1] << ' ' << turbAtPoint[2] << '\n';
              
            vx[0] = vx[1]; vy[0] = vy[1]; vz[0] = vz[1];                    //Setup for next iteration
            x[0] = x[1]; y[0] = y[1]; z[0] = z[1];
            
        }else{
            cout << "Something unexpected happened. Truncation error tests weren't hit. \n";
            return 1;
        }

    }//End for i

    file.close();
    fileTurb.close();
    
    cout << "Number of steps taken: " << i << endl;
    cout << "Number of bad steps taken: " << nbad << endl;
    cout << "Number of steps with too high accuracy: " << ntoolow << endl;
    cout << "Initial velocity: " << v_tot_init << " m/s" << endl;
    cout << "vx0: " << vx[0] << " ; vy0: " << vy[0] << " ; vz0: " << vz[0] << " m/s" << endl;
    cout << "Larmor radius: " << R_Larmor << " parsec" << endl;
    cout << "dt: " << dt << endl;

    return 0;
}

int Trajectory::RK4t_step_size_control_nw(){                //This version does not write anything to any files

/*
    The step size controls uses two RK4-solvers. The first one goes from t -> t+dt in one step.
    The second goes from t -> t + dt/2, then from t + dt/2 -> t + dt. This should in theory give 
    a more accurate approximation. The error is calculated as difference in the solution between the 
    half and whole-step approximations. In this equation the total speed is used for the truncation error.
    v(dt2) - v(dt1) = truncation error estimate, where
    v = sqrt( vx^2 + vy^2 + vz^2 )    
*/

    t = t_start;                                            
    double bx, by, bz;                                      //To hold field at the given point
    double dt2 = dt/2;                                      //Step size for size controlled part
    double truncation_error;                                //To hold truncation error
    double total_velocity_dt1, total_velocity_dt2;          //To calculate truncation error
    int nbad = 0, ntoolow = 0;                              //Counters to check number of steps where dt changes
    int count = 0;                                          //Counter to track when to send D_ij

    vector<double> kvxdt1(4), kvxdt2(4), lvydt1(4), lvydt2(4), hvzdt1(4), hvzdt2(4);    //Runge-kutta coefficients
    double vxdt1, vydt1, vzdt1, vxdt2, vydt2, vzdt2;        //to hold speed and check error

    //Generate the bfield at the starting point
    generate_bfield_at_point(Bx, By, Bz, t, bx, by, bz, x[0], y[0], z[0], do_gen_turb_at_point);
    
    int i;
    for (i = 0; t <= t_end; i++){
        dt2 = dt / 2;

        //Standard implementation of a RK4-method for 3 coupled ODEs.
        kvxdt1[0] = dt  * unit_coeff * (vy[0] * bz - vz[0] * by);
        lvydt1[0] = dt  * unit_coeff * (vz[0] * bx - vx[0] * bz);
        hvzdt1[0] = dt  * unit_coeff * (vx[0] * by - vy[0] * bx);
        kvxdt2[0] = dt2 * unit_coeff * (vy[0] * bz - vz[0] * by);
        lvydt2[0] = dt2 * unit_coeff * (vz[0] * bx - vx[0] * bz);
        hvzdt2[0] = dt2 * unit_coeff * (vx[0] * by - vy[0] * bx);

        kvxdt1[1] = dt  * unit_coeff * ((vy[0] + lvydt1[0] / 2) * bz - (vz[0] + hvzdt1[0] / 2) * by);
        lvydt1[1] = dt  * unit_coeff * ((vz[0] + hvzdt1[0] / 2) * bx - (vx[0] + kvxdt1[0] / 2) * bz);
        hvzdt1[1] = dt  * unit_coeff * ((vx[0] + kvxdt1[0] / 2) * by - (vy[0] + lvydt1[0] / 2) * bx);
        kvxdt2[1] = dt2 * unit_coeff * ((vy[0] + lvydt2[0] / 2) * bz - (vz[0] + hvzdt2[0] / 2) * by);
        lvydt2[1] = dt2 * unit_coeff * ((vz[0] + hvzdt2[0] / 2) * bx - (vx[0] + kvxdt2[0] / 2) * bz);
        hvzdt2[1] = dt2 * unit_coeff * ((vx[0] + kvxdt2[0] / 2) * by - (vy[0] + lvydt2[0] / 2) * bx);

        kvxdt1[2] = dt  * unit_coeff * ((vy[0] + lvydt1[1] / 2) * bz - (vz[0] + hvzdt1[1] / 2) * by);
        lvydt1[2] = dt  * unit_coeff * ((vz[0] + hvzdt1[1] / 2) * bx - (vx[0] + kvxdt1[1] / 2) * bz);
        hvzdt1[2] = dt  * unit_coeff * ((vx[0] + kvxdt1[1] / 2) * by - (vy[0] + lvydt1[1] / 2) * bx);
        kvxdt2[2] = dt2 * unit_coeff * ((vy[0] + lvydt2[1] / 2) * bz - (vz[0] + hvzdt2[1] / 2) * by);
        lvydt2[2] = dt2 * unit_coeff * ((vz[0] + hvzdt2[1] / 2) * bx - (vx[0] + kvxdt2[1] / 2) * bz);
        hvzdt2[2] = dt2 * unit_coeff * ((vx[0] + kvxdt2[1] / 2) * by - (vy[0] + lvydt2[1] / 2) * bx);

        kvxdt1[3] = dt  * unit_coeff * ((vy[0] + lvydt1[2]) * bz - (vz[0] + hvzdt1[2]) * by);
        lvydt1[3] = dt  * unit_coeff * ((vz[0] + hvzdt1[2]) * bx - (vx[0] + kvxdt1[2]) * bz);
        hvzdt1[3] = dt  * unit_coeff * ((vx[0] + kvxdt1[2]) * by - (vy[0] + lvydt1[2]) * bx);
        kvxdt2[3] = dt2 * unit_coeff * ((vy[0] + lvydt2[2]) * bz - (vz[0] + hvzdt2[2]) * by);
        lvydt2[3] = dt2 * unit_coeff * ((vz[0] + hvzdt2[2]) * bx - (vx[0] + kvxdt2[2]) * bz);
        hvzdt2[3] = dt2 * unit_coeff * ((vx[0] + kvxdt2[2]) * by - (vy[0] + lvydt2[2]) * bx);

/***** velocity calculation *****/
        vxdt1 = vx[0] + (1.0 / 6) * (kvxdt1[0] + 2 * (kvxdt1[1] + kvxdt1[2]) + kvxdt1[3]);
        vydt1 = vy[0] + (1.0 / 6) * (lvydt1[0] + 2 * (lvydt1[1] + lvydt1[2]) + lvydt1[3]);
        vzdt1 = vz[0] + (1.0 / 6) * (hvzdt1[0] + 2 * (hvzdt1[1] + hvzdt1[2]) + hvzdt1[3]);
        vxdt2 = vx[0] + (1.0 / 6) * (kvxdt2[0] + 2 * (kvxdt2[1] + kvxdt2[2]) + kvxdt2[3]);
        vydt2 = vy[0] + (1.0 / 6) * (lvydt2[0] + 2 * (lvydt2[1] + lvydt2[2]) + lvydt2[3]);
        vzdt2 = vz[0] + (1.0 / 6) * (hvzdt2[0] + 2 * (hvzdt2[1] + hvzdt2[2]) + hvzdt2[3]);
/********************************/
//Continue another halfstep for dt2

        kvxdt2[0] = dt2 * unit_coeff * (vydt2 * bz - vzdt2 * by);
        lvydt2[0] = dt2 * unit_coeff * (vzdt2 * bx - vxdt2 * bz);
        hvzdt2[0] = dt2 * unit_coeff * (vxdt2 * by - vydt2 * bx);

        kvxdt2[1] = dt2 * unit_coeff * ((vydt2 + lvydt2[0] / 2) * bz - (vzdt2 + hvzdt2[0] / 2) * by);
        lvydt2[1] = dt2 * unit_coeff * ((vzdt2 + hvzdt2[0] / 2) * bx - (vxdt2 + kvxdt2[0] / 2) * bz);
        hvzdt2[1] = dt2 * unit_coeff * ((vxdt2 + kvxdt2[0] / 2) * by - (vydt2 + lvydt2[0] / 2) * bx);

        kvxdt2[2] = dt2 * unit_coeff * ((vydt2 + lvydt2[1] / 2) * bz - (vzdt2 + hvzdt2[1] / 2) * by);
        lvydt2[2] = dt2 * unit_coeff * ((vzdt2 + hvzdt2[1] / 2) * bx - (vxdt2 + kvxdt2[1] / 2) * bz);
        hvzdt2[2] = dt2 * unit_coeff * ((vxdt2 + kvxdt2[1] / 2) * by - (vydt2 + lvydt2[1] / 2) * bx);

        kvxdt2[3] = dt2 * unit_coeff * ((vydt2 + lvydt2[2]) * bz - (vzdt2 + hvzdt2[2]) * by);
        lvydt2[3] = dt2 * unit_coeff * ((vzdt2 + hvzdt2[2]) * bx - (vxdt2 + kvxdt2[2]) * bz);
        hvzdt2[3] = dt2 * unit_coeff * ((vxdt2 + kvxdt2[2]) * by - (vydt2 + lvydt2[2]) * bx);

        vxdt2 = vxdt2 + (1.0 / 6) * (kvxdt2[0] + 2 * (kvxdt2[1] + kvxdt2[2]) + kvxdt2[3]);
        vydt2 = vydt2 + (1.0 / 6) * (lvydt2[0] + 2 * (lvydt2[1] + lvydt2[2]) + lvydt2[3]);
        vzdt2 = vzdt2 + (1.0 / 6) * (hvzdt2[0] + 2 * (hvzdt2[1] + hvzdt2[2]) + hvzdt2[3]);

/***********************************************************/
/********** CALCULATE AND CHECK TRUNCATION ERROR ***********/
        total_velocity_dt1 = sqrt(pow(vxdt1, 2) + pow(vydt1, 2) + pow(vzdt1, 2));
        total_velocity_dt2 = sqrt(pow(vxdt2, 2) + pow(vydt2, 2) + pow(vzdt2, 2));
        truncation_error = std::abs(total_velocity_dt2 - total_velocity_dt1);

        if (truncation_error <= max_err && truncation_error >= min_err){    //if error is within desired interval
            //Set next velocity
            vx[1] = vxdt2;
            vy[1] = vydt2;
            vz[1] = vzdt2;

            //Calculate the next position
            x[1] = x[0] + mtopc * vx[0] * dt;
            y[1] = y[0] + mtopc * vy[0] * dt;
            z[1] = z[0] + mtopc * vz[0] * dt;

            t += dt;       

            //Generate the field at the new point.
            generate_bfield_at_point(Bx, By, Bz, t, bx, by, bz, x[1], y[1], z[1], do_gen_turb_at_point);

            vx[0] = vx[1]; vy[0] = vy[1]; vz[0] = vz[1];                    //Setup for next iteration
            x[0] = x[1]; y[0] = y[1]; z[0] = z[1];

        }else if(truncation_error > max_err){           //If the truncation error is too large, we need to decrease the stepsize
            dt = dt2;
            i--;                                        //decrease i to redo the step
            nbad++;

        }else if(truncation_error < min_err){           //If the truncation error is too small, take a larger step next time.
            //Set next velocity
            vx[1] = vxdt2;
            vy[1] = vydt2;
            vz[1] = vzdt2;

            //Calculate the next position
            x[1] = x[0] + mtopc * vx[0] * dt;
            y[1] = y[0] + mtopc * vy[0] * dt;
            z[1] = z[0] + mtopc * vz[0] * dt;

            t += dt; 
            dt = 2*dt;                                                      //double the step size
            ntoolow++;

            //Generate the field at the new point.
            generate_bfield_at_point(Bx, By, Bz, t, bx, by, bz, x[1], y[1], z[1], do_gen_turb_at_point);
              
            vx[0] = vx[1]; vy[0] = vy[1]; vz[0] = vz[1];                    //Setup for next iteration
            x[0] = x[1]; y[0] = y[1]; z[0] = z[1];
            
        }else{
            cout << "Something unexpected happened. Truncation error tests weren't hit. \n";
            return 1;
        }

        if( log10(t/31557600) - count > __DBL_EPSILON__ ){                  //This should scompute D_ij for every whole integer on the 
                                                                            //log scale with time in years
            
            //add each particles coordinate to the sum in D_ij
            //D_ij holds the upper half, including the diagonal, of a symmetric matrix
            //stored as a 1D-array to comply with MPIs send and recieve functions.
            D_ij[6 * count + 0] += x[0] * x[0];                             //D_11
            D_ij[6 * count + 1] += x[0] * y[0];                             //D_12
            D_ij[6 * count + 2] += x[0] * z[0];                             //D_13
            D_ij[6 * count + 3] += y[0] * y[0];                             //D_22
            D_ij[6 * count + 4] += y[0] * z[0];                             //D_23
            D_ij[6 * count + 5] += z[0] * z[0];                             //D_33

            count++;
        }


    }//End for i
    
    return 0;
}//End RK4t_stepsizectrl_nw

/******************************************/
/********** RK4 NO STEP SIZE CTRL *********/
int Trajectory::RK4t(){

    t = t_start;
    double bx, by, bz;                                  //to hold field at the given point

    double kvx1, kvx2, kvx3, kvx4, lvy1, lvy2, lvy3, lvy4, hvz1, hvz2, hvz3, hvz4;


    ofstream file;

    
    file.open("cpp_runge-kutta/RK4_approx_t.dat");

    //file << setprecision(20);


    if(!file){
        cout << "Couldn't open RK4_approx_t.dat" << endl;
    }

    for (int i = 0; t<=t_end; i++){

        generate_bfield_at_point(Bx, By, Bz, t, bx, by, bz, x[0]/mtopc, y[0]/mtopc, z[0]/mtopc, do_gen_turb_at_point);
        

        //Standard implementation of a RK4-method for 3 coupled ODEs.
        kvx1 = dt * unit_coeff * (vy[0] * bz - vz[0] * by);
        lvy1 = dt * unit_coeff * (vz[0] * bx - vx[0] * bz);
        hvz1 = dt * unit_coeff * (vx[0] * by - vy[0] * bx);

        kvx2 = dt * unit_coeff * ((vy[0] + lvy1 / 2) * bz - (vz[0] + hvz1 / 2) * by);
        lvy2 = dt * unit_coeff * ((vz[0] + hvz1 / 2) * bx - (vx[0] + kvx1 / 2) * bz);
        hvz2 = dt * unit_coeff * ((vx[0] + kvx1 / 2) * by - (vy[0] + lvy1 / 2) * bx);

        kvx3 = dt * unit_coeff * ((vy[0] + lvy2 / 2) * bz - (vz[0] + hvz2 / 2) * by);
        lvy3 = dt * unit_coeff * ((vz[0] + hvz2 / 2) * bx - (vx[0] + kvx2 / 2) * bz);
        hvz3 = dt * unit_coeff * ((vx[0] + kvx2 / 2) * by - (vy[0] + lvy2 / 2) * bx);

        kvx4 = dt * unit_coeff * ((vy[0] + lvy3) * bz - (vz[0] + hvz3) * by);
        lvy4 = dt * unit_coeff * ((vz[0] + hvz3) * bx - (vx[0] + kvx3) * bz);
        hvz4 = dt * unit_coeff * ((vx[0] + kvx3) * by - (vy[0] + lvy3) * bx);

        //Calculate next position
        x[1] = x[0] + mtopc * vx[0] * dt;
        y[1] = y[0] + mtopc * vy[0] * dt;
        z[1] = z[0] + mtopc * vz[0] * dt;

        //Calculate next velocity
        vx[1] = vx[0] + (1.0 / 6) * (kvx1 + 2 * kvx2 + 2 * kvx3 + kvx4);
        vy[1] = vy[0] + (1.0 / 6) * (lvy1 + 2 * lvy2 + 2 * lvy3 + lvy4);
        vz[1] = vz[0] + (1.0 / 6) * (hvz1 + 2 * hvz2 + 2 * hvz3 + hvz4);

        file << x[0] << ' ' << y[0] << ' ' << z[0] << '\n';

        t += dt;
    }

    file.close();
    
    return 0;

}
/********  END RK4 NO STEP SIZE CTRL ******/

/******************************************/
/************** CONSTRUCTORS **************/

Trajectory::Trajectory() : Trajectory::Bfield(){
    qm = q/m;
    max_err = 1.0e-05;
    min_err = 1.0e-09;

    x.resize(2);  y.resize(2);  z.resize(2);
    vx.resize(2); vy.resize(2); vz.resize(2);

    x[0] = 0; y[0] = 0; z[0] = 0;                           //Initial conditions
    vx[0] = 1; vy[0] = 0; vz[0] = 1;

}

Trajectory::Trajectory(Trajectory_initializer &init) : Trajectory::Bfield(init){
    do_gen_turb_at_point = init.generate_turbulence;

    t_start = init.t_start; t_end = init.t_end;
    
    E = init.E;
    q = init.q; m = init.m;

    gamma_l = E/(1.0e6 * m);                                //Lorentz factor calculated from E = gamma_l * m * c^2, with [m] = [MeV/c^2]

    unit_coeff = 8.987 * q / (gamma_l * m);                 //Coefficient for the units when dv/dt = unit_coeff * V x B

    double square_root_in_v_tot = std::sqrt(1 - (1 / std::pow(gamma_l, 2)));
    v_tot_init = c * square_root_in_v_tot;                  //Find total initial velocity

    min_err = init.min_err;
    max_err = init.max_err;

    x.resize(2);  y.resize(2);  z.resize(2);
    vx.resize(2); vy.resize(2); vz.resize(2);
    
    D_ij.resize((log10(init.t_end_y) + 1) * 6);

    x[0] = init.x0; y[0] = init.y0; z[0] = init.z0;         
/*
    vx0, vy0 and vz0 are calculated from the energy and given percentages of total initial veloctiy
*/
    vx0 = init.vx0; vy0 = init.vy0; vz0 = init.vz0;
    vx[0] = init.vx0 * v_tot_init; vy[0] = init.vy0 * v_tot_init; vz[0] = init.vz0 * v_tot_init;

    R_Larmor = 1.0810076e-15 * (v_tot_init / c) * (E / init.B_0);
}

Trajectory::~Trajectory(){
    /* Do nothing */
}