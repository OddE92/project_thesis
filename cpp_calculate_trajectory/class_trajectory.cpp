#include "class_trajectory.h"
#include "cpp_bfield/class_bfield.h"

#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>

using std::cout; using std::endl; using std::ofstream;

/******************************************/
/************ RK4 STEP SIZE CTRL **********/
int Trajectory::RK4t_step_size_control(){

    t = t_start;
    double bx, by, bz;                                      //To hold field at the given point
    double dt2 = dt/2;                                      //Step size for size controlled part
    double truncation_error;                                //To hold truncation error
    double total_velocity_dt1, total_velocity_dt2;          //To calculate truncation error
    int nbad = 0, ntoolow = 0;                              //Counters to check number of steps where dt changes

    vector<double> kvxdt1(4), kvxdt2(4), lvydt1(4), lvydt2(4), hvzdt1(4), hvzdt2(4);    //Runge-kutta coefficients
    double vxdt1, vydt1, vzdt1, vxdt2, vydt2, vzdt2;        //to hold speed and check error

    x.reserve(2*n); y.reserve(2*n); z.reserve(2*n);         //Position and velocity vectors
    vx.reserve(2*n); vy.reserve(2*n); z.reserve(2*n);

    x.resize(1); y.resize(1); z.resize(1);                  //set the first element to zero, but the
    vx.resize(1); vy.resize(1); vz.resize(1);               //vector capacity is still kept as 2*n

    ofstream file;

    
    file.open("cpp_runge-kutta/RK4_approx_stepsizectrl.dat");

    //file << setprecision(20);

    if(!file){
        cout << "Couldn't open RK4_approx_stepsizectrl.dat" << endl;
    }

    //Generate the bfield at the starting point
    generate_bfield_at_point(Bx, By, Bz, t, bx, by, bz, x[0]/mtopc, y[0]/mtopc, z[0]/mtopc, do_gen_turb_at_point);

    for (int i = 0; t <= t_end; i++){
        dt2 = dt / 2;
        //Generate the field at the new point.
        //generate_bfield_at_point(Bx, By, Bz, t, bx, by, bz, x[i], y[i], z[i]);

        //Standard implementation of a RK4-method for 3 coupled ODEs.
        kvxdt1[0] = dt * unit_coeff * (vy[i] * bz - vz[i] * by);
        lvydt1[0] = dt * unit_coeff * (vz[i] * bx - vx[i] * bz);
        hvzdt1[0] = dt * unit_coeff * (vx[i] * by - vy[i] * bx);
        kvxdt2[0] = dt2 * unit_coeff * (vy[i] * bz - vz[i] * by);
        lvydt2[0] = dt2 * unit_coeff * (vz[i] * bx - vx[i] * bz);
        hvzdt2[0] = dt2 * unit_coeff * (vx[i] * by - vy[i] * bx);

        kvxdt1[1] = dt * unit_coeff * ((vy[i] + lvydt1[0] / 2) * bz - (vz[i] + hvzdt1[0] / 2) * by);
        lvydt1[1] = dt * unit_coeff * ((vz[i] + hvzdt1[0] / 2) * bx - (vx[i] + kvxdt1[0] / 2) * bz);
        hvzdt1[1] = dt * unit_coeff * ((vx[i] + kvxdt1[0] / 2) * by - (vy[i] + lvydt1[0] / 2) * bx);
        kvxdt2[1] = dt2 * unit_coeff * ((vy[i] + lvydt2[0] / 2) * bz - (vz[i] + hvzdt2[0] / 2) * by);
        lvydt2[1] = dt2 * unit_coeff * ((vz[i] + hvzdt2[0] / 2) * bx - (vx[i] + kvxdt2[0] / 2) * bz);
        hvzdt2[1] = dt2 * unit_coeff * ((vx[i] + kvxdt2[0] / 2) * by - (vy[i] + lvydt2[0] / 2) * bx);

        kvxdt1[2] = dt * unit_coeff * ((vy[i] + lvydt1[1] / 2) * bz - (vz[i] + hvzdt1[1] / 2) * by);
        lvydt1[2] = dt * unit_coeff * ((vz[i] + hvzdt1[1] / 2) * bx - (vx[i] + kvxdt1[1] / 2) * bz);
        hvzdt1[2] = dt * unit_coeff * ((vx[i] + kvxdt1[1] / 2) * by - (vy[i] + lvydt1[1] / 2) * bx);
        kvxdt2[2] = dt2 * unit_coeff * ((vy[i] + lvydt2[1] / 2) * bz - (vz[i] + hvzdt2[1] / 2) * by);
        lvydt2[2] = dt2 * unit_coeff * ((vz[i] + hvzdt2[1] / 2) * bx - (vx[i] + kvxdt2[1] / 2) * bz);
        hvzdt2[2] = dt2 * unit_coeff * ((vx[i] + kvxdt2[1] / 2) * by - (vy[i] + lvydt2[1] / 2) * bx);

        kvxdt1[3] = dt * unit_coeff * ((vy[i] + lvydt1[2]) * bz - (vz[i] + hvzdt1[2]) * by);
        lvydt1[3] = dt * unit_coeff * ((vz[i] + hvzdt1[2]) * bx - (vx[i] + kvxdt1[2]) * bz);
        hvzdt1[3] = dt * unit_coeff * ((vx[i] + kvxdt1[2]) * by - (vy[i] + lvydt1[2]) * bx);
        kvxdt2[3] = dt2 * unit_coeff * ((vy[i] + lvydt2[2]) * bz - (vz[i] + hvzdt2[2]) * by);
        lvydt2[3] = dt2 * unit_coeff * ((vz[i] + hvzdt2[2]) * bx - (vx[i] + kvxdt2[2]) * bz);
        hvzdt2[3] = dt2 * unit_coeff * ((vx[i] + kvxdt2[2]) * by - (vy[i] + lvydt2[2]) * bx);

/***** velocity calculation *****/
        vxdt1 = vx[i] + (1.0 / 6) * (kvxdt1[0] + 2 * (kvxdt1[1] + kvxdt1[2]) + kvxdt1[3]);
        vydt1 = vy[i] + (1.0 / 6) * (lvydt1[0] + 2 * (lvydt1[1] + lvydt1[2]) + lvydt1[3]);
        vzdt1 = vz[i] + (1.0 / 6) * (hvzdt1[0] + 2 * (hvzdt1[1] + hvzdt1[2]) + hvzdt1[3]);
        vxdt2 = vx[i] + (1.0 / 6) * (kvxdt2[0] + 2 * (kvxdt2[1] + kvxdt2[2]) + kvxdt2[3]);
        vydt2 = vy[i] + (1.0 / 6) * (lvydt2[0] + 2 * (lvydt2[1] + lvydt2[2]) + lvydt2[3]);
        vzdt2 = vz[i] + (1.0 / 6) * (hvzdt2[0] + 2 * (hvzdt2[1] + hvzdt2[2]) + hvzdt2[3]);
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
            //Calculate the next velocity
            vx.push_back(vxdt2);
            vy.push_back(vydt2);
            vz.push_back(vzdt2);

            //Calculate the next position
            x.push_back(x[i] + mtopc * vx[i] * dt);                         //Adjust for parsec instead of m
            y.push_back(y[i] + mtopc * vy[i] * dt);
            z.push_back(z[i] + mtopc * vz[i] * dt);

            file << x[i] << ' ' << y[i] << ' ' << z[i] << '\n';             //write position to file

            t += dt;       

            //Generate the field at the new point.
            generate_bfield_at_point(Bx, By, Bz, t, bx, by, bz, x[i+1]/mtopc, y[i+1]/mtopc, z[i+1]/mtopc, do_gen_turb_at_point);

        }else if(truncation_error > max_err){           //If the truncation error is too large, we need to decrease the stepsize
            dt = dt2;
            i--;                                        //decrease i to redo the step
            nbad++;

        }else if(truncation_error < min_err){           //If the truncation error is too small, take a larger step next time.
            //Calculate the next velocity
            vx.push_back(vxdt2);
            vy.push_back(vydt2);
            vz.push_back(vzdt2);

            //Calculate the next position
            x.push_back(x[i] + mtopc * vx[i] * dt);                         //Adjust for parsec instead of m
            y.push_back(y[i] + mtopc * vy[i] * dt);
            z.push_back(z[i] + mtopc * vz[i] * dt);

            file << x[i] << ' ' << y[i] << ' ' << z[i] << '\n';             //write position to file

            t += dt; 
            dt = 2*dt;                                                      //double the step size
            ntoolow++;

            //Generate the field at the new point.
            generate_bfield_at_point(Bx, By, Bz, t, bx, by, bz, x[i+1]/mtopc, y[i+1]/mtopc, z[i+1]/mtopc, do_gen_turb_at_point);
              
            
        }else{
            cout << "Something unexpected happened. Truncation error tests weren't hit. \n";
            return 1;
        }

        cout << truncation_error << '\n';

    }//End for i

    file.close();
    
    cout << "Number of steps taken: " << x.size() << endl;
    cout << "Number of bad steps taken: " << nbad << endl;
    cout << "Number of steps with too high accuracy: " << ntoolow << endl;
    cout << "Initial velocity: " << v_tot_init << " m/s" << endl;
    cout << "vx0: " << vx[0] << " ; vy0: " << vy[0] << " ; vz0: " << vz[0] << " m/s" << endl;
    cout << "Larmor radius: " << R_Larmor << " parsec" << endl;

    return 0;
}

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

        generate_bfield_at_point(Bx, By, Bz, t, bx, by, bz, x[i]/mtopc, y[i]/mtopc, z[i]/mtopc, do_gen_turb_at_point);
        

        //Standard implementation of a RK4-method for 3 coupled ODEs.
        kvx1 = dt * unit_coeff * (vy[i] * bz - vz[i] * by);
        lvy1 = dt * unit_coeff * (vz[i] * bx - vx[i] * bz);
        hvz1 = dt * unit_coeff * (vx[i] * by - vy[i] * bx);

        kvx2 = dt * unit_coeff * ((vy[i] + lvy1 / 2) * bz - (vz[i] + hvz1 / 2) * by);
        lvy2 = dt * unit_coeff * ((vz[i] + hvz1 / 2) * bx - (vx[i] + kvx1 / 2) * bz);
        hvz2 = dt * unit_coeff * ((vx[i] + kvx1 / 2) * by - (vy[i] + lvy1 / 2) * bx);

        kvx3 = dt * unit_coeff * ((vy[i] + lvy2 / 2) * bz - (vz[i] + hvz2 / 2) * by);
        lvy3 = dt * unit_coeff * ((vz[i] + hvz2 / 2) * bx - (vx[i] + kvx2 / 2) * bz);
        hvz3 = dt * unit_coeff * ((vx[i] + kvx2 / 2) * by - (vy[i] + lvy2 / 2) * bx);

        kvx4 = dt * unit_coeff * ((vy[i] + lvy3) * bz - (vz[i] + hvz3) * by);
        lvy4 = dt * unit_coeff * ((vz[i] + hvz3) * bx - (vx[i] + kvx3) * bz);
        hvz4 = dt * unit_coeff * ((vx[i] + kvx3) * by - (vy[i] + lvy3) * bx);

        //Calculate next position
        x[i + 1] = x[i] + mtopc * vx[i] * dt;
        y[i + 1] = y[i] + mtopc * vy[i] * dt;
        z[i + 1] = z[i] + mtopc * vz[i] * dt;

        //Calculate next velocity
        vx[i + 1] = vx[i] + (1.0 / 6) * (kvx1 + 2 * kvx2 + 2 * kvx3 + kvx4);
        vy[i + 1] = vy[i] + (1.0 / 6) * (lvy1 + 2 * lvy2 + 2 * lvy3 + lvy4);
        vz[i + 1] = vz[i] + (1.0 / 6) * (hvz1 + 2 * hvz2 + 2 * hvz3 + hvz4);

        file << x[i] << ' ' << y[i] << ' ' << z[i] << '\n';

        t += dt;
    }

    file.close();
    
    return 0;

}
/********  END RK4 NO STEP SIZE CTRL ******/

/******************************************/
/************** CONSTRUCTORS **************/

Trajectory::Trajectory() : Trajectory::Bfield(){
    n = ceil((t_end - t_start) / dt);
    unit_coeff = q/m;
    max_err = 1.0e-05;
    min_err = 1.0e-09;

    x.resize(n);  y.resize(n);  z.resize(n);
    vx.resize(n); vy.resize(n); vz.resize(n);

    x[0] = 0; y[0] = 0; z[0] = 0;                           //Initial conditions
    vx[0] = 1; vy[0] = 0; vz[0] = 1;

}

Trajectory::Trajectory(Trajectory_initializer &init) : Trajectory::Bfield(init){
    do_gen_turb_at_point = init.generate_turbulence;

    n = ceil((init.t_end - init.t_start) / init.dt);
    
    E = init.E;
    q = init.q; m = init.m;

    gamma_l = (m * c * c) / E;                              //Lorentz factor calculated from E = gamma_l * m * c^2

    unit_coeff = 8.987 * q / (gamma_l * m);                 //Coefficient for the units when dv/dt = unit_coeff * V x B

    double square_root_in_v_tot = std::sqrt(1 - (1 / std::pow(gamma_l, 2)));
    v_tot_init = c * square_root_in_v_tot;                  //Find total initial velocity

    min_err = init.min_err;
    max_err = init.max_err;

    x.resize(n);  y.resize(n);  z.resize(n);
    vx.resize(n); vy.resize(n); vz.resize(n);

    x[0] = init.x0; y[0] = init.y0; z[0] = init.z0;         
/*
    vx0, vy0 and vz0 are calculated from the energy and given percentages of total initial veloctiy
*/
    vx[0] = init.vx0 * v_tot_init; vy[0] = init.vy0 * v_tot_init; vz[0] = init.vz0 * v_tot_init;

    R_Larmor = 1.0810076e-15 * (v_tot_init / c) * (E / init.B_0);
}

Trajectory::~Trajectory(){
    /* Do nothing */
}