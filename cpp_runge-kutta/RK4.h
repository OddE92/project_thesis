#ifndef RUNGEKUTTA4
#define RUNGEKUTTA4

#include <vector>
#include <iostream>
#include <fstream>

#include "cpp_bfield/class_bfield.h"

using std::cout; using std::endl; using std::ofstream;

int RK4(void);

template<class Bfield_func_x, class Bfield_func_y, class Bfield_func_z>
int RK4(Bfield& bfield, Bfield_func_x& Bx_func, Bfield_func_y& By_func, Bfield_func_z& Bz_func, 
        double t_start, double t_end, double dt, double x0, double y0, double z0, double q, double m,
        double vx0, double vy0, double vz0, bool bool_gen_turb)
{

    int n = ceil((t_end - t_start) / dt);               //calculate to AT LEAST the given time
    double t = t_start;
    double qm = q/m;
    double bx, by, bz;                                  //to hold field at the given point
    //double xn1, yn1, zn1, vxn1, vyn1, vz1n;             //to hold current next step value while testing error
    //double dt2;

    vector<double> vx(n), vy(n), vz(n), x(n), y(n), z(n), err(2, 0);

    double kvx1, kvx2, kvx3, kvx4, lvy1, lvy2, lvy3, lvy4, hvz1, hvz2, hvz3, hvz4;
    //vector<double> kvxdt1(4), kvxdt2(7), lvydt1(4), lvydt2(7), hvzdt1(4), hvzdt2(7);

    x[0] = x0; y[0] = y0; z[0] = z0;                    //Initial conditions
    vx[0] = vx0; vy[0] = vy0; vz[0] = vz0;

    for (int i = 0; i < n; i++){

        //Generate the field at the new point.
        bfield.generate_bfield_at_point(Bx_func, By_func, Bz_func, t, bx, by, bz, x[i], y[i], z[i]);
 

        //Standard implementation of a RK4-method for 3 coupled ODEs.
        kvx1 = dt * qm * (vy[i] * bz - vz[i] * by);
        lvy1 = dt * qm * (vz[i] * bx - vx[i] * bz);
        hvz1 = dt * qm * (vx[i] * by - vy[i] * bx);

        kvx2 = dt * qm * ((vy[i] + lvy1 / 2) * bz - (vz[i] + hvz1 / 2) * by);
        lvy2 = dt * qm * ((vz[i] + hvz1 / 2) * bx - (vx[i] + kvx1 / 2) * bz);
        hvz2 = dt * qm * ((vx[i] + kvx1 / 2) * by - (vy[i] + lvy1 / 2) * bx);

        kvx3 = dt * qm * ((vy[i] + lvy2 / 2) * bz - (vz[i] + hvz2 / 2) * by);
        lvy3 = dt * qm * ((vz[i] + hvz2 / 2) * bx - (vx[i] + kvx2 / 2) * bz);
        hvz3 = dt * qm * ((vx[i] + kvx2 / 2) * by - (vy[i] + lvy2 / 2) * bx);

        kvx4 = dt * qm * ((vy[i] + lvy3) * bz - (vz[i] + hvz3) * by);
        lvy4 = dt * qm * ((vz[i] + hvz3) * bx - (vx[i] + kvx3) * bz);
        hvz4 = dt * qm * ((vx[i] + kvx3) * by - (vy[i] + lvy3) * bx);

        x[i + 1] = x[i] + vx[i] * dt;
        y[i + 1] = y[i] + vy[i] * dt;
        z[i + 1] = z[i] + vz[i] * dt;

        vx[i + 1] = vx[i] + (1.0 / 6) * (kvx1 + 2 * kvx2 + 2 * kvx3 + kvx4);
        vy[i + 1] = vy[i] + (1.0 / 6) * (lvy1 + 2 * lvy2 + 2 * lvy3 + lvy4);
        vz[i + 1] = vz[i] + (1.0 / 6) * (hvz1 + 2 * hvz2 + 2 * hvz3 + hvz4);

        t += dt;
    }

    ofstream file;

    
    file.open("cpp_runge-kutta/RK4_approx.dat");

    //file << setprecision(20);

    if(!file){
        cout << "Couldn't open RK4_approx.dat" << endl;
    }

    for(int i = 0; i < n; i++){
        file << x[i] << ' ' << y[i] << ' ' << z[i] << '\n';
    }


    file.close();
    
    return 0;

}

#endif