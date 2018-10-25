#include <vector>
#include <iostream>
#include <fstream>

#include "RK4.h"

using namespace std;

/********************************/
/* No input-version for testing */

int RK4(void){
    int n = 1000000;
    double dt = 0.1, t = 0.0;
    double q = 1.0, m = 1.0;
    double qm = q/m;

    vector<double> vx(n), vy(n), vz(n), x(n), y(n), z(n), bx(n, 0), by(n, 0), bz(n, 1);

    double kvx1, kvx2, kvx3, kvx4, lvy1, lvy2, lvy3, lvy4, hvz1, hvz2, hvz3, hvz4;
    double kx1, kx2, kx3, kx4, ly1, ly2, ly3, ly4, hz1, hz2, hz3, hz4;

    vx[0] = 1.0; vy[0] = 0.0; vz[0] = 10.0;

    for(int i = 0; i < n; i++){
        /*
        kx1 = x[i] / t;
        kx2 = (x[i] + kx1/2) / (t + dt/2);
        kx3 = (x[i] + kx2/2) / (t + dt/2);
        kx4 = (x[i] + kx3) / (t + dt);

        ly1 = y[i] / t;
        ly2 = (y[i] + ly1/2) / (t + dt/2);
        ly3 = (y[i] + ly2/2) / (t + dt/2);
        ly4 = (y[i] + ly3) / (t + dt);

        hz1 = z[i] / t;
        hz2 = (z[i] + hz1/2) / (t + dt/2);
        hz3 = (z[i] + hz2/2) / (t + dt/2);
        hz4 = (z[i] + hz3) / (t + dt);
        */

        kvx1 = dt * qm * (vy[i] * bz[i] - vz[i] * by[i]);
        lvy1 = dt * qm * (vz[i] * bx[i] - vx[i] * bz[i]);
        hvz1 = dt * qm * (vx[i] * by[i] - vy[i] * bx[i]);

        kvx2 = dt * qm * ((vy[i] + lvy1/2) * bz[i] - (vz[i] + hvz1/2) * by[i]);
        lvy2 = dt * qm * ((vz[i] + hvz1/2) * bx[i] - (vx[i] + kvx1/2) * bz[i]);
        hvz2 = dt * qm * ((vx[i] + kvx1/2) * by[i] - (vy[i] + lvy1/2) * bx[i]);

        kvx3 = dt * qm * ((vy[i] + lvy2/2) * bz[i] - (vz[i] + hvz2/2) * by[i]);
        lvy3 = dt * qm * ((vz[i] + hvz2/2) * bx[i] - (vx[i] + kvx2/2) * bz[i]);
        hvz3 = dt * qm * ((vx[i] + kvx2/2) * by[i] - (vy[i] + lvy2/2) * bx[i]);

        kvx4 = dt * qm * ((vy[i] + lvy3) * bz[i] - (vz[i] + hvz3) * by[i]);
        lvy4 = dt * qm * ((vz[i] + hvz3) * bx[i] - (vx[i] + kvx3) * bz[i]);
        hvz4 = dt * qm * ((vx[i] + kvx3) * by[i] - (vy[i] + lvy3) * bx[i]);

        /*
        x[i+1] = x[i] + (dt / 6) * (kx1 + 2*kx2 + 2*kx3 + kx4);
        y[i+1] = y[i] + (dt / 6) * (ly1 + 2*ly2 + 2*ly3 + ly4);
        z[i+1] = z[i] + (dt / 6) * (hz1 + 2*hz2 + 2*hz3 + hz4);
        */

        x[i+1] = x[i] + vx[i] * dt;
        y[i+1] = y[i] + vy[i] * dt;
        z[i+1] = z[i] + vz[i] * dt;

        vx[i+1] = vx[i] + (1.0/6) * (kvx1 + 2*kvx2 + 2*kvx3 + kvx4);
        vy[i+1] = vy[i] + (1.0/6) * (lvy1 + 2*lvy2 + 2*lvy3 + lvy4);
        vz[i+1] = vz[i] + (1.0/6) * (hvz1 + 2*hvz2 + 2*hvz3 + hvz4);

        t += dt;
    }

        ofstream file;

    
    file.open("RK4_approx.dat");

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