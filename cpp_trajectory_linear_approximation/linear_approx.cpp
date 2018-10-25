#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;

int main(void){
    int n = 1000000;
    double dt = 0.001;
    double ax, ay, az;
    double q = 1.0;
    double m = 1.0;


    vector<double> x(n), y(n), z(n), vx(n), vy(n), vz(n), bx(n), by(n), bz(n);

    for(int i = 0; i < n; i++){             //uniform field in z-direction
        bx[i] = 0;
        by[i] = 0;
        bz[i] = 1;
    }

    x[0] = 0; y[0] = 1; z[0] = 0;           //initial conditions
    vx[0] = 1.0; vy[0] = 0.0, vz[0] = 10.0;

    for(int i = 0; i < n-1; i++){
        ax = (q/m) * (vy[i] * bz[i] - vz[i] * by[i]);
        ay = (q/m) * (vz[i] * bx[i] - vx[i] * bz[i]);
        az = (q/m) * (vx[i] * by[i] - vy[i] * bx[i]);

        x[i+1] = x[i] + vx[i] * dt;
        y[i+1] = y[i] + vy[i] * dt;
        z[i+1] = z[i] + vz[i] * dt;

        vx[i+1] = vx[i] + ax * dt;
        vy[i+1] = vy[i] + ay * dt;
        vz[i+1] = vz[i] + az * dt;
    }

    ofstream file;

    
    file.open("linear_approx_traj.dat");

    //file << setprecision(20);

    if(!file){
        cout << "Couldn't open linear_approx_traj.dat" << endl;
    }

    for(int i = 0; i < n; i++){
        file << x[i] << ' ' << y[i] << ' ' << z[i] << '\n';
    }


    file.close();


}