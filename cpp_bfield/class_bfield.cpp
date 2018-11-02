#include <vector>
#include <complex>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <time.h>
#include <random>

#include "class_bfield.h"

using namespace std;

/********************************/
/********** GET BFIELD **********/
int Bfield::get_bfield(double& bx, double& by, double& bz, const double x, const double y, const double z){
    //Calculates the bfield at the given position, using the nearest neighbours.

    if(x > xyz_max || x < xyz_min || y > xyz_max || y < xyz_min || z > xyz_max || z < xyz_min ){
        cout << "Position outside of generated field. Please expand field." << endl;
        return 1;
    }

    double bx1, bx2, by1, by2, bz1, bz2, x1, x2, y1, y2, z1, z2, x1pos, x2pos, y1pos, y2pos, z1pos, z2pos;

    /* 
    if n*cscal = x/y/z:
    z: i = (z - xyz_min) / cscal
    y: i = c_size * (y - xyz_min) / cscal
    x: i = c2 * (x - xyz_min) / cscal
    
    Else we need to round off x/y/z (general case):
        n = x/cscal,
        floor(n) * cscal = x_min
        ceil (n) * cscal = x_max

    This should guarantee a whole number
    */

    z1 = floor(z/cscal) * cscal;                        // Finds the closest z in the array, given an arbitrary z-value
    z2 = ceil (z/cscal) * cscal;
    y1 = floor(y/cscal) * cscal;
    y2 = ceil (y/cscal) * cscal;
    x1 = floor(x/cscal) * cscal;
    x2 = ceil (x/cscal) * cscal; 

    z1pos = (z1 - xyz_min) / cscal;                       // Finds the position in the array of the given z-value
    z2pos = (z2 - xyz_min) / cscal;
    y1pos = c_size * (y1 - xyz_min) / cscal;
    y2pos = c_size * (y2 - xyz_min) / cscal;
    x1pos = c2 * (x1 - xyz_min) / cscal;
    x2pos = c2 * (x2 - xyz_min) / cscal;
    
    bx1 = B[x1pos + y1pos + z1pos][0];                  // Gets the field at the 6 nearest neighbours to the given x/y/z
    by1 = B[x1pos + y1pos + z1pos][1];
    bz1 = B[x1pos + y1pos + z1pos][2];
    bx2 = B[x2pos + y2pos + z1pos][0];
    by2 = B[x2pos + y2pos + z1pos][1];
    bz2 = B[x2pos + y2pos + z1pos][2];

    bx = (bx2 - bx1)/(x2 - x1) * (x - x1) + bx1;        // linear approximation of the field at x/y/z
    by = (by2 - bx1)/(y2 - y1) * (y - y1) + by1;
    bz = (bz2 - bx1)/(z2 - z1) * (z - z1) + bz1;

    return 0;
}
/******** END GET BFIELD ********/

/********************************/
/****** GENERATE TURBULENCE *****/
void Bfield::generate_turbulence(int x){
    vector<double> e_k(n_k);
    
    vector< complex<double> > delta_B(3);
    vector< complex<double> > F(n_k), epsilon_x(n_k), epsilon_y(n_k), epsilon_z(n_k), B_x_k(n_k), B_y_k(n_k), B_z_k(n_k);

    // This function takes a position x, and generate a "random" turbulence to the B-field

    // Follow the articles calculation of delta_B (or delta_Omega, where delta_B = delta_Omega * mc/q)

    for(int i = 0; i < n_k; i++){

        // e_k is the exponent of delta_B
        e_k[i] = k[i] * (st[i]*cp[i]*coord[x][0] + st[i]*sp[i]*coord[x][1] + ct[i]*coord[x][2]) + b[i];

        // F = B(k)*e^i(kz+beta), B_k = B(k)
        F[i].real(B_k[i] * cos(e_k[i]));  F[i].imag(B_k[i] * sin(e_k[i]));

        //Epsilon is a vector, so it is split into epsilon_(x/y/z), where (x/y/z) is the direction in the non-rotated coord.sys.
        epsilon_x[i].real( ca[i] * ct[i] * cp[i]);          epsilon_x[i].imag(-s[i] * sa[i] * sp[i]);
        epsilon_y[i].real( ca[i] * ct[i] * sp[i]);          epsilon_y[i].imag( s[i] * sa[i] * cp[i]);
        epsilon_z[i].real(-ca[i] * st[i]);                  epsilon_z[i].imag(0.0);

        //Calculate B(k) in the x-, y- and z-directions
        B_x_k[i] = F[i] * epsilon_x[i];                      // Complex number multiplication!
        B_y_k[i] = F[i] * epsilon_y[i];
        B_z_k[i] = F[i] * epsilon_z[i];
    }

    for(int i = 0; i < delta_B.size(); i++){                 // Need to set delta_B to 0 before we start summing
        delta_B[i].imag(0.0); delta_B[i].real(0.0);
    }

        //Calculates the final delta_B in the x-, y- and z-directions.
    for(int i = 0; i < n_k; i++){
        delta_B[0] += B_x_k[i];
        delta_B[1] += B_y_k[i];
        delta_B[2] += B_z_k[i];
    }

// As we are concerned about the real value of the B-field we need to calculate the real part of delta_B

    for(int i = 0; i < delta_B.size(); i++){
        turbulencevec[x][i] = real(delta_B[i]);
    }

// B = m*c*Omega/q, Omega = Omega_0 + dOmega. Check units!
}

int Bfield::generate_turbulence_at_point(double x, double y, double z){
    vector<double> e_k(n_k);
    
    vector< complex<double> > delta_B(3);
    vector< complex<double> > F(n_k), epsilon_x(n_k), epsilon_y(n_k), epsilon_z(n_k), B_x_k(n_k), B_y_k(n_k), B_z_k(n_k);

    // This function takes a position x, and generate a "random" turbulence to the B-field

    // Follow the articles calculation of delta_B (or delta_Omega, where delta_B = delta_Omega * mc/q)

    for(int i = 0; i < n_k; i++){

        // e_k is the exponent of delta_B
        e_k[i] = k[i] * (st[i]*cp[i]*x + st[i]*sp[i]*y + ct[i]*z) + b[i];

        // F = B(k)*e^i(kz+beta), B_k = B(k)
        F[i].real(B_k[i] * cos(e_k[i]));  F[i].imag(B_k[i] * sin(e_k[i]));

        //Epsilon is a vector, so it is split into epsilon_(x/y/z), where (x/y/z) is the direction in the non-rotated coord.sys.
        epsilon_x[i].real( ca[i] * ct[i] * cp[i]);          epsilon_x[i].imag(-s[i] * sa[i] * sp[i]);
        epsilon_y[i].real( ca[i] * ct[i] * sp[i]);          epsilon_y[i].imag( s[i] * sa[i] * cp[i]);
        epsilon_z[i].real(-ca[i] * st[i]);                  epsilon_z[i].imag(0.0);

        //Calculate B(k) in the x-, y- and z-directions
        B_x_k[i] = F[i] * epsilon_x[i];                      // Complex number multiplication!
        B_y_k[i] = F[i] * epsilon_y[i];
        B_z_k[i] = F[i] * epsilon_z[i];
    }

    for(int i = 0; i < delta_B.size(); i++){                 // Need to set delta_B to 0 before we start summing
        delta_B[i].imag(0.0); delta_B[i].real(0.0);
    }

        //Calculates the final delta_B in the x-, y- and z-directions.
    for(int i = 0; i < n_k; i++){
        delta_B[0] += B_x_k[i];
        delta_B[1] += B_y_k[i];
        delta_B[2] += B_z_k[i];
    }

// As we are concerned about the real value of the B-field we need to calculate the real part of delta_B

    for(int i = 0; i < delta_B.size(); i++){
        turbAtPoint[i] = real(delta_B[i]);
    }

// B = m*c*Omega/q, Omega = Omega_0 + dOmega. Check units!
    return 0;
}
/**** END GENERATE TURBULENCE ***/


/********************************/
/******** WRITE TO FILES ********/

void Bfield::write_turbulence_to_files(){

    string fileName;

    //cout << "Enter filename (remember \".dat\"): ";
    //cin >> fileName;

    fileName = "turbulence.dat";

    fileName = "cpp_bfield/" + fileName;

    ofstream file1, file2;

    file1.open(fileName);
    file2.open("cpp_bfield/turbulence_scaled.dat");

    file1 << setprecision(20);
    file2 << setprecision(20);

    if(!file1){
        cout << "Couldn't open " << fileName << endl;
    }
    if(!file2){
        cout << "Couldn't open turbulence_scaled.dat" << endl;
    }

    for(int i = 0; i < c3; i++){
        generate_turbulence(i);                                         // Generates the next line of the Bfield

        for(int j = 0; j < turbulencevec[i].size(); j++){               // Scales the Bfield for plotting
            Bvec_scal[j] = scal * turbulencevec[i][j];
        }


        for(int j = 0; j < coord[i].size(); j++){                       // First write the 3 coordinates
            file1 << coord[i][j] << ' ';
            file2 << coord[i][j] << ' ';
        }
        
        for(int j = 0; j < turbulencevec[i].size(); j++){               // Then write the Bfield components
            file1 << turbulencevec[i][j] << ' ';
            file2 << Bvec_scal[j] << ' ';
        }

        file1 << '\n';                                              
        file2 << '\n';

    }

    file1.close();
    file2.close();
}

void Bfield::write_B_to_file(){

    if(!generate_bfield_is_run){
        cout << "Field has not been generated. Run generate_bfield." << endl;
        return;
    }

    ofstream file;

    file.open("Total_B_field.dat");

    file << setprecision(20);

    if(!file){
        cout << "Couldn't open Total_B_field.dat" << endl;
    }

    for(int i = 0; i < c3; i++){

        for(int j = 0; j < coord[i].size(); j++){                       // First write the 3 coordinates
            file << coord[i][j] << ' ';
        }
        
        for(int j = 0; j < turbulencevec[i].size(); j++){               // Then write the Bfield components
            file << B[i][j] << ' ';
        }

        file << '\n';

    }

    file.close();
}


/********** INITIALIZE **********/

void Bfield::initialize_turbulence(){
    if (!turbulence_is_initialized){
        initialize_phases();
        initialize_normalization();
        initialize_coords();

        turbulence_is_initialized = true;
    }
}

/********************************/
/****** INITIALIZE PHASES *******/

void Bfield::initialize_phases(){
    double r;

    for(int i = 0; i < n_k; i++){                           // From the article:
        a[i] = two_pi * ran.doub();                // a = alpha
        b[i] = two_pi * ran.doub();                // b = beta
        p[i] = two_pi * ran.doub();                // p = phi
        t[i] = M_PI * ran.doub();                  // t = theta

        s[i] = 1.0;                                         // s = sign in epsilon (+-)
        r = ran.doub();
        if(r > 0.5) s[i] = -1.0;

        cp[i] = cos(p[i]);                                  //sines and cosines of theta, phi and alpha
        sp[i] = sin(p[i]);
        ct[i] = cos(t[i]);
        st[i] = sin(t[i]);
        ca[i] = cos(a[i]);
        sa[i] = cos(a[i]);
    }

}
/***** END INITIALIZE PHASES ****/
/********************************/

/********************************/
/*** INITIALIZE NORMALIZATION ***/

void Bfield::initialize_normalization(){

    double logk_diff, dlog_k;

    logk_diff = log10(k_max) - log10(k_min);
    dlog_k = logk_diff/(n_k - 1);


    //The for loop spreads the values of k equally along the logarithmic scale.
    for(int i = 0; i < n_k; i++){
        k[i] = pow(10, log10(k_min) + (i-1) * dlog_k);
        //B_k[i] = 1; to be used while testing.
    }

/*
    Normalization: The energy density S should stay the same for any number of wavevectors (n_k) used.
    Energy density is given by:
    S = sum( B(k)^2 / 8pi )

    The total energy density can be normalized by setting S = (B_0)^2 / 8pi, which gives

    (B_0)^2 / 8pi = B(k_min)^2 / 8pi * sum[ (k / k_min)^-gamma ]; B(k_min) = dB_min

    B(k_min) = sqrt{ (B_0)^2 / sum[ (k/k_min)^-gamma ]}

    Remember to check units against giacalone & jokipii 1994, where they use Omega instead of B.
*/
  
  double sum = 0;
  for(int i = 0; i < n_k; i++){
      sum += pow(k[i]/k_min, -gamma);
  }
  dB_min = sqrt(B_rms_turb / sum);

  for(int i = 0; i < n_k; i++){
      B_k[i] = dB_min * pow(k[i] / k_min, -gamma/2);
  }

  //Normalization is now done so B_0 = 1, not Omega(k) as in Giacalone & jokipii 1994.  
}
/****** END INITIALIZE NORM *****/
/********************************/

/********************************/
/**** INITIALIZE COORDINATES ****/

void Bfield::initialize_coords(){

    for(int i = 0; i < c_size; i++){                    //Sets the coordinates such that z increases from 0-(c_size-1)
        for(int j = 0; j < c_size; j++){                //y increases with 1 after each z-loop, and
            for(int m = 0; m < c_size; m++){            //x increases after y reaches c_size-1
                
                coord[i*c2 + j*c_size + m][0] =  i * cscal + xyz_min;
                coord[i*c2 + j*c_size + m][1] =  j * cscal + xyz_min;
                coord[i*c2 + j*c_size + m][2] =  m * cscal + xyz_min;

            }
        }
    }

}
/*** END INITIALIZE COORDINATES ***/
/**********************************/

Bfield::Bfield() : ran(15321){
    
    // See class declaration for vector descriptions.
    c2 = pow(c_size, 2), c3 = pow(c_size, 3);

    B.resize(c3);
    B0vec.resize(c3);
    turbulencevec.resize(c3);
    Bvec_scal.resize(3);
    turbAtPoint.resize(3);
    
    coord.resize(c3);

    for(int i = 0; i < c3; i++){
        B[i].resize(3);
        B0vec[i].resize(3);
        turbulencevec[i].resize(3);
        coord[i].resize(3);
    }

    im.imag(1.0); im.real(0.0);

    a.resize(n_k); b.resize(n_k); p.resize(n_k); t.resize(n_k); s.resize(n_k);              
    ca.resize(n_k); sa.resize(n_k); cp.resize(n_k); sp.resize(n_k); ct.resize(n_k); 
    st.resize(n_k);                                                                         

    B_k.resize(n_k); k.resize(n_k);
}

Bfield::Bfield( double xyz_max, double xyz_min, int numGridPoints, bool Bgenerate_turbulence) : xyz_max(xyz_max), xyz_min(xyz_min),
                c_size(numGridPoints + 1), ran(15321){

    // See class declaration for vector descriptions.
    c2 = pow(c_size, 2), c3 = pow(c_size, 3);

    B.resize(c3);
    B0vec.resize(c3);
    turbulencevec.resize(c3);
    Bvec_scal.resize(3);  
    coord.resize(c3);
    turbAtPoint.resize(3);

    for(int i = 0; i < c3; i++){
        B[i].resize(3);
        B0vec[i].resize(3);
        turbulencevec[i].resize(3);
        coord[i].resize(3);
    }

    im.imag(1.0); im.real(0.0);

    a.resize(n_k); b.resize(n_k); p.resize(n_k); t.resize(n_k); s.resize(n_k);              
    ca.resize(n_k); sa.resize(n_k); cp.resize(n_k); sp.resize(n_k); ct.resize(n_k); 
    st.resize(n_k);                                                                         

    B_k.resize(n_k); k.resize(n_k);

    cscal = (xyz_max - xyz_min) / numGridPoints;

}

Bfield::Bfield( Trajectory_initializer &init) : xyz_max(init.xyz_max), xyz_min(xyz_min),
                c_size(init.numGridPoints + 1), B_0(init.B_0), B_rms_turb(init.B_rms_turb), ran(15321 + init.seed){

    // See class declaration for vector descriptions.
    c2 = pow(c_size, 2), c3 = pow(c_size, 3);

    B.resize(c3);
    B0vec.resize(c3);
    turbulencevec.resize(c3);
    Bvec_scal.resize(3);  
    coord.resize(c3);
    turbAtPoint.resize(3);

    for(int i = 0; i < c3; i++){
        B[i].resize(3);
        B0vec[i].resize(3);
        turbulencevec[i].resize(3);
        coord[i].resize(3);
    }

    im.imag(1.0); im.real(0.0);

    a.resize(n_k); b.resize(n_k); p.resize(n_k); t.resize(n_k); s.resize(n_k);              
    ca.resize(n_k); sa.resize(n_k); cp.resize(n_k); sp.resize(n_k); ct.resize(n_k); 
    st.resize(n_k);                                                                         

    B_k.resize(n_k); k.resize(n_k);

    cscal = (init.xyz_max - init.xyz_min) / init.numGridPoints;

}

Bfield::~Bfield(){
    //do nothing
}
