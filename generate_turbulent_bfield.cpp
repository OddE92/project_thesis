/**********************************************

 * This program is a translation of the structure.f90 code written by Kristian Joten Andersen
 * for his thesis "Charged Particle Trajectories in the Local Superbubble".

 * The program was translated to C++ by Odd-Einar C. Nervik for his project thesis, under the guidance of
 * Michael Kachlerieß.
 
 * Some of the comments in the code are taken directly from the structure.f90-program.
 * The program follows the method shown in:

 *      Charged-particle motion in multidimensional magnetic-field turbulence
 *      Giacalone, J. and Jokipii, J. R.
 *      Department of Planetary Sciences, university of Arizona, Tucson AZ
 *      May 1994

 * The final calculation of the turbulent field follows the equation in:
 
 *      The transport of cosmic rays across a turbulent magnetic field
 *      Giacalone, J. and Jokipii, J. R.
 *      Department of Planetary Sciences, university of Arizona, Tucson AZ
 *      February 1999

 * Which has a different rotation and translation direction than the 1994 article.

 * The output of the program are two files. The first three numbers are the spatial coordinates x, y and z.
 * The last three numbers are Bx, By and Bz respectively.
 
 * The second file is scaled for improved readability when plotted.

**********************************************/
#include "class_bfield.h"

using namespace std;

int main(void){
    double q_max = 10.0;
    double q_min = 0.0;
    int numGridPoints = 10;
    
    Bfield myBfield(q_max, q_min, numGridPoints, true);

    myBfield.initialize_turbulence();
    
    myBfield.write_turbulence_to_files();                  //Write to files also generates the Bfield

    return 0;
}