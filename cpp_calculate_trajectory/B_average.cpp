#include "class_trajectory.h"
#include "cpp_bfield/class_bfield.h"
#include "cpp_bfield/ran.h"

#include <boost/mpi.hpp>

#include <iostream>
#include <fstream>

#define DX                0.01
#define INSTANCES         2000

#define N_TEST_PARTICLES  4000
#define N_RANDOM_MODES    500
#define T_RUN_FOR_YEARS   1e5
#define B_REGULAR_COMP    0.0                                         //microGauss
#define B_TURBULENT_COMP  4.0                                         //microGauss
#define E_TOTAL           1e15                                        //eV
#define LAMBDA_MAX        10.0                                        //pc
#define LAMBDA_MIN        0.27                                        //pc    (0.27 = Rl/10 for B = 4, E = e16)
#define Q_CHARGE          1                                           //# electron charges
#define M_MASS            938.2720813                                 //MeV/c^2
#define ERROR_MAX         1.0e-01
#define ERROR_MIN         1.0e-08

int initialize_init(Trajectory_initializer &init, int procID);

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  int numProcs;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  int procID;
  MPI_Comm_rank(MPI_COMM_WORLD, &procID);


/******************************************/ 
/************* MAIN PROCESSES *************/
  if(procID == 0){                                                      //Process to finalize calculations
    clock_t begin = clock();

    MPI_Status stat;

    int length_sample = (50*LAMBDA_MAX)/(2 * DX) - 2;

    vector< vector<double> > sumAtPoint(numProcs-1, vector<double>(length_sample, 0));
    vector<double> sumAtPoint_current(length_sample, 0);

    /* Recieve D_ij from all processes */
    for(int i = 1; i < numProcs; i++){                                  //Start from process 1, as this is process 0

      
      MPI_Recv(sumAtPoint_current.data(), length_sample, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &stat);
      
      sumAtPoint[stat.MPI_SOURCE - 1] = sumAtPoint_current;

      std::cout << "I recieved D_ij from process " << stat.MPI_SOURCE << std::endl;
    }//End for i

    for(int i = 0; i < sumAtPoint_current.size(); i++){
        sumAtPoint_current[i] = 0;
    }

    for(int i = 0; i < sumAtPoint.size(); i++){
        for(int j = 0; j < sumAtPoint[i].size(); j++){
            sumAtPoint_current[j] += sumAtPoint[i][j];
        }
    }

    for(int i = 0; i < sumAtPoint_current.size(); i++){
        sumAtPoint_current[i] = sumAtPoint_current[i]/((numProcs - 1)*INSTANCES);
    }

  /* Write <B^2> to a file */
  std::ofstream file;
  file.open("data/B_average.dat");

  if(!file){
    std::cout << "Couldn't open data/eigenvalues.dat. \n";
  }

  for(int i = 0; i < sumAtPoint_current.size(); i++){
    file << sumAtPoint_current[i] <<  '\n';
  }
  std::cout << "<B^2> has been written to data/B_average.dat. \n";
  file.close();

  //This should now have calculated <B^2>.

  clock_t end = clock();
  double elapsed_secs = double(end - begin)/CLOCKS_PER_SEC;

  std::cout << "Rank 0 ended in " << elapsed_secs << " seconds \n";

/*************************************************/    
/************* CALCULATION PROCESSES *************/
  }else{            

    //Settings are now done globally at the top of the program.

    Trajectory_initializer init;

    initialize_init(init, procID);

    Trajectory trajectory(init);
    
    double dx = DX;
    double sum = 0;
    int count = 1;
    int sample_length = 50 * init.lambda_max;
    int instances = INSTANCES;
    double x = 0, y = 0, z = 0;

    vector<double> sumAtPoint(sample_length/(2 * dx) - 2,0);
    
    trajectory.initialize_turbulence();

    for (int i = 0; i < instances; i++){
        int x_add = init.lambda_max * trajectory.ran.doub();
            y = init.lambda_max * trajectory.ran.doub();
            z = init.lambda_max * trajectory.ran.doub();
        for (x = 0; x <= sample_length; x += 2 * dx){

            trajectory.generate_turbulence_at_point(x + x_add, y, z);
            sum += pow(trajectory.turbAtPoint[0], 2) + pow(trajectory.turbAtPoint[1], 2) + pow(trajectory.turbAtPoint[2], 2);
            
            trajectory.generate_turbulence_at_point(x + x_add + dx, y, z);
            sum += 4 * pow(trajectory.turbAtPoint[0], 2) + pow(trajectory.turbAtPoint[1], 2) + pow(trajectory.turbAtPoint[2], 2);
            
            trajectory.generate_turbulence_at_point(x + x_add + 2 * dx, y, z);
            sum += pow(trajectory.turbAtPoint[0], 2) + pow(trajectory.turbAtPoint[1], 2) + pow(trajectory.turbAtPoint[2], 2);

            if(x > 3*dx){
            sumAtPoint[count - 1] += sum * dx / (3 * (x+2*dx));
            count++;
            }
            

        }
        trajectory.reinitialize_turbulence();
        count = 1;
        sum = 0;
        //cout << "Instance " << i << " complete \n";
    }
    
    MPI_Send(sumAtPoint.data(), sumAtPoint.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    
  }
/************* END CALCULATION PROCESSES *************/

  MPI_Finalize();
  return 0;
}




int initialize_init(Trajectory_initializer &init, int procID){
  init.procID = procID;

  std::srand(time(0));  //rand() +                                    //Add this line to seed for "more" randomness
  init.seed = rand() + 1000 * procID;                                 //Sets the seed for the RNG. Set to const 
                                                                      //to generate equal results

  init.n_k = N_RANDOM_MODES;                                          //#modes used to generate TMF
  init.N = N_TEST_PARTICLES;                                          //Number of test-particles
    
  init.t_end_y = T_RUN_FOR_YEARS;                                     //Set time in years
  init.t_start = 0.0;                                                 //Set time start and end here, and the timestep
  init.t_end = 31557600.0 * init.t_end_y; 
  init.dt = 0.1;        

  //(Beck, R. 2003), (Giacinti, Kachelriess & Semikoz, 2012)
  init.B_0 = B_REGULAR_COMP; init.B_rms_turb = B_TURBULENT_COMP;      //B_0 is the regular field, B_rms_turb is the RMS of the turb field

  init.E = E_TOTAL;                                                   //Energy in eV

  init.lambda_max = LAMBDA_MAX; init.lambda_min = LAMBDA_MIN;         //Wavelength in pc

  init.q = Q_CHARGE; init.m = M_MASS;                                 //Set the charge and mass
                                                                        
  init.max_err = ERROR_MAX;                                           //Set min and max error
  init.min_err = ERROR_MIN;
                              
  init.generate_turbulence = true;                                    //This does decide if you generate a turbulence or not
    
    

  init.x0 = 0.0; init.y0 = 0.0; init.z0 = 0.0;                        //Initial conditions
  init.vx0 = 1.0; init.vy0 = 0.0; init.vz0 = 0.0;
}