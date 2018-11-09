#include "cpp_calculate_trajectory/class_trajectory.h"
#include "cpp_bfield/class_bfield.h"
#include "cpp_bfield/ran.h"
#include "calculate_eigenvalues.h"

#include <boost/mpi.hpp>

#include <iostream>
#include <fstream>

#define N_TEST_PARTICLES  100
#define N_RANDOM_MODES    500
#define T_RUN_FOR_YEARS   1e5
#define B_REGULAR_COMP    0.0                                         //microGauss
#define B_TURBULENT_COMP  4.0                                         //microGauss
#define E_TOTAL           1e18                                        //eV
#define LAMBDA_MAX        10.0                                        //pc
#define LAMBDA_MIN        0.2                                         //pc
#define Q_CHARGE          1                                           //# electron charges (integer)
#define M_MASS            938.2720813                                 //MeV/c^2
#define ERROR_MAX         1.0e-03
#define ERROR_MIN         1.0e-08

int initialize_init(Trajectory_initializer &init, int procID);

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  int numProcs;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  int procID;
  MPI_Comm_rank(MPI_COMM_WORLD, &procID);
  //procID++;

/******************************************/ 
/************* MAIN PROCESSES *************/
  if(procID == 0){                                                      //Process to finalize calculations
    clock_t begin = clock();

    MPI_Status stat;

    int log_years = log10(T_RUN_FOR_YEARS) + 1;
    int length_Da_ij = log_years * 6; 


    //Vector of vectors to hold D_ij for each instane of the magnetic field.
    //i.e. D_ij[a] = D^(a)_ij
    vector< vector<double> > D_ij(numProcs - 1, vector<double>(length_Da_ij, 0));

    //Vector to hold eigenvalues
    //each outer vector holds the eigenvector for one logarithmic year 
    //(i.e eigenvalue_summed[instance][year] holds eigenvalue for instance i in year 10^i)
    vector< vector< vector<double> > > eigenvalues_summed(numProcs - 1, vector< vector<double> >(log_years, vector<double>(3, 0)));
    vector< vector<double> > eigenvalues_by_year(log_years, vector<double>(3, 0));
    vector<double> eigenvalues_current(3, 0);

    //vector to hold the recieved D^(a)_ij
    vector<double> Da_ij(length_Da_ij, 0);


    /* Recieve D_ij from all processes */
    for(int i = 1; i < numProcs; i++){                                  //Start from process 1, as this is process 0

      
      MPI_Recv(Da_ij.data(), length_Da_ij, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &stat);
      
      D_ij[stat.MPI_SOURCE - 1] = Da_ij;
      //D_ij[i] will now hold the matrix for D^(a)_ij

      std::cout << "I recieved D_ij from process " << stat.MPI_SOURCE << std::endl;
    }//End for i



    /* Calculate Eigenvalues */
    for(int instance = 0; instance < D_ij.size(); instance++){
      for(int i = 0; i < D_ij[instance].size() - 1; i += 6){
            
        calculate_eigenvalues_3x3_sym(D_ij[instance], i, eigenvalues_current);

        for(int d_i = 0; d_i < 3; d_i++){
          eigenvalues_summed[instance][i/6][d_i] = eigenvalues_current[d_i];    
          //eigenvalues_summed[instance][i/6] has the eigenvalues of D_ij[instance] for each logarithmic year (10^i/6)
        }

      }//End for i               
    }//End for instance



    /* sum up the eigenvalues for all instances */
    for(int instance = 0; instance < eigenvalues_summed.size(); instance++){
      for(int year = 0; year < eigenvalues_summed[instance].size(); year++){
        for(int d_i = 0; d_i < 3; d_i++){                                 //d_i is the eigenvalue d_1, d_2 or d_3 of the matrix D_ij[instance]
          eigenvalues_by_year[year][d_i] += eigenvalues_summed[instance][year][d_i];
        }
      }
    }//Eigenvalues_by_year now holds the sum in the calculation of yearly average eigenvalue


    for(int year = 0; year < eigenvalues_by_year.size(); year++){         //Divide sum by number of instances
      for(int d_i = 0; d_i < 3; d_i++){
        eigenvalues_by_year[year][d_i] = eigenvalues_by_year[year][d_i] / (numProcs - 1);
      }
    }//eigenvalue_by_year now holds the average eigenvalues for each logarithmic year.

  

  /* Write the eigenvalues to a file */
  std::ofstream file;
  file.open("data/eigenvalues.dat");

  if(!file){
    std::cout << "Couldn't open data/eigenvalues.dat. \n";
  }

  for(int year = 0; year < eigenvalues_by_year.size(); year++){
    file << eigenvalues_by_year[year][0] << ' ' << eigenvalues_by_year[year][1] << ' ' << eigenvalues_by_year[year][2] << '\n';
  }
  std::cout << "Eigenvalues has been written to data/eigenvalues.dat. \n";
  file.close();

  //This should now have calculated the eigenvalues.

  clock_t end = clock();
  double elapsed_secs = double(end - begin)/CLOCKS_PER_SEC;

  std::cout << "Rank 0 ended in " << elapsed_secs << " seconds \n";

/*************************************************/    
/************* CALCULATION PROCESSES *************/
  }else{            

    //Settings are now done globally at the top of the program.

    double ran1, ran2, ran3, ranTot;

    Trajectory_initializer init;

    initialize_init(init, procID);

    Ran rng(15321 + 100*procID + init.seed);                            //To generate directions
    Trajectory trajectory(init);
    
    for(int i = 0; i <= init.N; i++){                                   //i <= number of particles to test

      //Pick three random numbers to generate direction of particle
      ran1 = rng.doub(); ran2 = rng.doub(); ran3 = rng.doub();
      ranTot = ran1 + ran2 + ran3;

      trajectory.vx[0] = (ran1/ranTot) * trajectory.v_tot_init;
      trajectory.vy[0] = (ran2/ranTot) * trajectory.v_tot_init;
      trajectory.vz[0] = (ran3/ranTot) * trajectory.v_tot_init;

      trajectory.x[0] = 0.0; trajectory.y[0] = 0.0; trajectory.z[0] = 0.0;
   
      trajectory.RK4t_step_size_control_nw();                           //One step of the sum in D_ij is done per loop

      if( (static_cast<int>(((static_cast<double>(i)/init.N) * 1000))/10) % 10 == 0){
        std::cout << "Progress rank " << procID << ": " << (static_cast<double>(i)/init.N) * 100 << "\% \n";
      }


      //After this loop, D_ij holds the sum in the equation for D_ij
    }//End for i
  
    for(int i = 0; i < trajectory.D_ij.size(); i += 6){
      for(int j = 0; j < 6; j++){

        //i will represent the time in logarithmic years (i.e t_y = 10^(i/6))
        trajectory.D_ij[i + j] = trajectory.D_ij[i + j] / (2.0 * init.N * (static_cast<double>(std::pow(10,i/6))));
        //D_ij has now been calculated for one instance of the magnetic field.

      }
    } 
    
    MPI_Send(trajectory.D_ij.data(), trajectory.D_ij.size(), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    
  }
/************* END CALCULATION PROCESSES *************/

  MPI_Finalize();
  return 0;
}

/************* END GENERATE SAMPLES *************/
/************************************************/


int initialize_init(Trajectory_initializer &init, int procID){
  init.procID = procID;

  std::srand(time(0));  //rand() +                                    //Add this line to seed for "more" randomness
  init.seed = 1000 * procID;                                                 //Sets the seed for the RNG. Set to const 
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

