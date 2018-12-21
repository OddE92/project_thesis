#include "cpp_calculate_trajectory/class_trajectory.h"
#include "cpp_bfield/class_bfield.h"
#include "cpp_bfield/ran.h"
#include "calculate_eigenvalues.h"

#include <boost/mpi.hpp>

#include <iostream>
#include <fstream>
#include <string>

constexpr int N_TEST_PARTICLES      =   100;
constexpr int N_RANDOM_MODES        =   500;
constexpr int T_RUN_FOR_YEARS       =   1e5;
constexpr double B_REGULAR_COMP     =   0.0;                                         //microGauss
constexpr double B_TURBULENT_COMP   =   4.0;                                         //microGauss
constexpr double E_TOTAL            =   1e16;                                        //eV
constexpr double LAMBDA_MAX         =   10.0;                                        //pc
constexpr double LAMBDA_MIN         =   0.0000027;                                        //pc    (0.27 = Rl/10 for B = 4, E = e16)
constexpr double Q_CHARGE           =   1;                                           //# electron charges
constexpr double M_MASS             =   938.2720813;                                 //MeV/c^2
constexpr double ERROR_MAX          =   1.0e05;
constexpr double ERROR_MIN          =   1.0e01;
const int D_IJ_LENGTH               =   log10(T_RUN_FOR_YEARS) * 9 + 1;
constexpr bool GEN_TURB_GLOB        =   true;


int initialize_init(Trajectory_initializer &init, int procID);

int main(int argc, char* argv[])
{
  MPI_Init(&argc, &argv);

  int numProcs;
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  int procID;
  MPI_Comm_rank(MPI_COMM_WORLD, &procID);
  //procID = 1;       //For testing Comment for testing git

/******************************************/ 
/************* MAIN PROCESSES *************/
  if(procID == 0){                                                      //Process to finalize calculations
    clock_t begin = clock();

    MPI_Status stat;


    //This is the length needed to hold: 0.1, 0.2, 0.3, ... , 1, 2, 3, ... , 10, 20 , 30 and so on
    int log_years = D_IJ_LENGTH;               
    int length_Da_ij = log_years * 7;                                   //6 entries to hold a symmetric matrix


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
    vector<double> D_average(log_years, 0);


    /* Recieve D_ij from all processes */
    for(int i = 1; i < numProcs; i++){                                  //Start from process 1, as this is process 0

      
      MPI_Recv(Da_ij.data(), length_Da_ij, MPI_DOUBLE, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &stat);
      
      D_ij[stat.MPI_SOURCE - 1] = Da_ij;
      //D_ij[i] will now hold the matrix for D^(a)_ij

      std::cout << "I recieved D_ij from process " << stat.MPI_SOURCE << std::endl;
    }//End for i

    /* Calculate Eigenvalues */
    for(int instance = 0; instance < D_ij.size(); instance++){
      for(int i = 0; i < D_ij[instance].size() - 1; i += 7){

        calculate_eigenvalues_3x3_sym(D_ij[instance], i, eigenvalues_current);

        for(int d_i = 0; d_i < 3; d_i++){
          eigenvalues_summed[instance][i/7][d_i] = eigenvalues_current[d_i];    
          //eigenvalues_summed[instance][i/7] has the eigenvalues of D_ij[instance] for each logarithmic year (10^i/7)
        }

        //Creates the sum for the average diffusion coefficient, over all instances
        D_average[i/7] += D_ij[instance][i+6];


      }//End for i               
    }//End for instance



    /* sum up the eigenvalues for all instances */
    for(int instance = 0; instance < eigenvalues_summed.size(); instance++){
      for(int year = 0; year < eigenvalues_summed[instance].size(); year++){
        
        for(int d_i = 0; d_i < 3; d_i++){                  //d_i is the eigenvalue d_1, d_2 or d_3 of the matrix D_ij[instance]
          eigenvalues_by_year[year][d_i] += eigenvalues_summed[instance][year][d_i];
        }
      
      }
    }//Eigenvalues_by_year now holds the sum in the calculation of yearly average eigenvalue


    for(int year = 0; year < eigenvalues_by_year.size(); year++){         //Divide sum by number of instances
      for(int d_i = 0; d_i < 3; d_i++){
        eigenvalues_by_year[year][d_i] = eigenvalues_by_year[year][d_i] / (numProcs - 1);
      }

      D_average[year] = D_average[year] / (numProcs - 1);

    }//eigenvalue_by_year now holds the average eigenvalues for each logarithmic year.

  

  /* Write the eigenvalues to a file */
  std::ofstream file;
  file.open("data/eigenvalues.dat");

  if(!file){
    std::cout << "Couldn't open data/eigenvalues.dat. \n";
  }

  for(int year = 0; year < eigenvalues_by_year.size(); year++){
    file << eigenvalues_by_year[year][0] << ' ' << eigenvalues_by_year[year][1] << ' ' << eigenvalues_by_year[year][2];
    file << ' ' << D_average[year] << '\n';  
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
    double ran4, ran5, ran6;
    double percentCounter = 0.099999;

    Trajectory_initializer init;

    initialize_init(init, procID);

    Ran rng(15321 + 100*procID + init.seed);                            //To generate directions
    Trajectory trajectory(init);
    
    std::ofstream file;                                                 //To write x,y-coordinates of particles
    std::string filename;
    filename = "data/r_rank" + std::to_string(procID) + ".dat";
    file.open(filename);
    if(!file){
      std::cout << "Couldn't open " << filename << std::endl;
    }
      

    for(int i = 0; i < init.N; i++){                                    //i < number of particles to test
  
      //Pick three random numbers to generate direction of particle
      ran1 = rng.doub(); ran2 = rng.doub(); ran3 = rng.doub();
      ranTot = ran1 + ran2 + ran3;

    
      //Pick three random numbers to set the signs
      ran4 = rng.doub(); ran5 = rng.doub(); ran6 = rng.doub();
      if(ran4 > 0.5){ ran4 = -1.0; }else{ ran4 = 1.0; }
      if(ran5 > 0.5){ ran5 = -1.0; }else{ ran5 = 1.0; }
      if(ran6 > 0.5){ ran6 = -1.0; }else{ ran6 = 1.0; }

      trajectory.vx[0] = ran4 * (ran1/ranTot) * trajectory.v_tot_init;
      trajectory.vy[0] = ran5 * (ran2/ranTot) * trajectory.v_tot_init;
      trajectory.vz[0] = ran6 * (ran3/ranTot) * trajectory.v_tot_init;


      //trajectory.vx[0] = 1 * trajectory.v_tot_init; trajectory.vy[0] = 0; trajectory.vz[0] = 0;
      trajectory.x[0] = 0.0; trajectory.y[0] = 0.0; trajectory.z[0] = 0.0;


      trajectory.RK4t_step_size_control_nw();                           //One step of the sum in D_ij is done per loop
      
      //Print progress in % (static_cast<int>(((static_cast<double>(i)/init.N) * 1000))/10) % 10 == 0
      if( static_cast<double>(i+1)/init.N - percentCounter > __DBL_EPSILON__ ){
        std::cout << "Progress rank " << procID << ": " << (static_cast<double>(i+1)/init.N) * 100 << "\% \n";
        percentCounter += 0.1;
      }


      //Write position of particles to file
      for(int j = 0; j < trajectory.r_vect.size(); j += 3){
        file << trajectory.r_vect[j + 0] << ' ' << trajectory.r_vect[j + 1] << ' ' << trajectory.r_vect[j+2] << '\n';
      }

     //After this loop, D_ij holds the sum in the equation for D_ij
    }//End for particle i


    for(int i = 0; i < trajectory.D_ij.size(); i += 7){
      for(int j = 0; j < 7; j++){

        trajectory.D_ij[i + j] = trajectory.D_ij[i + j] / (2.0 * init.N * trajectory.D_ij_time[i/7]);      //std::pow(10, i/6)
        //D_ij has now been calculated for one instance of the magnetic field.

        if(j == 6) trajectory.D_ij[i + j] = trajectory.D_ij[i + j] / 3;     //To calculate the average diffusion coefficient

      }
    }//End for i

    
    
    std::cout << "Number of coordinate sets in r_vect: " << init.N * trajectory.r_vect.size() / 3.0 << std::endl;
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
  init.seed = rand() + 1000 * procID;                                 //Sets the seed for the RNG. Set to const 
                                                                      //to generate equal results

  init.n_k = N_RANDOM_MODES;                                          //#modes used to generate TMF
  init.N = N_TEST_PARTICLES;                                          //Number of test-particles
    
  init.t_end_y = T_RUN_FOR_YEARS;                                     //Set time in years
  init.t_start = 0.0;                                                 //Set time start and end here, and the timestep
  init.t_end = 31557600.0 * init.t_end_y; 
  init.dt = 0.1;        

  init.B_0 = B_REGULAR_COMP; init.B_rms_turb = B_TURBULENT_COMP;      //B_0 is the regular field, B_rms_turb is the RMS of the turb field

  init.E = E_TOTAL;                                                   //Energy in eV

  init.lambda_max = LAMBDA_MAX; init.lambda_min = LAMBDA_MIN;         //Wavelength in pc

  init.q = Q_CHARGE; init.m = M_MASS;                                 //Set the charge and mass
                                                                        
  init.max_err = ERROR_MAX;                                           //Set min and max error
  init.min_err = ERROR_MIN;
                              
  init.generate_turbulence = GEN_TURB_GLOB;                           //This decides if you generate a turbulence or not
    
  init.D_ij_size = D_IJ_LENGTH;                                       //Sets length of diffusion tensor-vector  

  init.x0 = 0.0; init.y0 = 0.0; init.z0 = 0.0;                        //Initial conditions
  init.vx0 = 1.0; init.vy0 = 0.0; init.vz0 = 0.0;
}

