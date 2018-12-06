#ifndef TRAJECTORY_INITIALIZER
#define TRAJECTORY_INITIALIZER

/******* INITIALIZING LIST ******/
struct Trajectory_initializer{
    double t_start, t_end, dt;
    int t_end_y;
    double B_0, B_rms_turb;
    double lambda_max, lambda_min;
    double x0, y0, z0;
    double vx0, vy0, vz0;
    double E;
    double q, m;
    double min_err, max_err;
    int seed;
    int procID;
    int n_k;
    int N;
    int D_ij_size;
    
    bool generate_turbulence;
};
/***** END INITIALIZING LIST ****/


#endif