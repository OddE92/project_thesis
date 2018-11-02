#ifndef TRAJECTORY_INITIALIZER
#define TRAJECTORY_INITIALIZER

/******* INITIALIZING LIST ******/
struct Trajectory_initializer{
    double t_start, t_end, dt;
    double B_0, B_rms_turb;
    double x0, y0, z0;
    double vx0, vy0, vz0;
    double E;
    double q, m;
    double min_err, max_err;
    double xyz_max, xyz_min;
    int numGridPoints;
    int seed;

    bool generate_discrete_field;
    bool generate_turbulence;
};
/***** END INITIALIZING LIST ****/


#endif