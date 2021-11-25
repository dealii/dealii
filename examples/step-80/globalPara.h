
// # physics constants
double global_PI = 3.1415927;

// # Laser
double global_Pow_laser = 0.4;      // power [W]
double global_spotsize_at_e_2 = 20e-6;  // laser spot size at e^(-2) [m]
double global_c_laser = global_spotsize_at_e_2 / 4.0;   // C parameter in Gaussian func, [m]
double global_c_hwhm = global_c_laser * 2.35482 / 2;    // HWHM, [m]
double global_V_scan_x = 10e-3;   // scan speed, [m/s]

double global_init_position_x0 = -50e-6;    // initial spot center position

// # material
// thin film
double global_rho_Tio2 = 4200;
double global_C_Tio2 = 690;
double global_k_Tio2 = 4.8;
// substrate
double global_rho_glass = 2200;
double global_C_glass = 700;
double global_k_glass = 1.8;

double global_film_thickness = 400e-9;  // film thickness, [m]

// # simulation time
double global_simulation_time_step = 1e-5;
double global_simulation_end_time = 100e-6 / global_V_scan_x; //

// # about the MESH
#define BOUNDARY_NUM 11
