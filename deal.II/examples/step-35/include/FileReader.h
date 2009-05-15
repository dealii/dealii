// This class is supposed to aid in the reading of parameter files
#ifndef _FILE_READER_H_
#define _FILE_READER_H_



#include <base/parameter_handler.h>



#include <fstream>
#include <iostream>



using namespace dealii;

enum Method_Formulation{
  METHOD_STANDARD,
  METHOD_ROTATIONAL
};

class Data_Storage{
  public:
    Data_Storage();
    ~Data_Storage();
    void read_data( char *filename );
    void print_usage();
    void print_status() const;
    // The data itself
    //// The type of method we want to use
    Method_Formulation form;
    //// physical data
    double initial_time,
           final_time,
           Reynolds;
    //// Time stepping data
    double initial_dt,
           final_dt,
           dt_decrement;
    //// Space discretization data
    unsigned int n_of_global_refines,
                 pressure_degree;
    //// Data to solve the velocity
    unsigned int vel_max_iterations,
                 vel_Krylov_size,
                 vel_off_diagonals,
                 vel_update_prec;
    double vel_eps,
           vel_diag_strength;
    //// Data to solve the projection
    unsigned int proj_max_iterations,
                 proj_off_diagonals;
    double proj_eps,
           proj_diag_strength;
    //// Data to do the pressure update step
    unsigned int pres_max_iterations,
                 pres_off_diagonals;
    double pres_eps,
           pres_diag_strength;
    //// Verbosity
    bool verbose;
    //// Frequency of the outputted data
    unsigned int output;

  protected:
    ParameterHandler prm;
};

#endif
