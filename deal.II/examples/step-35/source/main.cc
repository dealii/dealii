#include "../include/NavierStokes.h"


#include <iostream>
#include <sys/time.h>



int main( int argc, char **argv ){
  try{
    Data_Storage data;
    if( argc<2 ){
      std::cout<<std::endl<<std::endl
               <<" Missing parameter file!!"<<std::endl
               <<" Aborting."<<std::endl<<std::endl;
      data.print_usage();
      exit(1);
    }
    data.read_data( argv[1] );
    deallog.depth_console( data.verbose?2:0 );
    Navier_Stokes_Projection<deal_II_dimension> test( data );
    test.Create_Triangulation( data.n_of_global_refines );
    timeval init_time, end_time;
    gettimeofday( &init_time, 0 );
    for( double dt = data.initial_dt; dt >= data.final_dt; dt /= data.dt_decrement ){
      std::cout<<" dt = "<<dt<<std::endl;
      test.set_dt( dt );
      test.Initialize();
      test.run( data.verbose, data.output );
      test.Post_Process();
      std::cout<<"====================================="<<std::endl<<std::endl;
    }
    gettimeofday( &end_time, 0 );
    std::cout<<" Global time = "
             <<( double( end_time.tv_sec - init_time.tv_sec )+1e-6*double( end_time.tv_usec - init_time.tv_usec ) )
             <<" sec"<<std::endl;
  }
  catch (std::exception &exc){
    std::cerr << std::endl << std::endl
            << "----------------------------------------------------"
            << std::endl;
    std::cerr << "Exception on processing: " << std::endl
            << exc.what() << std::endl
            << "Aborting!" << std::endl
            << "----------------------------------------------------"
            << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
    std::cerr << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
    return 1;
  }
  std::cout<<"----------------------------------------------------"
           <<std::endl
           <<"Apparently everything went fine!"
           <<std::endl
           <<"Don't forget to brush your teeth :-)"
           <<std::endl<<std::endl;
  return 0;
}
