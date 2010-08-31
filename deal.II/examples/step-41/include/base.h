#ifndef BASE_H
#define BASE_H

#include <base/function.h>
#include <base/logstream.h>
#include <base/utilities.h>
#include <base/parsed_function.h>
#include <base/timer.h>
#include <grid/grid_tools.h>
#include <grid/grid_generator.h>
#include <numerics/data_out.h>

#include <grid/filtered_iterator.h>

#include "vector_space.h"
#include "domain.h"
#include "error_handler.h"
#include "linear_elastic.h"
#include "camclay.h"



using namespace dealii;
using namespace dealii::Functions;

template <int dim>
class Base 
{
 public:
  Base ();
  ~Base ();

  //the routine called from int main()
  void run ();
  
 private:

  //function with actual work, called from run()
  void run_cc (const unsigned int &cc);

  //function with actual work, called from run()
  void run_step (const unsigned int &cc,
		 const unsigned int &step);
  
  //parse the parameters and get things rollin
  void parameters ();

  //elastic predictor
  void elastic_predictor (double &time);

  //plastic corrector
  void plastic_corrector ();
  
  //writes the files
  void write_files (const unsigned int &cc,
		    const unsigned int &step,
		    const Vector<double> &my_elas,
		    const Vector<double> &my_plas,
		    const Vector<double> &my_other);

  void write_stress_strain (const unsigned int &cc,
			    const unsigned int &step,
			    const Vector<double> &my_elas);

  void write_plot_values(const unsigned int &cc,
			 const unsigned int &step,
			 const Vector<double> &my_elas,
			 const Vector<double> &my_hardening,
			 const Vector<double> &my_plastic,
			 const ParsedSymmetricTensorFunction<4,dim> &C);
  
  //A parameter handler object
  ParameterHandler prm;

  //ErrorHandler Object
  ErrorHandler<dim, Vector<double> >  error_handler;
  
  //This holds the triangulation information
  Domain<dim> domain;

  //A vectorspace for the thermoelasticity
  VectorSpace<dim> vspace;

  //Elasticity Class - to change models, just change the class
  LinearElastic<dim> elastic;
  //HypoElastic<dim> elastic;
  //HyperElastic<dim> elastic;
  //LinearThermoElastic<dim> elastic;

  //Plasticity Class - to change models, just change the class name
  CamClay<dim> plastic;
  //CamClayExplicit<dim> plastic;
  //CamClayInvariant<dim> plastic;
  //CamClayShearBand<dim> plastic;
  //J2Flow<dim> plastic;
  //J2FlowThermal<dim> plastic;
  
  //General Parameters
  double lin_red_tol;
  unsigned int num_steps;
  double end_time;
  unsigned int num_cc;
  unsigned int num_threads;
  unsigned int console_depth;

  //various functions for the laplace equation
  ParsedFunction<dim> exact_solution;


};

  
#endif
