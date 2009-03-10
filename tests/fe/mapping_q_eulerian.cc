//----------------------------  mapping_q_eulerian.cc  ---------------------------
//    mapping_q_eulerian.cc,v 1.3 2003/06/09 16:00:38 wolf Exp
//    Version: 
//
//    Copyright (C) 2008 by the deal.II authors and Joshua White
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mapping_q_eulerian.cc  ---------------------------

// compute some convergence results from computing pi on a mesh that
// is deformed to represent a quarter of a ring

#include "../tests.h"
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/numbers.h>
#include <base/convergence_table.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/grid_refinement.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_out.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_values.h>
#include <fe/mapping.h>
#include <fe/mapping_q1.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <fe/fe_system.h>
#include <fe/fe_q.h>
#include <fstream>
#include <iostream>

#include <fe/mapping_q_eulerian.h>


// .... IMPOSED DISPLACEMENT

template <int dim>
class ImposedDisplacement : public Function<dim>
{
  public:
    ImposedDisplacement() : Function<dim> (dim) { }
    virtual void vector_value(const Point<dim> &p,
                              Vector<double> &value) const;
};

template <>
void ImposedDisplacement<2>::vector_value(const Point<2> &p,
                                          Vector<double> &value) const
{
  double radius = 1 + (sqrt(5)-1)*p(0);
  double angle  = 0.5*numbers::PI*(1-p(1));
  value(0) = radius*sin(angle)-p(0);
  value(1) = radius*cos(angle)-p(1);
}


// .... MAPPING TEST CLASS

template <int dim>
class MappingTest
{
  public:
    MappingTest (unsigned int degree);
    ~MappingTest ();

    void run_test();
    void graphical_output();
    
  private:
    double compute_area();
    void explicitly_move_mesh();
    void write_tria_to_eps(std::string id);

    Triangulation<dim>     triangulation;
    DoFHandler<dim>        dof_handler;
    FESystem<dim>          fe;

    unsigned int           degree;

    ImposedDisplacement<dim> imposed_displacement;
    Vector<double>           displacements;
};


// .... CONSTRUCTOR

template <int dim>
MappingTest<dim>::MappingTest (unsigned int degree)
		:
		dof_handler (triangulation),
		fe (FE_Q<dim>(degree),dim),
                degree(degree)
{ }


// .... DESTRUCTOR

template <int dim>
MappingTest<dim>::~MappingTest () 
{
  dof_handler.clear ();
}


// .... COMPUTE AREA

template <int dim>
double MappingTest<dim>::compute_area () 
{  
  QGauss<dim>  quadrature_formula(degree+1);

  MappingQEulerian<dim> mapping(degree,displacements,dof_handler);

  FEValues<dim> fe_values (mapping, fe, quadrature_formula, 
			   update_JxW_values);

  const unsigned int   n_q_points = quadrature_formula.n_quadrature_points;

  long double area = 0.;

  typename DoFHandler<dim>::active_cell_iterator 
                              cell = dof_handler.begin_active(),
                              endc = dof_handler.end();

  for (; cell!=endc; ++cell) {
    fe_values.reinit (cell);
    for(unsigned int q=0; q<n_q_points; ++q) area += fe_values.JxW(q);
  }      

  return area;
}


// .... RUN TEST

template <int dim>
void MappingTest<dim>::run_test () 
{
  GridGenerator::hyper_cube (triangulation,0, 1);

  ConvergenceTable table;

  for(unsigned int ref_level = 0; 
                   ref_level < 5; 
                 ++ref_level, triangulation.refine_global(1)){

    dof_handler.distribute_dofs (fe);
    displacements.reinit (dof_handler.n_dofs());

    VectorTools::interpolate(MappingQ1<dim>(),dof_handler,
                             imposed_displacement,displacements);


    table.add_value("cells",triangulation.n_active_cells());
    table.add_value("dofs",dof_handler.n_dofs());

    long double area  = compute_area();
    long double error = std::fabs(numbers::PI-area)/numbers::PI;

    table.add_value("area",  static_cast<double> (area));
    table.add_value("error", static_cast<double> (error));
  }    

  table.set_precision("area", 8);
  table.set_precision("error", 4);
  table.set_scientific("error", true);
  table.evaluate_convergence_rates("error",
    ConvergenceTable::reduction_rate_log2); 
  table.write_text(deallog.get_file_stream());
  deallog << std::endl;

}


// .... EXPLICITLY MOVE MESH

template <int dim>
void MappingTest<dim>::explicitly_move_mesh () 
{  
  std::vector<bool> moved (triangulation.n_vertices(),false);
  unsigned int vpc = GeometryInfo<dim>::vertices_per_cell;

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active (),
    endc = dof_handler.end();

  for (; cell != endc; cell++) {
  for (unsigned int v=0; v < vpc; v++) {
  if (moved[cell->vertex_index(v)] == false){
      moved[cell->vertex_index(v)] =  true;
      Point<dim> vertex_disp;
      for (unsigned int d=0; d<dim; d++) {
        vertex_disp[d] = displacements(cell->vertex_dof_index(v,d));
      }
      cell->vertex(v) += vertex_disp;
  }}}
}



// .... GRAPHICAL OUTPUT

template <int dim>
void MappingTest<dim>::graphical_output () 
{
  GridGenerator::hyper_cube (triangulation,0, 1);
  triangulation.refine_global(4);

  dof_handler.distribute_dofs (fe);
  displacements.reinit (dof_handler.n_dofs());

  VectorTools::interpolate(MappingQ1<dim>(),dof_handler,
                           imposed_displacement,displacements);

  explicitly_move_mesh();
}


// .... MAIN

int main () 
{
  std::ofstream logfile ("mapping_q_eulerian/output");
  deallog << std::setprecision(2);
  deallog << std::fixed;  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

				   // convergence studies

  for(unsigned int degree = 1; degree <=4; ++degree)
  {
    deallog << ".... Q" << degree << " Mapping ...." << std::endl;
    MappingTest<2> test_one(degree);
    test_one.run_test();
  }

			// graphical output

  MappingTest<2> test_two(1);
  test_two.graphical_output();

  return 0;
}

