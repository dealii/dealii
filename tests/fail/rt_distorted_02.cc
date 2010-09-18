//----------------------------  rt_distorted_02.cc  ---------------------------
//    rt_distorted_02.cc,v 1.3 2003/06/09 16:00:38 wolf Exp
//    Version: 
//
//    Copyright (C) 2003, 2005, 2006, 2007, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  rt_distorted_02.cc  ---------------------------

/*
 * Snippet to demonstrate some properties of the deal.II implementation of
 * the RT spaces.
 */



#include "../tests.h"
#include <base/logstream.h>

#define PRECISION 2

#include <fstream>

std::ofstream logfile ("rt_distorted_02/output");

char buf[1000];

#include <fstream>

#include <grid/grid_generator.h>
#include <grid/grid_out.h>
#include <fe/mapping_q.h>
#include <fe/mapping_q1_eulerian.h>

#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <lac/vector.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <lac/vector_memory.h>
#include <lac/sparsity_pattern.h>
#include <lac/sparse_matrix.h>

#include <base/function.h>
#include <base/quadrature_lib.h>
#include <lac/constraint_matrix.h>
#include <dofs/dof_tools.h>

#include <fe/fe_system.h>
#include <fe/fe.h>
#include <fe/fe_system.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/fe_dgq.h>
#include <fe/mapping_q1_eulerian.h>

#include <fstream>


template <int dim>
class TestMap1 : public Function<dim>
{
  public:
    TestMap1 (const unsigned int n_components) :
		    Function<dim> (n_components)
      {}
  
    virtual ~TestMap1 () {}
  
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
  
    void vector_value (const Point<dim> &p,
		       Vector<double>   &return_value) const;
};



template <int dim>
double
TestMap1<dim>::value (const Point<dim>   &p,
		      const unsigned int  component) const 
{
				   // u = x^2, v = y^2, dudx = 2 x, dvdy = 2 y, div u = 2x + 2y
				   // I.e. \int div u = 2 (for unit square)

  if (component == 0)
    return (p(0) * p(0));
  if (component == 1)
    return (p(1) * p(1));

  return (0);
}


template <int dim>
void TestMap1<dim>::vector_value (const Point<dim> &p,
				  Vector<double>   &return_value) const
{
  Assert (return_value.size() == this->n_components,
	  ExcDimensionMismatch (return_value.size(), this->n_components));
  
				   // Just fill the vector with the appropriate components
  for (unsigned int iCount = 0; iCount < this->n_components; iCount++)
    return_value (iCount) = value (p, iCount);
}



///-----------------------------------------------------------------------

template <int dim>
class TestDef1 : public Function<dim>
{
  private:
    const double phi;

  public:
    TestDef1 (const unsigned int n_components, const double ph) :
		    Function<dim> (n_components),
		    phi (ph)
      {}
    
    virtual ~TestDef1 () {}
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
    
    void vector_value (const Point<dim> &p,
		       Vector<double>   &return_value) const;
};



template <int dim>
double
TestDef1<dim>::value (const Point<dim>   &p,
		      const unsigned int  component) const 
{
  Point<2> center;
  center(0) = 0.5;
  center(1) = 0.5;
  double rad = p.distance (center),
       phi_p = atan2 (p(0) - center(0), p(1) - center(1));

  if (component == 0)
    return rad * (sin (phi + phi_p) - sin (phi_p));
  else
    return rad * (cos (phi + phi_p) - cos (phi_p));
}


template <int dim>
void TestDef1<dim>::vector_value (const Point<dim> &p,
				  Vector<double>   &return_value) const
{
  Assert (return_value.size() == this->n_components,
	  ExcDimensionMismatch (return_value.size(), this->n_components));
  for (unsigned int iCount = 0; iCount < this->n_components; iCount++)
    return_value (iCount) = value (p, iCount);
}


///-----------------------------------------------------------------------


template <int dim>
class TestDef2 : public Function<dim>
{
  private:
    const double scale;

  public:
    TestDef2 (const unsigned int n_components, const double sc) :
		    Function<dim> (n_components),
		    scale (sc)
      {}
    
    virtual ~TestDef2 () {}
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
    
    void vector_value (const Point<dim> &p,
		       Vector<double>   &return_value) const;
};


template <int dim>
double
TestDef2<dim>::value (const Point<dim>   &p,
		      const unsigned int  component) const 
{
  double x = p(0),
	 y = p(1);

  if (component == 0)
    return scale * x;
  else
    return scale * y;
}


template <int dim>
void TestDef2<dim>::vector_value (const Point<dim> &p,
				  Vector<double>   &return_value) const
{
  Assert (return_value.size() == this->n_components,
	  ExcDimensionMismatch (return_value.size(), this->n_components));
  for (unsigned int iCount = 0; iCount < this->n_components; iCount++)
    return_value (iCount) = value (p, iCount);
}


///-----------------------------------------------------------------------
// testDef3 implements parallelograms ...


template <int dim>
class TestDef3 : public Function<dim>
{
  private:
    const double scale;

  public:
    TestDef3 (const unsigned int n_components, const double sc) :
		    Function<dim> (n_components),
		    scale (sc)
      {}
    
    virtual ~TestDef3 () {}
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
    
    void vector_value (const Point<dim> &p,
		       Vector<double>   &return_value) const;
};


template <int dim>
double
TestDef3<dim>::value (const Point<dim>   &p,
		      const unsigned int  component) const 
{
  double y = p(1);

  if (component == 0)
    return scale * y;
  else
    return 0;
}


template <int dim>
void TestDef3<dim>::vector_value (const Point<dim> &p,
				  Vector<double>   &return_value) const
{
  Assert (return_value.size() == this->n_components,
	  ExcDimensionMismatch (return_value.size(), this->n_components));
  for (unsigned int iCount = 0; iCount < this->n_components; iCount++)
    return_value (iCount) = value (p, iCount);
}



/*
 * Check the value of the derivative field.
 */

double EvaluateDiver (Mapping<2> &mapping,
		      DoFHandler<2> &dof_handler,
		      Vector<double> &solution)
{
				   // Use a high order quadrature. 
  QGauss<2> quad (6);
  FEValues<2> fe_values (mapping, dof_handler.get_fe (), quad, 
			 UpdateFlags(update_values    |
				     update_q_points  |
				     update_gradients |
				     update_JxW_values|
				     update_contravariant_transformation));

  const unsigned int   n_q_points    = quad.n_quadrature_points;
  const unsigned int   n_components   = dof_handler.get_fe().n_components();

				   // Cell iterators
  DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active(),
				      endc = dof_handler.end();
  double result = 0;

  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);

				       // Get values from solution vector (For Trap.Rule)
      std::vector<Vector<double> > this_value
	(n_q_points, Vector<double>(n_components));
      fe_values.get_function_values (solution, this_value);

      std::vector<std::vector<Tensor<1,2> > >
	grads_here (n_q_points, 
		    std::vector<Tensor<1,2> > (n_components));
      fe_values.get_function_grads (solution, grads_here);

      for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
	{
	  double JxW = fe_values.JxW (q_point);
	  double dudx = grads_here[q_point][0][0];
	  double dvdy = grads_here[q_point][1][1];
	  result += (dudx + dvdy) * JxW;
	}
    }

  return (result);
}


int main ()
{
  logfile.precision (PRECISION);
  logfile.setf(std::ios::fixed);  
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Triangulation<2> tria_test;
  Point<2> p1 (0,0),
    p2 (1, 1);

  GridGenerator::hyper_rectangle (tria_test, p1, p2);
  tria_test.refine_global (1);
  tria_test.distort_random (0.4);

				   // Create a DoFHandler for the RT space
  FE_RaviartThomas<2> fe (2);
  DoFHandler<2> dof_handler (tria_test);
  dof_handler.distribute_dofs (fe);

  QGauss<2> quad_temp(6);
  sprintf (buf,"DoFs per Quad: %i per Line %i per Vert %i\n", fe.dofs_per_quad, fe.dofs_per_line, fe.dofs_per_vertex);
  deallog << buf;
  
  sprintf (buf,"QPoints %i\n", quad_temp.n_quadrature_points);
  deallog << buf;

				   // Create an deformation object for the Eulerian mapping
  FESystem<2> fe_def (FE_Q<2>(1), 2);
  DoFHandler<2> dof_handler_def (tria_test);
  dof_handler_def.distribute_dofs (fe_def);

				   // Alloc some DoFs
  Vector<double> solution,
    solution_q,
    deformation;
  solution.reinit (dof_handler.n_dofs ());
  solution_q.reinit (dof_handler_def.n_dofs ());
  deformation.reinit (dof_handler_def.n_dofs ());

				   // Project solution onto RT FE field
  ConstraintMatrix     hn_constraints;
  hn_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler, 
					   hn_constraints);
  hn_constraints.close ();
  VectorTools::project (dof_handler, hn_constraints,
			QGauss<2> (6), TestMap1<2>(2),
			solution);

				   // Project reference solution onto RT FE field
  ConstraintMatrix     hn_constraints_def;
  hn_constraints_def.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler_def, 
					   hn_constraints_def);
  hn_constraints_def.close ();

  VectorTools::project (dof_handler_def, hn_constraints_def,
			QGauss<2> (6), TestMap1<2>(2),
			solution_q);
  
  MappingQ1Eulerian<2> mapping_euler (deformation, dof_handler_def);

  sprintf (buf,"FE_RT Area %e  FE_Q Area %e\n",
	   EvaluateDiver (mapping_euler, dof_handler, solution),
	   EvaluateDiver (mapping_euler, dof_handler_def, solution_q));
  deallog << buf;

				   // Try rotating the elements
  for (double rotat = 0; rotat < 2 * M_PI; rotat += 0.25 * M_PI)
    {
				       // Rotate element
      VectorTools::project (dof_handler_def, hn_constraints_def,
			    QGauss<2> (6), TestDef1<2>(2, rotat),
			    deformation);
      sprintf (buf,"phi = %e FE_RT Area %e  FE_Q Area %e\n",
	       rotat,
	       EvaluateDiver (mapping_euler, dof_handler, solution),
	       EvaluateDiver (mapping_euler, dof_handler_def, solution_q));
      deallog << buf;
    }

				   // Try resizing the elements
  for (double scale = -0.75; scale < 4.0; scale += 0.25)
    {
      VectorTools::project (dof_handler_def, hn_constraints_def,
			    QGauss<2> (6), TestDef2<2>(2, scale),
			    deformation);
      sprintf (buf,"Scale = %e FE_RT Area %e  FE_Q Area %e\n",
	       scale,
	       EvaluateDiver (mapping_euler, dof_handler, solution),
	       EvaluateDiver (mapping_euler, dof_handler_def, solution_q));
      deallog << buf;
    }

				   // Try paralellogramming the elements
  for (double scale = -1.0; scale < 1.0; scale += 0.25)
    {
      VectorTools::project (dof_handler_def, hn_constraints_def,
			    QGauss<2> (6), TestDef3<2>(2, scale),
			    deformation);
      sprintf (buf,"Scale = %e FE_RT Area %e  FE_Q Area %e\n",
	       scale,
	       EvaluateDiver (mapping_euler, dof_handler, solution),
	       EvaluateDiver (mapping_euler, dof_handler_def, solution_q));
      deallog << buf;
    }


  return (0);
}
