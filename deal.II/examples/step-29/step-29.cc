/*    $Id: step-28.cc 14713 2007-05-27 04:05:41Z bangerth $       */
/*    Version: $Name:  $                                          */
/*                                                                */
/*    Copyright (C) 2007 by the deal.II authors and Moritz Allmaras     */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>
#include <base/parameter_handler.h>
#include <base/subscriptor.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/sparse_direct.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_system.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <numerics/vectors.h>
#include <fstream>

#define DIM 2


using namespace dealii;

template <int dim>
class DirichletBoundaryValues : public Function<dim>
{
  public:
    DirichletBoundaryValues() : Function<dim> (2) {};
    
    virtual void vector_value (	const Point<dim> &p,
                                Vector<double>   &values) const;
 
    virtual void vector_value_list (const std::vector<Point<dim> > &	points,
                                    std::vector<Vector<double> > &		value_list) const;
};


template <int dim>
inline
void DirichletBoundaryValues<dim>::vector_value (	const Point<dim> &	/*p*/,
							Vector<double> &		values) const 
{
  Assert (values.size() == 2, ExcDimensionMismatch (values.size(), 2));
	
  values(0) = 1;
  values(1) = 0;
}


template <int dim>
void DirichletBoundaryValues<dim>::vector_value_list (const std::vector<Point<dim> > &	points,
						      std::vector<Vector<double> > &		value_list) const 
{
  Assert (value_list.size() == points.size(), 
          ExcDimensionMismatch (value_list.size(), points.size()));

  for (unsigned int p=0; p<points.size(); ++p)
    DirichletBoundaryValues<dim>::vector_value (points[p], value_list[p]);
}


class ParameterReader : public Subscriptor
{
  public:
    ParameterReader(ParameterHandler &);
    void read_parameters();
 
  private:
    void declare_parameters();
    ParameterHandler *prm;
};


ParameterReader::ParameterReader(ParameterHandler &paramhandler)
		:
		prm(&paramhandler)
{}


void ParameterReader::declare_parameters()
{
  prm->declare_entry(	"Number of refinements", "5",
                    	Patterns::Integer(1,10),
                    	"Number of global mesh refinement steps "
			"applied to initial coarse grid");

  prm->declare_entry(	"Focal distance", "0.3",
                    	Patterns::Double(0),
                    	"Distance of the focal point of the lens "
			"to the x-axis (or xy-plane in 3D)");

  prm->declare_entry(	"c", "1.5e5",
                    	Patterns::Double(),
                    	"Wave speed");

  prm->declare_entry(	"omega", "1.5e7",
                    	Patterns::Double(),
                    	"Frequency");

  prm->declare_entry(	"Output file", "solution",
			Patterns::Anything(),
			"Name of the output file (without extension)");

  DataOutInterface<1>::declare_parameters (*prm);
}


void ParameterReader::read_parameters()
{
  declare_parameters();

  const std::string parameter_file = "step-29.prm";
  prm->read_input (parameter_file);
}


template <int dim>
class UltrasoundProblem 
{
  public:
    UltrasoundProblem (ParameterHandler &);
    ~UltrasoundProblem ();
    void run ();
    
  private:
    void make_grid ();
    void setup_system ();
    void assemble_system ();
    void solve ();
    void postprocess ();
    void output_results () const;

    ParameterHandler &prm; 

    Triangulation<dim>			triangulation;
    DoFHandler<dim>         dof_handler;
    FESystem<dim>						fe;
    
    SparsityPattern         sparsity_pattern;
    SparseMatrix<double>    system_matrix;
		
    Vector<double>          solution, system_rhs, intensity;

    const double 						c, omega;
};


template <int dim>
UltrasoundProblem<dim>::UltrasoundProblem (ParameterHandler&	param) 
		:
		prm(param),
		dof_handler(triangulation),
		fe(FE_Q<dim>(1), 2),
		c(prm.get_double("c")), 
		omega(prm.get_double("omega"))
{}


template <int dim>
UltrasoundProblem<dim>::~UltrasoundProblem () 
{
  dof_handler.clear();
}


template <int dim>
void UltrasoundProblem<dim>::make_grid ()
{
  GridGenerator::subdivided_hyper_cube (triangulation, 5, 0, 1);
	
  const Point<dim> 	transducer = (dim == 2) ? 
				     Point<dim> (0.5, 0.0) :
				     Point<dim> (0.5, 0.5, 0.0), 
		       focal_point = (dim == 2) ?
				     Point<dim> (0.5, prm.get_double("Focal distance")) :
				     Point<dim> (0.5, 0.5, prm.get_double("Focal distance"));

  double radius = sqrt(	(focal_point.distance(transducer) * 
			 focal_point.distance(transducer)) + 
			((dim==2) ? 0.01 : 0.02));

  typename Triangulation<dim>::cell_iterator
    cell = triangulation.begin (),
    endc = triangulation.end();
	  
  for (; cell!=endc; ++cell)
    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
      if ( cell->face(face)->at_boundary() &&
	   ((cell->face(face)->center() - transducer).square() < 0.01) )
	cell->face(face)->set_boundary_indicator (1);

  const HyperBallBoundary<dim> boundary(focal_point, radius);
  triangulation.set_boundary(1, boundary);

  triangulation.refine_global (prm.get_integer("Number of refinements"));

  deallog << " Number of active cells:  "
	  << triangulation.n_active_cells()
	  << std::endl;

  triangulation.set_boundary(1);
} 


template <int dim>
void UltrasoundProblem<dim>::setup_system ()
{
  dof_handler.distribute_dofs (fe);

  deallog << " Number of degrees of freedom: "
          << dof_handler.n_dofs()
          << std::endl;

  sparsity_pattern.reinit (	dof_handler.n_dofs(),
				dof_handler.n_dofs(),
				dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);
  system_rhs.reinit (dof_handler.n_dofs());
  solution.reinit (dof_handler.n_dofs());
  intensity.reinit(triangulation.n_active_cells());
}


template <int dim>
void UltrasoundProblem<dim>::assemble_system () 
{ 
  const double om2 = omega * omega;
  const double c2 = c * c;
	
  QGauss<dim>		quadrature_formula(3);
  QGauss<dim-1> face_quadrature_formula(3);

  const unsigned int 	n_q_points    		= quadrature_formula.n_quadrature_points;
  const unsigned int 	n_face_q_points	= face_quadrature_formula.n_quadrature_points;

  const unsigned int 	dofs_per_cell = fe.dofs_per_cell;

  FullMatrix<double>	cell_matrix (dofs_per_cell, dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  
  FEValues<dim>  fe_values (fe, quadrature_formula, 
			    update_values | update_gradients |
			    update_JxW_values);

  FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, 
				    update_values | update_JxW_values);

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      cell_matrix = 0;
      fe_values.reinit (cell);

      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    {
	      if (fe.system_to_component_index(i).first == 
		  fe.system_to_component_index(j).first)
		{
		  for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
		    cell_matrix(i,j) += (((fe_values.shape_value(i,q_point) *
					   fe_values.shape_value(j,q_point)) *
					  (- om2)
					  +
					  (fe_values.shape_grad(i,q_point) *
					   fe_values.shape_grad(j,q_point)) *
					  c2) *
					 fe_values.JxW(q_point));
		}
	    }
	}

      if (cell->at_boundary())
	for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	  if (	cell->face(face)->at_boundary() &&
		(cell->face(face)->boundary_indicator() == 0) )
	    {
	      fe_face_values.reinit (cell, face);

	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		for (unsigned int j=0; j<dofs_per_cell; ++j)
		  if ((fe.system_to_component_index(i).first != 
		       fe.system_to_component_index(j).first) &&
		      fe.has_support_on_face(i, face) &&
		      fe.has_support_on_face(j, face))

		    for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
		      cell_matrix(i,j) += ((fe.system_to_component_index(i).first) ? 1 : (-1)) * 
					  fe_face_values.shape_value(i,q_point) *
					  fe_face_values.shape_value(j,q_point) *
					  c * omega *
					  fe_face_values.JxW(q_point);
	    }

      cell->get_dof_indices (local_dof_indices);

      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  system_matrix.add (	local_dof_indices[i],
				local_dof_indices[j],
				cell_matrix(i,j));
    }

  std::map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
					    1,
					    DirichletBoundaryValues<dim>(),
					    boundary_values);
  MatrixTools::apply_boundary_values (boundary_values,
				      system_matrix,
				      solution,
				      system_rhs);
}


template <int dim>
void UltrasoundProblem<dim>::solve () 
{
  SparseDirectUMFPACK			A_direct;

  A_direct.initialize(system_matrix);
  A_direct.vmult(solution,system_rhs);
}


template <int dim>
void UltrasoundProblem<dim>::postprocess () 
{
  QMidpoint<dim> 	midpoint_rule;
  FEValues<dim> fe_values(fe, midpoint_rule, update_values);

  std::vector<Vector<double> > 	values;
  values.resize(1);
  values[0].reinit(2);

  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (unsigned int	i=0; cell!=endc; ++cell, ++i)
    {
      fe_values.reinit(cell);
		
      fe_values.get_function_values (solution, values);

      intensity(i) = sqrt(	values[0](0)*values[0](0) 
				+ values[0](1)*values[0](1));
    }
}


template <int dim>
void UltrasoundProblem<dim>::output_results () const
{
  const std::string 	output_file(prm.get("Output file")), 
    output_format(prm.get("Output format"));
	
  DataOutBase::OutputFormat		format = DataOutBase::parse_output_format(output_format);

  const std::string		filename = 	output_file + 
						DataOutBase::default_suffix(format);

  std::ofstream output (filename.c_str());

  DataOut<dim> data_out;
  data_out.parse_parameters(prm);
  data_out.attach_dof_handler (dof_handler);

  std::vector<std::string> solution_names;
  solution_names.push_back ("Re_u");
  solution_names.push_back ("Im_u");
  data_out.add_data_vector (solution, solution_names);
	
  data_out.add_data_vector (intensity, "intensity");

  data_out.build_patches ();

  data_out.write (output, format);
}


template <int dim>
void UltrasoundProblem<dim>::run () 
{
  make_grid ();
  setup_system ();
  assemble_system ();
  solve ();
  postprocess ();
  output_results ();
}


int main () 
{
  try
    {
      ParameterHandler 	prm;
      ParameterReader		param(prm);
      param.read_parameters();

      Assert (DIM > 1, ExcNotImplemented());
		
      UltrasoundProblem<DIM>	ultrasound_problem (prm);
      ultrasound_problem.run ();
    }
  catch (std::exception &exc)
    {
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
  catch (...) 
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    }

  return 0;
}
