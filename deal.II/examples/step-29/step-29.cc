				// @sect3{Include files}

				// The following header files are unchanged 
				// from step-7 and have been discussed before:

#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>

#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>

#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>

#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>

#include <fe/fe_q.h>
#include <fe/fe_values.h>

#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <numerics/vectors.h>

#include <fstream>

				// This header file is needed for  
				// the ParameterHandler class that we will
				// use to read parameters from a 
				// configuration file during runtime:
#include <base/parameter_handler.h>

#include <lac/sparse_direct.h>
#include <fe/fe_system.h>

#define DIM 2

using namespace dealii;


template <int dim>
class DirichletBoundaryValues : public Function<dim>
{
  public:
    DirichletBoundaryValues() : Function<dim> (2) {};

    virtual void vector_value (const Point<dim> &p,
			       Vector<double>   &values) const;

    virtual void vector_value_list (const std::vector<Point<dim> > &points,
				    std::vector<Vector<double> >   &value_list) const;
};


template <int dim>
inline
void DirichletBoundaryValues<dim>::vector_value (const Point<dim> &/*p*/,
						 Vector<double>   &values) const 
{
  Assert (values.size() == 2, ExcDimensionMismatch (values.size(), 2));

  values(0) = 1;
  values(1) = 0;
}


template <int dim>
void DirichletBoundaryValues<dim>::vector_value_list (const std::vector<Point<dim> > &points,
						      std::vector<Vector<double> >   &value_list) const 
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
    ParameterHandler &prm;
};


ParameterReader::ParameterReader(ParameterHandler &paramhandler)
  :
  prm(paramhandler)
{}


void ParameterReader::read_parameters()
{
  declare_parameters();

  const std::string parameter_file = "step-29.prm";
  prm.read_input (parameter_file);
}


void ParameterReader::declare_parameters()
{
  prm.enter_subsection ("Mesh & geometry parameters");

    prm.declare_entry("Number of refinements", "6",
		      Patterns::Integer(1,10),
		      "Number of global mesh refinement steps "
		      "applied to initial coarse grid");

     prm.declare_entry("Focal distance", "0.3",
		       Patterns::Double(0),
		       "Distance of the focal point of the lens "
		       "to the x-axis (or xy-plane in 3D)");

  prm.leave_subsection ();

  prm.enter_subsection ("Physical constants");

    prm.declare_entry("c", "1.5e5",
		      Patterns::Double(),
		      "Wave speed");

    prm.declare_entry("omega", "5.0e7",
		      Patterns::Double(),
		      "Frequency");

  prm.leave_subsection ();

  prm.enter_subsection ("Output parameters");

    prm.declare_entry("Output file", "solution",
		      Patterns::Anything(),
		      "Name of the output file (without extension)");

    DataOutInterface<1>::declare_parameters (prm);

  prm.leave_subsection ();
}


template <int dim>
class Postprocessor : public DataPostprocessor<dim>
{
  public:

    void compute_derived_quantities_vector (
			const std::vector< Vector< double > > &, 
			const std::vector< std::vector< Tensor< 1, dim > > > &, 
			const std::vector< std::vector< Tensor< 2, dim > > > &, 
			const std::vector< Point< dim > >                    &,
			std::vector< Vector< double > > &
			) const;

    std::vector<std::string> get_names () const;
    UpdateFlags              get_needed_update_flags () const;
    unsigned int             n_output_variables () const;
};


template <int dim>
std::vector<std::string>
Postprocessor<dim>::get_names() const
{
  std::vector<std::string> field_names;

  field_names.push_back("Re_u");
  field_names.push_back("Im_u");
  field_names.push_back("Intensity");

  return field_names;
}


template <int dim>
UpdateFlags
Postprocessor<dim>::get_needed_update_flags () const
{
  return update_values;
}


template <int dim>
unsigned int
Postprocessor<dim>::n_output_variables () const
{
  return 3;
}


template <int dim>
void
Postprocessor<dim>::compute_derived_quantities_vector (
			const std::vector< Vector< double > >                 &uh,
			const std::vector< std::vector< Tensor< 1, dim > > >  &/*duh*/,
			const std::vector< std::vector< Tensor< 2, dim > > >  &/*dduh*/,
			const std::vector< Point< dim > >                     &/*normals*/,
			std::vector< Vector< double > >                       &computed_quantities
			) const
{
  Assert(computed_quantities.size() == uh.size(), 
         ExcDimensionMismatch (computed_quantities.size(), uh.size()));

  for (unsigned int i=0; i<computed_quantities.size(); i++)
  {
    Assert(computed_quantities[i].size() == 3, 
           ExcDimensionMismatch (computed_quantities[i].size(), 3));
    Assert(uh[i].size() == 2, ExcDimensionMismatch (uh[i].size(), 2));

    computed_quantities[i](0) = uh[i](0);
    computed_quantities[i](1) = uh[i](1);
    computed_quantities[i](2) = sqrt(uh[i](0)*uh[i](0) + uh[i](1)*uh[i](1));
  }
}


template <int dim>
class UltrasoundProblem 
{
  public:
    UltrasoundProblem (ParameterHandler &);
    ~UltrasoundProblem ();
    void run ();

  private:
    static double get_omega (ParameterHandler &);
    static double get_c (ParameterHandler &);

    void make_grid ();
    void setup_system ();
    void assemble_system ();
    void solve ();
    void output_results () const;

    ParameterHandler       &prm; 

    Triangulation<dim>     triangulation;
    DoFHandler<dim>        dof_handler;
    FESystem<dim>          fe;

    SparsityPattern        sparsity_pattern;
    SparseMatrix<double>   system_matrix;	
    Vector<double>         solution, system_rhs;

    const double           c, omega;
};


template <int dim>
double
UltrasoundProblem<dim>::get_omega (ParameterHandler &prm)
{
  prm.enter_subsection ("Physical constants");
    double omega_tmp = prm.get_double("omega");
  prm.leave_subsection ();

  return omega_tmp;
}


template <int dim>
double
UltrasoundProblem<dim>::get_c (ParameterHandler &prm)
{
  prm.enter_subsection ("Physical constants");
    double c_tmp = prm.get_double("c");
  prm.leave_subsection ();

  return c_tmp;
}


template <int dim>
UltrasoundProblem<dim>::UltrasoundProblem (ParameterHandler&  param) 
  :
  prm(param),
  dof_handler(triangulation),
  fe(FE_Q<dim>(1), 2),
  c(get_c(prm)),
  omega(get_omega(prm))
{}


template <int dim>
UltrasoundProblem<dim>::~UltrasoundProblem () 
{
  dof_handler.clear();
}


template <int dim>
void UltrasoundProblem<dim>::make_grid ()
{
  prm.enter_subsection ("Mesh & geometry parameters");

  const double		focal_distance = prm.get_double("Focal distance");
  const unsigned int	N_ref          = prm.get_integer("Number of refinements");

  prm.leave_subsection ();

  GridGenerator::subdivided_hyper_cube (triangulation, 5, 0, 1);

  const Point<dim> 	transducer = (dim == 2) ? 
				     Point<dim> (0.5, 0.0) :
				     Point<dim> (0.5, 0.5, 0.0), 
			focal_point = (dim == 2) ?
				      Point<dim> (0.5, focal_distance) :
				      Point<dim> (0.5, 0.5, focal_distance);

  double radius = sqrt( (focal_point.distance(transducer) * 
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

  triangulation.refine_global (N_ref);

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

  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());

  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);
  system_rhs.reinit (dof_handler.n_dofs());
  solution.reinit (dof_handler.n_dofs());
}


template <int dim>
void UltrasoundProblem<dim>::assemble_system () 
{
  const double om2 = omega * omega;
  const double c2 = c * c;

  QGauss<dim>    quadrature_formula(2);
  QGauss<dim-1>  face_quadrature_formula(2);

  const unsigned int n_q_points	      = quadrature_formula.n_quadrature_points,
		     n_face_q_points  = face_quadrature_formula.n_quadrature_points,
		     dofs_per_cell    = fe.dofs_per_cell;

  FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);

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
        if (cell->face(face)->at_boundary() &&
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
				      c *
				      omega *
				      fe_face_values.JxW(q_point);
        }

        cell->get_dof_indices (local_dof_indices);

    for (unsigned int i=0; i<dofs_per_cell; ++i)
      for (unsigned int j=0; j<dofs_per_cell; ++j)
        system_matrix.add (local_dof_indices[i],
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
  SparseDirectUMFPACK  A_direct;

  A_direct.initialize(system_matrix);
  A_direct.vmult(solution,system_rhs);
}


template <int dim>
void UltrasoundProblem<dim>::output_results () const
{
  Postprocessor<dim> pproc;
  DataOut<dim> data_out;

  data_out.attach_dof_handler (dof_handler);

  prm.enter_subsection("Output parameters");

  const std::string output_file    = prm.get("Output file"),
		    output_format  = prm.get("Output format");

  data_out.parse_parameters(prm);

  prm.leave_subsection ();

  DataOutBase::OutputFormat  format = DataOutBase::parse_output_format(output_format);

  const std::string filename = output_file +
			       DataOutBase::default_suffix(format);

  std::ofstream output (filename.c_str());

  data_out.add_data_vector (solution, pproc);
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
  output_results ();
}


int main () 
{
  try
  {
    ParameterHandler  prm;

    ParameterReader   param(prm);
    param.read_parameters();

    Assert (DIM > 1, ExcNotImplemented());

    UltrasoundProblem<DIM>  ultrasound_problem (prm);
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
