/* $Id$ */
/* Author: Ralf Hartmann, University of Heidelberg, 2000 */

				 // The first few files have already
				 // been covered in previous examples
				 // and will thus not be further
				 // commented on.
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <lac/vector.h>
#include <lac/sparse_matrix.h>
#include <lac/vector_memory.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_values.h>
#include <numerics/data_out.h>
#include <grid/grid_out.h>
#include <grid/grid_refinement.h>
#include <numerics/error_estimator.h>

#include <fe/mapping_q1.h>
#include <fe/fe_dgq.h>
#include <lac/precondition_block.h>
#include <lac/solver_richardson.h>


#include <fstream>


template <int dim>
class Beta
{
  public:
    Beta () {};

    void value_list (const std::vector<Point<dim> > &points,
		     std::vector<Point<dim> > &values) const;
};


template <int dim>
class RHS:  public Function<dim>
{
  public:
    RHS() {};
    
    virtual void value_list (const std::vector<Point<dim> > &points,
			     std::vector<double> &values,
			     const unsigned int component=0) const;
};


template <int dim>
class BoundaryFunction:  public Function<dim>
{
  public:
    BoundaryFunction() {};
    
    virtual void value_list (const std::vector<Point<dim> > &points,
			     std::vector<double> &values,
			     const unsigned int component=0) const;
};

template <int dim>
class DGAssembler
{
  public:
    DGAssembler() {};

    void assemble_cell_term(const FEValuesBase<dim>& fe_v,
			    FullMatrix<double> &cell_matrix,
			    Vector<double> &cell_vector);
    
    void assemble_face_term(const FEFaceValuesBase<dim>& fe_v,
			    const FEFaceValuesBase<dim>& fe_v_neighbor,
			    FullMatrix<double> &cell_matrix,
			    FullMatrix<double> &cell_inflow_matrix,
			    Vector<double> &cell_vector);
    
  private:
    Beta<dim> beta_function;
    RHS<dim> rhs_function;
    BoundaryFunction<dim> boundary_function;
};

				 // The main class is again almost
				 // unchanged. Two additions, however,
				 // are made: we have added the
				 // ``refine'' function, which is used
				 // to adaptively refine the grid
				 // (instead of the global refinement
				 // in the previous examples), and a
				 // variable which will hold the
				 // constraints associated to the
				 // hanging nodes.
template <int dim>
class TransportProblem
{
  public:
    TransportProblem ();
    ~TransportProblem ();

    void run ();
    
  private:
    void setup_system ();
    void assemble_system ();
    void solve ();
    void refine_grid ();
    void output_results (const unsigned int cycle) const;

    Triangulation<dim>   triangulation;
    MappingQ1<dim>       mapping;
    
				     // We need a finite element
				     // again. This time, we will want
				     // to use quadratic polynomials
				     // (but this is only specified in
				     // the constructor):
    FE_DGQ<dim>          fe;
    DoFHandler<dim>      dof_handler;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;
    
    Vector<double>       solution;
    Vector<double>       right_hand_side;

    DGAssembler<dim>     dg_assembler;
};


template <>
void Beta<2>::value_list(const std::vector<Point<2> > &points,
			 std::vector<Point<2> > &values) const
{
  Assert(values.size()==points.size(), ExcDimensionMismatch(values.size(),points.size()));
  for (unsigned int i=0; i<points.size(); ++i)
    {
      const Point<2> &p=points[i];
      Point<2> &beta=values[i];

      beta(0)=-p(1);
      beta(1)=p(0);
      beta/=sqrt(beta*beta);
    }
}




template <int dim>
void RHS<dim>::value_list(const std::vector<Point<dim> > &,
			  std::vector<double> &values,
			  const unsigned int) const
{
  for (unsigned int i=0; i<values.size(); ++i)
    values[i]=0;
}




template <int dim>
void BoundaryFunction<dim>::value_list(const std::vector<Point<dim> > &points,
				       std::vector<double> &values,
				       const unsigned int) const
{
  Assert(values.size()==points.size(), ExcDimensionMismatch(values.size(),points.size()));
  for (unsigned int i=0; i<values.size(); ++i)
    {
      if (points[i](0)<0.5)
	values[i]=1.;
      else
	values[i]=0.;
    }
}




template <int dim>
void DGAssembler<dim>::assemble_cell_term(const FEValuesBase<dim>& fe_v,
					  FullMatrix<double> &cell_matrix,
					  Vector<double> &cell_vector)
{
  const vector<vector<Tensor<1,2> > > &grad_v = fe_v.get_shape_grads ();
  const FullMatrix<double> &v = fe_v.get_shape_values ();
  const vector<double> &JxW = fe_v.get_JxW_values ();

  vector<Point<dim> > beta (fe_v.n_quadrature_points);
  vector<double> rhs (fe_v.n_quadrature_points);
  
  beta_function.value_list (fe_v.get_quadrature_points(), beta);
  rhs_function.value_list (fe_v.get_quadrature_points(), rhs);
  
  for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
    for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i) 
      {
	for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
	  cell_matrix(i,j) += beta[point]*grad_v[j][point]*
			      v(i,point) *
			      JxW[point];
	
	cell_vector(i) += rhs[point] *v(i,point) * JxW[point];
      }
}


template <int dim>
void DGAssembler<dim>::assemble_face_term(const FEFaceValuesBase<dim>& fe_v,
					  const FEFaceValuesBase<dim>& fe_v_neighbor,      
					  FullMatrix<double> &cell_matrix,
					  FullMatrix<double> &cell_inflow_matrix,
					  Vector<double> &cell_vector)
{
  DoFHandler<dim>::face_iterator face=fe_v.get_face();
  
  const FullMatrix<double> &v = fe_v.get_shape_values ();
  const FullMatrix<double> &v_neighbor = fe_v_neighbor.get_shape_values ();  
  const vector<double> &JxW = fe_v.get_JxW_values ();
  const vector<Point<dim> > &normals = fe_v.get_normal_vectors ();

  vector<Point<dim> > beta (fe_v.n_quadrature_points);
  vector<double> g(fe_v.n_quadrature_points);
  
  beta_function.value_list (fe_v.get_quadrature_points(), beta);

  if (face->at_boundary())
    boundary_function.value_list (fe_v.get_quadrature_points(), g);

  for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
    {
      double beta_n=beta[point] * normals[point];
      if (beta_n<0)
	{
	  for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
	    {
	      for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
		cell_matrix(i,j) -= beta_n *
				    v(j,point) *
				    v(i,point) *
				    JxW[point];

	      if (!face->at_boundary())
		for (unsigned int k=0; k<fe_v_neighbor.dofs_per_cell; ++k)
		  cell_inflow_matrix(i,k) += beta_n *
					     v_neighbor(k,point) *
					     v(i,point) *
					     JxW[point];
	      else
		cell_vector(i) -= beta_n *
				  g[point] *
				  v(i,point) *
				  JxW[point];
	    }
	}
    }
}

  
template <int dim>
TransportProblem<dim>::TransportProblem () :
                fe (1),
		dof_handler (triangulation)
{}


template <int dim>
TransportProblem<dim>::~TransportProblem () 
{
  dof_handler.clear ();
};



template <int dim>
void TransportProblem<dim>::setup_system ()
{
				   // To distribute degrees of
				   // freedom, the ``dof_handler''
				   // variable takes only the finite
				   // element object. In this case, it
				   // will distribute four degrees of
				   // freedom per cell.
  dof_handler.distribute_dofs (fe);

  sparsity_pattern.reinit (dof_handler.n_dofs(),
			   dof_handler.n_dofs(),
			   dof_handler.max_couplings_between_dofs());
  DoFTools::make_flux_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();
  
  system_matrix.reinit (sparsity_pattern);

  solution.reinit (dof_handler.n_dofs());
  right_hand_side.reinit (dof_handler.n_dofs());
};



template <int dim>
void TransportProblem<dim>::assemble_system () 
{
				   // See Cockburn paper for the proper quadrature.
  QGauss2<dim>  quadrature;
  QGauss2<dim-1>  face_quadrature;
  
  const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
  vector<unsigned int> dofs (dofs_per_cell);
  vector<unsigned int> dofs_neighbor (dofs_per_cell);

  UpdateFlags update_flags = UpdateFlags(update_values
					 | update_gradients
					 | update_q_points
					 | update_JxW_values);
  
  UpdateFlags face_update_flags = UpdateFlags(update_values
					      | update_q_points
					      | update_JxW_values
					      | update_normal_vectors);
   

  FEValues<dim> fe_v (
    mapping, fe, quadrature, update_flags);
  FEFaceValues<dim> fe_v_face (
    mapping, fe, face_quadrature, face_update_flags);
  FESubfaceValues<dim> fe_v_subface (
    mapping, fe, face_quadrature, face_update_flags);
  FEFaceValues<dim> fe_v_face_neighbor (
    mapping, fe, face_quadrature, UpdateFlags(update_values | update_default));
  FESubfaceValues<dim> fe_v_subface_neighbor (
    mapping, fe, face_quadrature, UpdateFlags(update_values | update_default));

				   // includes the u and v terms
  FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
				   // includes u_hat and v terms
  FullMatrix<double> cell_inflow_matrix (dofs_per_cell, dofs_per_cell);

  Vector<double>  cell_vector (dofs_per_cell);

  DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
					endc = dof_handler.end();
  DoFHandler<dim>::face_iterator face;
  DoFHandler<dim>::cell_iterator neighbor;
  DoFHandler<dim>::cell_iterator neighbor_child;

  for (;cell!=endc; ++cell) 
    {
				       // re-init fe values for this cell
      fe_v.reinit (cell);

      cell_matrix.clear ();
      cell_vector.clear ();

      dg_assembler.assemble_cell_term(fe_v,
				      cell_matrix,
				      cell_vector);
      
      cell->get_dof_indices (dofs);

      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
	{
	  face = cell->face(face_no);
	  
	  cell_inflow_matrix.clear();

	  if (face->at_boundary())
	    {
	      fe_v_face.reinit (cell, face_no);

	      dg_assembler.assemble_face_term(fe_v_face,
					      fe_v_face,
					      cell_matrix,
					      cell_inflow_matrix,
					      cell_vector);
	    }
	  else // if (!face->at_boundary())
	    {
	      Assert (cell->neighbor(face_no).state() == valid, ExcInternalError());
	      neighbor = cell->neighbor(face_no);
	      
	      if (face->has_children())  // i.e. neighbor is one level more refined than cell
		{
						   // store which number #cell# has in the
						   // list of neighbors of #neighbor#
		  const unsigned int neighbor2=cell->neighbor_of_neighbor(face_no);
		  
		  
						   // loop over all subfaces
		  for (unsigned int subface_no=0; subface_no<GeometryInfo<dim>::subfaces_per_face;
		       ++subface_no)
		    {
						       // get an iterator pointing to the
						       // cell behind the present subface
		      neighbor_child = neighbor->child(GeometryInfo<dim>::
						       child_cell_on_face(neighbor2,subface_no));
		      Assert (neighbor_child->face(neighbor2) == face->child(subface_no),
			      ExcInternalError());
		      Assert (!neighbor_child->has_children(), ExcInternalError());

		      fe_v_subface.reinit (cell, face_no, subface_no);
		      fe_v_face_neighbor.reinit (neighbor_child, neighbor2);
	
		      cell_inflow_matrix.clear();
	    
		      dg_assembler.assemble_face_term(fe_v_subface,
						      fe_v_face_neighbor,
						      cell_matrix,
						      cell_inflow_matrix,
						      cell_vector);
		  
						       // get indices of dofs of neighbor_child cell
		      neighbor_child->get_dof_indices (dofs_neighbor);
		      						
						       // distribute cell matrix
		      for (unsigned int i=0; i<dofs_per_cell; ++i)
			for (unsigned int k=0; k<dofs_per_cell; ++k)
			  system_matrix.add(dofs[i], dofs_neighbor[k],
					    cell_inflow_matrix(i,k));
		    }
		}
	      else // if (!face->has_children())
		{
		  if (neighbor->level() == cell->level()) 
		    {
						       // store which number #cell# has in the
						       // list of neighbors of #neighbor#
		      const unsigned int neighbor2=cell->neighbor_of_neighbor(face_no);

		      fe_v_face.reinit (cell, face_no);
		      fe_v_face_neighbor.reinit (neighbor, neighbor2);
		      
		      dg_assembler.assemble_face_term(fe_v_face,
						      fe_v_face_neighbor,
						      cell_matrix,
						      cell_inflow_matrix,
						      cell_vector);
		    }
		  else // if (neighbor->level() < cell->level()) i.e. neighbor is one level coarser than cell
		    {
		      Assert(neighbor->level() < cell->level(), ExcInternalError());

		      const std::pair<unsigned int, unsigned int> faceno_subfaceno=
			cell->neighbor_of_coarser_neighbor(face_no);
		      const unsigned int neighbor_face_no=faceno_subfaceno.first,
				      neighbor_subface_no=faceno_subfaceno.second;

		      Assert (neighbor->neighbor(neighbor_face_no)
			      ->child(GeometryInfo<dim>::child_cell_on_face(
				face_no,neighbor_subface_no)) == cell, ExcInternalError());
			
						       // now 'neighbor_face_no' stores the number
						       // of a face in the list of faces of 'neighbor'.
						       // This face has got a subface that is 
						       // between 'cell' and 'neighbor'.
						       // 'neighbor_subface_no' stores the number
						       // of this subface in the list of subfaces of this
						       // face 'neighbor->face(neighbor_face_no)'
						       // that is between 'cell' and 'neighbor'
		      fe_v_face.reinit (cell, face_no);
		      fe_v_subface_neighbor.reinit (neighbor, neighbor_face_no,
						    neighbor_subface_no);
		      
		      dg_assembler.assemble_face_term(fe_v_face,
						      fe_v_subface_neighbor,
						      cell_matrix,
						      cell_inflow_matrix,
						      cell_vector);
		    } // else // if (neighbor->level() < cell->level())

						   // get indices of dofs of neighbor_child cell
		  neighbor->get_dof_indices (dofs_neighbor);
		 						
						   // distribute cell_inflow_matrix
		  for (unsigned int i=0; i<dofs_per_cell; ++i)
		    for (unsigned int k=0; k<dofs_per_cell; ++k)
		      system_matrix.add(dofs[i], dofs_neighbor[k],
					cell_inflow_matrix(i,k));
		} // else // if (!face->has_children())
	    }  // else // if (!face->at_boundary())
	} //for (face_no...)
      
				       // distribute cell matrix
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  system_matrix.add(dofs[i], dofs[j], cell_matrix(i,j));
      
				       // distribute cell vector
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	right_hand_side(dofs[i]) += cell_vector(i);
    }  // for (cell...)
};



template <int dim>
void TransportProblem<dim>::solve () 
{  
  SolverControl           solver_control (1000, 1e-12);
  PrimitiveVectorMemory<> vector_memory;
  SolverRichardson<>      solver (solver_control, vector_memory);

  PreconditionBlockSSOR<double> preconditioner;
  preconditioner.initialize(system_matrix, fe.dofs_per_cell);
  preconditioner.invert_diagblocks();
  
  solver.solve (system_matrix, solution, right_hand_side,
		preconditioner);
};


template <int dim>
void TransportProblem<dim>::refine_grid ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

  FunctionMap<dim>::type neumann_boundary;

  KellyErrorEstimator<dim>::estimate (dof_handler,
				      QGauss3<dim-1>(),
				      neumann_boundary,
				      solution,
				      estimated_error_per_cell);

  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
						   estimated_error_per_cell,
						   0.3, 0.03);

  triangulation.execute_coarsening_and_refinement ();
};



template <int dim>
void TransportProblem<dim>::output_results (const unsigned int cycle) const
{
				   // We want to write the grid in
				   // each cycle. Here is another way
				   // to quickly produce a filename
				   // based on the cycle number. It
				   // assumes that the numbers `0'
				   // through `9' are represented
				   // consecutively in the character
				   // set (which is the case in all
				   // known character sets). However,
				   // this will only work if the cycle
				   // number is less than ten, which
				   // we check by an assertion.
  std::string filename = "grid-";
  filename += ('0' + cycle);
  Assert (cycle < 10, ExcInternalError());
  
  filename += ".eps";
  std::ofstream eps_output (filename.c_str());

				   // Using this filename, we write
				   // each grid as a postscript file.
  GridOut grid_out;
  grid_out.write_eps (triangulation, eps_output);

				   // output of the solution
  filename = "sol-";
  filename += ('0' + cycle);
  Assert (cycle < 10, ExcInternalError());
  
  filename += ".gnuplot";
  std::ofstream gnuplot_output (filename.c_str());
  
  DataOut<dim> data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "u");

  data_out.build_patches ();
  
  data_out.write_gnuplot(gnuplot_output);
};



template <int dim>
void TransportProblem<dim>::run () 
{
  for (unsigned int cycle=0; cycle<3; ++cycle)
    {
      std::cout << "Cycle " << cycle << ':' << std::endl;

      if (cycle == 0)
	{
	  GridGenerator::hyper_cube (triangulation);

	  triangulation.refine_global (3);
	}
      else
					 // In case this is not the
					 // first cycle, we want to
					 // refine the grid. Unlike
					 // the global refinement
					 // employed in the last
					 // example, we now use the
					 // adaptive procedure
					 // described in the function
					 // which we now call:
	{
	  refine_grid ();
	};
      

      std::cout << "   Number of active cells:       "
		<< triangulation.n_active_cells()
		<< std::endl;

      setup_system ();

      std::cout << "   Number of degrees of freedom: "
		<< dof_handler.n_dofs()
		<< std::endl;
      
      assemble_system ();
      solve ();
      output_results (cycle);
    }
}

int main () 
{

				   // The general idea behind the
				   // layout of this function is as
				   // follows: let's try to run the
				   // program as we did before...
  try
    {
      TransportProblem<2> Transport_problem_2d;
      Transport_problem_2d.run ();
    }
				   // ...and if this should fail, try
				   // to gather as much information as
				   // possible. Specifically, if the
				   // exception that was thrown is an
				   // object of a class that is
				   // derived from the C++ standard
				   // class ``exception'', then we can
				   // use the ``what'' member function
				   // to get a string which describes
				   // the reason why the exception was
				   // thrown. 
				   //
				   // The deal.II exception classes
				   // are all derived from the
				   // standard class, and in
				   // particular, the ``exc.what()''
				   // function will return
				   // approximately the same string as
				   // would be generated if the
				   // exception was thrown using the
				   // ``Assert'' macro. You have seen
				   // the output of such an exception
				   // in the previous example, and you
				   // then know that it contains the
				   // file and line number of where
				   // the exception occured, and some
				   // other information. This is also
				   // what would be printed in the
				   // following.
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
				       // We can't do much more than
				       // printing as much information
				       // as we can get to, so abort
				       // with error:
      return 1;
    }
				   // If the exception that was thrown
				   // somewhere was not an object of a
				   // class derived from the standard
				   // ``exception'' class, then we
				   // can't do anything at all. We
				   // then simply print an error
				   // message and exit.
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
    };

				   // If we got to this point, there
				   // was no exception which
				   // propagated up to the main
				   // function (maybe there were some,
				   // but they were caught somewhere
				   // in the program or the
				   // library). Therefore, the program
				   // performed as was expected and we
				   // can return without error.
  return 0;
};
