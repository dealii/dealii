				 // @sect3{Include files}

				 // This program is an daptation of step-20
				 // and includes some technique of DG method from step-12
				 // We list include files in the order
				 // base-lac-grid-dofs-fe-numerics.
				 
#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <base/function.h>
#include <lac/block_vector.h>
#include <lac/full_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_tools.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_renumbering.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_constraints.h>
#include <fe/fe_q.h>
				 //The Discontinuous Galerkin finite element is declared: 
#include <fe/fe_dgq.h>

#include <fe/fe_system.h>
#include <fe/fe_values.h>
#include <fe/mapping_q1.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <fstream>
#include <iostream>
#include <sstream>
				 // The Raviart-Thomas finite element is declared:
#include <fe/fe_raviart_thomas.h>

				 // In this program, we use a tensorial
				 // coefficient. Since it may have a
				 // spatial dependence, we consider it
				 // a tensor-valued function. The
				 // following include file provides
				 // the ``TensorFunction'' class that
				 // offers such functionality:
#include <base/tensor_function.h>


                                 // @sect3{The ``TwoPhaseFlowProblem'' class template}
                                 

template <int dim>
class TwoPhaseFlowProblem 
{
  public:
    TwoPhaseFlowProblem (const unsigned int degree);
    void run ();
    
  private:
    void make_grid_and_dofs ();
    void assemble_system ();
    void solve ();
    void compute_errors () const;
    void output_results (const unsigned int timestep_number) const;

    Vector<double> evaluate_solution (const Point<dim> &point) const;
    
    const unsigned int   degree;
    
    Triangulation<dim>   triangulation;
    FESystem<dim>        fe;
    DoFHandler<dim>      dof_handler;


    BlockSparsityPattern      sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;

    const unsigned int n_refinement_steps;
    
    double time_step;
    double vis;    
    double vfs_out;
    double v_out;
 
    BlockVector<double>       solution;
    BlockVector<double>       old_solution;
    BlockVector<double>       system_rhs;
    
    
};


				 //{Right hand side, boundary values and initial values}
                                
				 // we define the template for pressure right-hand side(source function)
                                 //and boundary values for pressure and saturation
                                 // initial values for saturation.

template <int dim>
class RightHandSide : public Function<dim> 
{
  public:
    RightHandSide () : Function<dim>(1) {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};



template <int dim>
class PressureBoundaryValues : public Function<dim> 
{
  public:
    PressureBoundaryValues () : Function<dim>(1) {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};


template <int dim>
class SaturationBoundaryValues : public Function<dim> 
{
  public:
    SaturationBoundaryValues () : Function<dim>(dim+2) {};
    
    virtual void vector_value (const Point<dim> &p, 
			       Vector<double>   &value) const;
};


template <int dim>
class InitialValues : public Function<dim> 
{
  public:
    InitialValues () : Function<dim>(dim+2) {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p, 
			       Vector<double>   &value) const;

};




				 // And then we also have to define
				 // these respective functions, of
				 // course. Given our discussion in
				 // the introduction of how the
				 // solution should look like, the
				 // following computations should be
				 // straightforward:
template <int dim>
double RightHandSide<dim>::value (const Point<dim>  &/*p*/,
				  const unsigned int /*component*/) const 
{
  return 0;
}



template <int dim>
double PressureBoundaryValues<dim>::value (const Point<dim>  &p,
					   const unsigned int /*component*/) const 
{
  return 1-p[0];
}


template <int dim>
void
SaturationBoundaryValues<dim>::vector_value (const Point<dim> &p,
					     Vector<double>   &values) const 
{
  Assert (values.size() == dim+2,
	  ExcDimensionMismatch (values.size(), dim+2));

  for (unsigned int d=0; d<dim+1; ++d)
    values(d) = 0;

  if (p[0] == 0)
    values(dim+1) = 1;
  else
    values(dim+1) = 0;
}



template <int dim>
double InitialValues<dim>::value (const Point<dim>  &p,
				  const unsigned int component) const 
{
  if(component<dim+1)
    return 0;
  else 
    { 
      if(p[0]==0)return 1;
      else return 0;
    }
  
}


template <int dim>
void
InitialValues<dim>::vector_value (const Point<dim> &p,
				  Vector<double>   &values) const 
{
  Assert (values.size() == dim+2,
	  ExcDimensionMismatch (values.size(), dim+2));

  for (unsigned int d=0; d<dim+1; ++d)
    values(d) = 0;
  values(dim+1) = InitialValues::value(p,dim+1);
}




				 // @sect3{The inverse permeability tensor, Coefficient and inverse mobility scalar}

                                
                                 //For the inverse  permeability tensor,
                                 // ``KInverse''.As in introduction, '
                                 // assume the heterogeneous is isotropic,
                                 // so it is a scalar multipy the identity matrix.
				 //DealII has a base class not only for
                                 // scalar and generally vector-valued
                                 // functions (the ``Function'' base
                                 // class) but also for functions that
                                 // return tensors of fixed dimension
                                 // and rank, the ``TensorFunction''
                                 // template. Here, the function under
                                 // consideration returns a dim-by-dim
                                 // matrix, i.e. a tensor of rank 2
                                 // and dimension ``dim''. We then
                                 // choose the template arguments of
                                 // the base class appropriately.
                                 //
                                 // The interface that the
                                 // ``TensorFunction'' class provides
                                 // is essentially equivalent to the
                                 // ``Function'' class. In particular,
                                 // there exists a ``value_list''
                                 // function that takes a list of
                                 // points at which to evaluate the
                                 // function, and returns the values
                                 // of the function in the second
                                 // argument, a list of tensors:
template <int dim>
class KInverse : public TensorFunction<2,dim>
{
  public:
    KInverse ();
    
    virtual void value_list (const std::vector<Point<dim> > &points,
			     std::vector<Tensor<2,dim> >    &values) const;

  private:
    std::vector<Point<dim> > centers;
};


template <int dim>
KInverse<dim>::KInverse () 
{
  const unsigned int N = 40;
  centers.resize (N);
  for (unsigned int i=0; i<N; ++i)
    for (unsigned int d=0; d<dim; ++d)
      centers[i][d] = 2.*rand()/RAND_MAX-1;
}



template <int dim>
void
KInverse<dim>::value_list (const std::vector<Point<dim> > &points,
                           std::vector<Tensor<2,dim> >    &values) const
{
  Assert (points.size() == values.size(),
	  ExcDimensionMismatch (points.size(), values.size()));

  for (unsigned int p=0; p<points.size(); ++p)
    {
      values[p].clear ();

      double permeability = 0;
      for (unsigned int i=0; i<centers.size(); ++i)
        permeability += std::exp(-(points[p]-centers[i]).square()
                                 / (0.1 * 0.1));
      
      const double normalized_permeability
        = std::max(permeability, 0.005);
      
      for (unsigned int d=0; d<dim; ++d)
	values[p][d][d] = 10./normalized_permeability;
    }
}



double mobility_inverse (const double S, const double vis)
{ 
  return 1.0 /(1.0/vis * S * S + (1-S) * (1-S));
}

double f_saturation(const double S, const double vis)
{   

  return S*S /( S * S +vis * (1-S) * (1-S));
}





                                 // @sect4{extract_u and friends}

                                 // The next five functions are
                                 // needed for matrix and right hand
                                 // side assembly. They are described
                                 // in detail in step-20:
template <int dim>
Tensor<1,dim>
extract_u (const FEValuesBase<dim> &fe_values,
	   const unsigned int i,
	   const unsigned int q)
{
  Tensor<1,dim> tmp;

  for (unsigned int d=0; d<dim; ++d)
    tmp[d] = fe_values.shape_value_component (i,q,d);

  return tmp;
}



template <int dim>
double
extract_div_u (const FEValuesBase<dim> &fe_values,
	       const unsigned int i,
	       const unsigned int q)
{
  double divergence = 0;
  for (unsigned int d=0; d<dim; ++d)
    divergence += fe_values.shape_grad_component (i,q,d)[d];

  return divergence;
}


  
template <int dim>
double extract_p (const FEValuesBase<dim> &fe_values,
                  const unsigned int i,
                  const unsigned int q)
{
  return fe_values.shape_value_component (i,q,dim);
}

template <int dim>
double extract_s (const FEValuesBase<dim> &fe_values,
                  const unsigned int i,
                  const unsigned int q)
{
  return fe_values.shape_value_component (i,q,dim+1);
}

template <int dim>
Tensor<1,dim>
extract_grad_s(const FEValuesBase<dim> &fe_values,
	       const unsigned int i,
	       const unsigned int q)
{
  Tensor<1,dim> tmp;
  for (unsigned int d=0; d<dim; ++d)
    tmp[d] = fe_values.shape_grad_component (i,q,dim+1)[d];

  return tmp;
}



                                 // @sect3{TwoPhaseFlowProblem class implementation}

                                 // @sect4{TwoPhaseFlowProblem::TwoPhaseFlowProblem}
                                 //  we use RT(k) X DG(k),DG(k) spaces.
                                 // time_step is small enough to make the solution 
                                 // converges stably. 

                               
                                 
template <int dim>
TwoPhaseFlowProblem<dim>::TwoPhaseFlowProblem (const unsigned int degree)
		:
		degree (degree),
                fe (FE_RaviartThomas<dim>(degree), 1,
                    FE_DGQ<dim>(degree), 1,
		    FE_DGQ<dim>(degree), 1),
		dof_handler (triangulation),
		n_refinement_steps (5),
		time_step (10.0/std::pow(2.0, double(n_refinement_steps))/6),
                vis (0.2)
                
{}



                                 // @sect4{TwoPhaseFlowProblem::make_grid_and_dofs}

                                 // This next function starts out with
                                 // well-known functions calls that
                                 // create and refine a mesh, and then
                                 // associate degrees of freedom with
                                 // it:
template <int dim>
void TwoPhaseFlowProblem<dim>::make_grid_and_dofs ()
{
  GridGenerator::hyper_cube (triangulation, 0, 1);
  
  for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
    { if (triangulation.begin()->face(f)->center()[0] == 0)
      triangulation.begin()->face(f)->set_boundary_indicator (1);
      if (triangulation.begin()->face(f)->center()[0] == 1)
	triangulation.begin()->face(f)->set_boundary_indicator (2);
    }

  triangulation.refine_global (n_refinement_steps);
  
  dof_handler.distribute_dofs (fe);

 
  DoFRenumbering::component_wise (dof_handler);

                                  
  std::vector<unsigned int> dofs_per_component (dim+2);
  DoFTools::count_dofs_per_component (dof_handler, dofs_per_component);  
  const unsigned int n_u = dofs_per_component[0],
                     n_p = dofs_per_component[dim],
		     n_s = dofs_per_component[dim+1];

  std::cout << "Number of active cells: "
	    << triangulation.n_active_cells()
	    << std::endl
	    << "Total number of cells: "
	    << triangulation.n_cells()
	    << std::endl
            << "Number of degrees of freedom: "
	    << dof_handler.n_dofs()
            << " (" << n_u << '+' << n_p << '+'<< n_s <<')'
	    << std::endl;

  
  const unsigned int
    n_couplings = dof_handler.max_couplings_between_dofs();
  
  sparsity_pattern.reinit (3,3);
  sparsity_pattern.block(0,0).reinit (n_u, n_u, n_couplings);
  sparsity_pattern.block(1,0).reinit (n_p, n_u, n_couplings);
  sparsity_pattern.block(2,0).reinit (n_s, n_u, n_couplings);
  sparsity_pattern.block(0,1).reinit (n_u, n_p, n_couplings);
  sparsity_pattern.block(1,1).reinit (n_p, n_p, n_couplings);
  sparsity_pattern.block(2,1).reinit (n_s, n_p, n_couplings);
  sparsity_pattern.block(0,2).reinit (n_u, n_s, n_couplings);
  sparsity_pattern.block(1,2).reinit (n_p, n_s, n_couplings);
  sparsity_pattern.block(2,2).reinit (n_s, n_s, n_couplings);
  
  sparsity_pattern.collect_sizes();

  
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);

                                   
  solution.reinit (3);
  solution.block(0).reinit (n_u);
  solution.block(1).reinit (n_p);
  solution.block(2).reinit (n_s);
  solution.collect_sizes ();
  
  old_solution.reinit (3);
  old_solution.block(0).reinit (n_u);
  old_solution.block(1).reinit (n_p);
  old_solution.block(2).reinit (n_s);
  old_solution.collect_sizes ();
  
  system_rhs.reinit (3);
  system_rhs.block(0).reinit (n_u);
  system_rhs.block(1).reinit (n_p);
  system_rhs.block(2).reinit (n_s);
  system_rhs.collect_sizes ();


}


                                 // @sect4{TwoPhaseFlowProblem::assemble_system}
                                 // The function that
                                 // assembles the linear system has
                                 // mostly been discussed already in
                                 // the introduction to this
                                 // test case. We want to emphasize that
                                 // we assemble the first two equations
                                 // for velocity and pressure, but 
                                 // for saturation we only assemble 
                                 // the Matrixblock(2,2), for Matrixblock(0,2)
                                 // we will assemble it in "solve()", because
                                 //at that time, we have the new velocity solved
                                 // we can use it to assemble Matrixblock(0,2)
                    

template <int dim>
void TwoPhaseFlowProblem<dim>::assemble_system () 
{  
  QGauss<dim>   quadrature_formula(degree+2); 
  QGauss<dim-1> face_quadrature_formula(degree+2);

  FEValues<dim> fe_values (fe, quadrature_formula, 
			   update_values    | update_gradients |
                           update_q_points  | update_JxW_values);
  FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, 
				    update_values    | update_normal_vectors |
				    update_q_points  | update_JxW_values);

  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
  
  const unsigned int   n_q_points      = quadrature_formula.n_quadrature_points;
  const unsigned int   n_face_q_points = face_quadrature_formula.n_quadrature_points;

  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  
                                   // The next step is to declare
                                   // objects that represent the
                                   // source term, pressure boundary
                                   // value, and coefficient in the
                                   // equation. In addition to these
                                   // objects that represent
                                   // continuous functions, we also
                                   // need arrays to hold their values
                                   // at the quadrature points of
                                   // individual cells (or faces, for
                                   // the boundary values). Note that
                                   // in the case of the coefficient,
                                   // the array has to be one of
                                   // matrices.
  const RightHandSide<dim>          right_hand_side;
  const PressureBoundaryValues<dim> pressure_boundary_values;
  const KInverse<dim>               k_inverse;
   
  
  std::vector<double>               rhs_values (n_q_points);
  std::vector<double>               boundary_values (n_face_q_points);
  std::vector<Tensor<2,dim> >       k_inverse_values (n_q_points);
  
  std::vector<Vector<double> >      old_solution_values(n_q_points, Vector<double>(dim+2));
  std::vector<std::vector<Tensor<1,dim> > >  old_solution_grads(n_q_points,
                                                                std::vector<Tensor<1,dim> > (dim+2));
  
  

                                   // With all this in place, we can
                                   // go on with the loop over all
                                   // cells. The body of this loop has
                                   // been discussed in the
                                   // introduction, and will not be
                                   // commented any further here:
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  unsigned int cellnum=0;
  system_matrix=0;
  system_rhs=0;
  for (; cell!=endc; ++cell)
    { cellnum++;
      fe_values.reinit (cell);
      local_matrix = 0;
      local_rhs = 0;

      fe_values.get_function_values (old_solution, old_solution_values);
      right_hand_side.value_list (fe_values.get_quadrature_points(),
                                  rhs_values);
      k_inverse.value_list (fe_values.get_quadrature_points(),
                            k_inverse_values);
      
      for (unsigned int q=0; q<n_q_points; ++q)            
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
	    const double old_s = old_solution_values[q](dim+1);

            const Tensor<1,dim> phi_i_u = extract_u (fe_values, i, q);
	    const double div_phi_i_u = extract_div_u (fe_values, i, q);
            const double phi_i_p = extract_p (fe_values, i, q);
	    const double phi_i_s = extract_s (fe_values, i, q); 
	    const Tensor<1,dim> grad_phi_i_s = extract_grad_s(fe_values, i, q);
	    	     
            
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              {
                const Tensor<1,dim> phi_j_u = extract_u (fe_values, j, q);
		const double div_phi_j_u = extract_div_u (fe_values, j, q);		
                const double phi_j_p = extract_p (fe_values, j, q);
                const double phi_j_s = extract_s (fe_values, j, q);  
		
                local_matrix(i,j) += (phi_i_u * k_inverse_values[q] *
				      mobility_inverse(old_s,vis) * phi_j_u            
                                      - div_phi_i_u * phi_j_p
                                      - phi_i_p * div_phi_j_u
				      + phi_i_s * phi_j_s
		)
                                     * fe_values.JxW(q);     
              }

            local_rhs(i) += (-phi_i_p * rhs_values[q])*
                            fe_values.JxW(q);
          }
      
				       //here, we compute the boundary values for pressure 

      for (unsigned int face_no=0;
	   face_no<GeometryInfo<dim>::faces_per_cell;
	   ++face_no)
	if (cell->at_boundary(face_no))
	  {
	    fe_face_values.reinit (cell, face_no);
	    
	    pressure_boundary_values
	      .value_list (fe_face_values.get_quadrature_points(),
			   boundary_values);

	    for (unsigned int q=0; q<n_face_q_points; ++q) 
	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
		  const Tensor<1,dim>
		    phi_i_u = extract_u (fe_face_values, i, q);

		  local_rhs(i) += -(phi_i_u *
				    fe_face_values.normal_vector(q) *
				    boundary_values[q] *
				    fe_face_values.JxW(q));
		}
	  }

                                       // The final step in the loop
                                       // over all cells is to
                                       // transfer local contributions
                                       // into the global matrix and
                                       // right hand side vector. Note
                                       // that we use exactly the same
                                       // interface as in previous
                                       // examples, although we now
                                       // use block matrices and
                                       // vectors instead of the
                                       // regular ones. In other
                                       // words, to the outside world,
                                       // block objects have the same
                                       // interface as matrices and
                                       // vectors, but they
                                       // additionally allow to access
                                       // individual blocks.
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
      
        for (unsigned int j=0; j<dofs_per_cell; ++j)
	  {    system_matrix.add (local_dof_indices[i],
				  local_dof_indices[j],
				  local_matrix(i,j));
	  }
      
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        system_rhs(local_dof_indices[i]) += local_rhs(i);	
	
    }
}


                                 // @sect3{Linear solvers and preconditioners}

                                 // @sect4{The ``InverseMatrix'' class template}
                                 
				 // Everything here is completely same with step-20
                                 


template <class Matrix>
class InverseMatrix : public Subscriptor
{
  public:
    InverseMatrix (const Matrix &m);

    void vmult (Vector<double>       &dst,
                const Vector<double> &src) const;

  private:
    const SmartPointer<const Matrix> matrix;

    mutable GrowingVectorMemory<> vector_memory;    
};


template <class Matrix>
InverseMatrix<Matrix>::InverseMatrix (const Matrix &m)
                :
                matrix (&m)
{}


                                 
template <class Matrix>
void InverseMatrix<Matrix>::vmult (Vector<double>       &dst,
                                   const Vector<double> &src) const
{
  SolverControl solver_control (src.size(), 1e-8*src.l2_norm());
  SolverCG<> cg (solver_control, vector_memory);

  dst = 0;
  
  cg.solve (*matrix, dst, src, PreconditionIdentity());        
}


                                 // @sect4{The ``SchurComplement'' class template}
                                                                 
 
class SchurComplement : public Subscriptor
{
  public:
    SchurComplement (const BlockSparseMatrix<double> &A,
                     const InverseMatrix<SparseMatrix<double> > &Minv);

    void vmult (Vector<double>       &dst,
                const Vector<double> &src) const;

  private:
    const SmartPointer<const BlockSparseMatrix<double> > system_matrix;
    const SmartPointer<const InverseMatrix<SparseMatrix<double> > > m_inverse;
    
    mutable Vector<double> tmp1, tmp2;
};


SchurComplement::SchurComplement (const BlockSparseMatrix<double> &A,
                                  const InverseMatrix<SparseMatrix<double> > &Minv)
                :
                system_matrix (&A),
                m_inverse (&Minv),
                tmp1 (A.block(0,0).m()),
                tmp2 (A.block(0,0).m())
{}


void SchurComplement::vmult (Vector<double>       &dst,
                             const Vector<double> &src) const
{
  system_matrix->block(0,1).vmult (tmp1, src);
  m_inverse->vmult (tmp2, tmp1);
  system_matrix->block(1,0).vmult (dst, tmp2);
}


                                 // @sect4{The ``ApproximateSchurComplement'' class template}

class ApproximateSchurComplement : public Subscriptor
{
  public:
    ApproximateSchurComplement (const BlockSparseMatrix<double> &A);

    void vmult (Vector<double>       &dst,
                const Vector<double> &src) const;

  private:
    const SmartPointer<const BlockSparseMatrix<double> > system_matrix;
    
    mutable Vector<double> tmp1, tmp2;
};


ApproximateSchurComplement::ApproximateSchurComplement (const BlockSparseMatrix<double> &A)
                :
                system_matrix (&A),
                tmp1 (A.block(0,0).m()),
                tmp2 (A.block(0,0).m())
{}


void ApproximateSchurComplement::vmult (Vector<double>       &dst,
                                        const Vector<double> &src) const
{
  system_matrix->block(0,1).vmult (tmp1, src);
  system_matrix->block(0,0).precondition_Jacobi (tmp2, tmp1);
  system_matrix->block(1,0).vmult (dst, tmp2);
}



                                 // @sect4{TwoPhaseFlowProblem::solve}

                                 // After all these preparations,
                                 // we finally solves the linear
                                 // system for velocity and pressure.
                                 // And remember, we still have to assemble 
                                 // the Matirxbloc(2,0) after velocity is computed
                                 // , then use it to solve saturation.
template <int dim>
void TwoPhaseFlowProblem<dim>::solve () 
{
  const InverseMatrix<SparseMatrix<double> >
    m_inverse (system_matrix.block(0,0));
  Vector<double> tmp (solution.block(0).size());
  Vector<double> schur_rhs (solution.block(1).size());
  Vector<double> tmp2 (solution.block(2).size());
  

				   // this part is for pressure
  {
    m_inverse.vmult (tmp, system_rhs.block(0));
    system_matrix.block(1,0).vmult (schur_rhs, tmp);
    schur_rhs -= system_rhs.block(1);

    
    SchurComplement
      schur_complement (system_matrix, m_inverse);
    
    ApproximateSchurComplement
      approximate_schur_complement (system_matrix);
      
    InverseMatrix<ApproximateSchurComplement>
      preconditioner (approximate_schur_complement);

    
    SolverControl solver_control (system_matrix.block(0,0).m(),
				  1e-12*schur_rhs.l2_norm());
    SolverCG<>    cg (solver_control);

    cg.solve (schur_complement, solution.block(1), schur_rhs,
              preconditioner);
  
    std::cout << solver_control.last_step()
              << " CG Schur complement iterations to obtain convergence for pressure."
              << std::endl;
  }

                                   //  this part is for velocity. The
                                   // equation reads MU=-B^TP+F, and
                                   // we solve it by first computing
                                   // the right hand side, and then
                                   // multiplying it with the object
                                   // that represents the inverse of
                                   // the mass matrix:
  {
    system_matrix.block(0,1).vmult (tmp, solution.block(1));
    tmp *= -1;
    tmp += system_rhs.block(0);
    
    m_inverse.vmult (solution.block(0), tmp);
  }

				   //This part is for saturation.
				   // Here are many complicated functions
				   //which are very similiar with the
				   //assemble_system() part.
				   // For DG(0), we have to consider the discontinuty
				   // of the solution, then as in Introduction,
				   // compute numerical flux and judge it is in-flow or out-flow.
				   // After assemble Matrixbloc(2,0)
				   // , we could compute saturation directly. 
 
  { 
    QGauss<dim>   quadrature_formula(degree+2); 
    QGauss<dim-1> face_quadrature_formula(degree+2);  
    FEValues<dim> fe_values (fe, quadrature_formula, 
			     update_values    | update_gradients |
			     update_q_points  | update_JxW_values);
    FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, 
				      update_values    | update_normal_vectors |
				      update_q_points  | update_JxW_values);
    FEFaceValues<dim> fe_face_values_neighbor (fe, face_quadrature_formula, 
					       update_values);
  
 
    const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
    const unsigned int   n_q_points      = quadrature_formula.n_quadrature_points;
    const unsigned int   n_face_q_points = face_quadrature_formula.n_quadrature_points;
  
    vfs_out = 0.0;
    v_out = 0.0;  
  
    Vector<double>       local_rhs (dofs_per_cell);
    std::vector<Vector<double> > old_solution_values(n_q_points, Vector<double>(dim+2));
    std::vector<Vector<double> > old_solution_values_face(n_face_q_points, Vector<double>(dim+2));
    std::vector<Vector<double> > old_solution_values_face_neighbor(n_face_q_points, Vector<double>(dim+2));
    std::vector<Vector<double> > present_solution_values(n_q_points, Vector<double>(dim+2));
    std::vector<Vector<double> > present_solution_values_face(n_face_q_points, Vector<double>(dim+2));

    std::vector<double> neighbor_saturation (n_face_q_points);
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  
    typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();

    for (; cell!=endc; ++cell)
      {
	local_rhs = 0;
	fe_values.reinit (cell);

	fe_values.get_function_values (old_solution, old_solution_values);
	fe_values.get_function_values (solution, present_solution_values);
    
	for (unsigned int q=0; q<n_q_points; ++q) 
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    {
	      const double old_s = old_solution_values[q](dim+1);
	      Tensor<1,dim> present_u;
	      for (unsigned int d=0; d<dim; ++d)
		present_u[d] = present_solution_values[q](d);

	      const double phi_i_s = extract_s(fe_values, i, q);
	      const Tensor<1,dim> grad_phi_i_s = extract_grad_s(fe_values, i, q);
	    	     
	      local_rhs(i) += (
		time_step *(f_saturation(old_s,vis) * present_u * grad_phi_i_s)+
		old_s * phi_i_s)
			      * fe_values.JxW(q);
	    }
					 //Here is our numerical flux computation
					 // Finding neighbor as step-12
     		     		  
	for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell;++face_no)
	  {
	    fe_face_values.reinit (cell, face_no);

	    fe_face_values.get_function_values (old_solution, old_solution_values_face);
	    fe_face_values.get_function_values (solution, present_solution_values_face);

	    if (cell->at_boundary(face_no))
	      {
		if (cell->face(face_no)->boundary_indicator() == 1)
		  for (unsigned int q=0;q<n_face_q_points;++q)
		    neighbor_saturation[q] = 1;
		else
		  for (unsigned int q=0;q<n_face_q_points;++q)
		    neighbor_saturation[q] = 0;	                 
	      }
	    else
					       // there is a neighbor behind this face
	      {
		const typename DoFHandler<dim>::active_cell_iterator
		  neighbor = cell->neighbor(face_no);
		const unsigned int
		  neighbor_face = cell->neighbor_of_neighbor(face_no);

		fe_face_values_neighbor.reinit (neighbor, neighbor_face);
	     
		fe_face_values_neighbor.get_function_values (old_solution,
							     old_solution_values_face_neighbor);
	     
		for (unsigned int q=0;q<n_face_q_points;++q)
		  neighbor_saturation[q] = old_solution_values_face_neighbor[q](dim+1);
	      }
          

	    if (cell->at_boundary(face_no))
	      {	
		if (cell->face(face_no)->boundary_indicator() ==2 )
		  {for (unsigned int q=0;q<n_face_q_points;++q)
		    {
		      vfs_out += present_solution_values_face[q](0)
				 *f_saturation(present_solution_values_face[q](dim+1),vis)
				 *fe_face_values.JxW(q);
		      v_out += present_solution_values_face[q](0)
			       *fe_face_values.JxW(q);
		    }     	                     
		  }
	      }
	    for (unsigned int q=0;q<n_face_q_points;++q)
	      {
		Tensor<1,dim> present_u_face;
		for (unsigned int d=0; d<dim; ++d)
		  { present_u_face[d] = present_solution_values_face[q](d);
		  }
		const double normal_flux = present_u_face *
					   fe_face_values.normal_vector(q);

		const bool is_outflow_q_point = (normal_flux >= 0);
            	     	     	     
		if (is_outflow_q_point == true)
		  {
		    for (unsigned int i=0; i<dofs_per_cell; ++i)
		      { 
			const double outflow = -time_step * normal_flux 
					       * f_saturation(old_solution_values_face[q](dim+1),vis)
					       * extract_s(fe_face_values,i,q)
					       * fe_face_values.JxW(q);
			local_rhs(i) += outflow;
		      } 
		  }
             
		else
		  {
		    for (unsigned int i=0; i<dofs_per_cell; ++i)
		      {
			const double inflow = -time_step * normal_flux 
					      * f_saturation( neighbor_saturation[q],vis)
					      * extract_s(fe_face_values,i,q)
					      * fe_face_values.JxW(q);
			local_rhs(i) += inflow;
		      }
               
		  }
       
	      }
	      
	  }
  
	cell->get_dof_indices (local_dof_indices);
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	    system_rhs(local_dof_indices[i]) += local_rhs(i);
	  }
        	
      }	
    SolverControl solver_control (system_matrix.block(2,2).m(),
				  1e-12*system_rhs.block(2).l2_norm());
    SolverCG<>   cg (solver_control);
    cg.solve (system_matrix.block(2,2), solution.block(2), system_rhs.block(2),
	      PreconditionIdentity());
		
	
    std::cout << solver_control.last_step()
              << " CG iterations to obtain convergence for saturation."
              << std::endl;		
  } 

   
  old_solution = solution; 

   
 
  
}
                                 
                                 // @sect4{TwoPhaseFlow::compute_errors}

                                 // After we have dealt with the
                                 // linear solver and preconditioners,
                                 // we continue with the
                                 // implementation of our main
                                 // class. In particular, the next
                                 // task is to compute the errors in
                                 // our numerical solution, in both
                                 // the pressures velocities as well as
                                 // saturations.
                                 //
                                 // To compute errors in the solution,
                                 // we will not use ``VectorTools::integrate_difference''
                                 // as  step-20,  since we don't have exact solutions.
                                 // What we will do is to give some points
                                 // and evaluate the values on these points.
                                 //For every solution, we get values on those points,
                                 // then we can compare the values as an error.
   

template <int dim>
Vector<double>
TwoPhaseFlowProblem<dim>::evaluate_solution (const Point<dim> &point) const
{
  static const MappingQ1<dim> mapping;
				   // first find the cell in which this point
				   // is, initialize a quadrature rule with
				   // it, and then a FEValues object
  const typename DoFHandler<dim>::active_cell_iterator
    cell = GridTools::find_active_cell_around_point (dof_handler, point);

  const Point<dim> unit_point
    = mapping.transform_real_to_unit_cell(cell, point);
  Assert (GeometryInfo<dim>::is_inside_unit_cell (unit_point),
          ExcInternalError());

  const Quadrature<dim> quadrature (unit_point);
  FEValues<dim> fe_values(mapping, fe, quadrature, update_values);
  fe_values.reinit(cell);
                                   // then use this to get at the values of
                                   // the given fe_function at this point
  std::vector<Vector<double> > u_value(1, Vector<double>(dim+2));
  fe_values.get_function_values(solution, u_value);

  return u_value[0];
}

				 //{TwoPhaseFlowProblem::compute_errors}

				 // The compute_errors function is to compute
				 // error on some euqally spaced fixed points
				 // use evaluation function to interpret 
				 // solution value at the point
				 // then output those fixed points' value
				 // For each mesh, we can compare the output
				 // to estimate errors.
   
template <int dim>
void TwoPhaseFlowProblem<dim>::compute_errors () const
{
  std::ofstream sampled_solution ("sampled_solution");

  const double dx = 0.01;
  const double dy = 0.01;

  for (double x=0; x<=1; x+=dx)
    for (double y=0; y<=1; y+=dy)
      {
	const Point<dim> point(x,y);

	Vector<double> solution_at_point(dim+2);

	solution_at_point = evaluate_solution (point);

	sampled_solution << point << " ";
	for (unsigned int c=0; c<dim+2; ++c)
	  sampled_solution << solution_at_point(c) << " ";
	sampled_solution << std::endl;
      }
}

                                 // @sect4{TwoPhaseFlowProblem::output_results}

                                 // The output_results function is
                                 // the one in which we generate
                                 // graphical output.
template <int dim>
void TwoPhaseFlowProblem<dim>::output_results 
(const unsigned int timestep_number)  const
{  
  std::vector<std::string> solution_names;
  switch (dim)
    {
      case 2:
            solution_names.push_back ("u");
            solution_names.push_back ("v");
            solution_names.push_back ("p");
	    solution_names.push_back ("S");
            break;
            
      case 3:
            solution_names.push_back ("u");
            solution_names.push_back ("v");
            solution_names.push_back ("w");
            solution_names.push_back ("p");
	    solution_names.push_back ("S");
            break;
            
      default:
            Assert (false, ExcNotImplemented());
    }
  
  DataOut<dim> data_out;

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, solution_names);

  data_out.build_patches (degree+1);
  
  std::ostringstream filename;
  filename << "solution-"<< timestep_number;

  std::ofstream output (filename.str().c_str());
  data_out.write_gnuplot (output);

				   //data_out.write_vtk (output);
}


                                 // @sect4{TwoPhaseFlowProblem::run}

                                 // This is the final function of our
                                 // main class. It's only job is to
                                 // call the other functions in their order:
template <int dim>
void TwoPhaseFlowProblem<dim>::run () 
{
  std::cout<<"Solving problem in " <<dim << " space dimensions." << std::endl;
  
  make_grid_and_dofs();
  
  ConstraintMatrix constraints;
  constraints.close();

  std::list<double> production_rate;
  std::list<double> production_time;

  Vector<double> tmp (old_solution.size());
  VectorTools::project (dof_handler, constraints, QGauss<dim>(degree+2),InitialValues<dim>(),tmp);
  std::copy (tmp.begin(), tmp.end(), old_solution.begin());
  
  unsigned int timestep_number = 1;
  
  for ( double time = time_step; time <=1; time+=time_step,  timestep_number++)
    { 
      std::cout<< "Timestep_number = "<< timestep_number<<std::endl; 
      assemble_system ();
      solve ();
      output_results(timestep_number);

      production_time.push_back (time);
      production_rate.push_back (1.0 - vfs_out/v_out);
      std::cout<<"production_rate="<<production_rate.back()<<std::endl;       
    }

  std::ofstream production_history ("production_history");
  std::list<double>::iterator
    list_element = production_rate.begin(),
    time_element = production_time.begin();
  for (; list_element != production_rate.end(); ++list_element, ++time_element)
    production_history << *time_element << " " << *list_element << std::endl;
  
   
  compute_errors ();
}

    
                                 // @sect3{The ``main'' function}

				 // In the main function, we pass the
				 // degree of the finite element space
				 // to the constructor of the TwoPhaseFlowProblem
				 // (here, we use zero-th order elements).
int main () 
{
  try
    {
      deallog.depth_console (0);

      TwoPhaseFlowProblem<2> two_phase_flow_problem(0);
      two_phase_flow_problem.run ();
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
