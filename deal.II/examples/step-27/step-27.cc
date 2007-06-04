/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 2000 */

/*    $Id$       */
/*    Version: $Name$                                          */
/*                                                                */
/*    Copyright (C) 2000, 2001, 2002, 2003, 2004, 2006, 2007 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

                                 // @sect3{Include files}

				 // The first few files have already
				 // been covered in previous examples
				 // and will thus not be further
				 // commented on.
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>
#include <base/utilities.h>
#include <base/timer.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/compressed_sparsity_pattern.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <hp/dof_handler.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_reordering.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <hp/fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>

#include <fstream>
#include <iostream>

#include <fe/fe_q.h>
#include <grid/grid_out.h>
#include <dofs/dof_constraints.h>
#include <grid/grid_refinement.h>
#include <numerics/error_estimator.h>

#include <complex>


				 // Finally, this is as in previous
				 // programs:
using namespace dealii;

template <int dim>
class LaplaceProblem 
{
  public:
    LaplaceProblem (const bool condense_glob = true);
    ~LaplaceProblem ();

    void run ();
    
  private:
    void setup_system ();
    void assemble_system ();
    void solve ();
    void create_coarse_grid ();
    void refine_grid ();
    void estimate_smoothness (Vector<float> &smoothness_indicators) const;
    void output_results (const unsigned int cycle) const;

    Triangulation<dim>   triangulation;

    hp::DoFHandler<dim>      dof_handler;
    hp::FECollection<dim>    fe_collection;
    hp::QCollection<dim>     quadrature_collection;
    hp::QCollection<dim-1>   face_quadrature_collection;

    ConstraintMatrix     hanging_node_constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       system_rhs;

    Timer distr, condense, hang, assemble, solver;

    const unsigned int max_degree;
    const bool condense_global;
};



template <int dim>
class RightHandSide : public Function<dim>
{
  public:
    RightHandSide () : Function<dim> () {}
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component) const;
};


template <int dim>
double
RightHandSide<dim>::value (const Point<dim>   &p,
			   const unsigned int  /*component*/) const
{
  switch (dim)
    {
      case 2:
      {
	double product = 1;
	for (unsigned int d=0; d<dim; ++d)
	  product *= (p[d]+1);
	return product;
      }
      
      case 3:
	    return (p[0]>std::fabs(p[1]) ? 1 : 0);
	    
      default:
	    Assert (false, ExcNotImplemented());
    }
  return 0.;
}




template <int dim>
LaplaceProblem<dim>::LaplaceProblem (bool condense_glob) :
  dof_handler (triangulation),
  max_degree (dim == 2 ? 7 : 5),
  condense_global (condense_glob)
{
  for (unsigned int degree=2; degree<=max_degree; ++degree)
    {
      fe_collection.push_back (FE_Q<dim>(degree));
      quadrature_collection.push_back (QGauss<dim>(degree+2));
      face_quadrature_collection.push_back (QGauss<dim-1>(degree+2));
    }
}


template <int dim>
LaplaceProblem<dim>::~LaplaceProblem () 
{
  dof_handler.clear ();
}

template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  distr.reset();
  distr.start();
  dof_handler.distribute_dofs (fe_collection);
  distr.stop();

  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());

  hanging_node_constraints.clear ();
  hang.reset();
  hang.start();
  DoFTools::make_hanging_node_constraints (dof_handler,
					   hanging_node_constraints);

  hanging_node_constraints.close ();
  hang.stop();

  if (dim < 3)
    {
      sparsity_pattern.reinit (dof_handler.n_dofs(),
			       dof_handler.n_dofs(),
			       dof_handler.max_couplings_between_dofs());
      condense.reset();
      if (condense_global)
	{
	  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
	  condense.start();
	  hanging_node_constraints.condense (sparsity_pattern);
	  condense.stop();
	}
      else
	DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern,
					 hanging_node_constraints);

      sparsity_pattern.compress();      
    }
  else
    {
      CompressedSparsityPattern csp (dof_handler.n_dofs(),
				     dof_handler.n_dofs());
      condense.reset();
      if (condense_global)
	{
	  DoFTools::make_sparsity_pattern (dof_handler, csp);

	  condense.start();
	  hanging_node_constraints.condense (csp);
	  condense.stop();
	}
      else
	DoFTools::make_sparsity_pattern (dof_handler, csp,
					 hanging_node_constraints);
	
      sparsity_pattern.copy_from (csp);
    }

  system_matrix.reinit (sparsity_pattern);
}

template <int dim>
void LaplaceProblem<dim>::assemble_system () 
{
  assemble.reset ();
  assemble.start ();
  
  hp::FEValues<dim> hp_fe_values (fe_collection,
				  quadrature_collection, 
				  update_values    |  update_gradients |
				  update_q_points  |  update_JxW_values);

  const RightHandSide<dim> rhs_function;
  
  typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      const unsigned int   dofs_per_cell = cell->get_fe().dofs_per_cell;
      FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
      Vector<double>       cell_rhs (dofs_per_cell);

      std::vector<unsigned int> local_dof_indices (dofs_per_cell);

      cell_matrix = 0;
      cell_rhs = 0;

      hp_fe_values.reinit (cell);
      
      const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();

      std::vector<double>  rhs_values (fe_values.n_quadrature_points);
      rhs_function.value_list (fe_values.get_quadrature_points(),
			       rhs_values);
      
      for (unsigned int q_point=0; q_point<fe_values.n_quadrature_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      cell_matrix(i,j) += (fe_values.shape_grad(i,q_point) *
				   fe_values.shape_grad(j,q_point) *
				   fe_values.JxW(q_point));

	    cell_rhs(i) += (fe_values.shape_value(i,q_point) *
			    rhs_values[q_point] *
			    fe_values.JxW(q_point));
	  }

      cell->get_dof_indices (local_dof_indices);

      if (condense_global)
	{
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    {
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		system_matrix.add (local_dof_indices[i],
				   local_dof_indices[j],
				   cell_matrix(i,j));
	      
	      system_rhs(local_dof_indices[i]) += cell_rhs(i);
	    }
	}
      else
	{
	  hanging_node_constraints
	    .distribute_local_to_global (cell_matrix,
					 local_dof_indices,
					 system_matrix);
	  
	  hanging_node_constraints
	    .distribute_local_to_global (cell_rhs,
					 local_dof_indices,
					 system_rhs);
	}
    }

  assemble.stop();
  
  if (condense_global)
    {
      condense.start();
      hanging_node_constraints.condense (system_matrix);
      hanging_node_constraints.condense (system_rhs);
      condense.stop();
    }

  std::map<unsigned int,double> boundary_values;
  VectorTools::interpolate_boundary_values (dof_handler,
					    0,
					    ZeroFunction<dim>(),
					    boundary_values);
  MatrixTools::apply_boundary_values (boundary_values,
				      system_matrix,
				      solution,
				      system_rhs);
}

template <int dim>
void LaplaceProblem<dim>::solve () 
{
  SolverControl           solver_control (system_rhs.size(),
					  1e-8*system_rhs.l2_norm());
  SolverCG<>              cg (solver_control);

  PreconditionSSOR<> preconditioner;
  preconditioner.initialize(system_matrix, 1.2);

  cg.solve (system_matrix, solution, system_rhs,
	    preconditioner);

  hanging_node_constraints.distribute (solution);
}


unsigned int
int_pow (const unsigned int x,
	 const unsigned int n)
{
  unsigned int p=1;
  for (unsigned int i=0; i<n; ++i)
    p *= x;
  return p;
}


template <int dim>
void
LaplaceProblem<dim>::
estimate_smoothness (Vector<float> &smoothness_indicators) const
{
  const unsigned int N = max_degree;

				   // form all the Fourier vectors
				   // that we want to
				   // consider. exclude k=0 to avoid
				   // problems with |k|^{-mu} and also
				   // logarithms of |k|
  std::vector<Tensor<1,dim> > k_vectors;
  std::vector<unsigned int>   k_vectors_magnitude;
  switch (dim)
    {
      case 2:
      {
	for (unsigned int i=0; i<N; ++i)
	  for (unsigned int j=0; j<N; ++j)
	    if (!((i==0) && (j==0))
		&&
		(i*i + j*j < N*N))
	      {
		k_vectors.push_back (Point<dim>(deal_II_numbers::PI * i,
						deal_II_numbers::PI * j));
		k_vectors_magnitude.push_back (i*i+j*j);
	      }
	
	break;
      }

      case 3:
      {
	for (unsigned int i=0; i<N; ++i)
	  for (unsigned int j=0; j<N; ++j)
	    for (unsigned int k=0; k<N; ++k)
	      if (!((i==0) && (j==0) && (k==0))
		  &&
		  (i*i + j*j + k*k < N*N))
		{
		  k_vectors.push_back (Point<dim>(deal_II_numbers::PI * i,
						  deal_II_numbers::PI * j,
						  deal_II_numbers::PI * k));
		  k_vectors_magnitude.push_back (i*i+j*j+k*k);
	      }
	
	break;
      }
      
      default:
	    Assert (false, ExcNotImplemented());
    }

  const unsigned n_fourier_modes = k_vectors.size();
  std::vector<double> ln_k (n_fourier_modes);
  for (unsigned int i=0; i<n_fourier_modes; ++i)
    ln_k[i] = std::log (k_vectors[i].norm());
  

				   // assemble the matrices that do
				   // the Fourier transforms for each
				   // of the finite elements we deal
				   // with. note that these matrices
				   // are complex-valued, so we can't
				   // use FullMatrix
  QGauss<1>      base_quadrature (2);
  QIterated<dim> quadrature (base_quadrature, N);
  
  std::vector<Table<2,std::complex<double> > >
    fourier_transform_matrices (fe_collection.size());
  for (unsigned int fe=0; fe<fe_collection.size(); ++fe)
    {
      fourier_transform_matrices[fe].reinit (n_fourier_modes,
					     fe_collection[fe].dofs_per_cell);

      for (unsigned int k=0; k<n_fourier_modes; ++k)
	for (unsigned int i=0; i<fe_collection[fe].dofs_per_cell; ++i)
	  {
	    std::complex<double> sum = 0;
	    for (unsigned int q=0; q<quadrature.n_quadrature_points; ++q)
	      {
		const Point<dim> x_q = quadrature.point(q);
		sum += std::exp(std::complex<double>(0,1) *
				(k_vectors[k] * x_q)) *
		       fe_collection[fe].shape_value(i,x_q) *
		       quadrature.weight(q);
	      }
	    fourier_transform_matrices[fe](k,i)
	      = sum / std::pow(2*deal_II_numbers::PI, 1.*dim/2);
	  }
    }

				   // the next thing is to loop over
				   // all cells and do our work there,
				   // i.e. to locally do the Fourier
				   // transform and estimate the decay
				   // coefficient
  std::vector<std::complex<double> > fourier_coefficients (n_fourier_modes);
  Vector<double>                     local_dof_values;
  
  typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (unsigned int index=0; cell!=endc; ++cell, ++index)
    {
      local_dof_values.reinit (cell->get_fe().dofs_per_cell);
      cell->get_dof_values (solution, local_dof_values);

				       // first compute the Fourier
				       // transform of the local
				       // solution
      std::fill (fourier_coefficients.begin(), fourier_coefficients.end(), 0);
      for (unsigned int f=0; f<n_fourier_modes; ++f)
	for (unsigned int i=0; i<cell->get_fe().dofs_per_cell; ++i)
	  fourier_coefficients[f] += 
	    fourier_transform_matrices[cell->active_fe_index()](f,i)
	    *
	    local_dof_values(i);

				       // enter the Fourier
				       // coefficients into a map with
				       // the k-magnitudes, to make
				       // sure that we get only the
				       // largest magnitude for each
				       // value of |k|
      std::map<unsigned int, double> k_to_max_U_map;
      for (unsigned int f=0; f<n_fourier_modes; ++f)
	if ((k_to_max_U_map.find (k_vectors_magnitude[f]) ==
	     k_to_max_U_map.end())
	    ||
	    (k_to_max_U_map[k_vectors_magnitude[f]] <
	     std::abs (fourier_coefficients[f])))
	  k_to_max_U_map[k_vectors_magnitude[f]]
	    = std::abs (fourier_coefficients[f]);
      
				       // now we have to calculate the
				       // various contributions to the
				       // formula for mu. we'll only
				       // take those fourier
				       // coefficients with the
				       // largest value for a given
				       // |k|
      double  sum_1           = 0,
	      sum_ln_k        = 0,
	      sum_ln_k_square = 0,
	      sum_ln_U        = 0,
	      sum_ln_U_ln_k   = 0;
      for (unsigned int f=0; f<n_fourier_modes; ++f)
	if (k_to_max_U_map[k_vectors_magnitude[f]] ==
	    std::abs (fourier_coefficients[f]))
	  {
	    sum_1 += 1;
	    sum_ln_k += ln_k[f];
	    sum_ln_k_square += ln_k[f]*ln_k[f];
	    sum_ln_U += std::log (std::abs (fourier_coefficients[f]));
	    sum_ln_U_ln_k += std::log (std::abs (fourier_coefficients[f])) * ln_k[f];
	  }

      const double mu
	= (1./(sum_1*sum_ln_k_square - sum_ln_k*sum_ln_k)
	   *
	   (sum_ln_k*sum_ln_U - sum_1*sum_ln_U_ln_k));
      
      smoothness_indicators(index) = mu - 1.*dim/2;
    }
}



  
template <int dim>
void LaplaceProblem<dim>::refine_grid ()
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
  KellyErrorEstimator<dim>::estimate (dof_handler,
				      face_quadrature_collection,
				      typename FunctionMap<dim>::type(),
				      solution,
				      estimated_error_per_cell);

  Vector<float> smoothness_indicators (triangulation.n_active_cells());
  estimate_smoothness (smoothness_indicators);
  
  GridRefinement::refine_and_coarsen_fixed_number (triangulation,
						   estimated_error_per_cell,
						   0.3, 0.03);

  float max_smoothness = 0,
	 min_smoothness = smoothness_indicators.linfty_norm();
  {
    typename hp::DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
    for (unsigned int index=0; cell!=endc; ++cell, ++index)
      if (cell->refine_flag_set())
	{
	  max_smoothness = std::max (max_smoothness,
				     smoothness_indicators(index));
	  min_smoothness = std::min (min_smoothness,
				     smoothness_indicators(index));
	}
  }
  const float cutoff_smoothness = (max_smoothness + min_smoothness) / 2;
  {
    typename hp::DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
    for (unsigned int index=0; cell!=endc; ++cell, ++index)
      if (cell->refine_flag_set()
	  &&
	  (smoothness_indicators(index) > cutoff_smoothness)
	  &&
	  !(cell->active_fe_index() == fe_collection.size() - 1))
	{
	  cell->clear_refine_flag();
	  cell->set_active_fe_index (std::min (cell->active_fe_index() + 1,
					       fe_collection.size() - 1));
	}
  } 
  
  triangulation.execute_coarsening_and_refinement ();
}



template <int dim>
void LaplaceProblem<dim>::output_results (const unsigned int cycle) const
{
  Assert (cycle < 10, ExcNotImplemented());
  
  {
    const std::string filename = "grid-" +
				 Utilities::int_to_string (cycle, 2) +
				 ".eps";
    std::ofstream output (filename.c_str());
    
    GridOut grid_out;
    grid_out.write_eps (triangulation, output);
  }
  
  {
    Vector<float> smoothness_indicators (triangulation.n_active_cells());
    estimate_smoothness (smoothness_indicators);

    Vector<double> smoothness_field (dof_handler.n_dofs());
    DoFTools::distribute_cell_to_dof_vector (dof_handler,
					     smoothness_indicators,
					     smoothness_field);

    Vector<float> fe_indices (triangulation.n_active_cells());
    {
      typename hp::DoFHandler<dim>::active_cell_iterator
	cell = dof_handler.begin_active(),
	endc = dof_handler.end();
      for (unsigned int index=0; cell!=endc; ++cell, ++index)
	{
	  
	  fe_indices(index) = cell->active_fe_index();
//	  smoothness_indicators(index) *= std::sqrt(cell->diameter());
	}
    }
    
    const std::string filename = "solution-" +
				 Utilities::int_to_string (cycle, 2) +
				 ".gmv";
    DataOut<dim,hp::DoFHandler<dim> > data_out;

    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, "solution");
    data_out.add_data_vector (smoothness_indicators, "smoothness1");
    data_out.add_data_vector (smoothness_field, "smoothness2");
    data_out.add_data_vector (fe_indices, "fe_index");
    data_out.build_patches ();
  
    std::ofstream output (filename.c_str());
    data_out.write_gmv (output);
  }
}


template <>
void LaplaceProblem<2>::create_coarse_grid ()
{
  const unsigned int dim = 2;
  
  static const Point<2> vertices_1[]
    = {  Point<2> (-1.,   -1.),
         Point<2> (-1./2, -1.),
         Point<2> (0.,    -1.),
         Point<2> (+1./2, -1.),
         Point<2> (+1,    -1.),
	     
         Point<2> (-1.,   -1./2.),
         Point<2> (-1./2, -1./2.),
         Point<2> (0.,    -1./2.),
         Point<2> (+1./2, -1./2.),
         Point<2> (+1,    -1./2.),
	     
         Point<2> (-1.,   0.),
         Point<2> (-1./2, 0.),
         Point<2> (+1./2, 0.),
         Point<2> (+1,    0.),
	     
         Point<2> (-1.,   1./2.),
         Point<2> (-1./2, 1./2.),
         Point<2> (0.,    1./2.),
         Point<2> (+1./2, 1./2.),
         Point<2> (+1,    1./2.),
	     
         Point<2> (-1.,   1.),
         Point<2> (-1./2, 1.),
         Point<2> (0.,    1.),			  
         Point<2> (+1./2, 1.),
         Point<2> (+1,    1.)    };
  const unsigned int
    n_vertices = sizeof(vertices_1) / sizeof(vertices_1[0]);
  const std::vector<Point<dim> > vertices (&vertices_1[0],
                                           &vertices_1[n_vertices]);
  static const int cell_vertices[][GeometryInfo<dim>::vertices_per_cell]
    = {{0, 1, 5, 6},
       {1, 2, 6, 7},
       {2, 3, 7, 8},
       {3, 4, 8, 9},
       {5, 6, 10, 11},
       {8, 9, 12, 13},
       {10, 11, 14, 15},
       {12, 13, 17, 18},
       {14, 15, 19, 20},
       {15, 16, 20, 21},
       {16, 17, 21, 22},
       {17, 18, 22, 23}};
  const unsigned int
    n_cells = sizeof(cell_vertices) / sizeof(cell_vertices[0]);

  std::vector<CellData<dim> > cells (n_cells, CellData<dim>());
  for (unsigned int i=0; i<n_cells; ++i) 
    {
      for (unsigned int j=0;
           j<GeometryInfo<dim>::vertices_per_cell;
           ++j)
        cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    }

  triangulation.create_triangulation (vertices,
                                    cells,
                                    SubCellData());
  triangulation.refine_global (3);
}



namespace BreastPhantom
{
  

				   // Radius of the sphere of the breast
				   // phantom geometry
  static const double hemisphere_radius = 5;

				   // Radius of the disk underneath
  static const double bottom_disk_radius = 10;

				   // Bottom z-coordinate of the disk
				   // underneath
  const double bottom_disk_floor = -3;
				   // Top z-coordinate of the disk
				   // underneath
  const double bottom_disk_ceil = -.5;

				   // radius of the inner set of cells of
				   // the sphere geometry
  const double interior_hemisphere_radius
  = hemisphere_radius/(1.+std::sqrt(2.0));
      
  template <int dim>
  class SphereBoundary : public HyperBallBoundary<dim> 
  {
    public:
      SphereBoundary () 
		      :
		      HyperBallBoundary<dim> (Point<dim>(), hemisphere_radius) 
	{}
  };
  
  
  template <int dim>
  class CylinderBoundary : public StraightBoundary<dim>
  {
    public:
      typedef
      typename Triangulation<dim>::line_iterator
      line_iterator;
      
      typedef
      typename Triangulation<dim>::quad_iterator
      quad_iterator;
      
      typedef
      typename Triangulation<dim>::face_iterator
      face_iterator;

				       /**
					* Constructor.
					*/
      CylinderBoundary (const double radius)
		      :
		      radius (radius)
	{}
          
          
      virtual Point<dim>
      get_new_point_on_line (const line_iterator &line) const;

      virtual Point<dim>
      get_new_point_on_quad (const quad_iterator &quad) const;

      virtual void
      get_intermediate_points_on_line (const line_iterator &line,
				       std::vector<Point<dim> > &points) const;

      virtual void
      get_intermediate_points_on_quad (const quad_iterator &quad,
				       std::vector<Point<dim> > &points) const;

      virtual void
      get_normals_at_vertices (const face_iterator &face,
			       typename Boundary<dim>::FaceVertexNormals &face_vertex_normals) const;

    private:
      const double radius;
          
      void
      get_intermediate_points_between_points (const Point<dim> &p0, const Point<dim> &p1,
					      std::vector<Point<dim> > &points) const;    
  };

  
  template <>
  Point<3>
  CylinderBoundary<3>::
  get_new_point_on_line (const line_iterator &line) const
  {
    const Point<3> middle = StraightBoundary<3>::get_new_point_on_line (line);
				     // project to boundary
    Point<3> p(middle[0], middle[1], 0);
    p *= radius/std::sqrt(p.square());

    return Point<3> (p[0], p[1], middle[2]);
  }


  template<>
  Point<3>
  CylinderBoundary<3>::
  get_new_point_on_quad (const quad_iterator &quad) const
  {
    Point<3> middle = StraightBoundary<3>::get_new_point_on_quad (quad);
  
				     // project to boundary
    Point<3> p(middle[0], middle[1], 0);
    p *= radius/std::sqrt(p.square());

    return Point<3> (p[0], p[1], middle[2]);
  }


  template <int dim>
  void
  CylinderBoundary<dim>::
  get_intermediate_points_on_line (const line_iterator &line,
				   std::vector<Point<dim> > &points) const
  {
    if (points.size()==1)
      points[0] = get_new_point_on_line(line);
    else
      get_intermediate_points_between_points(line->vertex(0), line->vertex(1), points);
  }

  
  template <int dim>
  void
  CylinderBoundary<dim>::
  get_intermediate_points_between_points (const Point<dim> &,
					  const Point<dim> &,
					  std::vector<Point<dim> > &) const
  {
    Assert (false, ExcNotImplemented());
  }


  template <>
  void
  CylinderBoundary<3>::
  get_intermediate_points_on_quad (const Triangulation<3>::quad_iterator &,
				   std::vector<Point<3> > &) const
  {
    Assert (false, ExcNotImplemented());
  }


  template <int dim>
  void
  CylinderBoundary<dim>::
  get_normals_at_vertices (const typename Triangulation<dim>::face_iterator &,
			   typename Boundary<dim>::FaceVertexNormals &) const
  {
    Assert (false, ExcNotImplemented());
  }


  void
  create_coarse_grid (Triangulation<3> &coarse_grid)
  {
    const unsigned int dim = 3;
  
    std::vector<Point<dim> >    vertices;
    std::vector<CellData<dim> > cells;
    SubCellData                 sub_cell_data;

    const unsigned char
      bottom_cylinder_boundary_id  = 10,
      middle_cylinder_boundary_id  = 11,
      spherical_boundary_id = 12,
      straight_nondirichlet_boundary = 14;
        
        
				     // first build up the cells of the
				     // cylinder
    {
				       // the vertices in each plane of
				       // the cylinder are located on
				       // three concentric rings of radii
				       // interior_hemisphere_radius,
				       // hemisphere_radius, and
				       // bottom_disk_radius,
				       // respectively. first generate
				       // these three rings
      const Point<3> ring_points[8] = { Point<3>(-1,0,0),
					Point<3>(-1,-1,0) / std::sqrt(2.),
					Point<3>(0,-1,0),
					Point<3>(+1,-1,0) / std::sqrt(2.),
					Point<3>(+1,0,0),
					Point<3>(+1,+1,0) / std::sqrt(2.),
					Point<3>(0,+1,0),
					Point<3>(-1,+1,0) / std::sqrt(2.) };

				       // first the point in the middle
				       // and the rest of those on the
				       // upper surface
      vertices.push_back (Point<3>(0,0,bottom_disk_ceil));
      for (unsigned int ring=0; ring<3; ++ring)
	for (unsigned int i=0; i<8; ++i)
	  vertices.push_back (ring_points[i] * (ring == 0 ?
						interior_hemisphere_radius :
						(ring == 1 ? hemisphere_radius :
						 bottom_disk_radius))
			      +
			      Point<3>(0,0,bottom_disk_ceil));

				       // then points on lower surface
      vertices.push_back (Point<3>(0,0,bottom_disk_floor));
      for (unsigned int ring=0; ring<3; ++ring)
	for (unsigned int i=0; i<8; ++i)
	  vertices.push_back (ring_points[i] * (ring == 0 ?
						interior_hemisphere_radius :
						(ring == 1 ?
						 hemisphere_radius :
						 bottom_disk_radius))
			      +
			      Point<3>(0,0,bottom_disk_floor));

      const unsigned int n_vertices_per_surface = 25;
      Assert (vertices.size() == n_vertices_per_surface*2,
	      ExcInternalError());
    
				       // next create cells from these
				       // vertices. only store the
				       // vertices of the upper surface,
				       // the lower ones are the same
				       // +12
      {
	const unsigned int connectivity[20][4]
	  = { { 1, 2, 3, 0 },  // four cells in the center
	      { 3, 4, 5, 0 },
	      { 0, 5, 6, 7 },
	      { 1, 0, 7, 8 },            
          
	      { 9, 10, 2, 1 },  // eight cells of inner ring
	      { 10, 11, 3, 2 },
	      { 11, 12, 4, 3 },
	      { 4, 12, 13, 5 },
	      { 5, 13, 14, 6 },
	      { 6, 14, 15, 7 },
	      { 8, 7, 15, 16 },
	      { 9, 1, 8, 16 },

	      { 17, 18, 10, 9 },  // eight cells of outer ring
	      { 18, 19, 11, 10 },
	      { 19, 20, 12, 11 },
	      { 12, 20, 21, 13 },
	      { 13, 21, 22, 14 },
	      { 14, 22, 23, 15 },
	      { 16, 15, 23, 24 },
	      { 17, 9, 16, 24 }   };

					 // now create cells out of this
	for (unsigned int i=0; i<20; ++i)
	  {
	    CellData<3> cell;
	    for (unsigned int j=0; j<4; ++j)
	      {
		cell.vertices[j]   = connectivity[i][j];
		cell.vertices[j+4] = connectivity[i][j]+n_vertices_per_surface;
	      }
	    cell.material_id = 0;
	    cells.push_back (cell);
	  }
      }
    
				       // associate edges and faces on the
				       // outer boundary with boundary
				       // indicator of the cylinder
				       // boundary indicator. do this the
				       // same way as above, just this
				       // time with faces (edges follow
				       // from this immediately. some
				       // edges are duplicated since they
				       // belong to more than one cell,
				       // but that doesn't harm us here)
      {
	const unsigned int connectivity[8][2]
	  = { { 17,18 }, { 18, 19 }, { 19, 20 }, { 20, 21 },
	      { 21,22 }, { 22, 23 }, { 23, 24 }, { 24, 17 }};

	for (unsigned int i=0; i<8; ++i)
	  {
	    const CellData<2> face = 
	      { { connectivity[i][0]+n_vertices_per_surface,
		  connectivity[i][1]+n_vertices_per_surface,
		  connectivity[i][1],
		  connectivity[i][0] },
		bottom_cylinder_boundary_id };
	    sub_cell_data.boundary_quads.push_back (face);

	    const CellData<1> edges[4] = 
	      { { { connectivity[i][0],    connectivity[i][1]    },
		  bottom_cylinder_boundary_id },
		{ { connectivity[i][0]+n_vertices_per_surface,
		    connectivity[i][1]+n_vertices_per_surface },
		  bottom_cylinder_boundary_id },
		{ { connectivity[i][0]+n_vertices_per_surface,
		    connectivity[i][0]    },
		  bottom_cylinder_boundary_id },
		{ { connectivity[i][1]+n_vertices_per_surface,
		    connectivity[i][1]    },
		  bottom_cylinder_boundary_id } };
	    for (unsigned int i=0; i<4; ++i)
	      sub_cell_data.boundary_lines.push_back (edges[i]);
	  }
      }
    }

				     // next build up the middle ring. for
				     // this, copy the first 17 vertices
				     // up to z=0
    {
      const unsigned int first_upper_vertex = vertices.size();
          
      for (unsigned int i=0; i<17; ++i)
	vertices.push_back (Point<3>(vertices[i][0], vertices[i][1], 0));

				       // next create cells from these
				       // vertices. only store the
				       // vertices of the lower surface,
				       // the lower ones are the same
				       // +first_upper_vertex
      const unsigned int connectivity[12][4]
	= { { 1, 2, 3, 0 },  // four cells in the center
	    { 3, 4, 5, 0 },
	    { 0, 5, 6, 7 },
	    { 1, 0, 7, 8 },            
          
	    { 9, 10, 2, 1 },  // eight cells of ring
	    { 10, 11, 3, 2 },
	    { 11, 12, 4, 3 },
	    { 4, 12, 13, 5 },
	    { 5, 13, 14, 6 },
	    { 6, 14, 15, 7 },
	    { 8, 7, 15, 16 },
	    { 9, 1, 8, 16 }};
				       // now create cells out of this
      for (unsigned int i=0; i<12; ++i)
	{
	  CellData<3> cell;
	  for (unsigned int j=0; j<4; ++j)
	    {
	      cell.vertices[j]   = connectivity[i][j]+first_upper_vertex;
	      cell.vertices[j+4] = connectivity[i][j];
	    }
	  cell.material_id = 0;
	  cells.push_back (cell);
	}

				       // mark the 8 vertical edges with
				       // the correct boundary indicator
      for (unsigned int i=0; i<8; ++i)
	{
	  const CellData<1> edge = { { 9, 9+first_upper_vertex },
				     middle_cylinder_boundary_id };
	  sub_cell_data.boundary_lines.push_back (edge);
	}
				       // likewise with the 8 tangential
				       // edges on the lower disk. the
				       // edges at the interface between
				       // the middle disk and the
				       // hemisphere are handled by the
				       // hemisphere boundary
      for (unsigned int i=0; i<8; ++i)
	{
	  const CellData<1> edge = { { 9+i, 9+(i+1)%8},
				     middle_cylinder_boundary_id };
	  sub_cell_data.boundary_lines.push_back (edge);
	}

				       // then assign face indicators
      for (unsigned int i=0; i<8; ++i)
	{
	  const CellData<2> face = { { 9+i,
				       9+(i+1)%8,
				       9+(i+1)%8+first_upper_vertex,
				       9+i+first_upper_vertex},
				     middle_cylinder_boundary_id };
	  sub_cell_data.boundary_quads.push_back (face);
	}
    }
        
				     // the final part is setting the
				     // half-sphere on top of this
    {
				       // add four cubes to the top of
				       // the inner four cells, as well
				       // as 8 to their outside
      {
					 // mirror the first nine vertices
					 // above the surface, and scale
					 // them to a certain distance
					 // outward
	const double rx = hemisphere_radius / (1+std::sqrt(3.0));
	for (unsigned int i=0; i<9; ++i)
	  {
	    const Point<3> p (vertices[i][0],
			      vertices[i][1],
			      i == 0 ?
			      1
			      :
			      std::max(std::fabs(vertices[i][0]),
				       std::fabs(vertices[i][1])));
	    vertices.push_back (p / std::sqrt(p.square()) * rx);
	  }
	Assert (vertices.size() == 76, ExcInternalError());

					 // same with the next ring of
					 // vertices, except that they
					 // go to hemisphere_radius
	for (unsigned int i=9; i<17; ++i)
	  {
	    Point<3> p (vertices[i][0],
			vertices[i][1],
			std::max(std::fabs(vertices[i][0]),
				 std::fabs(vertices[i][1])));
	    vertices.push_back (p / std::sqrt(p.square()) *
				hemisphere_radius);
	  }
	Assert (vertices.size() == 84, ExcInternalError());
      
					 // make 12 cells out of this
	const unsigned int connectivity[12][4]
	  = { { 1, 2, 3, 0 },  // four cells in the center
	      { 3, 4, 5, 0 },
	      { 0, 5, 6, 7 },
	      { 1, 0, 7, 8 },

	      { 9, 10, 2, 1 },  // eight cells of inner ring
	      { 10, 11, 3, 2 },
	      { 11, 12, 4, 3 },
	      { 4, 12, 13, 5 },
	      { 5, 13, 14, 6 },
	      { 6, 14, 15, 7 },
	      { 8, 7, 15, 16 },
	      { 9, 1, 8, 16 },
	  };

	for (unsigned int i=0; i<12; ++i)
	  {
	    CellData<3> cell;
	    for (unsigned int j=0; j<4; ++j)
	      {
		cell.vertices[j]   = connectivity[i][j]+67;
		cell.vertices[j+4] = connectivity[i][j]+50;
	      }
	    cell.material_id = 0;
	    cells.push_back (cell);
	  }
      }

				       // assign boundary indicators to
				       // the faces and edges of these
				       // cells
      {
					 // these are the numbers of the
					 // vertices on the top surface
					 // of the cylinder, with one
					 // "wrap-around":
	const unsigned int vertices[9] = 
	  { 9, 10, 11, 12, 13, 14, 15, 16, 9 };
					 // their counter-parts are the
					 // same +67
	for (unsigned int i=0; i<8; ++i)
	  {
					     // generate a face
	    const CellData<2> face = 
	      { { vertices[i]+50,   vertices[i+1]+50,
		  vertices[i+1]+67, vertices[i]+67 },
		spherical_boundary_id };
	    sub_cell_data.boundary_quads.push_back (face);

					     // same for the faces
	    const CellData<1> edges[4] =
	      { { { vertices[i]+50,   vertices[i+1]+50 },
		  spherical_boundary_id },
		{ { vertices[i]+67,   vertices[i+1]+67 },
		  spherical_boundary_id },
		{ { vertices[i]+50,   vertices[i]+67   },
		  spherical_boundary_id },
		{ { vertices[i+1]+50, vertices[i+1]+67 },
		  spherical_boundary_id } };
	    for (unsigned int j=0; j<4; ++j)
	      sub_cell_data.boundary_lines.push_back (edges[j]);
	  }
      }  


				       // finally top the building
				       // with four closing cells and
				       // the vertex at the top
      {
	vertices.push_back (Point<3> (0,0,hemisphere_radius));

	const unsigned int connectivity[4][8]
	  = { { 59, 60, 61, 67,   51, 52, 53, 50 },
	      { 61, 62, 63, 67,   53, 54, 55, 50 },
	      { 67, 63, 64, 65,   50, 55, 56, 57 },
	      { 59, 67, 65, 66,   51, 50, 57, 58 }};
      
	for (unsigned int i=0; i<4; ++i)
	  {
	    CellData<3> cell;
	    for (unsigned int j=0; j<8; ++j)
	      cell.vertices[j]   = connectivity[i][j]+17;
	    cell.material_id   = 0;
	    cells.push_back (cell);
	  }

					 // generate boundary
					 // information for these cells,
					 // too
	for (unsigned int i=0; i<4; ++i)
	  {
	    const CellData<2> face = 
	      { { connectivity[i][0]+17, connectivity[i][1]+17,
		  connectivity[i][2]+17, connectivity[i][3]+17 },
		spherical_boundary_id };
	    sub_cell_data.boundary_quads.push_back (face);

	    const CellData<1> edges[4] =
	      { { { connectivity[i][0]+17, connectivity[i][1]+17 },
		  spherical_boundary_id },
		{ { connectivity[i][1]+17, connectivity[i][2]+17 },
		  spherical_boundary_id },
		{ { connectivity[i][2]+17, connectivity[i][3]+17 },
		  spherical_boundary_id },
		{ { connectivity[i][3]+17, connectivity[i][0]+17 },
		  spherical_boundary_id } };
	    for (unsigned int j=0; j<4; ++j)
	      sub_cell_data.boundary_lines.push_back (edges[j]);
	  }
      }
    }
  

				     // finally generate a triangulation
				     // out of this
    GridReordering<3>::reorder_cells (cells);
    coarse_grid.create_triangulation_compatibility (vertices, cells,
						    sub_cell_data);

				     // then associate boundary objects
				     // with the different boundary
				     // indicators
    static const CylinderBoundary<3>
      bottom_cylinder_boundary (bottom_disk_radius);
    static const CylinderBoundary<3>
      middle_cylinder_boundary (hemisphere_radius);
    static const SphereBoundary<3> sphere_boundary;

    coarse_grid.set_boundary (bottom_cylinder_boundary_id,
			      bottom_cylinder_boundary);
    coarse_grid.set_boundary (middle_cylinder_boundary_id,
			      middle_cylinder_boundary);
    coarse_grid.set_boundary (spherical_boundary_id,
			      sphere_boundary);

    for (Triangulation<dim>::active_cell_iterator cell=coarse_grid.begin_active();
	 cell != coarse_grid.end(); ++cell)
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
	if ((cell->face(f)->boundary_indicator() == 0)
	    &&
	    (cell->face(f)->center()[2] >= (bottom_disk_floor+bottom_disk_ceil)/2))
	  cell->face(f)->set_boundary_indicator(straight_nondirichlet_boundary);
  }
}



template <>
void LaplaceProblem<3>::create_coarse_grid ()
{
  BreastPhantom::create_coarse_grid (triangulation);
}



template <int dim>
void LaplaceProblem<dim>::run () 
{
  for (unsigned int cycle=0; cycle<30; ++cycle)
    {
      std::cout << "Cycle " << cycle << ':' << std::endl;

      if (cycle == 0)
	create_coarse_grid ();
      else
	refine_grid ();
      

      std::cout << "   Number of active cells:       "
		<< triangulation.n_active_cells()
		<< std::endl;

      Timer all;
      all.reset();
      all.start();
      setup_system ();

      std::cout << "   Number of degrees of freedom: "
		<< dof_handler.n_dofs()
		<< std::endl;
      std::cout << "   Number of constraints       : "
		<< hanging_node_constraints.n_constraints()
		<< std::endl;

      assemble_system ();
      

      solver.reset();
      solver.start();
      solve ();
      solver.stop();
      
      all.stop();

      std::cout << "   All: " << all()
		<< ", distr: " << distr()
		<< ", hang: " << hang()
		<< ", condense: " << condense()
		<< ", assemble: " << assemble()
		<< ", solver: " << solver()
		<< std::endl;
      
      output_results (cycle);
    }
}

int main () 
{
  try
    {
      deallog.depth_console (0);

      LaplaceProblem<3> laplace_problem (false);
      laplace_problem.run ();
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
