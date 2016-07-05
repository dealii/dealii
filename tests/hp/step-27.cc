/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2006 - 2015 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * --------------------------------------------------------------------
 */


// A combination of step-27 from 8.4 with corrected k-vectors, that is 2\pi*k instead of \pi*k
// and a new step-27 from 8.5 which use FESeries namespace. By default, the new
// version is used, but the blessed output file is obtained using the
// modified 8.4 version.


#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_series.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_values.h>

#include <fstream>
#include <iostream>
#include <complex>


namespace Step27
{
  using namespace dealii;



  template <int dim>
  class LaplaceProblem
  {
  public:
    LaplaceProblem ();
    ~LaplaceProblem ();

    void run ();

  private:
    void setup_system ();
    void assemble_system ();
    void solve ();
    void create_coarse_grid ();
    void estimate_smoothness (Vector<float> &smoothness_indicators);
    void postprocess (const unsigned int cycle);

    Triangulation<dim>   triangulation;

    hp::DoFHandler<dim>      dof_handler;
    hp::FECollection<dim>    fe_collection;
    hp::QCollection<dim>     quadrature_collection;
    hp::QCollection<dim-1>   face_quadrature_collection;

    hp::QCollection<dim> fourier_q_collection;
    std_cxx11::shared_ptr<FESeries::Fourier<dim> > fourier;
    std::vector<double> ln_k;
    Table<dim,std::complex<double> > fourier_coefficients;

    ConstraintMatrix     constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double>       solution;
    Vector<double>       system_rhs;

    const unsigned int max_degree;
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
    double product = 1;
    for (unsigned int d=0; d<dim; ++d)
      product *= (p[d]+1);
    return product;
  }


  template <typename T>
  void resize(Table<2,T> &coeff, const unsigned int N)
  {
    coeff.reinit(N,N);
  }


  template <int dim>
  LaplaceProblem<dim>::LaplaceProblem ()
    :
    dof_handler (triangulation),
    max_degree (dim <= 2 ? 7 : 5)
  {
    for (unsigned int degree=2; degree<=max_degree; ++degree)
      {
        fe_collection.push_back (FE_Q<dim>(degree));
        quadrature_collection.push_back (QGauss<dim>(degree+1));
        face_quadrature_collection.push_back (QGauss<dim-1>(degree+1));
      }

    const unsigned int N = max_degree;

    QGauss<1>            base_quadrature (2);
    QIterated<dim>       quadrature (base_quadrature, N);
    for (unsigned int i = 0; i < fe_collection.size(); i++)
      fourier_q_collection.push_back(quadrature);

    fourier = std_cxx11::make_shared<FESeries::Fourier<dim> >(N,
                                                              fe_collection,
                                                              fourier_q_collection);
    resize(fourier_coefficients,N);
  }



  template <int dim>
  LaplaceProblem<dim>::~LaplaceProblem ()
  {
    dof_handler.clear ();
  }


  template <int dim>
  void LaplaceProblem<dim>::setup_system ()
  {
    dof_handler.distribute_dofs (fe_collection);

    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());

    constraints.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler,
                                             constraints);
    VectorTools::interpolate_boundary_values (dof_handler,
                                              0,
                                              ZeroFunction<dim>(),
                                              constraints);
    constraints.close ();

    DynamicSparsityPattern dsp (dof_handler.n_dofs(),
                                dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints, false);
    sparsity_pattern.copy_from (dsp);

    system_matrix.reinit (sparsity_pattern);
  }




  template <int dim>
  void LaplaceProblem<dim>::assemble_system ()
  {
    hp::FEValues<dim> hp_fe_values (fe_collection,
                                    quadrature_collection,
                                    update_values    |  update_gradients |
                                    update_quadrature_points  |  update_JxW_values);

    const RightHandSide<dim> rhs_function;

    FullMatrix<double>   cell_matrix;
    Vector<double>       cell_rhs;

    std::vector<types::global_dof_index> local_dof_indices;

    typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        const unsigned int   dofs_per_cell = cell->get_fe().dofs_per_cell;

        cell_matrix.reinit (dofs_per_cell, dofs_per_cell);
        cell_matrix = 0;

        cell_rhs.reinit (dofs_per_cell);
        cell_rhs = 0;

        hp_fe_values.reinit (cell);

        const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();

        std::vector<double>  rhs_values (fe_values.n_quadrature_points);
        rhs_function.value_list (fe_values.get_quadrature_points(),
                                 rhs_values);

        for (unsigned int q_point=0;
             q_point<fe_values.n_quadrature_points;
             ++q_point)
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

        local_dof_indices.resize (dofs_per_cell);
        cell->get_dof_indices (local_dof_indices);

        constraints.distribute_local_to_global (cell_matrix, cell_rhs,
                                                local_dof_indices,
                                                system_matrix, system_rhs);
      }
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

    constraints.distribute (solution);
  }




  template <int dim>
  void LaplaceProblem<dim>::postprocess (const unsigned int cycle)
  {
    Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
    KellyErrorEstimator<dim>::estimate (dof_handler,
                                        face_quadrature_collection,
                                        typename FunctionMap<dim>::type(),
                                        solution,
                                        estimated_error_per_cell);


    Vector<float> smoothness_indicators (triangulation.n_active_cells());
    estimate_smoothness (smoothness_indicators);

    // Output to VTK
    if (false)
      {
        Vector<float> fe_degrees (triangulation.n_active_cells());
        {
          typename hp::DoFHandler<dim>::active_cell_iterator
          cell = dof_handler.begin_active(),
          endc = dof_handler.end();
          for (; cell!=endc; ++cell)
            fe_degrees(cell->active_cell_index())
              = fe_collection[cell->active_fe_index()].degree;
        }

        DataOut<dim,hp::DoFHandler<dim> > data_out;

        data_out.attach_dof_handler (dof_handler);
        data_out.add_data_vector (solution, "solution");
        data_out.add_data_vector (estimated_error_per_cell, "error");
        data_out.add_data_vector (smoothness_indicators, "smoothness");
        data_out.add_data_vector (fe_degrees, "fe_degree");
        data_out.build_patches ();

        const std::string filename = "solution-" +
                                     Utilities::int_to_string (cycle, 2) +
                                     ".vtk";
        std::ofstream output (filename.c_str());
        data_out.write_vtk (output);
      }

    {
      GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                       estimated_error_per_cell,
                                                       0.3, 0.03);

      float max_smoothness = *std::min_element (smoothness_indicators.begin(),
                                                smoothness_indicators.end()),
                             min_smoothness = *std::max_element (smoothness_indicators.begin(),
                                                                 smoothness_indicators.end());
      {
        typename hp::DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
        for (; cell!=endc; ++cell)
          if (cell->refine_flag_set())
            {
              max_smoothness = std::max (max_smoothness,
                                         smoothness_indicators(cell->active_cell_index()));
              min_smoothness = std::min (min_smoothness,
                                         smoothness_indicators(cell->active_cell_index()));
            }
      }
      const float threshold_smoothness = (max_smoothness + min_smoothness) / 2;

      {
        typename hp::DoFHandler<dim>::active_cell_iterator
        cell = dof_handler.begin_active(),
        endc = dof_handler.end();
        for (; cell!=endc; ++cell)
          if (cell->refine_flag_set()
              &&
              (smoothness_indicators(cell->active_cell_index()) > threshold_smoothness)
              &&
              (cell->active_fe_index()+1 < fe_collection.size()))
            {
              cell->clear_refine_flag();
              cell->set_active_fe_index (cell->active_fe_index() + 1);
            }
      }

      triangulation.execute_coarsening_and_refinement ();
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
           Point<2> (+1,    1.)
        };
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
      {17, 18, 22, 23}
    };
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





  template <int dim>
  void LaplaceProblem<dim>::run ()
  {
    for (unsigned int cycle=0; cycle<6; ++cycle)
      {
        std::cout << "Cycle " << cycle << ':' << std::endl;

        if (cycle == 0)
          create_coarse_grid ();

        setup_system ();

        std::cout << "   Number of active cells:       "
                  << triangulation.n_active_cells()
                  << std::endl
                  << "   Number of degrees of freedom: "
                  << dof_handler.n_dofs()
                  << std::endl
                  << "   Number of constraints       : "
                  << constraints.n_constraints()
                  << std::endl;

        assemble_system ();
        solve ();
        postprocess (cycle);
      }
  }

  template <int dim>
  std::pair<bool,unsigned int>
  predicate_ind(const TableIndices<dim> &ind);

  template<>
  std::pair<bool,unsigned int>
  predicate_ind<2>(const TableIndices<2> &ind)
  {
    const unsigned int v = ind[0]*ind[0]+ind[1]*ind[1];
    if (v>0 && v<7*7)
      return std::make_pair(true,v);
    else
      return std::make_pair(false,v);
  }

  template <int dim>
  void
  LaplaceProblem<dim>::
  estimate_smoothness (Vector<float> &smoothness_indicators)
  {
#ifdef OLD
    const unsigned int N = max_degree;

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
                k_vectors.push_back (Point<dim>(2.*numbers::PI * i,
                                                2.*numbers::PI * j));
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
                  k_vectors.push_back (Point<dim>(2.*numbers::PI * i,
                                                  2.*numbers::PI * j,
                                                  2.*numbers::PI * k));
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


    std::vector<Table<2,std::complex<double> > >
    fourier_transform_matrices (fe_collection.size());

    QGauss<1>      base_quadrature (2);
    QIterated<dim> quadrature (base_quadrature, N);


    for (unsigned int fe=0; fe<fe_collection.size(); ++fe)
      {
        fourier_transform_matrices[fe].reinit (n_fourier_modes,
                                               fe_collection[fe].dofs_per_cell);

        for (unsigned int k=0; k<n_fourier_modes; ++k)
          for (unsigned int j=0; j<fe_collection[fe].dofs_per_cell; ++j)
            {
              std::complex<double> sum = 0;
              for (unsigned int q=0; q<quadrature.size(); ++q)
                {
                  const Point<dim> x_q = quadrature.point(q);
                  sum += std::exp(std::complex<double>(0,1) *
                                  (k_vectors[k] * x_q)) *
                         fe_collection[fe].shape_value(j,x_q) *
                         quadrature.weight(q);
                }
              fourier_transform_matrices[fe](k,j)
                = sum;
            }
      }

    std::vector<std::complex<double> > fourier_coefficients (n_fourier_modes);
    Vector<double>                     local_dof_values;

    typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        local_dof_values.reinit (cell->get_fe().dofs_per_cell);
        cell->get_dof_values (solution, local_dof_values);

        for (unsigned int f=0; f<n_fourier_modes; ++f)
          {
            fourier_coefficients[f] = 0;

            for (unsigned int i=0; i<cell->get_fe().dofs_per_cell; ++i)
              fourier_coefficients[f] +=
                fourier_transform_matrices[cell->active_fe_index()](f,i)
                *
                local_dof_values(i);
          }

        std::map<unsigned int, double> k_to_max_U_map;
        for (unsigned int f=0; f<n_fourier_modes; ++f)
          if ((k_to_max_U_map.find (k_vectors_magnitude[f]) ==
               k_to_max_U_map.end())
              ||
              (k_to_max_U_map[k_vectors_magnitude[f]] <
               std::abs (fourier_coefficients[f])))
            k_to_max_U_map[k_vectors_magnitude[f]]
              = std::abs (fourier_coefficients[f]);

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
              sum_ln_U_ln_k += std::log (std::abs (fourier_coefficients[f])) *
                               ln_k[f];
            }

        const double mu
          = (1./(sum_1*sum_ln_k_square - sum_ln_k*sum_ln_k)
             *
             (sum_ln_k*sum_ln_U - sum_1*sum_ln_U_ln_k));

        smoothness_indicators(cell->active_cell_index()) = mu - 1.*dim/2;
      }
#else
    Vector<double>                     local_dof_values;

    typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        local_dof_values.reinit (cell->get_fe().dofs_per_cell);
        cell->get_dof_values (solution, local_dof_values);

        fourier->calculate(local_dof_values,
                           cell->active_fe_index(),
                           fourier_coefficients);

        std::pair<std::vector<unsigned int>, std::vector<double> > res =
          FESeries::process_coefficients<dim>(fourier_coefficients,
                                              predicate_ind<dim>,
                                              VectorTools::Linfty_norm);

        Assert (res.first.size() == res.second.size(),
                ExcInternalError());

        if (ln_k.size() == 0)
          {
            ln_k.resize(res.first.size(),0);
            for (unsigned int f = 0; f < ln_k.size(); f++)
              ln_k[f] = std::log (2.0*numbers::PI*std::sqrt(1.*res.first[f]));
          }

        for (unsigned int f = 0; f < res.second.size(); f++)
          res.second[f] = std::log(res.second[f]);

        std::pair<double,double> fit = FESeries::linear_regression(ln_k,res.second);
        smoothness_indicators(cell->active_cell_index()) = -fit.first - 1.*dim/2;
      }
#endif
  }
}



int main ()
{
  try
    {
      using namespace dealii;
      using namespace Step27;

      LaplaceProblem<2> laplace_problem;
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
