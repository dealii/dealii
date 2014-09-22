// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// calculation of the solution of laplace equation in the external
// domain the Boundary Element Method (implemented only for 2d
// problems).

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>

// all include files you need here

#include <deal.II/base/convergence_table.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/vector.h>
#include <deal.II/numerics/data_out.h>

#include <cmath>
#include <iostream>
#include <math.h>
#include <string>

std::ofstream logfile("output");

template<int spacedim>
class BEM
{

public:

  BEM();

  void run();

private:

  void make_dofs();
  void assemble_system();
  void solve();
  void output_results();

  void read_grid(std::string filename);
  void write_grid(std::string filename);


  Triangulation<spacedim-1,spacedim>    tria;
  FE_DGQ<spacedim-1, spacedim>          fe;
  DoFHandler<spacedim-1, spacedim>      dof_handler;

  // finite elements used to smoothen
  // the solution (from piecewise
  // constant to continuous piecewise
  // quadratic)
  FE_Q<spacedim-1, spacedim>            fe_q;
  DoFHandler<spacedim-1, spacedim>      dof_handler_q;

  FullMatrix<double>                    system_matrix;

  Vector<double>                        solution,
         smooth_solution,
         system_rhs,
         tangential_derivative,
         tangential_velocity,
         error;

  Point<spacedim>                       velocity;



};

template <int spacedim>
BEM<spacedim>::BEM()
  :
  fe(0),
  dof_handler(tria),
  fe_q(1),
  dof_handler_q(tria)
{
  velocity[0] = 1.;
}

template <int spacedim>
void
BEM<spacedim>::run()
{
  deallog.attach(logfile);
  deallog.depth_console(0);

  ConvergenceTable table;

  if (spacedim == 2)
    {


      read_grid(SOURCE_DIR "/grids/circle_R10.inp");

      Point<spacedim> p;
      HyperBallBoundary<spacedim-1, spacedim> boundary(p,10.);
      tria.set_boundary(1, boundary);

      // works up to cycle<9, but for testin purpose, we stop at 4
      for (unsigned int cycle=0; cycle<4; ++cycle)
        {

          tria.set_boundary(1, boundary);
          tria.refine_global(1);

          double side_length  =
            tria.begin_active()->vertex(0).norm()  // CIRCLE RADIUS
            *2*sin( numbers::PI / tria.n_active_cells() ),
            global_error = 0.;

          table.add_value("cycle", cycle);
          table.add_value("cells", tria.n_active_cells());

          make_dofs();
          assemble_system();
          solve();
          output_results();

          tria.set_boundary(1);


//  std::cout
//      << "solution"
//      << std::endl;
//  for(unsigned int i=0; i<dof_handler.n_dofs(); ++i)
//      std::cout
//    << i<< " - "
//    <<solution(i)<< std::endl;

//  std::cout
//      << "smooth solution"
//      << std::endl;
//  for(unsigned int i=0; i<dof_handler_q.n_dofs(); ++i)
//      std::cout
//    << i<< " - "
//    <<smooth_solution(i)<< std::endl;

//  std::cout
//      << "tangential velocity"
//      << std::endl;
//  for(unsigned int i=0; i<dof_handler_q.n_dofs(); ++i)
//      std::cout
//    << i<< " - "
//    <<tangential_derivative(i)<< std::endl;

//  std::cout<<tria.n_active_cells()<<std::endl;

//  std::cout
//      << "cell error" << std::endl
//      << "solution // smooth_sol // tang_deriv // tang_vel // exact_sol // error"
//      << std::endl;
//  for(unsigned int i=0; i<dof_handler_q.n_dofs(); ++i)
//      std::cout
//    << i<< " - "
//    << solution(i) << " // "
//    << smooth_solution(i)<< " // "
//    << tangential_derivative(i) << " // "
//    << tangential_velocity(i) << " // "
//    // exact solution
//    << 2. * velocity[0]
//       * sin (
//           numbers::PI/2.
//           - (1./2.+i) * 2.*numbers::PI/tria.n_active_cells()
//             )
//    << " // "
//    << error(i)
//    << std::endl;

//      std::cout
//    << "global error"
//    <<std::endl;

          for (unsigned int i=0; i<dof_handler_q.n_dofs(); ++i)
            global_error +=
              pow(fabs(error(i)),2) * side_length;
//      std::cout
//    << global_error
//    << std::endl;

          table.add_value("L^2 norm error", global_error);

        }
    }
  else
    Assert(false, ExcNotImplemented());

  table.set_scientific("L^2 norm error", true);
  table.evaluate_convergence_rates(
    "L^2 norm error",ConvergenceTable::reduction_rate_log2);
//     table.write_text(std::cout);
}

template <int spacedim>
void
BEM<spacedim>::make_dofs()
{
  dof_handler.distribute_dofs(fe);
  dof_handler_q.distribute_dofs(fe_q);

  system_matrix.reinit(dof_handler.n_dofs(), dof_handler.n_dofs());
  solution.reinit (dof_handler.n_dofs());
  system_rhs.reinit (dof_handler.n_dofs());
}

template <int spacedim>
void
BEM<spacedim>::assemble_system()
{
  const QGauss<spacedim-1>        quadrature_formula(2);
  const QGaussLog<1>              qlog(1);
  FEValues<spacedim-1, spacedim>
  fe_values_i (fe, quadrature_formula,
               update_JxW_values |
               update_cell_normal_vectors |
               update_quadrature_points ),
                                        fe_values_j (fe, quadrature_formula,
                                                     update_JxW_values |
                                                     update_cell_normal_vectors |
                                                     update_quadrature_points );

  const unsigned int        dofs_per_cell = fe.dofs_per_cell;
  const unsigned int        n_q_points    = quadrature_formula.size();

  // The following two matrices will
  // be subtracted to obtain the
  // system matrix.
  FullMatrix<double>        DLP_matrix(dof_handler.n_dofs()); // DLP stands for double layer potential
  FullMatrix<double>        mass_matrix(dof_handler.n_dofs());

  FullMatrix<double>        cell_mass_matrix (dofs_per_cell, dofs_per_cell);
  FullMatrix<double>        cell_DLP_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>            cell_rhs (dofs_per_cell);

  std::vector< Point<spacedim> > cell_normals_i, cell_normals_j;

  std::vector<types::global_dof_index> local_dof_indices_i (dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices_j (dofs_per_cell);


//     std::cout
//  << "no. of cells: "<< tria.n_active_cells()<< std::endl
//  << "no. of dofs: "<< dof_handler.n_dofs()<< std::endl
//  << "dofs per cell: "<< dofs_per_cell<< std::endl
//  << "side length (r=10): "<< 10 * 2*sin(numbers::PI/tria.n_cells())<<std::endl ;


  typename DoFHandler<spacedim-1, spacedim>::active_cell_iterator
  cell_i = dof_handler.begin_active(),
  cell_j = dof_handler.begin_active(),
  endc   = dof_handler.end();

  // cycle on i index
  for (; cell_i!=endc; ++cell_i)
    {

      fe_values_i.reinit (cell_i);
      cell_normals_i = fe_values_i.get_cell_normal_vectors();
      cell_i->get_dof_indices (local_dof_indices_i);

//  if (cell_i->index()%100==0)
//    std::cout
//        << cell_i->index()<< " || ";
//      << cell_i->vertex(0)<< " || "
//      << velocity*cell_normals_i[0] << "; "
//      << local_dof_indices_i[0]
//      << std::endl;

      cell_DLP_matrix = 0.;
      cell_mass_matrix = 0.;




      // assembling of the right hand side
      Point<spacedim-1>
      a_unit(0.),
             b_unit(1.);
      Point<spacedim>
      A = fe_values_i.get_mapping().
          transform_unit_to_real_cell(cell_i,
                                      a_unit),
          B = fe_values_i.get_mapping().
              transform_unit_to_real_cell(cell_i,
                                          b_unit);

      cell_rhs = 0.;


      // In order to obtain the Gram
      // determinant it is used JxW and
      // then it is divided by the weight.
      // Both the jacobian and the normals
      // are constant on the cell, so they
      // are taken in the first quadrature
      // point.
      double constant_factor =
        - pow(fe_values_i.JxW (0)/fe_values_i.get_quadrature().weight(0), 2)
        / numbers::PI
        * velocity*cell_normals_i[0];

      // These are constant on the cell so
      // there is no need to loop on
      // quadrature points. Besides, the
      // cell_rhs vector element index is
      // chosen = 0 because it is actually
      // a 1x1 matrix. There should be a
      // loop over cell dofs, but this is
      // necessary just for higher degree
      // elements.
      cell_rhs(0) +=
        constant_factor
        *(
          log( A.distance(B) ) + 8.*qlog.weight(0)
        );
      // A Gauss integration is performed to compute the
      // integrals on K_i X K_j in the case when i!=j.
      // cycle on j index
      for (cell_j=dof_handler.begin_active(); cell_j!=endc; ++cell_j)
        {
          fe_values_j.reinit (cell_j);
          cell_normals_j = fe_values_j.get_cell_normal_vectors();
          cell_j->get_dof_indices (local_dof_indices_j);

          if (cell_j != cell_i)
            {
              // The mass matrix has only diagonal
              // elements.
              mass_matrix(cell_i->index(), cell_j->index())= 0.;

              // with constant elements there is
              // only 1 dof per cell, so there is
              // no real cycle over cell dofs
              for (unsigned int a=0; a<dofs_per_cell; ++a)
                for (unsigned int q_point_i=0; q_point_i<n_q_points; q_point_i++)
                  for (unsigned int q_point_j=0; q_point_j<n_q_points; q_point_j++)
                    {
                      // If the integration is performed
                      // on two different elements there
                      // are no singularities in the
                      // domain of integration, so usual
                      // gauss formulas can be used.

                      for (unsigned int b=0; b<dofs_per_cell; ++b)
                        {
                          // assembling of double layer potential

                          cell_DLP_matrix(a,b) +=
                            1./numbers::PI
                            * contract(
                              cell_normals_i[q_point_i],
                              (fe_values_j.quadrature_point(q_point_j)
                               -
                               fe_values_i.quadrature_point(q_point_i))
                            )
                            / pow(
                              (fe_values_j.quadrature_point(q_point_j)
                               -
                               fe_values_i.quadrature_point(q_point_i))
                              .norm()
                              ,2
                            )
                            * fe_values_i.JxW(q_point_i)
                            * fe_values_j.JxW(q_point_j);

//            if ( (cell_i->index()==0) && (cell_j->index()==5) )
//          std::cout
//              << "("<< cell_i->index()<<","<<cell_j->index()<<"; "
//              <<q_point_i<<","<<q_point_j<<") "
//              << cell_DLP_matrix(a,b)<< " | "
//              << fe_values_i.quadrature_point(q_point_i)<< "||"
//              << fe_values_j.quadrature_point(q_point_j)<< "|<"
//              << contract(
//            cell_normals_i[q_point_i],
//            (fe_values_j.quadrature_point(q_point_j)
//             -
//             fe_values_i.quadrature_point(q_point_i))
//                   )<< "<<"
//              ;

                        }



                      // assembling of right hand side
                      cell_rhs(a) +=
                        - 1./numbers::PI
                        * velocity*cell_normals_j[q_point_j]
                        * log(
                          fe_values_i.quadrature_point(q_point_i)
                          .distance(
                            fe_values_j.quadrature_point(q_point_j)
                          )
                        )
                        * fe_values_i.JxW(q_point_i)
                        * fe_values_j.JxW(q_point_j);
//    std::cout
//        << cell_j->index()<< " -- "
//        << cell_rhs(i);
                    }
              for (unsigned int a=0; a<dofs_per_cell; ++a)
                for (unsigned int b=0; b<dofs_per_cell; ++b)
                  DLP_matrix(local_dof_indices_i[a],local_dof_indices_j[b])
                  +=
                    cell_DLP_matrix(a,b);
            }
          else      // case when cell_i=cell_j
            {
              // The mass matrix is simply a
              // diagonal matrix with the area of
              // each element as entries.
              for (unsigned q_point_i=0; q_point_i<n_q_points; ++q_point_i)
                mass_matrix(cell_i->index(), cell_i->index())
                +=
                  fe_values_i.JxW(q_point_i);

              // The double layer potential matrix
              // has no diagonal terms since the
              // scalar product between the cell
              // normal and the difference
              // <tt>x-y<\tt> is always zero. This
              // is because, if <tt>x<\tt> and
              // <tt>y<\tt> are on the same
              // element and the element is a
              // straight segment, their
              // difference is always orthogonal
              // to the cell normal.
              DLP_matrix(cell_i->index(), cell_i->index()) = 0.;
            }

          cell_DLP_matrix = 0.;

//    std::cout
//        << cell_DLP_matrix(0,0)<< " - ";

        }

//      std::cout
//    << cell_rhs(0)<<std::endl;

      for (unsigned int a=0; a<dofs_per_cell; ++a)
        system_rhs(local_dof_indices_i[a]) += cell_rhs(a);

      // end of assembling of the right hand side

    }

//     TableIndices<2> table_indices(dof_handler.n_dofs(), dof_handler.n_dofs());

  system_matrix.add(
    1, mass_matrix,
    -1, DLP_matrix
  );

//     std::cout
//  << "rhs"<< std::endl;
//     for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
//  std::cout
//      << system_rhs(i)<< " | "
//      << std::endl;

//     std::cout
//  << "DLP matrix"<< std::endl;
//     for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
//     {
//  for (unsigned int j=0; j<dof_handler.n_dofs(); ++j)
//      std::cout
//    << DLP_matrix(i,j)<< " | ";
//  std::cout
//      << std::endl;
//     }

//     std::cout
//  << "mass matrix"<< std::endl;
//     for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
//     {
//  for (unsigned int j=0; j<dof_handler.n_dofs(); ++j)
//      std::cout
//    << mass_matrix(i,j)<< " | ";
//  std::cout
//      << std::endl;
//     }

//     std::cout
//  << "system matrix"<< std::endl;
//     for (unsigned int i=0; i<dof_handler.n_dofs(); ++i)
//     {
//  for (unsigned int j=0; j<dof_handler.n_dofs(); ++j)
//      std::cout
//    << system_matrix(i,j)<< " | ";
//  std::cout
//      << std::endl;
//     }


}

template <int spacedim>
void
BEM<spacedim>::solve()
{

  SolverControl      solver_control(1000, 1e-12);
  SolverCG<>         cg(solver_control);

  cg.solve( system_matrix,
            solution,
            system_rhs,
            PreconditionIdentity() );

  // Smoothen the piecewise constant
  // solution to a continuous piecewise
  // quadratic function.
  smooth_solution.reinit(dof_handler_q.n_dofs());

  FETools::interpolate( dof_handler,
                        solution,
                        dof_handler_q,
                        smooth_solution);

  // Calculate the tangential derivative
  // of the smoothened potential, that
  // is the tangential component of the
  // perturbation induced in the
  // velocity.
  tangential_velocity.reinit(tria.n_active_cells());
  tangential_derivative.reinit(tria.n_active_cells());
  error.reinit(tria.n_active_cells());
  QMidpoint<spacedim-1>             q_midpoint;
  QTrapez<spacedim-1>               q_trapez;
  const QIterated<spacedim-1>       q_iterated(q_midpoint, 1);

  FEValues<spacedim-1, spacedim>    fe_values_q(fe_q,
                                                q_iterated,
                                                update_values |
                                                update_gradients |
                                                update_cell_normal_vectors );

  std::vector< Point<spacedim> > cell_normals(q_iterated.size());
  std::vector< Point<spacedim> > cell_tangentials(q_iterated.size());
  std::vector<types::global_dof_index> local_dof_indices (fe_q.dofs_per_cell);

  typename DoFHandler<spacedim-1, spacedim>::active_cell_iterator
  cell = dof_handler_q.begin_active(),
  endc = dof_handler_q.end();
  for (; cell!=endc; ++cell)
    {
      fe_values_q.reinit(cell);
      cell->get_dof_indices (local_dof_indices);

//QUADRATIC INTERPOLATION OF THE POTENTIAL
//  std::cout
//      << "indices - ";
//  for(unsigned int i=0; i<fe_q.dofs_per_cell; ++i)
//      std::cout
//    << local_dof_indices[(i==0? 0
//              :
//              i==1 ?  2
//              :
//              1
//             )]<< ", ";
//  std::cout<<std::endl;


      cell_normals = fe_values_q.get_cell_normal_vectors();
      for (unsigned int i=0; i<q_iterated.size(); ++i)
        {
          cell_tangentials[i][0] = cell_normals[i][1];
          cell_tangentials[i][1] = -cell_normals[i][0];
        }

      // Create a vector where gradients at
      // quadrature points are
      // stored. Notice that the first
      // factor (smooth_solution..) is taken
      // so that it is the coefficient of
      // the fun_th shape function of the cell.
      std::vector< Tensor<1,spacedim> > gradient(q_iterated.size());
      for (unsigned int pnt=0; pnt<q_iterated.size(); ++pnt)
        for (unsigned int fun=0; fun<fe_q.dofs_per_cell; ++fun)
          {
            gradient[pnt]+=
              smooth_solution(local_dof_indices[fun])
//QUADRATIC INTERPOLATION OF THE POTENTIAL
//    smooth_solution(local_dof_indices[(fun==0? 0
//               :
//               fun==1 ? 2
//               :
//               1
//             )]
//                   )
              *
              fe_values_q.shape_grad(fun, pnt);

//      std::cout
//    << smooth_solution(local_dof_indices[(fun==0? 0
//                  :
//                  fun==1 ? 2
//                  :
//                  1
//                )])<< " | "
//    << fe_values_q.shape_grad(fun, pnt)
//    << std::endl;
          }

//  std::cout
//      << "gradienti "<<std::endl;
//  for(unsigned int pnt=0; pnt<q_iterated.size(); ++pnt)
//      std::cout
//    << pnt << " - "
//    << gradient[pnt]
//    << std::endl;

//    std::cout
//      << "derivata tangenziale "
//      << cell->index()
//      << std::endl;


      for (unsigned int pnt=0; pnt<q_iterated.size(); ++pnt)
        {
          tangential_derivative(cell->index())=
            contract(
              gradient[pnt],
              cell_tangentials[pnt]
            )
            +
            contract(velocity,
                     cell_tangentials[0]
                    );
          error(cell->index())=
            tangential_derivative(cell->index())
            -
            2. * velocity[0]
            * sin (numbers::PI/2. - (1./2.+cell->index()) * 2.*numbers::PI/tria.n_active_cells() );
          tangential_velocity(cell->index())=
            contract(velocity,
                     cell_tangentials[0]
                    );

        }

// //   for(unsigned int pnt=0; pnt<q_iterated.size(); ++pnt)
//      std::cout
// //     << pnt <<" - "
//    << tangential_derivative(cell->index())
//    << std::endl;
    }


//     DataOut<spacedim-1, DoFHandler<spacedim-1,spacedim> > dataout;
//     dataout.attach_dof_handler(dof_handler_q);
//     dataout.add_data_vector(smooth_solution, "quadratic_potential");
//     dataout.build_patches(fe_q.degree);
//     char outname[50];
//     sprintf(outname, "bem_gradient.vtk");
//     std::ofstream file(outname);
//     dataout.write_vtk(file);
//     dataout.write_vtk(logfile);



}

template <int spacedim>
void
BEM<spacedim>::output_results()
{

  DataOut<spacedim-1, DoFHandler<spacedim-1,spacedim> > dataout;
  dataout.attach_dof_handler(dof_handler_q);
  dataout.add_data_vector(dof_handler_q, smooth_solution, "linear_potential");
  dataout.add_data_vector(tangential_derivative, "tangential_velocity", DataOut<spacedim-1, DoFHandler<spacedim-1,spacedim> >::type_cell_data);
  dataout.add_data_vector(error, "error", DataOut<spacedim-1, DoFHandler<spacedim-1,spacedim> >::type_cell_data);
  dataout.build_patches();
//    dataout.build_patches(fe_q.degree);
  dataout.write_vtk(logfile);

}

template <int spacedim>
void
BEM<spacedim>::read_grid(std::string filename)
{
  GridIn<spacedim-1, spacedim> gi;
  gi.attach_triangulation (tria);
  std::ifstream in (filename.c_str());
  gi.read_ucd (in);
}


template <int spacedim>
void
BEM<spacedim>::write_grid(std::string filename)
{
  GridOut grid_out;
  grid_out.set_flags (GridOutFlags::Ucd(true));
  grid_out.write_ucd (tria, logfile);
  grid_out.write_msh (tria, logfile);
}

int main ()
{

  BEM<2> bem;
  bem.run();

  return 0;
}

