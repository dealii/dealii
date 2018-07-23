/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2016 by the deal.II authors
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
 * ---------------------------------------------------------------------

 * Author: Ryan Grove, Clemson University
 *         Timo Heister, Clemson University
 */

// @sect3{Include files}

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/flow_function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_gmres.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/manifold_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/block_info.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_bdm.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/lac/sparse_direct.h>

#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/grid/grid_out.h>

// We need to include the following file to do timings:
#include <deal.II/base/timer.h>

// This includes the files necessary for us to use geometric Multigrid
#include <deal.II/multigrid/multigrid.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_matrix.h>

#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/loop.h>

#include <deal.II/integrators/laplace.h>
#include <deal.II/algorithms/any_data.h>

#include <deal.II/fe/fe_facet.h>

#include <fstream>
#include <sstream>

namespace Step59
{
  using namespace dealii;

  // In order to make it easy to switch between the different solvers that are
  // being used, we declare an enum that can be passed as an argument to the
  // constructor of the main class.
  struct SolverType
  {
    enum type {FGMRES_ILU, FGMRES_GMG, UMFPACK};
  };

  struct GeoType
  {
  enum type {Cube, Circle, L_shape, L_shape_diagonal};
  };


  // @sect3{Functions for Solution and Righthand side}
  //
  // The class Solution is used to define the boundary conditions and to
  // compute errors of the numerical solution. Note that we need to define the
  // values and gradients in order to compute L2 and H1 errors. Here we
  // decided to separate the implementations for 2d and 3d using template
  // specialization.
  //
  // Note that the first dim components are the velocity components
  // and the last is the pressure.
  template <int dim>
  class Solution : public Function<dim>
  {
  public:
    Solution () : Function<dim>(dim+1) {}
    virtual double value (const Point<dim> &p,
                          const unsigned int component = 0) const;
    virtual Tensor<1,dim> gradient (const Point<dim> &p,
                                    const unsigned int component = 0) const;
  };

  template <>
  double
  Solution<2>::value (const Point<2> &p,
                      const unsigned int component) const
  {
    Assert (component <= 2+1, ExcIndexRange(component,0,2+1));

    // smooth corner

//    const double p_int = - 2.0*exp(1.0)*cos(1.0)-2.0+2.0*cos(1.0)+2.0*exp(1.0);

//    if (component == 0)
//      return -exp(x)*(-y*cos(-y)+sin(-y));
//    if (component == 1)
//      return exp(x)*y*sin(-y);
//    if (component == 2)
//      return 2*exp(x)*sin(-y) - p_int;

    if (true)
    {
        using numbers::PI;
        const double x = p(0);
        const double y = p(1);
    // zero on BD's
    	if (component == 0)
    		return PI*sin(PI*x)*sin(PI*x)*sin(2.0*PI*y);
    	if (component == 1)
    		return -PI*sin(PI*y)*sin(PI*y)*sin(2.0*PI*x);
    	if (component == 2)
        return cos(PI*x)*sin(PI*y);
    }
    else
    {
    // l singular:
    	Functions::StokesLSingularity s;
    	return s.value(p, component);
    }
    return 0;
  }

  template <>
  double
  Solution<3>::value (const Point<3> &p,
                      const unsigned int component) const
  {
    Assert (component <= 3+1, ExcIndexRange(component,0,3+1));

    using numbers::PI;
    const double x = p(0);
    const double y = p(1);
    const double z = p(2);

    if (component == 0)
      return 2.*PI*sin(PI*x)*sin(PI*x)*sin(2.0*PI*y)*sin(2.0*PI*z);
    if (component == 1)
      return   -PI*sin(PI*y)*sin(PI*y)*sin(2.0*PI*x)*sin(2.0*PI*z);
    if (component == 2)
      return   -PI*sin(PI*z)*sin(PI*z)*sin(2.0*PI*x)*sin(2.0*PI*y);
    if (component == 3)
      return sin (PI * x) * cos (PI * y) * sin (PI * z);

    return 0;
  }

  // Note that for the gradient we need to return a Tensor<1,dim>
  template <>
  Tensor<1,2>
  Solution<2>::gradient (const Point<2> &p,
                         const unsigned int component) const
  {
    Assert (component <= 2, ExcIndexRange(component,0,2+1));

    if (true)
    {
    using numbers::PI;
    const double x = p(0);
    const double y = p(1);

    Tensor<1,2> return_value;
    if (component == 0)
      {
//        return_value[0] = PI * cos (PI * x);
//        return_value[1] = 0.0;

    	return_value[0] = PI*PI*sin(2.0*PI*y)*sin(2.0*PI*x);
    	return_value[1] = 2.0*PI*PI*sin(PI*x)*sin(PI*x)*cos(2.0*PI*y);

      }
    else if (component == 1)
      {
//        return_value[0] = y * PI * PI * sin( PI * x);
//        return_value[1] = - PI * cos (PI * x);

    	return_value[0] = -2.0*PI*PI*sin(PI*y)*sin(PI*y)*cos(2.0*PI*x);
    	return_value[1] = -PI*PI*sin(2.0*PI*y)*sin(2.0*PI*x);

      }
    else if (component == 2)
      {
        return_value[0] = PI * cos (PI * x) * cos (PI * y);
        return_value[1] =  - PI * sin (PI * x) * sin(PI * y);
      }

    	return return_value;
    }

    else
    {
        Functions::StokesLSingularity s;
        std::vector<Point<2> > points(1, p);
        std::vector<std::vector<Tensor<1,2> > > gradients(1);
        gradients[0].resize(2+1);
        s.vector_gradient_list (points, gradients);
        return gradients[0][component];
    }
  }

  template <>
  Tensor<1,3>
  Solution<3>::gradient (const Point<3> &p,
                         const unsigned int component) const
  {
    Assert (component <= 3, ExcIndexRange(component,0,3+1));

    using numbers::PI;
    const double x = p(0);
    const double y = p(1);
    const double z = p(2);

    Tensor<1,3> return_value;
    if (component == 0)
      {
        return_value[0] = 2*PI*PI*sin(2*PI*x)*sin(2*PI*y)*sin(2*PI*z);
        return_value[1] = 4*PI*PI*sin(PI*x)*sin(PI*x)*cos(2*PI*y)*sin(2*PI*z);
        return_value[2] = 4*PI*PI*sin(PI*x)*sin(PI*x)*cos(2*PI*z)*sin(2*PI*y);
      }
    else if (component == 1)
      {
        return_value[0] = -2*PI*PI*sin(PI*y)*sin(PI*y)*cos(2*PI*x)*sin(2*PI*z);
        return_value[1] = -PI*PI*sin(2*PI*x)*sin(2*PI*y)*sin(2*PI*z);
        return_value[2] = -2*PI*PI*sin(PI*y)*sin(PI*y)*cos(2*PI*z)*sin(2*PI*x);
      }
    else if (component == 2)
      {
        return_value[0] = -2*PI*PI*sin(PI*z)*sin(PI*z)*cos(2*PI*x)*sin(2*PI*y);
        return_value[1] = -2*PI*PI*sin(PI*z)*sin(PI*z)*cos(2*PI*y)*sin(2*PI*x);
        return_value[2] = -PI*PI*sin(2*PI*x)*sin(2*PI*y)*sin(2*PI*z);
      }
    else if (component == 3)
      {
        return_value[0] = PI * cos (PI * x) * cos (PI * y) * sin (PI * z);
        return_value[1] =  - PI * sin (PI * x) * sin(PI * y) * sin (PI * z);
        return_value[2] = PI * sin (PI * x) * cos (PI * y) * cos (PI * z);
      }

    return return_value;
  }

  template <int dim>
  class VelocitySolution : public Function<dim>
  {
  public:
    VelocitySolution () : Function<dim>(dim) {}
    virtual double value (const Point<dim> &p,
                          const unsigned int component = 0) const;
  };

  template <int dim>
  double
  VelocitySolution<dim>::value (const Point<dim> &p,
                      const unsigned int component) const
  {
	  Solution<dim> solution;
	  return solution.value(p, component);
  }




  // Implementation of $f$. See the introduction for more information.
  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    RightHandSide () : Function<dim>(dim+1) {}

    virtual double value (const Point<dim> &p,
                          const unsigned int component = 0) const;
  };

  template <>
  double
  RightHandSide<2>::value (const Point<2> &p,
                           const unsigned int component) const
  {
    Assert (component <= 2, ExcIndexRange(component,0,2+1));

    using numbers::PI;
    double x = p(0);
    double y = p(1);
    double nu = 1.0;

    if (true)
    {
    // RHS for 0 BD's
    if (component == 0)
    	return -nu*2.0*PI*PI*PI*(-2.0*sin(PI*x)*sin(PI*x)+cos(2.*PI*x))*sin(2.0*PI*y)-PI*sin(PI*x)*sin(PI*y);
    if (component == 1)
    	return nu*2.0*PI*PI*PI*(2.0*cos(2.0*PI*y)-1)*sin(2.0*PI*x)+PI*cos(PI*x)*cos(PI*y);
    if (component == 2)
    	return 0.0;
    }
    else
    {
    	if (component == 0)
    		return 0;
    	if (component == 1)
    		return 0;
    	if (component == 2)
    		return 0;
    }

    return 0;
  }

  template <>
  double
  RightHandSide<3>::value (const Point<3>   &p,
                           const unsigned int component) const
  {
    Assert (component <= 3, ExcIndexRange(component,0,3+1));

    using numbers::PI;
    double x = p(0);
    double y = p(1);
    double z = p(2);

//    if (component == 0)
//      return 2 * PI * PI * sin(PI * x) + PI * cos(PI * x) * cos(PI * y) * sin(PI * z);
//    if (component == 1)
//      return  - PI * PI * PI * y * cos (PI * x) + PI * (-1) * sin(PI * y)*sin(PI * x)*sin(PI * z);
//    if (component == 2)
//      return - PI * PI * PI * z * cos (PI * x) + PI * cos(PI * z)*sin(PI * x)*cos(PI * y);
//    if (component == 3)
//      return 0;
    if(component == 0)
    	return 4.*PI*PI*PI*(4.*sin(PI*x)*sin(PI*x)-cos(2.*PI*x))*sin(2.*PI*y)*sin(2.*PI*z)
    			+ PI * cos(PI * x) * cos(PI * y) * sin(PI * z);
    if(component == 1)
    	return -2.*PI*PI*PI*(4.*sin(PI*y)*sin(PI*y)-cos(2.*PI*y))*sin(2.*PI*x)*sin(2.*PI*z)
    			+ PI * (-1) * sin(PI * y)*sin(PI * x)*sin(PI * z);
    if(component == 2)
    	return -2.*PI*PI*PI*(4.*sin(PI*z)*sin(PI*z)-cos(2.*PI*z))*sin(2.*PI*x)*sin(2.*PI*y)
    			+ PI * cos(PI * z)*sin(PI * x)*cos(PI * y);
    if (component == 3)
        return 0;


    return 0;
  }

  template <int dim>
    class StokesMatrixIntegrator : public MeshWorker::LocalIntegrator<dim>
    {
    public:
      void cell(MeshWorker::DoFInfo<dim> &dinfo,
                typename MeshWorker::IntegrationInfo<dim> &info) const;
      void boundary(MeshWorker::DoFInfo<dim> &dinfo,
                    typename MeshWorker::IntegrationInfo<dim> &info) const;
      void face(MeshWorker::DoFInfo<dim> &dinfo1,
                MeshWorker::DoFInfo<dim> &dinfo2,
                typename MeshWorker::IntegrationInfo<dim> &info1,
                typename MeshWorker::IntegrationInfo<dim> &info2) const;
    };

    template <int dim>
    void StokesMatrixIntegrator<dim>::cell(
      MeshWorker::DoFInfo<dim> &dinfo,
      typename MeshWorker::IntegrationInfo<dim> &info) const
    {
  //    LocalIntegrators::Laplace::cell_matrix(dinfo.matrix(0,false).matrix, info.fe_values());

      const FEValuesBase<dim> &fe           = info.fe_values();
      const unsigned int      n_q_points    = fe.n_quadrature_points;
      const unsigned int      n_dofs        = fe.dofs_per_cell;
      const unsigned int      n_components  = fe.get_fe().n_components();
      const unsigned int      v_components  = n_components-1;
      const unsigned int      p_components  = n_components-1;

      RightHandSide<dim> right_hand_side;
//      ZeroFunction<dim> right_hand_side(dim+1);
      std::vector<Vector<double> >   rhs_values (n_q_points, Vector<double>(dim+1));
      right_hand_side.vector_value_list(fe.get_quadrature_points(), rhs_values);

      for (unsigned int k=0; k<n_q_points; ++k)
      {
      	for (unsigned int i=0; i<n_dofs; ++i)
      	{
      		// matrix
      		for (unsigned int j=0; j<n_dofs; ++j)
      		{
      			for (unsigned int d=0; d< v_components;++d)
      			{
      				dinfo.matrix(0,false).matrix(i, j) +=
      						(fe.shape_grad_component(j,k,d) * fe.shape_grad_component(i,k,d))*fe.JxW(k);
      				dinfo.matrix(0,false).matrix(i, j) +=
      						-(fe.shape_value_component(j,k,p_components)*fe.shape_grad_component(i,k,d)[d])*fe.JxW(k); // p*div_v
      				dinfo.matrix(0,false).matrix(i, j) +=
      				      	-(fe.shape_value_component(i,k,p_components)*fe.shape_grad_component(j,k,d)[d])*fe.JxW(k); // q*div_u
      			}
  				dinfo.matrix(0,false).matrix(i, j) +=
  					(fe.shape_value_component(j,k,p_components) * fe.shape_value_component(i,k,p_components))*fe.JxW(k); // p*q

      		}

      		// rhs
      		for (unsigned int d=0; d< v_components;++d)
      			 dinfo.vector(0).block(0)(i) += rhs_values[k](d)*fe.shape_value_component(i,k,d)*fe.JxW(k);
      	}
      }
    }


    template <int dim>
    void StokesMatrixIntegrator<dim>::boundary(
      MeshWorker::DoFInfo<dim> &dinfo,
      typename MeshWorker::IntegrationInfo<dim> &info) const
    {
//    	return;
      const FEValuesBase<dim> &fe           = info.fe_values();
      const unsigned int deg = fe.get_fe().tensor_degree();
     // std::cout << deg << std::endl;

      std::vector<Vector<double> >  boundary_values (fe.n_quadrature_points, Vector<double>(dim+1));
      Solution<dim> exact_solution;
      exact_solution.vector_value_list(fe.get_quadrature_points(), boundary_values);

      const double penalty_bd = 2. * deg * (deg+1) * dinfo.face->measure() / dinfo.cell->measure();
      double penalty = 4.0 * LocalIntegrators::Laplace::compute_penalty(dinfo, dinfo, deg, deg);
      const unsigned int      n_q_points    = fe.n_quadrature_points;
      const unsigned int      n_dofs        = fe.dofs_per_cell;
      const unsigned int      n_components  = fe.get_fe().n_components();
      const unsigned int      v_components  = n_components-1;
      const unsigned int      p_components  = n_components-1;

      for (unsigned int k=0; k<n_q_points; ++k)
      {
      	const Tensor<1,dim> n = fe.normal_vector(k);
      	for (unsigned int i=0; i<n_dofs; ++i)
      	{
      		for (unsigned int j=0; j<n_dofs; ++j)
      			for (unsigned int d=0; d<v_components; ++d)
      			{
      				dinfo.matrix(0,false).matrix(i,j) +=
				  (/*2. * */ fe.shape_value_component(i,k,d) * penalty * fe.shape_value_component(j,k,d)
      				           - (n * fe.shape_grad_component(i,k,d)) * fe.shape_value_component(j,k,d)
						    - (n * fe.shape_grad_component(j,k,d)) * fe.shape_value_component(i,k,d))*fe.JxW(k);

      		/*  b(p, v)_bd = {p}[v]n = p_j*u_i*n */
      				dinfo.matrix(0,false).matrix(i,j) +=
      						    fe.shape_value_component(j,k,p_components)*(fe.shape_value_component(i,k,d)*n[d])*fe.JxW(k);
      		/*  b(q, u)_bd = {q}[u]n = p_i*u_j*n */
      				dinfo.matrix(0,false).matrix(i,j) +=
      						    fe.shape_value_component(i,k,p_components)*(fe.shape_value_component(j,k,d)*n[d])*fe.JxW(k);

      			}

      		for (unsigned int d=0; d< v_components; ++d)
      		{
      		    dinfo.vector(0).block(0)(i) += ( fe.shape_value_component(i,k,d) * penalty * boundary_values[k][d]
												 -(n * fe.shape_grad_component(i,k,d)) * boundary_values[k][d])
												 * fe.JxW(k);

			    dinfo.vector(0).block(0)(i) += fe.shape_value_component(i,k,p_components)*boundary_values[k][d]*n[d]*fe.JxW(k);
      		}
      	}
      }
    }

    // Interior faces use the interior penalty method
    template <int dim>
    void StokesMatrixIntegrator<dim>::face(
      MeshWorker::DoFInfo<dim> &dinfo1,
      MeshWorker::DoFInfo<dim> &dinfo2,
      typename MeshWorker::IntegrationInfo<dim> &info1,
      typename MeshWorker::IntegrationInfo<dim> &info2) const
    {
//    	return;
      const unsigned int deg = info1.fe_values(0).get_fe().tensor_degree();
  //    LocalIntegrators::Laplace::ip_matrix(
  //      dinfo1.matrix(0,false).matrix, dinfo1.matrix(0,true).matrix,
  //      dinfo2.matrix(0,true).matrix, dinfo2.matrix(0,false).matrix,
  //      info1.fe_values(0), info2.fe_values(0),
  //      LocalIntegrators::Laplace::compute_penalty(dinfo1, dinfo2, deg, deg));

      double penalty = 4.0 * LocalIntegrators::Laplace::compute_penalty(dinfo1, dinfo2, deg, deg);
      const FEValuesBase<dim> &fe1          = info1.fe_values(0);
      const FEValuesBase<dim> &fe2          = info2.fe_values(0);
      const unsigned int      n_q_points    = fe1.n_quadrature_points;
      const unsigned int      n_dofs        = fe1.dofs_per_cell;
      const unsigned int      n_components  = fe1.get_fe().n_components();
      const unsigned int      v_components  = n_components-1;
      const unsigned int      p_components  = n_components-1;

      for (unsigned int k=0; k<n_q_points; ++k)
      {
      	const Tensor<1,dim> n = fe1.normal_vector(k);
      	for (unsigned int d=0; d<v_components; ++d)
      	{
      		for (unsigned int i=0; i<n_dofs; ++i)
      		{
      			for (unsigned int j=0; j<n_dofs; ++j)
      			{
      				const double u1_i = fe1.shape_value_component(i,k,d);
      				const double u1_j = fe1.shape_value_component(j,k,d);
      				const double dn_u1_i = n * fe1.shape_grad_component(i,k,d);
      				const double dn_u1_j = n * fe1.shape_grad_component(j,k,d);
      				const double u2_i = fe2.shape_value_component(i,k,d);
      				const double u2_j = fe2.shape_value_component(j,k,d);
      				const double dn_u2_i = n * fe2.shape_grad_component(i,k,d);
      				const double dn_u2_j = n * fe2.shape_grad_component(j,k,d);
      				const double p1_i    = fe1.shape_value_component(i,k,p_components);
      				const double p1_j    = fe1.shape_value_component(j,k,p_components);
      				const double p2_i    = fe2.shape_value_component(i,k,p_components);
      				const double p2_j    = fe2.shape_value_component(j,k,p_components);

              /*     gamma [U][V] = (U1-U2)(V1-V2) = U1V1 + U2V2 - U1V2 - U2V1
               *  ==> gamma*( u1_j*v1_i + u2_i*v2_j - u1_j*v2_i - u2_j*v1_i )  */
      				dinfo1.matrix(0,false).matrix(i, j) +=  penalty*u1_j*u1_i*fe1.JxW(k);
      				dinfo2.matrix(0,false).matrix(i, j) +=  penalty*u2_j*u2_i*fe1.JxW(k);
      				dinfo1.matrix(0,true).matrix(i, j) += -penalty*u1_i*u2_j*fe1.JxW(k);
      				dinfo2.matrix(0,true).matrix(i, j) += -penalty*u1_j*u2_i*fe1.JxW(k);

      	    /*    - {grad_U}[Vn] = -0.5(grad_U1 + grad_U2)n(V1-V2)
      	     *  ==>  0.5*(grad_n_u1_j*v2_i - grad_n_u1_j*v1_i + grad_n_u2_j*v2_i - grad_n_u2_j*v1_i    ) */
      				dinfo1.matrix(0,false).matrix(i, j) += -0.5*dn_u1_j*u1_i*fe1.JxW(k);
      				dinfo2.matrix(0,false).matrix(i, j) +=  0.5*dn_u2_j*u2_i*fe1.JxW(k);
      				dinfo1.matrix(0,true).matrix(i, j) += -0.5*dn_u2_j*u1_i*fe1.JxW(k);
      				dinfo2.matrix(0,true).matrix(i, j) +=  0.5*dn_u1_j*u2_i*fe1.JxW(k);

      	    /*    - [Un]{grad_V} = -0.5(U1-U2)(grad_V1 + grad_V2)
      	     *  ==>  0.5*(-u1_j*grad_v1_i +u2_j*grad_v2_i - u1_j*grad_v2_i + u2_j*grad_v1_i)*/
      				dinfo1.matrix(0,false).matrix(i, j) += -0.5*u1_j*dn_u1_i*fe1.JxW(k);
      				dinfo2.matrix(0,false).matrix(i, j) +=  0.5*u2_j*dn_u2_i*fe1.JxW(k);
      				dinfo1.matrix(0,true).matrix(i, j) +=  0.5*u2_j*dn_u1_i*fe1.JxW(k);
      				dinfo2.matrix(0,true).matrix(i, j) += -0.5*u1_j*dn_u2_i*fe1.JxW(k);

      		/*     b(q, u)_face = {q}[u]n = 0.5*(q1+q2)*[u1-u2]n
      		 *  ==>  0.5*(q1_i*u1_j - q2_i*u2_j + q2_i*u1_j -q1_i*u2_j)n   */
      				dinfo1.matrix(0,false).matrix(i, j) +=  0.5*p1_i*n[d]*u1_j*fe1.JxW(k);
      				dinfo2.matrix(0,false).matrix(i, j) += -0.5*p2_i*n[d]*u2_j*fe1.JxW(k);
      				dinfo1.matrix(0,true).matrix(i, j)  += -0.5*p1_i*n[d]*u2_j*fe1.JxW(k);
      				dinfo2.matrix(0,true).matrix(i, j)  +=  0.5*p2_i*n[d]*u1_j*fe1.JxW(k);

      	    /* 	   b(p, v)_face = {p}[v]n = 0.5*(p1+p2)*[v1-v2]n
      	     *  ==>  0.5*(p1_j*v1_i - p2_j*v2_i + p2_j*v1_i -p1_j*v2_i )n	*/
      				dinfo1.matrix(0,false).matrix(i, j) +=  0.5*p1_j*n[d]*u1_i*fe1.JxW(k);
      				dinfo2.matrix(0,false).matrix(i, j) += -0.5*p2_j*n[d]*u2_i*fe1.JxW(k);
      				dinfo1.matrix(0,true).matrix(i, j)  +=  0.5*p2_j*n[d]*u1_i*fe1.JxW(k);
      				dinfo2.matrix(0,true).matrix(i, j)  += -0.5*p1_j*n[d]*u2_i*fe1.JxW(k);
      			}
      		}
      	}
      }
    }

    template <int dim>
      class RHSIntegrator : public MeshWorker::LocalIntegrator<dim>
      {
      public:
        void cell(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const;
        void boundary(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const;
        void face(MeshWorker::DoFInfo<dim> &dinfo1,
                  MeshWorker::DoFInfo<dim> &dinfo2,
                  typename MeshWorker::IntegrationInfo<dim> &info1,
                  typename MeshWorker::IntegrationInfo<dim> &info2) const;
      };


      template <int dim>
      void RHSIntegrator<dim>::cell(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const
      {
          const FEValuesBase<dim> &fe = info.fe_values();
          const unsigned int      n_q_points    = fe.n_quadrature_points;
          const unsigned int      n_dofs        = fe.dofs_per_cell;
          const unsigned int      n_components  = fe.get_fe().n_components();

          RightHandSide<dim> right_hand_side;
          std::vector<Vector<double> >   rhs_values (n_q_points, Vector<double>(dim+1));

          right_hand_side.vector_value_list(fe.get_quadrature_points(), rhs_values);

          for (unsigned int k=0; k<n_q_points;++k)
        	  for(unsigned int d=0; d<n_components-1; ++d)
        		  for(unsigned int i=0;i<n_dofs; ++i)
        			  dinfo.vector(0).block(0)(i) +=  rhs_values[k](d)*fe.shape_value_component(i,k,d)*fe.JxW(k);
      }

      template <int dim>
      void RHSIntegrator<dim>::boundary(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const
      {
//    	  return;
        const FEValuesBase<dim> &fe = info.fe_values();
        const unsigned int      n_components  = fe.get_fe().n_components();
        const unsigned int      v_components  = n_components-1;
        const unsigned int      p_components  = n_components-1;

        std::vector<Vector<double>> boundary_values(fe.n_quadrature_points);
        Solution<dim> exact_solution;
        exact_solution.value_list(fe.get_quadrature_points(), boundary_values);

        const unsigned int deg = fe.get_fe().tensor_degree();

        const double penalty = 2. * deg * (deg+1) * dinfo.face->measure() / dinfo.cell->measure();

        for (unsigned k=0; k<fe.n_quadrature_points; ++k)
        {
          const Tensor<1,dim> n = fe.normal_vector(k);
          for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        	  for (unsigned int d = 0; d<v_components; ++d)
        	  {
        		  dinfo.vector(0).block(0)(i) += ( fe.shape_value_component(i,k,d) * penalty * boundary_values[k][d]
												 -(n * fe.shape_grad_component(i,k,d)) * boundary_values[k][d])
												 * fe.JxW(k);

				  dinfo.vector(0).block(0)(i) += fe.shape_value_component(i,k,p_components)*boundary_values[k][d]*n[d]*fe.JxW(k);
        	  }
        }
      }


      template <int dim>
      void RHSIntegrator<dim>::face(MeshWorker::DoFInfo<dim> &,
                                    MeshWorker::DoFInfo<dim> &,
                                    typename MeshWorker::IntegrationInfo<dim> &,
                                    typename MeshWorker::IntegrationInfo<dim> &) const
      {}


      // The third local integrator is responsible for the contributions to the
      // error estimate. This is the standard energy estimator due to Karakashian
      // and Pascal (2003).

  // @sect3{ASPECT BlockSchurPreconditioner}

  // In the following, we will implement a preconditioner that expands
  // on the ideas discussed in the Results section of step-22.
  // Specifically, we
  // 1. use an upper block-triangular preconditioner because we want to
  // use right preconditioning.
  // 2. optionally allow using an inner solver for the velocity block instead
  // of a single preconditioner application.
  // 3. do not use InverseMatrix but explicitly call SolverCG.
  // This approach is also used in the ASPECT code
  // (see http://aspect.dealii.org) that solves the Stokes equations in
  // the context of simulating convection in the earth mantle, and which
  // has been used to solve problems on many thousands of processors.
  //
  // The bool flag @p do_solve_A in the constructor allows us to either
  // apply the preconditioner for the velocity block once or use an inner
  // iterative solver for a more accurate approximation instead.
  //
  // Notice how we keep track of the sum of the inner iterations
  // (preconditioner applications).
  template <class PreconditionerAType, class PreconditionerSType>
  class BlockSchurPreconditioner : public Subscriptor
  {
  public:
    BlockSchurPreconditioner (const BlockSparseMatrix<double>  &system_matrix,
                              const SparseMatrix<double> &schur_complement_matrix,
                              const PreconditionerAType &preconditioner_A,
                              const PreconditionerSType &preconditioner_S,
                              const bool do_solve_A);

    void vmult (BlockVector<double>       &dst,
                const BlockVector<double> &src) const;

    mutable unsigned int n_iterations_A;
    mutable unsigned int n_iterations_S;

  private:
    const BlockSparseMatrix<double> &system_matrix;
    const SparseMatrix<double> &schur_complement_matrix;
    const PreconditionerAType &preconditioner_A;
    const PreconditionerSType &preconditioner_S;

    const bool do_solve_A;
  };

  template <class PreconditionerAType, class PreconditionerSType>
  BlockSchurPreconditioner<PreconditionerAType, PreconditionerSType>::
  BlockSchurPreconditioner (const BlockSparseMatrix<double>  &system_matrix,
                            const SparseMatrix<double> &schur_complement_matrix,
                            const PreconditionerAType &preconditioner_A,
                            const PreconditionerSType &preconditioner_S,
                            const bool do_solve_A)
    :
    n_iterations_A (0),
    n_iterations_S (0),
    system_matrix (system_matrix),
    schur_complement_matrix (schur_complement_matrix),
    preconditioner_A (preconditioner_A),
    preconditioner_S (preconditioner_S),
    do_solve_A (do_solve_A)
  {}



  template <class PreconditionerAType, class PreconditionerSType>
  void
  BlockSchurPreconditioner<PreconditionerAType, PreconditionerSType>::
  vmult (BlockVector<double>       &dst,
         const BlockVector<double> &src) const
  {
    Vector<double> utmp(src.block(0));

    // First solve with the approximation for S
    {
      SolverControl solver_control(1000, 1e-6 * src.block(1).l2_norm());
      SolverCG<>    cg (solver_control);

      dst.block(1) = 0.0;
      cg.solve(schur_complement_matrix,
               dst.block(1),
               src.block(1),
               preconditioner_S);

      n_iterations_S += solver_control.last_step();
      dst.block(1) *= -1.0;
    }

    // Second, apply the top right block (B^T)
    {
      system_matrix.block(0,1).vmult(utmp, dst.block(1));
      utmp *= -1.0;
      utmp += src.block(0);
    }

    // Finally, either solve with the top left block
    // or just apply one preconditioner sweep
    if (do_solve_A == true)
      {
        SolverControl solver_control(10000, utmp.l2_norm()*1e-4);
        SolverCG<>    cg (solver_control);

        dst.block(0) = 0.0;
        cg.solve(system_matrix.block(0,0),
                 dst.block(0),
                 utmp,
                 preconditioner_A);

        n_iterations_A += solver_control.last_step();
      }
    else
      {
        preconditioner_A.vmult (dst.block(0), utmp);
        n_iterations_A += 1;
      }
  }

  template <int dim>
    class Estimator : public MeshWorker::LocalIntegrator<dim>
    {
    public:
 	 Estimator(double viscosity) : MeshWorker::LocalIntegrator<dim>(), viscosity (viscosity) {}
      void cell(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const;
      void boundary(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const;
      void face(MeshWorker::DoFInfo<dim> &dinfo1,
                MeshWorker::DoFInfo<dim> &dinfo2,
                typename MeshWorker::IntegrationInfo<dim> &info1,
                typename MeshWorker::IntegrationInfo<dim> &info2) const;
    private:
      double viscosity;
    };

   // The cell contribution is the Laplacian of the discrete solution, since
   // the right hand side is zero.
   template <int dim>
   void Estimator<dim>::cell(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const
   {
     // compute h \| f + \triangle u -  \|_0

     const FEValuesBase<dim> &fe           = info.fe_values();
     const unsigned int      n_q_points    = fe.n_quadrature_points;
//     const unsigned int      dofs_per_cell = fe.dofs_per_cell;

     //  Check problem type
    const RightHandSide<dim>          right_hand_side;
 	std::vector<Vector<double> >      rhs_values (n_q_points, Vector<double>(dim+1));

     right_hand_side.vector_value_list(fe.get_quadrature_points(), rhs_values);


     const double h = dinfo.cell->diameter();

     const std::vector<Tensor<2,dim> > &DDuh1 = info.hessians[0][0];
     const std::vector<Tensor<2,dim> > &DDuh2 = info.hessians[0][1];

     const std::vector<Tensor<1,dim> > &Duh1  = info.gradients[0][0];
     const std::vector<Tensor<1,dim> > &Duh2  = info.gradients[0][1];

     const std::vector<double> &uh1 = info.values[0][0];
     const std::vector<double> &uh2 = info.values[0][1];

     const std::vector<Tensor<1,dim> > &Dph   = info.gradients[0][2];

 //    std::cout << " Dof_per_cell: " << fe.dofs_per_cell << std::endl;
 //    std::cout << " Quadrature_n: " << fe.n_quadrature_points << std::endl;
 //    std::cout << " Size: " << info.values[0].size() << std::endl;

     double R1 = 0.0;
     double R2 = 0.0;
     for (unsigned k=0; k<fe.n_quadrature_points; ++k)
       {
         Vector<double> R1_q(dim);
         double         R2_q;

         Tensor<1, dim> uh;
         uh[0] = uh1[k];
         uh[1] = uh2[k];

//         R1_q(0) = rhs_values[k][0] + viscosity*trace(DDuh1[k]) -(Duh1[k]*uh) - Dph[k][0];
//         R1_q(1) = rhs_values[k][1] + viscosity*trace(DDuh2[k]) -(Duh2[k]*uh) - Dph[k][1];

         R1_q(0) = rhs_values[k][0] + viscosity*trace(DDuh1[k]) - Dph[k][0];
         R1_q(1) = rhs_values[k][1] + viscosity*trace(DDuh2[k]) - Dph[k][1];

 		R2_q = Duh1[k][0] + Duh2[k][1];

 		R1 += R1_q * R1_q * fe.JxW(k);
        R2 += R2_q * R2_q * fe.JxW(k);
       }

     R1 = h*std::sqrt(R1); // not h*h
     R2 = std::sqrt(R2);

     dinfo.value(0) += R1;
     dinfo.value(0) += R2;
   }

   // At the boundary, we use simply a weighted form of the boundary residual,
   // namely the norm of the difference between the finite element solution and
   // the correct boundary condition.
   template <int dim>
   void Estimator<dim>::boundary(MeshWorker::DoFInfo<dim> &, typename MeshWorker::IntegrationInfo<dim> &) const
   {}


   // Finally, on interior faces, the estimator consists of the jumps of the
   // solution and its normal derivative, weighted appropriately.
   template <int dim>
   void Estimator<dim>::face(MeshWorker::DoFInfo<dim> &dinfo1,
                             MeshWorker::DoFInfo<dim> &dinfo2,
                             typename MeshWorker::IntegrationInfo<dim> &info1,
                             typename MeshWorker::IntegrationInfo<dim> &info2) const
   {

 //    // 1/2 * sqrt(h_e) * \| [n . \nabla u] \|_0

     const FEValuesBase<dim> &fe = info1.fe_values();
     const std::vector<double> &ph1 = info1.values[0][2];
     const std::vector<double> &ph2 = info2.values[0][2];

 //    std::cout << " Quadrature_n: " << fe.n_quadrature_points << std::endl;
 //    std::cout << " Size= " << info1.values[0][2].size() << std::endl;

     const std::vector<Tensor<1,dim> > &Du1_h1 = info1.gradients[0][0];
     const std::vector<Tensor<1,dim> > &Du2_h1 = info1.gradients[0][1];

     const std::vector<Tensor<1,dim> > &Du1_h2 = info2.gradients[0][0];
     const std::vector<Tensor<1,dim> > &Du2_h2 = info2.gradients[0][1];

     const double h = dinfo1.face->measure();

     for (unsigned k=0; k<fe.n_quadrature_points; ++k)
       {
     	Vector<double> diff(dim);

     	diff(0) =  fe.normal_vector(k)*Du1_h1[k] - fe.normal_vector(k)[0]*ph1[k]
 				  -fe.normal_vector(k)*Du1_h2[k] + fe.normal_vector(k)[0]*ph2[k];

       	diff(1) =  fe.normal_vector(k)*Du2_h1[k] - fe.normal_vector(k)[1]*ph1[k]
 				  -fe.normal_vector(k)*Du2_h2[k] + fe.normal_vector(k)[1]*ph2[k];

     	  dinfo1.value(0) += h*(diff*diff) * fe.JxW(k);
       }

     dinfo1.value(0) = 0.5*std::sqrt(dinfo1.value(0));
     dinfo2.value(0) = dinfo1.value(0);
   }

   template <int dim>
     class ErrorIntegrator : public MeshWorker::LocalIntegrator<dim>
     {
     public:
       void cell(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const;
       void boundary(MeshWorker::DoFInfo<dim> &dinfo, typename MeshWorker::IntegrationInfo<dim> &info) const;
       void face(MeshWorker::DoFInfo<dim> &dinfo1,
                 MeshWorker::DoFInfo<dim> &dinfo2,
                 typename MeshWorker::IntegrationInfo<dim> &info1,
                 typename MeshWorker::IntegrationInfo<dim> &info2) const;
     };

     // Here we have the integration on cells. There is currently no good
     // interface in MeshWorker that would allow us to access values of regular
     // functions in the quadrature points. Thus, we have to create the vectors
     // for the exact function's values and gradients inside the cell
     // integrator. After that, everything is as before and we just add up the
     // squares of the differences.

     // Additionally to computing the error in the energy norm, we use the
     // capability of the mesh worker to compute two functionals at the same time
     // and compute the <i>L<sup>2</sup></i>-error in the same loop. Obviously,
     // this one does not have any jump terms and only appears in the integration
     // on cells.
     template <int dim>
     void ErrorIntegrator<dim>::cell(
       MeshWorker::DoFInfo<dim> &dinfo,
       typename MeshWorker::IntegrationInfo<dim> &info) const
     {
       const FEValuesBase<dim> &fe = info.fe_values();

       const std::vector<Tensor<1,dim> > &Du1h = info.gradients[0][0];
       const std::vector<Tensor<1,dim> > &Du2h = info.gradients[0][1];

       for (unsigned k=0; k<fe.n_quadrature_points; ++k)
         {
           dinfo.value(0) += (Du1h[k][0]+Du2h[k][1])*(Du1h[k][0]+Du2h[k][1]) * fe.JxW(k);
         }
       dinfo.value(0) = std::sqrt(dinfo.value(0));
     }


     template <int dim>
     void ErrorIntegrator<dim>::boundary(
       MeshWorker::DoFInfo<dim> &,
       typename MeshWorker::IntegrationInfo<dim> &) const
     {
   //    const FEValuesBase<dim> &fe = info.fe_values();
   //
   //    std::vector<double> exact_values(fe.n_quadrature_points);
   //    exact_solution.value_list(fe.get_quadrature_points(), exact_values);
   //
   //    const std::vector<double> &uh = info.values[0][0];
   //
   //    const unsigned int deg = fe.get_fe().tensor_degree();
   //    const double penalty = 2. * deg * (deg+1) * dinfo.face->measure() / dinfo.cell->measure();
   //
   //    for (unsigned k=0; k<fe.n_quadrature_points; ++k)
   //      {
   //        const double diff = exact_values[k] - uh[k];
   //        dinfo.value(0) += penalty * diff * diff * fe.JxW(k);
   //      }
   //    dinfo.value(0) = std::sqrt(dinfo.value(0));
     }


     template <int dim>
     void ErrorIntegrator<dim>::face(
       MeshWorker::DoFInfo<dim> &,
       MeshWorker::DoFInfo<dim> &,
       typename MeshWorker::IntegrationInfo<dim> &,
       typename MeshWorker::IntegrationInfo<dim> &) const
     {
   //    const FEValuesBase<dim> &fe = info1.fe_values();
   //    const std::vector<double> &uh1 = info1.values[0][0];
   //    const std::vector<double> &uh2 = info2.values[0][0];
   //
   //    const unsigned int deg = fe.get_fe().tensor_degree();
   //    const double penalty1 = deg * (deg+1) * dinfo1.face->measure() / dinfo1.cell->measure();
   //    const double penalty2 = deg * (deg+1) * dinfo2.face->measure() / dinfo2.cell->measure();
   //    const double penalty = penalty1 + penalty2;
   //
   //    for (unsigned k=0; k<fe.n_quadrature_points; ++k)
   //      {
   //        double diff = uh1[k] - uh2[k];
   //        dinfo1.value(0) += (penalty * diff*diff)
   //                           * fe.JxW(k);
   //      }
   //    dinfo1.value(0) = std::sqrt(dinfo1.value(0));
   //    dinfo2.value(0) = dinfo1.value(0);
     }

  // @sect3{The StokesProblem class}
  //
  // This is the main class of the problem.
  template <int dim>
  class StokesProblem
  {
  public:
    StokesProblem (FiniteElement<dim> &fe,
    		       const unsigned int pressure_degree,
                   SolverType::type solver_type,
				   GeoType::type geo_type);
    void run ();

  private:
    void setup_dofs ();
    void assemble_system ();
    void assemble_matrix_rhs_DG();
    void assemble_system_DG ();
    void assemble_system_mesh_loop ();
    void assemble_rhs_DG ();
    void assemble_multigrid ();
    void solve ();
    void compute_errors (unsigned int k);
    void output_results (const unsigned int refinement_cycle) const;

    const unsigned int            pressure_degree;
    SolverType::type              solver_type;
    GeoType::type                 geo_type;

    Triangulation<dim>            triangulation;
    FiniteElement<dim>            &fe;
    DoFHandler<dim>               dof_handler;
    DoFHandler<dim>               velocity_dof_handler;

    ConstraintMatrix              constraints;

    BlockSparsityPattern          sparsity_pattern;
    BlockSparseMatrix<double>     system_matrix;
    SparseMatrix<double>          pressure_mass_matrix;

    BlockVector<double>           solution;
    BlockVector<double>           system_rhs;

    MGLevelObject<SparsityPattern>        mg_sparsity_patterns;
    MGLevelObject<SparseMatrix<double> >  mg_matrices;
    MGLevelObject<SparseMatrix<double> >  mg_interface_matrices;
    MGConstrainedDoFs                     mg_constrained_dofs;

    TimerOutput                           computing_timer;

    double                        last_l2_error;
    double                        last_H1_error;
    double                        last_Hdiv_error1;
    double                        last_Hdiv_error2;

    Vector<double>    global_div_diff;
    Vector<double>    cell_div_diff;
    Vector<double>    global_l2_diff;
    Vector<double>    global_h1_diff;
  };



  template <int dim>
  StokesProblem<dim>::StokesProblem (FiniteElement<dim> &fe,
		     	 	 	 	 	 	 const unsigned int pressure_degree,
                                     SolverType::type solver_type,
									 GeoType::type geo_type)
    :
    pressure_degree (pressure_degree),
    solver_type (solver_type),
	geo_type(geo_type),
//    triangulation (Triangulation<dim>::limit_level_difference_at_vertices),
    triangulation (Triangulation<dim>::maximum_smoothing),
    // Finite element for the velocity only:
//    velocity_fe (FE_RaviartThomas<dim>(pressure_degree), 1),
////    velocity_fe (FE_BDM<dim>(pressure_degree+1), 1),
//    // Finite element for the whole system:
//    fe (velocity_fe, 1, FE_DGQ<dim> (pressure_degree), 1),
	fe(fe),
    dof_handler (triangulation),
    velocity_dof_handler (triangulation),
    computing_timer (std::cout, TimerOutput::never,
                     TimerOutput::wall_times)
  {}

// @sect4{StokesProblem::setup_dofs}

// This function sets up the DoFHandler, matrices, vectors, and Multigrid
// structures (if needed).
  template <int dim>
  void StokesProblem<dim>::setup_dofs ()
  {
    TimerOutput::Scope scope(computing_timer, "Setup");

    system_matrix.clear ();
    pressure_mass_matrix.clear ();

    // The main DoFHandler only needs active DoFs, so we are not calling
    // distribute_mg_dofs() here
    dof_handler.distribute_dofs(fe);

    // This block structure separates the dim velocity components from
    // the pressure component (used for reordering). Note that we have
    // 2 instead of dim+1 blocks like in step-22, because our FESystem
    // is nested and the dim velocity components appear as one block.
    std::vector<unsigned int> block_component (2);
    block_component[0] = 0;
    block_component[1] = 1;

    // Velocities start at component 0:
    const FEValuesExtractors::Vector velocities(0);

    // ILU behaves better if we apply a reordering to reduce fillin. There
    // is no advantage in doing this for the other solvers.
    if (solver_type == SolverType::FGMRES_ILU)
      {
        TimerOutput::Scope ilu_specific(computing_timer, "(ILU specific)");
        DoFRenumbering::Cuthill_McKee (dof_handler);
      }

    // This ensures that all velocities DoFs are enumerated before the
    // pressure unknowns. This allows us to use blocks for vectors and
    // matrices and allows us to get the same DoF numbering for
    // dof_handler and velocity_dof_handler.
    DoFRenumbering::block_wise (dof_handler);

    std::vector<types::global_dof_index> dofs_per_block (2);
    DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);
    const unsigned int n_u = dofs_per_block[0],
                       n_p = dofs_per_block[1];

    std::cout << "\tNumber of active cells: "
              << triangulation.n_active_cells()
              << std::endl
              << "\tNumber of degrees of freedom: "
              << dof_handler.n_dofs()
              << " (" << n_u << '+' << n_p << ')'
              << std::endl;

    constraints.reinit ();
    if (false)
    {
        FEValuesExtractors::Vector velocities(0);
        DoFTools::make_hanging_node_constraints (dof_handler, constraints);
        std::set<types::boundary_id> no_normal_flux_boundaries;
        no_normal_flux_boundaries.insert (0);
        VectorTools::compute_no_normal_flux_constraints(dof_handler, 0, no_normal_flux_boundaries, constraints);
    }

    if (fe.base_element(0).conforms(FiniteElementData<dim>::H1))
    {
    FEValuesExtractors::Vector velocities(0);
    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    VectorTools::interpolate_boundary_values (dof_handler,
                                              0,
                                              Solution<dim>(),
                                              constraints,
                                              fe.component_mask(velocities));
    }
    else if (fe.base_element(0).conforms(FiniteElementData<dim>::Hdiv))
    {
        FEValuesExtractors::Vector velocities(0);
        DoFTools::make_hanging_node_constraints (dof_handler, constraints);
      VectorTools::project_boundary_values_div_conforming(dof_handler, 0, VelocitySolution<dim>(), 0, constraints);
    }

    constraints.close ();

    {
      BlockDynamicSparsityPattern csp (dofs_per_block, dofs_per_block);
//      DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);
      DoFTools::make_flux_sparsity_pattern(dof_handler, csp, constraints);
      sparsity_pattern.copy_from (csp);
    }
    system_matrix.reinit (sparsity_pattern);

    solution.reinit (dofs_per_block);
    system_rhs.reinit (dofs_per_block);
  }


// @sect4{StokesProblem::assemble_system}

// In this function, the system matrix is assembled. We assemble the pressure
// mass matrix in the (1,1) block (if needed) and move it out of this location
// at the end of this function.
  template <int dim>
  void StokesProblem<dim>::assemble_system ()
  {
    TimerOutput::Scope assemble(computing_timer, "Assemble");
    system_matrix=0;
    system_rhs=0;

    // If true, we will assemble the pressure mass matrix in the (1,1) block:
    const bool assemble_pressure_mass_matrix = (solver_type == SolverType::UMFPACK) ? false : true;

    QGauss<dim>   quadrature_formula(pressure_degree+2);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    |
                             update_quadrature_points  |
                             update_JxW_values |
                             update_gradients);

    const unsigned int   dofs_per_cell   = fe.dofs_per_cell;

    const unsigned int   n_q_points      = quadrature_formula.size();

    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       local_rhs (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    const RightHandSide<dim>          right_hand_side;
    std::vector<Vector<double> >      rhs_values (n_q_points,
                                                  Vector<double>(dim+1));

    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);

    std::vector<SymmetricTensor<2,dim> > symgrad_phi_u (dofs_per_cell);
    std::vector<double>                  div_phi_u   (dofs_per_cell);
    std::vector<double>                  phi_p       (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        fe_values.reinit (cell);
        local_matrix = 0;
        local_rhs = 0;

        right_hand_side.vector_value_list(fe_values.get_quadrature_points(),
                                          rhs_values);

        for (unsigned int q=0; q<n_q_points; ++q)
          {
            for (unsigned int k=0; k<dofs_per_cell; ++k)
              {
                symgrad_phi_u[k] = fe_values[velocities].symmetric_gradient (k, q);
                div_phi_u[k]     = fe_values[velocities].divergence (k, q);
                phi_p[k]         = fe_values[pressure].value (k, q);
              }

            for (unsigned int i=0; i<dofs_per_cell; ++i)
              {
                for (unsigned int j=0; j<=i; ++j)
                  {
                    local_matrix(i,j) += (//2 * (symgrad_phi_u[i] * symgrad_phi_u[j])
                                          (scalar_product( fe_values[velocities].gradient (i,q)
                                           , fe_values[velocities].gradient (j,q)
                                                          )
                                           )
                                          - div_phi_u[i] * phi_p[j]
                                          - phi_p[i] * div_phi_u[j]
                                          + (assemble_pressure_mass_matrix ? phi_p[i] * phi_p[j] : 0))
                                         * fe_values.JxW(q);

                  }

                const unsigned int component_i =
                  fe.system_to_component_index(i).first;
                local_rhs(i) += fe_values.shape_value(i,q) *
                                rhs_values[q](component_i) *
                                fe_values.JxW(q);
              }
          }

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=i+1; j<dofs_per_cell; ++j)
            local_matrix(i,j) = local_matrix(j,i);

        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global (local_matrix, local_rhs,
                                                local_dof_indices,
                                                system_matrix, system_rhs);
      }

    if (solver_type != SolverType::UMFPACK)
      {
        pressure_mass_matrix.reinit(sparsity_pattern.block(1,1));
        pressure_mass_matrix.copy_from(system_matrix.block(1,1));
        system_matrix.block(1,1) = 0;
      }
  }

  template <int dim>
    void
   StokesProblem<dim>::assemble_system_mesh_loop()
    {
    system_matrix=0;
    system_rhs=0;

    typedef decltype(dof_handler.begin_active()) Iterator;
    const RightHandSide<dim> rhs_function;
    const Solution<dim> boundary_function;
    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);

    auto penalty_parameter = [] (const double degree,
        const double extent1,
        const double extent2) -> double
    {
        return 4.0 * degree * (degree+1.0) * 0.5 * (1.0/extent1 + 1.0/extent2);
      };

    auto cell_worker = [&] (const Iterator &cell, ScratchData<dim> &scratch_data, CopyData &copy_data)
    {
      scratch_data.fe_values.reinit (cell);

      const unsigned int dofs_per_cell   = scratch_data.fe_values.get_fe().dofs_per_cell;
      const unsigned int n_q_points      = scratch_data.fe_values.get_quadrature().size();

      copy_data.reinit(cell, dofs_per_cell);

      const FEValues<dim> &fe_v = scratch_data.fe_values;
      const std::vector<double> &JxW = fe_v.get_JxW_values ();

      const double nu = 1.0;
      std::vector<Vector<double> >      rhs_values (n_q_points,
                                                    Vector<double>(dim+1));
      rhs_function.vector_value_list (fe_v.get_quadrature_points(), rhs_values);
      Tensor<1,dim> force_f;

      for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
        {
          for (unsigned int d=0;d<dim;++d)
          force_f[d] = rhs_values[point](d);
        for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
          {
            for (unsigned int j=0; j<fe_v.dofs_per_cell; ++j)
              copy_data.cell_matrix(i,j) +=
                  (
                  // nu \nabla v : \nabla u
                  nu
                  * scalar_product(fe_v[velocities].gradient (i,point),
                                    fe_v[velocities].gradient (j,point))
                  // -q, div u
                  - fe_v[pressure].value (i, point)
                    * fe_v[velocities].divergence (j, point)
                    // -p, div v
                    - fe_v[pressure].value (j, point)
                    * fe_v[velocities].divergence (i, point)

                    // p,q
                    + fe_v[pressure].value (j, point)
                    * fe_v[pressure].value (i, point)
                  ) * JxW[point];

            copy_data.cell_rhs(i) +=
                // f,v
                (force_f * fe_v[velocities].value(i, point))
                * JxW[point];
          }
        }
    };

    auto boundary_worker = [&] (const Iterator &cell, const unsigned int &face_no, ScratchData<dim> &scratch_data, CopyData &copy_data)
    {
        /*
      FEFaceValues<dim> &fe_fv = scratch_data.internal_fe_face_values;
      fe_fv.reinit(cell, face_no);

      const auto &q_points = fe_fv.get_quadrature_points();

      const std::vector<double> &JxW = fe_fv.get_JxW_values ();
      const std::vector<Tensor<1,dim> > &normals = fe_fv.get_normal_vectors ();

      const double nu = 1.0;

      std::vector<Vector<double> >      g_values (q_points.size(),
                                                    Vector<double>(dim+1));
      boundary_function.vector_value_list (q_points, g_values);
      Tensor<1,dim> g;

      const double degree = std::max(1.0, static_cast<double>(fe_fv.get_fe().degree));
      const double extent1 = cell->extent_in_direction(GeometryInfo<dim>::unit_normal_direction[face_no]);
      const double penalty = penalty_parameter(degree, extent1, extent1);

      for (unsigned int point=0; point<q_points.size(); ++point)
        {
          for (unsigned int d=0;d<dim;++d)
            g[d] = g_values[point](d);

          for (unsigned int i=0; i<fe_fv.dofs_per_cell; ++i)
            for (unsigned int j=0; j<fe_fv.dofs_per_cell; ++j)
              copy_data.cell_matrix(i,j) +=
                  (
                    // - nu (\nabla u n) . v
                    - nu
                    * ((fe_fv[velocities].gradient(j,point) * normals[point])
                    * fe_fv[velocities].value(i,point))

                    // - nu u . (\nabla v n)  // NIPG: use +
                    - nu
                    * (fe_fv[velocities].value(j,point)
                    * (fe_fv[velocities].gradient(i,point) * normals[point]))

                    // + nu * penalty u . v
                    + nu
                    * penalty
                    * (fe_fv[velocities].value(j,point)
                    * fe_fv[velocities].value(i,point))

                    // p (v.n)
                    + fe_fv[pressure].value(j, point)
                    * scalar_product(fe_fv[velocities].value(i,point),
                                     normals[point])

                    // q (u.n)
                    + fe_fv[pressure].value(i, point)
                    * scalar_product(fe_fv[velocities].value(j,point),
                                     normals[point])

                    ) * JxW[point];

          for (unsigned int i=0; i<fe_fv.dofs_per_cell; ++i)
            copy_data.cell_rhs(i) +=
                (
                  // -nu g . (\nabla v n) // NIPG: use +
                  - nu
                  * scalar_product(g, (fe_fv[velocities].gradient(i,point) * normals[point]))

                  // +nu penalty g . v
                  + nu
                  * penalty
                  * scalar_product(g, fe_fv[velocities].value(i,point))

                  // q (g.n) (weak normal component of boundary condition)
                  + fe_fv[pressure].value(i, point)
                  * scalar_product(g, normals[point])
                  ) * JxW[point];

        }*/
    };

    auto face_worker = [&]
                       (const Iterator &cell, const unsigned int &f, const unsigned int &sf,
                        const Iterator &ncell, const unsigned int &nf, const unsigned int &nsf,
                        ScratchData<dim> &scratch_data, CopyData &copy_data)
    {
      auto & fe_fv = scratch_data.fe_facet_values;
      fe_fv.reinit(cell, f, sf, ncell, nf, nsf);

      copy_data.face_data.emplace_back();
      CopyDataFace &copy_data_face = copy_data.face_data.back();
      const unsigned int dofs_per_cell   = scratch_data.fe_values.get_fe().dofs_per_cell;

      copy_data_face.joint_dof_indices = scratch_data.fe_facet_values.get_joint_dof_indices();


      copy_data_face.cell_matrix.reinit(2*dofs_per_cell, 2*dofs_per_cell);

      const FEFaceValuesBase<dim> &fe_v = *scratch_data.fe_face_values;
      const FEFaceValuesBase<dim> &fe_v_neighbor = *scratch_data.fe_face_values_neighbor;

      const std::vector<double> &JxW = fe_v.get_JxW_values ();
      const std::vector<Tensor<1,dim> > &normals = fe_v.get_normal_vectors ();

      double nu = 1.0;

      const double degree = std::max(1.0, static_cast<double>(fe_v.get_fe().degree));
      const double extent1 = cell->extent_in_direction(GeometryInfo<dim>::unit_normal_direction[f])
          * (cell->has_children() ? 2.0 : 1.0);
      const double extent2 = ncell->extent_in_direction(GeometryInfo<dim>::unit_normal_direction[nf])
          * (ncell->has_children() ? 2.0 : 1.0);
      const double penalty = penalty_parameter(degree, extent1, extent2);

      for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
        {
          for (unsigned int i=0; i<2*fe_v.dofs_per_cell; ++i)
            for (unsigned int j=0; j<2*fe_v.dofs_per_cell; ++j)
              copy_data_face.cell_matrix(i,j) +=
                  (
                  // - nu {\nabla u}n . [v] (consistency)
                  - nu
                  * (scratch_data.gradient_n_avg(j, point, velocities)
                  * scratch_data.jump(i, point, velocities))

                    // - nu [u] . {\nabla v}n  (symmetry) // NIPG: use +
                    - nu
                    * (scratch_data.jump(j, point, velocities)
                    * scratch_data.gradient_n_avg(i, point, velocities))

                    // nu sigma [u].[v] (penalty)
                    + nu * penalty
                    * scalar_product(scratch_data.jump(j, point, velocities),
                    scratch_data.jump(i, point, velocities))

                    // {p} ([v].n)
                    + scratch_data.avg(j, point, pressure)
                    * scalar_product(scratch_data.jump(i, point, velocities),
                                     normals[point])

                    // {q} ([u].n)
                    + scratch_data.avg(i, point, pressure)
                    * scalar_product(scratch_data.jump(j, point, velocities),
                                     normals[point])

                  ) * JxW[point];
        }

    };

    auto copier = [&] (const CopyData &c)
    {
        copy(c, constraints, system_matrix, system_rhs);
      };

    const unsigned int n_gauss_points = pressure_degree+2;

    static MappingQ1<dim> mapping;
    ScratchData<dim> scratch_data(mapping, fe, n_gauss_points);
    CopyData cd;
    MeshWorker::mesh_loop(dof_handler.begin_active(),
                          dof_handler.end(),
                          cell_worker,
                          copier,
                          scratch_data,
                          cd,
                          MeshWorker::assemble_own_cells | MeshWorker::assemble_boundary_faces | MeshWorker::assemble_own_interior_faces_once,
                          boundary_worker,
                          face_worker);

    pressure_mass_matrix.reinit(sparsity_pattern.block(1,1));
    pressure_mass_matrix.copy_from(system_matrix.block(1,1));
    system_matrix.block(1,1) = 0;
  }

  template <int dim>
    void
	 StokesProblem<dim>::assemble_system_DG()
    {
	  system_matrix=0;

	  static MappingQ1<dim> mapping;
      MeshWorker::IntegrationInfoBox<dim> info_box;
      UpdateFlags update_flags = update_quadrature_points | update_values | update_gradients;
      info_box.add_update_flags_all(update_flags);
      info_box.initialize(fe, mapping);

      MeshWorker::DoFInfo<dim> dof_info(dof_handler);

      MeshWorker::Assembler::MatrixSimple<BlockSparseMatrix<double> > assembler;
      assembler.initialize(system_matrix);
//      assembler.initialize(constraints);

      StokesMatrixIntegrator<dim> integrator;
      MeshWorker::integration_loop<dim, dim>(
        dof_handler.begin_active(), dof_handler.end(),
        dof_info, info_box,
        integrator, assembler);

      pressure_mass_matrix.reinit(sparsity_pattern.block(1,1));
      pressure_mass_matrix.copy_from(system_matrix.block(1,1));
      system_matrix.block(1,1) = 0;

    }

  template <int dim>
    void
	StokesProblem<dim>::assemble_rhs_DG()
    {
	  system_rhs=0;

	  static MappingQ1<dim> mapping;
      MeshWorker::IntegrationInfoBox<dim> info_box;
      UpdateFlags update_flags = update_quadrature_points | update_values | update_gradients;
      info_box.add_update_flags_all(update_flags);
      info_box.initialize(fe, mapping);

      MeshWorker::DoFInfo<dim> dof_info(dof_handler);

      // Since this assembler allows us to fill several vectors, the interface is
      // a little more complicated as above. The pointers to the vectors have to
      // be stored in a AnyData object. While this seems to cause two extra
      // lines of code here, it actually comes handy in more complex
      // applications.
      MeshWorker::Assembler::ResidualSimple<BlockVector<double> > assembler;
      AnyData data;
      data.add<BlockVector<double>*>(&system_rhs, "RHS");
      assembler.initialize(data);

      RHSIntegrator<dim> integrator;
      MeshWorker::integration_loop<dim, dim>(
        dof_handler.begin_active(), dof_handler.end(),
        dof_info, info_box,
        integrator, assembler);
    }

  template <int dim>
    void
	StokesProblem<dim>::assemble_matrix_rhs_DG()
    {
	  system_matrix=0;
	  system_rhs=0;

//	  static MappingQ1<dim> mapping;
	  static MappingQ<dim> mapping(1);
      MeshWorker::IntegrationInfoBox<dim> info_box;
      UpdateFlags update_flags = update_quadrature_points | update_values | update_gradients;
      info_box.add_update_flags_all(update_flags);
      info_box.initialize(fe, mapping);

      MeshWorker::DoFInfo<dim> dof_info(dof_handler);

      // Since this assembler allows us to fill several vectors, the interface is
      // a little more complicated as above. The pointers to the vectors have to
      // be stored in a AnyData object. While this seems to cause two extra
      // lines of code here, it actually comes handy in more complex
      // applications.
      MeshWorker::Assembler::SystemSimple<BlockSparseMatrix<double>, BlockVector<double> > assembler;

      assembler.initialize(system_matrix, system_rhs);
      assembler.initialize(constraints);

      StokesMatrixIntegrator<dim> integrator;
      MeshWorker::integration_loop<dim, dim>(
        dof_handler.begin_active(), dof_handler.end(),
        dof_info, info_box,
        integrator, assembler);

      pressure_mass_matrix.reinit(sparsity_pattern.block(1,1));
      pressure_mass_matrix.copy_from(system_matrix.block(1,1));
      system_matrix.block(1,1) = 0;

//      std::cout << " RHS block(1) norm: " << system_rhs.block(1).l2_norm() << std::endl;
//      std::cout << " Matrix block(1, 0) norm: " << system_matrix.block(1,0).linfty_norm() << std::endl;
//      std::cout << " Matrix block(1, 1) norm: " << system_matrix.block(1,1).linfty_norm() << std::endl;

    }



  // @sect4{StokesProblem::assemble_multigrid}

  // Here, like in step-16, we have a function that assembles the level
  // and interface matrices necessary for the multigrid preconditioner.
  template <int dim>
  void StokesProblem<dim>::assemble_multigrid ()
  {}

// @sect4{StokesProblem::solve}

// This function sets up things differently based on if you want to use ILU
// or GMG as a preconditioner.  Both methods share the same solver (FGMRES)
// but require a different preconditioner to be initialized. Here we time not
// only the entire solve function, but we separately time the setup of the
// preconditioner as well as the solve itself.
  template <int dim>
  void StokesProblem<dim>::solve ()
  {
    // Here we must make sure to solve for the residual with "good enough" accuracy
    SolverControl solver_control (system_matrix.m(),
                                  1e-12*system_rhs.l2_norm());

    // This is used to pass whether or not we want to solve for A inside
    // the preconditioner.  One could change this to false to see if
    // there is still convergence and if so does the program then run
    // faster or slower
    const bool use_expensive = true;

    SolverFGMRES<BlockVector<double> > solver (solver_control);

        computing_timer.enter_subsection ("(ILU specific)");
        computing_timer.enter_subsection ("Solve - Set-up Preconditioner");

//        std::cout << "   Computing preconditioner..." << std::endl << std::flush;

        SparseDirectUMFPACK  A_preconditioner;
        A_preconditioner.initialize(system_matrix.block(0,0));

//        SparseILU<double> A_preconditioner;
//        A_preconditioner.initialize (system_matrix.block(0,0));

        SparseILU<double> S_preconditioner;
        S_preconditioner.initialize (pressure_mass_matrix);

        const BlockSchurPreconditioner<SparseDirectUMFPACK, SparseILU<double> >
        preconditioner (system_matrix,
                        pressure_mass_matrix,
                        A_preconditioner,
                        S_preconditioner,
                        use_expensive);

        computing_timer.leave_subsection();
        computing_timer.leave_subsection();

        {
          TimerOutput::Scope solve_fmgres(computing_timer, "Solve - FGMRES");

          solution = 0;
          solver.solve (system_matrix,
                        solution,
                        system_rhs,
                        preconditioner);
        }

    constraints.distribute (solution);


//          unsigned int n_iterations_A = preconditioner.n_iterations_A;
//          unsigned int n_iterations_S = preconditioner.n_iterations_S;
//    std::cout << std::endl
//              << "\tNumber of FGMRES iterations: "
//              << solver_control.last_step() << std::endl
//              << "\tTotal number of iterations used for approximation of A inverse: "
//              << n_iterations_A << std::endl
//              << "\tTotal number of iterations used for approximation of S inverse: "
//              << n_iterations_S << std::endl
//              << std::endl;
  }


// @sect4{StokesProblem::process_solution}

// This function computes the L2 and H1 errors of the solution. For this,
// we need to make sure the pressure has mean zero.
  template <int dim>
  void StokesProblem<dim>::compute_errors (unsigned int k)
  {
    const ComponentSelectFunction<dim> velocity_mask(std::make_pair(0, dim), dim+1);
    const ComponentSelectFunction<dim> pressure_mask(dim, 1.0, dim+1);


    Vector<float> difference_per_cell (triangulation.n_active_cells());
    VectorTools::integrate_difference (dof_handler,
                                       solution,
                                       Solution<dim>(),
                                       difference_per_cell,
                                       QGauss<dim>(pressure_degree+3),
                                       VectorTools::L2_norm,
                                       &velocity_mask);

//    const double Velocity_L2_error = VectorTools::compute_global_error(triangulation,
//  		  	  	  	  	  	  	  	  	  	  	  	  	  	difference_per_cell,
//                                                              VectorTools::L2_norm);

    const double Velocity_L2_error = difference_per_cell.l2_norm();

    global_l2_diff = difference_per_cell;


    VectorTools::integrate_difference (dof_handler,
                                       solution,
                                       Solution<dim>(),
                                       difference_per_cell,
                                       QGauss<dim>(pressure_degree+3),
                                       VectorTools::H1_norm,
                                       &velocity_mask);

//    const double Velocity_H1_error = VectorTools::compute_global_error(triangulation,
//  		  	  	  	  	  	  	  	  	  	  	  	  	  	difference_per_cell,
//                                                              VectorTools::H1_norm);

    const double Velocity_H1_error = difference_per_cell.l2_norm();

    global_h1_diff = difference_per_cell;

    VectorTools::integrate_difference (dof_handler,
                                       solution,
                                       ZeroFunction<dim>(dim+1),
                                       difference_per_cell,
                                       QGauss<dim>(pressure_degree+3),
                                       VectorTools::Hdiv_seminorm,
                                       &velocity_mask);

//    const double Velocity_Hdiv_error1 = VectorTools::compute_global_error(triangulation,
//  		  	  	  	  	  	  	  	  	  	  	  	  	  	difference_per_cell,
//                                                              VectorTools::Hdiv_seminorm);

    const double Velocity_Hdiv_error1 = difference_per_cell.l2_norm();

    global_div_diff = difference_per_cell;

    VectorTools::integrate_difference (dof_handler,
                                       solution,
                                       Solution<dim>(),
                                       difference_per_cell,
                                       QGauss<dim>(pressure_degree+3),
                                       VectorTools::L2_norm,
                                       &pressure_mask);
    const double Pressure_L2_error = VectorTools::compute_global_error(triangulation,
                                                            difference_per_cell,
                                                                 VectorTools::L2_norm);

    static double last_Pressure_L2_error = 0;

    VectorTools::integrate_difference (dof_handler,
                                       solution,
                                       ZeroFunction<dim>(dim+1),
                                       difference_per_cell,
                                       QGauss<dim>(pressure_degree+3),
                                       VectorTools::mean,
                                       &pressure_mask);
    const double Pressure_mean = VectorTools::compute_global_error(triangulation,
                                                            difference_per_cell,
                                                                 VectorTools::mean);



    std::vector<unsigned int> old_user_indices;
    triangulation.save_user_indices(old_user_indices);
    unsigned int i=0;
    for (typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell,++i)
      cell->set_user_index(i);

     	static MappingQ1<dim> mapping;
    	BlockVector<double> errors(1);
		errors.block(0).reinit(triangulation.n_active_cells());
		errors.collect_sizes();
   	    MeshWorker::IntegrationInfoBox<dim> info_box;
   	    const unsigned int n_gauss_points = dof_handler.get_fe().tensor_degree()+1;
   	    info_box.initialize_gauss_quadrature(n_gauss_points, n_gauss_points+1, n_gauss_points);
   	    AnyData solution_data;
   	    solution_data.add<BlockVector<double>*>(&solution, "solution");
   	    info_box.cell_selector.add("solution", true, true, false);
//   	    info_box.boundary_selector.add("solution", true, false, false);
//   	    info_box.face_selector.add("solution", true, false, false);
   	    info_box.add_update_flags_cell(update_quadrature_points);
   	    info_box.add_update_flags_boundary(update_quadrature_points);
   	    info_box.initialize(fe, mapping, solution_data, solution);
   	    MeshWorker::DoFInfo<dim> dof_info(dof_handler);
   	    MeshWorker::Assembler::CellsAndFaces<double> assembler;
   	    AnyData out_data;
   	    out_data.add<BlockVector<double>* >(&errors, "cells");
   	    assembler.initialize(out_data, false);
   	    ErrorIntegrator<dim> integrator;
   	    MeshWorker::integration_loop<dim, dim> (
   	      dof_handler.begin_active(), dof_handler.end(),
   	      dof_info, info_box,
   	      integrator, assembler);
   	    triangulation.load_user_indices(old_user_indices);
     	cell_div_diff = errors.block(0);
//   	    double Velocity_Hdiv_error2 = errors.block(0).l2_norm();

	std::cout   << " At " << k+1 << "th mesh" << std::endl
//	            << " DoFs: " << dof_handler.n_dofs() << std::endl
	            << " L2 error:     " << std::setw(12) << Velocity_L2_error  << std::setw(0)
             	<< " L2_Conv_rate: " << std::setw(6)<< (k==0? 0:last_l2_error/Velocity_L2_error) << std::endl
				<< " H1 error:     " << std::setw(12) << Velocity_H1_error << std::setw(0)
				<< " H1_Conv_rate: " << std::setw(6)<< (k==0? 0:last_H1_error/Velocity_H1_error) << std::endl
				<< " Hdiv error1:  " << std::setw(12) << Velocity_Hdiv_error1 << std::setw(0)
				<< " Hdiv_Conv_rate1: " << std::setw(6)<< (k==0? 0:last_Hdiv_error1/Velocity_Hdiv_error1) << std::endl
//				<< " Hdiv error2:  " << std::setw(6) << Velocity_Hdiv_error2 << std::setw(0)
//				<< " Hdiv_Conv_rate2: " << std::setw(6)<< (k==0? 0:last_Hdiv_error2/Velocity_Hdiv_error2) << std::endl
        << " L2 pressure: " << std::setw(12) << Pressure_L2_error << std::setw(0)
        << " rate: " << std::setw(6)<< (k==0? 0:last_Pressure_L2_error/Pressure_L2_error) << std::endl
        << " pressure mean: " << std::setw(12) << Pressure_mean << std::setw(0) << std::endl
				<< "         *          " << std::endl;
	last_l2_error = Velocity_L2_error;
	last_H1_error = Velocity_H1_error;
	last_Hdiv_error1 = Velocity_Hdiv_error1;
//	last_Hdiv_error2 = Velocity_Hdiv_error2;
  last_Pressure_L2_error = Pressure_L2_error;


  }


// @sect4{StokesProblem::output_results}

// This function generates graphical output like it is done in step-22.
  template <int dim>
  void
  StokesProblem<dim>::output_results (const unsigned int refinement_cycle)  const
  {
    std::vector<std::string> solution_names (dim, "velocity");
    solution_names.push_back ("pressure");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation
    (dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation
    .push_back (DataComponentInterpretation::component_is_scalar);

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (solution, solution_names,
                              DataOut<dim>::type_dof_data,
                              data_component_interpretation);

//    Vector<double> ref(dof_handler.n_dofs());
//    VectorTools::interpolate(dof_handler, Solution<dim>(), ref);
//
//    {
//      std::vector<std::string> solution_names (dim, "ref_u");
//    solution_names.push_back ("ref_p");
//
//    data_out.add_data_vector (ref, solution_names,
//                              DataOut<dim>::type_dof_data,
//                              data_component_interpretation);
//    }

//    Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
//    FEValuesExtractors::Vector velocity(0);
//    KellyErrorEstimator<dim>::estimate (dof_handler,
//                                        QGauss<dim-1>(pressure_degree+2),
//                                        typename FunctionMap<dim>::type(),
//                                        solution,
//                                        estimated_error_per_cell,
//                                        fe.component_mask(velocity));
//
    data_out.add_data_vector (global_div_diff, "global_div_diff",
                              DataOut<dim>::type_cell_data);
    data_out.add_data_vector (cell_div_diff, "cell_div_diff",
                              DataOut<dim>::type_cell_data);
    data_out.add_data_vector (global_l2_diff, "global_l2_diff",
                              DataOut<dim>::type_cell_data);
    data_out.add_data_vector (global_h1_diff, "global_h1_diff",
                              DataOut<dim>::type_cell_data);

    data_out.build_patches (3);

    std::ostringstream filename;
    filename << "solution-"
             << Utilities::int_to_string (refinement_cycle, 2)
             << ".vtk";

    std::ofstream output (filename.str().c_str());
    data_out.write_vtk (output);
  }



// @sect4{StokesProblem::run}

// The last step in the Stokes class is, as usual, the function that
// generates the initial grid and calls the other functions in the
// respective order.
  template <int dim>
  void StokesProblem<dim>::run ()
  {
	if (geo_type == GeoType::Cube)
	  {
	    
	    	GridGenerator::hyper_cube (triangulation);
	    /*	      GridIn<dim> gridin;
	      gridin.attach_triangulation(triangulation);
	      std::ifstream f("square.msh");
	      gridin.read_msh(f);*/
	    
	  }
	

  if (geo_type == GeoType::Circle)
	{
		GridGenerator::hyper_ball (triangulation);
		static const SphericalManifold<dim> boundary;
		triangulation.set_all_manifold_ids_on_boundary (0);
		triangulation.set_manifold (0, boundary);
	}
  if (geo_type == GeoType::L_shape)
	  {
      GridGenerator::hyper_L (triangulation);
    }
  if (geo_type == GeoType::L_shape_diagonal)
    {
      std::vector<unsigned int> repetitions(dim, 1);
      repetitions[0] = 2;
      Point<dim> p1 (-1.0, 0.0);
      Point<dim> p2 (1.0, -1.0);
      GridGenerator::subdivided_hyper_rectangle (triangulation,
                                                 repetitions,
                                                 p1,
                                                 p2, false);
      auto cell=triangulation.begin();
      cell->vertex(0)(1)=1.0;
      cell->vertex(2)(1)=1.0;
      cell->vertex(2)(0)=0.0;
      cell->vertex(1)(0)=-1.0;
      //cell->vertex(3)(0)+=0.2;

      //cell->vertex(2)(0)=0.0;

      std::ofstream out ("grid.eps");
       GridOut grid_out;
       grid_out.write_eps (triangulation, out);
	    
	  }
	

	if (false)
		GridTools::distort_random(0.2, triangulation);

    triangulation.refine_global (1);

    std::cout << "  Now running with "<< fe.get_name() << std::endl;

    for (unsigned int refinement_cycle = 0; refinement_cycle<6;
         ++refinement_cycle)
      {
//        std::cout << "Refinement cycle " << refinement_cycle << std::endl;

        if (refinement_cycle > 0)
          {
            int ref_type = 0;
            switch (ref_type)
              {
              case 0: // global
                {
//                  std::cout << " global..." << std::endl;
                  triangulation.refine_global();
                  break;
                }
              case 1: // kelly u
                {
//                  std::cout << " Kelly..." << std::endl;
              	  Vector<float> estimated_error_per_cell (triangulation.n_global_active_cells());
              	  FEValuesExtractors::Vector velocity(0);

                  KellyErrorEstimator<dim>::estimate(dof_handler,
                                                       QGauss<dim-1>(fe.degree+2),
                                                       typename FunctionMap<dim>::type(),
                                                       solution,
                                                       estimated_error_per_cell,
                                                       fe.component_mask(velocity)
                                                      );

                  GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                      estimated_error_per_cell,
                                                      0.3, 0.0);

                  triangulation.execute_coarsening_and_refinement ();

                  break;
                }
              case 2: // res_based
              {
//            	  std::cout << " Residual based..." << std::endl;
              	std::vector<unsigned int> old_user_indices;
              	triangulation.save_user_indices(old_user_indices);
            	BlockVector<double> estimates(1);
        		estimates.block(0).reinit(triangulation.n_active_cells());
        		estimates.collect_sizes();
        		unsigned int i=0;
        		for (typename Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active();
        		 	   cell != triangulation.end(); ++cell,++i)
        		 	   cell->set_user_index(i);

        		MeshWorker::IntegrationInfoBox<dim> info_box;
        		const unsigned int n_gauss_points = dof_handler.get_fe().tensor_degree()+1;
        		info_box.initialize_gauss_quadrature(n_gauss_points, n_gauss_points+1, n_gauss_points);

        		AnyData solution_data;
        		solution_data.add<const BlockVector<double>*>(&solution, "solution");

        		info_box.cell_selector.add("solution", true, true, true);
        		info_box.boundary_selector.add("solution", true, true, false);
        		info_box.face_selector.add("solution", true, true, false);

        		info_box.add_update_flags_cell(update_quadrature_points);
        		info_box.add_update_flags_boundary(update_quadrature_points);
        		static MappingQ1<dim> mapping;
        		info_box.initialize(fe, mapping , solution_data, solution);
        		MeshWorker::DoFInfo<dim> dof_info(dof_handler);

        		MeshWorker::Assembler::CellsAndFaces<double> assembler;
        		AnyData out_data;
        		out_data.add<BlockVector<double>*>(&estimates, "cells");
        		assembler.initialize(out_data, false);
                Estimator<dim> integrator(1.0);
        	    MeshWorker::integration_loop<dim, dim> (dof_handler.begin_active(),
        	    		  	  	  	  	  	  	  	  	 dof_handler.end(),
        	 											 dof_info, info_box,
        												 integrator, assembler);

        	    GridRefinement::refine_and_coarsen_fixed_number (triangulation,
        	                                                      estimates.block(0),
        	                                                      0.3, 0.0);

        	    triangulation.load_user_indices(old_user_indices);

        	 	triangulation.execute_coarsening_and_refinement ();
//        	    std::cout << " estimates: " << estimates.block(0).l2_norm() << std::endl;
        		break;
              }

              }
          }

//        std::cout << "   Set-up..." << std::endl;
        setup_dofs();

//        std::cout << "   Assembling..." << std::endl;
        int assemble_type = 2;
        switch (assemble_type)
        {
        case 0:
        {
        	assemble_system ();
        	break;
        }
        case 1:
        {
//        	assemble_system_DG();
//        	assemble_rhs_DG();
        	assemble_matrix_rhs_DG();
        	break;
        }
          case 2:
          {
              assemble_system_mesh_loop ();
              break;
            }
        }

//        std::cout << "   Solving..." << std::flush;
        solve ();

        compute_errors (refinement_cycle);

        output_results (refinement_cycle);

        Utilities::System::MemoryStats mem;
        Utilities::System::get_memory_stats(mem);
//        std::cout << "   VM Peak: " << mem.VmPeak << std::endl;
//
//        computing_timer.print_summary ();
        computing_timer.reset ();
      }
  }
}

// @sect3{The main function}
int main ()
{
  try
    {
      using namespace dealii;
      using namespace Step59;

      deallog.depth_console(0);


      const int dim = 2;

      for (unsigned int d=0;d<3;++d)
      {
          FE_RaviartThomas<dim> test(d);

      std::cout << test.get_name()
                << ": "
                << test.dofs_per_face
                << " "
                << test.dofs_per_cell
                << std::endl;
        }
      for (unsigned int d=1;d<3;++d)
      {
          FE_BDM<dim> test(d);

      std::cout << test.get_name()
                << ": "
                << test.dofs_per_face
                << " "
                << test.dofs_per_cell
                << std::endl;
        }

      const int degree = 2;

      //FESystem<dim> fe(FESystem<dim>(FE_Q<dim>(degree), dim), 1, FE_Q<dim>(degree-1), 1);
      //FESystem<dim> fe(FESystem<dim>(FE_DGQ<dim>(degree), dim), 1, FE_DGQ<dim>(degree-1), 1);
      FESystem<dim> fe(FE_RaviartThomas<dim>(degree), 1, FE_DGQ<dim>(degree), 1);
      //FESystem<dim> fe(FE_BDM<dim>(degree), 1, FE_DGP<dim>(degree-1), 1);
            std::cout << fe.degree
                      << " " << fe.tensor_degree()
                      << std::endl;
      StokesProblem<dim> flow_problem(fe,
    		                          degree,
									  SolverType::FGMRES_ILU,
                    GeoType::Circle);

      flow_problem.run ();
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

