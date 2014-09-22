/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2009 - 2013 by the deal.II authors
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

 *
 * Author: Yaqi Wang, Texas A&M University, 2009, 2010
 */


// @sect3{Include files}

// We start with a bunch of include files that have already been explained in
// previous tutorial programs:
#include <deal.II/base/timer.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <fstream>
#include <iostream>

#include <deal.II/base/utilities.h>

// We use the next include file to access block vectors which provide us a
// convenient way to manage solution and right hand side vectors of all energy
// groups:
#include <deal.II/lac/block_vector.h>

// This include file is for transferring solutions from one mesh to another
// different mesh. We use it when we are initializing solutions after each
// mesh iteration:
#include <deal.II/numerics/solution_transfer.h>

// When integrating functions defined on one mesh against shape functions
// defined on a different mesh, we need a function @p get_finest_common_cells
// (as discussed in the introduction) which is defined in the following header
// file:
#include <deal.II/grid/grid_tools.h>

// Here are two more C++ standard headers that we use to define list data
// types as well as to fine-tune the output we generate:
#include <list>
#include <iomanip>

// The last step is as in all previous programs:
namespace Step28
{
  using namespace dealii;


  // @sect3{Material data}

  // First up, we need to define a class that provides material data
  // (including diffusion coefficients, removal cross sections, scattering
  // cross sections, fission cross sections and fission spectra) to the main
  // class.
  //
  // The parameter to the constructor determines for how many energy groups we
  // set up the relevant tables. At present, this program only includes data
  // for 2 energy groups, but a more sophisticated program may be able to
  // initialize the data structures for more groups as well, depending on how
  // many energy groups are selected in the parameter file.
  //
  // For each of the different coefficient types, there is one function that
  // returns the value of this coefficient for a particular energy group (or
  // combination of energy groups, as for the distribution cross section
  // $\chi_g\nu\Sigma_{f,g'}$ or scattering cross section $\Sigma_{s,g'\to
  // g}$). In addition to the energy group or groups, these coefficients
  // depend on the type of fuel or control rod, as explained in the
  // introduction. The functions therefore take an additional parameter, @p
  // material_id, that identifies the particular kind of rod. Within this
  // program, we use <code>n_materials=8</code> different kinds of rods.
  //
  // Except for the scattering cross section, each of the coefficients
  // therefore can be represented as an entry in a two-dimensional array of
  // floating point values indexed by the energy group number as well as the
  // material ID. The Table class template is the ideal way to store such
  // data. Finally, the scattering coefficient depends on both two energy
  // group indices and therefore needs to be stored in a three-dimensional
  // array, for which we again use the Table class, where this time the first
  // template argument (denoting the dimensionality of the array) of course
  // needs to be three:
  class MaterialData
  {
  public:
    MaterialData (const unsigned int n_groups);

    double get_diffusion_coefficient (const unsigned int group,
                                      const unsigned int material_id) const;
    double get_removal_XS (const unsigned int group,
                           const unsigned int material_id) const;
    double get_fission_XS (const unsigned int group,
                           const unsigned int material_id) const;
    double get_fission_dist_XS (const unsigned int group_1,
                                const unsigned int group_2,
                                const unsigned int material_id) const;
    double get_scattering_XS (const unsigned int group_1,
                              const unsigned int group_2,
                              const unsigned int material_id) const;
    double get_fission_spectrum (const unsigned int group,
                                 const unsigned int material_id) const;

  private:
    const unsigned int n_groups;
    const unsigned int n_materials;

    Table<2,double> diffusion;
    Table<2,double> sigma_r;
    Table<2,double> nu_sigma_f;
    Table<3,double> sigma_s;
    Table<2,double> chi;
  };

  // The constructor of the class is used to initialize all the material data
  // arrays. It takes the number of energy groups as an argument (an throws an
  // error if that value is not equal to two, since at presently only data for
  // two energy groups is implemented; however, using this, the function
  // remains flexible and extendable into the future). In the member
  // initialization part at the beginning, it also resizes the arrays to their
  // correct sizes.
  //
  // At present, material data is stored for 8 different types of
  // material. This, as well, may easily be extended in the future.
  MaterialData::MaterialData (const unsigned int n_groups)
    :
    n_groups (n_groups),
    n_materials (8),
    diffusion (n_materials, n_groups),
    sigma_r (n_materials, n_groups),
    nu_sigma_f (n_materials, n_groups),
    sigma_s (n_materials, n_groups, n_groups),
    chi (n_materials, n_groups)
  {
    switch (n_groups)
      {
      case 2:
      {
        for (unsigned int m=0; m<n_materials; ++m)
          {
            diffusion[m][0] = 1.2;
            diffusion[m][1] = 0.4;
            chi[m][0]       = 1.0;
            chi[m][1]       = 0.0;
            sigma_r[m][0]   = 0.03;
            for (unsigned int group_1=0; group_1<n_groups; ++group_1)
              for (unsigned int group_2=0; group_2<n_groups; ++ group_2)
                sigma_s[m][group_1][group_2]   = 0.0;
          }


        diffusion[5][1]  = 0.2;

        sigma_r[4][0]    = 0.026;
        sigma_r[5][0]    = 0.051;
        sigma_r[6][0]    = 0.026;
        sigma_r[7][0]    = 0.050;

        sigma_r[0][1]    = 0.100;
        sigma_r[1][1]    = 0.200;
        sigma_r[2][1]    = 0.250;
        sigma_r[3][1]    = 0.300;
        sigma_r[4][1]    = 0.020;
        sigma_r[5][1]    = 0.040;
        sigma_r[6][1]    = 0.020;
        sigma_r[7][1]    = 0.800;

        nu_sigma_f[0][0] = 0.0050;
        nu_sigma_f[1][0] = 0.0075;
        nu_sigma_f[2][0] = 0.0075;
        nu_sigma_f[3][0] = 0.0075;
        nu_sigma_f[4][0] = 0.000;
        nu_sigma_f[5][0] = 0.000;
        nu_sigma_f[6][0] = 1e-7;
        nu_sigma_f[7][0] = 0.00;

        nu_sigma_f[0][1] = 0.125;
        nu_sigma_f[1][1] = 0.300;
        nu_sigma_f[2][1] = 0.375;
        nu_sigma_f[3][1] = 0.450;
        nu_sigma_f[4][1] = 0.000;
        nu_sigma_f[5][1] = 0.000;
        nu_sigma_f[6][1] = 3e-6;
        nu_sigma_f[7][1] = 0.00;

        sigma_s[0][0][1] = 0.020;
        sigma_s[1][0][1] = 0.015;
        sigma_s[2][0][1] = 0.015;
        sigma_s[3][0][1] = 0.015;
        sigma_s[4][0][1] = 0.025;
        sigma_s[5][0][1] = 0.050;
        sigma_s[6][0][1] = 0.025;
        sigma_s[7][0][1] = 0.010;

        break;
      }


      default:
        Assert (false,
                ExcMessage ("Presently, only data for 2 groups is implemented"));
      }
  }


  // Next are the functions that return the coefficient values for given
  // materials and energy groups. All they do is to make sure that the given
  // arguments are within the allowed ranges, and then look the respective
  // value up in the corresponding tables:
  double
  MaterialData::get_diffusion_coefficient (const unsigned int group,
                                           const unsigned int material_id) const
  {
    Assert (group < n_groups,
            ExcIndexRange (group, 0, n_groups));
    Assert (material_id < n_materials,
            ExcIndexRange (material_id, 0, n_materials));

    return diffusion[material_id][group];
  }



  double
  MaterialData::get_removal_XS (const unsigned int group,
                                const unsigned int material_id) const
  {
    Assert (group < n_groups,
            ExcIndexRange (group, 0, n_groups));
    Assert (material_id < n_materials,
            ExcIndexRange (material_id, 0, n_materials));

    return sigma_r[material_id][group];
  }


  double
  MaterialData::get_fission_XS (const unsigned int group,
                                const unsigned int material_id) const
  {
    Assert (group < n_groups,
            ExcIndexRange (group, 0, n_groups));
    Assert (material_id < n_materials,
            ExcIndexRange (material_id, 0, n_materials));

    return nu_sigma_f[material_id][group];
  }



  double
  MaterialData::get_scattering_XS (const unsigned int group_1,
                                   const unsigned int group_2,
                                   const unsigned int material_id) const
  {
    Assert (group_1 < n_groups,
            ExcIndexRange (group_1, 0, n_groups));
    Assert (group_2 < n_groups,
            ExcIndexRange (group_2, 0, n_groups));
    Assert (material_id < n_materials,
            ExcIndexRange (material_id, 0, n_materials));

    return sigma_s[material_id][group_1][group_2];
  }



  double
  MaterialData::get_fission_spectrum (const unsigned int group,
                                      const unsigned int material_id) const
  {
    Assert (group < n_groups,
            ExcIndexRange (group, 0, n_groups));
    Assert (material_id < n_materials,
            ExcIndexRange (material_id, 0, n_materials));

    return chi[material_id][group];
  }


  // The function computing the fission distribution cross section is slightly
  // different, since it computes its value as the product of two other
  // coefficients. We don't need to check arguments here, since this already
  // happens when we call the two other functions involved, even though it
  // would probably not hurt either:
  double
  MaterialData::get_fission_dist_XS (const unsigned int group_1,
                                     const unsigned int group_2,
                                     const unsigned int material_id) const
  {
    return (get_fission_spectrum(group_1, material_id) *
            get_fission_XS(group_2, material_id));
  }



  // @sect3{The <code>EnergyGroup</code> class}

  // The first interesting class is the one that contains everything that is
  // specific to a single energy group. To group things that belong together
  // into individual objects, we declare a structure that holds the
  // Triangulation and DoFHandler objects for the mesh used for a single
  // energy group, and a number of other objects and member functions that we
  // will discuss in the following sections.
  //
  // The main reason for this class is as follows: for both the forward
  // problem (with a specified right hand side) as well as for the eigenvalue
  // problem, one typically solves a sequence of problems for a single energy
  // group each, rather than the fully coupled problem. This becomes
  // understandable once one realizes that the system matrix for a single
  // energy group is symmetric and positive definite (it is simply a diffusion
  // operator), whereas the matrix for the fully coupled problem is generally
  // nonsymmetric and not definite. It is also very large and quite full if
  // more than a few energy groups are involved.
  //
  // Let us first look at the equation to solve in the case of an external
  // right hand side (for the time independent case): @f{eqnarray*} -\nabla
  // \cdot(D_g(x) \nabla \phi_g(x)) + \Sigma_{r,g}(x)\phi_g(x) =
  // \chi_g\sum_{g'=1}^G\nu\Sigma_{f,g'}(x)\phi_{g'}(x) + \sum_{g'\ne
  // g}\Sigma_{s,g'\to g}(x)\phi_{g'}(x) + s_{\mathrm{ext},g}(x) @f}
  //
  // We would typically solve this equation by moving all the terms on the
  // right hand side with $g'=g$ to the left hand side, and solving for
  // $\phi_g$. Of course, we don't know $\phi_{g'}$ yet, since the equations
  // for those variables include right hand side terms involving
  // $\phi_g$. What one typically does in such situations is to iterate:
  // compute @f{eqnarray*} -\nabla \cdot(D_g(x) \nabla \phi^{(n)}_g(x)) &+&
  // \Sigma_{r,g}(x)\phi^{(n)}_g(x) \\ &=&
  // \chi_g\sum_{g'=1}^{g-1}\nu\Sigma_{f,g'}(x)\phi^{(n)}_{g'}(x) +
  // \chi_g\sum_{g'=g}^G\nu\Sigma_{f,g'}(x)\phi^{(n-1)}_{g'}(x) + \sum_{g'\ne
  // g, g'<g}\Sigma_{s,g'\to g}(x)\phi^{(n)}_{g'}(x) + \sum_{g'\ne g,
  // g'>g}\Sigma_{s,g'\to g}(x)\phi^{(n-1)}_{g'}(x) + s_{\mathrm{ext},g}(x)
  // @f}
  //
  // In other words, we solve the equation one by one, using values for
  // $\phi_{g'}$ from the previous iteration $n-1$ if $g'\ge g$ and already
  // computed values for $\phi_{g'}$ from the present iteration if $g'<g$.
  //
  // When computing the eigenvalue, we do a very similar iteration, except
  // that we have no external right hand side and that the solution is scaled
  // after each iteration as explained in the introduction.
  //
  // In either case, these two cases can be treated jointly if all we do is to
  // equip the following class with these abilities: (i) form the left hand
  // side matrix, (ii) form the in-group right hand side contribution,
  // i.e. involving the extraneous source, and (iii) form that contribution to
  // the right hand side that stems from group $g'$. This class does exactly
  // these tasks (as well as some book-keeping, such as mesh refinement,
  // setting up matrices and vectors, etc). On the other hand, the class
  // itself has no idea how many energy groups there are, and in particular
  // how they interact, i.e. the decision of how the outer iteration looks
  // (and consequently whether we solve an eigenvalue or a direct problem) is
  // left to the NeutronDiffusionProblem class further down below in this
  // program.
  //
  // So let us go through the class and its interface:
  template <int dim>
  class EnergyGroup
  {
  public:

    // @sect5{Public member functions}
    //
    // The class has a good number of public member functions, since its the
    // way it operates is controlled from the outside, and therefore all
    // functions that do something significant need to be called from another
    // class. Let's start off with book-keeping: the class obviously needs to
    // know which energy group it represents, which material data to use, and
    // from what coarse grid to start. The constructor takes this information
    // and initializes the relevant member variables with that (see below).
    //
    // Then we also need functions that set up the linear system,
    // i.e. correctly size the matrix and its sparsity pattern, etc, given a
    // finite element object to use. The <code>setup_linear_system</code>
    // function does that. Finally, for this initial block, there are two
    // functions that return the number of active cells and degrees of freedom
    // used in this object -- using this, we can make the triangulation and
    // DoF handler member variables private, and do not have to grant external
    // use to it, enhancing encapsulation:
    EnergyGroup (const unsigned int        group,
                 const MaterialData       &material_data,
                 const Triangulation<dim> &coarse_grid,
                 const FiniteElement<dim> &fe);

    void setup_linear_system ();

    unsigned int n_active_cells () const;
    unsigned int n_dofs () const;

    // Then there are functions that assemble the linear system for each
    // iteration and the present energy group. Note that the matrix is
    // independent of the iteration number, so only has to be computed once
    // for each refinement cycle. The situation is a bit more involved for the
    // right hand side that has to be updated in each inverse power iteration,
    // and that is further complicated by the fact that computing it may
    // involve several different meshes as explained in the introduction. To
    // make things more flexible with regard to solving the forward or the
    // eigenvalue problem, we split the computation of the right hand side
    // into a function that assembles the extraneous source and in-group
    // contributions (which we will call with a zero function as source terms
    // for the eigenvalue problem) and one that computes contributions to the
    // right hand side from another energy group:
    void assemble_system_matrix ();
    void assemble_ingroup_rhs (const Function<dim> &extraneous_source);
    void assemble_cross_group_rhs (const EnergyGroup<dim> &g_prime);

    // Next we need a set of functions that actually compute the solution of a
    // linear system, and do something with it (such as computing the fission
    // source contribution mentioned in the introduction, writing graphical
    // information to an output file, computing error indicators, or actually
    // refining the grid based on these criteria and thresholds for refinement
    // and coarsening). All these functions will later be called from the
    // driver class <code>NeutronDiffusionProblem</code>, or any other class
    // you may want to implement to solve a problem involving the neutron flux
    // equations:
    void   solve ();

    double get_fission_source () const;

    void   output_results (const unsigned int cycle) const;

    void   estimate_errors (Vector<float> &error_indicators) const;

    void   refine_grid (const Vector<float> &error_indicators,
                        const double         refine_threshold,
                        const double         coarsen_threshold);

    // @sect5{Public data members}
    //
    // As is good practice in object oriented programming, we hide most data
    // members by making them private. However, we have to grant the class
    // that drives the process access to the solution vector as well as the
    // solution of the previous iteration, since in the power iteration, the
    // solution vector is scaled in every iteration by the present guess of
    // the eigenvalue we are looking for:
  public:

    Vector<double> solution;
    Vector<double> solution_old;


    // @sect5{Private data members}
    //
    // The rest of the data members are private. Compared to all the previous
    // tutorial programs, the only new data members are an integer storing
    // which energy group this object represents, and a reference to the
    // material data object that this object's constructor gets passed from
    // the driver class. Likewise, the constructor gets a reference to the
    // finite element object we are to use.
    //
    // Finally, we have to apply boundary values to the linear system in each
    // iteration, i.e. quite frequently. Rather than interpolating them every
    // time, we interpolate them once on each new mesh and then store them
    // along with all the other data of this class:
  private:

    const unsigned int            group;
    const MaterialData           &material_data;

    Triangulation<dim>            triangulation;
    const FiniteElement<dim>     &fe;
    DoFHandler<dim>               dof_handler;

    SparsityPattern               sparsity_pattern;
    SparseMatrix<double>          system_matrix;

    Vector<double>                system_rhs;

    std::map<types::global_dof_index,double> boundary_values;
    ConstraintMatrix              hanging_node_constraints;


    // @sect5{Private member functions}
    //
    // There is one private member function in this class. It recursively
    // walks over cells of two meshes to compute the cross-group right hand
    // side terms. The algorithm for this is explained in the introduction to
    // this program. The arguments to this function are a reference to an
    // object representing the energy group against which we want to integrate
    // a right hand side term, an iterator to a cell of the mesh used for the
    // present energy group, an iterator to a corresponding cell on the other
    // mesh, and the matrix that interpolates the degrees of freedom from the
    // coarser of the two cells to the finer one:
  private:

    void
    assemble_cross_group_rhs_recursive (const EnergyGroup<dim>                        &g_prime,
                                        const typename DoFHandler<dim>::cell_iterator &cell_g,
                                        const typename DoFHandler<dim>::cell_iterator &cell_g_prime,
                                        const FullMatrix<double>                       prolongation_matrix);
  };


  // @sect4{Implementation of the <code>EnergyGroup</code> class}

  // The first few functions of this class are mostly self-explanatory. The
  // constructor only sets a few data members and creates a copy of the given
  // triangulation as the base for the triangulation used for this energy
  // group. The next two functions simply return data from private data
  // members, thereby enabling us to make these data members private.
  template <int dim>
  EnergyGroup<dim>::EnergyGroup (const unsigned int        group,
                                 const MaterialData       &material_data,
                                 const Triangulation<dim> &coarse_grid,
                                 const FiniteElement<dim> &fe)
    :
    group (group),
    material_data (material_data),
    fe (fe),
    dof_handler (triangulation)
  {
    triangulation.copy_triangulation (coarse_grid);
    dof_handler.distribute_dofs (fe);
  }



  template <int dim>
  unsigned int
  EnergyGroup<dim>::n_active_cells () const
  {
    return triangulation.n_active_cells ();
  }



  template <int dim>
  unsigned int
  EnergyGroup<dim>::n_dofs () const
  {
    return dof_handler.n_dofs ();
  }



  // @sect5{<code>EnergyGroup::setup_linear_system</code>}
  //
  // The first "real" function is the one that sets up the mesh, matrices,
  // etc, on the new mesh or after mesh refinement. We use this function to
  // initialize sparse system matrices, and the right hand side vector. If the
  // solution vector has never been set before (as indicated by a zero size),
  // we also initialize it and set it to a default value. We don't do that if
  // it already has a non-zero size (i.e. this function is called after mesh
  // refinement) since in that case we want to preserve the solution across
  // mesh refinement (something we do in the
  // <code>EnergyGroup::refine_grid</code> function).
  template <int dim>
  void
  EnergyGroup<dim>::setup_linear_system ()
  {
    const unsigned int n_dofs = dof_handler.n_dofs();

    hanging_node_constraints.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler,
                                             hanging_node_constraints);
    hanging_node_constraints.close ();

    system_matrix.clear ();

    sparsity_pattern.reinit (n_dofs, n_dofs,
                             dof_handler.max_couplings_between_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
    hanging_node_constraints.condense (sparsity_pattern);
    sparsity_pattern.compress ();

    system_matrix.reinit (sparsity_pattern);

    system_rhs.reinit (n_dofs);

    if (solution.size() == 0)
      {
        solution.reinit (n_dofs);
        solution_old.reinit(n_dofs);
        solution_old = 1.0;
        solution = solution_old;
      }


    // At the end of this function, we update the list of boundary nodes and
    // their values, by first clearing this list and the re-interpolating
    // boundary values (remember that this function is called after first
    // setting up the mesh, and each time after mesh refinement).
    //
    // To understand the code, it is necessary to realize that we create the
    // mesh using the <code>GridGenerator::subdivided_hyper_rectangle</code>
    // function (in <code>NeutronDiffusionProblem::initialize_problem</code>)
    // where we set the last parameter to <code>true</code>. This means that
    // boundaries of the domain are "colored", i.e. the four (or six, in 3d)
    // sides of the domain are assigned different boundary indicators. As it
    // turns out, the bottom boundary gets indicator zero, the top one
    // boundary indicator one, and left and right boundaries get indicators
    // two and three, respectively.
    //
    // In this program, we simulate only one, namely the top right, quarter of
    // a reactor. That is, we want to interpolate boundary conditions only on
    // the top and right boundaries, while do nothing on the bottom and left
    // boundaries (i.e. impose natural, no-flux Neumann boundary
    // conditions). This is most easily generalized to arbitrary dimension by
    // saying that we want to interpolate on those boundaries with indicators
    // 1, 3, ..., which we do in the following loop (note that calls to
    // <code>VectorTools::interpolate_boundary_values</code> are additive,
    // i.e. they do not first clear the boundary value map):
    boundary_values.clear();

    for (unsigned int i=0; i<dim; ++i)
      VectorTools::interpolate_boundary_values (dof_handler,
                                                2*i+1,
                                                ZeroFunction<dim>(),
                                                boundary_values);
  }



  // @sect5{<code>EnergyGroup::assemble_system_matrix</code>}
  //
  // Next we need functions assembling the system matrix and right hand
  // sides. Assembling the matrix is straightforward given the equations
  // outlined in the introduction as well as what we've seen in previous
  // example programs. Note the use of <code>cell->material_id()</code> to get
  // at the kind of material from which a cell is made up of. Note also how we
  // set the order of the quadrature formula so that it is always appropriate
  // for the finite element in use.
  //
  // Finally, note that since we only assemble the system matrix here, we
  // can't yet eliminate boundary values (we need the right hand side vector
  // for this). We defer this to the <code>EnergyGroup::solve</code> function,
  // at which point all the information is available.
  template <int dim>
  void
  EnergyGroup<dim>::assemble_system_matrix ()
  {
    const QGauss<dim>  quadrature_formula(fe.degree + 1);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    |  update_gradients |
                             update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points    = quadrature_formula.size();

    FullMatrix<double> cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    for (; cell!=endc; ++cell)
      {
        cell_matrix = 0;

        fe_values.reinit (cell);

        const double diffusion_coefficient
          = material_data.get_diffusion_coefficient (group, cell->material_id());
        const double removal_XS
          = material_data.get_removal_XS (group,cell->material_id());

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              cell_matrix(i,j) += ((diffusion_coefficient *
                                    fe_values.shape_grad(i,q_point) *
                                    fe_values.shape_grad(j,q_point)
                                    +
                                    removal_XS *
                                    fe_values.shape_value(i,q_point) *
                                    fe_values.shape_value(j,q_point))
                                   *
                                   fe_values.JxW(q_point));

        cell->get_dof_indices (local_dof_indices);

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            system_matrix.add (local_dof_indices[i],
                               local_dof_indices[j],
                               cell_matrix(i,j));
      }

    hanging_node_constraints.condense (system_matrix);
  }



  // @sect5{<code>EnergyGroup::assemble_ingroup_rhs</code>}
  //
  // As explained in the documentation of the <code>EnergyGroup</code> class,
  // we split assembling the right hand side into two parts: the ingroup and
  // the cross-group couplings. First, we need a function to assemble the
  // right hand side of one specific group here, i.e. including an extraneous
  // source (that we will set to zero for the eigenvalue problem) as well as
  // the ingroup fission contributions.  (In-group scattering has already been
  // accounted for with the definition of removal cross section.) The
  // function's workings are pretty standard as far as assembling right hand
  // sides go, and therefore does not require more comments except that we
  // mention that the right hand side vector is set to zero at the beginning
  // of the function -- something we are not going to do for the cross-group
  // terms that simply add to the right hand side vector.
  template <int dim>
  void EnergyGroup<dim>::assemble_ingroup_rhs (const Function<dim> &extraneous_source)
  {
    system_rhs.reinit (dof_handler.n_dofs());

    const QGauss<dim>  quadrature_formula (fe.degree + 1);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values    |  update_quadrature_points  |
                             update_JxW_values);

    Vector<double>            cell_rhs (dofs_per_cell);
    std::vector<double>       extraneous_source_values (n_q_points);
    std::vector<double>       solution_old_values (n_q_points);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    for (; cell!=endc; ++cell)
      {
        cell_rhs = 0;

        fe_values.reinit (cell);

        const double fission_dist_XS
          = material_data.get_fission_dist_XS (group, group, cell->material_id());

        extraneous_source.value_list (fe_values.get_quadrature_points(),
                                      extraneous_source_values);

        fe_values.get_function_values (solution_old, solution_old_values);

        cell->get_dof_indices (local_dof_indices);

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            cell_rhs(i) += ((extraneous_source_values[q_point]
                             +
                             fission_dist_XS *
                             solution_old_values[q_point]) *
                            fe_values.shape_value(i,q_point) *
                            fe_values.JxW(q_point));

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          system_rhs(local_dof_indices[i]) += cell_rhs(i);
      }
  }



  // @sect5{<code>EnergyGroup::assemble_cross_group_rhs</code>}
  //
  // The more interesting function for assembling the right hand side vector
  // for the equation of a single energy group is the one that couples energy
  // group $g$ and $g'$. As explained in the introduction, we first have to
  // find the set of cells common to the meshes of the two energy
  // groups. First we call <code>get_finest_common_cells</code> to obtain this
  // list of pairs of common cells from both meshes. Both cells in a pair may
  // not be active but at least one of them is. We then hand each of these
  // cell pairs off to a function that computes the right hand side terms
  // recursively.
  //
  // Note that ingroup coupling is handled already before, so we exit the
  // function early if $g=g'$.
  template <int dim>
  void EnergyGroup<dim>::assemble_cross_group_rhs (const EnergyGroup<dim> &g_prime)
  {
    if (group == g_prime.group)
      return;

    const std::list<std::pair<typename DoFHandler<dim>::cell_iterator,
          typename DoFHandler<dim>::cell_iterator> >
          cell_list
          = GridTools::get_finest_common_cells (dof_handler,
                                                g_prime.dof_handler);

    typename std::list<std::pair<typename DoFHandler<dim>::cell_iterator,
             typename DoFHandler<dim>::cell_iterator> >
             ::const_iterator
             cell_iter = cell_list.begin();

    for (; cell_iter!=cell_list.end(); ++cell_iter)
      {
        FullMatrix<double> unit_matrix (fe.dofs_per_cell);
        for (unsigned int i=0; i<unit_matrix.m(); ++i)
          unit_matrix(i,i) = 1;
        assemble_cross_group_rhs_recursive (g_prime,
                                            cell_iter->first,
                                            cell_iter->second,
                                            unit_matrix);
      }
  }



  // @sect5{<code>EnergyGroup::assemble_cross_group_rhs_recursive</code>}
  //
  // This is finally the function that handles assembling right hand side
  // terms on potentially different meshes recursively, using the algorithm
  // described in the introduction. The function takes a reference to the
  // object representing energy group $g'$, as well as iterators to
  // corresponding cells in the meshes for energy groups $g$ and $g'$. At
  // first, i.e. when this function is called from the one above, these two
  // cells will be matching cells on two meshes; however, one of the two may
  // be further refined, and we will call the function recursively with one of
  // the two iterators replaced by one of the children of the original cell.
  //
  // The last argument is the matrix product matrix $B_{c^{(k)}}^T \cdots
  // B_{c'}^T B_c^T$ from the introduction that interpolates from the coarser
  // of the two cells to the finer one. If the two cells match, then this is
  // the identity matrix -- exactly what we pass to this function initially.
  //
  // The function has to consider two cases: that both of the two cells are
  // not further refined, i.e. have no children, in which case we can finally
  // assemble the right hand side contributions of this pair of cells; and
  // that one of the two cells is further refined, in which case we have to
  // keep recursing by looping over the children of the one cell that is not
  // active. These two cases will be discussed below:
  template <int dim>
  void
  EnergyGroup<dim>::
  assemble_cross_group_rhs_recursive (const EnergyGroup<dim>                        &g_prime,
                                      const typename DoFHandler<dim>::cell_iterator &cell_g,
                                      const typename DoFHandler<dim>::cell_iterator &cell_g_prime,
                                      const FullMatrix<double>                       prolongation_matrix)
  {
    // The first case is that both cells are no further refined. In that case,
    // we can assemble the relevant terms (see the introduction). This
    // involves assembling the mass matrix on the finer of the two cells (in
    // fact there are two mass matrices with different coefficients, one for
    // the fission distribution cross section $\chi_g\nu\Sigma_{f,g'}$ and one
    // for the scattering cross section $\Sigma_{s,g'\to g}$). This is
    // straight forward, but note how we determine which of the two cells is
    // the finer one by looking at the refinement level of the two cells:
    if (!cell_g->has_children() && !cell_g_prime->has_children())
      {
        const QGauss<dim>  quadrature_formula (fe.degree+1);
        const unsigned int n_q_points = quadrature_formula.size();

        FEValues<dim> fe_values (fe, quadrature_formula,
                                 update_values  |  update_JxW_values);

        if (cell_g->level() > cell_g_prime->level())
          fe_values.reinit (cell_g);
        else
          fe_values.reinit (cell_g_prime);

        const double fission_dist_XS
          = material_data.get_fission_dist_XS (group, g_prime.group,
                                               cell_g_prime->material_id());

        const double scattering_XS
          = material_data.get_scattering_XS (g_prime.group, group,
                                             cell_g_prime->material_id());

        FullMatrix<double>    local_mass_matrix_f (fe.dofs_per_cell,
                                                   fe.dofs_per_cell);
        FullMatrix<double>    local_mass_matrix_g (fe.dofs_per_cell,
                                                   fe.dofs_per_cell);

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
            for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
              {
                local_mass_matrix_f(i,j) += (fission_dist_XS *
                                             fe_values.shape_value(i,q_point) *
                                             fe_values.shape_value(j,q_point) *
                                             fe_values.JxW(q_point));
                local_mass_matrix_g(i,j) += (scattering_XS *
                                             fe_values.shape_value(i,q_point) *
                                             fe_values.shape_value(j,q_point) *
                                             fe_values.JxW(q_point));
              }

        // Now we have all the interpolation (prolongation) matrices as well
        // as local mass matrices, so we only have to form the product @f[
        // F_i|_{K_{cc'\cdots c^{(k)}}} = [B_c B_{c'} \cdots B_{c^{(k)}}
        // M_{K_{cc'\cdots c^{(k)}}}]^{ij} \phi_{g'}^j, @f] or @f[
        // F_i|_{K_{cc'\cdots c^{(k)}}} = [(B_c B_{c'} \cdots B_{c^{(k)}}
        // M_{K_{cc'\cdots c^{(k)}}})^T]^{ij} \phi_{g'}^j, @f] depending on
        // which of the two cells is the finer. We do this using either the
        // matrix-vector product provided by the <code>vmult</code> function,
        // or the product with the transpose matrix using <code>Tvmult</code>.
        // After doing so, we transfer the result into the global right hand
        // side vector of energy group $g$.
        Vector<double>       g_prime_new_values (fe.dofs_per_cell);
        Vector<double>       g_prime_old_values (fe.dofs_per_cell);
        cell_g_prime->get_dof_values (g_prime.solution_old, g_prime_old_values);
        cell_g_prime->get_dof_values (g_prime.solution,     g_prime_new_values);

        Vector<double>       cell_rhs (fe.dofs_per_cell);
        Vector<double>       tmp (fe.dofs_per_cell);

        if (cell_g->level() > cell_g_prime->level())
          {
            prolongation_matrix.vmult (tmp, g_prime_old_values);
            local_mass_matrix_f.vmult (cell_rhs, tmp);

            prolongation_matrix.vmult (tmp, g_prime_new_values);
            local_mass_matrix_g.vmult_add (cell_rhs, tmp);
          }
        else
          {
            local_mass_matrix_f.vmult (tmp, g_prime_old_values);
            prolongation_matrix.Tvmult (cell_rhs, tmp);

            local_mass_matrix_g.vmult (tmp, g_prime_new_values);
            prolongation_matrix.Tvmult_add (cell_rhs, tmp);
          }

        std::vector<types::global_dof_index> local_dof_indices (fe.dofs_per_cell);
        cell_g->get_dof_indices (local_dof_indices);

        for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
          system_rhs(local_dof_indices[i]) += cell_rhs(i);
      }

    // The alternative is that one of the two cells is further refined. In
    // that case, we have to loop over all the children, multiply the existing
    // interpolation (prolongation) product of matrices from the left with the
    // interpolation from the present cell to its child (using the
    // matrix-matrix multiplication function <code>mmult</code>), and then
    // hand the result off to this very same function again, but with the cell
    // that has children replaced by one of its children:
    else
      for (unsigned int child=0; child<GeometryInfo<dim>::max_children_per_cell; ++child)
        {
          FullMatrix<double>   new_matrix (fe.dofs_per_cell, fe.dofs_per_cell);
          fe.get_prolongation_matrix(child).mmult (new_matrix,
                                                   prolongation_matrix);

          if (cell_g->has_children())
            assemble_cross_group_rhs_recursive (g_prime,
                                                cell_g->child(child), cell_g_prime,
                                                new_matrix);
          else
            assemble_cross_group_rhs_recursive (g_prime,
                                                cell_g, cell_g_prime->child(child),
                                                new_matrix);
        }
  }


  // @sect5{<code>EnergyGroup::get_fission_source</code>}
  //
  // In the (inverse) power iteration, we use the integrated fission source to
  // update the $k$-eigenvalue. Given its definition, the following function
  // is essentially self-explanatory:
  template <int dim>
  double EnergyGroup<dim>::get_fission_source () const
  {
    const QGauss<dim>  quadrature_formula (fe.degree + 1);
    const unsigned int n_q_points    = quadrature_formula.size();

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values  |  update_JxW_values);

    std::vector<double>       solution_values (n_q_points);

    double fission_source = 0;

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        fe_values.reinit (cell);

        const double fission_XS
          = material_data.get_fission_XS(group, cell->material_id());

        fe_values.get_function_values (solution, solution_values);

        for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
          fission_source += (fission_XS *
                             solution_values[q_point] *
                             fe_values.JxW(q_point));
      }

    return fission_source;
  }


  // @sect5{<code>EnergyGroup::solve</code>}
  //
  // Next a function that solves the linear system assembled before. Things
  // are pretty much standard, except that we delayed applying boundary values
  // until we get here, since in all the previous functions we were still
  // adding up contributions the right hand side vector.
  template <int dim>
  void
  EnergyGroup<dim>::solve ()
  {
    hanging_node_constraints.condense (system_rhs);
    MatrixTools::apply_boundary_values (boundary_values,
                                        system_matrix,
                                        solution,
                                        system_rhs);

    SolverControl           solver_control (system_matrix.m(),
                                            1e-12*system_rhs.l2_norm());
    SolverCG<>              cg (solver_control);

    PreconditionSSOR<> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve (system_matrix, solution, system_rhs, preconditioner);

    hanging_node_constraints.distribute (solution);
  }



  // @sect5{<code>EnergyGroup::estimate_errors</code>}
  //
  // Mesh refinement is split into two functions. The first estimates the
  // error for each cell, normalizes it by the magnitude of the solution, and
  // returns it in the vector given as an argument. The calling function
  // collects all error indicators from all energy groups, and computes
  // thresholds for refining and coarsening cells.
  template <int dim>
  void EnergyGroup<dim>::estimate_errors (Vector<float> &error_indicators) const
  {
    KellyErrorEstimator<dim>::estimate (dof_handler,
                                        QGauss<dim-1> (fe.degree + 1),
                                        typename FunctionMap<dim>::type(),
                                        solution,
                                        error_indicators);
    error_indicators /= solution.linfty_norm();
  }



  // @sect5{<code>EnergyGroup::refine_grid</code>}
  //
  // The second part is to refine the grid given the error indicators compute
  // in the previous function and error thresholds above which cells shall be
  // refined or below which cells shall be coarsened. Note that we do not use
  // any of the functions in <code>GridRefinement</code> here, but rather set
  // refinement flags ourselves.
  //
  // After setting these flags, we use the SolutionTransfer class to move the
  // solution vector from the old to the new mesh. The procedure used here is
  // described in detail in the documentation of that class:
  template <int dim>
  void EnergyGroup<dim>::refine_grid (const Vector<float> &error_indicators,
                                      const double         refine_threshold,
                                      const double         coarsen_threshold)
  {
    typename Triangulation<dim>::active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();

    for (unsigned int cell_index=0; cell!=endc; ++cell, ++cell_index)
      if (error_indicators(cell_index) > refine_threshold)
        cell->set_refine_flag ();
      else if (error_indicators(cell_index) < coarsen_threshold)
        cell->set_coarsen_flag ();

    SolutionTransfer<dim> soltrans(dof_handler);

    triangulation.prepare_coarsening_and_refinement();
    soltrans.prepare_for_coarsening_and_refinement(solution);

    triangulation.execute_coarsening_and_refinement ();
    dof_handler.distribute_dofs (fe);

    solution.reinit (dof_handler.n_dofs());
    soltrans.interpolate(solution_old, solution);

    solution_old.reinit (dof_handler.n_dofs());
    solution_old = solution;
  }


  // @sect5{<code>EnergyGroup::output_results</code>}
  //
  // The last function of this class outputs meshes and solutions after each
  // mesh iteration. This has been shown many times before. The only thing
  // worth pointing out is the use of the
  // <code>Utilities::int_to_string</code> function to convert an integer into
  // its string representation. The second argument of that function denotes
  // how many digits we shall use -- if this value was larger than one, then
  // the number would be padded by leading zeros.
  template <int dim>
  void
  EnergyGroup<dim>::output_results (const unsigned int cycle) const
  {
    {
      const std::string filename = std::string("grid-") +
                                   Utilities::int_to_string(group,1) +
                                   "." +
                                   Utilities::int_to_string(cycle,1) +
                                   ".eps";
      std::ofstream output (filename.c_str());

      GridOut grid_out;
      grid_out.write_eps (triangulation, output);
    }

    {
      const std::string filename = std::string("solution-") +
                                   Utilities::int_to_string(group,1) +
                                   "." +
                                   Utilities::int_to_string(cycle,1) +
                                   ".gmv";

      DataOut<dim> data_out;

      data_out.attach_dof_handler (dof_handler);
      data_out.add_data_vector (solution, "solution");
      data_out.build_patches ();

      std::ofstream output (filename.c_str());
      data_out.write_gmv (output);
    }
  }



  // @sect3{The <code>NeutronDiffusionProblem</code> class template}

  // This is the main class of the program, not because it implements all the
  // functionality (in fact, most of it is implemented in the
  // <code>EnergyGroup</code> class) but because it contains the driving
  // algorithm that determines what to compute and when. It is mostly as shown
  // in many of the other tutorial programs in that it has a public
  // <code>run</code> function and private functions doing all the rest. In
  // several places, we have to do something for all energy groups, in which
  // case we will start threads for each group to let these things run in
  // parallel if deal.II was configured for multithreading.  For strategies of
  // parallelization, take a look at the @ref threads module.
  //
  // The biggest difference to previous example programs is that we also
  // declare a nested class that has member variables for all the run-time
  // parameters that can be passed to the program in an input file. Right now,
  // these are the number of energy groups, the number of refinement cycles,
  // the polynomial degree of the finite element to be used, and the tolerance
  // used to determine when convergence of the inverse power iteration has
  // occurred. In addition, we have a constructor of this class that sets all
  // these values to their default values, a function
  // <code>declare_parameters</code> that described to the ParameterHandler
  // class already used in step-19 what parameters are accepted in the input
  // file, and a function <code>get_parameters</code> that can extract the
  // values of these parameters from a ParameterHandler object.
  template <int dim>
  class NeutronDiffusionProblem
  {
  public:
    class Parameters
    {
    public:
      Parameters ();

      static void declare_parameters (ParameterHandler &prm);
      void get_parameters (ParameterHandler &prm);

      unsigned int n_groups;
      unsigned int n_refinement_cycles;

      unsigned int fe_degree;

      double convergence_tolerance;
    };



    NeutronDiffusionProblem (const Parameters &parameters);
    ~NeutronDiffusionProblem ();

    void run ();

  private:
    // @sect5{Private member functions}

    // There are not that many member functions in this class since most of
    // the functionality has been moved into the <code>EnergyGroup</code>
    // class and is simply called from the <code>run()</code> member function
    // of this class. The ones that remain have self-explanatory names:
    void initialize_problem();

    void refine_grid ();

    double get_total_fission_source () const;


    // @sect5{Private member variables}

    // Next, we have a few member variables. In particular, these are (i) a
    // reference to the parameter object (owned by the main function of this
    // program, and passed to the constructor of this class), (ii) an object
    // describing the material parameters for the number of energy groups
    // requested in the input file, and (iii) the finite element to be used by
    // all energy groups:
    const Parameters  &parameters;
    const MaterialData material_data;
    FE_Q<dim>          fe;

    // Furthermore, we have (iv) the value of the computed eigenvalue at the
    // present iteration. This is, in fact, the only part of the solution that
    // is shared between all energy groups -- all other parts of the solution,
    // such as neutron fluxes are particular to one or the other energy group,
    // and are therefore stored in objects that describe a single energy
    // group:
    double k_eff;

    // Finally, (v), we have an array of pointers to the energy group
    // objects. The length of this array is, of course, equal to the number of
    // energy groups specified in the parameter file.
    std::vector<EnergyGroup<dim>*> energy_groups;
  };


  // @sect4{Implementation of the <code>NeutronDiffusionProblem::Parameters</code> class}

  // Before going on to the implementation of the outer class, we have to
  // implement the functions of the parameters structure. This is pretty
  // straightforward and, in fact, looks pretty much the same for all such
  // parameters classes using the ParameterHandler capabilities. We will
  // therefore not comment further on this:
  template <int dim>
  NeutronDiffusionProblem<dim>::Parameters::Parameters ()
    :
    n_groups (2),
    n_refinement_cycles (5),
    fe_degree (2),
    convergence_tolerance (1e-12)
  {}



  template <int dim>
  void
  NeutronDiffusionProblem<dim>::Parameters::
  declare_parameters (ParameterHandler &prm)
  {
    prm.declare_entry ("Number of energy groups", "2",
                       Patterns::Integer (),
                       "The number of energy different groups considered");
    prm.declare_entry ("Refinement cycles", "5",
                       Patterns::Integer (),
                       "Number of refinement cycles to be performed");
    prm.declare_entry ("Finite element degree", "2",
                       Patterns::Integer (),
                       "Polynomial degree of the finite element to be used");
    prm.declare_entry ("Power iteration tolerance", "1e-12",
                       Patterns::Double (),
                       "Inner power iterations are stopped when the change in k_eff falls "
                       "below this tolerance");
  }



  template <int dim>
  void
  NeutronDiffusionProblem<dim>::Parameters::
  get_parameters (ParameterHandler &prm)
  {
    n_groups              = prm.get_integer ("Number of energy groups");
    n_refinement_cycles   = prm.get_integer ("Refinement cycles");
    fe_degree             = prm.get_integer ("Finite element degree");
    convergence_tolerance = prm.get_double ("Power iteration tolerance");
  }




  // @sect4{Implementation of the <code>NeutronDiffusionProblem</code> class}

  // Now for the <code>NeutronDiffusionProblem</code> class. The constructor
  // and destructor have nothing of much interest:
  template <int dim>
  NeutronDiffusionProblem<dim>::
  NeutronDiffusionProblem (const Parameters &parameters)
    :
    parameters (parameters),
    material_data (parameters.n_groups),
    fe (parameters.fe_degree)
  {}



  template <int dim>
  NeutronDiffusionProblem<dim>::~NeutronDiffusionProblem ()
  {
    for (unsigned int group=0; group<energy_groups.size(); ++group)
      delete energy_groups[group];

    energy_groups.resize (0);
  }

  // @sect5{<code>NeutronDiffusionProblem::initialize_problem</code>}
  //
  // The first function of interest is the one that sets up the geometry of
  // the reactor core. This is described in more detail in the introduction.
  //
  // The first part of the function defines geometry data, and then creates a
  // coarse mesh that has as many cells as there are fuel rods (or pin cells,
  // for that matter) in that part of the reactor core that we simulate. As
  // mentioned when interpolating boundary values above, the last parameter to
  // the <code>GridGenerator::subdivided_hyper_rectangle</code> function
  // specifies that sides of the domain shall have unique boundary indicators
  // that will later allow us to determine in a simple way which of the
  // boundaries have Neumann and which have Dirichlet conditions attached to
  // them.
  template <int dim>
  void NeutronDiffusionProblem<dim>::initialize_problem()
  {
    const unsigned int rods_per_assembly_x = 17,
                       rods_per_assembly_y = 17;
    const double pin_pitch_x = 1.26,
                 pin_pitch_y = 1.26;
    const double assembly_height = 200;

    const unsigned int assemblies_x = 2,
                       assemblies_y = 2,
                       assemblies_z = 1;

    const Point<dim> bottom_left = Point<dim>();
    const Point<dim> upper_right = (dim == 2
                                    ?
                                    Point<dim> (assemblies_x*rods_per_assembly_x*pin_pitch_x,
                                                assemblies_y*rods_per_assembly_y*pin_pitch_y)
                                    :
                                    Point<dim> (assemblies_x*rods_per_assembly_x*pin_pitch_x,
                                                assemblies_y*rods_per_assembly_y*pin_pitch_y,
                                                assemblies_z*assembly_height));

    std::vector<unsigned int> n_subdivisions;
    n_subdivisions.push_back (assemblies_x*rods_per_assembly_x);
    if (dim >= 2)
      n_subdivisions.push_back (assemblies_y*rods_per_assembly_y);
    if (dim >= 3)
      n_subdivisions.push_back (assemblies_z);

    Triangulation<dim> coarse_grid;
    GridGenerator::subdivided_hyper_rectangle (coarse_grid,
                                               n_subdivisions,
                                               bottom_left,
                                               upper_right,
                                               true);


    // The second part of the function deals with material numbers of pin
    // cells of each type of assembly. Here, we define four different types of
    // assembly, for which we describe the arrangement of fuel rods in the
    // following tables.
    //
    // The assemblies described here are taken from the benchmark mentioned in
    // the introduction and are (in this order): <ol> <li>'UX' Assembly: UO2
    // fuel assembly with 24 guide tubes and a central Moveable Fission
    // Chamber <li>'UA' Assembly: UO2 fuel assembly with 24 AIC and a central
    // Moveable Fission Chamber <li>'PX' Assembly: MOX fuel assembly with 24
    // guide tubes and a central Moveable Fission Chamber <li>'R' Assembly: a
    // reflector.  </ol>
    //
    // Note that the numbers listed here and taken from the benchmark
    // description are, in good old Fortran fashion, one-based. We will later
    // subtract one from each number when assigning materials to individual
    // cells to convert things into the C-style zero-based indexing.
    const unsigned int n_assemblies=4;
    const unsigned int
    assembly_materials[n_assemblies][rods_per_assembly_x][rods_per_assembly_y]
    =
    {
      {
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 5, 1, 1, 5, 1, 1, 7, 1, 1, 5, 1, 1, 5, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 5, 1, 1, 1, 1, 1, 1, 1, 1, 1, 5, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 5, 1, 1, 5, 1, 1, 5, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }
      },
      {
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 8, 1, 1, 8, 1, 1, 7, 1, 1, 8, 1, 1, 8, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 8, 1, 1, 1, 1, 1, 1, 1, 1, 1, 8, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 8, 1, 1, 8, 1, 1, 8, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 },
        { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 }
      },
      {
        { 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 },
        { 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2 },
        { 2, 3, 3, 3, 3, 5, 3, 3, 5, 3, 3, 5, 3, 3, 3, 3, 2 },
        { 2, 3, 3, 5, 3, 4, 4, 4, 4, 4, 4, 4, 3, 5, 3, 3, 2 },
        { 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 2 },
        { 2, 3, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 3, 2 },
        { 2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2 },
        { 2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2 },
        { 2, 3, 5, 4, 4, 5, 4, 4, 7, 4, 4, 5, 4, 4, 5, 3, 2 },
        { 2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2 },
        { 2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 2 },
        { 2, 3, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 4, 4, 5, 3, 2 },
        { 2, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 2 },
        { 2, 3, 3, 5, 3, 4, 4, 4, 4, 4, 4, 4, 3, 5, 3, 3, 2 },
        { 2, 3, 3, 3, 3, 5, 3, 3, 5, 3, 3, 5, 3, 3, 3, 3, 2 },
        { 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2 },
        { 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 }
      },
      {
        { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 },
        { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 },
        { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 },
        { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 },
        { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 },
        { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 },
        { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 },
        { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 },
        { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 },
        { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 },
        { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 },
        { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 },
        { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 },
        { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 },
        { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 },
        { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 },
        { 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6 }
      }
    };

    // After the description of the materials that make up an assembly, we
    // have to specify the arrangement of assemblies within the core. We use a
    // symmetric pattern that in fact only uses the 'UX' and 'PX' assemblies:
    const unsigned int core[assemblies_x][assemblies_y][assemblies_z]
    =  {{{0}, {2}}, {{2}, {0}}};

    // We are now in a position to actually set material IDs for each cell. To
    // this end, we loop over all cells, look at the location of the cell's
    // center, and determine which assembly and fuel rod this would be in. (We
    // add a few checks to see that the locations we compute are within the
    // bounds of the arrays in which we have to look up materials.) At the end
    // of the loop, we set material identifiers accordingly:
    for (typename Triangulation<dim>::active_cell_iterator
         cell = coarse_grid.begin_active();
         cell!=coarse_grid.end();
         ++cell)
      {
        const Point<dim> cell_center = cell->center();

        const unsigned int tmp_x = int(cell_center[0]/pin_pitch_x);
        const unsigned int ax = tmp_x/rods_per_assembly_x;
        const unsigned int cx = tmp_x - ax * rods_per_assembly_x;

        const unsigned tmp_y = int(cell_center[1]/pin_pitch_y);
        const unsigned int ay = tmp_y/rods_per_assembly_y;
        const unsigned int cy = tmp_y - ay * rods_per_assembly_y;

        const unsigned int az = (dim == 2
                                 ?
                                 0
                                 :
                                 int (cell_center[dim-1]/assembly_height));

        Assert (ax < assemblies_x, ExcInternalError());
        Assert (ay < assemblies_y, ExcInternalError());
        Assert (az < assemblies_z, ExcInternalError());

        Assert (core[ax][ay][az] < n_assemblies, ExcInternalError());

        Assert (cx < rods_per_assembly_x, ExcInternalError());
        Assert (cy < rods_per_assembly_y, ExcInternalError());

        cell->set_material_id(assembly_materials[core[ax][ay][az]][cx][cy] - 1);
      }

    // With the coarse mesh so initialized, we create the appropriate number
    // of energy group objects and let them initialize their individual meshes
    // with the coarse mesh generated above:
    energy_groups.resize (parameters.n_groups);
    for (unsigned int group=0; group<parameters.n_groups; ++group)
      energy_groups[group] = new EnergyGroup<dim> (group, material_data,
                                                   coarse_grid, fe);
  }


  // @sect5{<code>NeutronDiffusionProblem::get_total_fission_source</code>}
  //
  // In the eigenvalue computation, we need to calculate total fission neutron
  // source after each power iteration. The total power then is used to renew
  // k-effective.
  //
  // Since the total fission source is a sum over all the energy groups, and
  // since each of these sums can be computed independently, we actually do
  // this in parallel. One of the problems is that the function in the
  // <code>EnergyGroup</code> class that computes the fission source returns a
  // value. If we now simply spin off a new thread, we have to later capture
  // the return value of the function run on that thread. The way this can be
  // done is to use the return value of the Threads::new_thread function,
  // which returns an object of type Threads::Thread@<double@> if the function
  // spawned returns a double. We can then later ask this object for the
  // returned value (when doing so, the Threads::Thread::return_value function
  // first waits for the thread to finish if it hasn't done so already).
  //
  // The way this function then works is to first spawn one thread for each
  // energy group we work with, then one-by-one collecting the returned values
  // of each thread and return the sum.
  template <int dim>
  double NeutronDiffusionProblem<dim>::get_total_fission_source () const
  {
    std::vector<Threads::Thread<double> > threads;
    for (unsigned int group=0; group<parameters.n_groups; ++group)
      threads.push_back (Threads::new_thread (&EnergyGroup<dim>::get_fission_source,
                                              *energy_groups[group]));

    double fission_source = 0;
    for (unsigned int group=0; group<parameters.n_groups; ++group)
      fission_source += threads[group].return_value ();

    return fission_source;
  }




  // @sect5{<code>NeutronDiffusionProblem::refine_grid</code>}
  //
  // The next function lets the individual energy group objects refine their
  // meshes. Much of this, again, is a task that can be done independently in
  // parallel: first, let all the energy group objects calculate their error
  // indicators in parallel, then compute the maximum error indicator over all
  // energy groups and determine thresholds for refinement and coarsening of
  // cells, and then ask all the energy groups to refine their meshes
  // accordingly, again in parallel.
  template <int dim>
  void NeutronDiffusionProblem<dim>::refine_grid ()
  {
    std::vector<types::global_dof_index> n_cells (parameters.n_groups);
    for (unsigned int group=0; group<parameters.n_groups; ++group)
      n_cells[group] = energy_groups[group]->n_active_cells();

    BlockVector<float>  group_error_indicators(n_cells);

    {
      Threads::ThreadGroup<> threads;
      for (unsigned int group=0; group<parameters.n_groups; ++group)
        threads += Threads::new_thread (&EnergyGroup<dim>::estimate_errors,
                                        *energy_groups[group],
                                        group_error_indicators.block(group));
      threads.join_all ();
    }

    const float max_error         = group_error_indicators.linfty_norm();
    const float refine_threshold  = 0.3*max_error;
    const float coarsen_threshold = 0.01*max_error;

    {
      Threads::ThreadGroup<> threads;
      for (unsigned int group=0; group<parameters.n_groups; ++group)
        threads += Threads::new_thread (&EnergyGroup<dim>::refine_grid,
                                        *energy_groups[group],
                                        group_error_indicators.block(group),
                                        refine_threshold,
                                        coarsen_threshold);
      threads.join_all ();
    }
  }


  // @sect5{<code>NeutronDiffusionProblem::run</code>}
  //
  // Finally, this is the function where the meat is: iterate on a sequence of
  // meshes, and on each of them do a power iteration to compute the
  // eigenvalue.
  //
  // Given the description of the algorithm in the introduction, there is
  // actually not much to comment on:
  template <int dim>
  void NeutronDiffusionProblem<dim>::run ()
  {
    std::cout << std::setprecision (12) << std::fixed;

    double k_eff_old = k_eff;

    Timer timer;
    timer.start ();

    for (unsigned int cycle=0; cycle<parameters.n_refinement_cycles; ++cycle)
      {
        std::cout << "Cycle " << cycle << ':' << std::endl;

        if (cycle == 0)
          initialize_problem();
        else
          {
            refine_grid ();
            for (unsigned int group=0; group<parameters.n_groups; ++group)
              energy_groups[group]->solution *= k_eff;
          }

        for (unsigned int group=0; group<parameters.n_groups; ++group)
          energy_groups[group]->setup_linear_system ();

        std::cout << "   Numbers of active cells:       ";
        for (unsigned int group=0; group<parameters.n_groups; ++group)
          std::cout << energy_groups[group]->n_active_cells()
                    << ' ';
        std::cout << std::endl;
        std::cout << "   Numbers of degrees of freedom: ";
        for (unsigned int group=0; group<parameters.n_groups; ++group)
          std::cout << energy_groups[group]->n_dofs()
                    << ' ';
        std::cout << std::endl << std::endl;


        Threads::ThreadGroup<> threads;
        for (unsigned int group=0; group<parameters.n_groups; ++group)
          threads += Threads::new_thread
                     (&EnergyGroup<dim>::assemble_system_matrix,
                      *energy_groups[group]);
        threads.join_all ();

        double error;
        unsigned int iteration = 1;
        do
          {
            for (unsigned int group=0; group<parameters.n_groups; ++group)
              {
                energy_groups[group]->assemble_ingroup_rhs (ZeroFunction<dim>());

                for (unsigned int bgroup=0; bgroup<parameters.n_groups; ++bgroup)
                  energy_groups[group]->assemble_cross_group_rhs (*energy_groups[bgroup]);

                energy_groups[group]->solve ();
              }

            k_eff = get_total_fission_source();
            error = fabs(k_eff-k_eff_old)/fabs(k_eff);
            std::cout << "   Iteration " << iteration
                      << ": k_eff=" << k_eff
                      << std::endl;
            k_eff_old=k_eff;

            for (unsigned int group=0; group<parameters.n_groups; ++group)
              {
                energy_groups[group]->solution_old = energy_groups[group]->solution;
                energy_groups[group]->solution_old /= k_eff;
              }

            ++iteration;
          }
        while ((error > parameters.convergence_tolerance)
               &&
               (iteration < 500));

        for (unsigned int group=0; group<parameters.n_groups; ++group)
          energy_groups[group]->output_results (cycle);

        std::cout << std::endl;
        std::cout << "   Cycle=" << cycle
                  << ", n_dofs=" << energy_groups[0]->n_dofs() + energy_groups[1]->n_dofs()
                  << ",  k_eff=" << k_eff
                  << ", time=" << timer()
                  << std::endl;


        std::cout << std::endl << std::endl;
      }
  }
}



// @sect3{The <code>main()</code> function}
//
// The last thing in the program in the <code>main()</code> function. The
// structure is as in most other tutorial programs, with the only exception
// that we here handle a parameter file.  To this end, we first look at the
// command line arguments passed to this function: if no input file is
// specified on the command line, then use "project.prm", otherwise take the
// filename given as the first argument on the command line.
//
// With this, we create a ParameterHandler object, let the
// <code>NeutronDiffusionProblem::Parameters</code> class declare all the
// parameters it wants to see in the input file (or, take the default values,
// if nothing is listed in the parameter file), then read the input file, ask
// the parameters object to extract the values, and finally hand everything
// off to an object of type <code>NeutronDiffusionProblem</code> for
// computation of the eigenvalue:
int main (int argc, char **argv)
{
  try
    {
      using namespace dealii;
      using namespace Step28;

      deallog.depth_console (0);

      std::string filename;
      if (argc < 2)
        filename = "project.prm";
      else
        filename = argv[1];


      const unsigned int dim = 2;

      ParameterHandler parameter_handler;

      NeutronDiffusionProblem<dim>::Parameters parameters;
      parameters.declare_parameters (parameter_handler);

      parameter_handler.read_input (filename);

      parameters.get_parameters (parameter_handler);


      NeutronDiffusionProblem<dim> neutron_diffusion_problem (parameters);
      neutron_diffusion_problem.run ();
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
