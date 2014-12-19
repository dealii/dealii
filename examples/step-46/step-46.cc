/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2011 - 2013 by the deal.II authors
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
 * Author: Wolfgang Bangerth, Texas A&M University, 2011
 */


// @sect3{Include files}

// The include files for this program are the same as for many others
// before. The only new one is the one that declares FE_Nothing as discussed
// in the introduction. The ones in the hp directory have already been
// discussed in step-27.

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <fstream>
#include <sstream>


namespace Step46
{
  using namespace dealii;

  // @sect3{The <code>FluidStructureProblem</code> class template}

  // This is the main class. It is, if you want, a combination of step-8 and
  // step-22 in that it has member variables that either address the global
  // problem (the Triangulation and hp::DoFHandler objects, as well as the
  // hp::FECollection and various linear algebra objects) or that pertain to
  // either the elasticity or Stokes sub-problems. The general structure of
  // the class, however, is like that of most of the other programs
  // implementing stationary problems.
  //
  // There are a few helper functions (<code>cell_is_in_fluid_domain,
  // cell_is_in_solid_domain</code>) of self-explanatory nature (operating on
  // the symbolic names for the two subdomains that will be used as
  // material_ids for cells belonging to the subdomains, as explained in the
  // introduction) and a few functions (<code>make_grid,
  // set_active_fe_indices, assemble_interface_terms</code>) that have been
  // broken out of other functions that can be found in many of the other
  // tutorial programs and that will be discussed as we get to their
  // implementation.
  //
  // The final set of variables (<code>viscosity, lambda, eta</code>)
  // describes the material properties used for the two physics models.
  template <int dim>
  class FluidStructureProblem
  {
  public:
    FluidStructureProblem (const unsigned int stokes_degree,
                           const unsigned int elasticity_degree);
    void run ();

  private:
    enum
    {
      fluid_domain_id,
      solid_domain_id
    };

    static bool
    cell_is_in_fluid_domain (const typename hp::DoFHandler<dim>::cell_iterator &cell);

    static bool
    cell_is_in_solid_domain (const typename hp::DoFHandler<dim>::cell_iterator &cell);


    void make_grid ();
    void set_active_fe_indices ();
    void setup_dofs ();
    void assemble_system ();
    void assemble_interface_term (const FEFaceValuesBase<dim>          &elasticity_fe_face_values,
                                  const FEFaceValuesBase<dim>          &stokes_fe_face_values,
                                  std::vector<Tensor<1,dim> >          &elasticity_phi,
                                  std::vector<SymmetricTensor<2,dim> > &stokes_symgrad_phi_u,
                                  std::vector<double>                  &stokes_phi_p,
                                  FullMatrix<double>                   &local_interface_matrix) const;
    void solve ();
    void output_results (const unsigned int refinement_cycle) const;
    void refine_mesh ();

    const unsigned int    stokes_degree;
    const unsigned int    elasticity_degree;

    Triangulation<dim>    triangulation;
    FESystem<dim>         stokes_fe;
    FESystem<dim>         elasticity_fe;
    hp::FECollection<dim> fe_collection;
    hp::DoFHandler<dim>   dof_handler;

    ConstraintMatrix      constraints;

    SparsityPattern       sparsity_pattern;
    SparseMatrix<double>  system_matrix;

    Vector<double>        solution;
    Vector<double>        system_rhs;

    const double          viscosity;
    const double          lambda;
    const double          mu;
  };


  // @sect3{Boundary values and right hand side}

  // The following classes do as their names suggest. The boundary values for
  // the velocity are $\mathbf u=(0, \sin(\pi x))^T$ in 2d and $\mathbf u=(0,
  // 0, \sin(\pi x)\sin(\pi y))^T$ in 3d, respectively. The remaining boundary
  // conditions for this problem are all homogeneous and have been discussed in
  // the introduction. The right hand side forcing term is zero for both the
  // fluid and the solid.
  template <int dim>
  class StokesBoundaryValues : public Function<dim>
  {
  public:
    StokesBoundaryValues () : Function<dim>(dim+1+dim) {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;
  };


  template <int dim>
  double
  StokesBoundaryValues<dim>::value (const Point<dim>  &p,
                                    const unsigned int component) const
  {
    Assert (component < this->n_components,
            ExcIndexRange (component, 0, this->n_components));

    if (component == dim-1)
      switch (dim)
        {
        case 2:
          return std::sin(numbers::PI*p[0]);
        case 3:
          return std::sin(numbers::PI*p[0]) * std::sin(numbers::PI*p[1]);
        default:
          Assert (false, ExcNotImplemented());
        }

    return 0;
  }


  template <int dim>
  void
  StokesBoundaryValues<dim>::vector_value (const Point<dim> &p,
                                           Vector<double>   &values) const
  {
    for (unsigned int c=0; c<this->n_components; ++c)
      values(c) = StokesBoundaryValues<dim>::value (p, c);
  }



  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    RightHandSide () : Function<dim>(dim+1) {}

    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p,
                               Vector<double>   &value) const;

  };


  template <int dim>
  double
  RightHandSide<dim>::value (const Point<dim>  &/*p*/,
                             const unsigned int /*component*/) const
  {
    return 0;
  }


  template <int dim>
  void
  RightHandSide<dim>::vector_value (const Point<dim> &p,
                                    Vector<double>   &values) const
  {
    for (unsigned int c=0; c<this->n_components; ++c)
      values(c) = RightHandSide<dim>::value (p, c);
  }



  // @sect3{The <code>FluidStructureProblem</code> implementation}

  // @sect4{Constructors and helper functions}

  // Let's now get to the implementation of the primary class of this
  // program. The first few functions are the constructor and the helper
  // functions that can be used to determine which part of the domain a cell
  // is in. Given the discussion of these topics in the introduction, their
  // implementation is rather obvious. In the constructor, note that we have
  // to construct the hp::FECollection object from the base elements for
  // Stokes and elasticity; using the hp::FECollection::push_back function
  // assigns them spots zero and one in this collection, an order that we have
  // to remember and use consistently in the rest of the program.
  template <int dim>
  FluidStructureProblem<dim>::
  FluidStructureProblem (const unsigned int stokes_degree,
                         const unsigned int elasticity_degree)
    :
    stokes_degree (stokes_degree),
    elasticity_degree (elasticity_degree),
    triangulation (Triangulation<dim>::maximum_smoothing),
    stokes_fe (FE_Q<dim>(stokes_degree+1), dim,
               FE_Q<dim>(stokes_degree), 1,
               FE_Nothing<dim>(), dim),
    elasticity_fe (FE_Nothing<dim>(), dim,
                   FE_Nothing<dim>(), 1,
                   FE_Q<dim>(elasticity_degree), dim),
    dof_handler (triangulation),
    viscosity (2),
    lambda (1),
    mu (1)
  {
    fe_collection.push_back (stokes_fe);
    fe_collection.push_back (elasticity_fe);
  }




  template <int dim>
  bool
  FluidStructureProblem<dim>::
  cell_is_in_fluid_domain (const typename hp::DoFHandler<dim>::cell_iterator &cell)
  {
    return (cell->material_id() == fluid_domain_id);
  }


  template <int dim>
  bool
  FluidStructureProblem<dim>::
  cell_is_in_solid_domain (const typename hp::DoFHandler<dim>::cell_iterator &cell)
  {
    return (cell->material_id() == solid_domain_id);
  }


  // @sect4{Meshes and assigning subdomains}

  // The next pair of functions deals with generating a mesh and making sure
  // all flags that denote subdomains are correct. <code>make_grid</code>, as
  // discussed in the introduction, generates an $8\times 8$ mesh (or an
  // $8\times 8\times 8$ mesh in 3d) to make sure that each coarse mesh cell
  // is completely within one of the subdomains. After generating this mesh,
  // we loop over its boundary and set the boundary indicator to one at the
  // top boundary, the only place where we set nonzero Dirichlet boundary
  // conditions. After this, we loop again over all cells to set the material
  // indicator &mdash; used to denote which part of the domain we are in, to
  // either the fluid or solid indicator.
  template <int dim>
  void
  FluidStructureProblem<dim>::make_grid ()
  {
    GridGenerator::subdivided_hyper_cube (triangulation, 8, -1, 1);

    for (typename Triangulation<dim>::active_cell_iterator
         cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell)
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        if (cell->face(f)->at_boundary()
            &&
            (cell->face(f)->center()[dim-1] == 1))
          cell->face(f)->set_all_boundary_indicators(1);


    for (typename Triangulation<dim>::active_cell_iterator
         cell = dof_handler.begin_active();
         cell != dof_handler.end(); ++cell)
      if (((std::fabs(cell->center()[0]) < 0.25)
           &&
           (cell->center()[dim-1] > 0.5))
          ||
          ((std::fabs(cell->center()[0]) >= 0.25)
           &&
           (cell->center()[dim-1] > -0.5)))
        cell->set_material_id (fluid_domain_id);
      else
        cell->set_material_id (solid_domain_id);
  }


  // The second part of this pair of functions determines which finite element
  // to use on each cell. Above we have set the material indicator for each
  // coarse mesh cell, and as mentioned in the introduction, this information
  // is inherited from mother to child cell upon mesh refinement.
  //
  // In other words, whenever we have refined (or created) the mesh, we can
  // rely on the material indicators to be a correct description of which part
  // of the domain a cell is in. We then use this to set the active FE index
  // of the cell to the corresponding element of the hp::FECollection member
  // variable of this class: zero for fluid cells, one for solid cells.
  template <int dim>
  void
  FluidStructureProblem<dim>::set_active_fe_indices ()
  {
    for (typename hp::DoFHandler<dim>::active_cell_iterator
         cell = dof_handler.begin_active();
         cell != dof_handler.end(); ++cell)
      {
        if (cell_is_in_fluid_domain(cell))
          cell->set_active_fe_index (0);
        else if (cell_is_in_solid_domain(cell))
          cell->set_active_fe_index (1);
        else
          Assert (false, ExcNotImplemented());
      }
  }


  // @sect4{<code>FluidStructureProblem::setup_dofs</code>}

  // The next step is to setup the data structures for the linear system. To
  // this end, we first have to set the active FE indices with the function
  // immediately above, then distribute degrees of freedom, and then determine
  // constraints on the linear system. The latter includes hanging node
  // constraints as usual, but also the inhomogeneous boundary values at the
  // top fluid boundary, and zero boundary values along the perimeter of the
  // solid subdomain.
  template <int dim>
  void
  FluidStructureProblem<dim>::setup_dofs ()
  {
    set_active_fe_indices ();
    dof_handler.distribute_dofs (fe_collection);

    {
      constraints.clear ();
      DoFTools::make_hanging_node_constraints (dof_handler,
                                               constraints);

      const FEValuesExtractors::Vector velocities(0);
      VectorTools::interpolate_boundary_values (dof_handler,
                                                1,
                                                StokesBoundaryValues<dim>(),
                                                constraints,
                                                fe_collection.component_mask(velocities));

      const FEValuesExtractors::Vector displacements(dim+1);
      VectorTools::interpolate_boundary_values (dof_handler,
                                                0,
                                                ZeroFunction<dim>(dim+1+dim),
                                                constraints,
                                                fe_collection.component_mask(displacements));
    }

    // There are more constraints we have to handle, though: we have to make
    // sure that the velocity is zero at the interface between fluid and
    // solid. The following piece of code was already presented in the
    // introduction:
    {
      std::vector<types::global_dof_index> local_face_dof_indices (stokes_fe.dofs_per_face);
      for (typename hp::DoFHandler<dim>::active_cell_iterator
           cell = dof_handler.begin_active();
           cell != dof_handler.end(); ++cell)
        if (cell_is_in_fluid_domain (cell))
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
            if (!cell->at_boundary(f))
              {
                bool face_is_on_interface = false;

                if ((cell->neighbor(f)->has_children() == false)
                    &&
                    (cell_is_in_solid_domain (cell->neighbor(f))))
                  face_is_on_interface = true;
                else if (cell->neighbor(f)->has_children() == true)
                  {
                    for (unsigned int sf=0; sf<cell->face(f)->n_children(); ++sf)
                      if (cell_is_in_solid_domain (cell->neighbor_child_on_subface
                                                   (f, sf)))
                        {
                          face_is_on_interface = true;
                          break;
                        }
                  }

                if (face_is_on_interface)
                  {
                    cell->face(f)->get_dof_indices (local_face_dof_indices, 0);
                    for (unsigned int i=0; i<local_face_dof_indices.size(); ++i)
                      if (stokes_fe.face_system_to_component_index(i).first < dim)
                        constraints.add_line (local_face_dof_indices[i]);
                  }
              }
    }

    // At the end of all this, we can declare to the constraints object that
    // we now have all constraints ready to go and that the object can rebuild
    // its internal data structures for better efficiency:
    constraints.close ();

    std::cout << "   Number of active cells: "
              << triangulation.n_active_cells()
              << std::endl
              << "   Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << std::endl;

    // In the rest of this function we create a sparsity pattern as discussed
    // extensively in the introduction, and use it to initialize the matrix;
    // then also set vectors to their correct sizes:
    {
      CompressedSimpleSparsityPattern csp (dof_handler.n_dofs(),
                                           dof_handler.n_dofs());

      Table<2,DoFTools::Coupling> cell_coupling (fe_collection.n_components(),
                                                 fe_collection.n_components());
      Table<2,DoFTools::Coupling> face_coupling (fe_collection.n_components(),
                                                 fe_collection.n_components());

      for (unsigned int c=0; c<fe_collection.n_components(); ++c)
        for (unsigned int d=0; d<fe_collection.n_components(); ++d)
          {
            if (((c<dim+1) && (d<dim+1)
                 && !((c==dim) && (d==dim)))
                ||
                ((c>=dim+1) && (d>=dim+1)))
              cell_coupling[c][d] = DoFTools::always;

            if ((c>=dim+1) && (d<dim+1))
              face_coupling[c][d] = DoFTools::always;
          }

      DoFTools::make_flux_sparsity_pattern (dof_handler, csp,
                                            cell_coupling, face_coupling);
      constraints.condense (csp);
      sparsity_pattern.copy_from (csp);
    }

    system_matrix.reinit (sparsity_pattern);

    solution.reinit (dof_handler.n_dofs());
    system_rhs.reinit (dof_handler.n_dofs());
  }



  // @sect4{<code>FluidStructureProblem::assemble_system</code>}

  // Following is the central function of this program: the one that assembles
  // the linear system. It has a long section of setting up auxiliary
  // functions at the beginning: from creating the quadrature formulas and
  // setting up the FEValues, FEFaceValues and FESubfaceValues objects
  // necessary to integrate the cell terms as well as the interface terms for
  // the case where cells along the interface come together at same size or
  // with differing levels of refinement...
  template <int dim>
  void FluidStructureProblem<dim>::assemble_system ()
  {
    system_matrix=0;
    system_rhs=0;

    const QGauss<dim> stokes_quadrature(stokes_degree+2);
    const QGauss<dim> elasticity_quadrature(elasticity_degree+2);

    hp::QCollection<dim>  q_collection;
    q_collection.push_back (stokes_quadrature);
    q_collection.push_back (elasticity_quadrature);

    hp::FEValues<dim> hp_fe_values (fe_collection, q_collection,
                                    update_values    |
                                    update_quadrature_points  |
                                    update_JxW_values |
                                    update_gradients);

    const QGauss<dim-1> common_face_quadrature(std::max (stokes_degree+2,
                                                         elasticity_degree+2));

    FEFaceValues<dim>    stokes_fe_face_values (stokes_fe,
                                                common_face_quadrature,
                                                update_JxW_values |
                                                update_normal_vectors |
                                                update_gradients);
    FEFaceValues<dim>    elasticity_fe_face_values (elasticity_fe,
                                                    common_face_quadrature,
                                                    update_values);
    FESubfaceValues<dim> stokes_fe_subface_values (stokes_fe,
                                                   common_face_quadrature,
                                                   update_JxW_values |
                                                   update_normal_vectors |
                                                   update_gradients);
    FESubfaceValues<dim> elasticity_fe_subface_values (elasticity_fe,
                                                       common_face_quadrature,
                                                       update_values);

    // ...to objects that are needed to describe the local contributions to
    // the global linear system...
    const unsigned int        stokes_dofs_per_cell     = stokes_fe.dofs_per_cell;
    const unsigned int        elasticity_dofs_per_cell = elasticity_fe.dofs_per_cell;

    FullMatrix<double>        local_matrix;
    FullMatrix<double>        local_interface_matrix (elasticity_dofs_per_cell,
                                                      stokes_dofs_per_cell);
    Vector<double>            local_rhs;

    std::vector<types::global_dof_index> local_dof_indices;
    std::vector<types::global_dof_index> neighbor_dof_indices (stokes_dofs_per_cell);

    const RightHandSide<dim>  right_hand_side;

    // ...to variables that allow us to extract certain components of the
    // shape functions and cache their values rather than having to recompute
    // them at every quadrature point:
    const FEValuesExtractors::Vector     velocities (0);
    const FEValuesExtractors::Scalar     pressure (dim);
    const FEValuesExtractors::Vector     displacements (dim+1);

    std::vector<SymmetricTensor<2,dim> > stokes_symgrad_phi_u (stokes_dofs_per_cell);
    std::vector<double>                  stokes_div_phi_u     (stokes_dofs_per_cell);
    std::vector<double>                  stokes_phi_p         (stokes_dofs_per_cell);

    std::vector<Tensor<2,dim> >          elasticity_grad_phi (elasticity_dofs_per_cell);
    std::vector<double>                  elasticity_div_phi  (elasticity_dofs_per_cell);
    std::vector<Tensor<1,dim> >          elasticity_phi      (elasticity_dofs_per_cell);

    // Then comes the main loop over all cells and, as in step-27, the
    // initialization of the hp::FEValues object for the current cell and the
    // extraction of a FEValues object that is appropriate for the current
    // cell:
    typename hp::DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        hp_fe_values.reinit (cell);

        const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values();

        local_matrix.reinit (cell->get_fe().dofs_per_cell,
                             cell->get_fe().dofs_per_cell);
        local_rhs.reinit (cell->get_fe().dofs_per_cell);

        // With all of this done, we continue to assemble the cell terms for
        // cells that are part of the Stokes and elastic regions. While we
        // could in principle do this in one formula, in effect implementing
        // the one bilinear form stated in the introduction, we realize that
        // our finite element spaces are chosen in such a way that on each
        // cell, one set of variables (either velocities and pressure, or
        // displacements) are always zero, and consequently a more efficient
        // way of computing local integrals is to do only what's necessary
        // based on an <code>if</code> clause that tests which part of the
        // domain we are in.
        //
        // The actual computation of the local matrix is the same as in
        // step-22 as well as that given in the @ref vector_valued
        // documentation module for the elasticity equations:
        if (cell_is_in_fluid_domain (cell))
          {
            const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
            Assert (dofs_per_cell == stokes_dofs_per_cell,
                    ExcInternalError());

            for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q)
              {
                for (unsigned int k=0; k<dofs_per_cell; ++k)
                  {
                    stokes_symgrad_phi_u[k] = fe_values[velocities].symmetric_gradient (k, q);
                    stokes_div_phi_u[k]     = fe_values[velocities].divergence (k, q);
                    stokes_phi_p[k]         = fe_values[pressure].value (k, q);
                  }

                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  for (unsigned int j=0; j<dofs_per_cell; ++j)
                    local_matrix(i,j) += (2 * viscosity * stokes_symgrad_phi_u[i] * stokes_symgrad_phi_u[j]
                                          - stokes_div_phi_u[i] * stokes_phi_p[j]
                                          - stokes_phi_p[i] * stokes_div_phi_u[j])
                                         * fe_values.JxW(q);
              }
          }
        else
          {
            const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
            Assert (dofs_per_cell == elasticity_dofs_per_cell,
                    ExcInternalError());

            for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q)
              {
                for (unsigned int k=0; k<dofs_per_cell; ++k)
                  {
                    elasticity_grad_phi[k] = fe_values[displacements].gradient (k, q);
                    elasticity_div_phi[k]  = fe_values[displacements].divergence (k, q);
                  }

                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  for (unsigned int j=0; j<dofs_per_cell; ++j)
                    {
                      local_matrix(i,j)
                      +=  (lambda *
                           elasticity_div_phi[i] * elasticity_div_phi[j]
                           +
                           mu *
                           scalar_product(elasticity_grad_phi[i], elasticity_grad_phi[j])
                           +
                           mu *
                           scalar_product(elasticity_grad_phi[i], transpose(elasticity_grad_phi[j]))
                          )
                          *
                          fe_values.JxW(q);
                    }
              }
          }

        // Once we have the contributions from cell integrals, we copy them
        // into the global matrix (taking care of constraints right away,
        // through the ConstraintMatrix::distribute_local_to_global
        // function). Note that we have not written anything into the
        // <code>local_rhs</code> variable, though we still need to pass it
        // along since the elimination of nonzero boundary values requires the
        // modification of local and consequently also global right hand side
        // values:
        local_dof_indices.resize (cell->get_fe().dofs_per_cell);
        cell->get_dof_indices (local_dof_indices);
        constraints.distribute_local_to_global (local_matrix, local_rhs,
                                                local_dof_indices,
                                                system_matrix, system_rhs);

        // The more interesting part of this function is where we see about
        // face terms along the interface between the two subdomains. To this
        // end, we first have to make sure that we only assemble them once
        // even though a loop over all faces of all cells would encounter each
        // part of the interface twice. We arbitrarily make the decision that
        // we will only evaluate interface terms if the current cell is part
        // of the solid subdomain and if, consequently, a face is not at the
        // boundary and the potential neighbor behind it is part of the fluid
        // domain. Let's start with these conditions:
        if (cell_is_in_solid_domain (cell))
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
            if (cell->at_boundary(f) == false)
              {
                // At this point we know that the current cell is a candidate
                // for integration and that a neighbor behind face
                // <code>f</code> exists. There are now three possibilities:
                //
                // - The neighbor is at the same refinement level and has no
                //   children.
                // - The neighbor has children.
                // - The neighbor is coarser.
                //
                // In all three cases, we are only interested in it if it is
                // part of the fluid subdomain. So let us start with the first
                // and simplest case: if the neighbor is at the same level,
                // has no children, and is a fluid cell, then the two cells
                // share a boundary that is part of the interface along which
                // we want to integrate interface terms. All we have to do is
                // initialize two FEFaceValues object with the current face
                // and the face of the neighboring cell (note how we find out
                // which face of the neighboring cell borders on the current
                // cell) and pass things off to the function that evaluates
                // the interface terms (the third through fifth arguments to
                // this function provide it with scratch arrays). The result
                // is then again copied into the global matrix, using a
                // function that knows that the DoF indices of rows and
                // columns of the local matrix result from different cells:
                if ((cell->neighbor(f)->level() == cell->level())
                    &&
                    (cell->neighbor(f)->has_children() == false)
                    &&
                    cell_is_in_fluid_domain (cell->neighbor(f)))
                  {
                    elasticity_fe_face_values.reinit (cell, f);
                    stokes_fe_face_values.reinit (cell->neighbor(f),
                                                  cell->neighbor_of_neighbor(f));

                    assemble_interface_term (elasticity_fe_face_values, stokes_fe_face_values,
                                             elasticity_phi, stokes_symgrad_phi_u, stokes_phi_p,
                                             local_interface_matrix);

                    cell->neighbor(f)->get_dof_indices (neighbor_dof_indices);
                    constraints.distribute_local_to_global(local_interface_matrix,
                                                           local_dof_indices,
                                                           neighbor_dof_indices,
                                                           system_matrix);
                  }

                // The second case is if the neighbor has further children. In
                // that case, we have to loop over all the children of the
                // neighbor to see if they are part of the fluid subdomain. If
                // they are, then we integrate over the common interface,
                // which is a face for the neighbor and a subface of the
                // current cell, requiring us to use an FEFaceValues for the
                // neighbor and an FESubfaceValues for the current cell:
                else if ((cell->neighbor(f)->level() == cell->level())
                         &&
                         (cell->neighbor(f)->has_children() == true))
                  {
                    for (unsigned int subface=0;
                         subface<cell->face(f)->n_children();
                         ++subface)
                      if (cell_is_in_fluid_domain (cell->neighbor_child_on_subface
                                                   (f, subface)))
                        {
                          elasticity_fe_subface_values.reinit (cell,
                                                               f,
                                                               subface);
                          stokes_fe_face_values.reinit (cell->neighbor_child_on_subface (f, subface),
                                                        cell->neighbor_of_neighbor(f));

                          assemble_interface_term (elasticity_fe_subface_values,
                                                   stokes_fe_face_values,
                                                   elasticity_phi,
                                                   stokes_symgrad_phi_u, stokes_phi_p,
                                                   local_interface_matrix);

                          cell->neighbor_child_on_subface (f, subface)
                          ->get_dof_indices (neighbor_dof_indices);
                          constraints.distribute_local_to_global(local_interface_matrix,
                                                                 local_dof_indices,
                                                                 neighbor_dof_indices,
                                                                 system_matrix);
                        }
                  }

                // The last option is that the neighbor is coarser. In that
                // case we have to use an FESubfaceValues object for the
                // neighbor and a FEFaceValues for the current cell; the rest
                // is the same as before:
                else if (cell->neighbor_is_coarser(f)
                         &&
                         cell_is_in_fluid_domain(cell->neighbor(f)))
                  {
                    elasticity_fe_face_values.reinit (cell, f);
                    stokes_fe_subface_values.reinit (cell->neighbor(f),
                                                     cell->neighbor_of_coarser_neighbor(f).first,
                                                     cell->neighbor_of_coarser_neighbor(f).second);

                    assemble_interface_term (elasticity_fe_face_values,
                                             stokes_fe_subface_values,
                                             elasticity_phi,
                                             stokes_symgrad_phi_u, stokes_phi_p,
                                             local_interface_matrix);

                    cell->neighbor(f)->get_dof_indices (neighbor_dof_indices);
                    constraints.distribute_local_to_global(local_interface_matrix,
                                                           local_dof_indices,
                                                           neighbor_dof_indices,
                                                           system_matrix);

                  }
              }
      }
  }



  // In the function that assembles the global system, we passed computing
  // interface terms to a separate function we discuss here. The key is that
  // even though we can't predict the combination of FEFaceValues and
  // FESubfaceValues objects, they are both derived from the FEFaceValuesBase
  // class and consequently we don't have to care: the function is simply
  // called with two such objects denoting the values of the shape functions
  // on the quadrature points of the two sides of the face. We then do what we
  // always do: we fill the scratch arrays with the values of shape functions
  // and their derivatives, and then loop over all entries of the matrix to
  // compute the local integrals. The details of the bilinear form we evaluate
  // here are given in the introduction.
  template <int dim>
  void
  FluidStructureProblem<dim>::
  assemble_interface_term (const FEFaceValuesBase<dim>          &elasticity_fe_face_values,
                           const FEFaceValuesBase<dim>          &stokes_fe_face_values,
                           std::vector<Tensor<1,dim> >          &elasticity_phi,
                           std::vector<SymmetricTensor<2,dim> > &stokes_symgrad_phi_u,
                           std::vector<double>                  &stokes_phi_p,
                           FullMatrix<double>                   &local_interface_matrix) const
  {
    Assert (stokes_fe_face_values.n_quadrature_points ==
            elasticity_fe_face_values.n_quadrature_points,
            ExcInternalError());
    const unsigned int n_face_quadrature_points
      = elasticity_fe_face_values.n_quadrature_points;

    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);
    const FEValuesExtractors::Vector displacements (dim+1);

    local_interface_matrix = 0;
    for (unsigned int q=0; q<n_face_quadrature_points; ++q)
      {
        const Tensor<1,dim> normal_vector = stokes_fe_face_values.normal_vector(q);

        for (unsigned int k=0; k<stokes_fe_face_values.dofs_per_cell; ++k)
          stokes_symgrad_phi_u[k] = stokes_fe_face_values[velocities].symmetric_gradient (k, q);
        for (unsigned int k=0; k<elasticity_fe_face_values.dofs_per_cell; ++k)
          elasticity_phi[k] = elasticity_fe_face_values[displacements].value (k,q);

        for (unsigned int i=0; i<elasticity_fe_face_values.dofs_per_cell; ++i)
          for (unsigned int j=0; j<stokes_fe_face_values.dofs_per_cell; ++j)
            local_interface_matrix(i,j) += -((2 * viscosity *
                                              (stokes_symgrad_phi_u[j] *
                                               normal_vector)
                                              +
                                              stokes_phi_p[j] *
                                              normal_vector) *
                                             elasticity_phi[i] *
                                             stokes_fe_face_values.JxW(q));
      }
  }


  // @sect4{<code>FluidStructureProblem::solve</code>}

  // As discussed in the introduction, we use a rather trivial solver here: we
  // just pass the linear system off to the SparseDirectUMFPACK direct solver
  // (see, for example, step-29). The only thing we have to do after solving
  // is ensure that hanging node and boundary value constraints are correct.
  template <int dim>
  void
  FluidStructureProblem<dim>::solve ()
  {
    SparseDirectUMFPACK direct_solver;
    direct_solver.initialize (system_matrix);
    direct_solver.vmult (solution, system_rhs);

    constraints.distribute (solution);
  }



  // @sect4{<code>FluidStructureProblem::output_results</code>}

  // Generating graphical output is rather trivial here: all we have to do is
  // identify which components of the solution vector belong to scalars and/or
  // vectors (see, for example, step-22 for a previous example), and then pass
  // it all on to the DataOut class (with the second template argument equal
  // to hp::DoFHandler instead of the usual default DoFHandler):
  template <int dim>
  void
  FluidStructureProblem<dim>::
  output_results (const unsigned int refinement_cycle)  const
  {
    std::vector<std::string> solution_names (dim, "velocity");
    solution_names.push_back ("pressure");
    for (unsigned int d=0; d<dim; ++d)
      solution_names.push_back ("displacement");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation
    (dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation
    .push_back (DataComponentInterpretation::component_is_scalar);
    for (unsigned int d=0; d<dim; ++d)
      data_component_interpretation
      .push_back (DataComponentInterpretation::component_is_part_of_vector);

    DataOut<dim,hp::DoFHandler<dim> > data_out;
    data_out.attach_dof_handler (dof_handler);

    data_out.add_data_vector (solution, solution_names,
                              DataOut<dim,hp::DoFHandler<dim> >::type_dof_data,
                              data_component_interpretation);
    data_out.build_patches ();

    std::ostringstream filename;
    filename << "solution-"
             << Utilities::int_to_string (refinement_cycle, 2)
             << ".vtk";

    std::ofstream output (filename.str().c_str());
    data_out.write_vtk (output);
  }


  // @sect4{<code>FluidStructureProblem::refine_mesh</code>}

  // The next step is to refine the mesh. As was discussed in the
  // introduction, this is a bit tricky primarily because the fluid and the
  // solid subdomains use variables that have different physical dimensions
  // and for which the absolute magnitude of error estimates is consequently
  // not directly comparable. We will therefore have to scale them. At the top
  // of the function, we therefore first compute error estimates for the
  // different variables separately (using the velocities but not the pressure
  // for the fluid domain, and the displacements in the solid domain):
  template <int dim>
  void
  FluidStructureProblem<dim>::refine_mesh ()
  {
    Vector<float>
    stokes_estimated_error_per_cell (triangulation.n_active_cells());
    Vector<float>
    elasticity_estimated_error_per_cell (triangulation.n_active_cells());

    const QGauss<dim-1> stokes_face_quadrature(stokes_degree+2);
    const QGauss<dim-1> elasticity_face_quadrature(elasticity_degree+2);

    hp::QCollection<dim-1> face_q_collection;
    face_q_collection.push_back (stokes_face_quadrature);
    face_q_collection.push_back (elasticity_face_quadrature);

    const FEValuesExtractors::Vector velocities(0);
    KellyErrorEstimator<dim>::estimate (dof_handler,
                                        face_q_collection,
                                        typename FunctionMap<dim>::type(),
                                        solution,
                                        stokes_estimated_error_per_cell,
                                        fe_collection.component_mask(velocities));

    const FEValuesExtractors::Vector displacements(dim+1);
    KellyErrorEstimator<dim>::estimate (dof_handler,
                                        face_q_collection,
                                        typename FunctionMap<dim>::type(),
                                        solution,
                                        elasticity_estimated_error_per_cell,
                                        fe_collection.component_mask(displacements));

    // We then normalize error estimates by dividing by their norm and scale
    // the fluid error indicators by a factor of 4 as discussed in the
    // introduction. The results are then added together into a vector that
    // contains error indicators for all cells:
    stokes_estimated_error_per_cell
    *= 4. / stokes_estimated_error_per_cell.l2_norm();
    elasticity_estimated_error_per_cell
    *= 1. / elasticity_estimated_error_per_cell.l2_norm();

    Vector<float>
    estimated_error_per_cell (triangulation.n_active_cells());

    estimated_error_per_cell += stokes_estimated_error_per_cell;
    estimated_error_per_cell += elasticity_estimated_error_per_cell;

    // The second to last part of the function, before actually refining the
    // mesh, involves a heuristic that we have already mentioned in the
    // introduction: because the solution is discontinuous, the
    // KellyErrorEstimator class gets all confused about cells that sit at the
    // boundary between subdomains: it believes that the error is large there
    // because the jump in the gradient is large, even though this is entirely
    // expected and a feature that is in fact present in the exact solution as
    // well and therefore not indicative of any numerical error.
    //
    // Consequently, we set the error indicators to zero for all cells at the
    // interface; the conditions determining which cells this affects are
    // slightly awkward because we have to account for the possibility of
    // adaptively refined meshes, meaning that the neighboring cell can be
    // coarser than the current one, or could in fact be refined some
    // more. The structure of these nested conditions is much the same as we
    // encountered when assembling interface terms in
    // <code>assemble_system</code>.
    {
      unsigned int cell_index = 0;
      for (typename hp::DoFHandler<dim>::active_cell_iterator
           cell = dof_handler.begin_active();
           cell != dof_handler.end(); ++cell, ++cell_index)
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell_is_in_solid_domain (cell))
            {
              if ((cell->at_boundary(f) == false)
                  &&
                  (((cell->neighbor(f)->level() == cell->level())
                    &&
                    (cell->neighbor(f)->has_children() == false)
                    &&
                    cell_is_in_fluid_domain (cell->neighbor(f)))
                   ||
                   ((cell->neighbor(f)->level() == cell->level())
                    &&
                    (cell->neighbor(f)->has_children() == true)
                    &&
                    (cell_is_in_fluid_domain (cell->neighbor_child_on_subface
                                              (f, 0))))
                   ||
                   (cell->neighbor_is_coarser(f)
                    &&
                    cell_is_in_fluid_domain(cell->neighbor(f)))
                  ))
                estimated_error_per_cell(cell_index) = 0;
            }
          else
            {
              if ((cell->at_boundary(f) == false)
                  &&
                  (((cell->neighbor(f)->level() == cell->level())
                    &&
                    (cell->neighbor(f)->has_children() == false)
                    &&
                    cell_is_in_solid_domain (cell->neighbor(f)))
                   ||
                   ((cell->neighbor(f)->level() == cell->level())
                    &&
                    (cell->neighbor(f)->has_children() == true)
                    &&
                    (cell_is_in_solid_domain (cell->neighbor_child_on_subface
                                              (f, 0))))
                   ||
                   (cell->neighbor_is_coarser(f)
                    &&
                    cell_is_in_solid_domain(cell->neighbor(f)))
                  ))
                estimated_error_per_cell(cell_index) = 0;
            }
    }

    GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                     estimated_error_per_cell,
                                                     0.3, 0.0);
    triangulation.execute_coarsening_and_refinement ();
  }



  // @sect4{<code>FluidStructureProblem::run</code>}

  // This is, as usual, the function that controls the overall flow of
  // operation. If you've read through tutorial programs step-1 through
  // step-6, for example, then you are already quite familiar with the
  // following structure:
  template <int dim>
  void FluidStructureProblem<dim>::run ()
  {
    make_grid ();

    for (unsigned int refinement_cycle = 0; refinement_cycle<10-2*dim;
         ++refinement_cycle)
      {
        std::cout << "Refinement cycle " << refinement_cycle << std::endl;

        if (refinement_cycle > 0)
          refine_mesh ();

        setup_dofs ();

        std::cout << "   Assembling..." << std::endl;
        assemble_system ();

        std::cout << "   Solving..." << std::endl;
        solve ();

        std::cout << "   Writing output..." << std::endl;
        output_results (refinement_cycle);

        std::cout << std::endl;
      }
  }
}



// @sect4{The <code>main()</code> function}

// This, final, function contains pretty much exactly what most of the other
// tutorial programs have:
int main ()
{
  try
    {
      using namespace dealii;
      using namespace Step46;

      deallog.depth_console (0);

      FluidStructureProblem<2> flow_problem(1, 1);
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
