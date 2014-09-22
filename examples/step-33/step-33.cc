/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2007 - 2013 by the deal.II authors
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
 * Author: David Neckels, Boulder, Colorado, 2007, 2008
 */


// @sect3{Include files}

// First a standard set of deal.II includes. Nothing special to comment on
// here:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/solution_transfer.h>

// Then, as mentioned in the introduction, we use various Trilinos packages as
// linear solvers as well as for automatic differentiation. These are in the
// following include files.
//
// Since deal.II provides interfaces to the basic Trilinos matrices, vectors,
// preconditioners and solvers, we include them similarly as deal.II linear
// algebra structures.
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>


// Sacado is the automatic differentiation package within Trilinos, which is
// used to find the Jacobian for a fully implicit Newton iteration:
#include <Sacado.hpp>


// And this again is C++:
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>

// To end this section, introduce everything in the dealii library into the
// namespace into which the contents of this program will go:
namespace Step33
{
  using namespace dealii;


  // @sect3{Euler equation specifics}

  // Here we define the flux function for this particular system of
  // conservation laws, as well as pretty much everything else that's specific
  // to the Euler equations for gas dynamics, for reasons discussed in the
  // introduction. We group all this into a structure that defines everything
  // that has to do with the flux. All members of this structure are static,
  // i.e. the structure has no actual state specified by instance member
  // variables. The better way to do this, rather than a structure with all
  // static members would be to use a namespace -- but namespaces can't be
  // templatized and we want some of the member variables of the structure to
  // depend on the space dimension, which we in our usual way introduce using
  // a template parameter.
  template <int dim>
  struct EulerEquations
  {
    // @sect4{Component description}

    // First a few variables that describe the various components of our
    // solution vector in a generic way. This includes the number of
    // components in the system (Euler's equations have one entry for momenta
    // in each spatial direction, plus the energy and density components, for
    // a total of <code>dim+2</code> components), as well as functions that
    // describe the index within the solution vector of the first momentum
    // component, the density component, and the energy density
    // component. Note that all these %numbers depend on the space dimension;
    // defining them in a generic way (rather than by implicit convention)
    // makes our code more flexible and makes it easier to later extend it,
    // for example by adding more components to the equations.
    static const unsigned int n_components             = dim + 2;
    static const unsigned int first_momentum_component = 0;
    static const unsigned int density_component        = dim;
    static const unsigned int energy_component         = dim+1;

    // When generating graphical output way down in this program, we need to
    // specify the names of the solution variables as well as how the various
    // components group into vector and scalar fields. We could describe this
    // there, but in order to keep things that have to do with the Euler
    // equation localized here and the rest of the program as generic as
    // possible, we provide this sort of information in the following two
    // functions:
    static
    std::vector<std::string>
    component_names ()
    {
      std::vector<std::string> names (dim, "momentum");
      names.push_back ("density");
      names.push_back ("energy_density");

      return names;
    }


    static
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    component_interpretation ()
    {
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation
      (dim, DataComponentInterpretation::component_is_part_of_vector);
      data_component_interpretation
      .push_back (DataComponentInterpretation::component_is_scalar);
      data_component_interpretation
      .push_back (DataComponentInterpretation::component_is_scalar);

      return data_component_interpretation;
    }


    // @sect4{Transformations between variables}

    // Next, we define the gas constant. We will set it to 1.4 in its
    // definition immediately following the declaration of this class (unlike
    // integer variables, like the ones above, static const floating point
    // member variables cannot be initialized within the class declaration in
    // C++). This value of 1.4 is representative of a gas that consists of
    // molecules composed of two atoms, such as air which consists up to small
    // traces almost entirely of $N_2$ and $O_2$.
    static const double gas_gamma;


    // In the following, we will need to compute the kinetic energy and the
    // pressure from a vector of conserved variables. This we can do based on
    // the energy density and the kinetic energy $\frac 12 \rho |\mathbf v|^2
    // = \frac{|\rho \mathbf v|^2}{2\rho}$ (note that the independent
    // variables contain the momentum components $\rho v_i$, not the
    // velocities $v_i$).
    //
    // There is one slight problem: We will need to call the following
    // functions with input arguments of type
    // <code>std::vector@<number@></code> and
    // <code>Vector@<number@></code>. The problem is that the former has an
    // access operator <code>operator[]</code> whereas the latter, for
    // historical reasons, has <code>operator()</code>. We wouldn't be able to
    // write the function in a generic way if we were to use one or the other
    // of these. Fortunately, we can use the following trick: instead of
    // writing <code>v[i]</code> or <code>v(i)</code>, we can use
    // <code>*(v.begin() + i)</code>, i.e. we generate an iterator that points
    // to the <code>i</code>th element, and then dereference it. This works
    // for both kinds of vectors -- not the prettiest solution, but one that
    // works.
    template <typename number, typename InputVector>
    static
    number
    compute_kinetic_energy (const InputVector &W)
    {
      number kinetic_energy = 0;
      for (unsigned int d=0; d<dim; ++d)
        kinetic_energy += *(W.begin()+first_momentum_component+d) *
                          *(W.begin()+first_momentum_component+d);
      kinetic_energy *= 1./(2 * *(W.begin() + density_component));

      return kinetic_energy;
    }


    template <typename number, typename InputVector>
    static
    number
    compute_pressure (const InputVector &W)
    {
      return ((gas_gamma-1.0) *
              (*(W.begin() + energy_component) -
               compute_kinetic_energy<number>(W)));
    }


    // @sect4{EulerEquations::compute_flux_matrix}

    // We define the flux function $F(W)$ as one large matrix.  Each row of
    // this matrix represents a scalar conservation law for the component in
    // that row.  The exact form of this matrix is given in the
    // introduction. Note that we know the size of the matrix: it has as many
    // rows as the system has components, and <code>dim</code> columns; rather
    // than using a FullMatrix object for such a matrix (which has a variable
    // number of rows and columns and must therefore allocate memory on the
    // heap each time such a matrix is created), we use a rectangular array of
    // numbers right away.
    //
    // We templatize the numerical type of the flux function so that we may
    // use the automatic differentiation type here.  Similarly, we will call
    // the function with different input vector data types, so we templatize
    // on it as well:
    template <typename InputVector, typename number>
    static
    void compute_flux_matrix (const InputVector &W,
                              number (&flux)[n_components][dim])
    {
      // First compute the pressure that appears in the flux matrix, and then
      // compute the first <code>dim</code> columns of the matrix that
      // correspond to the momentum terms:
      const number pressure = compute_pressure<number> (W);

      for (unsigned int d=0; d<dim; ++d)
        {
          for (unsigned int e=0; e<dim; ++e)
            flux[first_momentum_component+d][e]
              = W[first_momentum_component+d] *
                W[first_momentum_component+e] /
                W[density_component];

          flux[first_momentum_component+d][d] += pressure;
        }

      // Then the terms for the density (i.e. mass conservation), and, lastly,
      // conservation of energy:
      for (unsigned int d=0; d<dim; ++d)
        flux[density_component][d] = W[first_momentum_component+d];

      for (unsigned int d=0; d<dim; ++d)
        flux[energy_component][d] = W[first_momentum_component+d] /
                                    W[density_component] *
                                    (W[energy_component] + pressure);
    }


    // @sect4{EulerEquations::compute_normal_flux}

    // On the boundaries of the domain and across hanging nodes we use a
    // numerical flux function to enforce boundary conditions.  This routine
    // is the basic Lax-Friedrich's flux with a stabilization parameter
    // $\alpha$. It's form has also been given already in the introduction:
    template <typename InputVector>
    static
    void numerical_normal_flux (const Point<dim>          &normal,
                                const InputVector         &Wplus,
                                const InputVector         &Wminus,
                                const double               alpha,
                                Sacado::Fad::DFad<double> (&normal_flux)[n_components])
    {
      Sacado::Fad::DFad<double> iflux[n_components][dim];
      Sacado::Fad::DFad<double> oflux[n_components][dim];

      compute_flux_matrix (Wplus, iflux);
      compute_flux_matrix (Wminus, oflux);

      for (unsigned int di=0; di<n_components; ++di)
        {
          normal_flux[di] = 0;
          for (unsigned int d=0; d<dim; ++d)
            normal_flux[di] += 0.5*(iflux[di][d] + oflux[di][d]) * normal[d];

          normal_flux[di] += 0.5*alpha*(Wplus[di] - Wminus[di]);
        }
    }

    // @sect4{EulerEquations::compute_forcing_vector}

    // In the same way as describing the flux function $\mathbf F(\mathbf w)$,
    // we also need to have a way to describe the right hand side forcing
    // term. As mentioned in the introduction, we consider only gravity here,
    // which leads to the specific form $\mathbf G(\mathbf w) = \left(
    // g_1\rho, g_2\rho, g_3\rho, 0, \rho \mathbf g \cdot \mathbf v
    // \right)^T$, shown here for the 3d case. More specifically, we will
    // consider only $\mathbf g=(0,0,-1)^T$ in 3d, or $\mathbf g=(0,-1)^T$ in
    // 2d. This naturally leads to the following function:
    template <typename InputVector, typename number>
    static
    void compute_forcing_vector (const InputVector &W,
                                 number (&forcing)[n_components])
    {
      const double gravity = -1.0;

      for (unsigned int c=0; c<n_components; ++c)
        switch (c)
          {
          case first_momentum_component+dim-1:
            forcing[c] = gravity * W[density_component];
            break;
          case energy_component:
            forcing[c] = gravity *
                         W[density_component] *
                         W[first_momentum_component+dim-1];
            break;
          default:
            forcing[c] = 0;
          }
    }


    // @sect4{Dealing with boundary conditions}

    // Another thing we have to deal with is boundary conditions. To this end,
    // let us first define the kinds of boundary conditions we currently know
    // how to deal with:
    enum BoundaryKind
    {
      inflow_boundary,
      outflow_boundary,
      no_penetration_boundary,
      pressure_boundary
    };


    // The next part is to actually decide what to do at each kind of
    // boundary. To this end, remember from the introduction that boundary
    // conditions are specified by choosing a value $\mathbf w^-$ on the
    // outside of a boundary given an inhomogeneity $\mathbf j$ and possibly
    // the solution's value $\mathbf w^+$ on the inside. Both are then passed
    // to the numerical flux $\mathbf H(\mathbf{w}^+, \mathbf{w}^-,
    // \mathbf{n})$ to define boundary contributions to the bilinear form.
    //
    // Boundary conditions can in some cases be specified for each component
    // of the solution vector independently. For example, if component $c$ is
    // marked for inflow, then $w^-_c = j_c$. If it is an outflow, then $w^-_c
    // = w^+_c$. These two simple cases are handled first in the function
    // below.
    //
    // There is a little snag that makes this function unpleasant from a C++
    // language viewpoint: The output vector <code>Wminus</code> will of
    // course be modified, so it shouldn't be a <code>const</code>
    // argument. Yet it is in the implementation below, and needs to be in
    // order to allow the code to compile. The reason is that we call this
    // function at a place where <code>Wminus</code> is of type
    // <code>Table@<2,Sacado::Fad::DFad@<double@> @></code>, this being 2d
    // table with indices representing the quadrature point and the vector
    // component, respectively. We call this function with
    // <code>Wminus[q]</code> as last argument; subscripting a 2d table yields
    // a temporary accessor object representing a 1d vector, just what we want
    // here. The problem is that a temporary accessor object can't be bound to
    // a non-const reference argument of a function, as we would like here,
    // according to the C++ 1998 and 2003 standards (something that will be
    // fixed with the next standard in the form of rvalue references).  We get
    // away with making the output argument here a constant because it is the
    // <i>accessor</i> object that's constant, not the table it points to:
    // that one can still be written to. The hack is unpleasant nevertheless
    // because it restricts the kind of data types that may be used as
    // template argument to this function: a regular vector isn't going to do
    // because that one can not be written to when marked
    // <code>const</code>. With no good solution around at the moment, we'll
    // go with the pragmatic, even if not pretty, solution shown here:
    template <typename DataVector>
    static
    void
    compute_Wminus (const BoundaryKind  (&boundary_kind)[n_components],
                    const Point<dim>     &normal_vector,
                    const DataVector     &Wplus,
                    const Vector<double> &boundary_values,
                    const DataVector     &Wminus)
    {
      for (unsigned int c = 0; c < n_components; c++)
        switch (boundary_kind[c])
          {
          case inflow_boundary:
          {
            Wminus[c] = boundary_values(c);
            break;
          }

          case outflow_boundary:
          {
            Wminus[c] = Wplus[c];
            break;
          }

          // Prescribed pressure boundary conditions are a bit more
          // complicated by the fact that even though the pressure is
          // prescribed, we really are setting the energy component here,
          // which will depend on velocity and pressure. So even though this
          // seems like a Dirichlet type boundary condition, we get
          // sensitivities of energy to velocity and density (unless these are
          // also prescribed):
          case pressure_boundary:
          {
            const typename DataVector::value_type
            density = (boundary_kind[density_component] ==
                       inflow_boundary
                       ?
                       boundary_values(density_component)
                       :
                       Wplus[density_component]);

            typename DataVector::value_type kinetic_energy = 0;
            for (unsigned int d=0; d<dim; ++d)
              if (boundary_kind[d] == inflow_boundary)
                kinetic_energy += boundary_values(d)*boundary_values(d);
              else
                kinetic_energy += Wplus[d]*Wplus[d];
            kinetic_energy *= 1./2./density;

            Wminus[c] = boundary_values(c) / (gas_gamma-1.0) +
                        kinetic_energy;

            break;
          }

          case no_penetration_boundary:
          {
            // We prescribe the velocity (we are dealing with a particular
            // component here so that the average of the velocities is
            // orthogonal to the surface normal.  This creates sensitivities of
            // across the velocity components.
            Sacado::Fad::DFad<double> vdotn = 0;
            for (unsigned int d = 0; d < dim; d++)
              {
                vdotn += Wplus[d]*normal_vector[d];
              }

            Wminus[c] = Wplus[c] - 2.0*vdotn*normal_vector[c];
            break;
          }

          default:
            Assert (false, ExcNotImplemented());
          }
    }


    // @sect4{EulerEquations::compute_refinement_indicators}

    // In this class, we also want to specify how to refine the mesh. The
    // class <code>ConservationLaw</code> that will use all the information we
    // provide here in the <code>EulerEquation</code> class is pretty agnostic
    // about the particular conservation law it solves: as doesn't even really
    // care how many components a solution vector has. Consequently, it can't
    // know what a reasonable refinement indicator would be. On the other
    // hand, here we do, or at least we can come up with a reasonable choice:
    // we simply look at the gradient of the density, and compute
    // $\eta_K=\log\left(1+|\nabla\rho(x_K)|\right)$, where $x_K$ is the
    // center of cell $K$.
    //
    // There are certainly a number of equally reasonable refinement
    // indicators, but this one does, and it is easy to compute:
    static
    void
    compute_refinement_indicators (const DoFHandler<dim> &dof_handler,
                                   const Mapping<dim>    &mapping,
                                   const Vector<double>  &solution,
                                   Vector<double>        &refinement_indicators)
    {
      const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;
      std::vector<unsigned int> dofs (dofs_per_cell);

      const QMidpoint<dim>  quadrature_formula;
      const UpdateFlags update_flags = update_gradients;
      FEValues<dim> fe_v (mapping, dof_handler.get_fe(),
                          quadrature_formula, update_flags);

      std::vector<std::vector<Tensor<1,dim> > >
      dU (1, std::vector<Tensor<1,dim> >(n_components));

      typename DoFHandler<dim>::active_cell_iterator
      cell = dof_handler.begin_active(),
      endc = dof_handler.end();
      for (unsigned int cell_no=0; cell!=endc; ++cell, ++cell_no)
        {
          fe_v.reinit(cell);
          fe_v.get_function_grads (solution, dU);

          refinement_indicators(cell_no)
            = std::log(1+
                       std::sqrt(dU[0][density_component] *
                                 dU[0][density_component]));
        }
    }



    // @sect4{EulerEquations::Postprocessor}

    // Finally, we declare a class that implements a postprocessing of data
    // components. The problem this class solves is that the variables in the
    // formulation of the Euler equations we use are in conservative rather
    // than physical form: they are momentum densities $\mathbf m=\rho\mathbf
    // v$, density $\rho$, and energy density $E$. What we would like to also
    // put into our output file are velocities $\mathbf v=\frac{\mathbf
    // m}{\rho}$ and pressure $p=(\gamma-1)(E-\frac{1}{2} \rho |\mathbf
    // v|^2)$.
    //
    // In addition, we would like to add the possibility to generate schlieren
    // plots. Schlieren plots are a way to visualize shocks and other sharp
    // interfaces. The word "schlieren" is a German word that may be
    // translated as "striae" -- it may be simpler to explain it by an
    // example, however: schlieren is what you see when you, for example, pour
    // highly concentrated alcohol, or a transparent saline solution, into
    // water; the two have the same color, but they have different refractive
    // indices and so before they are fully mixed light goes through the
    // mixture along bent rays that lead to brightness variations if you look
    // at it. That's "schlieren". A similar effect happens in compressible
    // flow because the refractive index depends on the pressure (and
    // therefore the density) of the gas.
    //
    // The origin of the word refers to two-dimensional projections of a
    // three-dimensional volume (we see a 2d picture of the 3d fluid). In
    // computational fluid dynamics, we can get an idea of this effect by
    // considering what causes it: density variations. Schlieren plots are
    // therefore produced by plotting $s=|\nabla \rho|^2$; obviously, $s$ is
    // large in shocks and at other highly dynamic places. If so desired by
    // the user (by specifying this in the input file), we would like to
    // generate these schlieren plots in addition to the other derived
    // quantities listed above.
    //
    // The implementation of the algorithms to compute derived quantities from
    // the ones that solve our problem, and to output them into data file,
    // rests on the DataPostprocessor class. It has extensive documentation,
    // and other uses of the class can also be found in step-29. We therefore
    // refrain from extensive comments.
    class Postprocessor : public DataPostprocessor<dim>
    {
    public:
      Postprocessor (const bool do_schlieren_plot);

      virtual
      void
      compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                         const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                         const std::vector<std::vector<Tensor<2,dim> > > &dduh,
                                         const std::vector<Point<dim> >                  &normals,
                                         const std::vector<Point<dim> >                  &evaluation_points,
                                         std::vector<Vector<double> >                    &computed_quantities) const;

      virtual std::vector<std::string> get_names () const;

      virtual
      std::vector<DataComponentInterpretation::DataComponentInterpretation>
      get_data_component_interpretation () const;

      virtual UpdateFlags get_needed_update_flags () const;

    private:
      const bool do_schlieren_plot;
    };
  };


  template <int dim>
  const double EulerEquations<dim>::gas_gamma = 1.4;



  template <int dim>
  EulerEquations<dim>::Postprocessor::
  Postprocessor (const bool do_schlieren_plot)
    :
    do_schlieren_plot (do_schlieren_plot)
  {}


  // This is the only function worth commenting on. When generating graphical
  // output, the DataOut and related classes will call this function on each
  // cell, with values, gradients, Hessians, and normal vectors (in case we're
  // working on faces) at each quadrature point. Note that the data at each
  // quadrature point is itself vector-valued, namely the conserved
  // variables. What we're going to do here is to compute the quantities we're
  // interested in at each quadrature point. Note that for this we can ignore
  // the Hessians ("dduh") and normal vectors; to avoid compiler warnings
  // about unused variables, we comment out their names.
  template <int dim>
  void
  EulerEquations<dim>::Postprocessor::
  compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                     const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                     const std::vector<std::vector<Tensor<2,dim> > > &/*dduh*/,
                                     const std::vector<Point<dim> >                  &/*normals*/,
                                     const std::vector<Point<dim> >                  &/*evaluation_points*/,
                                     std::vector<Vector<double> >                    &computed_quantities) const
  {
    // At the beginning of the function, let us make sure that all variables
    // have the correct sizes, so that we can access individual vector
    // elements without having to wonder whether we might read or write
    // invalid elements; we also check that the <code>duh</code> vector only
    // contains data if we really need it (the system knows about this because
    // we say so in the <code>get_needed_update_flags()</code> function
    // below). For the inner vectors, we check that at least the first element
    // of the outer vector has the correct inner size:
    const unsigned int n_quadrature_points = uh.size();

    if (do_schlieren_plot == true)
      Assert (duh.size() == n_quadrature_points,
              ExcInternalError())
      else
        Assert (duh.size() == 0,
                ExcInternalError());

    Assert (computed_quantities.size() == n_quadrature_points,
            ExcInternalError());

    Assert (uh[0].size() == n_components,
            ExcInternalError());

    if (do_schlieren_plot == true)
      Assert (computed_quantities[0].size() == dim+2, ExcInternalError())
      else
        Assert (computed_quantities[0].size() == dim+1, ExcInternalError());

    // Then loop over all quadrature points and do our work there. The code
    // should be pretty self-explanatory. The order of output variables is
    // first <code>dim</code> velocities, then the pressure, and if so desired
    // the schlieren plot. Note that we try to be generic about the order of
    // variables in the input vector, using the
    // <code>first_momentum_component</code> and
    // <code>density_component</code> information:
    for (unsigned int q=0; q<n_quadrature_points; ++q)
      {
        const double density = uh[q](density_component);

        for (unsigned int d=0; d<dim; ++d)
          computed_quantities[q](d)
            = uh[q](first_momentum_component+d) / density;

        computed_quantities[q](dim) = compute_pressure<double> (uh[q]);

        if (do_schlieren_plot == true)
          computed_quantities[q](dim+1) = duh[q][density_component] *
                                          duh[q][density_component];
      }
  }


  template <int dim>
  std::vector<std::string>
  EulerEquations<dim>::Postprocessor::
  get_names () const
  {
    std::vector<std::string> names;
    for (unsigned int d=0; d<dim; ++d)
      names.push_back ("velocity");
    names.push_back ("pressure");

    if (do_schlieren_plot == true)
      names.push_back ("schlieren_plot");

    return names;
  }


  template <int dim>
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
  EulerEquations<dim>::Postprocessor::
  get_data_component_interpretation () const
  {
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation (dim,
                    DataComponentInterpretation::component_is_part_of_vector);

    interpretation.push_back (DataComponentInterpretation::
                              component_is_scalar);

    if (do_schlieren_plot == true)
      interpretation.push_back (DataComponentInterpretation::
                                component_is_scalar);

    return interpretation;
  }



  template <int dim>
  UpdateFlags
  EulerEquations<dim>::Postprocessor::
  get_needed_update_flags () const
  {
    if (do_schlieren_plot == true)
      return update_values | update_gradients;
    else
      return update_values;
  }


  // @sect3{Run time parameter handling}

  // Our next job is to define a few classes that will contain run-time
  // parameters (for example solver tolerances, number of iterations,
  // stabilization parameter, and the like). One could do this in the main
  // class, but we separate it from that one to make the program more modular
  // and easier to read: Everything that has to do with run-time parameters
  // will be in the following namespace, whereas the program logic is in the
  // main class.
  //
  // We will split the run-time parameters into a few separate structures,
  // which we will all put into a namespace <code>Parameters</code>. Of these
  // classes, there are a few that group the parameters for individual groups,
  // such as for solvers, mesh refinement, or output. Each of these classes
  // have functions <code>declare_parameters()</code> and
  // <code>parse_parameters()</code> that declare parameter subsections and
  // entries in a ParameterHandler object, and retrieve actual parameter
  // values from such an object, respectively. These classes declare all their
  // parameters in subsections of the ParameterHandler.
  //
  // The final class of the following namespace combines all the previous
  // classes by deriving from them and taking care of a few more entries at
  // the top level of the input file, as well as a few odd other entries in
  // subsections that are too short to warrant a structure by themselves.
  //
  // It is worth pointing out one thing here: None of the classes below have a
  // constructor that would initialize the various member variables. This
  // isn't a problem, however, since we will read all variables declared in
  // these classes from the input file (or indirectly: a ParameterHandler
  // object will read it from there, and we will get the values from this
  // object), and they will be initialized this way. In case a certain
  // variable is not specified at all in the input file, this isn't a problem
  // either: The ParameterHandler class will in this case simply take the
  // default value that was specified when declaring an entry in the
  // <code>declare_parameters()</code> functions of the classes below.
  namespace Parameters
  {

    // @sect4{Parameters::Solver}
    //
    // The first of these classes deals with parameters for the linear inner
    // solver. It offers parameters that indicate which solver to use (GMRES
    // as a solver for general non-symmetric indefinite systems, or a sparse
    // direct solver), the amount of output to be produced, as well as various
    // parameters that tweak the thresholded incomplete LU decomposition
    // (ILUT) that we use as a preconditioner for GMRES.
    //
    // In particular, the ILUT takes the following parameters:
    // - ilut_fill: the number of extra entries to add when forming the ILU
    //   decomposition
    // - ilut_atol, ilut_rtol: When forming the preconditioner, for certain
    //   problems bad conditioning (or just bad luck) can cause the
    //   preconditioner to be very poorly conditioned.  Hence it can help to
    //   add diagonal perturbations to the original matrix and form the
    //   preconditioner for this slightly better matrix.  ATOL is an absolute
    //   perturbation that is added to the diagonal before forming the prec,
    //   and RTOL is a scaling factor $rtol \geq 1$.
    // - ilut_drop: The ILUT will drop any values that have magnitude less
    //   than this value.  This is a way to manage the amount of memory used
    //   by this preconditioner.
    //
    // The meaning of each parameter is also briefly described in the third
    // argument of the ParameterHandler::declare_entry call in
    // <code>declare_parameters()</code>.
    struct Solver
    {
      enum SolverType { gmres, direct };
      SolverType solver;

      enum  OutputType { quiet, verbose };
      OutputType output;

      double linear_residual;
      int max_iterations;

      double ilut_fill;
      double ilut_atol;
      double ilut_rtol;
      double ilut_drop;

      static void declare_parameters (ParameterHandler &prm);
      void parse_parameters (ParameterHandler &prm);
    };



    void Solver::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("linear solver");
      {
        prm.declare_entry("output", "quiet",
                          Patterns::Selection("quiet|verbose"),
                          "State whether output from solver runs should be printed. "
                          "Choices are <quiet|verbose>.");
        prm.declare_entry("method", "gmres",
                          Patterns::Selection("gmres|direct"),
                          "The kind of solver for the linear system. "
                          "Choices are <gmres|direct>.");
        prm.declare_entry("residual", "1e-10",
                          Patterns::Double(),
                          "Linear solver residual");
        prm.declare_entry("max iters", "300",
                          Patterns::Integer(),
                          "Maximum solver iterations");
        prm.declare_entry("ilut fill", "2",
                          Patterns::Double(),
                          "Ilut preconditioner fill");
        prm.declare_entry("ilut absolute tolerance", "1e-9",
                          Patterns::Double(),
                          "Ilut preconditioner tolerance");
        prm.declare_entry("ilut relative tolerance", "1.1",
                          Patterns::Double(),
                          "Ilut relative tolerance");
        prm.declare_entry("ilut drop tolerance", "1e-10",
                          Patterns::Double(),
                          "Ilut drop tolerance");
      }
      prm.leave_subsection();
    }




    void Solver::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("linear solver");
      {
        const std::string op = prm.get("output");
        if (op == "verbose")
          output = verbose;
        if (op == "quiet")
          output = quiet;

        const std::string sv = prm.get("method");
        if (sv == "direct")
          solver = direct;
        else if (sv == "gmres")
          solver = gmres;

        linear_residual = prm.get_double("residual");
        max_iterations  = prm.get_integer("max iters");
        ilut_fill       = prm.get_double("ilut fill");
        ilut_atol       = prm.get_double("ilut absolute tolerance");
        ilut_rtol       = prm.get_double("ilut relative tolerance");
        ilut_drop       = prm.get_double("ilut drop tolerance");
      }
      prm.leave_subsection();
    }



    // @sect4{Parameters::Refinement}
    //
    // Similarly, here are a few parameters that determine how the mesh is to
    // be refined (and if it is to be refined at all). For what exactly the
    // shock parameters do, see the mesh refinement functions further down.
    struct Refinement
    {
      bool do_refine;
      double shock_val;
      double shock_levels;

      static void declare_parameters (ParameterHandler &prm);
      void parse_parameters (ParameterHandler &prm);
    };



    void Refinement::declare_parameters (ParameterHandler &prm)
    {

      prm.enter_subsection("refinement");
      {
        prm.declare_entry("refinement", "true",
                          Patterns::Bool(),
                          "Whether to perform mesh refinement or not");
        prm.declare_entry("refinement fraction", "0.1",
                          Patterns::Double(),
                          "Fraction of high refinement");
        prm.declare_entry("unrefinement fraction", "0.1",
                          Patterns::Double(),
                          "Fraction of low unrefinement");
        prm.declare_entry("max elements", "1000000",
                          Patterns::Double(),
                          "maximum number of elements");
        prm.declare_entry("shock value", "4.0",
                          Patterns::Double(),
                          "value for shock indicator");
        prm.declare_entry("shock levels", "3.0",
                          Patterns::Double(),
                          "number of shock refinement levels");
      }
      prm.leave_subsection();
    }


    void Refinement::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("refinement");
      {
        do_refine     = prm.get_bool ("refinement");
        shock_val     = prm.get_double("shock value");
        shock_levels  = prm.get_double("shock levels");
      }
      prm.leave_subsection();
    }



    // @sect4{Parameters::Flux}
    //
    // Next a section on flux modifications to make it more stable. In
    // particular, two options are offered to stabilize the Lax-Friedrichs
    // flux: either choose $\mathbf{H}(\mathbf{a},\mathbf{b},\mathbf{n}) =
    // \frac{1}{2}(\mathbf{F}(\mathbf{a})\cdot \mathbf{n} +
    // \mathbf{F}(\mathbf{b})\cdot \mathbf{n} + \alpha (\mathbf{a} -
    // \mathbf{b}))$ where $\alpha$ is either a fixed number specified in the
    // input file, or where $\alpha$ is a mesh dependent value. In the latter
    // case, it is chosen as $\frac{h}{2\delta T}$ with $h$ the diameter of
    // the face to which the flux is applied, and $\delta T$ the current time
    // step.
    struct Flux
    {
      enum StabilizationKind { constant, mesh_dependent };
      StabilizationKind stabilization_kind;

      double stabilization_value;

      static void declare_parameters (ParameterHandler &prm);
      void parse_parameters (ParameterHandler &prm);
    };


    void Flux::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("flux");
      {
        prm.declare_entry("stab", "mesh",
                          Patterns::Selection("constant|mesh"),
                          "Whether to use a constant stabilization parameter or "
                          "a mesh-dependent one");
        prm.declare_entry("stab value", "1",
                          Patterns::Double(),
                          "alpha stabilization");
      }
      prm.leave_subsection();
    }


    void Flux::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("flux");
      {
        const std::string stab = prm.get("stab");
        if (stab == "constant")
          stabilization_kind = constant;
        else if (stab == "mesh")
          stabilization_kind = mesh_dependent;
        else
          AssertThrow (false, ExcNotImplemented());

        stabilization_value = prm.get_double("stab value");
      }
      prm.leave_subsection();
    }



    // @sect4{Parameters::Output}
    //
    // Then a section on output parameters. We offer to produce Schlieren
    // plots (the squared gradient of the density, a tool to visualize shock
    // fronts), and a time interval between graphical output in case we don't
    // want an output file every time step.
    struct Output
    {
      bool schlieren_plot;
      double output_step;

      static void declare_parameters (ParameterHandler &prm);
      void parse_parameters (ParameterHandler &prm);
    };



    void Output::declare_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("output");
      {
        prm.declare_entry("schlieren plot", "true",
                          Patterns::Bool (),
                          "Whether or not to produce schlieren plots");
        prm.declare_entry("step", "-1",
                          Patterns::Double(),
                          "Output once per this period");
      }
      prm.leave_subsection();
    }



    void Output::parse_parameters (ParameterHandler &prm)
    {
      prm.enter_subsection("output");
      {
        schlieren_plot = prm.get_bool("schlieren plot");
        output_step = prm.get_double("step");
      }
      prm.leave_subsection();
    }



    // @sect4{Parameters::AllParameters}
    //
    // Finally the class that brings it all together. It declares a number of
    // parameters itself, mostly ones at the top level of the parameter file
    // as well as several in section too small to warrant their own
    // classes. It also contains everything that is actually space dimension
    // dependent, like initial or boundary conditions.
    //
    // Since this class is derived from all the ones above, the
    // <code>declare_parameters()</code> and <code>parse_parameters()</code>
    // functions call the respective functions of the base classes as well.
    //
    // Note that this class also handles the declaration of initial and
    // boundary conditions specified in the input file. To this end, in both
    // cases, there are entries like "w_0 value" which represent an expression
    // in terms of $x,y,z$ that describe the initial or boundary condition as
    // a formula that will later be parsed by the FunctionParser
    // class. Similar expressions exist for "w_1", "w_2", etc, denoting the
    // <code>dim+2</code> conserved variables of the Euler system. Similarly,
    // we allow up to <code>max_n_boundaries</code> boundary indicators to be
    // used in the input file, and each of these boundary indicators can be
    // associated with an inflow, outflow, or pressure boundary condition,
    // with homogeneous boundary conditions being specified for each
    // component and each boundary indicator separately.
    //
    // The data structure used to store the boundary indicators is a bit
    // complicated. It is an array of <code>max_n_boundaries</code> elements
    // indicating the range of boundary indicators that will be accepted. For
    // each entry in this array, we store a pair of data in the
    // <code>BoundaryCondition</code> structure: first, an array of size
    // <code>n_components</code> that for each component of the solution
    // vector indicates whether it is an inflow, outflow, or other kind of
    // boundary, and second a FunctionParser object that describes all
    // components of the solution vector for this boundary id at once.
    //
    // The <code>BoundaryCondition</code> structure requires a constructor
    // since we need to tell the function parser object at construction time
    // how many vector components it is to describe. This initialization can
    // therefore not wait till we actually set the formulas the FunctionParser
    // object represents later in
    // <code>AllParameters::parse_parameters()</code>
    //
    // For the same reason of having to tell Function objects their vector
    // size at construction time, we have to have a constructor of the
    // <code>AllParameters</code> class that at least initializes the other
    // FunctionParser object, i.e. the one describing initial conditions.
    template <int dim>
    struct AllParameters : public Solver,
      public Refinement,
      public Flux,
      public Output
    {
      static const unsigned int max_n_boundaries = 10;

      struct BoundaryConditions
      {
        typename EulerEquations<dim>::BoundaryKind
        kind[EulerEquations<dim>::n_components];

        FunctionParser<dim> values;

        BoundaryConditions ();
      };


      AllParameters ();

      double diffusion_power;

      double time_step, final_time;
      double theta;
      bool is_stationary;

      std::string mesh_filename;

      FunctionParser<dim> initial_conditions;
      BoundaryConditions  boundary_conditions[max_n_boundaries];

      static void declare_parameters (ParameterHandler &prm);
      void parse_parameters (ParameterHandler &prm);
    };



    template <int dim>
    AllParameters<dim>::BoundaryConditions::BoundaryConditions ()
      :
      values (EulerEquations<dim>::n_components)
    {}


    template <int dim>
    AllParameters<dim>::AllParameters ()
      :
      initial_conditions (EulerEquations<dim>::n_components)
    {}


    template <int dim>
    void
    AllParameters<dim>::declare_parameters (ParameterHandler &prm)
    {
      prm.declare_entry("mesh", "grid.inp",
                        Patterns::Anything(),
                        "intput file name");

      prm.declare_entry("diffusion power", "2.0",
                        Patterns::Double(),
                        "power of mesh size for diffusion");

      prm.enter_subsection("time stepping");
      {
        prm.declare_entry("time step", "0.1",
                          Patterns::Double(0),
                          "simulation time step");
        prm.declare_entry("final time", "10.0",
                          Patterns::Double(0),
                          "simulation end time");
        prm.declare_entry("theta scheme value", "0.5",
                          Patterns::Double(0,1),
                          "value for theta that interpolated between explicit "
                          "Euler (theta=0), Crank-Nicolson (theta=0.5), and "
                          "implicit Euler (theta=1).");
      }
      prm.leave_subsection();


      for (unsigned int b=0; b<max_n_boundaries; ++b)
        {
          prm.enter_subsection("boundary_" +
                               Utilities::int_to_string(b));
          {
            prm.declare_entry("no penetration", "false",
                              Patterns::Bool(),
                              "whether the named boundary allows gas to "
                              "penetrate or is a rigid wall");

            for (unsigned int di=0; di<EulerEquations<dim>::n_components; ++di)
              {
                prm.declare_entry("w_" + Utilities::int_to_string(di),
                                  "outflow",
                                  Patterns::Selection("inflow|outflow|pressure"),
                                  "<inflow|outflow|pressure>");

                prm.declare_entry("w_" + Utilities::int_to_string(di) +
                                  " value", "0.0",
                                  Patterns::Anything(),
                                  "expression in x,y,z");
              }
          }
          prm.leave_subsection();
        }

      prm.enter_subsection("initial condition");
      {
        for (unsigned int di=0; di<EulerEquations<dim>::n_components; ++di)
          prm.declare_entry("w_" + Utilities::int_to_string(di) + " value",
                            "0.0",
                            Patterns::Anything(),
                            "expression in x,y,z");
      }
      prm.leave_subsection();

      Parameters::Solver::declare_parameters (prm);
      Parameters::Refinement::declare_parameters (prm);
      Parameters::Flux::declare_parameters (prm);
      Parameters::Output::declare_parameters (prm);
    }


    template <int dim>
    void
    AllParameters<dim>::parse_parameters (ParameterHandler &prm)
    {
      mesh_filename = prm.get("mesh");
      diffusion_power = prm.get_double("diffusion power");

      prm.enter_subsection("time stepping");
      {
        time_step = prm.get_double("time step");
        if (time_step == 0)
          {
            is_stationary = true;
            time_step = 1.0;
            final_time = 1.0;
          }
        else
          is_stationary = false;

        final_time = prm.get_double("final time");
        theta = prm.get_double("theta scheme value");
      }
      prm.leave_subsection();

      for (unsigned int boundary_id=0; boundary_id<max_n_boundaries;
           ++boundary_id)
        {
          prm.enter_subsection("boundary_" +
                               Utilities::int_to_string(boundary_id));
          {
            std::vector<std::string>
            expressions(EulerEquations<dim>::n_components, "0.0");

            const bool no_penetration = prm.get_bool("no penetration");

            for (unsigned int di=0; di<EulerEquations<dim>::n_components; ++di)
              {
                const std::string boundary_type
                  = prm.get("w_" + Utilities::int_to_string(di));

                if ((di < dim) && (no_penetration == true))
                  boundary_conditions[boundary_id].kind[di]
                    = EulerEquations<dim>::no_penetration_boundary;
                else if (boundary_type == "inflow")
                  boundary_conditions[boundary_id].kind[di]
                    = EulerEquations<dim>::inflow_boundary;
                else if (boundary_type == "pressure")
                  boundary_conditions[boundary_id].kind[di]
                    = EulerEquations<dim>::pressure_boundary;
                else if (boundary_type == "outflow")
                  boundary_conditions[boundary_id].kind[di]
                    = EulerEquations<dim>::outflow_boundary;
                else
                  AssertThrow (false, ExcNotImplemented());

                expressions[di] = prm.get("w_" + Utilities::int_to_string(di) +
                                          " value");
              }

            boundary_conditions[boundary_id].values
            .initialize (FunctionParser<dim>::default_variable_names(),
                         expressions,
                         std::map<std::string, double>());
          }
          prm.leave_subsection();
        }

      prm.enter_subsection("initial condition");
      {
        std::vector<std::string> expressions (EulerEquations<dim>::n_components,
                                              "0.0");
        for (unsigned int di = 0; di < EulerEquations<dim>::n_components; di++)
          expressions[di] = prm.get("w_" + Utilities::int_to_string(di) +
                                    " value");
        initial_conditions.initialize (FunctionParser<dim>::default_variable_names(),
                                       expressions,
                                       std::map<std::string, double>());
      }
      prm.leave_subsection();

      Parameters::Solver::parse_parameters (prm);
      Parameters::Refinement::parse_parameters (prm);
      Parameters::Flux::parse_parameters (prm);
      Parameters::Output::parse_parameters (prm);
    }
  }




  // @sect3{Conservation law class}

  // Here finally comes the class that actually does something with all the
  // Euler equation and parameter specifics we've defined above. The public
  // interface is pretty much the same as always (the constructor now takes
  // the name of a file from which to read parameters, which is passed on the
  // command line). The private function interface is also pretty similar to
  // the usual arrangement, with the <code>assemble_system</code> function
  // split into three parts: one that contains the main loop over all cells
  // and that then calls the other two for integrals over cells and faces,
  // respectively.
  template <int dim>
  class ConservationLaw
  {
  public:
    ConservationLaw (const char *input_filename);
    void run ();

  private:
    void setup_system ();

    void assemble_system ();
    void assemble_cell_term (const FEValues<dim>             &fe_v,
                             const std::vector<types::global_dof_index> &dofs);
    void assemble_face_term (const unsigned int               face_no,
                             const FEFaceValuesBase<dim>     &fe_v,
                             const FEFaceValuesBase<dim>     &fe_v_neighbor,
                             const std::vector<types::global_dof_index> &dofs,
                             const std::vector<types::global_dof_index> &dofs_neighbor,
                             const bool                       external_face,
                             const unsigned int               boundary_id,
                             const double                     face_diameter);

    std::pair<unsigned int, double> solve (Vector<double> &solution);

    void compute_refinement_indicators (Vector<double> &indicator) const;
    void refine_grid (const Vector<double> &indicator);

    void output_results () const;



    // The first few member variables are also rather standard. Note that we
    // define a mapping object to be used throughout the program when
    // assembling terms (we will hand it to every FEValues and FEFaceValues
    // object); the mapping we use is just the standard $Q_1$ mapping --
    // nothing fancy, in other words -- but declaring one here and using it
    // throughout the program will make it simpler later on to change it if
    // that should become necessary. This is, in fact, rather pertinent: it is
    // known that for transsonic simulations with the Euler equations,
    // computations do not converge even as $h\rightarrow 0$ if the boundary
    // approximation is not of sufficiently high order.
    Triangulation<dim>   triangulation;
    const MappingQ1<dim> mapping;

    const FESystem<dim>  fe;
    DoFHandler<dim>      dof_handler;

    const QGauss<dim>    quadrature;
    const QGauss<dim-1>  face_quadrature;

    // Next come a number of data vectors that correspond to the solution of
    // the previous time step (<code>old_solution</code>), the best guess of
    // the current solution (<code>current_solution</code>; we say
    // <i>guess</i> because the Newton iteration to compute it may not have
    // converged yet, whereas <code>old_solution</code> refers to the fully
    // converged final result of the previous time step), and a predictor for
    // the solution at the next time step, computed by extrapolating the
    // current and previous solution one time step into the future:
    Vector<double>       old_solution;
    Vector<double>       current_solution;
    Vector<double>       predictor;

    Vector<double>       right_hand_side;

    // This final set of member variables (except for the object holding all
    // run-time parameters at the very bottom and a screen output stream that
    // only prints something if verbose output has been requested) deals with
    // the interface we have in this program to the Trilinos library that
    // provides us with linear solvers. Similarly to including PETSc matrices
    // in step-17, step-18, and step-19, all we need to do is to create a
    // Trilinos sparse matrix instead of the standard deal.II class. The
    // system matrix is used for the Jacobian in each Newton step. Since we do
    // not intend to run this program in parallel (which wouldn't be too hard
    // with Trilinos data structures, though), we don't have to think about
    // anything else like distributing the degrees of freedom.
    TrilinosWrappers::SparseMatrix system_matrix;

    Parameters::AllParameters<dim>  parameters;
    ConditionalOStream              verbose_cout;
  };


  // @sect4{ConservationLaw::ConservationLaw}
  //
  // There is nothing much to say about the constructor. Essentially, it reads
  // the input file and fills the parameter object with the parsed values:
  template <int dim>
  ConservationLaw<dim>::ConservationLaw (const char *input_filename)
    :
    mapping (),
    fe (FE_Q<dim>(1), EulerEquations<dim>::n_components),
    dof_handler (triangulation),
    quadrature (2),
    face_quadrature (2),
    verbose_cout (std::cout, false)
  {
    ParameterHandler prm;
    Parameters::AllParameters<dim>::declare_parameters (prm);

    prm.read_input (input_filename);
    parameters.parse_parameters (prm);

    verbose_cout.set_condition (parameters.output == Parameters::Solver::verbose);
  }



  // @sect4{ConservationLaw::setup_system}
  //
  // The following (easy) function is called each time the mesh is
  // changed. All it does is to resize the Trilinos matrix according to a
  // sparsity pattern that we generate as in all the previous tutorial
  // programs.
  template <int dim>
  void ConservationLaw<dim>::setup_system ()
  {
    CompressedSparsityPattern sparsity_pattern (dof_handler.n_dofs(),
                                                dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);

    system_matrix.reinit (sparsity_pattern);
  }


  // @sect4{ConservationLaw::assemble_system}
  //
  // This and the following two functions are the meat of this program: They
  // assemble the linear system that results from applying Newton's method to
  // the nonlinear system of conservation equations.
  //
  // This first function puts all of the assembly pieces together in a routine
  // that dispatches the correct piece for each cell/face.  The actual
  // implementation of the assembly on these objects is done in the following
  // functions.
  //
  // At the top of the function we do the usual housekeeping: allocate
  // FEValues, FEFaceValues, and FESubfaceValues objects necessary to do the
  // integrations on cells, faces, and subfaces (in case of adjoining cells on
  // different refinement levels). Note that we don't need all information
  // (like values, gradients, or real locations of quadrature points) for all
  // of these objects, so we only let the FEValues classes whatever is
  // actually necessary by specifying the minimal set of UpdateFlags. For
  // example, when using a FEFaceValues object for the neighboring cell we
  // only need the shape values: Given a specific face, the quadrature points
  // and <code>JxW</code> values are the same as for the current cells, and
  // the normal vectors are known to be the negative of the normal vectors of
  // the current cell.
  template <int dim>
  void ConservationLaw<dim>::assemble_system ()
  {
    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

    std::vector<types::global_dof_index> dof_indices (dofs_per_cell);
    std::vector<types::global_dof_index> dof_indices_neighbor (dofs_per_cell);

    const UpdateFlags update_flags               = update_values
                                                   | update_gradients
                                                   | update_q_points
                                                   | update_JxW_values,
                                                   face_update_flags          = update_values
                                                       | update_q_points
                                                       | update_JxW_values
                                                       | update_normal_vectors,
                                                       neighbor_face_update_flags = update_values;

    FEValues<dim>        fe_v                  (mapping, fe, quadrature,
                                                update_flags);
    FEFaceValues<dim>    fe_v_face             (mapping, fe, face_quadrature,
                                                face_update_flags);
    FESubfaceValues<dim> fe_v_subface          (mapping, fe, face_quadrature,
                                                face_update_flags);
    FEFaceValues<dim>    fe_v_face_neighbor    (mapping, fe, face_quadrature,
                                                neighbor_face_update_flags);
    FESubfaceValues<dim> fe_v_subface_neighbor (mapping, fe, face_quadrature,
                                                neighbor_face_update_flags);

    // Then loop over all cells, initialize the FEValues object for the
    // current cell and call the function that assembles the problem on this
    // cell.
    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        fe_v.reinit (cell);
        cell->get_dof_indices (dof_indices);

        assemble_cell_term(fe_v, dof_indices);

        // Then loop over all the faces of this cell.  If a face is part of
        // the external boundary, then assemble boundary conditions there (the
        // fifth argument to <code>assemble_face_terms</code> indicates
        // whether we are working on an external or internal face; if it is an
        // external face, the fourth argument denoting the degrees of freedom
        // indices of the neighbor is ignored, so we pass an empty vector):
        for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell;
             ++face_no)
          if (cell->at_boundary(face_no))
            {
              fe_v_face.reinit (cell, face_no);
              assemble_face_term (face_no, fe_v_face,
                                  fe_v_face,
                                  dof_indices,
                                  std::vector<types::global_dof_index>(),
                                  true,
                                  cell->face(face_no)->boundary_indicator(),
                                  cell->face(face_no)->diameter());
            }

        // The alternative is that we are dealing with an internal face. There
        // are two cases that we need to distinguish: that this is a normal
        // face between two cells at the same refinement level, and that it is
        // a face between two cells of the different refinement levels.
        //
        // In the first case, there is nothing we need to do: we are using a
        // continuous finite element, and face terms do not appear in the
        // bilinear form in this case. The second case usually does not lead
        // to face terms either if we enforce hanging node constraints
        // strongly (as in all previous tutorial programs so far whenever we
        // used continuous finite elements -- this enforcement is done by the
        // ConstraintMatrix class together with
        // DoFTools::make_hanging_node_constraints). In the current program,
        // however, we opt to enforce continuity weakly at faces between cells
        // of different refinement level, for two reasons: (i) because we can,
        // and more importantly (ii) because we would have to thread the
        // automatic differentiation we use to compute the elements of the
        // Newton matrix from the residual through the operations of the
        // ConstraintMatrix class. This would be possible, but is not trivial,
        // and so we choose this alternative approach.
        //
        // What needs to be decided is which side of an interface between two
        // cells of different refinement level we are sitting on.
        //
        // Let's take the case where the neighbor is more refined first. We
        // then have to loop over the children of the face of the current cell
        // and integrate on each of them. We sprinkle a couple of assertions
        // into the code to ensure that our reasoning trying to figure out
        // which of the neighbor's children's faces coincides with a given
        // subface of the current cell's faces is correct -- a bit of
        // defensive programming never hurts.
        //
        // We then call the function that integrates over faces; since this is
        // an internal face, the fifth argument is false, and the sixth one is
        // ignored so we pass an invalid value again:
          else
            {
              if (cell->neighbor(face_no)->has_children())
                {
                  const unsigned int neighbor2=
                    cell->neighbor_of_neighbor(face_no);

                  for (unsigned int subface_no=0;
                       subface_no < cell->face(face_no)->n_children();
                       ++subface_no)
                    {
                      const typename DoFHandler<dim>::active_cell_iterator
                      neighbor_child
                        = cell->neighbor_child_on_subface (face_no, subface_no);

                      Assert (neighbor_child->face(neighbor2) ==
                              cell->face(face_no)->child(subface_no),
                              ExcInternalError());
                      Assert (neighbor_child->has_children() == false,
                              ExcInternalError());

                      fe_v_subface.reinit (cell, face_no, subface_no);
                      fe_v_face_neighbor.reinit (neighbor_child, neighbor2);

                      neighbor_child->get_dof_indices (dof_indices_neighbor);

                      assemble_face_term (face_no, fe_v_subface,
                                          fe_v_face_neighbor,
                                          dof_indices,
                                          dof_indices_neighbor,
                                          false,
                                          numbers::invalid_unsigned_int,
                                          neighbor_child->face(neighbor2)->diameter());
                    }
                }

              // The other possibility we have to care for is if the neighbor
              // is coarser than the current cell (in particular, because of
              // the usual restriction of only one hanging node per face, the
              // neighbor must be exactly one level coarser than the current
              // cell, something that we check with an assertion). Again, we
              // then integrate over this interface:
              else if (cell->neighbor(face_no)->level() != cell->level())
                {
                  const typename DoFHandler<dim>::cell_iterator
                  neighbor = cell->neighbor(face_no);
                  Assert(neighbor->level() == cell->level()-1,
                         ExcInternalError());

                  neighbor->get_dof_indices (dof_indices_neighbor);

                  const std::pair<unsigned int, unsigned int>
                  faceno_subfaceno = cell->neighbor_of_coarser_neighbor(face_no);
                  const unsigned int neighbor_face_no    = faceno_subfaceno.first,
                                     neighbor_subface_no = faceno_subfaceno.second;

                  Assert (neighbor->neighbor_child_on_subface (neighbor_face_no,
                                                               neighbor_subface_no)
                          == cell,
                          ExcInternalError());

                  fe_v_face.reinit (cell, face_no);
                  fe_v_subface_neighbor.reinit (neighbor,
                                                neighbor_face_no,
                                                neighbor_subface_no);

                  assemble_face_term (face_no, fe_v_face,
                                      fe_v_subface_neighbor,
                                      dof_indices,
                                      dof_indices_neighbor,
                                      false,
                                      numbers::invalid_unsigned_int,
                                      cell->face(face_no)->diameter());
                }
            }
      }

    // After all this assembling, notify the Trilinos matrix object that the
    // matrix is done:
    system_matrix.compress(VectorOperation::add);
  }


  // @sect4{ConservationLaw::assemble_cell_term}
  //
  // This function assembles the cell term by computing the cell part of the
  // residual, adding its negative to the right hand side vector, and adding
  // its derivative with respect to the local variables to the Jacobian
  // (i.e. the Newton matrix). Recall that the cell contributions to the
  // residual read $F_i = \left(\frac{\mathbf{w}_{n+1} - \mathbf{w}_n}{\delta
  // t},\mathbf{z}_i\right)_K - \left(\mathbf{F}(\tilde{\mathbf{w}}),
  // \nabla\mathbf{z}_i\right)_K + h^{\eta}(\nabla \mathbf{w} , \nabla
  // \mathbf{z}_i)_K - (\mathbf{G}(\tilde{\mathbf w}), \mathbf{z}_i)_K$ where
  // $\tilde{\mathbf w}$ is represented by the variable <code>W_theta</code>,
  // $\mathbf{z}_i$ is the $i$th test function, and the scalar product
  // $\left(\mathbf{F}(\tilde{\mathbf{w}}), \nabla\mathbf{z}\right)_K$ is
  // understood as $\int_K \sum_{c=1}^{\text{n\_components}}
  // \sum_{d=1}^{\text{dim}} \mathbf{F}(\tilde{\mathbf{w}})_{cd}
  // \frac{\partial z_c}{x_d}$.
  //
  // At the top of this function, we do the usual housekeeping in terms of
  // allocating some local variables that we will need later. In particular,
  // we will allocate variables that will hold the values of the current
  // solution $W_{n+1}^k$ after the $k$th Newton iteration (variable
  // <code>W</code>), the previous time step's solution $W_{n}$ (variable
  // <code>W_old</code>), as well as the linear combination $\theta W_{n+1}^k
  // + (1-\theta)W_n$ that results from choosing different time stepping
  // schemes (variable <code>W_theta</code>).
  //
  // In addition to these, we need the gradients of the current variables.  It
  // is a bit of a shame that we have to compute these; we almost don't.  The
  // nice thing about a simple conservation law is that the flux doesn't
  // generally involve any gradients.  We do need these, however, for the
  // diffusion stabilization.
  //
  // The actual format in which we store these variables requires some
  // explanation. First, we need values at each quadrature point for each of
  // the <code>EulerEquations::n_components</code> components of the solution
  // vector. This makes for a two-dimensional table for which we use deal.II's
  // Table class (this is more efficient than
  // <code>std::vector@<std::vector@<T@> @></code> because it only needs to
  // allocate memory once, rather than once for each element of the outer
  // vector). Similarly, the gradient is a three-dimensional table, which the
  // Table class also supports.
  //
  // Secondly, we want to use automatic differentiation. To this end, we use
  // the Sacado::Fad::DFad template for everything that is a computed from the
  // variables with respect to which we would like to compute
  // derivatives. This includes the current solution and gradient at the
  // quadrature points (which are linear combinations of the degrees of
  // freedom) as well as everything that is computed from them such as the
  // residual, but not the previous time step's solution. These variables are
  // all found in the first part of the function, along with a variable that
  // we will use to store the derivatives of a single component of the
  // residual:
  template <int dim>
  void
  ConservationLaw<dim>::
  assemble_cell_term (const FEValues<dim>             &fe_v,
                      const std::vector<types::global_dof_index> &dof_indices)
  {
    const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
    const unsigned int n_q_points    = fe_v.n_quadrature_points;

    Table<2,Sacado::Fad::DFad<double> >
    W (n_q_points, EulerEquations<dim>::n_components);

    Table<2,double>
    W_old (n_q_points, EulerEquations<dim>::n_components);

    Table<2,Sacado::Fad::DFad<double> >
    W_theta (n_q_points, EulerEquations<dim>::n_components);

    Table<3,Sacado::Fad::DFad<double> >
    grad_W (n_q_points, EulerEquations<dim>::n_components, dim);

    std::vector<double> residual_derivatives (dofs_per_cell);

    // Next, we have to define the independent variables that we will try to
    // determine by solving a Newton step. These independent variables are the
    // values of the local degrees of freedom which we extract here:
    std::vector<Sacado::Fad::DFad<double> > independent_local_dof_values(dofs_per_cell);
    for (unsigned int i=0; i<dofs_per_cell; ++i)
      independent_local_dof_values[i] = current_solution(dof_indices[i]);

    // The next step incorporates all the magic: we declare a subset of the
    // autodifferentiation variables as independent degrees of freedom,
    // whereas all the other ones remain dependent functions. These are
    // precisely the local degrees of freedom just extracted. All calculations
    // that reference them (either directly or indirectly) will accumulate
    // sensitivities with respect to these variables.
    //
    // In order to mark the variables as independent, the following does the
    // trick, marking <code>independent_local_dof_values[i]</code> as the
    // $i$th independent variable out of a total of
    // <code>dofs_per_cell</code>:
    for (unsigned int i=0; i<dofs_per_cell; ++i)
      independent_local_dof_values[i].diff (i, dofs_per_cell);

    // After all these declarations, let us actually compute something. First,
    // the values of <code>W</code>, <code>W_old</code>, <code>W_theta</code>,
    // and <code>grad_W</code>, which we can compute from the local DoF values
    // by using the formula $W(x_q)=\sum_i \mathbf W_i \Phi_i(x_q)$, where
    // $\mathbf W_i$ is the $i$th entry of the (local part of the) solution
    // vector, and $\Phi_i(x_q)$ the value of the $i$th vector-valued shape
    // function evaluated at quadrature point $x_q$. The gradient can be
    // computed in a similar way.
    //
    // Ideally, we could compute this information using a call into something
    // like FEValues::get_function_values and FEValues::get_function_grads,
    // but since (i) we would have to extend the FEValues class for this, and
    // (ii) we don't want to make the entire <code>old_solution</code> vector
    // fad types, only the local cell variables, we explicitly code the loop
    // above. Before this, we add another loop that initializes all the fad
    // variables to zero:
    for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned int c=0; c<EulerEquations<dim>::n_components; ++c)
        {
          W[q][c]       = 0;
          W_old[q][c]   = 0;
          W_theta[q][c] = 0;
          for (unsigned int d=0; d<dim; ++d)
            grad_W[q][c][d] = 0;
        }

    for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          const unsigned int c = fe_v.get_fe().system_to_component_index(i).first;

          W[q][c] += independent_local_dof_values[i] *
                     fe_v.shape_value_component(i, q, c);
          W_old[q][c] += old_solution(dof_indices[i]) *
                         fe_v.shape_value_component(i, q, c);
          W_theta[q][c] += (parameters.theta *
                            independent_local_dof_values[i]
                            +
                            (1-parameters.theta) *
                            old_solution(dof_indices[i])) *
                           fe_v.shape_value_component(i, q, c);

          for (unsigned int d = 0; d < dim; d++)
            grad_W[q][c][d] += independent_local_dof_values[i] *
                               fe_v.shape_grad_component(i, q, c)[d];
        }


    // Next, in order to compute the cell contributions, we need to evaluate
    // $F(\tilde{\mathbf w})$ and $G(\tilde{\mathbf w})$ at all quadrature
    // points. To store these, we also need to allocate a bit of memory. Note
    // that we compute the flux matrices and right hand sides in terms of
    // autodifferentiation variables, so that the Jacobian contributions can
    // later easily be computed from it:
    typedef Sacado::Fad::DFad<double> FluxMatrix[EulerEquations<dim>::n_components][dim];
    FluxMatrix *flux = new FluxMatrix[n_q_points];

    typedef Sacado::Fad::DFad<double> ForcingVector[EulerEquations<dim>::n_components];
    ForcingVector *forcing = new ForcingVector[n_q_points];

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        EulerEquations<dim>::compute_flux_matrix (W_theta[q], flux[q]);
        EulerEquations<dim>::compute_forcing_vector (W_theta[q], forcing[q]);
      }


    // We now have all of the pieces in place, so perform the assembly.  We
    // have an outer loop through the components of the system, and an inner
    // loop over the quadrature points, where we accumulate contributions to
    // the $i$th residual $F_i$. The general formula for this residual is
    // given in the introduction and at the top of this function. We can,
    // however, simplify it a bit taking into account that the $i$th
    // (vector-valued) test function $\mathbf{z}_i$ has in reality only a
    // single nonzero component (more on this topic can be found in the @ref
    // vector_valued module). It will be represented by the variable
    // <code>component_i</code> below. With this, the residual term can be
    // re-written as $F_i = \left(\frac{(\mathbf{w}_{n+1} -
    // \mathbf{w}_n)_{\text{component\_i}}}{\delta
    // t},(\mathbf{z}_i)_{\text{component\_i}}\right)_K$ $-
    // \sum_{d=1}^{\text{dim}} \left(\mathbf{F}
    // (\tilde{\mathbf{w}})_{\text{component\_i},d},
    // \frac{\partial(\mathbf{z}_i)_{\text{component\_i}}} {\partial
    // x_d}\right)_K$ $+ \sum_{d=1}^{\text{dim}} h^{\eta} \left(\frac{\partial
    // \mathbf{w}_{\text{component\_i}}}{\partial x_d} , \frac{\partial
    // (\mathbf{z}_i)_{\text{component\_i}}}{\partial x_d} \right)_K$
    // $-(\mathbf{G}(\tilde{\mathbf{w}} )_{\text{component\_i}},
    // (\mathbf{z}_i)_{\text{component\_i}})_K$, where integrals are
    // understood to be evaluated through summation over quadrature points.
    //
    // We initially sum all contributions of the residual in the positive
    // sense, so that we don't need to negative the Jacobian entries.  Then,
    // when we sum into the <code>right_hand_side</code> vector, we negate
    // this residual.
    for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
      {
        Sacado::Fad::DFad<double> F_i = 0;

        const unsigned int
        component_i = fe_v.get_fe().system_to_component_index(i).first;

        // The residual for each row (i) will be accumulating into this fad
        // variable.  At the end of the assembly for this row, we will query
        // for the sensitivities to this variable and add them into the
        // Jacobian.

        for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
          {
            if (parameters.is_stationary == false)
              F_i += 1.0 / parameters.time_step *
                     (W[point][component_i] - W_old[point][component_i]) *
                     fe_v.shape_value_component(i, point, component_i) *
                     fe_v.JxW(point);

            for (unsigned int d=0; d<dim; d++)
              F_i -= flux[point][component_i][d] *
                     fe_v.shape_grad_component(i, point, component_i)[d] *
                     fe_v.JxW(point);

            for (unsigned int d=0; d<dim; d++)
              F_i += 1.0*std::pow(fe_v.get_cell()->diameter(),
                                  parameters.diffusion_power) *
                     grad_W[point][component_i][d] *
                     fe_v.shape_grad_component(i, point, component_i)[d] *
                     fe_v.JxW(point);

            F_i -= forcing[point][component_i] *
                   fe_v.shape_value_component(i, point, component_i) *
                   fe_v.JxW(point);
          }

        // At the end of the loop, we have to add the sensitivities to the
        // matrix and subtract the residual from the right hand side. Trilinos
        // FAD data type gives us access to the derivatives using
        // <code>F_i.fastAccessDx(k)</code>, so we store the data in a
        // temporary array. This information about the whole row of local dofs
        // is then added to the Trilinos matrix at once (which supports the
        // data types we have chosen).
        for (unsigned int k=0; k<dofs_per_cell; ++k)
          residual_derivatives[k] = F_i.fastAccessDx(k);
        system_matrix.add(dof_indices[i], dof_indices, residual_derivatives);
        right_hand_side(dof_indices[i]) -= F_i.val();
      }

    delete[] forcing;
    delete[] flux;
  }


  // @sect4{ConservationLaw::assemble_face_term}
  //
  // Here, we do essentially the same as in the previous function. t the top,
  // we introduce the independent variables. Because the current function is
  // also used if we are working on an internal face between two cells, the
  // independent variables are not only the degrees of freedom on the current
  // cell but in the case of an interior face also the ones on the neighbor.
  template <int dim>
  void
  ConservationLaw<dim>::assemble_face_term(const unsigned int           face_no,
                                           const FEFaceValuesBase<dim> &fe_v,
                                           const FEFaceValuesBase<dim> &fe_v_neighbor,
                                           const std::vector<types::global_dof_index>   &dof_indices,
                                           const std::vector<types::global_dof_index>   &dof_indices_neighbor,
                                           const bool                   external_face,
                                           const unsigned int           boundary_id,
                                           const double                 face_diameter)
  {
    const unsigned int n_q_points = fe_v.n_quadrature_points;
    const unsigned int dofs_per_cell = fe_v.dofs_per_cell;

    std::vector<Sacado::Fad::DFad<double> >
    independent_local_dof_values (dofs_per_cell),
                                 independent_neighbor_dof_values (external_face == false ?
                                     dofs_per_cell :
                                     0);

    const unsigned int n_independent_variables = (external_face == false ?
                                                  2 * dofs_per_cell :
                                                  dofs_per_cell);

    for (unsigned int i = 0; i < dofs_per_cell; i++)
      {
        independent_local_dof_values[i] = current_solution(dof_indices[i]);
        independent_local_dof_values[i].diff(i, n_independent_variables);
      }

    if (external_face == false)
      for (unsigned int i = 0; i < dofs_per_cell; i++)
        {
          independent_neighbor_dof_values[i]
            = current_solution(dof_indices_neighbor[i]);
          independent_neighbor_dof_values[i]
          .diff(i+dofs_per_cell, n_independent_variables);
        }


    // Next, we need to define the values of the conservative variables
    // $\tilde {\mathbf W}$ on this side of the face ($\tilde {\mathbf W}^+$)
    // and on the opposite side ($\tilde {\mathbf W}^-$). The former can be
    // computed in exactly the same way as in the previous function, but note
    // that the <code>fe_v</code> variable now is of type FEFaceValues or
    // FESubfaceValues:
    Table<2,Sacado::Fad::DFad<double> >
    Wplus (n_q_points, EulerEquations<dim>::n_components),
          Wminus (n_q_points, EulerEquations<dim>::n_components);

    for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          const unsigned int component_i = fe_v.get_fe().system_to_component_index(i).first;
          Wplus[q][component_i] += (parameters.theta *
                                    independent_local_dof_values[i]
                                    +
                                    (1.0-parameters.theta) *
                                    old_solution(dof_indices[i])) *
                                   fe_v.shape_value_component(i, q, component_i);
        }

    // Computing $\tilde {\mathbf W}^-$ is a bit more complicated. If this is
    // an internal face, we can compute it as above by simply using the
    // independent variables from the neighbor:
    if (external_face == false)
      {
        for (unsigned int q=0; q<n_q_points; ++q)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              const unsigned int component_i = fe_v_neighbor.get_fe().
                                               system_to_component_index(i).first;
              Wminus[q][component_i] += (parameters.theta *
                                         independent_neighbor_dof_values[i]
                                         +
                                         (1.0-parameters.theta) *
                                         old_solution(dof_indices_neighbor[i]))*
                                        fe_v_neighbor.shape_value_component(i, q, component_i);
            }
      }
    // On the other hand, if this is an external boundary face, then the
    // values of $W^-$ will be either functions of $W^+$, or they will be
    // prescribed, depending on the kind of boundary condition imposed here.
    //
    // To start the evaluation, let us ensure that the boundary id specified
    // for this boundary is one for which we actually have data in the
    // parameters object. Next, we evaluate the function object for the
    // inhomogeneity.  This is a bit tricky: a given boundary might have both
    // prescribed and implicit values.  If a particular component is not
    // prescribed, the values evaluate to zero and are ignored below.
    //
    // The rest is done by a function that actually knows the specifics of
    // Euler equation boundary conditions. Note that since we are using fad
    // variables here, sensitivities will be updated appropriately, a process
    // that would otherwise be tremendously complicated.
    else
      {
        Assert (boundary_id < Parameters::AllParameters<dim>::max_n_boundaries,
                ExcIndexRange (boundary_id, 0,
                               Parameters::AllParameters<dim>::max_n_boundaries));

        std::vector<Vector<double> >
        boundary_values(n_q_points, Vector<double>(EulerEquations<dim>::n_components));
        parameters.boundary_conditions[boundary_id]
        .values.vector_value_list(fe_v.get_quadrature_points(),
                                  boundary_values);

        for (unsigned int q = 0; q < n_q_points; q++)
          EulerEquations<dim>::compute_Wminus (parameters.boundary_conditions[boundary_id].kind,
                                               fe_v.normal_vector(q),
                                               Wplus[q],
                                               boundary_values[q],
                                               Wminus[q]);
      }


    // Now that we have $\mathbf w^+$ and $\mathbf w^-$, we can go about
    // computing the numerical flux function $\mathbf H(\mathbf w^+,\mathbf
    // w^-, \mathbf n)$ for each quadrature point. Before calling the function
    // that does so, we also need to determine the Lax-Friedrich's stability
    // parameter:
    typedef Sacado::Fad::DFad<double> NormalFlux[EulerEquations<dim>::n_components];
    NormalFlux *normal_fluxes = new NormalFlux[n_q_points];

    double alpha;

    switch (parameters.stabilization_kind)
      {
      case Parameters::Flux::constant:
        alpha = parameters.stabilization_value;
        break;
      case Parameters::Flux::mesh_dependent:
        alpha = face_diameter/(2.0*parameters.time_step);
        break;
      default:
        Assert (false, ExcNotImplemented());
        alpha = 1;
      }

    for (unsigned int q=0; q<n_q_points; ++q)
      EulerEquations<dim>::numerical_normal_flux(fe_v.normal_vector(q),
                                                 Wplus[q], Wminus[q], alpha,
                                                 normal_fluxes[q]);

    // Now assemble the face term in exactly the same way as for the cell
    // contributions in the previous function. The only difference is that if
    // this is an internal face, we also have to take into account the
    // sensitivities of the residual contributions to the degrees of freedom on
    // the neighboring cell:
    std::vector<double> residual_derivatives (dofs_per_cell);
    for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
      if (fe_v.get_fe().has_support_on_face(i, face_no) == true)
        {
          Sacado::Fad::DFad<double> F_i = 0;

          for (unsigned int point=0; point<n_q_points; ++point)
            {
              const unsigned int
              component_i = fe_v.get_fe().system_to_component_index(i).first;

              F_i += normal_fluxes[point][component_i] *
                     fe_v.shape_value_component(i, point, component_i) *
                     fe_v.JxW(point);
            }

          for (unsigned int k=0; k<dofs_per_cell; ++k)
            residual_derivatives[k] = F_i.fastAccessDx(k);
          system_matrix.add(dof_indices[i], dof_indices, residual_derivatives);

          if (external_face == false)
            {
              for (unsigned int k=0; k<dofs_per_cell; ++k)
                residual_derivatives[k] = F_i.fastAccessDx(dofs_per_cell+k);
              system_matrix.add (dof_indices[i], dof_indices_neighbor,
                                 residual_derivatives);
            }

          right_hand_side(dof_indices[i]) -= F_i.val();
        }

    delete[] normal_fluxes;
  }


  // @sect4{ConservationLaw::solve}
  //
  // Here, we actually solve the linear system, using either of Trilinos'
  // Aztec or Amesos linear solvers. The result of the computation will be
  // written into the argument vector passed to this function. The result is a
  // pair of number of iterations and the final linear residual.

  template <int dim>
  std::pair<unsigned int, double>
  ConservationLaw<dim>::solve (Vector<double> &newton_update)
  {
    switch (parameters.solver)
      {
      // If the parameter file specified that a direct solver shall be used,
      // then we'll get here. The process is straightforward, since deal.II
      // provides a wrapper class to the Amesos direct solver within
      // Trilinos. All we have to do is to create a solver control object
      // (which is just a dummy object here, since we won't perform any
      // iterations), and then create the direct solver object. When
      // actually doing the solve, note that we don't pass a
      // preconditioner. That wouldn't make much sense for a direct solver
      // anyway.  At the end we return the solver control statistics &mdash;
      // which will tell that no iterations have been performed and that the
      // final linear residual is zero, absent any better information that
      // may be provided here:
      case Parameters::Solver::direct:
      {
        SolverControl solver_control (1,0);
        TrilinosWrappers::SolverDirect direct (solver_control,
                                               parameters.output ==
                                               Parameters::Solver::verbose);

        direct.solve (system_matrix, newton_update, right_hand_side);

        return std::pair<unsigned int, double> (solver_control.last_step(),
                                                solver_control.last_value());
      }

      // Likewise, if we are to use an iterative solver, we use Aztec's GMRES
      // solver. We could use the Trilinos wrapper classes for iterative
      // solvers and preconditioners here as well, but we choose to use an
      // Aztec solver directly. For the given problem, Aztec's internal
      // preconditioner implementations are superior over the ones deal.II has
      // wrapper classes to, so we use ILU-T preconditioning within the
      // AztecOO solver and set a bunch of options that can be changed from
      // the parameter file.
      //
      // There are two more practicalities: Since we have built our right hand
      // side and solution vector as deal.II Vector objects (as opposed to the
      // matrix, which is a Trilinos object), we must hand the solvers
      // Trilinos Epetra vectors.  Luckily, they support the concept of a
      // 'view', so we just send in a pointer to our deal.II vectors. We have
      // to provide an Epetra_Map for the vector that sets the parallel
      // distribution, which is just a dummy object in serial. The easiest way
      // is to ask the matrix for its map, and we're going to be ready for
      // matrix-vector products with it.
      //
      // Secondly, the Aztec solver wants us to pass a Trilinos
      // Epetra_CrsMatrix in, not the deal.II wrapper class itself. So we
      // access to the actual Trilinos matrix in the Trilinos wrapper class by
      // the command trilinos_matrix(). Trilinos wants the matrix to be
      // non-constant, so we have to manually remove the constantness using a
      // const_cast.
      case Parameters::Solver::gmres:
      {
        Epetra_Vector x(View, system_matrix.domain_partitioner(),
                        newton_update.begin());
        Epetra_Vector b(View, system_matrix.range_partitioner(),
                        right_hand_side.begin());

        AztecOO solver;
        solver.SetAztecOption(AZ_output,
                              (parameters.output ==
                               Parameters::Solver::quiet
                               ?
                               AZ_none
                               :
                               AZ_all));
        solver.SetAztecOption(AZ_solver, AZ_gmres);
        solver.SetRHS(&b);
        solver.SetLHS(&x);

        solver.SetAztecOption(AZ_precond,         AZ_dom_decomp);
        solver.SetAztecOption(AZ_subdomain_solve, AZ_ilut);
        solver.SetAztecOption(AZ_overlap,         0);
        solver.SetAztecOption(AZ_reorder,         0);

        solver.SetAztecParam(AZ_drop,      parameters.ilut_drop);
        solver.SetAztecParam(AZ_ilut_fill, parameters.ilut_fill);
        solver.SetAztecParam(AZ_athresh,   parameters.ilut_atol);
        solver.SetAztecParam(AZ_rthresh,   parameters.ilut_rtol);

        solver.SetUserMatrix(const_cast<Epetra_CrsMatrix *>
                             (&system_matrix.trilinos_matrix()));

        solver.Iterate(parameters.max_iterations, parameters.linear_residual);

        return std::pair<unsigned int, double> (solver.NumIters(),
                                                solver.TrueResidual());
      }
      }

    Assert (false, ExcNotImplemented());
    return std::pair<unsigned int, double> (0,0);
  }


  // @sect4{ConservationLaw::compute_refinement_indicators}

  // This function is real simple: We don't pretend that we know here what a
  // good refinement indicator would be. Rather, we assume that the
  // <code>EulerEquation</code> class would know about this, and so we simply
  // defer to the respective function we've implemented there:
  template <int dim>
  void
  ConservationLaw<dim>::
  compute_refinement_indicators (Vector<double> &refinement_indicators) const
  {
    EulerEquations<dim>::compute_refinement_indicators (dof_handler,
                                                        mapping,
                                                        predictor,
                                                        refinement_indicators);
  }



  // @sect4{ConservationLaw::refine_grid}

  // Here, we use the refinement indicators computed before and refine the
  // mesh. At the beginning, we loop over all cells and mark those that we
  // think should be refined:
  template <int dim>
  void
  ConservationLaw<dim>::refine_grid (const Vector<double> &refinement_indicators)
  {
    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    for (unsigned int cell_no=0; cell!=endc; ++cell, ++cell_no)
      {
        cell->clear_coarsen_flag();
        cell->clear_refine_flag();

        if ((cell->level() < parameters.shock_levels) &&
            (std::fabs(refinement_indicators(cell_no)) > parameters.shock_val))
          cell->set_refine_flag();
        else if ((cell->level() > 0) &&
                 (std::fabs(refinement_indicators(cell_no)) < 0.75*parameters.shock_val))
          cell->set_coarsen_flag();
      }

    // Then we need to transfer the various solution vectors from the old to
    // the new grid while we do the refinement. The SolutionTransfer class is
    // our friend here; it has a fairly extensive documentation, including
    // examples, so we won't comment much on the following code. The last
    // three lines simply re-set the sizes of some other vectors to the now
    // correct size:
    std::vector<Vector<double> > transfer_in;
    std::vector<Vector<double> > transfer_out;

    transfer_in.push_back(old_solution);
    transfer_in.push_back(predictor);

    triangulation.prepare_coarsening_and_refinement();

    SolutionTransfer<dim> soltrans(dof_handler);
    soltrans.prepare_for_coarsening_and_refinement(transfer_in);

    triangulation.execute_coarsening_and_refinement ();

    dof_handler.clear();
    dof_handler.distribute_dofs (fe);

    {
      Vector<double> new_old_solution(1);
      Vector<double> new_predictor(1);

      transfer_out.push_back(new_old_solution);
      transfer_out.push_back(new_predictor);
      transfer_out[0].reinit(dof_handler.n_dofs());
      transfer_out[1].reinit(dof_handler.n_dofs());
    }

    soltrans.interpolate(transfer_in, transfer_out);

    old_solution.reinit (transfer_out[0].size());
    old_solution = transfer_out[0];

    predictor.reinit (transfer_out[1].size());
    predictor = transfer_out[1];

    current_solution.reinit(dof_handler.n_dofs());
    current_solution = old_solution;
    right_hand_side.reinit (dof_handler.n_dofs());
  }


  // @sect4{ConservationLaw::output_results}

  // This function now is rather straightforward. All the magic, including
  // transforming data from conservative variables to physical ones has been
  // abstracted and moved into the EulerEquations class so that it can be
  // replaced in case we want to solve some other hyperbolic conservation law.
  //
  // Note that the number of the output file is determined by keeping a
  // counter in the form of a static variable that is set to zero the first
  // time we come to this function and is incremented by one at the end of
  // each invocation.
  template <int dim>
  void ConservationLaw<dim>::output_results () const
  {
    typename EulerEquations<dim>::Postprocessor
    postprocessor (parameters.schlieren_plot);

    DataOut<dim> data_out;
    data_out.attach_dof_handler (dof_handler);

    data_out.add_data_vector (current_solution,
                              EulerEquations<dim>::component_names (),
                              DataOut<dim>::type_dof_data,
                              EulerEquations<dim>::component_interpretation ());

    data_out.add_data_vector (current_solution, postprocessor);

    data_out.build_patches ();

    static unsigned int output_file_number = 0;
    std::string filename = "solution-" +
                           Utilities::int_to_string (output_file_number, 3) +
                           ".vtk";
    std::ofstream output (filename.c_str());
    data_out.write_vtk (output);

    ++output_file_number;
  }




  // @sect4{ConservationLaw::run}

  // This function contains the top-level logic of this program:
  // initialization, the time loop, and the inner Newton iteration.
  //
  // At the beginning, we read the mesh file specified by the parameter file,
  // setup the DoFHandler and various vectors, and then interpolate the given
  // initial conditions on this mesh. We then perform a number of mesh
  // refinements, based on the initial conditions, to obtain a mesh that is
  // already well adapted to the starting solution. At the end of this
  // process, we output the initial solution.
  template <int dim>
  void ConservationLaw<dim>::run ()
  {
    {
      GridIn<dim> grid_in;
      grid_in.attach_triangulation(triangulation);

      std::ifstream input_file(parameters.mesh_filename.c_str());
      Assert (input_file, ExcFileNotOpen(parameters.mesh_filename.c_str()));

      grid_in.read_ucd(input_file);
    }

    dof_handler.clear();
    dof_handler.distribute_dofs (fe);

    // Size all of the fields.
    old_solution.reinit (dof_handler.n_dofs());
    current_solution.reinit (dof_handler.n_dofs());
    predictor.reinit (dof_handler.n_dofs());
    right_hand_side.reinit (dof_handler.n_dofs());

    setup_system();

    VectorTools::interpolate(dof_handler,
                             parameters.initial_conditions, old_solution);
    current_solution = old_solution;
    predictor = old_solution;

    if (parameters.do_refine == true)
      for (unsigned int i=0; i<parameters.shock_levels; ++i)
        {
          Vector<double> refinement_indicators (triangulation.n_active_cells());

          compute_refinement_indicators(refinement_indicators);
          refine_grid(refinement_indicators);

          setup_system();

          VectorTools::interpolate(dof_handler,
                                   parameters.initial_conditions, old_solution);
          current_solution = old_solution;
          predictor = old_solution;
        }

    output_results ();

    // We then enter into the main time stepping loop. At the top we simply
    // output some status information so one can keep track of where a
    // computation is, as well as the header for a table that indicates
    // progress of the nonlinear inner iteration:
    Vector<double> newton_update (dof_handler.n_dofs());

    double time = 0;
    double next_output = time + parameters.output_step;

    predictor = old_solution;
    while (time < parameters.final_time)
      {
        std::cout << "T=" << time << std::endl
                  << "   Number of active cells:       "
                  << triangulation.n_active_cells()
                  << std::endl
                  << "   Number of degrees of freedom: "
                  << dof_handler.n_dofs()
                  << std::endl
                  << std::endl;

        std::cout << "   NonLin Res     Lin Iter       Lin Res" << std::endl
                  << "   _____________________________________" << std::endl;

        // Then comes the inner Newton iteration to solve the nonlinear
        // problem in each time step. The way it works is to reset matrix and
        // right hand side to zero, then assemble the linear system. If the
        // norm of the right hand side is small enough, then we declare that
        // the Newton iteration has converged. Otherwise, we solve the linear
        // system, update the current solution with the Newton increment, and
        // output convergence information. At the end, we check that the
        // number of Newton iterations is not beyond a limit of 10 -- if it
        // is, it appears likely that iterations are diverging and further
        // iterations would do no good. If that happens, we throw an exception
        // that will be caught in <code>main()</code> with status information
        // being displayed before the program aborts.
        //
        // Note that the way we write the AssertThrow macro below is by and
        // large equivalent to writing something like <code>if (!(nonlin_iter
        // @<= 10)) throw ExcMessage ("No convergence in nonlinear
        // solver");</code>. The only significant difference is that
        // AssertThrow also makes sure that the exception being thrown carries
        // with it information about the location (file name and line number)
        // where it was generated. This is not overly critical here, because
        // there is only a single place where this sort of exception can
        // happen; however, it is generally a very useful tool when one wants
        // to find out where an error occurred.
        unsigned int nonlin_iter = 0;
        current_solution = predictor;
        while (true)
          {
            system_matrix = 0;

            right_hand_side = 0;
            assemble_system ();

            const double res_norm = right_hand_side.l2_norm();
            if (std::fabs(res_norm) < 1e-10)
              {
                std::printf("   %-16.3e (converged)\n\n", res_norm);
                break;
              }
            else
              {
                newton_update = 0;

                std::pair<unsigned int, double> convergence
                  = solve (newton_update);

                current_solution += newton_update;

                std::printf("   %-16.3e %04d        %-5.2e\n",
                            res_norm, convergence.first, convergence.second);
              }

            ++nonlin_iter;
            AssertThrow (nonlin_iter <= 10,
                         ExcMessage ("No convergence in nonlinear solver"));
          }

        // We only get to this point if the Newton iteration has converged, so
        // do various post convergence tasks here:
        //
        // First, we update the time and produce graphical output if so
        // desired. Then we update a predictor for the solution at the next
        // time step by approximating $\mathbf w^{n+1}\approx \mathbf w^n +
        // \delta t \frac{\partial \mathbf w}{\partial t} \approx \mathbf w^n
        // + \delta t \; \frac{\mathbf w^n-\mathbf w^{n-1}}{\delta t} = 2
        // \mathbf w^n - \mathbf w^{n-1}$ to try and make adaptivity work
        // better.  The idea is to try and refine ahead of a front, rather
        // than stepping into a coarse set of elements and smearing the
        // old_solution.  This simple time extrapolator does the job. With
        // this, we then refine the mesh if so desired by the user, and
        // finally continue on with the next time step:
        time += parameters.time_step;

        if (parameters.output_step < 0)
          output_results ();
        else if (time >= next_output)
          {
            output_results ();
            next_output += parameters.output_step;
          }

        predictor = current_solution;
        predictor.sadd (2.0, -1.0, old_solution);

        old_solution = current_solution;

        if (parameters.do_refine == true)
          {
            Vector<double> refinement_indicators (triangulation.n_active_cells());
            compute_refinement_indicators(refinement_indicators);

            refine_grid(refinement_indicators);
            setup_system();

            newton_update.reinit (dof_handler.n_dofs());
          }
      }
  }
}

// @sect3{main()}

// The following ``main'' function is similar to previous examples and need
// not to be commented on. Note that the program aborts if no input file name
// is given on the command line.
int main (int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      using namespace Step33;

      deallog.depth_console(0);
      if (argc != 2)
        {
          std::cout << "Usage:" << argv[0] << " input_file" << std::endl;
          std::exit(1);
        }

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, dealii::numbers::invalid_unsigned_int);

      ConservationLaw<2> cons (argv[1]);
      cons.run ();
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
    };

  return 0;
}
