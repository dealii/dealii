/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2012 - 2014 by the deal.II authors
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
 * Authors: Joerg Frohne, Texas A&M University and
 *                        University of Siegen, 2012, 2013
 *          Wolfgang Bangerth, Texas A&M University, 2012, 2013
 *          Timo Heister, Texas A&M University, 2013
 */

// @sect3{Include files}
// The set of include files is not much of a surprise any more at this time:
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/timer.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/solution_transfer.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/fe_field_function.h>

#include <fstream>
#include <iostream>

// This final include file provides the <code>mkdir</code> function
// that we will use to create a directory for output files, if necessary:
#include <sys/stat.h>

namespace Step42
{
  using namespace dealii;

  // @sect3{The <code>ConstitutiveLaw</code> class template}

  // This class provides an interface for a constitutive law, i.e., for the
  // relationship between strain $\varepsilon(\mathbf u)$ and stress
  // $\sigma$. In this example we are using an elastoplastic material behavior
  // with linear, isotropic hardening. Such materials are characterized by
  // Young's modulus $E$, Poisson's ratio $\nu$, the initial yield stress
  // $\sigma_0$ and the isotropic hardening parameter $\gamma$.  For $\gamma =
  // 0$ we obtain perfect elastoplastic behavior.
  //
  // As explained in the paper that describes this program, the first Newton
  // steps are solved with a completely elastic material model to avoid having
  // to deal with both nonlinearities (plasticity and contact) at once. To this
  // end, this class has a function <code>set_sigma_0()</code> that we use later
  // on to simply set $\sigma_0$ to a very large value -- essentially
  // guaranteeing that the actual stress will not exceed it, and thereby
  // producing an elastic material. When we are ready to use a plastic model, we
  // set $\sigma_0$ back to its proper value, using the same function.  As a
  // result of this approach, we need to leave <code>sigma_0</code> as the only
  // non-const member variable of this class.
  template <int dim>
  class ConstitutiveLaw
  {
  public:
    ConstitutiveLaw (const double E,
                     const double nu,
                     const double sigma_0,
                     const double gamma);

    void
    set_sigma_0 (double sigma_zero);

    bool
    get_stress_strain_tensor (const SymmetricTensor<2, dim> &strain_tensor,
                              SymmetricTensor<4, dim> &stress_strain_tensor) const;

    void
    get_linearized_stress_strain_tensors (const SymmetricTensor<2, dim> &strain_tensor,
                                          SymmetricTensor<4, dim> &stress_strain_tensor_linearized,
                                          SymmetricTensor<4, dim> &stress_strain_tensor) const;

  private:
    const double kappa;
    const double mu;
    double       sigma_0;
    const double gamma;

    const SymmetricTensor<4, dim> stress_strain_tensor_kappa;
    const SymmetricTensor<4, dim> stress_strain_tensor_mu;
  };

  // The constructor of the ConstitutiveLaw class sets the required material
  // parameter for our deformable body. Material parameters for elastic
  // isotropic media can be defined in a variety of ways, such as the pair $E,
  // \nu$ (elastic modulus and Poisson's number), using the Lame parameters
  // $\lambda,mu$ or several other commonly used conventions. Here, the
  // constructor takes a description of material parameters in the form of
  // $E,\nu$, but since this turns out to these are not the coefficients that
  // appear in the equations of the plastic projector, we immediately convert
  // them into the more suitable set $\kappa,\mu$ of bulk and shear moduli.  In
  // addition, the constructor takes $\sigma_0$ (the yield stress absent any
  // plastic strain) and $\gamma$ (the hardening parameter) as arguments. In
  // this constructor, we also compute the two principal components of the
  // stress-strain relation and its linearization.
  template <int dim>
  ConstitutiveLaw<dim>::ConstitutiveLaw (double E,
                                         double nu,
                                         double sigma_0,
                                         double gamma)
    :
    kappa (E / (3 * (1 - 2 * nu))),
    mu (E / (2 * (1 + nu))),
    sigma_0(sigma_0),
    gamma(gamma),
    stress_strain_tensor_kappa (kappa
                                * outer_product(unit_symmetric_tensor<dim>(),
                                                unit_symmetric_tensor<dim>())),
    stress_strain_tensor_mu (2 * mu
                             * (identity_tensor<dim>()
                                - outer_product(unit_symmetric_tensor<dim>(),
                                                unit_symmetric_tensor<dim>()) / 3.0))
  {}


  template <int dim>
  void
  ConstitutiveLaw<dim>::set_sigma_0 (double sigma_zero)
  {
    sigma_0 = sigma_zero;
  }


  // @sect4{ConstitutiveLaw::get_stress_strain_tensor}

  // This is the principal component of the constitutive law. It projects the
  // deviatoric part of the stresses in a quadrature point back to the yield
  // stress (i.e., the original yield stress $\sigma_0$ plus the term that
  // describes linear isotropic hardening).  We need this function to calculate
  // the nonlinear residual in PlasticityContactProblem::residual_nl_system. The
  // computations follow the formulas laid out in the introduction.
  //
  // The function returns whether the quadrature point is plastic to allow for
  // some statistics downstream on how many of the quadrature points are
  // plastic and how many are elastic.
  template <int dim>
  bool
  ConstitutiveLaw<dim>::
  get_stress_strain_tensor (const SymmetricTensor<2, dim> &strain_tensor,
                            SymmetricTensor<4, dim> &stress_strain_tensor) const
  {
    Assert (dim == 3, ExcNotImplemented());

    SymmetricTensor<2, dim> stress_tensor;
    stress_tensor = (stress_strain_tensor_kappa + stress_strain_tensor_mu)
                    * strain_tensor;

    const SymmetricTensor<2, dim> deviator_stress_tensor = deviator(stress_tensor);
    const double deviator_stress_tensor_norm = deviator_stress_tensor.norm();

    stress_strain_tensor = stress_strain_tensor_mu;
    if (deviator_stress_tensor_norm > sigma_0)
      {
        const double beta = sigma_0 / deviator_stress_tensor_norm;
        stress_strain_tensor *= (gamma + (1 - gamma) * beta);
      }

    stress_strain_tensor += stress_strain_tensor_kappa;

    return (deviator_stress_tensor_norm > sigma_0);
  }


  // @sect4{ConstitutiveLaw::get_linearized_stress_strain_tensors}

  // This function returns the linearized stress strain tensor, linearized
  // around the solution $u^{i-1}$ of the previous Newton step $i-1$.  The
  // parameter <code>strain_tensor</code> (commonly denoted
  // $\varepsilon(u^{i-1})$) must be passed as an argument, and serves as the
  // linearization point. The function returns the derivative of the nonlinear
  // constitutive law in the variable stress_strain_tensor, as well as the
  // stress-strain tensor of the linearized problem in
  // stress_strain_tensor_linearized.  See
  // PlasticityContactProblem::assemble_nl_system where this function is used.
  template <int dim>
  void
  ConstitutiveLaw<dim>::
  get_linearized_stress_strain_tensors (const SymmetricTensor<2, dim> &strain_tensor,
                                        SymmetricTensor<4, dim> &stress_strain_tensor_linearized,
                                        SymmetricTensor<4, dim> &stress_strain_tensor) const
  {
    Assert (dim == 3, ExcNotImplemented());

    SymmetricTensor<2, dim> stress_tensor;
    stress_tensor = (stress_strain_tensor_kappa + stress_strain_tensor_mu)
                    * strain_tensor;

    stress_strain_tensor = stress_strain_tensor_mu;
    stress_strain_tensor_linearized = stress_strain_tensor_mu;

    SymmetricTensor<2, dim> deviator_stress_tensor = deviator(stress_tensor);
    const double deviator_stress_tensor_norm = deviator_stress_tensor.norm();

    if (deviator_stress_tensor_norm > sigma_0)
      {
        const double beta = sigma_0 / deviator_stress_tensor_norm;
        stress_strain_tensor *= (gamma + (1 - gamma) * beta);
        stress_strain_tensor_linearized *= (gamma + (1 - gamma) * beta);
        deviator_stress_tensor /= deviator_stress_tensor_norm;
        stress_strain_tensor_linearized -= (1 - gamma) * beta * 2 * mu
                                           * outer_product(deviator_stress_tensor,
                                                           deviator_stress_tensor);
      }

    stress_strain_tensor += stress_strain_tensor_kappa;
    stress_strain_tensor_linearized += stress_strain_tensor_kappa;
  }

  // <h3>Equation data: boundary forces, boundary values, obstacles</h3>
  //
  // The following should be relatively standard. We need classes for
  // the boundary forcing term (which we here choose to be zero)
  // and boundary values on those part of the boundary that are not part
  // of the contact surface (also chosen to be zero here).
  namespace EquationData
  {
    template <int dim>
    class BoundaryForce : public Function<dim>
    {
    public:
      BoundaryForce ();

      virtual
      double value (const Point<dim> &p,
                    const unsigned int component = 0) const;

      virtual
      void vector_value (const Point<dim> &p,
                         Vector<double> &values) const;
    };

    template <int dim>
    BoundaryForce<dim>::BoundaryForce ()
      :
      Function<dim>(dim)
    {}


    template <int dim>
    double
    BoundaryForce<dim>::value (const Point<dim> &,
                               const unsigned int) const
    {
      return 0.;
    }

    template <int dim>
    void
    BoundaryForce<dim>::vector_value (const Point<dim> &p,
                                      Vector<double> &values) const
    {
      for (unsigned int c = 0; c < this->n_components; ++c)
        values(c) = BoundaryForce<dim>::value(p, c);
    }



    template <int dim>
    class BoundaryValues : public Function<dim>
    {
    public:
      BoundaryValues ();

      virtual double value (const Point<dim> &p,
                            const unsigned int component = 0) const;

      virtual
      void vector_value (const Point<dim> &p,
                         Vector<double> &values) const;
    };


    template <int dim>
    BoundaryValues<dim>::BoundaryValues ()
      :
      Function<dim>(dim)
    {}


    template <int dim>
    double
    BoundaryValues<dim>::value (const Point<dim> &,
                                const unsigned int) const
    {
      return 0.;
    }

    template <int dim>
    void
    BoundaryValues<dim>::vector_value (const Point<dim> &p,
                                       Vector<double> &values) const
    {
      for (unsigned int c = 0; c < this->n_components; ++c)
        values(c) = BoundaryValues<dim>::value(p, c);
    }


    // @sect4{The <code>SphereObstacle</code> class}

    // The following class is the first of two obstacles that can be
    // selected from the input file. It describes a sphere centered
    // at position $x=y=0.5, z=z_{\text{surface}}+0.59$ and radius $r=0.6$,
    // where $z_{\text{surface}}$ is the vertical position of the (flat)
    // surface of the deformable body. The function's <code>value</code>
    // returns the location of the obstacle for a given $x,y$ value if the
    // point actually lies below the sphere, or a large positive value that
    // can't possibly interfere with the deformation if it lies outside
    // the "shadow" of the sphere.
    template <int dim>
    class SphereObstacle : public Function<dim>
    {
    public:
      SphereObstacle (const double z_surface);

      virtual
      double value (const Point<dim> &p,
                    const unsigned int component = 0) const;

      virtual
      void vector_value (const Point<dim> &p,
                         Vector<double> &values) const;

    private:
      const double z_surface;
    };


    template <int dim>
    SphereObstacle<dim>::SphereObstacle (const double z_surface)
      :
      Function<dim>(dim),
      z_surface(z_surface)
    {}


    template <int dim>
    double
    SphereObstacle<dim>::value (
      const Point<dim> &p, const unsigned int component) const
    {
      if (component == 0)
        return p(0);
      else if (component == 1)
        return p(1);
      else if (component == 2)
        {
          if ((p(0) - 0.5) * (p(0) - 0.5) + (p(1) - 0.5) * (p(1) - 0.5)
              < 0.36)
            return (-std::sqrt(
                      0.36 - (p(0) - 0.5) * (p(0) - 0.5)
                      - (p(1) - 0.5) * (p(1) - 0.5)) + z_surface + 0.59);
          else
            return 1000;
        }

      Assert(false, ExcNotImplemented());
      return 1e9; // an unreasonable value; ignored in debug mode because of the preceding Assert
    }


    template <int dim>
    void
    SphereObstacle<dim>::vector_value (const Point<dim> &p,
                                       Vector<double> &values) const
    {
      for (unsigned int c = 0; c < this->n_components; ++c)
        values(c) = SphereObstacle<dim>::value(p, c);
    }

    // @sect4{The <code>BitmapFile</code> and <code>ChineseObstacle</code> classes}

    // The following two classes describe the obstacle outlined in the introduction,
    // i.e., the Chinese character. The first of the two, <code>BitmapFile</code>
    // is responsible for reading in data from a picture file
    // stored in pbm ascii format. This data will be bilinearly interpolated and
    // thereby provides a function that describes the obstacle. (The code below
    // shows how one can construct a function by interpolating between given
    // data points. One could use the Functions::InterpolatedUniformGridData,
    // introduced after this tutorial program was written, which does exactly
    // what we want here, but it is instructive to see how to do it by hand.)
    //
    // The data which we read from the file will be stored in a double std::vector
    // named obstacle_data.  This vector composes the base to calculate a
    // piecewise bilinear function as a polynomial interpolation. The data we will
    // read from a file consists of zeros (white) and ones (black).
    //
    // The <code>hx,hy</code> variables denote the spacing between pixels in $x$
    // and $y$ directions. <code>nx,ny</code> are the numbers of pixels in each of
    // these directions.  <code>get_value()</code> returns the value of the image
    // at a given location, interpolated from the adjacent pixel values.
    template <int dim>
    class BitmapFile
    {
    public:
      BitmapFile(const std::string &name);

      double
      get_value(const double x, const double y) const;

    private:
      std::vector<double> obstacle_data;
      double hx, hy;
      int nx, ny;

      double
      get_pixel_value(const int i, const int j) const;
    };

    // The constructor of this class reads in the data that describes
    // the obstacle from the given file name.
    template <int dim>
    BitmapFile<dim>::BitmapFile(const std::string &name)
      :
      obstacle_data(0),
      hx(0),
      hy(0),
      nx(0),
      ny(0)
    {
      std::ifstream f(name.c_str());
      AssertThrow (f, ExcMessage (std::string("Can't read from file <") +
                                  name + ">!"));

      std::string temp;
      f >> temp >> nx >> ny;

      AssertThrow(nx > 0 && ny > 0, ExcMessage("Invalid file format."));

      for (int k = 0; k < nx * ny; k++)
        {
          double val;
          f >> val;
          obstacle_data.push_back(val);
        }

      hx = 1.0 / (nx - 1);
      hy = 1.0 / (ny - 1);

      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        std::cout << "Read obstacle from file <" << name << ">" << std::endl
                  << "Resolution of the scanned obstacle picture: " << nx
                  << " x " << ny << std::endl;
    }

    // The following two functions return the value of a given pixel with
    // coordinates $i,j$, which we identify with the values of a function
    // defined at positions <code>i*hx, j*hy</code>, and at arbitrary
    // coordinates $x,y$ where we do a bilinear interpolation between
    // point values returned by the first of the two functions. In the
    // second function, for each $x,y$, we first compute the (integer)
    // location of the nearest pixel coordinate to the bottom left of
    // $x,y$, and then compute the coordinates $\xi,\eta$ within this
    // pixel. We truncate both kinds of variables from both below
    // and above to avoid problems when evaluating the function outside
    // of its defined range as may happen due to roundoff errors.
    template <int dim>
    double
    BitmapFile<dim>::get_pixel_value(const int i,
                                     const int j) const
    {
      assert(i >= 0 && i < nx);
      assert(j >= 0 && j < ny);
      return obstacle_data[nx * (ny - 1 - j) + i];
    }

    template <int dim>
    double
    BitmapFile<dim>::get_value(const double x,
                               const double y) const
    {
      const int ix = std::min(std::max((int) (x / hx), 0), nx - 2);
      const int iy = std::min(std::max((int) (y / hy), 0), ny - 2);

      const double xi  = std::min(std::max((x-ix*hx)/hx, 1.), 0.);
      const double eta = std::min(std::max((y-iy*hy)/hy, 1.), 0.);

      return ((1-xi)*(1-eta)*get_pixel_value(ix,iy)
              +
              xi*(1-eta)*get_pixel_value(ix+1,iy)
              +
              (1-xi)*eta*get_pixel_value(ix,iy+1)
              +
              xi*eta*get_pixel_value(ix+1,iy+1));
    }

    // Finally, this is the class that actually uses the class above. It
    // has a BitmapFile object as a member that describes the height of the
    // obstacle. As mentioned above, the BitmapFile class will provide us
    // with a mask, i.e., values that are either zero or one (and, if you
    // ask for locations between pixels, values that are interpolated between
    // zero and one). This class translates this to heights that are either
    // 0.001 below the surface of the deformable body (if the BitmapFile
    // class reports a one at this location) or 0.999 above the obstacle (if
    // the BitmapFile class reports a zero). The following function should then
    // be self-explanatory.
    template <int dim>
    class ChineseObstacle : public Function<dim>
    {
    public:
      ChineseObstacle(const std::string &filename,
                      const double z_surface);

      virtual
      double value (const Point<dim> &p,
                    const unsigned int component = 0) const;

      virtual
      void vector_value (const Point<dim> &p,
                         Vector<double> &values) const;

    private:
      const BitmapFile<dim> input_obstacle;
      double z_surface;
    };


    template <int dim>
    ChineseObstacle<dim>::ChineseObstacle(const std::string &filename,
                                          const double z_surface)
      :
      Function<dim>(dim),
      input_obstacle(filename),
      z_surface(z_surface)
    {}


    template <int dim>
    double
    ChineseObstacle<dim>::value (const Point<dim> &p,
                                 const unsigned int component) const
    {
      if (component == 0)
        return p(0);
      if (component == 1)
        return p(1);
      else if (component==2)
        {
          if (p(0) >= 0.0 && p(0) <= 1.0 && p(1) >= 0.0 && p(1) <= 1.0)
            return z_surface + 0.999 - input_obstacle.get_value(p(0), p(1));
        }

      Assert (false, ExcNotImplemented());
      return 1e9; // an unreasonable value; ignored in debug mode because of the preceding Assert
    }

    template <int dim>
    void
    ChineseObstacle<dim>::vector_value (const Point<dim> &p,
                                        Vector<double> &values) const
    {
      for (unsigned int c = 0; c < this->n_components; ++c)
        values(c) = ChineseObstacle<dim>::value(p, c);
    }
  }

  // @sect3{The <code>PlasticityContactProblem</code> class template}

  // This is the main class of this program and supplies all functions
  // and variables needed to describe
  // the nonlinear contact problem. It is
  // close to step-41 but with some additional
  // features like handling hanging nodes,
  // a Newton method, using Trilinos and p4est
  // for parallel distributed computing.
  // To deal with hanging nodes makes
  // life a bit more complicated since
  // we need another ConstraintMatrix now.
  // We create a Newton method for the
  // active set method for the contact
  // situation and to handle the nonlinear
  // operator for the constitutive law.
  //
  // The general layout of this class is very much like for most other tutorial programs.
  // To make our life a bit easier, this class reads a set of input parameters from an input file. These
  // parameters, using the ParameterHandler class, are declared in the <code>declare_parameters</code>
  // function (which is static so that it can be called before we even create an object of the current
  // type), and a ParameterHandler object that has been used to read an input file will then be passed
  // to the constructor of this class.
  //
  // The remaining member functions are by and large as we have seen in several of the other tutorial
  // programs, though with additions for the current nonlinear system. We will comment on their purpose
  // as we get to them further below.
  template <int dim>
  class PlasticityContactProblem
  {
  public:
    PlasticityContactProblem (const ParameterHandler &prm);

    void run ();

    static void declare_parameters (ParameterHandler &prm);

  private:
    void make_grid ();
    void setup_system ();
    void compute_dirichlet_constraints ();
    void update_solution_and_constraints ();
    void assemble_mass_matrix_diagonal (TrilinosWrappers::SparseMatrix &mass_matrix);
    void assemble_newton_system (const TrilinosWrappers::MPI::Vector &linearization_point);
    void compute_nonlinear_residual (const TrilinosWrappers::MPI::Vector &linearization_point);
    void solve_newton_system ();
    void solve_newton ();
    void refine_grid ();
    void move_mesh (const TrilinosWrappers::MPI::Vector &displacement) const;
    void output_results (const std::string &filename_base);
    void output_contact_force () const;

    // As far as member variables are concerned, we start with ones that we use to
    // indicate the MPI universe this program runs on, a stream we use to let
    // exactly one processor produce output to the console (see step-17) and
    // a variable that is used to time the various sections of the program:
    MPI_Comm           mpi_communicator;
    ConditionalOStream pcout;
    TimerOutput        computing_timer;

    // The next group describes the mesh and the finite element space.
    // In particular, for this parallel program, the finite element
    // space has associated with it variables that indicate which degrees
    // of freedom live on the current processor (the index sets, see
    // also step-40 and the @ref distributed documentation module) as
    // well as a variety of constraints: those imposed by hanging nodes,
    // by Dirichlet boundary conditions, and by the active set of
    // contact nodes. Of the three ConstraintMatrix variables defined
    // here, the first only contains hanging node constraints, the
    // second also those associated with Dirichlet boundary conditions,
    // and the third these plus the contact constraints.
    //
    // The variable <code>active_set</code> consists of those degrees
    // of freedom constrained by the contact, and we use
    // <code>fraction_of_plastic_q_points_per_cell</code> to keep
    // track of the fraction of quadrature points on each cell where
    // the stress equals the yield stress. The latter is only used to
    // create graphical output showing the plastic zone, but not for
    // any further computation; the variable is a member variable of
    // this class since the information is computed as a by-product
    // of computing the residual, but is used only much later. (Note
    // that the vector is a vector of length equal to the number of
    // active cells on the <i>local mesh</i>; it is never used to
    // exchange information between processors and can therefore be
    // a regular deal.II vector.)
    const unsigned int                        n_initial_global_refinements;
    parallel::distributed::Triangulation<dim> triangulation;

    const unsigned int fe_degree;
    FESystem<dim>      fe;
    DoFHandler<dim>    dof_handler;

    IndexSet           locally_owned_dofs;
    IndexSet           locally_relevant_dofs;

    ConstraintMatrix   constraints_hanging_nodes;
    ConstraintMatrix   constraints_dirichlet_and_hanging_nodes;
    ConstraintMatrix   all_constraints;

    IndexSet           active_set;
    Vector<float>      fraction_of_plastic_q_points_per_cell;


    // The next block of variables corresponds to the solution
    // and the linear systems we need to form. In particular, this
    // includes the Newton matrix and right hand side; the vector
    // that corresponds to the residual (i.e., the Newton right hand
    // side) but from which we have not eliminated the various
    // constraints and that is used to determine which degrees of
    // freedom need to be constrained in the next iteration; and
    // a vector that corresponds to the diagonal of the $B$ matrix
    // briefly mentioned in the introduction and discussed in the
    // accompanying paper.
    TrilinosWrappers::SparseMatrix    newton_matrix;

    TrilinosWrappers::MPI::Vector     solution;
    TrilinosWrappers::MPI::Vector     newton_rhs;
    TrilinosWrappers::MPI::Vector     newton_rhs_uncondensed;
    TrilinosWrappers::MPI::Vector     diag_mass_matrix_vector;

    // The next block contains the variables that describe the material
    // response:
    const double         e_modulus, nu, gamma, sigma_0;
    ConstitutiveLaw<dim> constitutive_law;

    // And then there is an assortment of other variables that are used
    // to identify the mesh we are asked to build as selected by the
    // parameter file, the obstacle that is being pushed into the
    // deformable body, the mesh refinement strategy, whether to transfer
    // the solution from one mesh to the next, and how many mesh
    // refinement cycles to perform. As possible, we mark these kinds
    // of variables as <code>const</code> to help the reader identify
    // which ones may or may not be modified later on (the output directory
    // being an exception -- it is never modified outside the constructor
    // but it is awkward to initialize in the member-initializer-list
    // following the colon in the constructor since there we have only
    // one shot at setting it; the same is true for the mesh refinement
    // criterion):
    const std::string                                  base_mesh;
    const std_cxx11::shared_ptr<const Function<dim> >  obstacle;

    struct RefinementStrategy
    {
      enum value
      {
        refine_global,
        refine_percentage,
        refine_fix_dofs
      };
    };
    typename RefinementStrategy::value                 refinement_strategy;

    const bool                                         transfer_solution;
    std::string                                        output_dir;
    const unsigned int                                 n_refinement_cycles;
    unsigned int                                       current_refinement_cycle;
  };


  // @sect3{Implementation of the <code>PlasticityContactProblem</code> class}

  // @sect4{PlasticityContactProblem::declare_parameters}

  // Let us start with the declaration of run-time parameters that can be
  // selected in the input file. These values will be read back in the
  // constructor of this class to initialize the member variables of this
  // class:
  template <int dim>
  void
  PlasticityContactProblem<dim>::declare_parameters (ParameterHandler &prm)
  {
    prm.declare_entry("polynomial degree", "1",
                      Patterns::Integer(),
                      "Polynomial degree of the FE_Q finite element space, typically 1 or 2.");
    prm.declare_entry("number of initial refinements", "2",
                      Patterns::Integer(),
                      "Number of initial global mesh refinement steps before "
                      "the first computation.");
    prm.declare_entry("refinement strategy", "percentage",
                      Patterns::Selection("global|percentage"),
                      "Mesh refinement strategy:\n"
                      " global: one global refinement\n"
                      " percentage: a fixed percentage of cells gets refined using the Kelly estimator.");
    prm.declare_entry("number of cycles", "5",
                      Patterns::Integer(),
                      "Number of adaptive mesh refinement cycles to run.");
    prm.declare_entry("obstacle", "sphere",
                      Patterns::Selection("sphere|read from file"),
                      "The name of the obstacle to use. This may either be 'sphere' if we should "
                      "use a spherical obstacle, or 'read from file' in which case the obstacle "
                      "will be read from a file named 'obstacle.pbm' that is supposed to be in "
                      "ASCII PBM format.");
    prm.declare_entry("output directory", "",
                      Patterns::Anything(),
                      "Directory for output files (graphical output and benchmark "
                      "statistics). If empty, use the current directory.");
    prm.declare_entry("transfer solution", "false",
                      Patterns::Bool(),
                      "Whether the solution should be used as a starting guess "
                      "for the next finer mesh. If false, then the iteration starts at "
                      "zero on every mesh.");
    prm.declare_entry("base mesh", "box",
                      Patterns::Selection("box|half sphere"),
                      "Select the shape of the domain: 'box' or 'half sphere'");
  }


  // @sect4{The <code>PlasticityContactProblem</code> constructor}

  // Given the declarations of member variables as well as the
  // declarations of run-time parameters that are read from the input
  // file, there is nothing surprising in this constructor. In the body
  // we initialize the mesh refinement strategy and the output directory,
  // creating such a directory if necessary.
  template <int dim>
  PlasticityContactProblem<dim>::
  PlasticityContactProblem (const ParameterHandler &prm)
    :
    mpi_communicator(MPI_COMM_WORLD),
    pcout(std::cout,
          (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
    computing_timer(MPI_COMM_WORLD, pcout, TimerOutput::never,
                    TimerOutput::wall_times),

    n_initial_global_refinements (prm.get_integer("number of initial refinements")),
    triangulation(mpi_communicator),
    fe_degree (prm.get_integer("polynomial degree")),
    fe(FE_Q<dim>(QGaussLobatto<1>(fe_degree+1)), dim),
    dof_handler(triangulation),

    e_modulus (200000),
    nu (0.3),
    gamma (0.01),
    sigma_0(400.0),
    constitutive_law (e_modulus,
                      nu,
                      sigma_0,
                      gamma),

    base_mesh (prm.get("base mesh")),
    obstacle (prm.get("obstacle") == "read from file"
              ?
              static_cast<const Function<dim>*>
              (new EquationData::ChineseObstacle<dim>("obstacle.pbm", (base_mesh == "box" ? 1.0 : 0.5)))
              :
              static_cast<const Function<dim>*>
              (new EquationData::SphereObstacle<dim>(base_mesh == "box" ? 1.0 : 0.5))),

    transfer_solution (prm.get_bool("transfer solution")),
    n_refinement_cycles (prm.get_integer("number of cycles"))
  {
    std::string strat = prm.get("refinement strategy");
    if (strat == "global")
      refinement_strategy = RefinementStrategy::refine_global;
    else if (strat == "percentage")
      refinement_strategy = RefinementStrategy::refine_percentage;
    else
      AssertThrow (false, ExcNotImplemented());

    output_dir = prm.get("output directory");
    if (output_dir != "" && *(output_dir.rbegin()) != '/')
      output_dir += "/";
    mkdir(output_dir.c_str(), 0777);

    pcout << "    Using output directory '" << output_dir << "'" << std::endl;
    pcout << "    FE degree " << fe_degree << std::endl;
    pcout << "    transfer solution "
          << (transfer_solution ? "true" : "false") << std::endl;
  }



  // @sect4{PlasticityContactProblem::make_grid}

  // The next block deals with constructing the starting mesh.
  // We will use the following helper function and the first
  // block of the <code>make_grid()</code> to construct a
  // mesh that corresponds to a half sphere. deal.II has a function
  // that creates such a mesh, but it is in the wrong location
  // and facing the wrong direction, so we need to shift and rotate
  // it a bit before using it.
  //
  // For later reference, as described in the documentation of
  // GridGenerator::half_hyper_ball(), the flat surface of the halfsphere
  // has boundary indicator zero, while the remainder has boundary
  // indicator one.
  Point<3>
  rotate_half_sphere (const Point<3> &in)
  {
    return Point<3>(in(2), in(1), -in(0));
  }

  template <int dim>
  void
  PlasticityContactProblem<dim>::make_grid ()
  {
    if (base_mesh == "half sphere")
      {
        const Point<dim> center(0, 0, 0);
        const double radius = 0.8;
        GridGenerator::half_hyper_ball(triangulation, center, radius);

        GridTools::transform(&rotate_half_sphere, triangulation);
        GridTools::shift(Point<dim>(0.5, 0.5, 0.5), triangulation);

        static HyperBallBoundary<dim> boundary_description(Point<dim>(0.5, 0.5, 0.5), radius);
        triangulation.set_boundary(0, boundary_description);
      }
    // Alternatively, create a hypercube mesh. After creating it,
    // assign boundary indicators as follows:
    // @code
    //    _______
    //   /  1    /|
    //  /______ / |
    // |       | 8|
    // |   8   | /
    // |_______|/
    //     6
    // @endcode
    // In other words, the boundary indicators of the sides of the cube are 8.
    // The boundary indicator of the bottom is 6 and the top has indicator 1.
    // We will make use of these indicators later when evaluating which
    // boundary will carry Dirichlet boundary conditions or will be
    // subject to potential contact.
    else
      {
        const Point<dim> p1(0, 0, 0);
        const Point<dim> p2(1.0, 1.0, 1.0);

        GridGenerator::hyper_rectangle(triangulation, p1, p2);

        /* boundary_indicators:

         */
        Triangulation<3>::active_cell_iterator
        cell = triangulation.begin_active(),
        endc = triangulation.end();
        for (; cell != endc; ++cell)
          for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
               ++face)
            {
              if (cell->face(face)->center()[2] == p2(2))
                cell->face(face)->set_boundary_indicator(1);
              if (cell->face(face)->center()[0] == p1(0)
                  || cell->face(face)->center()[0] == p2(0)
                  || cell->face(face)->center()[1] == p1(1)
                  || cell->face(face)->center()[1] == p2(1))
                cell->face(face)->set_boundary_indicator(8);
              if (cell->face(face)->center()[2] == p1(2))
                cell->face(face)->set_boundary_indicator(6);
            }
      }

    triangulation.refine_global(n_initial_global_refinements);
  }



  // @sect4{PlasticityContactProblem::setup_system}

  // The next piece in the puzzle is to set up the DoFHandler, resize
  // vectors and take care of various other status variables such as
  // index sets and constraint matrices.
  //
  // In the following, each group of operations is put into a brace-enclosed
  // block that is being timed by the variable declared at the top of the
  // block (the constructor of the TimerOutput::Scope variable starts the
  // timed section, the destructor that is called at the end of the block
  // stops it again).
  template <int dim>
  void
  PlasticityContactProblem<dim>::setup_system ()
  {
    /* setup dofs and get index sets for locally owned and relevant dofs */
    {
      TimerOutput::Scope t(computing_timer, "Setup: distribute DoFs");
      dof_handler.distribute_dofs(fe);

      locally_owned_dofs = dof_handler.locally_owned_dofs();
      locally_relevant_dofs.clear();
      DoFTools::extract_locally_relevant_dofs(dof_handler,
                                              locally_relevant_dofs);
    }

    /* setup hanging nodes and Dirichlet constraints */
    {
      TimerOutput::Scope t(computing_timer, "Setup: constraints");
      constraints_hanging_nodes.reinit(locally_relevant_dofs);
      DoFTools::make_hanging_node_constraints(dof_handler,
                                              constraints_hanging_nodes);
      constraints_hanging_nodes.close();

      pcout << "   Number of active cells: "
            << triangulation.n_global_active_cells() << std::endl
            << "   Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;

      compute_dirichlet_constraints();
    }

    /* initialization of vectors and the active set */
    {
      TimerOutput::Scope t(computing_timer, "Setup: vectors");
      solution.reinit(locally_relevant_dofs, mpi_communicator);
      newton_rhs.reinit(locally_owned_dofs, mpi_communicator);
      newton_rhs_uncondensed.reinit(locally_owned_dofs, mpi_communicator);
      diag_mass_matrix_vector.reinit(locally_owned_dofs, mpi_communicator);
      fraction_of_plastic_q_points_per_cell.reinit(triangulation.n_active_cells());

      active_set.clear();
      active_set.set_size(dof_handler.n_dofs());
    }

    // Finally, we set up sparsity patterns and matrices.
    // We temporarily (ab)use the system matrix to also build the (diagonal)
    // matrix that we use in eliminating degrees of freedom that are in contact
    // with the obstacle, but we then immediately set the Newton matrix back
    // to zero.
    {
      TimerOutput::Scope t(computing_timer, "Setup: matrix");
      TrilinosWrappers::SparsityPattern sp(locally_owned_dofs,
                                           mpi_communicator);

      DoFTools::make_sparsity_pattern(dof_handler, sp,
                                      constraints_dirichlet_and_hanging_nodes, false,
                                      Utilities::MPI::this_mpi_process(mpi_communicator));
      sp.compress();
      newton_matrix.reinit(sp);


      TrilinosWrappers::SparseMatrix &mass_matrix = newton_matrix;

      assemble_mass_matrix_diagonal(mass_matrix);

      const unsigned int start = (newton_rhs.local_range().first),
                         end = (newton_rhs.local_range().second);
      for (unsigned int j = start; j < end; j++)
        diag_mass_matrix_vector(j) = mass_matrix.diag_element(j);
      diag_mass_matrix_vector.compress(VectorOperation::insert);

      mass_matrix = 0;
    }
  }


  // @sect4{PlasticityContactProblem::compute_dirichlet_constraints}

  // This function, broken out of the preceding one, computes the constraints
  // associated with Dirichlet-type boundary conditions and puts them into the
  // <code>constraints_dirichlet_and_hanging_nodes</code> variable by merging
  // with the constraints that come from hanging nodes.
  //
  // As laid out in the introduction, we need to distinguish between two
  // cases:
  // - If the domain is a box, we set the displacement to zero at the bottom,
  //   and allow vertical movement in z-direction along the sides. As
  //   shown in the <code>make_grid()</code> function, the former corresponds
  //   to boundary indicator 6, the latter to 8.
  // - If the domain is a half sphere, then we impose zero displacement along
  //   the curved part of the boundary, associated with boundary indicator zero.
  template <int dim>
  void
  PlasticityContactProblem<dim>::compute_dirichlet_constraints ()
  {
    constraints_dirichlet_and_hanging_nodes.reinit(locally_relevant_dofs);
    constraints_dirichlet_and_hanging_nodes.merge(constraints_hanging_nodes);

    if (base_mesh == "box")
      {
        // interpolate all components of the solution
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 6,
                                                 EquationData::BoundaryValues<dim>(),
                                                 constraints_dirichlet_and_hanging_nodes,
                                                 ComponentMask());

        // interpolate x- and y-components of the
        // solution (this is a bit mask, so apply
        // operator| )
        const FEValuesExtractors::Scalar x_displacement(0);
        const FEValuesExtractors::Scalar y_displacement(1);
        VectorTools::interpolate_boundary_values(dof_handler,
                                                 8,
                                                 EquationData::BoundaryValues<dim>(),
                                                 constraints_dirichlet_and_hanging_nodes,
                                                 (fe.component_mask(x_displacement) | fe.component_mask(y_displacement)));
      }
    else
      VectorTools::interpolate_boundary_values(dof_handler,
                                               0,
                                               EquationData::BoundaryValues<dim>(),
                                               constraints_dirichlet_and_hanging_nodes,
                                               ComponentMask());

    constraints_dirichlet_and_hanging_nodes.close();
  }



  // @sect4{PlasticityContactProblem::assemble_mass_matrix_diagonal}

  // The next helper function computes the (diagonal) mass matrix that
  // is used to determine the active set of the active set method we use in
  // the contact algorithm. This matrix is of mass matrix type, but unlike
  // the standard mass matrix, we can make it diagonal (even in the case of
  // higher order elements) by using a quadrature formula that has its
  // quadrature points at exactly the same locations as the interpolation points
  // for the finite element are located. We achieve this by using a
  // QGaussLobatto quadrature formula here, along with initializing the finite
  // element with a set of interpolation points derived from the same quadrature
  // formula. The remainder of the function is relatively straightfoward: we
  // put the resulting matrix into the given argument; because we know the
  // matrix is diagonal, it is sufficient to have a loop over only $i$ not
  // not over $j$. Strictly speaking, we could even avoid multiplying the
  // shape function's values at quadrature point <code>q_point</code> by itself
  // because we know the shape value to be a vector with exactly one one which
  // when dotted with itself yields one. Since this function is not time
  // critical we add this term for clarity.
  template <int dim>
  void
  PlasticityContactProblem<dim>::
  assemble_mass_matrix_diagonal (TrilinosWrappers::SparseMatrix &mass_matrix)
  {
    QGaussLobatto<dim-1> face_quadrature_formula(fe.degree + 1);

    FEFaceValues<dim> fe_values_face(fe, face_quadrature_formula,
                                     update_values | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    const FEValuesExtractors::Vector displacement(0);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
        for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell;
             ++face)
          if (cell->face(face)->at_boundary()
              &&
              cell->face(face)->boundary_indicator() == 1)
            {
              fe_values_face.reinit(cell, face);
              cell_matrix = 0;

              for (unsigned int q_point = 0; q_point<n_face_q_points; ++q_point)
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  cell_matrix(i, i) += (fe_values_face[displacement].value(i, q_point) *
                                        fe_values_face[displacement].value(i, q_point) *
                                        fe_values_face.JxW(q_point));

              cell->get_dof_indices(local_dof_indices);

              for (unsigned int i = 0; i < dofs_per_cell; i++)
                mass_matrix.add(local_dof_indices[i],
                                local_dof_indices[i],
                                cell_matrix(i, i));
            }
    mass_matrix.compress(VectorOperation::add);
  }


  // @sect4{PlasticityContactProblem::update_solution_and_constraints}

  // The following function is the first function we call in each Newton
  // iteration in the <code>solve_newton()</code> function. What it does is
  // to project the solution onto the feasible set and update the active set
  // for the degrees of freedom that touch or penetrate the obstacle.
  //
  // In order to function, we first need to do some bookkeeping: We need
  // to write into the solution vector (which we can only do with fully
  // distributed vectors without ghost elements) and we need to read
  // the Lagrange multiplier and the elements of the diagonal mass matrix
  // from their respective vectors (which we can only do with vectors that
  // do have ghost elements), so we create the respective vectors. We then
  // also initialize the constraints object that will contain constraints
  // from contact and all other sources, as well as an object that contains
  // an index set of all locally owned degrees of freedom that are part of
  // the contact:
  template <int dim>
  void
  PlasticityContactProblem<dim>::update_solution_and_constraints ()
  {
    std::vector<bool> dof_touched(dof_handler.n_dofs(), false);

    TrilinosWrappers::MPI::Vector distributed_solution(locally_owned_dofs, mpi_communicator);
    distributed_solution = solution;

    TrilinosWrappers::MPI::Vector lambda(locally_relevant_dofs, mpi_communicator);
    lambda = newton_rhs_uncondensed;

    TrilinosWrappers::MPI::Vector diag_mass_matrix_vector_relevant(locally_relevant_dofs, mpi_communicator);
    diag_mass_matrix_vector_relevant = diag_mass_matrix_vector;


    all_constraints.reinit(locally_relevant_dofs);
    active_set.clear();

    // The second part is a loop over all cells in which we look at each
    // point where a degree of freedom is defined whether the active set
    // condition is true and we need to add this degree of freedom to
    // the active set of contact nodes. As we always do, if we want to
    // evaluate functions at individual points, we do this with an
    // FEValues object (or, here, an FEFaceValues object since we need to
    // check contact at the surface) with an appropriately chosen quadrature
    // object. We create this face quadrature object by choosing the
    // "support points" of the shape functions defined on the faces
    // of cells (for more on support points, see this
    // @ref GlossSupport "glossary entry"). As a consequence, we have as
    // many quadrature points as there are shape functions per face and
    // looping over quadrature points is equivalent to looping over shape
    // functions defined on a face. With this, the code looks as follows:
    Quadrature<dim-1> face_quadrature(fe.get_unit_face_support_points());
    FEFaceValues<dim> fe_values_face(fe, face_quadrature,
                                     update_quadrature_points);

    const unsigned int dofs_per_face = fe.dofs_per_face;
    const unsigned int n_face_q_points = face_quadrature.size();

    std::vector<types::global_dof_index> dof_indices(dofs_per_face);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    for (; cell != endc; ++cell)
      if (!cell->is_artificial())
        for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
          if (cell->face(face)->at_boundary()
              &&
              cell->face(face)->boundary_indicator() == 1)
            {
              fe_values_face.reinit(cell, face);
              cell->face(face)->get_dof_indices(dof_indices);

              for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
                {
                  // At each quadrature point (i.e., at each support point of a degree
                  // of freedom located on the contact boundary), we then ask whether
                  // it is part of the z-displacement degrees of freedom and if we
                  // haven't encountered this degree of freedom yet (which can happen
                  // for those on the edges between faces), we need to evaluate the gap
                  // between the deformed object and the obstacle. If the active set
                  // condition is true, then we add a constraint to the ConstraintMatrix
                  // object that the next Newton update needs to satisfy, set the solution
                  // vector's corresponding element to the correct value, and add the
                  // index to the IndexSet object that stores which degree of freedom is
                  // part of the contact:
                  const unsigned int
                  component = fe.face_system_to_component_index(q_point).first;

                  const unsigned int index_z = dof_indices[q_point];

                  if ((component == 2) && (dof_touched[index_z] == false))
                    {
                      dof_touched[index_z] = true;

                      const Point<dim> this_support_point = fe_values_face.quadrature_point(q_point);

                      const double obstacle_value = obstacle->value(this_support_point, 2);
                      const double solution_here = solution(index_z);
                      const double undeformed_gap = obstacle_value - this_support_point(2);

                      const double c = 100.0 * e_modulus;
                      if ((lambda(index_z) / diag_mass_matrix_vector_relevant(index_z)
                           +
                           c * (solution_here - undeformed_gap)
                           > 0)
                          &&
                          !constraints_hanging_nodes.is_constrained(index_z))
                        {
                          all_constraints.add_line(index_z);
                          all_constraints.set_inhomogeneity(index_z, undeformed_gap);
                          distributed_solution(index_z) = undeformed_gap;

                          active_set.add_index(index_z);
                        }
                    }
                }
            }

    // At the end of this function, we exchange data between processors updating
    // those ghost elements in the <code>solution</code> variable that have been
    // written by other processors. We then merge the Dirichlet constraints and
    // those from hanging nodes into the ConstraintMatrix object that already
    // contains the active set. We finish the function by outputting the total
    // number of actively constrained degrees of freedom for which we sum over
    // the number of actively constrained degrees of freedom owned by each
    // of the processors. This number of locally owned constrained degrees of
    // freedom is of course the number of elements of the intersection of the
    // active set and the set of locally owned degrees of freedom, which
    // we can get by using <code>operator&</code> on two IndexSets:
    distributed_solution.compress(VectorOperation::insert);
    solution = distributed_solution;

    all_constraints.close();
    all_constraints.merge(constraints_dirichlet_and_hanging_nodes);

    pcout << "         Size of active set: "
          << Utilities::MPI::sum((active_set & locally_owned_dofs).n_elements(),
                                 mpi_communicator)
          << std::endl;
  }


  // @sect4{PlasticityContactProblem::assemble_newton_system}

  // Given the complexity of the problem, it may come as a bit of a surprise
  // that assembling the linear system we have to solve in each Newton iteration
  // is actually fairly straightforward. The following function builds the Newton
  // right hand side and Newton matrix. It looks fairly innocent because the
  // heavy lifting happens in the call to
  // <code>ConstitutiveLaw::get_linearized_stress_strain_tensors()</code> and in
  // particular in ConstraintMatrix::distribute_local_to_global(), using the
  // constraints we have previously computed.
  template <int dim>
  void
  PlasticityContactProblem<dim>::
  assemble_newton_system (const TrilinosWrappers::MPI::Vector &linearization_point)
  {
    TimerOutput::Scope t(computing_timer, "Assembling");

    QGauss<dim> quadrature_formula(fe.degree + 1);
    QGauss<dim - 1> face_quadrature_formula(fe.degree + 1);

    FEValues<dim> fe_values(fe, quadrature_formula,
                            update_values | update_gradients | update_JxW_values);

    FEFaceValues<dim> fe_values_face(fe, face_quadrature_formula,
                                     update_values | update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell   = fe.dofs_per_cell;
    const unsigned int n_q_points      = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    const EquationData::BoundaryForce<dim> boundary_force;
    std::vector<Vector<double> >           boundary_force_values(n_face_q_points,
        Vector<double>(dim));

    FullMatrix<double>                     cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>                         cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index>   local_dof_indices(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    const FEValuesExtractors::Vector displacement(0);

    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          cell_matrix = 0;
          cell_rhs = 0;

          std::vector<SymmetricTensor<2, dim> > strain_tensor(n_q_points);
          fe_values[displacement].get_function_symmetric_gradients(linearization_point,
                                                                   strain_tensor);

          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            {
              SymmetricTensor<4, dim> stress_strain_tensor_linearized;
              SymmetricTensor<4, dim> stress_strain_tensor;
              constitutive_law.get_linearized_stress_strain_tensors(strain_tensor[q_point],
                                                                    stress_strain_tensor_linearized,
                                                                    stress_strain_tensor);

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  // Having computed the stress-strain tensor and its linearization,
                  // we can now put together the parts of the matrix and right hand side.
                  // In both, we need the linearized stress-strain tensor times the
                  // symmetric gradient of $\varphi_i$, i.e. the term $I_\Pi\varepsilon(\varphi_i)$,
                  // so we introduce an abbreviation of this term. Recall that the
                  // matrix corresponds to the bilinear form
                  // $A_{ij}=(I_\Pi\varepsilon(\varphi_i),\varepsilon(\varphi_j))$ in the
                  // notation of the accompanying publication, whereas the right
                  // hand side is $F_i=([I_\Pi-P_\Pi C]\varepsilon(\varphi_i),\varepsilon(\mathbf u))$
                  // where $u$ is the current linearization points (typically the last solution).
                  // This might suggest that the right hand side will be zero if the material
                  // is completely elastic (where $I_\Pi=P_\Pi$) but this ignores the fact
                  // that the right hand side will also contain contributions from
                  // non-homogeneous constraints due to the contact.
                  //
                  // The code block that follows this adds contributions that are due to
                  // boundary forces, should there be any.
                  const SymmetricTensor<2, dim>
                  stress_phi_i = stress_strain_tensor_linearized
                                 * fe_values[displacement].symmetric_gradient(i, q_point);

                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    cell_matrix(i, j) += (stress_phi_i
                                          * fe_values[displacement].symmetric_gradient(j, q_point)
                                          * fe_values.JxW(q_point));

                  cell_rhs(i) += ((stress_phi_i
                                   -
                                   stress_strain_tensor
                                   * fe_values[displacement].symmetric_gradient(i, q_point))
                                  * strain_tensor[q_point]
                                  * fe_values.JxW(q_point));
                }
            }

          for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
            if (cell->face(face)->at_boundary()
                &&
                cell->face(face)->boundary_indicator() == 1)
              {
                fe_values_face.reinit(cell, face);

                boundary_force.vector_value_list(fe_values_face.get_quadrature_points(),
                                                 boundary_force_values);

                for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
                  {
                    Tensor<1, dim> rhs_values;
                    rhs_values[2] = boundary_force_values[q_point][2];
                    for (unsigned int i = 0; i < dofs_per_cell; ++i)
                      cell_rhs(i) += (fe_values_face[displacement].value(i, q_point)
                                      * rhs_values
                                      * fe_values_face.JxW(q_point));
                  }
              }

          cell->get_dof_indices(local_dof_indices);
          all_constraints.distribute_local_to_global(cell_matrix, cell_rhs,
                                                     local_dof_indices,
                                                     newton_matrix,
                                                     newton_rhs,
                                                     true);

        }

    newton_matrix.compress(VectorOperation::add);
    newton_rhs.compress(VectorOperation::add);
  }



  // @sect4{PlasticityContactProblem::compute_nonlinear_residual}

  // The following function computes the nonlinear residual of the equation
  // given the current solution (or any other linearization point). This
  // is needed in the linear search algorithm where we need to try various
  // linear combinations of previous and current (trial) solution to
  // compute the (real, globalized) solution of the current Newton step.
  //
  // That said, in a slight abuse of the name of the function, it actually
  // does significantly more. For example, it also computes the vector
  // that corresponds to the Newton residual but without eliminating
  // constrained degrees of freedom. We need this vector to compute contact
  // forces and, ultimately, to compute the next active set. Likewise, by
  // keeping track of how many quadrature points we encounter on each cell
  // that show plastic yielding, we also compute the
  // <code>fraction_of_plastic_q_points_per_cell</code> vector that we
  // can later output to visualize the plastic zone. In both of these cases,
  // the results are not necessary as part of the line search, and so we may
  // be wasting a small amount of time computing them. At the same time, this
  // information appears as a natural by-product of what we need to do here
  // anyway, and we want to collect it once at the end of each Newton
  // step, so we may as well do it here.
  //
  // The actual implementation of this function should be rather obvious:
  template <int dim>
  void
  PlasticityContactProblem<dim>::
  compute_nonlinear_residual (const TrilinosWrappers::MPI::Vector &linearization_point)
  {
    QGauss<dim>   quadrature_formula(fe.degree + 1);
    QGauss<dim-1> face_quadrature_formula(fe.degree + 1);

    FEValues<dim> fe_values(fe, quadrature_formula,
                            update_values | update_gradients |
                            update_JxW_values);

    FEFaceValues<dim> fe_values_face(fe, face_quadrature_formula,
                                     update_values | update_quadrature_points |
                                     update_JxW_values);

    const unsigned int dofs_per_cell   = fe.dofs_per_cell;
    const unsigned int n_q_points      = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    const EquationData::BoundaryForce<dim> boundary_force;
    std::vector<Vector<double> >           boundary_force_values(n_face_q_points,
        Vector<double>(dim));

    Vector<double> cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    const FEValuesExtractors::Vector displacement(0);

    newton_rhs                            = 0;
    newton_rhs_uncondensed                = 0;

    fraction_of_plastic_q_points_per_cell = 0;

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    unsigned int cell_number = 0;
    for (; cell != endc; ++cell, ++cell_number)
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          cell_rhs = 0;

          std::vector<SymmetricTensor<2, dim> > strain_tensors(n_q_points);
          fe_values[displacement].get_function_symmetric_gradients(linearization_point,
                                                                   strain_tensors);

          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            {
              SymmetricTensor<4, dim> stress_strain_tensor;
              const bool q_point_is_plastic
                = constitutive_law.get_stress_strain_tensor(strain_tensors[q_point],
                                                            stress_strain_tensor);
              if (q_point_is_plastic)
                ++fraction_of_plastic_q_points_per_cell(cell_number);

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  cell_rhs(i) -= (strain_tensors[q_point]
                                  * stress_strain_tensor
                                  * fe_values[displacement].symmetric_gradient(i, q_point)
                                  * fe_values.JxW(q_point));

                  Tensor<1, dim> rhs_values;
                  rhs_values = 0;
                  cell_rhs(i) += (fe_values[displacement].value(i, q_point)
                                  * rhs_values
                                  * fe_values.JxW(q_point));
                }
            }

          for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell; ++face)
            if (cell->face(face)->at_boundary()
                && cell->face(face)->boundary_indicator() == 1)
              {
                fe_values_face.reinit(cell, face);

                boundary_force.vector_value_list(fe_values_face.get_quadrature_points(),
                                                 boundary_force_values);

                for (unsigned int q_point = 0; q_point < n_face_q_points;
                     ++q_point)
                  {
                    Tensor<1, dim> rhs_values;
                    rhs_values[2] = boundary_force_values[q_point][2];
                    for (unsigned int i = 0; i < dofs_per_cell; ++i)
                      cell_rhs(i) += (fe_values_face[displacement].value(i, q_point) * rhs_values
                                      * fe_values_face.JxW(q_point));
                  }
              }

          cell->get_dof_indices(local_dof_indices);
          constraints_dirichlet_and_hanging_nodes.distribute_local_to_global(cell_rhs,
              local_dof_indices,
              newton_rhs);

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            newton_rhs_uncondensed(local_dof_indices[i]) += cell_rhs(i);
        }

    fraction_of_plastic_q_points_per_cell /= quadrature_formula.size();
    newton_rhs.compress(VectorOperation::add);
    newton_rhs_uncondensed.compress(VectorOperation::add);
  }





  // @sect4{PlasticityContactProblem::solve_newton_system}

  // The last piece before we can discuss the actual Newton iteration
  // on a single mesh is the solver for the linear systems. There are
  // a couple of complications that slightly obscure the code, but
  // mostly it is just setup then solve. Among the complications are:
  //
  // - For the hanging nodes we have to apply
  //   the ConstraintMatrix::set_zero function to newton_rhs.
  //   This is necessary if a hanging node with solution value $x_0$
  //   has one neighbor with value $x_1$ which is in contact with the
  //   obstacle and one neighbor $x_2$ which is not in contact. Because
  //   the update for the former will be prescribed, the hanging node constraint
  //   will have an inhomogeneity and will look like $x_0 = x_1/2 + \text{gap}/2$.
  //   So the corresponding entries in the
  //   ride-hang-side are non-zero with a
  //   meaningless value. These values we have to
  //   to set to zero.
  // - Like in step-40, we need to shuffle between vectors that do and do
  //   do not have ghost elements when solving or using the solution.
  //
  // The rest of the function is similar to step-40 and
  // step-41 except that we use a BiCGStab solver
  // instead of CG. This is due to the fact that for very small hardening
  // parameters $\gamma$, the linear system becomes almost semidefinite though
  // still symmetric. BiCGStab appears to have an easier time with such linear
  // systems.
  template <int dim>
  void
  PlasticityContactProblem<dim>::solve_newton_system ()
  {
    TimerOutput::Scope t(computing_timer, "Solve");

    TrilinosWrappers::MPI::Vector distributed_solution(locally_owned_dofs, mpi_communicator);
    distributed_solution = solution;

    constraints_hanging_nodes.set_zero(distributed_solution);
    constraints_hanging_nodes.set_zero(newton_rhs);

    TrilinosWrappers::PreconditionAMG preconditioner;
    {
      TimerOutput::Scope t(computing_timer, "Solve: setup preconditioner");

      std::vector<std::vector<bool> > constant_modes;
      DoFTools::extract_constant_modes(dof_handler, ComponentMask(),
                                       constant_modes);

      TrilinosWrappers::PreconditionAMG::AdditionalData additional_data;
      additional_data.constant_modes = constant_modes;
      additional_data.elliptic = true;
      additional_data.n_cycles = 1;
      additional_data.w_cycle = false;
      additional_data.output_details = false;
      additional_data.smoother_sweeps = 2;
      additional_data.aggregation_threshold = 1e-2;

      preconditioner.initialize(newton_matrix, additional_data);
    }

    {
      TimerOutput::Scope t(computing_timer, "Solve: iterate");

      TrilinosWrappers::MPI::Vector tmp(locally_owned_dofs, mpi_communicator);

      const double relative_accuracy = 1e-8;
      const double solver_tolerance  = relative_accuracy
                                       * newton_matrix.residual(tmp, distributed_solution,
                                                                newton_rhs);

      SolverControl solver_control(newton_matrix.m(),
                                   solver_tolerance);
      SolverBicgstab<TrilinosWrappers::MPI::Vector> solver(solver_control);
      solver.solve(newton_matrix, distributed_solution,
                   newton_rhs, preconditioner);

      pcout << "         Error: " << solver_control.initial_value()
            << " -> " << solver_control.last_value() << " in "
            << solver_control.last_step() << " Bicgstab iterations."
            << std::endl;
    }

    all_constraints.distribute(distributed_solution);

    solution = distributed_solution;
  }


  // @sect4{PlasticityContactProblem::solve_newton}

  // This is, finally, the function that implements the damped Newton method
  // on the current mesh. There are two nested loops: the outer loop for the Newton
  // iteration and the inner loop for the line search which
  // will be used only if necessary. To obtain a good and reasonable
  // starting value we solve an elastic problem in the very first Newton step on each
  // mesh (or only on the first mesh if we transfer solutions between meshes). We
  // do so by setting the yield stress to an unreasonably large value in these
  // iterations and then setting it back to the correct value in subsequent
  // iterations.
  //
  // Other than this, the top part of this function should be
  // reasonably obvious. We initialize the variable
  // <code>previous_residual_norm</code> to the most negative value
  // representable with double precision numbers so that the
  // comparison whether the current residual is less than that of the
  // previous step will always fail in the first step.
  template <int dim>
  void
  PlasticityContactProblem<dim>::solve_newton ()
  {
    TrilinosWrappers::MPI::Vector old_solution(locally_owned_dofs, mpi_communicator);
    TrilinosWrappers::MPI::Vector residual(locally_owned_dofs, mpi_communicator);
    TrilinosWrappers::MPI::Vector tmp_vector(locally_owned_dofs, mpi_communicator);
    TrilinosWrappers::MPI::Vector locally_relevant_tmp_vector(locally_relevant_dofs, mpi_communicator);
    TrilinosWrappers::MPI::Vector distributed_solution(locally_owned_dofs, mpi_communicator);

    double residual_norm;
    double previous_residual_norm = -std::numeric_limits<double>::max();

    const double correct_sigma = sigma_0;

    IndexSet old_active_set(active_set);

    for (unsigned int newton_step = 1; newton_step <= 100; ++newton_step)
      {
        if (newton_step == 1
            &&
            ((transfer_solution && current_refinement_cycle == 0)
             ||
             !transfer_solution))
          constitutive_law.set_sigma_0(1e+10);
        else if (newton_step == 2
                 ||
                 current_refinement_cycle > 0
                 ||
                 !transfer_solution)
          constitutive_law.set_sigma_0(correct_sigma);

        pcout << " " << std::endl;
        pcout << "   Newton iteration " << newton_step << std::endl;
        pcout << "      Updating active set..." << std::endl;

        {
          TimerOutput::Scope t(computing_timer, "update active set");
          update_solution_and_constraints();
        }

        pcout << "      Assembling system... " << std::endl;
        newton_matrix = 0;
        newton_rhs = 0;
        assemble_newton_system(solution);

        pcout << "      Solving system... " << std::endl;
        solve_newton_system();

        // It gets a bit more hairy after we have computed the
        // trial solution $\tilde{\mathbf u}$ of the current Newton step.
        // We handle a highly nonlinear problem so we have to damp
        // Newton's method using a line search. To understand how we do this,
        // recall that in our formulation, we compute a trial solution
        // in each Newton step and not the update between old and new solution.
        // Since the solution set is a convex set, we will use a line
        // search that tries linear combinations of the
        // previous and the trial solution to guarantee that the
        // damped solution is in our solution set again.
        // At most we apply 5 damping steps.
        //
        // There are exceptions to when we use a line search. First,
        // if this is the first Newton step on any mesh, then we don't have
        // any point to compare the residual to, so we always accept a full
        // step. Likewise, if this is the second Newton step on the first mesh (or
        // the second on any mesh if we don't transfer solutions from
        // mesh to mesh), then we have computed the first of these steps using
        // just an elastic model (see how we set the yield stress sigma to
        // an unreasonably large value above). In this case, the first Newton
        // solution was a purely elastic one, the second one a plastic one,
        // and any linear combination would not necessarily be expected to
        // lie in the feasible set -- so we just accept the solution we just
        // got.
        //
        // In either of these two cases, we bypass the line search and just
        // update residual and other vectors as necessary.
        if ((newton_step==1)
            ||
            (transfer_solution && newton_step == 2 && current_refinement_cycle == 0)
            ||
            (!transfer_solution && newton_step == 2))
          {
            compute_nonlinear_residual(solution);
            old_solution = solution;

            residual = newton_rhs;
            const unsigned int start_res = (residual.local_range().first),
                               end_res = (residual.local_range().second);
            for (unsigned int n = start_res; n < end_res; ++n)
              if (all_constraints.is_inhomogeneously_constrained(n))
                residual(n) = 0;

            residual.compress(VectorOperation::insert);

            residual_norm = residual.l2_norm();

            pcout << "      Accepting Newton solution with residual: "
                  << residual_norm << std::endl;
          }
        else
          {
            for (unsigned int i = 0; i < 5; i++)
              {
                distributed_solution = solution;

                const double alpha = std::pow(0.5, static_cast<double>(i));
                tmp_vector = old_solution;
                tmp_vector.sadd(1 - alpha, alpha, distributed_solution);

                TimerOutput::Scope t(computing_timer, "Residual and lambda");

                locally_relevant_tmp_vector = tmp_vector;
                compute_nonlinear_residual(locally_relevant_tmp_vector);
                residual = newton_rhs;

                const unsigned int start_res = (residual.local_range().first),
                                   end_res = (residual.local_range().second);
                for (unsigned int n = start_res; n < end_res; ++n)
                  if (all_constraints.is_inhomogeneously_constrained(n))
                    residual(n) = 0;

                residual.compress(VectorOperation::insert);

                residual_norm = residual.l2_norm();

                pcout << "      Residual of the non-contact part of the system: "
                      << residual_norm << std::endl
                      << "         with a damping parameter alpha = " << alpha
                      << std::endl;

                if (residual_norm < previous_residual_norm)
                  break;
              }

            solution = tmp_vector;
            old_solution = solution;
          }

        old_active_set = active_set;
        previous_residual_norm = residual_norm;


        // The final step is to check for convergence. If the active set
        // has not changed across all processors and the residual is
        // less than a threshold of $10^{-10}$, then we terminate
        // the iteration on the current mesh:
        if (Utilities::MPI::sum((active_set == old_active_set) ? 0 : 1,
                                mpi_communicator) == 0)
          {
            pcout << "      Active set did not change!" << std::endl;
            if (residual_norm < 1e-10)
              break;
          }
      }
  }

  // @sect4{PlasticityContactProblem::refine_grid}

  // If you've made it this far into the deal.II tutorial, the following
  // function refining the mesh should not pose any challenges to you
  // any more. It refines the mesh, either globally or using the Kelly
  // error estimator, and if so asked also transfers the solution from
  // the previous to the next mesh. In the latter case, we also need
  // to compute the active set and other quantities again, for which we
  // need the information computed by <code>compute_nonlinear_residual()</code>.
  template <int dim>
  void
  PlasticityContactProblem<dim>::refine_grid ()
  {
    if (refinement_strategy == RefinementStrategy::refine_global)
      {
        for (typename Triangulation<dim>::active_cell_iterator
             cell = triangulation.begin_active();
             cell != triangulation.end(); ++cell)
          if (cell->is_locally_owned())
            cell->set_refine_flag ();
      }
    else
      {
        Vector<float> estimated_error_per_cell(triangulation.n_active_cells());
        KellyErrorEstimator<dim>::estimate(dof_handler,
                                           QGauss<dim - 1>(fe.degree + 2),
                                           typename FunctionMap<dim>::type(),
                                           solution,
                                           estimated_error_per_cell);

        parallel::distributed::GridRefinement
        ::refine_and_coarsen_fixed_number(triangulation,
                                          estimated_error_per_cell,
                                          0.3, 0.03);
      }

    triangulation.prepare_coarsening_and_refinement();

    parallel::distributed::SolutionTransfer<dim,
             TrilinosWrappers::MPI::Vector> solution_transfer(dof_handler);
    if (transfer_solution)
      solution_transfer.prepare_for_coarsening_and_refinement(solution);

    triangulation.execute_coarsening_and_refinement();

    setup_system();

    if (transfer_solution)
      {
        TrilinosWrappers::MPI::Vector distributed_solution(locally_owned_dofs, mpi_communicator);
        solution_transfer.interpolate(distributed_solution);
        solution = distributed_solution;
        compute_nonlinear_residual(solution);
      }
  }


  // @sect4{PlasticityContactProblem::move_mesh}

  // The remaining three functions before we get to <code>run()</code>
  // have to do with generating output. The following one is an attempt
  // at showing the deformed body in its deformed configuration. To this
  // end, this function takes a displacement vector field and moves every
  // vertex of the (local part) of the mesh by the previously computed
  // displacement. We will call this function with the current
  // displacement field before we generate graphical output, and we will
  // call it again after generating graphical output with the negative
  // displacement field to undo the changes to the mesh so made.
  //
  // The function itself is pretty straightforward. All we have to do
  // is keep track which vertices we have already touched, as we
  // encounter the same vertices multiple times as we loop over cells.
  template <int dim>
  void
  PlasticityContactProblem<dim>::
  move_mesh (const TrilinosWrappers::MPI::Vector &displacement) const
  {
    std::vector<bool> vertex_touched(triangulation.n_vertices(), false);

    for (typename DoFHandler<dim>::active_cell_iterator cell =
           dof_handler.begin_active();
         cell != dof_handler.end(); ++cell)
      if (cell->is_locally_owned())
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
          if (vertex_touched[cell->vertex_index(v)] == false)
            {
              vertex_touched[cell->vertex_index(v)] = true;

              Point<dim> vertex_displacement;
              for (unsigned int d = 0; d < dim; ++d)
                vertex_displacement[d] = displacement(cell->vertex_dof_index(v, d));

              cell->vertex(v) += vertex_displacement;
            }
  }



  // @sect4{PlasticityContactProblem::output_results}

  // Next is the function we use to actually generate graphical output. The
  // function is a bit tedious, but not actually particularly complicated.
  // It moves the mesh at the top (and moves it back at the end), then
  // computes the contact forces along the contact surface. We can do
  // so (as shown in the accompanying paper) by taking the untreated
  // residual vector and identifying which degrees of freedom
  // correspond to those with contact by asking whether they have an
  // inhomogeneous constraints associated with them. As always, we need
  // to be mindful that we can only write into completely distributed
  // vectors (i.e., vectors without ghost elements) but that when we
  // want to generate output, we need vectors that do indeed have
  // ghost entries for all locally relevant degrees of freedom.
  template <int dim>
  void
  PlasticityContactProblem<dim>::output_results (const std::string &filename_base)
  {
    TimerOutput::Scope t(computing_timer, "Graphical output");

    pcout << "      Writing graphical output... " << std::flush;

    move_mesh(solution);

    // Calculation of the contact forces
    TrilinosWrappers::MPI::Vector distributed_lambda(locally_owned_dofs, mpi_communicator);
    const unsigned int start_res = (newton_rhs_uncondensed.local_range().first),
                       end_res = (newton_rhs_uncondensed.local_range().second);
    for (unsigned int n = start_res; n < end_res; ++n)
      if (all_constraints.is_inhomogeneously_constrained(n))
        distributed_lambda(n) = newton_rhs_uncondensed(n) /
                                diag_mass_matrix_vector(n);
    distributed_lambda.compress(VectorOperation::insert);
    constraints_hanging_nodes.distribute(distributed_lambda);

    TrilinosWrappers::MPI::Vector lambda(locally_relevant_dofs, mpi_communicator);
    lambda = distributed_lambda;


    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);

    const std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(dim, DataComponentInterpretation::component_is_part_of_vector);
    data_out.add_data_vector(solution,
                             std::vector<std::string> (dim, "displacement"),
                             DataOut<dim>::type_dof_data, data_component_interpretation);
    data_out.add_data_vector(lambda,
                             std::vector<std::string> (dim, "contact_force"),
                             DataOut<dim>::type_dof_data, data_component_interpretation);
    data_out.add_data_vector(active_set,
                             std::vector<std::string> (dim, "active_set"),
                             DataOut<dim>::type_dof_data, data_component_interpretation);

    Vector<float> subdomain(triangulation.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector(subdomain, "subdomain");

    data_out.add_data_vector(fraction_of_plastic_q_points_per_cell,
                             "fraction_of_plastic_q_points");

    data_out.build_patches();

    // In the remainder of the function, we generate one VTU file on
    // every processor, indexed by the subdomain id of this processor.
    // On the first processor, we then also create a <code>.pvtu</code>
    // file that indexes <i>all</i> of the VTU files so that the entire
    // set of output files can be read at once. These <code>.pvtu</code>
    // are used by Paraview to describe an entire parallel computation's
    // output files. We then do the same again for the competitor of
    // Paraview, the Visit visualization program, by creating a matching
    // <code>.visit</code> file.
    const std::string filename =
      (output_dir + filename_base + "-"
       + Utilities::int_to_string(triangulation.locally_owned_subdomain(), 4));

    std::ofstream output_vtu((filename + ".vtu").c_str());
    data_out.write_vtu(output_vtu);
    pcout << output_dir + filename_base << ".pvtu" << std::endl;

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        std::vector<std::string> filenames;
        for (unsigned int i = 0;
             i < Utilities::MPI::n_mpi_processes(mpi_communicator); ++i)
          filenames.push_back(filename_base + "-" +
                              Utilities::int_to_string(i, 4) +
                              ".vtu");

        std::ofstream pvtu_master_output((output_dir + filename_base + ".pvtu").c_str());
        data_out.write_pvtu_record(pvtu_master_output, filenames);

        std::ofstream visit_master_output((output_dir + filename_base + ".visit").c_str());
        data_out.write_visit_record(visit_master_output, filenames);
      }

    TrilinosWrappers::MPI::Vector tmp(solution);
    tmp *= -1;
    move_mesh(tmp);
  }


  // @sect4{PlasticityContactProblem::output_contact_force}

  // This last auxiliary function computes the contact force by
  // calculating an integral over the contact pressure in z-direction
  // over the contact area. For this purpose we set the contact
  // pressure lambda to 0 for all inactive dofs (whether a degree
  // of freedom is part of the contact is determined just as
  // we did in the previous function). For all
  // active dofs, lambda contains the quotient of the nonlinear
  // residual (newton_rhs_uncondensed) and corresponding diagonal entry
  // of the mass matrix (diag_mass_matrix_vector). Because it is
  // not unlikely that hanging nodes show up in the contact area
  // it is important to apply contraints_hanging_nodes.distribute
  // to the distributed_lambda vector.
  template <int dim>
  void
  PlasticityContactProblem<dim>::output_contact_force () const
  {
    TrilinosWrappers::MPI::Vector distributed_lambda(locally_owned_dofs, mpi_communicator);
    const unsigned int start_res = (newton_rhs_uncondensed.local_range().first),
                       end_res = (newton_rhs_uncondensed.local_range().second);
    for (unsigned int n = start_res; n < end_res; ++n)
      if (all_constraints.is_inhomogeneously_constrained(n))
        distributed_lambda(n) = newton_rhs_uncondensed(n) / diag_mass_matrix_vector(n);
      else
        distributed_lambda(n) = 0;
    distributed_lambda.compress(VectorOperation::insert);
    constraints_hanging_nodes.distribute(distributed_lambda);

    TrilinosWrappers::MPI::Vector lambda(locally_relevant_dofs, mpi_communicator);
    lambda = distributed_lambda;

    double contact_force = 0.0;

    QGauss<dim-1> face_quadrature_formula(fe.degree + 1);
    FEFaceValues<dim> fe_values_face(fe, face_quadrature_formula,
                                     update_values | update_JxW_values);

    const unsigned int n_face_q_points = face_quadrature_formula.size();

    const FEValuesExtractors::Vector displacement(0);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
        for (unsigned int face = 0; face<GeometryInfo<dim>::faces_per_cell; ++face)
          if (cell->face(face)->at_boundary()
              &&
              cell->face(face)->boundary_indicator() == 1)
            {
              fe_values_face.reinit(cell, face);

              std::vector<Tensor<1, dim> > lambda_values(n_face_q_points);
              fe_values_face[displacement].get_function_values(lambda,
                                                               lambda_values);

              for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
                contact_force += lambda_values[q_point][2]
                                 * fe_values_face.JxW(q_point);
            }
    contact_force = Utilities::MPI::sum(contact_force, MPI_COMM_WORLD);

    pcout << "Contact force = " << contact_force << std::endl;
  }


  // @sect4{PlasticityContactProblem::run}

  // As in all other tutorial programs, the <code>run()</code> function contains
  // the overall logic. There is not very much to it here: in essence, it
  // performs the loops over all mesh refinement cycles, and within each, hands
  // things over to the Newton solver in <code>solve_newton()</code> on the
  // current mesh and calls the function that creates graphical output for
  // the so-computed solution. It then outputs some statistics concerning both
  // run times and memory consumption that has been collected over the course of
  // computations on this mesh.
  template <int dim>
  void
  PlasticityContactProblem<dim>::run ()
  {
    computing_timer.reset();
    for (current_refinement_cycle = 0;
         current_refinement_cycle < n_refinement_cycles;
         ++current_refinement_cycle)
      {
        {
          TimerOutput::Scope t(computing_timer, "Setup");

          pcout << std::endl;
          pcout << "Cycle " << current_refinement_cycle << ':' << std::endl;

          if (current_refinement_cycle == 0)
            {
              make_grid();
              setup_system();
            }
          else
            {
              TimerOutput::Scope t(computing_timer, "Setup: refine mesh");
              refine_grid();
            }
        }

        solve_newton();

        output_results((std::string("solution-") +
                        Utilities::int_to_string(current_refinement_cycle, 2)).c_str());

        computing_timer.print_summary();
        computing_timer.reset();

        Utilities::System::MemoryStats stats;
        Utilities::System::get_memory_stats(stats);
        pcout << "Peak virtual memory used, resident in kB: " << stats.VmSize << " "
              << stats.VmRSS << std::endl;

        if (base_mesh == "box")
          output_contact_force();
      }
  }
}

// @sect3{The <code>main</code> function}

// There really isn't much to the <code>main()</code> function. It looks
// like they always do:
int main (int argc, char *argv[])
{
  using namespace dealii;
  using namespace Step42;

  try
    {
      deallog.depth_console(0);
      ParameterHandler prm;
      PlasticityContactProblem<3>::declare_parameters(prm);
      if (argc != 2)
        {
          std::cerr << "*** Call this program as <./step-42 input.prm>" << std::endl;
          return 1;
        }

      prm.read_input(argv[1]);
      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
      {
        PlasticityContactProblem<3> problem(prm);
        problem.run();
      }
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
