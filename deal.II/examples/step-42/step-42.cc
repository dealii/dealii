/* ---------------------------------------------------------------------
 * $Id$
 *
 * Copyright (C) 2012 - 2013 by the deal.II authors
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
    // provides in this way a function which describes an obstacle.
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

    template <int dim>
    double
    BitmapFile<dim>::get_pixel_value(const int i, const int j) const
    {
      assert(i >= 0 && i < nx);
      assert(j >= 0 && j < ny);
      return obstacle_data[nx * (ny - 1 - j) + i];
    }

    template <int dim>
    double
    BitmapFile<dim>::get_value(const double x, const double y) const
    {
      const int ix = std::min(std::max((int) (x / hx), 0), nx - 2);
      const int iy = std::min(std::max((int) (y / hy), 0), ny - 2);

      FullMatrix<double> H(4, 4);
      Vector<double> X(4);
      Vector<double> b(4);

      double xx, yy;

      xx = ix * hx;
      yy = iy * hy;
      H(0, 0) = xx;
      H(0, 1) = yy;
      H(0, 2) = xx * yy;
      H(0, 3) = 1.0;
      b(0) = get_pixel_value(ix, iy);

      xx = (ix + 1) * hx;
      yy = iy * hy;
      H(1, 0) = xx;
      H(1, 1) = yy;
      H(1, 2) = xx * yy;
      H(1, 3) = 1.0;
      b(1) = get_pixel_value(ix + 1, iy);

      xx = (ix + 1) * hx;
      yy = (iy + 1) * hy;
      H(2, 0) = xx;
      H(2, 1) = yy;
      H(2, 2) = xx * yy;
      H(2, 3) = 1.0;
      b(2) = get_pixel_value(ix + 1, iy + 1);

      xx = ix * hx;
      yy = (iy + 1) * hy;
      H(3, 0) = xx;
      H(3, 1) = yy;
      H(3, 2) = xx * yy;
      H(3, 3) = 1.0;
      b(3) = get_pixel_value(ix, iy + 1);

      H.gauss_jordan();
      H.vmult(X, b);

      return (X(0) * x + X(1) * y + X(2) * x * y + X(3));
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
    void assemble_nl_system (const TrilinosWrappers::MPI::Vector &u);
    void compute_nonlinear_residual (const TrilinosWrappers::MPI::Vector &current_solution);
    void assemble_mass_matrix_diagonal (TrilinosWrappers::SparseMatrix &mass_matrix);
    void update_solution_and_constraints ();
    void dirichlet_constraints ();
    void solve ();
    void solve_newton ();
    void refine_grid ();
    void move_mesh (const TrilinosWrappers::MPI::Vector &_complete_displacement) const;
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
    TrilinosWrappers::SparseMatrix    system_matrix_newton;

    TrilinosWrappers::MPI::Vector     solution;
    TrilinosWrappers::MPI::Vector     system_rhs_newton;
    TrilinosWrappers::MPI::Vector     system_rhs_lambda;
    TrilinosWrappers::MPI::Vector     resid_vector;
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
    const std_cxx1x::shared_ptr<const Function<dim> >  obstacle;

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
    prm.declare_entry("polynomial degree", "1", Patterns::Integer(),
                      "Polynomial degree of the FE_Q finite element space, typically 1 or 2.");
    prm.declare_entry("number of initial refinements", "2",
                      Patterns::Integer(),
                      "Number of initial global mesh refinement steps before "
                      "the first computation.");
    prm.declare_entry("refinement strategy", "percentage",
                      Patterns::Selection("global|percentage|fix dofs"),
                      "Mesh refinement strategy:\n"
                      " global: one global refinement\n"
                      " percentage: fixed percentage gets refined using the Kelly estimator\n"
                      " fix dofs: tries to achieve 2^initial_refinement*300 dofs after "
                      "refinement cycle 1 (only use 2 cycles!). This requires the code to "
                      "make some changes to the coarse mesh as well.");
    prm.declare_entry("number of cycles", "5", Patterns::Integer(),
                      "Number of adaptive mesh refinement cycles to run.");
    prm.declare_entry("obstacle filename", "", Patterns::Anything(),
                      "Obstacle file to read, use 'obstacle_file.pbm' or leave empty to use a sphere.");
    prm.declare_entry("output directory", "", Patterns::Anything(),
                      "Directory for output files (graphical output and benchmark "
                      "statistics). If empty, use the current directory.");
    prm.declare_entry("transfer solution", "false", Patterns::Bool(),
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
    obstacle (prm.get("obstacle filename") != ""
              ?
              static_cast<const Function<dim>*>
              (new EquationData::ChineseObstacle<dim>(prm.get("obstacle filename"), (base_mesh == "box" ? 1.0 : 0.5)))
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
  // it a bit before using it:
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



  template <int dim>
  void
  PlasticityContactProblem<dim>::setup_system ()
  {
    // setup dofs
    {
      TimerOutput::Scope t(computing_timer, "Setup: distribute DoFs");
      dof_handler.distribute_dofs(fe);

      locally_owned_dofs = dof_handler.locally_owned_dofs();
      locally_relevant_dofs.clear();
      DoFTools::extract_locally_relevant_dofs(dof_handler,
                                              locally_relevant_dofs);
    }

    // setup hanging nodes and Dirichlet constraints
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

      dirichlet_constraints();
    }

    // Initialization for matrices and vectors
    {
      TimerOutput::Scope t(computing_timer, "Setup: vectors");
      solution.reinit(locally_relevant_dofs, mpi_communicator);
      system_rhs_newton.reinit(locally_owned_dofs, mpi_communicator);
      system_rhs_lambda.reinit(system_rhs_newton);
      resid_vector.reinit(system_rhs_newton);
      diag_mass_matrix_vector.reinit(system_rhs_newton);
      fraction_of_plastic_q_points_per_cell.reinit(triangulation.n_active_cells());
      active_set.clear();
      active_set.set_size(locally_relevant_dofs.size());
    }

    // setup sparsity pattern
    {
      TimerOutput::Scope t(computing_timer, "Setup: matrix");
      TrilinosWrappers::SparsityPattern sp(locally_owned_dofs,
                                           mpi_communicator);

      DoFTools::make_sparsity_pattern(dof_handler, sp,
                                      constraints_dirichlet_and_hanging_nodes, false,
                                      Utilities::MPI::this_mpi_process(mpi_communicator));

      sp.compress();

      system_matrix_newton.reinit(sp);

      // we are going to reuse the system
      // matrix for assembling the diagonal
      // of the mass matrix so that we do not
      // need to allocate two sparse matrices
      // at the same time:
      TrilinosWrappers::SparseMatrix &mass_matrix = system_matrix_newton;
      assemble_mass_matrix_diagonal(mass_matrix);
      const unsigned int start = (system_rhs_newton.local_range().first),
                         end = (system_rhs_newton.local_range().second);
      for (unsigned int j = start; j < end; j++)
        diag_mass_matrix_vector(j) = mass_matrix.diag_element(j);

      diag_mass_matrix_vector.compress(VectorOperation::insert);

      // remove the mass matrix entries from the matrix:
      mass_matrix = 0;
    }
  }

  template <int dim>
  void
  PlasticityContactProblem<dim>::assemble_nl_system (const TrilinosWrappers::MPI::Vector &u)
  {
    TimerOutput::Scope t(computing_timer, "Assembling");

    QGauss<dim> quadrature_formula(fe.degree + 1);
    QGauss<dim - 1> face_quadrature_formula(fe.degree + 1);

    FEValues<dim> fe_values(fe, quadrature_formula,
                            UpdateFlags(
                              update_values | update_gradients | update_q_points
                              | update_JxW_values));

    FEFaceValues<dim> fe_values_face(fe, face_quadrature_formula,
                                     update_values | update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    const EquationData::BoundaryForce<dim> boundary_force;
    std::vector<Vector<double> > boundary_force_values(n_face_q_points,
                                                       Vector<double>(dim));

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double> cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell =
      dof_handler.begin_active(), endc = dof_handler.end();

    const FEValuesExtractors::Vector displacement(0);

    const double kappa = 1.0;
    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          cell_matrix = 0;
          cell_rhs = 0;

          std::vector<SymmetricTensor<2, dim> > strain_tensor(n_q_points);
          fe_values[displacement].get_function_symmetric_gradients(u,
                                                                   strain_tensor);

          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            {
              SymmetricTensor<4, dim> stress_strain_tensor_linearized;
              SymmetricTensor<4, dim> stress_strain_tensor;
              SymmetricTensor<2, dim> stress_tensor;

              constitutive_law.get_linearized_stress_strain_tensors(strain_tensor[q_point],
                                                                    stress_strain_tensor_linearized,
                                                                    stress_strain_tensor);

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  stress_tensor = stress_strain_tensor_linearized
                                  * fe_values[displacement].symmetric_gradient(i, q_point);

                  for (unsigned int j = 0; j < dofs_per_cell; ++j)
                    cell_matrix(i, j) += (stress_tensor
                                          * fe_values[displacement].symmetric_gradient(j, q_point)
                                          * fe_values.JxW(q_point));

                  // the linearized part a(v^i;v^i,v) of the rhs
                  cell_rhs(i) += (stress_tensor * strain_tensor[q_point]
                                  * fe_values.JxW(q_point));

                  // the residual part a(v^i;v) of the rhs
                  cell_rhs(i) -= (strain_tensor[q_point]
                                  * stress_strain_tensor
                                  * fe_values[displacement].symmetric_gradient(i, q_point)
                                  * fe_values.JxW(q_point));

                  // the residual part F(v) of the rhs
                  Tensor<1, dim> rhs_values;
                  rhs_values = 0;
                  cell_rhs(i) += (fe_values[displacement].value(i, q_point)
                                  * rhs_values * fe_values.JxW(q_point));
                }
            }

          for (unsigned int face = 0;
               face < GeometryInfo<dim>::faces_per_cell; ++face)
            {
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
                        cell_rhs(i) += (fe_values_face[displacement].value(i,
                                                                           q_point) * rhs_values
                                        * fe_values_face.JxW(q_point));
                    }
                }
            }

          cell->get_dof_indices(local_dof_indices);
          all_constraints.distribute_local_to_global(cell_matrix, cell_rhs,
                                                     local_dof_indices,
                                                     system_matrix_newton,
                                                     system_rhs_newton,
                                                     true);

        }

    system_matrix_newton.compress(VectorOperation::add);
    system_rhs_newton.compress(VectorOperation::add);
  }



  template <int dim>
  void
  PlasticityContactProblem<dim>::compute_nonlinear_residual (const TrilinosWrappers::MPI::Vector &current_solution)
  {
    QGauss<dim> quadrature_formula(fe.degree + 1);
    QGauss<dim-1> face_quadrature_formula(fe.degree + 1);

    FEValues<dim> fe_values(fe, quadrature_formula,
                            update_values | update_gradients |
                            update_q_points  | update_JxW_values);

    FEFaceValues<dim> fe_values_face(fe, face_quadrature_formula,
                                     update_values | update_quadrature_points |
                                     update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_q_points = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    const EquationData::BoundaryForce<dim> boundary_force;
    std::vector<Vector<double> > boundary_force_values(n_face_q_points,
                                                       Vector<double>(dim));

    Vector<double> cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    const FEValuesExtractors::Vector displacement(0);

    typename DoFHandler<dim>::active_cell_iterator cell =
      dof_handler.begin_active(), endc = dof_handler.end();

    unsigned int elast_points = 0;
    unsigned int plast_points = 0;
    double yield = 0;
    unsigned int cell_number = 0;
    fraction_of_plastic_q_points_per_cell = 0;

    for (; cell != endc; ++cell, ++cell_number)
      if (cell->is_locally_owned())
        {
          fe_values.reinit(cell);
          cell_rhs = 0;

          std::vector<SymmetricTensor<2, dim> > strain_tensor(n_q_points);
          fe_values[displacement].get_function_symmetric_gradients(current_solution,
                                                                   strain_tensor);

          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            {
              SymmetricTensor<4, dim> stress_strain_tensor;
              const bool q_point_is_plastic
                = constitutive_law.get_stress_strain_tensor(strain_tensor[q_point],
                                                            stress_strain_tensor);
              if (q_point_is_plastic)
                {
                  ++plast_points;
                  ++fraction_of_plastic_q_points_per_cell(cell_number);
                }
              else
                ++elast_points;

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  cell_rhs(i) -= (strain_tensor[q_point]
                                  * stress_strain_tensor
                                  * fe_values[displacement].symmetric_gradient(i, q_point)
                                  * fe_values.JxW(q_point));

                  Tensor<1, dim> rhs_values;
                  rhs_values = 0;
                  cell_rhs(i) += (fe_values[displacement].value(i, q_point)
                                  * rhs_values * fe_values.JxW(q_point));
                }
            }

          for (unsigned int face = 0;
               face < GeometryInfo<dim>::faces_per_cell; ++face)
            {
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
            }

          cell->get_dof_indices(local_dof_indices);
          constraints_dirichlet_and_hanging_nodes.distribute_local_to_global(cell_rhs,
              local_dof_indices,
              system_rhs_newton);

          for (unsigned int i = 0; i < dofs_per_cell; i++)
            system_rhs_lambda(local_dof_indices[i]) += cell_rhs(i);
        }

    fraction_of_plastic_q_points_per_cell /= quadrature_formula.size();
    system_rhs_newton.compress(VectorOperation::add);
    system_rhs_lambda.compress(VectorOperation::add);

    const unsigned int sum_elast_points = Utilities::MPI::sum(elast_points,
                                                              mpi_communicator);
    const unsigned int sum_plast_points = Utilities::MPI::sum(plast_points,
                                                              mpi_communicator);
    pcout << "      Number of elastic quadrature points: " << sum_elast_points
          << " and plastic quadrature points: " << sum_plast_points
          << std::endl;
  }



  template <int dim>
  void
  PlasticityContactProblem<dim>::assemble_mass_matrix_diagonal (TrilinosWrappers::SparseMatrix &mass_matrix)
  {
    QGaussLobatto<dim - 1> face_quadrature_formula(fe.degree + 1);

    FEFaceValues<dim> fe_values_face(fe, face_quadrature_formula,
                                     update_values | update_quadrature_points | update_JxW_values);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Tensor<1, dim, double> ones(dim);
    for (unsigned i = 0; i < dim; i++)
      ones[i] = 1.0;

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    const FEValuesExtractors::Vector displacement(0);

    typename DoFHandler<dim>::active_cell_iterator cell =
      dof_handler.begin_active(), endc = dof_handler.end();

    for (; cell != endc; ++cell)
      if (cell->is_locally_owned())
        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
             ++face)
          if (cell->face(face)->at_boundary()
              && cell->face(face)->boundary_indicator() == 1)
            {
              fe_values_face.reinit(cell, face);
              cell_matrix = 0;

              for (unsigned int q_point = 0; q_point < n_face_q_points;
                   ++q_point)
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  cell_matrix(i, i) += (fe_values_face[displacement].value(i,
                                                                           q_point) * ones * fe_values_face.JxW(q_point));

              cell->get_dof_indices(local_dof_indices);

              for (unsigned int i = 0; i < dofs_per_cell; i++)
                mass_matrix.add(local_dof_indices[i], local_dof_indices[i],
                                cell_matrix(i, i));
            }
    mass_matrix.compress(VectorOperation::add);
  }

// @sect4{PlasticityContactProblem::update_solution_and_constraints}

// Projection and updating of the active set
// for the dofs which penetrates the obstacle.
  template <int dim>
  void
  PlasticityContactProblem<dim>::update_solution_and_constraints ()
  {
    std::vector<bool> vertex_touched(dof_handler.n_dofs(), false);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    TrilinosWrappers::MPI::Vector distributed_solution(system_rhs_newton);
    distributed_solution = solution;
    TrilinosWrappers::MPI::Vector lambda(solution);
    lambda = resid_vector;
    TrilinosWrappers::MPI::Vector diag_mass_matrix_vector_relevant(solution);
    diag_mass_matrix_vector_relevant = diag_mass_matrix_vector;

    all_constraints.reinit(locally_relevant_dofs);
    active_set.clear();
    IndexSet active_set_locally_owned;
    active_set_locally_owned.set_size(locally_owned_dofs.size());
    const double c = 100.0 * e_modulus;

    Quadrature<dim - 1> face_quadrature(fe.get_unit_face_support_points());
    FEFaceValues<dim> fe_values_face(fe, face_quadrature,
                                     update_quadrature_points);

    const unsigned int dofs_per_face = fe.dofs_per_face;
    const unsigned int n_face_q_points = face_quadrature.size();

    std::vector<types::global_dof_index> dof_indices(dofs_per_face);

    // pcout<< "dofs_per_face = " << dofs_per_face
    //      << "n_face_q_points = " << n_face_q_points
    //      <<std::endl;
    unsigned int counter_hanging_nodes = 0;
    for (; cell != endc; ++cell)
      if (!cell->is_artificial())
        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
             ++face)
          if (cell->face(face)->at_boundary()
              && cell->face(face)->boundary_indicator() == 1)
            {
              fe_values_face.reinit(cell, face);
              cell->face(face)->get_dof_indices(dof_indices);

              for (unsigned int q_point = 0; q_point < n_face_q_points;
                   ++q_point)
                {
                  unsigned int component = fe.face_system_to_component_index(
                                             q_point).first;

                  if (component == 2)
                    {
                      unsigned int index_z = dof_indices[q_point];

                      if (vertex_touched[index_z] == false)
                        vertex_touched[index_z] = true;
                      else
                        continue;

                      // the local row where
                      Point<dim> point(
                        fe_values_face.quadrature_point(q_point));

                      double obstacle_value = obstacle->value(point, 2);
                      double solution_index_z = solution(index_z);
                      double gap = obstacle_value - point(2);

                      if (lambda(index_z)
                          / diag_mass_matrix_vector_relevant(index_z)
                          + c * (solution_index_z - gap) > 0
                          && !(constraints_hanging_nodes.is_constrained(index_z)))
                        {
                          all_constraints.add_line(index_z);
                          all_constraints.set_inhomogeneity(index_z, gap);
                          distributed_solution(index_z) = gap;

                          if (locally_owned_dofs.is_element(index_z))
                            {
                              active_set_locally_owned.add_index(index_z);
                              if (locally_relevant_dofs.is_element(index_z))
                                active_set.add_index(index_z);
                            }

                        }
                      else if (lambda(index_z)
                               / diag_mass_matrix_vector_relevant(index_z)
                               + c * (solution_index_z - gap) > 0
                               && constraints_hanging_nodes.is_constrained(index_z))
                        {
                          if (locally_owned_dofs.is_element(index_z))
                            counter_hanging_nodes += 1;
                        }
                    }
                }
            }
    distributed_solution.compress(VectorOperation::insert);

    const unsigned int sum_contact_constraints
      = Utilities::MPI::sum(active_set_locally_owned.n_elements(),
                            mpi_communicator);
    pcout << "         Size of active set: " << sum_contact_constraints
          << std::endl;
    const unsigned int sum_contact_hanging_nodes
      = Utilities::MPI::sum(counter_hanging_nodes,
                            mpi_communicator);
    pcout << "         Number of hanging nodes in contact: "
          << sum_contact_hanging_nodes << std::endl;

    solution = distributed_solution;

    all_constraints.close();
    all_constraints.merge(constraints_dirichlet_and_hanging_nodes);
  }

// @sect4{PlasticityContactProblem::dirichlet_constraints}

// This function defines the new ConstraintMatrix
// constraints_dirichlet_hanging_nodes. It contains
// the Dirichlet boundary values as well as the
// hanging nodes constraints.
  template <int dim>
  void
  PlasticityContactProblem<dim>::dirichlet_constraints ()
  {
    /* boundary_indicators:
     _______
     /  1    /|
     /______ / |
     8|       | 8|
     |   8   | /
     |_______|/
     6
     */

    constraints_dirichlet_and_hanging_nodes.reinit(locally_relevant_dofs);
    constraints_dirichlet_and_hanging_nodes.merge(constraints_hanging_nodes);

    // interpolate all components of the solution
    VectorTools::interpolate_boundary_values(dof_handler,
                                             base_mesh == "box" ? 6 : 0, EquationData::BoundaryValues<dim>(),
                                             constraints_dirichlet_and_hanging_nodes, ComponentMask());

    // interpolate x- and y-components of the
    // solution (this is a bit mask, so apply
    // operator| )
    const FEValuesExtractors::Scalar x_displacement(0);
    const FEValuesExtractors::Scalar y_displacement(1);
    VectorTools::interpolate_boundary_values(dof_handler, 8,
                                             EquationData::BoundaryValues<dim>(),
                                             constraints_dirichlet_and_hanging_nodes,
                                             (fe.component_mask(x_displacement) | fe.component_mask(y_displacement)));
    constraints_dirichlet_and_hanging_nodes.close();
  }

// @sect4{PlasticityContactProblem::solve}

// In addition to step-41 we have
// to deal with the hanging node
// constraints. Again we also consider
// the locally_owned_dofs only by
// creating the vector distributed_solution.
//
// For the hanging nodes we have to apply
// the set_zero function to system_rhs_newton.
// This is necessary if a hanging node value x_0
// has one neighbor which is in contact with
// value x_0 and one neighbor which is not with
// value x_1. This leads to an inhomogeneity
// constraint with value x_1/2 = gap/2 in the
// ConstraintMatrix.
// So the corresponding entries in the
// ride-hang-side are non-zero with a
// meaningless value. These values have to
// to set to zero.

// The rest of the function is similar to
// step-41 except that we use a FGMRES-solver
// instead of CG. For a very small hardening
// value gamma the linear system becomes
// almost semi definite but still symmetric.
  template <int dim>
  void
  PlasticityContactProblem<dim>::solve ()
  {
    TimerOutput::Scope t(computing_timer, "Solve");

    TrilinosWrappers::MPI::Vector distributed_solution(system_rhs_newton);
    distributed_solution = solution;

    constraints_hanging_nodes.set_zero(distributed_solution);
    constraints_hanging_nodes.set_zero(system_rhs_newton);
    distributed_solution.compress(VectorOperation::insert);
    system_rhs_newton.compress(VectorOperation::insert);

    TrilinosWrappers::PreconditionAMG preconditioner;
    {
      TimerOutput::Scope t(computing_timer, "Solve: setup preconditioner");

      std::vector < std::vector<bool> > constant_modes;
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

      preconditioner.initialize(system_matrix_newton, additional_data);
    }

    {
      TimerOutput::Scope t(computing_timer, "Solve: iterate");

      PrimitiveVectorMemory<TrilinosWrappers::MPI::Vector> mem;
      TrilinosWrappers::MPI::Vector tmp(system_rhs_newton);
      // 1e-4 seems to be the fasted option altogether, but to get more
      // reproducible parallel benchmark results, we use a small residual:
      double relative_accuracy = 1e-8;
      if (output_dir.compare("its/") == 0)
        relative_accuracy = 1e-4;

      const double solver_tolerance = relative_accuracy
                                      * system_matrix_newton.residual(tmp, distributed_solution,
                                          system_rhs_newton);

      SolverControl solver_control(system_matrix_newton.m(),
                                   solver_tolerance);
      SolverBicgstab<TrilinosWrappers::MPI::Vector> solver(solver_control,
                                                           mem/*,
               SolverFGMRES<TrilinosWrappers::MPI::Vector>::
               AdditionalData(30, true)*/);
      solver.solve(system_matrix_newton, distributed_solution,
                   system_rhs_newton, preconditioner);

      pcout << "         Error: " << solver_control.initial_value()
            << " -> " << solver_control.last_value() << " in "
            << solver_control.last_step() << " Bicgstab iterations."
            << std::endl;
    }

    all_constraints.distribute(distributed_solution);

    solution = distributed_solution;
  }

// @sect4{PlasticityContactProblem::solve_newton}

// In this function the damped Newton method is implemented.
// That means two nested loops: the outer loop for the newton
// iteration and the inner loop for the damping steps which
// will be used only if necessary. To obtain a good and reasonable
// starting value we solve an elastic problem in very first step (j=1).
  template <int dim>
  void
  PlasticityContactProblem<dim>::solve_newton ()
  {
    TimerOutput::Scope t(computing_timer, "solve newton setup");

    double resid = 0;
    double resid_old = 100000;
    TrilinosWrappers::MPI::Vector old_solution(system_rhs_newton);
    TrilinosWrappers::MPI::Vector res(system_rhs_newton);
    TrilinosWrappers::MPI::Vector tmp_vector(system_rhs_newton);

    double sigma_hlp = sigma_0;

    IndexSet active_set_old(active_set);

    t.stop(); // stop newton setup timer

    unsigned int j = 1;
    unsigned int number_assemble_system = 0;
    for (; j <= 100; j++)
      {
        if (transfer_solution)
          {
            if (transfer_solution && j == 1 && current_refinement_cycle == 0)
              constitutive_law.set_sigma_0(1e+10);
            else if (transfer_solution && (j == 2 || current_refinement_cycle > 0))
              constitutive_law.set_sigma_0(sigma_hlp);
          }
        else
          {
            if (j == 1)
              constitutive_law.set_sigma_0(1e+10);
            else
              constitutive_law.set_sigma_0(sigma_hlp);
          }

        pcout << " " << std::endl;
        pcout << "   Newton iteration " << j << std::endl;
        pcout << "      Updating active set..." << std::endl;

        {
          TimerOutput::Scope t(computing_timer, "update active set");
          update_solution_and_constraints();
        }

        pcout << "      Assembling system... " << std::endl;
        system_matrix_newton = 0;
        system_rhs_newton = 0;
        assemble_nl_system(solution); //compute Newton-Matrix

        number_assemble_system += 1;

        pcout << "      Solving system... " << std::endl;
        solve();

        TrilinosWrappers::MPI::Vector distributed_solution(system_rhs_newton);
        distributed_solution = solution;

        // We handle a highly nonlinear problem so we have to damp
        // the Newtons method. We refer that we iterate the new solution
        // in each Newton step and not only the solution update.
        // Since the solution set is a convex set and not a space we
        // compute for the damping a linear combination of the
        // previous and the current solution to guarantee that the
        // damped solution is in our solution set again.
        // At most we apply 10 damping steps.
        bool damped = false;
        tmp_vector = old_solution;
        double a = 0;
        for (unsigned int i = 0; (i < 5) && (!damped); i++)
          {
            a = std::pow(0.5, static_cast<double>(i));
            old_solution = tmp_vector;
            old_solution.sadd(1 - a, a, distributed_solution);
            old_solution.compress(VectorOperation::add);

            TimerOutput::Scope t(computing_timer, "Residual and lambda");

            system_rhs_newton = 0;
            system_rhs_lambda = 0;

            solution = old_solution;
            compute_nonlinear_residual(solution);
            res = system_rhs_newton;

            const unsigned int start_res = (res.local_range().first),
                               end_res = (res.local_range().second);
            for (unsigned int n = start_res; n < end_res; ++n)
              if (all_constraints.is_inhomogeneously_constrained(n))
                res(n) = 0;

            res.compress(VectorOperation::insert);

            resid = res.l2_norm();

            if (resid < resid_old)
              damped = true;

            pcout << "      Residual of the non-contact part of the system: "
                  << resid << std::endl
                  << "         with a damping parameter alpha = " << a
                  << std::endl;

            // The previous iteration of step 0 is the solution of an elastic problem.
            // So a linear combination of a plastic and an elastic solution makes no sense
            // since the elastic solution is not in the convex set of the plastic solution.
            if (!transfer_solution && j == 2)
              break;
            if (transfer_solution && j == 2 && current_refinement_cycle == 0)
              break;
          }

        resid_old = resid;

        resid_vector = system_rhs_lambda;

        if (Utilities::MPI::sum((active_set == active_set_old) ? 0 : 1,
                                mpi_communicator) == 0)
          {
            pcout << "      Active set did not change!" << std::endl;
            if (output_dir.compare("its/") != 0 && resid < 1e-7)
              break;
            else if (output_dir.compare("its/") == 0 && resid < 1e-10)
              break;
          }
        active_set_old = active_set;
      }

    pcout << "" << std::endl << "      Number of assembled systems = "
          << number_assemble_system << std::endl;
  }

// @sect3{The <code>refine_grid</code> function}

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

        parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
          triangulation, estimated_error_per_cell, 0.3, 0.03);
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
        TrilinosWrappers::MPI::Vector distributed_solution(system_rhs_newton);
        distributed_solution = solution;
        solution_transfer.interpolate(distributed_solution);
        solution = distributed_solution;
        compute_nonlinear_residual(solution);
        resid_vector = system_rhs_lambda;
      }
  }

// @sect3{The <code>move_mesh</code> function}

  template <int dim>
  void
  PlasticityContactProblem<dim>::move_mesh (
    const TrilinosWrappers::MPI::Vector &_complete_displacement) const
  {
    std::vector<bool> vertex_touched(triangulation.n_vertices(), false);

    for (typename DoFHandler<dim>::active_cell_iterator cell =
           dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
      if (cell->is_locally_owned())
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell;
             ++v)
          {
            if (vertex_touched[cell->vertex_index(v)] == false)
              {
                vertex_touched[cell->vertex_index(v)] = true;

                Point<dim> vertex_displacement;
                for (unsigned int d = 0; d < dim; ++d)
                  {
                    if (_complete_displacement(cell->vertex_dof_index(v, d))
                        != 0)
                      vertex_displacement[d] = _complete_displacement(
                                                 cell->vertex_dof_index(v, d));
                  }

                cell->vertex(v) += vertex_displacement;
              }
          }
  }

// @sect4{PlasticityContactProblem::output_results}

  template <int dim>
  void
  PlasticityContactProblem<dim>::output_results (const std::string &filename_base)
  {
    TimerOutput::Scope t(computing_timer, "Graphical output");

    pcout << "      Writing graphical output... " << std::flush;

    move_mesh(solution);

    // Calculation of the contact forces
    TrilinosWrappers::MPI::Vector lambda(solution);
    TrilinosWrappers::MPI::Vector distributed_lambda(system_rhs_newton);
    const unsigned int start_res = (resid_vector.local_range().first),
                       end_res = (resid_vector.local_range().second);
    for (unsigned int n = start_res; n < end_res; ++n)
      if (all_constraints.is_inhomogeneously_constrained(n))
        distributed_lambda(n) = resid_vector(n) / diag_mass_matrix_vector(n);
    distributed_lambda.compress(VectorOperation::insert);
    constraints_hanging_nodes.distribute(distributed_lambda);
    lambda = distributed_lambda;
    TrilinosWrappers::MPI::Vector resid_vector_relevant(solution);
    TrilinosWrappers::MPI::Vector distributed_resid_vector(resid_vector);
    constraints_hanging_nodes.distribute(distributed_resid_vector);
    resid_vector_relevant = distributed_resid_vector;

    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);

    const std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation(
      dim, DataComponentInterpretation::component_is_part_of_vector);
    data_out.add_data_vector(solution,
                             std::vector < std::string > (dim, "Displacement"),
                             DataOut<dim>::type_dof_data, data_component_interpretation);
    data_out.add_data_vector(lambda,
                             std::vector < std::string > (dim, "ContactForce"),
                             DataOut<dim>::type_dof_data, data_component_interpretation);
    data_out.add_data_vector(active_set,
                             std::vector < std::string > (dim, "ActiveSet"),
                             DataOut<dim>::type_dof_data, data_component_interpretation);
    data_out.add_data_vector(resid_vector_relevant,
                             std::vector < std::string > (dim, "Residual"),
                             DataOut<dim>::type_dof_data, data_component_interpretation);

    Vector<float> subdomain(triangulation.n_active_cells());
    for (unsigned int i = 0; i < subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector(subdomain, "subdomain");

    data_out.add_data_vector(fraction_of_plastic_q_points_per_cell, "FractionOfPlasticQPoints");

    data_out.build_patches();

    const std::string filename =
      (output_dir + filename_base + "-"
       + Utilities::int_to_string(triangulation.locally_owned_subdomain(), 4));

    std::ofstream output_vtu((filename + ".vtu").c_str());
    data_out.write_vtu(output_vtu);
    pcout << output_dir + filename_base << ".pvtu" << std::endl;

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        std::vector < std::string > filenames;
        for (unsigned int i = 0;
             i < Utilities::MPI::n_mpi_processes(mpi_communicator); ++i)
          filenames.push_back(
            filename_base + "-" + Utilities::int_to_string(i, 4) + ".vtu");

        std::ofstream master_output((output_dir + filename_base + ".pvtu").c_str());
        data_out.write_pvtu_record(master_output, filenames);
      }

    TrilinosWrappers::MPI::Vector tmp(solution);
    tmp *= -1;
    move_mesh(tmp);
  }

// @sect4{PlasticityContactProblem::output_contact_force}

// This function provides the contact force by calculating
// an integral over the contact pressure in z-directions
// over the contact area. For this purpose we set the contact
// pressure lambda to 0 for all inactive dofs. For all
// active dofs we lambda contains the quotient of the nonlinear
// residual (resid_vector) and corresponding diagonal entry
// of the mass matrix (diag_mass_matrix_vector). Because it is
// not unlikely that hanging nodes shows up in the contact area
// it is important to apply contraints_hanging_nodes.distribute
// to the distributed_lambda vector.
// To calculate the contact pressure in a certain point in the
// contact area, we apply the Functions::FEFieldFunction.
// In parallel this is a little tricky because we have to find the
// process with the right cell which contains this point. If
// a processor does not own the cell with the point we have to
// catch these cases.
  template <int dim>
  void
  PlasticityContactProblem<dim>::output_contact_force () const
  {
    Functions::FEFieldFunction<dim, DoFHandler<dim>,
              TrilinosWrappers::MPI::Vector> solution_function(dof_handler,
                                                               solution);
    std::cout.precision(10);

    Vector<double> solution_p1(dim);
    std::vector<Tensor<1, dim> > solution_gradient_p1(dim);

    // Here we calculate the contact pressure as a vector lambda.
    // If a dof is element of the active set lambda contains the
    // nonlinear residual this dof divided by the according entry
    // of the mass matrix. In all other dofs lambda will be set to
    // zero.
    TrilinosWrappers::MPI::Vector lambda(solution);
    TrilinosWrappers::MPI::Vector distributed_lambda(system_rhs_newton);
    const unsigned int start_res = (resid_vector.local_range().first),
                       end_res = (resid_vector.local_range().second);
    for (unsigned int n = start_res; n < end_res; ++n)
      if (all_constraints.is_inhomogeneously_constrained(n))
        distributed_lambda(n) = resid_vector(n) / diag_mass_matrix_vector(n);
      else
        distributed_lambda(n) = 0;
    distributed_lambda.compress(VectorOperation::insert);
    constraints_hanging_nodes.distribute(distributed_lambda);
    lambda = distributed_lambda;
    Functions::FEFieldFunction<dim, DoFHandler<dim>,
              TrilinosWrappers::MPI::Vector> lambda_function(dof_handler, lambda);

    // Here we try to find the MPI-process which owns the cell
    // with the point_of_interest. If it is the wrong MPI-process
    // we catch this case and set point_found to false.
    const Point<dim> point_of_interest(0.49, 0.5001, 1.0);
    Vector<double> contact_pressure_in_point(dim);
    bool point_found = true;

    MPI_Barrier(MPI_COMM_WORLD);
    try
      {
        lambda_function.vector_value(point_of_interest,
                                     contact_pressure_in_point);
      }
    catch (const typename Functions::FEFieldFunction<dim, DoFHandler<dim>,
             TrilinosWrappers::MPI::Vector>::ExcPointNotAvailableHere &)
      {
        point_found = false;
      }

    if (point_found == true)
      {
        std::cout << "PoI contact pressure: " << contact_pressure_in_point(2)
                  << std::endl;
      }

    // To obtain the contact force we have to compute an integral of the contact pressure
    // in z-direction over the whole contact area. To be accurate enough we use the
    // Gaussian quadrature rule with fe.degree + 1.
    double contact_force = 0.0;
    {
      QGauss<dim - 1> face_quadrature_formula(fe.degree + 1);

      FEFaceValues<dim> fe_values_face(fe, face_quadrature_formula,
                                       update_values | update_quadrature_points | update_JxW_values);

      const unsigned int n_face_q_points = face_quadrature_formula.size();

      const FEValuesExtractors::Vector displacement(0);

      typename DoFHandler<dim>::active_cell_iterator cell =
        dof_handler.begin_active(), endc = dof_handler.end();
      for (; cell != endc; ++cell)
        if (cell->is_locally_owned())
          for (unsigned int face = 0;
               face < GeometryInfo<dim>::faces_per_cell; ++face)
            if (cell->face(face)->at_boundary()
                && cell->face(face)->boundary_indicator() == 1)
              {
                fe_values_face.reinit(cell, face);

                std::vector<Tensor<1, dim> > lambda_values(n_face_q_points);
                fe_values_face[displacement].get_function_values(lambda,
                                                                 lambda_values);

                for (unsigned int q_point = 0; q_point < n_face_q_points;
                     ++q_point)
                  {
                    contact_force += lambda_values[q_point][2]
                                     * fe_values_face.JxW(q_point);
                  }
              }
      contact_force = Utilities::MPI::sum(contact_force, MPI_COMM_WORLD);
      pcout << "Contact force = " << contact_force << std::endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
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
