/* $Id$ */
/* Author: Wolfgang Bangerth, University of Heidelberg, 1999 */

/*    $Id$       */
/*                                                                */
/*    Copyright (C) 1999-2012 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyrightG and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

// @sect3{Include files}

// We are using the the same
// include files as in step-41:

#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/base/timer.h>
#include <fstream>
#include <iostream>
#include <list>
#include <time.h>


#include <deal.II/base/logstream.h>

namespace Step42
{
  using namespace dealii;

  // @sect3{The <code>PlasticityContactProblem</code> class template}

  template <int dim> class Input;

  // This class provides an interface
  // for a constitutive law. In this
  // example we are using an elastic
  // plastic material with linear,
  // isotropic hardening.

  template <int dim> class ConstitutiveLaw;

  // This class supplies all function
  // and variables needed to describe
  // the nonlinear contact problem. It is
  // close to step-41 but with some additional
  // features like: handling hanging nodes,
  // a newton method, using Trilinos and p4est
  // for parallel distributed computing.
  // To deal with hanging nodes makes
  // life a bit more complicated since
  // we need an other ConstraintMatrix now.
  // We create a newton method for the
  // active set method for the contact
  // situation and to handle the nonlinear
  // operator for the constitutive law.

  template <int dim>
  class PlasticityContactProblem
  {
  public:
    PlasticityContactProblem (int _n_refinements_global);
    void run ();

  private:
    void make_grid ();
    void setup_system();
    void assemble_nl_system (TrilinosWrappers::MPI::Vector &u);
    void residual_nl_system (TrilinosWrappers::MPI::Vector &u);
    void assemble_mass_matrix_diagonal (TrilinosWrappers::SparseMatrix &mass_matrix);
    void update_solution_and_constraints ();
    void dirichlet_constraints ();
    void solve ();
    void solve_newton ();
    void refine_grid ();
    void move_mesh (const TrilinosWrappers::MPI::Vector &_complete_displacement) const;
    void output_results (const std::string &title) const;

    int                  n_refinements_global;

    MPI_Comm             mpi_communicator;

    parallel::distributed::Triangulation<dim>   triangulation;

    FESystem<dim>        fe;
    DoFHandler<dim>      dof_handler;

    IndexSet             locally_owned_dofs;
    IndexSet             locally_relevant_dofs;

    unsigned int         number_iterations;
    std::vector<double>  run_time;

    ConstraintMatrix     constraints;
    ConstraintMatrix     constraints_hanging_nodes;
    ConstraintMatrix     constraints_dirichlet_hanging_nodes;

    TrilinosWrappers::SparseMatrix system_matrix_newton;

    TrilinosWrappers::MPI::Vector       solution;
    TrilinosWrappers::MPI::Vector       old_solution;
    TrilinosWrappers::MPI::Vector       system_rhs_newton;
    TrilinosWrappers::MPI::Vector       resid_vector;
    TrilinosWrappers::MPI::Vector       diag_mass_matrix_vector;
    IndexSet                            active_set;

    ConditionalOStream pcout;

    TrilinosWrappers::PreconditionAMG::AdditionalData additional_data;
    TrilinosWrappers::PreconditionAMG preconditioner_u;
    TrilinosWrappers::PreconditionAMG preconditioner_t;

    std::unique_ptr<Input<dim> >               input_obstacle;
    std::unique_ptr<ConstitutiveLaw<dim> >     plast_lin_hard;

    double sigma_0;    // Yield stress
    double gamma;      // Parameter for the linear isotropic hardening
    double e_modul;    // E-Modul
    double nu;         // Poisson ratio
  };

  template <int dim>
  class Input
  {
   public:
    Input (const char* _name) :
      name (_name),
      mpi_communicator (MPI_COMM_WORLD),
      pcout (std::cout,
               (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
      HV (NULL),
      lx (0),
      ly (0),
      nx (0),
      ny (0)
      {read_surface (name);}

      double hv(int i, int j);

      double& set_height(int i, int j);

      double mikro_height(double x,double y, double z);

      void read_surface(const char* name);

   private:
      const char*          name;
      MPI_Comm             mpi_communicator;
      ConditionalOStream   pcout;
      double*              HV;
      double               lx, ly;
      double               hx, hy;
      int                  nx, ny;
  };

  template <int dim>
  double Input<dim>::hv(int i, int j)
  {
    return HV[nx*j+i];  // i indiziert x-werte, j indiziert y-werte
  }

  template <int dim>
  double& Input<dim>::set_height(int i, int j)
  {
    return HV[nx*j+i];  // i indiziert x-werte, j indiziert y-werte
  }

  template <int dim>
  double Input<dim>::mikro_height(double x,double y, double z)
  {
    int ix = (int)(x/hx);
    int iy = (int)(y/hy);

    if (ix<0)
      ix = 0;

    if (iy<0)
      iy = 0;

    if (ix>=nx-1)
      ix = nx-2;

    if (iy>=ny-1)
      iy = ny-2;

    double val = 0.0;
    {
      FullMatrix<double> H(4,4);
      Vector<double>  X(4);
      Vector<double>  b(4);

      double xx = 0.0;
      double yy = 0.0;

      xx = ix*hx;
      yy = iy*hy;
      H(0,0) = xx;
      H(0,1) = yy;
      H(0,2) = xx*yy;
      H(0,3) = 1.0;
      b(0)   = hv (ix, iy);

      xx = (ix + 1)*hx;
      yy = iy*hy;
      H(1,0) = xx;
      H(1,1) = yy;
      H(1,2) = xx*yy;
      H(1,3) = 1.0;
      b(1)   = hv (ix + 1, iy);

      xx = (ix + 1)*hx;
      yy = (iy + 1)*hy;
      H(2,0) = xx;
      H(2,1) = yy;
      H(2,2) = xx*yy;
      H(2,3) = 1.0;
      b(2)   = hv (ix + 1, iy + 1);

      xx = ix*hx;
      yy = (iy + 1)*hy;
      H(3,0) = xx;
      H(3,1) = yy;
      H(3,2) = xx*yy;
      H(3,3) = 1.0;
      b(3)   = hv (ix, iy + 1);

      H.gauss_jordan ();
      H.vmult (X, b);

      val = X(0)*x + X(1)*y + X(2)*x*y + X(3);
    }

    return val;
  }

  template <int dim>
  void Input<dim>::read_surface(const char* name)
  {
    int SZ = 100000;
    FILE* fp = fopen (name, "r");
    char* zeile   = new char[SZ];
    char* hlp_str = new char[SZ];
    double hlp;
    int POS;

    fgets (zeile, SZ, fp);
    POS = strcspn (zeile, "=");
    for (int i=0; i<=POS; i++)
      zeile[i] = ' ';
    sscanf (zeile, "%d", &nx);

    fgets (zeile,SZ,fp);
    POS = strcspn (zeile, "=");
    for (int i=0; i<=POS; i++)
      zeile[i] = ' ';
    sscanf (zeile, "%d", &ny);

    fgets (zeile, SZ, fp);
    POS = strcspn (zeile, "=");
    for (int i=0; i<=POS; i++)
      zeile[i] = ' ';
    sscanf (zeile, "%lf", &lx);

    fgets (zeile, SZ, fp);
    POS = strcspn(zeile,"=");
    for (int i=0; i<=POS; i++)
      zeile[i] = ' ';
    sscanf(zeile,"%lf",&ly);

    hx = lx/(nx - 1);
    hy = ly/(ny - 1);

    pcout<< "Resolution of the scanned obstacle picture: " << nx << " x " << ny <<std::endl;

    if (HV) delete[] HV;
    HV = new double [nx*ny];

    int j = 0;
    double max_hlp = 0;
    double min_hlp = 1e+10;
    while (fgets (zeile, SZ, fp))
      {
        int reached = 0;
        for (int k=0; !reached; k++)
          {
            sscanf (zeile, "%lf", &hlp);

            if (hlp > max_hlp)
              max_hlp=hlp;
            if (hlp < min_hlp)
              min_hlp=hlp;

            set_height (k, ny - 1 - j) = hlp;
            int pos = strcspn (zeile, ",");
            if (!strpbrk (zeile, ","))
              {
                reached = 1;
                continue;
              }

            for (int i=0; i<=pos; i++)
              {
                zeile[i] = ' ';
              }
          }
        j++;
      }
    pcout<< "Highest point of the obstacle: " << max_hlp <<std::endl;
    pcout<< "Lowest point of the obstacle:  " << min_hlp <<std::endl;
  }

  template <int dim>
  class ConstitutiveLaw
  {
  public:
    ConstitutiveLaw (double _E,
                     double _nu,
                     double _sigma_0,
                     double _gamma,
                     MPI_Comm _mpi_communicator,
                     ConditionalOStream _pcout);

    void plast_linear_hardening (SymmetricTensor<4,dim>  &stress_strain_tensor,
                                 SymmetricTensor<2,dim>  &strain_tensor,
                                 unsigned int            &elast_points,
                                 unsigned int            &plast_points,
                                 double                  &yield);
    void linearized_plast_linear_hardening (SymmetricTensor<4,dim>  &stress_strain_tensor_linearized,
                                            SymmetricTensor<4,dim>  &stress_strain_tensor,
                                            SymmetricTensor<2,dim>  &strain_tensor);
    inline SymmetricTensor<2,dim> get_strain (const FEValues<dim> &fe_values,
                                              const unsigned int  shape_func,
                                              const unsigned int  q_point) const;

  private:
    SymmetricTensor<4,dim>  stress_strain_tensor_mu;
    SymmetricTensor<4,dim>  stress_strain_tensor_kappa;
    double E;
    double nu;
    double sigma_0;
    double gamma;
    double mu;
    double kappa;
    MPI_Comm mpi_communicator;
    ConditionalOStream pcout;
  };

  template <int dim>
  ConstitutiveLaw<dim>::ConstitutiveLaw(double _E, double _nu, double _sigma_0, double _gamma, MPI_Comm _mpi_communicator, ConditionalOStream _pcout)
    :E (_E),
     nu (_nu),
     sigma_0 (_sigma_0),
     gamma (_gamma),
     mpi_communicator (_mpi_communicator),
     pcout (_pcout)
  {
    mu = E/(2*(1+nu));
    kappa = E/(3*(1-2*nu));
    pcout<< "-----> mu = " << mu << ", kappa = " << kappa <<std::endl;
    stress_strain_tensor_kappa = kappa*outer_product(unit_symmetric_tensor<dim>(), unit_symmetric_tensor<dim>());
    stress_strain_tensor_mu = 2*mu*(identity_tensor<dim>() - outer_product(unit_symmetric_tensor<dim>(), unit_symmetric_tensor<dim>())/3.0);
  }

  template <int dim>
  inline
  SymmetricTensor<2,dim> ConstitutiveLaw<dim>::get_strain (const FEValues<dim> &fe_values,
                                                           const unsigned int   shape_func,
                                                           const unsigned int   q_point) const
  {
    const FEValuesExtractors::Vector displacement (0);
    SymmetricTensor<2,dim> tmp;

    tmp = fe_values[displacement].symmetric_gradient (shape_func,q_point);

    return tmp;
  }

  template <int dim>
  void ConstitutiveLaw<dim>::plast_linear_hardening (SymmetricTensor<4,dim>  &stress_strain_tensor,
                                                     SymmetricTensor<2,dim>  &strain_tensor,
                                                     unsigned int            &elast_points,
                                                     unsigned int            &plast_points,
                                                     double                  &yield)
  {
    if (dim == 3)
      {
        SymmetricTensor<2,dim> stress_tensor;
        stress_tensor = (stress_strain_tensor_kappa + stress_strain_tensor_mu)*strain_tensor;
        double tmp = E/((1+nu)*(1-2*nu));

        SymmetricTensor<2,dim> deviator_stress_tensor = deviator(stress_tensor);

        double deviator_stress_tensor_norm = deviator_stress_tensor.norm ();

        yield = 0;
        stress_strain_tensor = stress_strain_tensor_mu;
        double beta = 1.0;
        if (deviator_stress_tensor_norm >= sigma_0)
          {
            beta = (sigma_0 + gamma)/deviator_stress_tensor_norm;
            stress_strain_tensor *= beta;
            yield = 1;
            plast_points += 1;
          }
        else
          elast_points += 1;

        stress_strain_tensor += stress_strain_tensor_kappa;
      }
  }

  template <int dim>
  void ConstitutiveLaw<dim>::linearized_plast_linear_hardening (SymmetricTensor<4,dim>  &stress_strain_tensor_linearized,
      SymmetricTensor<4,dim>  &stress_strain_tensor,
      SymmetricTensor<2,dim>  &strain_tensor)
  {
    if (dim == 3)
      {
        SymmetricTensor<2,dim> stress_tensor;
        stress_tensor = (stress_strain_tensor_kappa + stress_strain_tensor_mu)*strain_tensor;
        double tmp = E/((1+nu)*(1-2*nu));

        stress_strain_tensor = stress_strain_tensor_mu;
        stress_strain_tensor_linearized = stress_strain_tensor_mu;

        SymmetricTensor<2,dim> deviator_stress_tensor = deviator(stress_tensor);

        double deviator_stress_tensor_norm = deviator_stress_tensor.norm ();

        double beta = 1.0;
        if (deviator_stress_tensor_norm >= sigma_0)
          {
            beta = (sigma_0 + gamma)/deviator_stress_tensor_norm;
            stress_strain_tensor *= beta;
            stress_strain_tensor_linearized *= beta;
            deviator_stress_tensor /= deviator_stress_tensor_norm;
            stress_strain_tensor_linearized -= beta*2*mu*outer_product(deviator_stress_tensor, deviator_stress_tensor);
          }

        stress_strain_tensor += stress_strain_tensor_kappa;
        stress_strain_tensor_linearized += stress_strain_tensor_kappa;
      }
  }

  namespace EquationData
  {
    template <int dim>
    class RightHandSide : public Function<dim>
    {
    public:
      RightHandSide () : Function<dim>(dim) {}

      virtual double value (const Point<dim>   &p,
                            const unsigned int  component = 0) const;

      virtual void vector_value (const Point<dim> &p,
                                 Vector<double>   &values) const;
    };

    template <int dim>
    double RightHandSide<dim>::value (const Point<dim> &p,
                                      const unsigned int component) const
    {
      double return_value = 0.0;

      if (component == 0)
        return_value = 0.0;
      if (component == 1)
        return_value = 0.0;
      if (component == 2)
        // if ((p(0)-0.5)*(p(0)-0.5)+(p(1)-0.5)*(p(1)-0.5) < 0.2)
        //  return_value = -5000;
        // else
        return_value = 0.0;
      // for (unsigned int i=0; i<dim; ++i)
      //   return_value += 4*std::pow(p(i), 4);

      return return_value;
    }

    template <int dim>
    void RightHandSide<dim>::vector_value (const Point<dim> &p,
                                           Vector<double>   &values) const
    {
      for (unsigned int c=0; c<this->n_components; ++c)
        values(c) = RightHandSide<dim>::value (p, c);
    }


    template <int dim>
    class BoundaryValues : public Function<dim>
    {
    public:
      BoundaryValues () : Function<dim>(dim) {};

      virtual double value (const Point<dim>   &p,
                            const unsigned int  component = 0) const;

      virtual void vector_value (const Point<dim> &p,
                                 Vector<double>   &values) const;
    };

    template <int dim>
    double BoundaryValues<dim>::value (const Point<dim> &p,
                                       const unsigned int component) const
    {
      double return_value = 0;

      if (component == 0)
        return_value = 0.0;
      if (component == 1)
        return_value = 0.0;
      if (component == 2)
        return_value = 0.0;

      return return_value;
    }

    template <int dim>
    void BoundaryValues<dim>::vector_value (const Point<dim> &p,
                                            Vector<double>   &values) const
    {
      for (unsigned int c=0; c<this->n_components; ++c)
        values(c) = BoundaryValues<dim>::value (p, c);
    }


    template <int dim>
    class Obstacle : public Function<dim>
    {
    public:
      Obstacle (std::unique_ptr<Input<dim> > const &_input) :
        Function<dim>(dim),
        input_obstacle_copy(std::move (_input)) {};

      virtual double value (const Point<dim>   &p,
                            const unsigned int  component = 0) const;

      virtual void vector_value (const Point<dim> &p,
                                 Vector<double>   &values) const;

    private:
      std::unique_ptr<Input<dim> >  const &input_obstacle_copy;
    };

    template <int dim>
    double Obstacle<dim>::value (const Point<dim> &p,
                                 const unsigned int component) const
    {
      double R = 0.03;
      double return_value = 100.0;
      if (component == 0)
        return_value = p(0);
      if (component == 1)
        return_value = p(1);
      if (component == 2)
        {
          // Hindernis Dortmund
//          double x1 = p(0);
//          double x2 = p(1);
//          if (((x2-0.5)*(x2-0.5)+(x1-0.5)*(x1-0.5)<=0.3*0.3)&&((x2-0.5)*(x2-0.5)+(x1-1.0)*(x1-1.0)>=0.4*0.4)&&((x2-0.5)*(x2-0.5)+x1*x1>=0.4*0.4))
//            return_value = 0.999;
//          else
//            return_value = 1e+10;

          // Hindernis Werkzeug TKSE
           return_value = 1.999 - input_obstacle_copy->mikro_height (p(0), p(1), p(2));
//           std::cout<< "Obstacle value: " << return_value
//               << " p(0) = " << p(0)
//               << " p(1) = " << p(1)
//               <<std::endl;

          // Ball with radius R
          // double R = 1.0;
          // if (std::pow ((p(0)-1.0/2.0), 2) + std::pow ((p(1)-1.0/2.0), 2) < R*R)
          //   return_value = 1.0 + R - 0.01 - sqrt (R*R  - std::pow ((p(0)-1.0/2.0), 2)
          //           - std::pow ((p(1)-1.0/2.0), 2));
          // else
          //   return_value = 1e+5;
        }
      return return_value;
    }

    template <int dim>
    void Obstacle<dim>::vector_value (const Point<dim> &p,
                                      Vector<double> &values) const
    {
      for (unsigned int c=0; c<this->n_components; ++c)
        values(c) = Obstacle<dim>::value (p, c);
    }
  }


  // @sect3{Implementation of the <code>PlasticityContactProblem</code> class}

  // Next for the implementation of the class
  // template that makes use of the functions
  // above. As before, we will write everything

  template <int dim>
  PlasticityContactProblem<dim>::PlasticityContactProblem (int _n_refinements_global)
    :
    n_refinements_global (_n_refinements_global),
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator),
    fe (FE_Q<dim>(1), dim),
    dof_handler (triangulation),
    pcout (std::cout,
           (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
    sigma_0 (400),
    gamma (1.e-2),
    e_modul (2.0e5),
    nu (0.3)
  {
    // double _E, double _nu, double _sigma_0, double _gamma
    plast_lin_hard.reset (new ConstitutiveLaw<dim> (e_modul, nu, sigma_0, gamma, mpi_communicator, pcout));
  }

  template <int dim>
  void PlasticityContactProblem<dim>::make_grid ()
  {
    std::vector<unsigned int> repet(3);
    repet[0] = 1;//20;
    repet[1] = 1;
    repet[2] = 1;

    Point<dim> p1 (0,0,0);
    Point<dim> p2 (1.0, 1.0, 1.0);
    GridGenerator::subdivided_hyper_rectangle (triangulation, repet, p1, p2);

    Triangulation<3>::active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();

    /* boundary_indicators:
              _______
             /  9    /|
            /______ / |
          8|       | 8|
           |   8   | /
           |_______|/
               6
     */

    for (; cell!=endc; ++cell)
      for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
        {
          if (cell->face (face)->center ()[2] == p2(2))
            cell->face (face)->set_boundary_indicator (9);
          if (cell->face (face)->center ()[0] == p1(0) ||
              cell->face (face)->center ()[0] == p2(0) ||
              cell->face (face)->center ()[1] == p1(1) ||
              cell->face (face)->center ()[1] == p2(1))
            cell->face (face)->set_boundary_indicator (8);
          if (cell->face (face)->center ()[2] == p1(2))
            cell->face (face)->set_boundary_indicator (6);
        }

    triangulation.refine_global (n_refinements_global);
  }

  template <int dim>
  void PlasticityContactProblem<dim>::setup_system ()
  {
    // setup dofs
    {
      dof_handler.distribute_dofs (fe);

      locally_owned_dofs = dof_handler.locally_owned_dofs ();
      locally_relevant_dofs.clear();
      DoFTools::extract_locally_relevant_dofs (dof_handler,
                                               locally_relevant_dofs);
    }

    // setup hanging nodes and dirichlet constraints
    {
      constraints_hanging_nodes.clear ();
      constraints_hanging_nodes.reinit (locally_relevant_dofs);
      DoFTools::make_hanging_node_constraints (dof_handler,
                                               constraints_hanging_nodes);
      constraints_hanging_nodes.close ();

      pcout << "   Number of active cells: "
            << triangulation.n_active_cells()
            << std::endl
            << "   Total number of cells: "
            << triangulation.n_cells()
            << std::endl
            << "   Number of degrees of freedom: "
            << dof_handler.n_dofs ()
            << std::endl;

      dirichlet_constraints ();
    }

    // Initialzation for matrices and vectors
    {
      solution.reinit (locally_relevant_dofs, mpi_communicator);
      system_rhs_newton.reinit (locally_owned_dofs, mpi_communicator);
      old_solution.reinit (system_rhs_newton);
      resid_vector.reinit (system_rhs_newton);
      diag_mass_matrix_vector.reinit (system_rhs_newton);
      active_set.clear ();
      active_set.set_size (locally_relevant_dofs.size ());
    }

    // setup sparsity pattern
    {
      TrilinosWrappers::SparsityPattern sp (locally_owned_dofs,
                                            mpi_communicator);

      DoFTools::make_sparsity_pattern (dof_handler, sp, constraints_dirichlet_hanging_nodes, false,
                                       Utilities::MPI::this_mpi_process(mpi_communicator));

      sp.compress();

      system_matrix_newton.reinit (sp);

      TrilinosWrappers::SparseMatrix mass_matrix;
      mass_matrix.reinit (sp);
      assemble_mass_matrix_diagonal (mass_matrix);
      const unsigned int
      start = (system_rhs_newton.local_range().first),
      end   = (system_rhs_newton.local_range().second);
      for (unsigned int j=start; j<end; j++)
        diag_mass_matrix_vector (j) = mass_matrix.diag_element (j);
      number_iterations = 0;

      diag_mass_matrix_vector.compress (VectorOperation::insert);
    }
  }

  template <int dim>
  void PlasticityContactProblem<dim>::assemble_nl_system (TrilinosWrappers::MPI::Vector &u)
  {
    QGauss<dim>  quadrature_formula(2);
    QGauss<dim-1>  face_quadrature_formula(2);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             UpdateFlags(update_values    |
                                         update_gradients |
                                         update_q_points  |
                                         update_JxW_values));

    FEFaceValues<dim> fe_values_face (fe, face_quadrature_formula,
                                      update_values   | update_quadrature_points |
                                      update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size ();
    const unsigned int   n_face_q_points = face_quadrature_formula.size();

    const EquationData::RightHandSide<dim> right_hand_side;
    std::vector<Vector<double> > right_hand_side_values (n_q_points,
                                                         Vector<double>(dim));
    std::vector<Vector<double> > right_hand_side_values_face (n_face_q_points,
                                                              Vector<double>(dim));

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       cell_rhs (dofs_per_cell);

    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();

    const FEValuesExtractors::Vector displacement (0);

    TrilinosWrappers::MPI::Vector   test_rhs(solution);
    const double kappa = 1.0;
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          cell_matrix = 0;
          cell_rhs = 0;

          right_hand_side.vector_value_list (fe_values.get_quadrature_points(),
                                             right_hand_side_values);

          std::vector<SymmetricTensor<2,dim> > strain_tensor (n_q_points);
          fe_values[displacement].get_function_symmetric_gradients (u, strain_tensor);

          for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
            {
              SymmetricTensor<4,dim> stress_strain_tensor_linearized;
              SymmetricTensor<4,dim> stress_strain_tensor;
              SymmetricTensor<2,dim> stress_tensor;

              plast_lin_hard->linearized_plast_linear_hardening (stress_strain_tensor_linearized,
                                                                 stress_strain_tensor,
                                                                 strain_tensor[q_point]);

              for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                  stress_tensor = stress_strain_tensor_linearized * plast_lin_hard->get_strain(fe_values, i, q_point);

                  for (unsigned int j=0; j<dofs_per_cell; ++j)
                    {
                      cell_matrix(i,j) += (stress_tensor *
                                           plast_lin_hard->get_strain(fe_values, j, q_point) *
                                           fe_values.JxW (q_point));
                    }

                  // the linearized part a(v^i;v^i,v) of the rhs
                  cell_rhs(i) += (stress_tensor *
                                  strain_tensor[q_point] *
                                  fe_values.JxW (q_point));

                  // the residual part a(v^i;v) of the rhs
                  cell_rhs(i) -= (strain_tensor[q_point] * stress_strain_tensor *
                                  plast_lin_hard->get_strain(fe_values, i, q_point) *
                                  fe_values.JxW (q_point));

                  // the residual part F(v) of the rhs
                  Tensor<1,dim> rhs_values;
                  rhs_values = 0;
                  cell_rhs(i) += (fe_values[displacement].value (i, q_point) *
                                  rhs_values *
                                  fe_values.JxW (q_point));
                }
            }

          for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
            {
              if (cell->face (face)->at_boundary()
                  && cell->face (face)->boundary_indicator () == 9)
                {
                  fe_values_face.reinit (cell, face);

                  right_hand_side.vector_value_list (fe_values_face.get_quadrature_points(),
                                                     right_hand_side_values_face);

                  for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
                    {
                      Tensor<1,dim> rhs_values;
                      rhs_values = 0;
                      for (unsigned int i=0; i<dofs_per_cell; ++i)
                        cell_rhs(i) += (fe_values_face[displacement].value (i, q_point) *
                                        rhs_values *
                                        fe_values_face.JxW (q_point));
                    }
                }
            }

          cell->get_dof_indices (local_dof_indices);
          constraints.distribute_local_to_global (cell_matrix, cell_rhs,
                                                  local_dof_indices,
                                                  system_matrix_newton, system_rhs_newton, true);
        };

    system_matrix_newton.compress (VectorOperation::add);
    system_rhs_newton.compress (VectorOperation::add);
  }

  template <int dim>
  void PlasticityContactProblem<dim>::residual_nl_system (TrilinosWrappers::MPI::Vector &u)
  {
    QGauss<dim>  quadrature_formula(2);
    QGauss<dim-1> face_quadrature_formula(2);

    FEValues<dim> fe_values (fe, quadrature_formula,
                             UpdateFlags(update_values    |
                                         update_gradients |
                                         update_q_points  |
                                         update_JxW_values));

    FEFaceValues<dim> fe_values_face (fe, face_quadrature_formula,
                                      update_values   | update_quadrature_points |
                                      update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size ();
    const unsigned int   n_face_q_points = face_quadrature_formula.size();

    const EquationData::RightHandSide<dim> right_hand_side;
    std::vector<Vector<double> > right_hand_side_values (n_q_points,
                                                         Vector<double>(dim));
    std::vector<Vector<double> > right_hand_side_values_face (n_face_q_points,
                                                              Vector<double>(dim));

    Vector<double>       cell_rhs (dofs_per_cell);

    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    const FEValuesExtractors::Vector displacement (0);

    typename DoFHandler<dim>::active_cell_iterator cell = dof_handler.begin_active(),
                                                   endc = dof_handler.end();

    unsigned int elast_points = 0;
    unsigned int plast_points = 0;
    double       yield = 0;
    unsigned int cell_number = 0;
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          cell_rhs = 0;

          right_hand_side.vector_value_list (fe_values.get_quadrature_points(),
                                             right_hand_side_values);

          std::vector<SymmetricTensor<2,dim> > strain_tensor (n_q_points);
          fe_values[displacement].get_function_symmetric_gradients (u, strain_tensor);

          for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
            {
              SymmetricTensor<4,dim> stress_strain_tensor;
              SymmetricTensor<2,dim> stress_tensor;

              plast_lin_hard->plast_linear_hardening (stress_strain_tensor, strain_tensor[q_point],
                                                      elast_points, plast_points, yield);

              for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                  cell_rhs(i) -= (strain_tensor[q_point] * stress_strain_tensor * //(stress_tensor) *
                                  plast_lin_hard->get_strain(fe_values, i, q_point) *
                                  fe_values.JxW (q_point));

                  Tensor<1,dim> rhs_values;
                  rhs_values = 0;
                  cell_rhs(i) += ((fe_values[displacement].value (i, q_point) *
                                   rhs_values) *
                                  fe_values.JxW (q_point));
                };
            };

          for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
            {
              if (cell->face (face)->at_boundary()
                  && cell->face (face)->boundary_indicator () == 9)
                {
                  fe_values_face.reinit (cell, face);

                  right_hand_side.vector_value_list (fe_values_face.get_quadrature_points(),
                                                     right_hand_side_values_face);

                  for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
                    {
                      Tensor<1,dim> rhs_values;
                      rhs_values = 0;
                      for (unsigned int i=0; i<dofs_per_cell; ++i)
                        cell_rhs(i) += (fe_values_face[displacement].value (i, q_point) *
                                        rhs_values *
                                        fe_values_face.JxW (q_point));
                    }
                }
            }

          cell->get_dof_indices (local_dof_indices);
          constraints_dirichlet_hanging_nodes.distribute_local_to_global (cell_rhs,
              local_dof_indices,
              system_rhs_newton);

          cell_number += 1;
        };

    system_rhs_newton.compress (VectorOperation::add);

    unsigned int sum_elast_points = Utilities::MPI::sum(elast_points, mpi_communicator);
    unsigned int sum_plast_points = Utilities::MPI::sum(plast_points, mpi_communicator);
    pcout << "      Number of elastic quadrature points: " << sum_elast_points
          << " and plastic quadrature points: " << sum_plast_points <<std::endl;
  }

  template <int dim>
  void PlasticityContactProblem<dim>::assemble_mass_matrix_diagonal (TrilinosWrappers::SparseMatrix &mass_matrix)
  {
    QTrapez<dim-1>  face_quadrature_formula;

    FEFaceValues<dim> fe_values_face (fe, face_quadrature_formula,
                                      update_values   |
                                      update_quadrature_points |
                                      update_JxW_values);

    const unsigned int   dofs_per_cell      = fe.dofs_per_cell;
    const unsigned int   n_face_q_points    = face_quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    const FEValuesExtractors::Vector displacement (0);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
          if (cell->face (face)->at_boundary()
              && cell->face (face)->boundary_indicator () == 9)
            {
              fe_values_face.reinit (cell, face);
              cell_matrix = 0;

              for (unsigned int q_point=0; q_point<n_face_q_points; ++q_point)
                for (unsigned int i=0; i<dofs_per_cell; ++i)
                  cell_matrix(i,i) += (fe_values_face[displacement].value (i, q_point) *
                                       fe_values_face[displacement].value (i, q_point) *
                                       fe_values_face.JxW (q_point));

              cell->get_dof_indices (local_dof_indices);

              constraints_dirichlet_hanging_nodes.distribute_local_to_global (cell_matrix,
                  local_dof_indices,
                  mass_matrix);
            }

    mass_matrix.compress (VectorOperation::add);
  }

  // @sect4{PlasticityContactProblem::update_solution_and_constraints}

  // Projection and updating of the active set
  // for the dofs which penetrates the obstacle.
  template <int dim>
  void PlasticityContactProblem<dim>::update_solution_and_constraints ()
  {
    clock_t                        start_proj, end_proj;

    const EquationData::Obstacle<dim>     obstacle (input_obstacle);
    std::vector<bool>                     vertex_touched (dof_handler.n_dofs (), false);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    TrilinosWrappers::MPI::Vector         distributed_solution (system_rhs_newton);
    distributed_solution = solution;
    TrilinosWrappers::MPI::Vector         lambda (solution);
    lambda = resid_vector;
    TrilinosWrappers::MPI::Vector         diag_mass_matrix_vector_relevant (solution);
    diag_mass_matrix_vector_relevant = diag_mass_matrix_vector;

    constraints.reinit(locally_relevant_dofs);
    active_set.clear ();
    IndexSet     active_set_locally_owned;
    active_set_locally_owned.set_size (locally_owned_dofs.size ());
    const double c = 100.0*e_modul;

    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
          if (cell->face (face)->at_boundary()
              && cell->face (face)->boundary_indicator () == 9)
            for (unsigned int v=0; v<GeometryInfo<dim-1>::vertices_per_cell; ++v)
              {
                unsigned int index_z = cell->face (face)->vertex_dof_index (v,2);

                if (vertex_touched[cell->face (face)->vertex_index(v)] == false)
                  vertex_touched[cell->face (face)->vertex_index(v)] = true;
                else
                  continue;

                // the local row where
                Point<dim> point (cell->face (face)->vertex (v)[0],
                                  cell->face (face)->vertex (v)[1],
                                  cell->face (face)->vertex (v)[2]);

                double obstacle_value = obstacle.value (point, 2);
                double solution_index_z = solution (index_z);
                double gap = obstacle_value - point (2);

                if (lambda (index_z) +
                    c *
                    diag_mass_matrix_vector_relevant (index_z) *
                    (solution_index_z - gap)
                    > 0 &&
                    !(constraints_hanging_nodes.is_constrained(index_z)))
                  {
                    constraints.add_line (index_z);
                    constraints.set_inhomogeneity (index_z, gap);

                    distributed_solution (index_z) = gap;

                    if (locally_relevant_dofs.is_element (index_z))
                      active_set.add_index (index_z);

                    if (locally_owned_dofs.is_element (index_z))
                      active_set_locally_owned.add_index (index_z);
                  }
              }
    distributed_solution.compress (VectorOperation::insert);

    unsigned int sum_contact_constraints = Utilities::MPI::sum(active_set_locally_owned.n_elements (),
                                                               mpi_communicator);
    pcout << "         Size of active set: " << sum_contact_constraints <<std::endl;

    solution = distributed_solution;

    constraints.close ();

    constraints.merge (constraints_dirichlet_hanging_nodes);
  }

  // @sect4{PlasticityContactProblem::dirichlet_constraints}

  // This function defines the new ConstraintMatrix
  // constraints_dirichlet_hanging_nodes. It contains
  // the dirichlet boundary values as well as the
  // hanging nodes constraints.
  template <int dim>
  void PlasticityContactProblem<dim>::dirichlet_constraints ()
  {
    /* boundary_indicators:
              _______
             /  9    /|
            /______ / |
          8|       | 8|
           |   8   | /
           |_______|/
               6
     */

    constraints_dirichlet_hanging_nodes.reinit (locally_relevant_dofs);
    constraints_dirichlet_hanging_nodes.merge (constraints_hanging_nodes);

    // interpolate all components of the solution
    VectorTools::interpolate_boundary_values (dof_handler,
                                              6,
                                              EquationData::BoundaryValues<dim>(),
                                              constraints_dirichlet_hanging_nodes,
                                              ComponentMask());

    // interpolate x- and y-components of the
    // solution (this is a bit mask, so apply
    // operator| )
    FEValuesExtractors::Scalar x_displacement(0);
    FEValuesExtractors::Scalar y_displacement(1);
    VectorTools::interpolate_boundary_values (dof_handler,
                                              8,
                                              EquationData::BoundaryValues<dim>(),
                                              constraints_dirichlet_hanging_nodes,
                                              (fe.component_mask(x_displacement)
                                               |
                                               fe.component_mask(y_displacement)));
    constraints_dirichlet_hanging_nodes.close ();
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

  // The rest of the funtion is smiliar to
  // step-41 except that we use a FGMRES-solver
  // instead of CG. For a very small hardening
  // value gamma the linear system becomes
  // almost semi definite but still symmetric.
  template <int dim>
  void PlasticityContactProblem<dim>::solve ()
  {
    Timer t;

    TrilinosWrappers::MPI::Vector    distributed_solution (system_rhs_newton);
    distributed_solution = solution;

    constraints_hanging_nodes.set_zero (distributed_solution);
    constraints_hanging_nodes.set_zero (system_rhs_newton);

    MPI_Barrier (mpi_communicator);
    t.restart();

    preconditioner_u.initialize (system_matrix_newton, additional_data);

    MPI_Barrier (mpi_communicator);
    t.stop();
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      run_time[6] += t.wall_time();

    MPI_Barrier (mpi_communicator);
    t.restart();

    PrimitiveVectorMemory<TrilinosWrappers::MPI::Vector> mem;
    TrilinosWrappers::MPI::Vector    tmp (system_rhs_newton);
    const double solver_tolerance = 1e-3 *
                                    system_matrix_newton.residual (tmp, distributed_solution, system_rhs_newton);

    SolverControl solver_control (system_matrix_newton.m(), solver_tolerance);
    SolverFGMRES<TrilinosWrappers::MPI::Vector>
    solver(solver_control, mem,
           SolverFGMRES<TrilinosWrappers::MPI::Vector>::
           AdditionalData(30, true));
    solver.solve(system_matrix_newton, distributed_solution, system_rhs_newton, preconditioner_u);

    pcout << "         Error: " << solver_control.initial_value()
          << " -> " << solver_control.last_value()
          << " in " << solver_control.last_step()
          << " FGMRES iterations."
          << std::endl;

    MPI_Barrier (mpi_communicator);
    t.stop();
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      run_time[7] += t.wall_time();

    number_iterations += solver_control.last_step();

    constraints.distribute (distributed_solution);

    solution = distributed_solution;
  }



  template <int dim>
  void PlasticityContactProblem<dim>::solve_newton ()
  {
    double                         resid=0;
    double                         resid_old=100000;
    TrilinosWrappers::MPI::Vector  res (system_rhs_newton);
    TrilinosWrappers::MPI::Vector  tmp_vector (system_rhs_newton);
    Timer                          t;

    std::vector<std::vector<bool> > constant_modes;
    DoFTools::extract_constant_modes (dof_handler,
                                      ComponentMask(),
                                      constant_modes);

    additional_data.elliptic = true;
    additional_data.n_cycles = 1;
    additional_data.w_cycle = false;
    additional_data.output_details = false;
    additional_data.smoother_sweeps = 2;
    additional_data.aggregation_threshold = 1e-2;

    IndexSet                        active_set_old (active_set);
    unsigned int j = 0;
    unsigned int number_assemble_system = 0;
    for (; j<=100; j++)
      {
        pcout<< " " <<std::endl;
        pcout<< "   Newton iteration " << j <<std::endl;
        pcout<< "      Updating active set..." <<std::endl;

        MPI_Barrier (mpi_communicator);
        t.restart();

        update_solution_and_constraints ();

        MPI_Barrier (mpi_communicator);
        t.stop();
        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
          run_time[5] += t.wall_time();

        pcout<< "      Assembling system... " <<std::endl;
        MPI_Barrier (mpi_communicator);
        t.restart();
        system_matrix_newton = 0;
        system_rhs_newton = 0;
        assemble_nl_system (solution);  //compute Newton-Matrix
        MPI_Barrier (mpi_communicator);
        t.stop();
        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
          run_time[1] += t.wall_time();

        number_assemble_system += 1;

        MPI_Barrier (mpi_communicator);
        t.restart();
        pcout<< "      Solving system... " <<std::endl;
        solve ();
        MPI_Barrier (mpi_communicator);
        t.stop();
        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
          run_time[2] += t.wall_time();

        TrilinosWrappers::MPI::Vector    distributed_solution (system_rhs_newton);
        distributed_solution = solution;

        int damped = 0;
        tmp_vector = old_solution;
        double a = 0;
        for (unsigned int i=0; (i<10)&&(!damped); i++)
          {
            a=std::pow(0.5, static_cast<double>(i));
            old_solution = tmp_vector;
            old_solution.sadd(1-a,a, distributed_solution);

            MPI_Barrier (mpi_communicator);
            t.restart();
            system_rhs_newton = 0;
            solution = old_solution;
            residual_nl_system (solution);
            res = system_rhs_newton;

            const unsigned int
            start_res     = (res.local_range().first),
            end_res       = (res.local_range().second);
            for (unsigned int n=start_res; n<end_res; ++n)
              if (constraints.is_inhomogeneously_constrained (n))
                res(n) = 0;

            resid = res.l2_norm ();

            if (resid<resid_old)
              damped=1;

            MPI_Barrier (mpi_communicator);
            t.stop();
            if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
              run_time[3] += t.wall_time();

            pcout << "      Residual of the non-contact part of the system: " << resid
                              << std::endl
                              << "         with a damping parameter alpha = " << a
                              << std::endl;
          }

        if (resid<1e-8)
          break;

        resid_old=resid;

        resid_vector = system_rhs_newton;

        if (active_set == active_set_old && resid < 1e-10)
          break;
        active_set_old = active_set;
      } // End of active-set-loop

    pcout << "" << std::endl
          << "      Number of assembled systems = " << number_assemble_system
          << std::endl
          << "      Number of Solver-Iterations = " << number_iterations << std::endl;
  }



  template <int dim>
  void PlasticityContactProblem<dim>::refine_grid ()
  {
    Vector<float> estimated_error_per_cell (triangulation.n_active_cells());
    KellyErrorEstimator<dim>::estimate (dof_handler,
                                        QGauss<dim-1>(3),
                                        typename FunctionMap<dim>::type(),
                                        solution,
                                        estimated_error_per_cell);
    parallel::distributed::GridRefinement::
    refine_and_coarsen_fixed_number (triangulation,
                                     estimated_error_per_cell,
                                     0.3, 0.03);
    triangulation.execute_coarsening_and_refinement ();

  }



  template <int dim>
  void PlasticityContactProblem<dim>::move_mesh (const TrilinosWrappers::MPI::Vector &_complete_displacement) const
  {
    std::vector<bool> vertex_touched (triangulation.n_vertices(),
                                      false);

    for (typename DoFHandler<dim>::active_cell_iterator
         cell = dof_handler.begin_active ();
         cell != dof_handler.end(); ++cell)
      if (cell->is_locally_owned())
        for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
          {
            if (vertex_touched[cell->vertex_index(v)] == false)
              {
                vertex_touched[cell->vertex_index(v)] = true;

                Point<dim> vertex_displacement;
                for (unsigned int d=0; d<dim; ++d)
                  {
                    if (_complete_displacement(cell->vertex_dof_index(v,d)) != 0)
                      vertex_displacement[d]
                        = _complete_displacement(cell->vertex_dof_index(v,d));
                  }

                cell->vertex(v) += vertex_displacement;
              }
          }
  }



  template <int dim>
  void PlasticityContactProblem<dim>::output_results (const std::string &title) const
  {
    move_mesh (solution);

    TrilinosWrappers::MPI::Vector         lambda (solution);
    lambda = resid_vector;

    DataOut<dim> data_out;

    data_out.attach_dof_handler (dof_handler);

    const std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation
    (dim, DataComponentInterpretation::component_is_part_of_vector);
    data_out.add_data_vector (solution, std::vector<std::string>(dim, "Displacement"),
                              DataOut<dim>::type_dof_data,
                              data_component_interpretation);
    data_out.add_data_vector (lambda, std::vector<std::string>(dim, "Residual"),
                              DataOut<dim>::type_dof_data,
                              data_component_interpretation);
    data_out.add_data_vector (active_set, std::vector<std::string>(dim, "ActiveSet"),
                              DataOut<dim>::type_dof_data,
                              data_component_interpretation);

    Vector<float> subdomain (triangulation.n_active_cells());
    for (unsigned int i=0; i<subdomain.size(); ++i)
      subdomain(i) = triangulation.locally_owned_subdomain();
    data_out.add_data_vector (subdomain, "subdomain");

    data_out.build_patches ();

    const std::string filename = (title + "-" +
                                  Utilities::int_to_string
                                  (triangulation.locally_owned_subdomain(), 4));

    std::ofstream output_vtu ((filename + ".vtu").c_str ());
    data_out.write_vtu (output_vtu);

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
      {
        std::vector<std::string> filenames;
        for (unsigned int i=0;
             i<Utilities::MPI::n_mpi_processes(mpi_communicator);
             ++i)
          filenames.push_back (title + "-" +
                               Utilities::int_to_string (i, 4) +
                               ".vtu");

        std::ofstream master_output ((filename + ".pvtu").c_str());
        data_out.write_pvtu_record (master_output, filenames);
      }

    TrilinosWrappers::MPI::Vector  tmp (solution);
    tmp *= -1;
    move_mesh (tmp);
  }



  template <int dim>
  void PlasticityContactProblem<dim>::run ()
  {
    pcout << "Read the obstacle from a file." << std::endl;
    input_obstacle.reset (new Input<dim>("obstacle_file.dat"));
    pcout << "Ostacle is available now." << std::endl;

    Timer             t;
    run_time.resize (8);

    const unsigned int n_cycles = 6;
    for (unsigned int cycle=0; cycle<n_cycles; ++cycle)
      {
        pcout << "" <<std::endl;
        pcout << "Cycle " << cycle << ':' << std::endl;

        MPI_Barrier (mpi_communicator);
        t.restart();

        if (cycle == 0)
          {
            make_grid();
          }
        else
          refine_grid ();

        setup_system ();

        MPI_Barrier (mpi_communicator);
        t.stop();
        run_time[0] += t.wall_time();;

        solve_newton ();

        pcout<< "      Writing graphical output..." <<std::endl;
        MPI_Barrier (mpi_communicator);
        t.restart();
        std::ostringstream filename_solution;
        filename_solution << "solution-";
        filename_solution << cycle;
        output_results (filename_solution.str ());
        MPI_Barrier (mpi_communicator);
        t.stop();
        if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
          run_time[4] += t.wall_time();

        pcout << "      Computing time for:" << std::endl
              << "         making grid and setup = " << run_time[0] << std::endl
              << "         updating active set = " << run_time[5] <<std::endl
              << "         assembling system = " << run_time[1] <<std::endl
              << "         solving system = " << run_time[2] <<std::endl
              << "         preconditioning = " << run_time[6] <<std::endl
              << "         solving with FGMRES = " << run_time[7] <<std::endl
              << "         computing error and lambda = " << run_time[3] <<std::endl
              << "         writing graphical output = " << run_time[4] <<std::endl;
      }
  }
}

// @sect3{The <code>main</code> function}

int main (int argc, char *argv[])
{
  using namespace dealii;
  using namespace Step42;

  deallog.depth_console (0);

  clock_t     start, end;

  start = clock();

  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);
  {
    int _n_refinements_global = 3;

    if (argc == 3)
      {
        _n_refinements_global = atoi(argv[1]);
      }

    PlasticityContactProblem<3> laplace_problem_3d (_n_refinements_global);
    laplace_problem_3d.run ();
  }

  end = clock();
  cout<< "%%%%%% Rechenzeit overall = " << (double)(end-start)/CLOCKS_PER_SEC <<std::endl;

  return 0;
}
