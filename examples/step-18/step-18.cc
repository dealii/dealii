/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2000 - 2013 by the deal.II authors
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
 * Author: Wolfgang Bangerth, University of Texas at Austin, 2000, 2004, 2005,
 * Timo Heister, 2013
 */


// First the usual list of header files that have already been used in
// previous example programs:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_solver.h>
#include <deal.II/lac/petsc_precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

// And here the only two new things among the header files: an include file in
// which symmetric tensors of rank 2 and 4 are implemented, as introduced in
// the introduction:
#include <deal.II/base/symmetric_tensor.h>

// And a header that implements filters for iterators looping over all
// cells. We will use this when selecting only those cells for output that are
// owned by the present process in a %parallel program:
#include <deal.II/grid/filtered_iterator.h>

// This is then simply C++ again:
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>

// The last step is as in all previous programs:
namespace Step18
{
  using namespace dealii;

  // @sect3{The <code>PointHistory</code> class}

  // As was mentioned in the introduction, we have to store the old stress in
  // quadrature point so that we can compute the residual forces at this point
  // during the next time step. This alone would not warrant a structure with
  // only one member, but in more complicated applications, we would have to
  // store more information in quadrature points as well, such as the history
  // variables of plasticity, etc. In essence, we have to store everything
  // that affects the present state of the material here, which in plasticity
  // is determined by the deformation history variables.
  //
  // We will not give this class any meaningful functionality beyond being
  // able to store data, i.e. there are no constructors, destructors, or other
  // member functions. In such cases of `dumb' classes, we usually opt to
  // declare them as <code>struct</code> rather than <code>class</code>, to
  // indicate that they are closer to C-style structures than C++-style
  // classes.
  template <int dim>
  struct PointHistory
  {
    SymmetricTensor<2,dim> old_stress;
  };


  // @sect3{The stress-strain tensor}

  // Next, we define the linear relationship between the stress and the strain
  // in elasticity. It is given by a tensor of rank 4 that is usually written
  // in the form $C_{ijkl} = \mu (\delta_{ik} \delta_{jl} + \delta_{il}
  // \delta_{jk}) + \lambda \delta_{ij} \delta_{kl}$. This tensor maps
  // symmetric tensor of rank 2 to symmetric tensors of rank 2. A function
  // implementing its creation for given values of the Lame constants $\lambda$
  // and $\mu$ is straightforward:
  template <int dim>
  SymmetricTensor<4,dim>
  get_stress_strain_tensor (const double lambda, const double mu)
  {
    SymmetricTensor<4,dim> tmp;
    for (unsigned int i=0; i<dim; ++i)
      for (unsigned int j=0; j<dim; ++j)
        for (unsigned int k=0; k<dim; ++k)
          for (unsigned int l=0; l<dim; ++l)
            tmp[i][j][k][l] = (((i==k) && (j==l) ? mu : 0.0) +
                               ((i==l) && (j==k) ? mu : 0.0) +
                               ((i==j) && (k==l) ? lambda : 0.0));
    return tmp;
  }

  // With this function, we will define a static member variable of the main
  // class below that will be used throughout the program as the stress-strain
  // tensor. Note that in more elaborate programs, this will probably be a
  // member variable of some class instead, or a function that returns the
  // stress-strain relationship depending on other input. For example in
  // damage theory models, the Lame constants are considered a function of the
  // prior stress/strain history of a point. Conversely, in plasticity the
  // form of the stress-strain tensor is modified if the material has reached
  // the yield stress in a certain point, and possibly also depending on its
  // prior history.
  //
  // In the present program, however, we assume that the material is
  // completely elastic and linear, and a constant stress-strain tensor is
  // sufficient for our present purposes.



  // @sect3{Auxiliary functions}

  // Before the rest of the program, here are a few functions that we need as
  // tools. These are small functions that are called in inner loops, so we
  // mark them as <code>inline</code>.
  //
  // The first one computes the symmetric strain tensor for shape function
  // <code>shape_func</code> at quadrature point <code>q_point</code> by
  // forming the symmetric gradient of this shape function. We need that when
  // we want to form the matrix, for example.
  //
  // We should note that in previous examples where we have treated
  // vector-valued problems, we have always asked the finite element object in
  // which of the vector component the shape function is actually non-zero,
  // and thereby avoided to compute any terms that we could prove were zero
  // anyway. For this, we used the <code>fe.system_to_component_index</code>
  // function that returns in which component a shape function was zero, and
  // also that the <code>fe_values.shape_value</code> and
  // <code>fe_values.shape_grad</code> functions only returned the value and
  // gradient of the single non-zero component of a shape function if this is
  // a vector-valued element.
  //
  // This was an optimization, and if it isn't terribly time critical, we can
  // get away with a simpler technique: just ask the <code>fe_values</code>
  // for the value or gradient of a given component of a given shape function
  // at a given quadrature point. This is what the
  // <code>fe_values.shape_grad_component(shape_func,q_point,i)</code> call
  // does: return the full gradient of the <code>i</code>th component of shape
  // function <code>shape_func</code> at quadrature point
  // <code>q_point</code>. If a certain component of a certain shape function
  // is always zero, then this will simply always return zero.
  //
  // As mentioned, using <code>fe_values.shape_grad_component</code> instead
  // of the combination of <code>fe.system_to_component_index</code> and
  // <code>fe_values.shape_grad</code> may be less efficient, but its
  // implementation is optimized for such cases and shouldn't be a big
  // slowdown. We demonstrate the technique here since it is so much simpler
  // and straightforward.
  template <int dim>
  inline
  SymmetricTensor<2,dim>
  get_strain (const FEValues<dim> &fe_values,
              const unsigned int   shape_func,
              const unsigned int   q_point)
  {
    // Declare a temporary that will hold the return value:
    SymmetricTensor<2,dim> tmp;

    // First, fill diagonal terms which are simply the derivatives in
    // direction <code>i</code> of the <code>i</code> component of the
    // vector-valued shape function:
    for (unsigned int i=0; i<dim; ++i)
      tmp[i][i] = fe_values.shape_grad_component (shape_func,q_point,i)[i];

    // Then fill the rest of the strain tensor. Note that since the tensor is
    // symmetric, we only have to compute one half (here: the upper right
    // corner) of the off-diagonal elements, and the implementation of the
    // <code>SymmetricTensor</code> class makes sure that at least to the
    // outside the symmetric entries are also filled (in practice, the class
    // of course stores only one copy). Here, we have picked the upper right
    // half of the tensor, but the lower left one would have been just as
    // good:
    for (unsigned int i=0; i<dim; ++i)
      for (unsigned int j=i+1; j<dim; ++j)
        tmp[i][j]
          = (fe_values.shape_grad_component (shape_func,q_point,i)[j] +
             fe_values.shape_grad_component (shape_func,q_point,j)[i]) / 2;

    return tmp;
  }


  // The second function does something very similar (and therefore is given
  // the same name): compute the symmetric strain tensor from the gradient of
  // a vector-valued field. If you already have a solution field, the
  // <code>fe_values.get_function_grads</code> function allows you to extract
  // the gradients of each component of your solution field at a quadrature
  // point. It returns this as a vector of rank-1 tensors: one rank-1 tensor
  // (gradient) per vector component of the solution. From this we have to
  // reconstruct the (symmetric) strain tensor by transforming the data
  // storage format and symmetrization. We do this in the same way as above,
  // i.e. we avoid a few computations by filling first the diagonal and then
  // only one half of the symmetric tensor (the <code>SymmetricTensor</code>
  // class makes sure that it is sufficient to write only one of the two
  // symmetric components).
  //
  // Before we do this, though, we make sure that the input has the kind of
  // structure we expect: that is that there are <code>dim</code> vector
  // components, i.e. one displacement component for each coordinate
  // direction. We test this with the <code>Assert</code> macro that will
  // simply abort our program if the condition is not met.
  template <int dim>
  inline
  SymmetricTensor<2,dim>
  get_strain (const std::vector<Tensor<1,dim> > &grad)
  {
    Assert (grad.size() == dim, ExcInternalError());

    SymmetricTensor<2,dim> strain;
    for (unsigned int i=0; i<dim; ++i)
      strain[i][i] = grad[i][i];

    for (unsigned int i=0; i<dim; ++i)
      for (unsigned int j=i+1; j<dim; ++j)
        strain[i][j] = (grad[i][j] + grad[j][i]) / 2;

    return strain;
  }


  // Finally, below we will need a function that computes the rotation matrix
  // induced by a displacement at a given point. In fact, of course, the
  // displacement at a single point only has a direction and a magnitude, it
  // is the change in direction and magnitude that induces rotations. In
  // effect, the rotation matrix can be computed from the gradients of a
  // displacement, or, more specifically, from the curl.
  //
  // The formulas by which the rotation matrices are determined are a little
  // awkward, especially in 3d. For 2d, there is a simpler way, so we
  // implement this function twice, once for 2d and once for 3d, so that we
  // can compile and use the program in both space dimensions if so desired --
  // after all, deal.II is all about dimension independent programming and
  // reuse of algorithm thoroughly tested with cheap computations in 2d, for
  // the more expensive computations in 3d. Here is one case, where we have to
  // implement different algorithms for 2d and 3d, but then can write the rest
  // of the program in a way that is independent of the space dimension.
  //
  // So, without further ado to the 2d implementation:
  Tensor<2,2>
  get_rotation_matrix (const std::vector<Tensor<1,2> > &grad_u)
  {
    // First, compute the curl of the velocity field from the gradients. Note
    // that we are in 2d, so the rotation is a scalar:
    const double curl = (grad_u[1][0] - grad_u[0][1]);

    // From this, compute the angle of rotation:
    const double angle = std::atan (curl);

    // And from this, build the antisymmetric rotation matrix:
    const double t[2][2] = {{ cos(angle), sin(angle) },
      {-sin(angle), cos(angle) }
    };
    return Tensor<2,2>(t);
  }


  // The 3d case is a little more contrived:
  Tensor<2,3>
  get_rotation_matrix (const std::vector<Tensor<1,3> > &grad_u)
  {
    // Again first compute the curl of the velocity field. This time, it is a
    // real vector:
    const Point<3> curl (grad_u[2][1] - grad_u[1][2],
                         grad_u[0][2] - grad_u[2][0],
                         grad_u[1][0] - grad_u[0][1]);

    // From this vector, using its magnitude, compute the tangent of the angle
    // of rotation, and from it the actual angle:
    const double tan_angle = std::sqrt(curl*curl);
    const double angle = std::atan (tan_angle);

    // Now, here's one problem: if the angle of rotation is too small, that
    // means that there is no rotation going on (for example a translational
    // motion). In that case, the rotation matrix is the identity matrix.
    //
    // The reason why we stress that is that in this case we have that
    // <code>tan_angle==0</code>. Further down, we need to divide by that
    // number in the computation of the axis of rotation, and we would get
    // into trouble when dividing doing so. Therefore, let's shortcut this and
    // simply return the identity matrix if the angle of rotation is really
    // small:
    if (angle < 1e-9)
      {
        static const double rotation[3][3]
        = {{ 1, 0, 0}, { 0, 1, 0 }, { 0, 0, 1 } };
        static const Tensor<2,3> rot(rotation);
        return rot;
      }

    // Otherwise compute the real rotation matrix. The algorithm for this is
    // not exactly obvious, but can be found in a number of books,
    // particularly on computer games where rotation is a very frequent
    // operation. Online, you can find a description at
    // http://www.makegames.com/3drotation/ and (this particular form, with
    // the signs as here) at
    // http://www.gamedev.net/reference/articles/article1199.asp:
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    const double t = 1-c;

    const Point<3> axis = curl/tan_angle;
    const double rotation[3][3]
    = {{
        t *axis[0] *axis[0]+c,
        t *axis[0] *axis[1]+s *axis[2],
        t *axis[0] *axis[2]-s *axis[1]
      },
      {
        t *axis[0] *axis[1]-s *axis[2],
        t *axis[1] *axis[1]+c,
        t *axis[1] *axis[2]+s *axis[0]
      },
      {
        t *axis[0] *axis[2]+s *axis[1],
        t *axis[1] *axis[1]-s *axis[0],
        t *axis[2] *axis[2]+c
      }
    };
    return Tensor<2,3>(rotation);
  }



  // @sect3{The <code>TopLevel</code> class}

  // This is the main class of the program. Since the namespace already
  // indicates what problem we are solving, let's call it by what it does: it
  // directs the flow of the program, i.e. it is the toplevel driver.
  //
  // The member variables of this class are essentially as before, i.e. it has
  // to have a triangulation, a DoF handler and associated objects such as
  // constraints, variables that describe the linear system, etc. There are a
  // good number of more member functions now, which we will explain below.
  //
  // The external interface of the class, however, is unchanged: it has a
  // public constructor and desctructor, and it has a <code>run</code>
  // function that initiated all the work.
  template <int dim>
  class TopLevel
  {
  public:
    TopLevel ();
    ~TopLevel ();
    void run ();

  private:
    // The private interface is more extensive than in step-17. First, we
    // obviously need functions that create the initial mesh, set up the
    // variables that describe the linear system on the present mesh
    // (i.e. matrices and vectors), and then functions that actually assemble
    // the system, direct what has to be solved in each time step, a function
    // that solves the linear system that arises in each timestep (and returns
    // the number of iterations it took), and finally output the solution
    // vector on the correct mesh:
    void create_coarse_grid ();

    void setup_system ();

    void assemble_system ();

    void solve_timestep ();

    unsigned int solve_linear_problem ();

    void output_results () const;

    // All, except for the first two, of these functions are called in each
    // timestep. Since the first time step is a little special, we have
    // separate functions that describe what has to happen in a timestep: one
    // for the first, and one for all following timesteps:
    void do_initial_timestep ();

    void do_timestep ();

    // Then we need a whole bunch of functions that do various things. The
    // first one refines the initial grid: we start on the coarse grid with a
    // pristine state, solve the problem, then look at it and refine the mesh
    // accordingly, and start the same process over again, again with a
    // pristine state. Thus, refining the initial mesh is somewhat simpler
    // than refining a grid between two successive time steps, since it does
    // not involve transferring data from the old to the new triangulation, in
    // particular the history data that is stored in each quadrature point.
    void refine_initial_grid ();

    // At the end of each time step, we want to move the mesh vertices around
    // according to the incremental displacement computed in this time
    // step. This is the function in which this is done:
    void move_mesh ();

    // Next are two functions that handle the history variables stored in each
    // quadrature point. The first one is called before the first timestep to
    // set up a pristine state for the history variables. It only works on
    // those quadrature points on cells that belong to the present processor:
    void setup_quadrature_point_history ();

    // The second one updates the history variables at the end of each
    // timestep:
    void update_quadrature_point_history ();

    // After the member functions, here are the member variables. The first
    // ones have all been discussed in more detail in previous example
    // programs:
    Triangulation<dim>   triangulation;

    FESystem<dim>        fe;

    DoFHandler<dim>      dof_handler;

    ConstraintMatrix     hanging_node_constraints;

    // One difference of this program is that we declare the quadrature
    // formula in the class declaration. The reason is that in all the other
    // programs, it didn't do much harm if we had used different quadrature
    // formulas when computing the matrix and the right hand side, for
    // example. However, in the present case it does: we store information in
    // the quadrature points, so we have to make sure all parts of the program
    // agree on where they are and how many there are on each cell. Thus, let
    // us first declare the quadrature formula that will be used throughout...
    const QGauss<dim>          quadrature_formula;

    // ... and then also have a vector of history objects, one per quadrature
    // point on those cells for which we are responsible (i.e. we don't store
    // history data for quadrature points on cells that are owned by other
    // processors).
    std::vector<PointHistory<dim> > quadrature_point_history;

    // The way this object is accessed is through a <code>user pointer</code>
    // that each cell, face, or edge holds: it is a <code>void*</code> pointer
    // that can be used by application programs to associate arbitrary data to
    // cells, faces, or edges. What the program actually does with this data
    // is within its own responsibility, the library just allocates some space
    // for these pointers, and application programs can set and read the
    // pointers for each of these objects.


    // Further: we need the objects of linear systems to be solved,
    // i.e. matrix, right hand side vector, and the solution vector. Since we
    // anticipate solving big problems, we use the same types as in step-17,
    // i.e. distributed %parallel matrices and vectors built on top of the
    // PETSc library. Conveniently, they can also be used when running on only
    // a single machine, in which case this machine happens to be the only one
    // in our %parallel universe.
    //
    // However, as a difference to step-17, we do not store the solution
    // vector -- which here is the incremental displacements computed in each
    // time step -- in a distributed fashion. I.e., of course it must be a
    // distributed vector when computing it, but immediately after that we
    // make sure each processor has a complete copy. The reason is that we had
    // already seen in step-17 that many functions needed a complete
    // copy. While it is not hard to get it, this requires communication on
    // the network, and is thus slow. In addition, these were repeatedly the
    // same operations, which is certainly undesirable unless the gains of not
    // always having to store the entire vector outweighs it. When writing
    // this program, it turned out that we need a complete copy of the
    // solution in so many places that it did not seem worthwhile to only get
    // it when necessary. Instead, we opted to obtain the complete copy once
    // and for all, and instead get rid of the distributed copy
    // immediately. Thus, note that the declaration of
    // <code>inremental_displacement</code> does not denote a distribute
    // vector as would be indicated by the middle namespace <code>MPI</code>:
    PETScWrappers::MPI::SparseMatrix system_matrix;

    PETScWrappers::MPI::Vector       system_rhs;

    PETScWrappers::Vector            incremental_displacement;

    // The next block of variables is then related to the time dependent
    // nature of the problem: they denote the length of the time interval
    // which we want to simulate, the present time and number of time step,
    // and length of present timestep:
    double       present_time;
    double       present_timestep;
    double       end_time;
    unsigned int timestep_no;

    // Then a few variables that have to do with %parallel processing: first,
    // a variable denoting the MPI communicator we use, and then two numbers
    // telling us how many participating processors there are, and where in
    // this world we are. Finally, a stream object that makes sure only one
    // processor is actually generating output to the console. This is all the
    // same as in step-17:
    MPI_Comm mpi_communicator;

    const unsigned int n_mpi_processes;

    const unsigned int this_mpi_process;

    ConditionalOStream pcout;

    // Here is a vector where each entry denotes the numbers of degrees of
    // freedom that are stored on the processor with that particular number:
    std::vector<types::global_dof_index> local_dofs_per_process;

    // Next, how many degrees of freedom the present processor stores. This
    // is, of course, an abbreviation to
    // <code>local_dofs_per_process[this_mpi_process]</code>.
    types::global_dof_index n_local_dofs;

    // In the same direction, also cache how many cells the present processor
    // owns. Note that the cells that belong to a processor are not
    // necessarily contiguously numbered (when iterating over them using
    // <code>active_cell_iterator</code>).
    unsigned int         n_local_cells;

    // Finally, we have a static variable that denotes the linear relationship
    // between the stress and strain. Since it is a constant object that does
    // not depend on any input (at least not in this program), we make it a
    // static variable and will initialize it in the same place where we
    // define the constructor of this class:
    static const SymmetricTensor<4,dim> stress_strain_tensor;
  };


  // @sect3{The <code>BodyForce</code> class}

  // Before we go on to the main functionality of this program, we have to
  // define what forces will act on the body whose deformation we want to
  // study. These may either be body forces or boundary forces. Body forces
  // are generally mediated by one of the four basic physical types of forces:
  // gravity, strong and weak interaction, and electromagnetism. Unless one
  // wants to consider subatomic objects (for which quasistatic deformation is
  // irrelevant and an inappropriate description anyway), only gravity and
  // electromagnetic forces need to be considered. Let us, for simplicity
  // assume that our body has a certain mass density, but is either
  // non-magnetic and not electrically conducting or that there are no
  // significant electromagnetic fields around. In that case, the body forces
  // are simply <code>rho g</code>, where <code>rho</code> is the material
  // density and <code>g</code> is a vector in negative z-direction with
  // magnitude 9.81 m/s^2.  Both the density and <code>g</code> are defined in
  // the function, and we take as the density 7700 kg/m^3, a value commonly
  // assumed for steel.
  //
  // To be a little more general and to be able to do computations in 2d as
  // well, we realize that the body force is always a function returning a
  // <code>dim</code> dimensional vector. We assume that gravity acts along
  // the negative direction of the last, i.e. <code>dim-1</code>th
  // coordinate. The rest of the implementation of this function should be
  // mostly self-explanatory given similar definitions in previous example
  // programs. Note that the body force is independent of the location; to
  // avoid compiler warnings about unused function arguments, we therefore
  // comment out the name of the first argument of the
  // <code>vector_value</code> function:
  template <int dim>
  class BodyForce :  public Function<dim>
  {
  public:
    BodyForce ();

    virtual
    void
    vector_value (const Point<dim> &p,
                  Vector<double>   &values) const;

    virtual
    void
    vector_value_list (const std::vector<Point<dim> > &points,
                       std::vector<Vector<double> >   &value_list) const;
  };


  template <int dim>
  BodyForce<dim>::BodyForce ()
    :
    Function<dim> (dim)
  {}


  template <int dim>
  inline
  void
  BodyForce<dim>::vector_value (const Point<dim> &/*p*/,
                                Vector<double>   &values) const
  {
    Assert (values.size() == dim,
            ExcDimensionMismatch (values.size(), dim));

    const double g   = 9.81;
    const double rho = 7700;

    values = 0;
    values(dim-1) = -rho * g;
  }



  template <int dim>
  void
  BodyForce<dim>::vector_value_list (const std::vector<Point<dim> > &points,
                                     std::vector<Vector<double> >   &value_list) const
  {
    const unsigned int n_points = points.size();

    Assert (value_list.size() == n_points,
            ExcDimensionMismatch (value_list.size(), n_points));

    for (unsigned int p=0; p<n_points; ++p)
      BodyForce<dim>::vector_value (points[p],
                                    value_list[p]);
  }



  // @sect3{The <code>IncrementalBoundaryValue</code> class}

  // In addition to body forces, movement can be induced by boundary forces
  // and forced boundary displacement. The latter case is equivalent to forces
  // being chosen in such a way that they induce certain displacement.
  //
  // For quasistatic displacement, typical boundary forces would be pressure
  // on a body, or tangential friction against another body. We chose a
  // somewhat simpler case here: we prescribe a certain movement of (parts of)
  // the boundary, or at least of certain components of the displacement
  // vector. We describe this by another vector-valued function that, for a
  // given point on the boundary, returns the prescribed displacement.
  //
  // Since we have a time-dependent problem, the displacement increment of the
  // boundary equals the displacement accumulated during the length of the
  // timestep. The class therefore has to know both the present time and the
  // length of the present time step, and can then approximate the incremental
  // displacement as the present velocity times the present timestep.
  //
  // For the purposes of this program, we choose a simple form of boundary
  // displacement: we displace the top boundary with constant velocity
  // downwards. The rest of the boundary is either going to be fixed (and is
  // then described using an object of type <code>ZeroFunction</code>) or free
  // (Neumann-type, in which case nothing special has to be done).  The
  // implementation of the class describing the constant downward motion
  // should then be obvious using the knowledge we gained through all the
  // previous example programs:
  template <int dim>
  class IncrementalBoundaryValues :  public Function<dim>
  {
  public:
    IncrementalBoundaryValues (const double present_time,
                               const double present_timestep);

    virtual
    void
    vector_value (const Point<dim> &p,
                  Vector<double>   &values) const;

    virtual
    void
    vector_value_list (const std::vector<Point<dim> > &points,
                       std::vector<Vector<double> >   &value_list) const;

  private:
    const double velocity;
    const double present_time;
    const double present_timestep;
  };


  template <int dim>
  IncrementalBoundaryValues<dim>::
  IncrementalBoundaryValues (const double present_time,
                             const double present_timestep)
    :
    Function<dim> (dim),
    velocity (.1),
    present_time (present_time),
    present_timestep (present_timestep)
  {}


  template <int dim>
  void
  IncrementalBoundaryValues<dim>::
  vector_value (const Point<dim> &/*p*/,
                Vector<double>   &values) const
  {
    Assert (values.size() == dim,
            ExcDimensionMismatch (values.size(), dim));

    values = 0;
    values(2) = -present_timestep * velocity;
  }



  template <int dim>
  void
  IncrementalBoundaryValues<dim>::
  vector_value_list (const std::vector<Point<dim> > &points,
                     std::vector<Vector<double> >   &value_list) const
  {
    const unsigned int n_points = points.size();

    Assert (value_list.size() == n_points,
            ExcDimensionMismatch (value_list.size(), n_points));

    for (unsigned int p=0; p<n_points; ++p)
      IncrementalBoundaryValues<dim>::vector_value (points[p],
                                                    value_list[p]);
  }



  // @sect3{Implementation of the <code>TopLevel</code> class}

  // Now for the implementation of the main class. First, we initialize the
  // stress-strain tensor, which we have declared as a static const
  // variable. We chose Lame constants that are appropriate for steel:
  template <int dim>
  const SymmetricTensor<4,dim>
  TopLevel<dim>::stress_strain_tensor
    = get_stress_strain_tensor<dim> (/*lambda = */ 9.695e10,
                                                   /*mu     = */ 7.617e10);



  // @sect4{The public interface}

  // The next step is the definition of constructors and destructors. There
  // are no surprises here: we choose linear and continuous finite elements
  // for each of the <code>dim</code> vector components of the solution, and a
  // Gaussian quadrature formula with 2 points in each coordinate
  // direction. The destructor should be obvious:
  template <int dim>
  TopLevel<dim>::TopLevel ()
    :
    fe (FE_Q<dim>(1), dim),
    dof_handler (triangulation),
    quadrature_formula (2),
    mpi_communicator (MPI_COMM_WORLD),
    n_mpi_processes (Utilities::MPI::n_mpi_processes(mpi_communicator)),
    this_mpi_process (Utilities::MPI::this_mpi_process(mpi_communicator)),
    pcout (std::cout, this_mpi_process == 0)
  {}



  template <int dim>
  TopLevel<dim>::~TopLevel ()
  {
    dof_handler.clear ();
  }



  // The last of the public functions is the one that directs all the work,
  // <code>run()</code>. It initializes the variables that describe where in
  // time we presently are, then runs the first time step, then loops over all
  // the other time steps. Note that for simplicity we use a fixed time step,
  // whereas a more sophisticated program would of course have to choose it in
  // some more reasonable way adaptively:
  template <int dim>
  void TopLevel<dim>::run ()
  {
    present_time = 0;
    present_timestep = 1;
    end_time = 10;
    timestep_no = 0;

    do_initial_timestep ();

    while (present_time < end_time)
      do_timestep ();
  }


  // @sect4{TopLevel::create_coarse_grid}

  // The next function in the order in which they were declared above is the
  // one that creates the coarse grid from which we start. For this example
  // program, we want to compute the deformation of a cylinder under axial
  // compression. The first step therefore is to generate a mesh for a
  // cylinder of length 3 and with inner and outer radii of 0.8 and 1,
  // respectively. Fortunately, there is a library function for such a mesh.
  //
  // In a second step, we have to associated boundary conditions with the
  // upper and lower faces of the cylinder. We choose a boundary indicator of
  // 0 for the boundary faces that are characterized by their midpoints having
  // z-coordinates of either 0 (bottom face), an indicator of 1 for z=3 (top
  // face); finally, we use boundary indicator 2 for all faces on the inside
  // of the cylinder shell, and 3 for the outside.
  template <int dim>
  void TopLevel<dim>::create_coarse_grid ()
  {
    const double inner_radius = 0.8,
                 outer_radius = 1;
    GridGenerator::cylinder_shell (triangulation,
                                   3, inner_radius, outer_radius);
    for (typename Triangulation<dim>::active_cell_iterator
         cell=triangulation.begin_active();
         cell!=triangulation.end(); ++cell)
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        if (cell->face(f)->at_boundary())
          {
            const Point<dim> face_center = cell->face(f)->center();

            if (face_center[2] == 0)
              cell->face(f)->set_boundary_indicator (0);
            else if (face_center[2] == 3)
              cell->face(f)->set_boundary_indicator (1);
            else if (std::sqrt(face_center[0]*face_center[0] +
                               face_center[1]*face_center[1])
                     <
                     (inner_radius + outer_radius) / 2)
              cell->face(f)->set_boundary_indicator (2);
            else
              cell->face(f)->set_boundary_indicator (3);
          }

    // In order to make sure that new vertices are placed correctly on mesh
    // refinement, we have to associate objects describing those parts of the
    // boundary that do not consist of straight parts. Corresponding to the
    // cylinder shell generator function used above, there are classes that
    // can be used to describe the geometry of cylinders. We need to use
    // different objects for the inner and outer parts of the cylinder, with
    // different radii; the second argument to the constructor indicates the
    // axis around which the cylinder revolves -- in this case the
    // z-axis. Note that the boundary objects need to live as long as the
    // triangulation does; we can achieve this by making the objects static,
    // which means that they live as long as the program runs:
    static const CylinderBoundary<dim> inner_cylinder (inner_radius, 2);
    static const CylinderBoundary<dim> outer_cylinder (outer_radius, 2);
    // We then attach these two objects to the triangulation, and make them
    // correspond to boundary indicators 2 and 3:
    triangulation.set_boundary (2, inner_cylinder);
    triangulation.set_boundary (3, outer_cylinder);

    // There's one more thing we have to take care of (we should have done so
    // above already, but for didactic reasons it was more appropriate to
    // handle it after discussing boundary objects). %Boundary indicators in
    // deal.II, for mostly historic reasons, serve a dual purpose: they
    // describe the type of a boundary for other places in a program where
    // different boundary conditions are implemented; and they describe which
    // boundary object (as the ones associated above) should be queried when
    // new boundary points need to be placed upon mesh refinement. In the
    // prefix to this function, we have discussed the boundary condition
    // issue, and the boundary geometry issue was mentioned just above. But
    // there is a case where we have to be careful with geometry: what happens
    // if a cell is refined that has two faces with different boundary
    // indicators? For example one at the edges of the cylinder? In that case,
    // the library wouldn't know where to put new points in the middle of
    // edges (one of the twelve lines of a hexahedron). In fact, the library
    // doesn't even care about the boundary indicator of adjacent faces when
    // refining edges: it considers the boundary indicators associated with
    // the edges themselves. So what do we want to happen with the edges of
    // the cylinder shell: they sit on both faces with boundary indicators 2
    // or 3 (inner or outer shell) and 0 or 1 (for which no boundary objects
    // have been specified, and for which the library therefore assumes
    // straight lines). Obviously, we want these lines to follow the curved
    // shells, so we have to assign all edges along faces with boundary
    // indicators 2 or 3 these same boundary indicators to make sure they are
    // refined using the appropriate geometry objects. This is easily done:
    for (typename Triangulation<dim>::active_face_iterator
         face=triangulation.begin_active_face();
         face!=triangulation.end_face(); ++face)
      if (face->at_boundary())
        if ((face->boundary_indicator() == 2)
            ||
            (face->boundary_indicator() == 3))
          for (unsigned int edge = 0; edge<GeometryInfo<dim>::lines_per_face;
               ++edge)
            face->line(edge)
            ->set_boundary_indicator (face->boundary_indicator());

    // Once all this is done, we can refine the mesh once globally:
    triangulation.refine_global (1);


    // As the final step, we need to set up a clean state of the data that we
    // store in the quadrature points on all cells that are treated on the
    // present processor. To do so, we also have to know which processors are
    // ours in the first place. This is done in the following two function
    // calls:
    GridTools::partition_triangulation (n_mpi_processes, triangulation);
    setup_quadrature_point_history ();
  }




  // @sect4{TopLevel::setup_system}

  // The next function is the one that sets up the data structures for a given
  // mesh. This is done in most the same way as in step-17: distribute the
  // degrees of freedom, then sort these degrees of freedom in such a way that
  // each processor gets a contiguous chunk of them. Note that subdivisions into
  // chunks for each processor is handled in the functions that create or
  // refine grids, unlike in the previous example program (the point where
  // this happens is mostly a matter of taste; here, we chose to do it when
  // grids are created since in the <code>do_initial_timestep</code> and
  // <code>do_timestep</code> functions we want to output the number of cells
  // on each processor at a point where we haven't called the present function
  // yet).
  template <int dim>
  void TopLevel<dim>::setup_system ()
  {
    dof_handler.distribute_dofs (fe);
    DoFRenumbering::subdomain_wise (dof_handler);

    // The next thing is to store some information for later use on how many
    // cells or degrees of freedom the present processor, or any of the
    // processors has to work on. First the cells local to this processor...
    n_local_cells
      = GridTools::count_cells_with_subdomain_association (triangulation,
                                                           this_mpi_process);

    // ...and then a list of numbers of how many degrees of freedom each
    // processor has to handle:
    local_dofs_per_process.resize (n_mpi_processes);
    for (unsigned int i=0; i<n_mpi_processes; ++i)
      local_dofs_per_process[i]
        = DoFTools::count_dofs_with_subdomain_association (dof_handler, i);

    // Finally, make it easier to denote how many degrees of freedom the
    // present process has to deal with, by introducing an abbreviation:
    n_local_dofs = local_dofs_per_process[this_mpi_process];

    // The next step is to set up constraints due to hanging nodes. This has
    // been handled many times before:
    hanging_node_constraints.clear ();
    DoFTools::make_hanging_node_constraints (dof_handler,
                                             hanging_node_constraints);
    hanging_node_constraints.close ();

    // And then we have to set up the matrix. Here we deviate from step-17, in
    // which we simply used PETSc's ability to just know about the size of the
    // matrix and later allocate those nonzero elements that are being written
    // to. While this works just fine from a correctness viewpoint, it is not
    // at all efficient: if we don't give PETSc a clue as to which elements
    // are written to, it is (at least at the time of this writing) unbearably
    // slow when we set the elements in the matrix for the first time (i.e. in
    // the first timestep). Later on, when the elements have been allocated,
    // everything is much faster. In experiments we made, the first timestep
    // can be accelerated by almost two orders of magnitude if we instruct
    // PETSc which elements will be used and which are not.
    //
    // To do so, we first generate the sparsity pattern of the matrix we are
    // going to work with, and make sure that the condensation of hanging node
    // constraints add the necessary additional entries in the sparsity
    // pattern:
    CompressedSparsityPattern sparsity_pattern (dof_handler.n_dofs(),
                                                dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
    hanging_node_constraints.condense (sparsity_pattern);
    // Note that we have used the <code>CompressedSparsityPattern</code> class
    // here that was already introduced in step-11, rather than the
    // <code>SparsityPattern</code> class that we have used in all other
    // cases. The reason for this is that for the latter class to work we have
    // to give an initial upper bound for the number of entries in each row, a
    // task that is traditionally done by
    // <code>DoFHandler::max_couplings_between_dofs()</code>. However, this
    // function suffers from a serious problem: it has to compute an upper
    // bound to the number of nonzero entries in each row, and this is a
    // rather complicated task, in particular in 3d. In effect, while it is
    // quite accurate in 2d, it often comes up with much too large a number in
    // 3d, and in that case the <code>SparsityPattern</code> allocates much
    // too much memory at first, often several 100 MBs. This is later
    // corrected when <code>DoFTools::make_sparsity_pattern</code> is called
    // and we realize that we don't need all that much memory, but at time it
    // is already too late: for large problems, the temporary allocation of
    // too much memory can lead to out-of-memory situations.
    //
    // In order to avoid this, we resort to the
    // <code>CompressedSparsityPattern</code> class that is slower but does
    // not require any up-front estimate on the number of nonzero entries per
    // row. It therefore only ever allocates as much memory as it needs at any
    // given time, and we can build it even for large 3d problems.
    //
    // It is also worth noting that the sparsity pattern we construct is
    // global, i.e. comprises all degrees of freedom whether they will be
    // owned by the processor we are on or another one (in case this program
    // is run in %parallel via MPI). This of course is not optimal -- it
    // limits the size of the problems we can solve, since storing the entire
    // sparsity pattern (even if only for a short time) on each processor does
    // not scale well. However, there are several more places in the program
    // in which we do this, for example we always keep the global
    // triangulation and DoF handler objects around, even if we only work on
    // part of them. At present, deal.II does not have the necessary
    // facilities to completely distribute these objects (a task that, indeed,
    // is very hard to achieve with adaptive meshes, since well-balanced
    // subdivisions of a domain tend to become unbalanced as the mesh is
    // adaptively refined).
    //
    // With this data structure, we can then go to the PETSc sparse matrix and
    // tell it to preallocate all the entries we will later want to write to:
    system_matrix.reinit (mpi_communicator,
                          sparsity_pattern,
                          local_dofs_per_process,
                          local_dofs_per_process,
                          this_mpi_process);
    // After this point, no further explicit knowledge of the sparsity pattern
    // is required any more and we can let the <code>sparsity_pattern</code>
    // variable go out of scope without any problem.

    // The last task in this function is then only to reset the right hand
    // side vector as well as the solution vector to its correct size;
    // remember that the solution vector is a local one, unlike the right hand
    // side that is a distributed %parallel one and therefore needs to know
    // the MPI communicator over which it is supposed to transmit messages:
    system_rhs.reinit (mpi_communicator, dof_handler.n_dofs(), n_local_dofs);
    incremental_displacement.reinit (dof_handler.n_dofs());
  }



  // @sect4{TopLevel::assemble_system}

  // Again, assembling the system matrix and right hand side follows the same
  // structure as in many example programs before. In particular, it is mostly
  // equivalent to step-17, except for the different right hand side that now
  // only has to take into account internal stresses. In addition, assembling
  // the matrix is made significantly more transparent by using the
  // <code>SymmetricTensor</code> class: note the elegance of forming the
  // scalar products of symmetric tensors of rank 2 and 4. The implementation
  // is also more general since it is independent of the fact that we may or
  // may not be using an isotropic elasticity tensor.
  //
  // The first part of the assembly routine is as always:
  template <int dim>
  void TopLevel<dim>::assemble_system ()
  {
    system_rhs = 0;
    system_matrix = 0;

    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values   | update_gradients |
                             update_quadrature_points | update_JxW_values);

    const unsigned int   dofs_per_cell = fe.dofs_per_cell;
    const unsigned int   n_q_points    = quadrature_formula.size();

    FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       cell_rhs (dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

    BodyForce<dim>      body_force;
    std::vector<Vector<double> > body_force_values (n_q_points,
                                                    Vector<double>(dim));

    // As in step-17, we only need to loop over all cells that belong to the
    // present processor:
    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->subdomain_id() == this_mpi_process)
        {
          cell_matrix = 0;
          cell_rhs = 0;

          fe_values.reinit (cell);

          // Then loop over all indices i,j and quadrature points and assemble
          // the system matrix contributions from this cell.  Note how we
          // extract the symmetric gradients (strains) of the shape functions
          // at a given quadrature point from the <code>FEValues</code>
          // object, and the elegance with which we form the triple
          // contraction <code>eps_phi_i : C : eps_phi_j</code>; the latter
          // needs to be compared to the clumsy computations needed in
          // step-17, both in the introduction as well as in the respective
          // place in the program:
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              for (unsigned int q_point=0; q_point<n_q_points;
                   ++q_point)
                {
                  const SymmetricTensor<2,dim>
                  eps_phi_i = get_strain (fe_values, i, q_point),
                  eps_phi_j = get_strain (fe_values, j, q_point);

                  cell_matrix(i,j)
                  += (eps_phi_i * stress_strain_tensor * eps_phi_j
                      *
                      fe_values.JxW (q_point));
                }


          // Then also assemble the local right hand side contributions. For
          // this, we need to access the prior stress value in this quadrature
          // point. To get it, we use the user pointer of this cell that
          // points into the global array to the quadrature point data
          // corresponding to the first quadrature point of the present cell,
          // and then add an offset corresponding to the index of the
          // quadrature point we presently consider:
          const PointHistory<dim> *local_quadrature_points_data
            = reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());
          // In addition, we need the values of the external body forces at
          // the quadrature points on this cell:
          body_force.vector_value_list (fe_values.get_quadrature_points(),
                                        body_force_values);
          // Then we can loop over all degrees of freedom on this cell and
          // compute local contributions to the right hand side:
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              const unsigned int
              component_i = fe.system_to_component_index(i).first;

              for (unsigned int q_point=0; q_point<n_q_points; ++q_point)
                {
                  const SymmetricTensor<2,dim> &old_stress
                    = local_quadrature_points_data[q_point].old_stress;

                  cell_rhs(i) += (body_force_values[q_point](component_i) *
                                  fe_values.shape_value (i,q_point)
                                  -
                                  old_stress *
                                  get_strain (fe_values,i,q_point))
                                 *
                                 fe_values.JxW (q_point);
                }
            }

          // Now that we have the local contributions to the linear system, we
          // need to transfer it into the global objects. This is done exactly
          // as in step-17:
          cell->get_dof_indices (local_dof_indices);

          hanging_node_constraints
          .distribute_local_to_global (cell_matrix, cell_rhs,
                                       local_dof_indices,
                                       system_matrix, system_rhs);
        }

    // Now compress the vector and the system matrix:
    system_matrix.compress(VectorOperation::add);
    system_rhs.compress(VectorOperation::add);


    // The last step is to again fix up boundary values, just as we already
    // did in previous programs. A slight complication is that the
    // <code>apply_boundary_values</code> function wants to have a solution
    // vector compatible with the matrix and right hand side (i.e. here a
    // distributed %parallel vector, rather than the sequential vector we use
    // in this program) in order to preset the entries of the solution vector
    // with the correct boundary values. We provide such a compatible vector
    // in the form of a temporary vector which we then copy into the
    // sequential one.

    // We make up for this complication by showing how boundary values can be
    // used flexibly: following the way we create the triangulation, there are
    // three distinct boundary indicators used to describe the domain,
    // corresponding to the bottom and top faces, as well as the inner/outer
    // surfaces. We would like to impose boundary conditions of the following
    // type: The inner and outer cylinder surfaces are free of external
    // forces, a fact that corresponds to natural (Neumann-type) boundary
    // conditions for which we don't have to do anything. At the bottom, we
    // want no movement at all, corresponding to the cylinder being clamped or
    // cemented in at this part of the boundary. At the top, however, we want
    // a prescribed vertical downward motion compressing the cylinder; in
    // addition, we only want to restrict the vertical movement, but not the
    // horizontal ones -- one can think of this situation as a well-greased
    // plate sitting on top of the cylinder pushing it downwards: the atoms of
    // the cylinder are forced to move downward, but they are free to slide
    // horizontally along the plate.

    // The way to describe this is as follows: for boundary indicator zero
    // (bottom face) we use a dim-dimensional zero function representing no
    // motion in any coordinate direction. For the boundary with indicator 1
    // (top surface), we use the <code>IncrementalBoundaryValues</code> class,
    // but we specify an additional argument to the
    // <code>VectorTools::interpolate_boundary_values</code> function denoting
    // which vector components it should apply to; this is a vector of bools
    // for each vector component and because we only want to restrict vertical
    // motion, it has only its last component set:
    FEValuesExtractors::Scalar z_component (dim-1);
    std::map<types::global_dof_index,double> boundary_values;
    VectorTools::
    interpolate_boundary_values (dof_handler,
                                 0,
                                 ZeroFunction<dim> (dim),
                                 boundary_values);
    VectorTools::
    interpolate_boundary_values (dof_handler,
                                 1,
                                 IncrementalBoundaryValues<dim>(present_time,
                                                                present_timestep),
                                 boundary_values,
                                 fe.component_mask(z_component));

    PETScWrappers::MPI::Vector tmp (mpi_communicator, dof_handler.n_dofs(),
                                    n_local_dofs);
    MatrixTools::apply_boundary_values (boundary_values,
                                        system_matrix, tmp,
                                        system_rhs, false);
    incremental_displacement = tmp;
  }



  // @sect4{TopLevel::solve_timestep}

  // The next function is the one that controls what all has to happen within
  // a timestep. The order of things should be relatively self-explanatory
  // from the function names:
  template <int dim>
  void TopLevel<dim>::solve_timestep ()
  {
    pcout << "    Assembling system..." << std::flush;
    assemble_system ();
    pcout << " norm of rhs is " << system_rhs.l2_norm()
          << std::endl;

    const unsigned int n_iterations = solve_linear_problem ();

    pcout << "    Solver converged in " << n_iterations
          << " iterations." << std::endl;

    pcout << "    Updating quadrature point data..." << std::flush;
    update_quadrature_point_history ();
    pcout << std::endl;
  }



  // @sect4{TopLevel::solve_linear_problem}

  // Solving the linear system again works mostly as before. The only
  // difference is that we want to only keep a complete local copy of the
  // solution vector instead of the distributed one that we get as output from
  // PETSc's solver routines. To this end, we declare a local temporary
  // variable for the distributed vector and initialize it with the contents
  // of the local variable (remember that the
  // <code>apply_boundary_values</code> function called in
  // <code>assemble_system</code> preset the values of boundary nodes in this
  // vector), solve with it, and at the end of the function copy it again into
  // the complete local vector that we declared as a member variable. Hanging
  // node constraints are then distributed only on the local copy,
  // i.e. independently of each other on each of the processors:
  template <int dim>
  unsigned int TopLevel<dim>::solve_linear_problem ()
  {
    PETScWrappers::MPI::Vector
    distributed_incremental_displacement (mpi_communicator,
                                          dof_handler.n_dofs(),
                                          n_local_dofs);
    distributed_incremental_displacement = incremental_displacement;

    SolverControl           solver_control (dof_handler.n_dofs(),
                                            1e-16*system_rhs.l2_norm());
    PETScWrappers::SolverCG cg (solver_control,
                                mpi_communicator);

    PETScWrappers::PreconditionBlockJacobi preconditioner(system_matrix);

    cg.solve (system_matrix, distributed_incremental_displacement, system_rhs,
              preconditioner);

    incremental_displacement = distributed_incremental_displacement;

    hanging_node_constraints.distribute (incremental_displacement);

    return solver_control.last_step();
  }



  // @sect4{TopLevel::output_results}

  // This function generates the graphical output in .vtu format as explained
  // in the introduction. Each process will only work on the cells it owns,
  // and then write the result into a file of its own. Additionally, processor
  // 0 will write the record files the reference all the .vtu files.
  //
  // The crucial part of this function is to give the <code>DataOut</code>
  // class a way to only work on the cells that the present process owns. This
  // class is already well-equipped for that: it has two virtual functions
  // <code>first_cell</code> and <code>next_cell</code> that return the first
  // cell to be worked on, and given one cell return the next cell to be
  // worked on. By default, these functions return the first active cell
  // (i.e. the first one that has no children) and the next active cell. What
  // we have to do here is derive a class from <code>DataOut</code> that
  // overloads these two functions to only iterate over those cells with the
  // right subdomain indicator.
  //
  // We do this at the beginning of this function. The <code>first_cell</code>
  // function just starts with the first active cell, and then iterates to the
  // next cells while the cell presently under consideration does not yet have
  // the correct subdomain id. The only thing that needs to be taken care of
  // is that we don't try to keep iterating when we have hit the end iterator.
  //
  // The <code>next_cell</code> function could be implemented in a similar
  // way. However, we use this occasion as a pretext to introduce one more
  // thing that the library offers: filtered iterators. These are wrappers for
  // the iterator classes that just skip all cells (or faces, lines, etc) that
  // do not satisfy a certain predicate (a predicate in computer-lingo is a
  // function that when applied to a data element either returns true or
  // false). In the present case, the predicate is that the cell has to have a
  // certain subdomain id, and the library already has this predicate built
  // in. If the cell iterator is not the end iterator, what we then have to do
  // is to initialize such a filtered iterator with the present cell and the
  // predicate, and then increase the iterator exactly once. While the more
  // conventional loop would probably not have been much longer, this is
  // definitely the more elegant way -- and then, these example programs also
  // serve the purpose of introducing what is available in deal.II.
  template<int dim>
  class FilteredDataOut : public DataOut<dim>
  {
  public:
    FilteredDataOut (const unsigned int subdomain_id)
      :
      subdomain_id (subdomain_id)
    {}

    virtual typename DataOut<dim>::cell_iterator
    first_cell ()
    {
      typename DataOut<dim>::active_cell_iterator
      cell = this->dofs->begin_active();
      while ((cell != this->dofs->end()) &&
             (cell->subdomain_id() != subdomain_id))
        ++cell;

      return cell;
    }

    virtual typename DataOut<dim>::cell_iterator
    next_cell (const typename DataOut<dim>::cell_iterator &old_cell)
    {
      if (old_cell != this->dofs->end())
        {
          const IteratorFilters::SubdomainEqualTo
          predicate(subdomain_id);

          return
            ++(FilteredIterator
               <typename DataOut<dim>::active_cell_iterator>
               (predicate,old_cell));
        }
      else
        return old_cell;
    }

  private:
    const unsigned int subdomain_id;
  };



  template <int dim>
  void TopLevel<dim>::output_results () const
  {
    // With this newly defined class, declare an object that is going to
    // generate the graphical output and attach the dof handler with it from
    // which to get the solution vector:
    FilteredDataOut<dim> data_out(this_mpi_process);
    data_out.attach_dof_handler (dof_handler);

    // Then, just as in step-17, define the names of solution variables (which
    // here are the displacement increments) and queue the solution vector for
    // output. Note in the following switch how we make sure that if the space
    // dimension should be unhandled that we throw an exception saying that we
    // haven't implemented this case yet (another case of defensive
    // programming):
    std::vector<std::string> solution_names;
    switch (dim)
      {
      case 1:
        solution_names.push_back ("delta_x");
        break;
      case 2:
        solution_names.push_back ("delta_x");
        solution_names.push_back ("delta_y");
        break;
      case 3:
        solution_names.push_back ("delta_x");
        solution_names.push_back ("delta_y");
        solution_names.push_back ("delta_z");
        break;
      default:
        Assert (false, ExcNotImplemented());
      }

    data_out.add_data_vector (incremental_displacement,
                              solution_names);


    // The next thing is that we wanted to output something like the average
    // norm of the stresses that we have stored in each cell. This may seem
    // complicated, since on the present processor we only store the stresses
    // in quadrature points on those cells that actually belong to the present
    // process. In other words, it seems as if we can't compute the average
    // stresses for all cells. However, remember that our class derived from
    // <code>DataOut</code> only iterates over those cells that actually do
    // belong to the present processor, i.e. we don't have to compute anything
    // for all the other cells as this information would not be touched. The
    // following little loop does this. We enclose the entire block into a
    // pair of braces to make sure that the iterator variables do not remain
    // accidentally visible beyond the end of the block in which they are
    // used:
    Vector<double> norm_of_stress (triangulation.n_active_cells());
    {
      // Loop over all the cells...
      typename Triangulation<dim>::active_cell_iterator
      cell = triangulation.begin_active(),
      endc = triangulation.end();
      for (unsigned int index=0; cell!=endc; ++cell, ++index)
        // ... and pick those that are relevant to us:
        if (cell->subdomain_id() == this_mpi_process)
          {
            // On these cells, add up the stresses over all quadrature
            // points...
            SymmetricTensor<2,dim> accumulated_stress;
            for (unsigned int q=0;
                 q<quadrature_formula.size();
                 ++q)
              accumulated_stress +=
                reinterpret_cast<PointHistory<dim>*>(cell->user_pointer())[q]
                .old_stress;

            // ...then write the norm of the average to their destination:
            norm_of_stress(index)
              = (accumulated_stress /
                 quadrature_formula.size()).norm();
          }
      // And on the cells that we are not interested in, set the respective
      // value in the vector to a bogus value (norms must be positive, and a
      // large negative value should catch your eye) in order to make sure
      // that if we were somehow wrong about our assumption that these
      // elements would not appear in the output file, that we would find out
      // by looking at the graphical output:
        else
          norm_of_stress(index) = -1e+20;
    }
    // Finally attach this vector as well to be treated for output:
    data_out.add_data_vector (norm_of_stress, "norm_of_stress");

    // As a last piece of data, let us also add the partitioning of the domain
    // into subdomains associated with the processors if this is a parallel
    // job. This works in the exact same way as in the step-17 program:
    std::vector<types::subdomain_id> partition_int (triangulation.n_active_cells());
    GridTools::get_subdomain_association (triangulation, partition_int);
    const Vector<double> partitioning(partition_int.begin(),
                                      partition_int.end());
    data_out.add_data_vector (partitioning, "partitioning");

    // Finally, with all this data, we can instruct deal.II to munge the
    // information and produce some intermediate data structures that contain
    // all these solution and other data vectors:
    data_out.build_patches ();


    // Let us determine the name of the file we will want to write it to. We
    // compose it of the prefix <code>solution-</code>, followed by the time
    // step number, and finally the processor id (encoded as a three digit
    // number):
    std::string filename = "solution-" + Utilities::int_to_string(timestep_no,4)
                           + "." + Utilities::int_to_string(this_mpi_process,3)
                           + ".vtu";

    // The following assertion makes sure that there are less than 1000
    // processes (a very conservative check, but worth having anyway) as our
    // scheme of generating process numbers would overflow if there were 1000
    // processes or more. Note that we choose to use <code>AssertThrow</code>
    // rather than <code>Assert</code> since the number of processes is a
    // variable that depends on input files or the way the process is started,
    // rather than static assumptions in the program code. Therefore, it is
    // inappropriate to use <code>Assert</code> that is optimized away in
    // optimized mode, whereas here we actually can assume that users will run
    // the largest computations with the most processors in optimized mode,
    // and we should check our assumptions in this particular case, and not
    // only when running in debug mode:
    AssertThrow (n_mpi_processes < 1000, ExcNotImplemented());

    // With the so-completed filename, let us open a file and write the data
    // we have generated into it:
    std::ofstream output (filename.c_str());
    data_out.write_vtu (output);

    // The record files must be written only once and not by each processor,
    // so we do this on processor 0:
    if (this_mpi_process==0)
      {
        // Here we collect all filenames of the current timestep (same format as above)
        std::vector<std::string> filenames;
        for (unsigned int i=0; i<n_mpi_processes; ++i)
          filenames.push_back ("solution-" + Utilities::int_to_string(timestep_no,4)
                               + "." + Utilities::int_to_string(i,3)
                               + ".vtu");

        // Now we write the .visit file. The naming is similar to the .vtu files, only
        // that the file obviously doesn't contain a processor id.
        const std::string
        visit_master_filename = ("solution-" +
                                 Utilities::int_to_string(timestep_no,4) +
                                 ".visit");
        std::ofstream visit_master (visit_master_filename.c_str());
        data_out.write_visit_record (visit_master, filenames);

        // Similarly, we write the paraview .pvtu:
        const std::string
        pvtu_master_filename = ("solution-" +
                                Utilities::int_to_string(timestep_no,4) +
                                ".pvtu");
        std::ofstream pvtu_master (pvtu_master_filename.c_str());
        data_out.write_pvtu_record (pvtu_master, filenames);

        // Finally, we write the paraview record, that references all .pvtu files and
        // their respective time. Note that the variable times_and_names is declared
        // static, so it will retain the entries from the pervious timesteps.
        static std::vector<std::pair<double,std::string> > times_and_names;
        times_and_names.push_back (std::pair<double,std::string> (present_time, pvtu_master_filename));
        std::ofstream pvd_output ("solution.pvd");
        data_out.write_pvd_record (pvd_output, times_and_names);
      }

  }



  // @sect4{TopLevel::do_initial_timestep}

  // This and the next function handle the overall structure of the first and
  // following timesteps, respectively. The first timestep is slightly more
  // involved because we want to compute it multiple times on successively
  // refined meshes, each time starting from a clean state. At the end of
  // these computations, in which we compute the incremental displacements
  // each time, we use the last results obtained for the incremental
  // displacements to compute the resulting stress updates and move the mesh
  // accordingly. On this new mesh, we then output the solution and any
  // additional data we consider important.
  //
  // All this is interspersed by generating output to the console to update
  // the person watching the screen on what is going on. As in step-17, the
  // use of <code>pcout</code> instead of <code>std::cout</code> makes sure
  // that only one of the parallel processes is actually writing to the
  // console, without having to explicitly code an if-statement in each place
  // where we generate output:
  template <int dim>
  void TopLevel<dim>::do_initial_timestep ()
  {
    present_time += present_timestep;
    ++timestep_no;
    pcout << "Timestep " << timestep_no << " at time " << present_time
          << std::endl;

    for (unsigned int cycle=0; cycle<2; ++cycle)
      {
        pcout << "  Cycle " << cycle << ':' << std::endl;

        if (cycle == 0)
          create_coarse_grid ();
        else
          refine_initial_grid ();

        pcout << "    Number of active cells:       "
              << triangulation.n_active_cells()
              << " (by partition:";
        for (unsigned int p=0; p<n_mpi_processes; ++p)
          pcout << (p==0 ? ' ' : '+')
                << (GridTools::
                    count_cells_with_subdomain_association (triangulation,p));
        pcout << ")" << std::endl;

        setup_system ();

        pcout << "    Number of degrees of freedom: "
              << dof_handler.n_dofs()
              << " (by partition:";
        for (unsigned int p=0; p<n_mpi_processes; ++p)
          pcout << (p==0 ? ' ' : '+')
                << (DoFTools::
                    count_dofs_with_subdomain_association (dof_handler,p));
        pcout << ")" << std::endl;

        solve_timestep ();
      }

    move_mesh ();
    output_results ();

    pcout << std::endl;
  }



  // @sect4{TopLevel::do_timestep}

  // Subsequent timesteps are simpler, and probably do not require any more
  // documentation given the explanations for the previous function above:
  template <int dim>
  void TopLevel<dim>::do_timestep ()
  {
    present_time += present_timestep;
    ++timestep_no;
    pcout << "Timestep " << timestep_no << " at time " << present_time
          << std::endl;
    if (present_time > end_time)
      {
        present_timestep -= (present_time - end_time);
        present_time = end_time;
      }


    solve_timestep ();

    move_mesh ();
    output_results ();

    pcout << std::endl;
  }


  // @sect4{TopLevel::refine_initial_grid}

  // The following function is called when solving the first time step on
  // successively refined meshes. After each iteration, it computes a
  // refinement criterion, refines the mesh, and sets up the history variables
  // in each quadrature point again to a clean state.
  template <int dim>
  void TopLevel<dim>::refine_initial_grid ()
  {
    // First, let each process compute error indicators for the cells it owns:
    Vector<float> error_per_cell (triangulation.n_active_cells());
    KellyErrorEstimator<dim>::estimate (dof_handler,
                                        QGauss<dim-1>(2),
                                        typename FunctionMap<dim>::type(),
                                        incremental_displacement,
                                        error_per_cell,
                                        ComponentMask(),
                                        0,
                                        multithread_info.n_threads(),
                                        this_mpi_process);

    // Then set up a global vector into which we merge the local indicators
    // from each of the %parallel processes:
    const unsigned int n_local_cells
      = GridTools::count_cells_with_subdomain_association (triangulation,
                                                           this_mpi_process);
    PETScWrappers::MPI::Vector
    distributed_error_per_cell (mpi_communicator,
                                triangulation.n_active_cells(),
                                n_local_cells);

    for (unsigned int i=0; i<error_per_cell.size(); ++i)
      if (error_per_cell(i) != 0)
        distributed_error_per_cell(i) = error_per_cell(i);
    distributed_error_per_cell.compress (VectorOperation::insert);

    // Once we have that, copy it back into local copies on all processors and
    // refine the mesh accordingly:
    error_per_cell = distributed_error_per_cell;
    GridRefinement::refine_and_coarsen_fixed_number (triangulation,
                                                     error_per_cell,
                                                     0.35, 0.03);
    triangulation.execute_coarsening_and_refinement ();

    // Finally, set up quadrature point data again on the new mesh, and only
    // on those cells that we have determined to be ours:
    GridTools::partition_triangulation (n_mpi_processes, triangulation);
    setup_quadrature_point_history ();
  }



  // @sect4{TopLevel::move_mesh}

  // At the end of each time step, we move the nodes of the mesh according to
  // the incremental displacements computed in this time step. To do this, we
  // keep a vector of flags that indicate for each vertex whether we have
  // already moved it around, and then loop over all cells and move those
  // vertices of the cell that have not been moved yet. It is worth noting
  // that it does not matter from which of the cells adjacent to a vertex we
  // move this vertex: since we compute the displacement using a continuous
  // finite element, the displacement field is continuous as well and we can
  // compute the displacement of a given vertex from each of the adjacent
  // cells. We only have to make sure that we move each node exactly once,
  // which is why we keep the vector of flags.
  //
  // There are two noteworthy things in this function. First, how we get the
  // displacement field at a given vertex using the
  // <code>cell-@>vertex_dof_index(v,d)</code> function that returns the index
  // of the <code>d</code>th degree of freedom at vertex <code>v</code> of the
  // given cell. In the present case, displacement in the k-th coordinate
  // direction corresponds to the k-th component of the finite element. Using a
  // function like this bears a certain risk, because it uses knowledge of the
  // order of elements that we have taken together for this program in the
  // <code>FESystem</code> element. If we decided to add an additional
  // variable, for example a pressure variable for stabilization, and happened
  // to insert it as the first variable of the element, then the computation
  // below will start to produce nonsensical results. In addition, this
  // computation rests on other assumptions: first, that the element we use
  // has, indeed, degrees of freedom that are associated with vertices. This
  // is indeed the case for the present Q1 element, as would be for all Qp
  // elements of polynomial order <code>p</code>. However, it would not hold
  // for discontinuous elements, or elements for mixed formulations. Secondly,
  // it also rests on the assumption that the displacement at a vertex is
  // determined solely by the value of the degree of freedom associated with
  // this vertex; in other words, all shape functions corresponding to other
  // degrees of freedom are zero at this particular vertex. Again, this is the
  // case for the present element, but is not so for all elements that are
  // presently available in deal.II. Despite its risks, we choose to use this
  // way in order to present a way to query individual degrees of freedom
  // associated with vertices.
  //
  // In this context, it is instructive to point out what a more general way
  // would be. For general finite elements, the way to go would be to take a
  // quadrature formula with the quadrature points in the vertices of a
  // cell. The <code>QTrapez</code> formula for the trapezoidal rule does
  // exactly this. With this quadrature formula, we would then initialize an
  // <code>FEValues</code> object in each cell, and use the
  // <code>FEValues::get_function_values</code> function to obtain the values
  // of the solution function in the quadrature points, i.e. the vertices of
  // the cell. These are the only values that we really need, i.e. we are not
  // at all interested in the weights (or the <code>JxW</code> values)
  // associated with this particular quadrature formula, and this can be
  // specified as the last argument in the constructor to
  // <code>FEValues</code>. The only point of minor inconvenience in this
  // scheme is that we have to figure out which quadrature point corresponds
  // to the vertex we consider at present, as they may or may not be ordered
  // in the same order.
  //
  // This inconvenience could be avoided if finite elements have support
  // points on vertices (which the one here has; for the concept of support
  // points, see @ref GlossSupport "support points"). For such a case, one
  // could construct a custom quadrature rule using
  // FiniteElement::get_unit_support_points(). The first
  // <code>GeometryInfo@<dim@>::vertices_per_cell*fe.dofs_per_vertex</code>
  // quadrature points will then correspond to the vertices of the cell and
  // are ordered consistent with <code>cell-@>vertex(i)</code>, taking into
  // account that support points for vector elements will be duplicated
  // <code>fe.dofs_per_vertex</code> times.
  //
  // Another point worth explaining about this short function is the way in
  // which the triangulation class exports information about its vertices:
  // through the <code>Triangulation::n_vertices</code> function, it
  // advertises how many vertices there are in the triangulation. Not all of
  // them are actually in use all the time -- some are left-overs from cells
  // that have been coarsened previously and remain in existence since deal.II
  // never changes the number of a vertex once it has come into existence,
  // even if vertices with lower number go away. Secondly, the location
  // returned by <code>cell-@>vertex(v)</code> is not only a read-only object
  // of type <code>Point@<dim@></code>, but in fact a reference that can be
  // written to. This allows to move around the nodes of a mesh with relative
  // ease, but it is worth pointing out that it is the responsibility of an
  // application program using this feature to make sure that the resulting
  // cells are still useful, i.e. are not distorted so much that the cell is
  // degenerated (indicated, for example, by negative Jacobians). Note that we
  // do not have any provisions in this function to actually ensure this, we
  // just have faith.
  //
  // After this lengthy introduction, here are the full 20 or so lines of
  // code:
  template <int dim>
  void TopLevel<dim>::move_mesh ()
  {
    pcout << "    Moving mesh..." << std::endl;

    std::vector<bool> vertex_touched (triangulation.n_vertices(),
                                      false);
    for (typename DoFHandler<dim>::active_cell_iterator
         cell = dof_handler.begin_active ();
         cell != dof_handler.end(); ++cell)
      for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
        if (vertex_touched[cell->vertex_index(v)] == false)
          {
            vertex_touched[cell->vertex_index(v)] = true;

            Point<dim> vertex_displacement;
            for (unsigned int d=0; d<dim; ++d)
              vertex_displacement[d]
                = incremental_displacement(cell->vertex_dof_index(v,d));

            cell->vertex(v) += vertex_displacement;
          }
  }


  // @sect4{TopLevel::setup_quadrature_point_history}

  // At the beginning of our computations, we needed to set up initial values
  // of the history variables, such as the existing stresses in the material,
  // that we store in each quadrature point. As mentioned above, we use the
  // <code>user_pointer</code> for this that is available in each cell.
  //
  // To put this into larger perspective, we note that if we had previously
  // available stresses in our model (which we assume do not exist for the
  // purpose of this program), then we would need to interpolate the field of
  // preexisting stresses to the quadrature points. Likewise, if we were to
  // simulate elasto-plastic materials with hardening/softening, then we would
  // have to store additional history variables like the present yield stress
  // of the accumulated plastic strains in each quadrature
  // points. Pre-existing hardening or weakening would then be implemented by
  // interpolating these variables in the present function as well.
  template <int dim>
  void TopLevel<dim>::setup_quadrature_point_history ()
  {
    // What we need to do here is to first count how many quadrature points
    // are within the responsibility of this processor. This, of course,
    // equals the number of cells that belong to this processor times the
    // number of quadrature points our quadrature formula has on each cell.
    //
    // For good measure, we also set all user pointers of all cells, whether
    // ours of not, to the null pointer. This way, if we ever access the user
    // pointer of a cell which we should not have accessed, a segmentation
    // fault will let us know that this should not have happened:
    unsigned int our_cells = 0;
    for (typename Triangulation<dim>::active_cell_iterator
         cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell)
      if (cell->subdomain_id() == this_mpi_process)
        ++our_cells;

    triangulation.clear_user_data();

    // Next, allocate as many quadrature objects as we need. Since the
    // <code>resize</code> function does not actually shrink the amount of
    // allocated memory if the requested new size is smaller than the old
    // size, we resort to a trick to first free all memory, and then
    // reallocate it: we declare an empty vector as a temporary variable and
    // then swap the contents of the old vector and this temporary
    // variable. This makes sure that the
    // <code>quadrature_point_history</code> is now really empty, and we can
    // let the temporary variable that now holds the previous contents of the
    // vector go out of scope and be destroyed. In the next step. we can then
    // re-allocate as many elements as we need, with the vector
    // default-initializing the <code>PointHistory</code> objects, which
    // includes setting the stress variables to zero.
    {
      std::vector<PointHistory<dim> > tmp;
      tmp.swap (quadrature_point_history);
    }
    quadrature_point_history.resize (our_cells *
                                     quadrature_formula.size());

    // Finally loop over all cells again and set the user pointers from the
    // cells that belong to the present processor to point to the first
    // quadrature point objects corresponding to this cell in the vector of
    // such objects:
    unsigned int history_index = 0;
    for (typename Triangulation<dim>::active_cell_iterator
         cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell)
      if (cell->subdomain_id() == this_mpi_process)
        {
          cell->set_user_pointer (&quadrature_point_history[history_index]);
          history_index += quadrature_formula.size();
        }

    // At the end, for good measure make sure that our count of elements was
    // correct and that we have both used up all objects we allocated
    // previously, and not point to any objects beyond the end of the
    // vector. Such defensive programming strategies are always good checks to
    // avoid accidental errors and to guard against future changes to this
    // function that forget to update all uses of a variable at the same
    // time. Recall that constructs using the <code>Assert</code> macro are
    // optimized away in optimized mode, so do not affect the run time of
    // optimized runs:
    Assert (history_index == quadrature_point_history.size(),
            ExcInternalError());
  }




  // @sect4{TopLevel::update_quadrature_point_history}

  // At the end of each time step, we should have computed an incremental
  // displacement update so that the material in its new configuration
  // accommodates for the difference between the external body and boundary
  // forces applied during this time step minus the forces exerted through
  // preexisting internal stresses. In order to have the preexisting
  // stresses available at the next time step, we therefore have to update the
  // preexisting stresses with the stresses due to the incremental
  // displacement computed during the present time step. Ideally, the
  // resulting sum of internal stresses would exactly counter all external
  // forces. Indeed, a simple experiment can make sure that this is so: if we
  // choose boundary conditions and body forces to be time independent, then
  // the forcing terms (the sum of external forces and internal stresses)
  // should be exactly zero. If you make this experiment, you will realize
  // from the output of the norm of the right hand side in each time step that
  // this is almost the case: it is not exactly zero, since in the first time
  // step the incremental displacement and stress updates were computed
  // relative to the undeformed mesh, which was then deformed. In the second
  // time step, we again compute displacement and stress updates, but this
  // time in the deformed mesh -- there, the resulting updates are very small
  // but not quite zero. This can be iterated, and in each such iteration the
  // residual, i.e. the norm of the right hand side vector, is reduced; if one
  // makes this little experiment, one realizes that the norm of this residual
  // decays exponentially with the number of iterations, and after an initial
  // very rapid decline is reduced by roughly a factor of about 3.5 in each
  // iteration (for one testcase I looked at, other testcases, and other
  // numbers of unknowns change the factor, but not the exponential decay).

  // In a sense, this can then be considered as a quasi-timestepping scheme to
  // resolve the nonlinear problem of solving large-deformation elasticity on
  // a mesh that is moved along in a Lagrangian manner.
  //
  // Another complication is that the existing (old) stresses are defined on
  // the old mesh, which we will move around after updating the stresses. If
  // this mesh update involves rotations of the cell, then we need to also
  // rotate the updated stress, since it was computed relative to the
  // coordinate system of the old cell.
  //
  // Thus, what we need is the following: on each cell which the present
  // processor owns, we need to extract the old stress from the data stored
  // with each quadrature point, compute the stress update, add the two
  // together, and then rotate the result together with the incremental
  // rotation computed from the incremental displacement at the present
  // quadrature point. We will detail these steps below:
  template <int dim>
  void TopLevel<dim>::update_quadrature_point_history ()
  {
    // First, set up an <code>FEValues</code> object by which we will evaluate
    // the incremental displacements and the gradients thereof at the
    // quadrature points, together with a vector that will hold this
    // information:
    FEValues<dim> fe_values (fe, quadrature_formula,
                             update_values | update_gradients);
    std::vector<std::vector<Tensor<1,dim> > >
    displacement_increment_grads (quadrature_formula.size(),
                                  std::vector<Tensor<1,dim> >(dim));

    // Then loop over all cells and do the job in the cells that belong to our
    // subdomain:
    for (typename DoFHandler<dim>::active_cell_iterator
         cell = dof_handler.begin_active();
         cell != dof_handler.end(); ++cell)
      if (cell->subdomain_id() == this_mpi_process)
        {
          // Next, get a pointer to the quadrature point history data local to
          // the present cell, and, as a defensive measure, make sure that
          // this pointer is within the bounds of the global array:
          PointHistory<dim> *local_quadrature_points_history
            = reinterpret_cast<PointHistory<dim> *>(cell->user_pointer());
          Assert (local_quadrature_points_history >=
                  &quadrature_point_history.front(),
                  ExcInternalError());
          Assert (local_quadrature_points_history <
                  &quadrature_point_history.back(),
                  ExcInternalError());

          // Then initialize the <code>FEValues</code> object on the present
          // cell, and extract the gradients of the displacement at the
          // quadrature points for later computation of the strains
          fe_values.reinit (cell);
          fe_values.get_function_grads (incremental_displacement,
                                        displacement_increment_grads);

          // Then loop over the quadrature points of this cell:
          for (unsigned int q=0; q<quadrature_formula.size(); ++q)
            {
              // On each quadrature point, compute the strain increment from
              // the gradients, and multiply it by the stress-strain tensor to
              // get the stress update. Then add this update to the already
              // existing strain at this point:
              const SymmetricTensor<2,dim> new_stress
                = (local_quadrature_points_history[q].old_stress
                   +
                   (stress_strain_tensor *
                    get_strain (displacement_increment_grads[q])));

              // Finally, we have to rotate the result. For this, we first
              // have to compute a rotation matrix at the present quadrature
              // point from the incremental displacements. In fact, it can be
              // computed from the gradients, and we already have a function
              // for that purpose:
              const Tensor<2,dim> rotation
                = get_rotation_matrix (displacement_increment_grads[q]);
              // Note that the result, a rotation matrix, is in general an
              // antisymmetric tensor of rank 2, so we must store it as a full
              // tensor.

              // With this rotation matrix, we can compute the rotated tensor
              // by contraction from the left and right, after we expand the
              // symmetric tensor <code>new_stress</code> into a full tensor:
              const SymmetricTensor<2,dim> rotated_new_stress
                = symmetrize(transpose(rotation) *
                             static_cast<Tensor<2,dim> >(new_stress) *
                             rotation);
              // Note that while the result of the multiplication of these
              // three matrices should be symmetric, it is not due to floating
              // point round off: we get an asymmetry on the order of 1e-16 of
              // the off-diagonal elements of the result. When assigning the
              // result to a <code>SymmetricTensor</code>, the constructor of
              // that class checks the symmetry and realizes that it isn't
              // exactly symmetric; it will then raise an exception. To avoid
              // that, we explicitly symmetrize the result to make it exactly
              // symmetric.

              // The result of all these operations is then written back into
              // the original place:
              local_quadrature_points_history[q].old_stress
                = rotated_new_stress;
            }
        }
  }

  // This ends the project specific namespace <code>Step18</code>. The rest is
  // as usual and as already shown in step-17: A <code>main()</code> function
  // that initializes and terminates PETSc, calls the classes that do the
  // actual work, and makes sure that we catch all exceptions that propagate
  // up to this point:
}


int main (int argc, char **argv)
{
  try
    {
      using namespace dealii;
      using namespace Step18;

      Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
      {
        deallog.depth_console (0);

        TopLevel<3> elastic_problem;
        elastic_problem.run ();
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
