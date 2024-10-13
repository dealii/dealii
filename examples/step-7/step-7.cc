/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2000 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */


// @sect3{Include files}

// These first include files have all been treated in previous examples, so we
// won't explain what is in them again.
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/data_out.h>

// In this example, we will not use the numeration scheme which is used per
// default by the DoFHandler class, but will renumber them using the
// Cuthill-McKee algorithm. As has already been explained in step-2, the
// necessary functions are declared in the following file:
#include <deal.II/dofs/dof_renumbering.h>
// Then we will show a little trick how we can make sure that objects are not
// deleted while they are still in use. For this purpose, deal.II has the
// ObserverPointer helper class, which is declared in this file:
#include <deal.II/base/observer_pointer.h>
// Next, we will want to use the function VectorTools::integrate_difference()
// mentioned in the introduction, and we are going to use a ConvergenceTable
// that collects all important data during a run and prints it at the end as a
// table. These comes from the following two files:
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/base/convergence_table.h>
// And finally, we need to use the FEFaceValues class, which is declared in
// the same file as the FEValues class:
#include <deal.II/fe/fe_values.h>

#include <array>
#include <fstream>
#include <iostream>

// The last step before we go on with the actual implementation is to open a
// namespace <code>Step7</code> into which we will put everything, as
// discussed at the end of the introduction, and to import the members of
// namespace <code>dealii</code> into it:
namespace Step7
{
  using namespace dealii;

  // @sect3{Equation data}

  // Before implementing the classes that actually solve something, we first
  // declare and define some function classes that represent right hand side
  // and solution classes. Since we want to compare the numerically obtained
  // solution to the exact continuous one, we need a function object that
  // represents the continuous solution. On the other hand, we need the right
  // hand side function, and that one of course shares some characteristics
  // with the solution. In order to reduce dependencies which arise if we have
  // to change something in both classes at the same time, we move the common
  // characteristics of both functions into a base class.
  //
  // The common characteristics for solution (as explained in the
  // introduction, we choose a sum of three exponentials) and right hand side,
  // are these: the number of exponentials, their centers, and their half
  // width. We declare them in the following class. Since the number of
  // exponentials is a compile-time constant we use a fixed-length
  // <code>std::array</code> to store the center points:
  template <int dim>
  class SolutionBase
  {
  protected:
    static const std::array<Point<dim>, 3> source_centers;
    static const double                    width;
  };


  // The variables which denote the centers and the width of the exponentials
  // have just been declared, now we still need to assign values to
  // them. Here, we can show another small piece of template sorcery, namely
  // how we can assign different values to these variables depending on the
  // dimension. We will only use the 2d case in the program, but we show the
  // 1d case for exposition of a useful technique.
  //
  // First we assign values to the centers for the 1d case, where we place the
  // centers equidistantly at -1/3, 0, and 1/3. The <code>template
  // &lt;&gt;</code> header for this definition indicates an explicit
  // specialization. This means, that the variable belongs to a template, but
  // that instead of providing the compiler with a template from which it can
  // specialize a concrete variable by substituting <code>dim</code> with some
  // concrete value, we provide a specialization ourselves, in this case for
  // <code>dim=1</code>. If the compiler then sees a reference to this
  // variable in a place where the template argument equals one, it knows that
  // it doesn't have to generate the variable from a template by substituting
  // <code>dim</code>, but can immediately use the following definition:
  template <>
  const std::array<Point<1>, 3> SolutionBase<1>::source_centers = {
    {Point<1>(-1.0 / 3.0), Point<1>(0.0), Point<1>(+1.0 / 3.0)}};

  // Likewise, we can provide an explicit specialization for
  // <code>dim=2</code>. We place the centers for the 2d case as follows:
  template <>
  const std::array<Point<2>, 3> SolutionBase<2>::source_centers = {
    {Point<2>(-0.5, +0.5), Point<2>(-0.5, -0.5), Point<2>(+0.5, -0.5)}};

  // There remains to assign a value to the half-width of the exponentials. We
  // would like to use the same value for all dimensions. In this case, we
  // simply provide the compiler with a template from which it can generate a
  // concrete instantiation by substituting <code>dim</code> with a concrete
  // value:
  template <int dim>
  const double SolutionBase<dim>::width = 1. / 8.;



  // After declaring and defining the characteristics of solution and right
  // hand side, we can declare the classes representing these two. They both
  // represent continuous functions, so they are derived from the
  // Function&lt;dim&gt; base class, and they also inherit the characteristics
  // defined in the SolutionBase class.
  //
  // The actual classes are declared in the following. Note that in order to
  // compute the error of the numerical solution against the continuous one in
  // the L2 and H1 (semi-)norms, we have to provide value and gradient of the
  // exact solution. This is more than we have done in previous examples, where
  // all we provided was the value at one or a list of points. Fortunately, the
  // Function class also has virtual functions for the gradient, so we can
  // simply overload the respective virtual member functions in the Function
  // base class. Note that the gradient of a function in <code>dim</code>
  // space dimensions is a vector of size <code>dim</code>, i.e. a tensor of
  // rank 1 and dimension <code>dim</code>. As for so many other things, the
  // library provides a suitable class for this. One new thing about this
  // class is that it explicitly uses the Tensor objects, which previously
  // appeared as intermediate terms in step-3 and step-4. A tensor is a
  // generalization of scalars (rank zero tensors), vectors (rank one
  // tensors), and matrices (rank two tensors), as well as higher dimensional
  // objects. The Tensor class requires two template arguments: the tensor
  // rank and tensor dimension. For example, here we use tensors of rank one
  // (vectors) with dimension <code>dim</code> (so they have <code>dim</code>
  // entries.) While this is a bit less flexible than using Vector, the
  // compiler can generate faster code when the length of the vector is known
  // at compile time. Additionally, specifying a Tensor of rank one and
  // dimension <code>dim</code> guarantees that the tensor will have the right
  // shape (since it is built into the type of the object itself), so the
  // compiler can catch most size-related mistakes for us.
  template <int dim>
  class Solution : public Function<dim>, protected SolutionBase<dim>
  {
  public:
    virtual double value(const Point<dim>  &p,
                         const unsigned int component = 0) const override;

    virtual Tensor<1, dim>
    gradient(const Point<dim>  &p,
             const unsigned int component = 0) const override;
  };


  // The actual definition of the values and gradients of the exact solution
  // class is according to their mathematical definition and does not need
  // much explanation.
  //
  // The only thing that is worth mentioning is that if we access
  // elements of a base class that is template dependent (in this case
  // the elements of SolutionBase&lt;dim&gt;), then the C++ language
  // forces us to write <code>this-&gt;source_centers</code>, and
  // similarly for other members of the base class. C++ does not
  // require the <code>this-&gt;</code> qualification if the base
  // class is not template dependent. The reason why this is necessary
  // is complicated; C++ books will explain under the phrase
  // <i>two-stage (name) lookup</i>, and there is also a lengthy
  // description in the deal.II FAQs.
  template <int dim>
  double Solution<dim>::value(const Point<dim> &p, const unsigned int) const
  {
    double return_value = 0;
    for (const auto &center : this->source_centers)
      {
        const Tensor<1, dim> x_minus_xi = p - center;
        return_value +=
          std::exp(-x_minus_xi.norm_square() / (this->width * this->width));
      }

    return return_value;
  }


  // Likewise, this is the computation of the gradient of the solution.  In
  // order to accumulate the gradient from the contributions of the
  // exponentials, we allocate an object <code>return_value</code> that
  // denotes the mathematical quantity of a tensor of rank <code>1</code> and
  // dimension <code>dim</code>. Its default constructor sets it to the vector
  // containing only zeroes, so we need not explicitly care for its
  // initialization.
  //
  // Note that we could as well have taken the type of the object to be
  // Point&lt;dim&gt; instead of Tensor&lt;1,dim&gt;. Tensors of rank 1 and
  // points are almost exchangeable, and have only very slightly different
  // mathematical meanings. In fact, the Point&lt;dim&gt; class is derived
  // from the Tensor&lt;1,dim&gt; class, which makes up for their mutual
  // exchange ability. Their main difference is in what they logically mean:
  // points are points in space, such as the location at which we want to
  // evaluate a function (see the type of the first argument of this function
  // for example). On the other hand, tensors of rank 1 share the same
  // transformation properties, for example that they need to be rotated in a
  // certain way when we change the coordinate system; however, they do not
  // share the same connotation that points have and are only objects in a
  // more abstract space than the one spanned by the coordinate
  // directions. (In fact, gradients live in `reciprocal' space, since the
  // dimension of their components is not that of a length, but of one over
  // length).
  template <int dim>
  Tensor<1, dim> Solution<dim>::gradient(const Point<dim> &p,
                                         const unsigned int) const
  {
    Tensor<1, dim> return_value;

    for (const auto &center : this->source_centers)
      {
        const Tensor<1, dim> x_minus_xi = p - center;

        // For the gradient, note that its direction is along (x-x_i), so we
        // add up multiples of this distance vector, where the factor is given
        // by the exponentials.
        return_value +=
          (-2. / (this->width * this->width) *
           std::exp(-x_minus_xi.norm_square() / (this->width * this->width)) *
           x_minus_xi);
      }

    return return_value;
  }



  // Besides the function that represents the exact solution, we also need a
  // function which we can use as right hand side when assembling the linear
  // system of discretized equations. This is accomplished using the following
  // class and the following definition of its function. Note that here we
  // only need the value of the function, not its gradients or higher
  // derivatives.
  template <int dim>
  class RightHandSide : public Function<dim>, protected SolutionBase<dim>
  {
  public:
    virtual double value(const Point<dim>  &p,
                         const unsigned int component = 0) const override;
  };


  // The value of the right hand side is given by the negative Laplacian of
  // the solution plus the solution itself, since we wanted to solve
  // Helmholtz's equation:
  template <int dim>
  double RightHandSide<dim>::value(const Point<dim> &p,
                                   const unsigned int) const
  {
    double return_value = 0;
    for (const auto &center : this->source_centers)
      {
        const Tensor<1, dim> x_minus_xi = p - center;

        // The first contribution is the Laplacian:
        return_value +=
          ((2. * dim -
            4. * x_minus_xi.norm_square() / (this->width * this->width)) /
           (this->width * this->width) *
           std::exp(-x_minus_xi.norm_square() / (this->width * this->width)));
        // And the second is the solution itself:
        return_value +=
          std::exp(-x_minus_xi.norm_square() / (this->width * this->width));
      }

    return return_value;
  }


  // @sect3{The Helmholtz solver class}

  // Then we need the class that does all the work. Except for its name, its
  // interface is mostly the same as in previous examples.
  //
  // One of the differences is that we will use this class in several modes:
  // for different finite elements, as well as for adaptive and global
  // refinement. The decision whether global or adaptive refinement shall be
  // used is communicated to the constructor of this class through an
  // enumeration type declared at the top of the class. The constructor then
  // takes a finite element object and the refinement mode as arguments.
  //
  // The rest of the member functions are as before except for the
  // <code>process_solution</code> function: After the solution has been
  // computed, we perform some analysis on it, such as computing the error in
  // various norms. To enable some output, it requires the number of the
  // refinement cycle, and consequently gets it as an argument.
  template <int dim>
  class HelmholtzProblem
  {
  public:
    enum RefinementMode
    {
      global_refinement,
      adaptive_refinement
    };

    HelmholtzProblem(const FiniteElement<dim> &fe,
                     const RefinementMode      refinement_mode);

    void run();

  private:
    void setup_system();
    void assemble_system();
    void solve();
    void refine_grid();
    void process_solution(const unsigned int cycle);

    // Now for the data elements of this class. Among the variables that we
    // have already used in previous examples, only the finite element object
    // differs: The finite elements which the objects of this class operate on
    // are passed to the constructor of this class. It has to store a pointer
    // to the finite element for the member functions to use. Now, for the
    // present class there is no big deal in that, but since we want to show
    // techniques rather than solutions in these programs, we will here point
    // out a problem that often occurs -- and of course the right solution as
    // well.
    //
    // Consider the following situation that occurs in all the example
    // programs: we have a triangulation object, and we have a finite element
    // object, and we also have an object of type DoFHandler that uses both of
    // the first two. These three objects all have a lifetime that is rather
    // long compared to most other objects: they are basically set at the
    // beginning of the program or an outer loop, and they are destroyed at
    // the very end. The question is: can we guarantee that the two objects
    // which the DoFHandler uses, live at least as long as they are in use?
    // This means that the DoFHandler must have some kind of knowledge on the
    // destruction of the other objects.
    //
    // We will show here how the library managed to find out that there are
    // still active references to an object and the object is still alive
    // from the point of view of a using object. Basically, the method is along
    // the following line: all objects that are subject to such potentially
    // dangerous pointers are derived from a class called
    // EnableObserverPointer. For example, the Triangulation,
    // DoFHandler, and a base class of the FiniteElement class are derived from
    // EnableObserverPointer. This latter class does not offer much
    // functionality, but it has a built-in counter which we can subscribe to,
    // thus the name of the class. Whenever we initialize a pointer to that
    // object, we can increase its use counter, and when we move away our
    // pointer or do not need it any more, we decrease the counter again. This
    // way, we can always check how many objects still use that object.
    // Additionally, the class requires to know about a pointer that it can use
    // to tell the subscribing object about its invalidation.
    //
    // If an object of a class that is derived from the
    // EnableObserverPointer class is destroyed, it also has to
    // call the destructor of the EnableObserverPointer class. In
    // this destructor, we tell all the subscribing objects about the
    // invalidation of the object using the stored pointers. The same happens
    // when the object appears on the right hand side of a move expression,
    // i.e., it will no longer contain valid content after the operation. The
    // subscribing class is expected to check the value stored in its
    // corresponding pointer before trying to access the object subscribed to.
    //
    // This is exactly what the ObserverPointer class is doing. It basically
    // acts just like a pointer, i.e. it can be dereferenced, can be assigned to
    // and from other pointers, and so on. On top of that it uses the mechanism
    // described above to find out if the pointer this class is representing is
    // dangling when we try to dereference it. In that case an exception is
    // thrown.
    //
    // In the present example program, we want to protect the finite element
    // object from the situation that for some reason the finite element
    // pointed to is destroyed while still in use. We therefore use a
    // ObserverPointer to the finite element object; since the finite element
    // object is actually never changed in our computations, we pass a const
    // FiniteElement&lt;dim&gt; as template argument to the ObserverPointer
    // class. Note that the pointer so declared is assigned at construction
    // time of the solve object, and destroyed upon destruction, so the lock
    // on the destruction of the finite element object extends throughout the
    // lifetime of this HelmholtzProblem object.
    Triangulation<dim> triangulation;
    DoFHandler<dim>    dof_handler;

    ObserverPointer<const FiniteElement<dim>> fe;

    AffineConstraints<double> hanging_node_constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> system_matrix;

    Vector<double> solution;
    Vector<double> system_rhs;

    // The second to last variable stores the refinement mode passed to the
    // constructor. Since it is only set in the constructor, we can declare
    // this variable constant, to avoid that someone sets it involuntarily
    // (e.g. in an `if'-statement where == was written as = by chance).
    const RefinementMode refinement_mode;

    // For each refinement level some data (like the number of cells, or the
    // L2 error of the numerical solution) will be generated and later
    // printed. The TableHandler can be used to collect all this data and to
    // output it at the end of the run as a table in a simple text or in LaTeX
    // format. Here we don't only use the TableHandler but we use the derived
    // class ConvergenceTable that additionally evaluates rates of
    // convergence:
    ConvergenceTable convergence_table;
  };


  // @sect3{The HelmholtzProblem class implementation}

  // @sect4{HelmholtzProblem::HelmholtzProblem constructor}

  // In the constructor of this class, we only set the variables passed as
  // arguments, and associate the DoF handler object with the triangulation
  // (which is empty at present, however).
  template <int dim>
  HelmholtzProblem<dim>::HelmholtzProblem(const FiniteElement<dim> &fe,
                                          const RefinementMode refinement_mode)
    : dof_handler(triangulation)
    , fe(&fe)
    , refinement_mode(refinement_mode)
  {}


  // @sect4{HelmholtzProblem::setup_system}

  // The following function sets up the degrees of freedom, sizes of matrices
  // and vectors, etc. Most of its functionality has been showed in previous
  // examples, the only difference being the renumbering step immediately
  // after first distributing degrees of freedom.
  //
  // Renumbering the degrees of freedom is not overly difficult, as long as
  // you use one of the algorithms included in the library. It requires only a
  // single line of code. Some more information on this can be found in
  // step-2.
  //
  // Note, however, that when you renumber the degrees of freedom, you must do
  // so immediately after distributing them, since such things as hanging
  // nodes, the sparsity pattern etc. depend on the absolute numbers which are
  // altered by renumbering.
  //
  // The reason why we introduce renumbering here is that it is a relatively
  // cheap operation but often has a beneficial effect: While the CG iteration
  // itself is independent of the actual ordering of degrees of freedom, we
  // will use SSOR as a preconditioner. SSOR goes through all degrees of
  // freedom and does some operations that depend on what happened before; the
  // SSOR operation is therefore not independent of the numbering of degrees
  // of freedom, and it is known that its performance improves by using
  // renumbering techniques. A little experiment shows that indeed, for
  // example, the number of CG iterations for the fifth refinement cycle of
  // adaptive refinement with the Q1 program used here is 40 without, but 36
  // with renumbering. Similar savings can generally be observed for all the
  // computations in this program.
  template <int dim>
  void HelmholtzProblem<dim>::setup_system()
  {
    dof_handler.distribute_dofs(*fe);
    DoFRenumbering::Cuthill_McKee(dof_handler);

    hanging_node_constraints.clear();
    DoFTools::make_hanging_node_constraints(dof_handler,
                                            hanging_node_constraints);
    hanging_node_constraints.close();

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    hanging_node_constraints.condense(dsp);
    sparsity_pattern.copy_from(dsp);

    system_matrix.reinit(sparsity_pattern);

    solution.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());
  }


  // @sect4{HelmholtzProblem::assemble_system}

  // Assembling the system of equations for the problem at hand is mostly as
  // for the example programs before. However, some things have changed
  // anyway, so we comment on this function fairly extensively.
  //
  // At the top of the function you will find the usual assortment of variable
  // declarations. Compared to previous programs, of importance is only that
  // we expect to solve problems also with bi-quadratic elements and therefore
  // have to use sufficiently accurate quadrature formula. In addition, we
  // need to compute integrals over faces, i.e. <code>dim-1</code> dimensional
  // objects. The declaration of a face quadrature formula is then
  // straightforward:
  template <int dim>
  void HelmholtzProblem<dim>::assemble_system()
  {
    const QGauss<dim>     quadrature_formula(fe->degree + 1);
    const QGauss<dim - 1> face_quadrature_formula(fe->degree + 1);

    const unsigned int n_q_points      = quadrature_formula.size();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    const unsigned int dofs_per_cell = fe->n_dofs_per_cell();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    Vector<double>     cell_rhs(dofs_per_cell);

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    // Then we need objects which can evaluate the values, gradients, etc of
    // the shape functions at the quadrature points. While it seems that it
    // should be feasible to do it with one object for both domain and face
    // integrals, there is a subtle difference since the weights in the domain
    // integrals include the measure of the cell in the domain, while the face
    // integral quadrature requires the measure of the face in a
    // lower-dimensional manifold. Internally these two classes are rooted in
    // a common base class which does most of the work and offers the same
    // interface to both domain and interface integrals.
    //
    // For the domain integrals in the bilinear form for Helmholtz's equation,
    // we need to compute the values and gradients, as well as the weights at
    // the quadrature points. Furthermore, we need the quadrature points on
    // the real cell (rather than on the unit cell) to evaluate the right hand
    // side function. The object we use to get at this information is the
    // FEValues class discussed previously.
    //
    // For the face integrals, we only need the values of the shape functions,
    // as well as the weights. We also need the normal vectors and quadrature
    // points on the real cell since we want to determine the Neumann values
    // from the exact solution object (see below). The class that gives us
    // this information is called FEFaceValues:
    FEValues<dim> fe_values(*fe,
                            quadrature_formula,
                            update_values | update_gradients |
                              update_quadrature_points | update_JxW_values);

    FEFaceValues<dim> fe_face_values(*fe,
                                     face_quadrature_formula,
                                     update_values | update_quadrature_points |
                                       update_normal_vectors |
                                       update_JxW_values);

    // Then we need some objects already known from previous examples: An
    // object denoting the right hand side function, its values at the
    // quadrature points on a cell, the cell matrix and right hand side, and
    // the indices of the degrees of freedom on a cell.
    //
    // Note that the operations we will do with the right hand side object are
    // only querying data, never changing the object. We can therefore declare
    // it <code>const</code>:
    RightHandSide<dim>  right_hand_side;
    std::vector<double> rhs_values(n_q_points);

    // Finally we define an object denoting the exact solution function. We
    // will use it to compute the Neumann values at the boundary from
    // it. Usually, one would of course do so using a separate object, in
    // particular since the exact solution is generally unknown while the
    // Neumann values are prescribed. We will, however, be a little bit lazy
    // and use what we already have in information. Real-life programs would
    // to go other ways here, of course.
    Solution<dim> exact_solution;

    // Now for the main loop over all cells. This is mostly unchanged from
    // previous examples, so we only comment on the things that have changed.
    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        fe_values.reinit(cell);

        cell_matrix = 0.;
        cell_rhs    = 0.;

        right_hand_side.value_list(fe_values.get_quadrature_points(),
                                   rhs_values);

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              for (unsigned int j = 0; j < dofs_per_cell; ++j)
                // The first thing that has changed is the bilinear form. It
                // now contains the additional term from the Helmholtz
                // equation:
                cell_matrix(i, j) +=
                  ((fe_values.shape_grad(i, q_point) *     // grad phi_i(x_q)
                      fe_values.shape_grad(j, q_point)     // grad phi_j(x_q)
                    +                                      //
                    fe_values.shape_value(i, q_point) *    // phi_i(x_q)
                      fe_values.shape_value(j, q_point)) * // phi_j(x_q)
                   fe_values.JxW(q_point));                // dx


              cell_rhs(i) += (fe_values.shape_value(i, q_point) * // phi_i(x_q)
                              rhs_values[q_point] *               // f(x_q)
                              fe_values.JxW(q_point));            // dx
            }

        // Then there is that second term on the right hand side, the contour
        // integral. First we have to find out whether the intersection of the
        // faces of this cell with the boundary part Gamma2 is nonzero. To
        // this end, we loop over all faces and check whether its boundary
        // indicator equals <code>1</code>, which is the value that we have
        // assigned to that portions of the boundary composing Gamma2 in the
        // <code>run()</code> function further below. (The default value of
        // boundary indicators is <code>0</code>, so faces can only have an
        // indicator equal to <code>1</code> if we have explicitly set it.)
        for (const auto &face : cell->face_iterators())
          if (face->at_boundary() && (face->boundary_id() == 1))
            {
              // If we came into here, then we have found an external face
              // belonging to Gamma2. Next, we have to compute the values of
              // the shape functions and the other quantities which we will
              // need for the computation of the contour integral. This is
              // done using the <code>reinit</code> function which we already
              // know from the FEValue class:
              fe_face_values.reinit(cell, face);

              // And we can then perform the integration by using a loop over
              // all quadrature points.
              //
              // On each quadrature point, we first compute the value of the
              // normal derivative. We do so using the gradient of the exact
              // solution and the normal vector to the face at the present
              // quadrature point obtained from the
              // <code>fe_face_values</code> object. This is then used to
              // compute the additional contribution of this face to the right
              // hand side:
              for (unsigned int q_point = 0; q_point < n_face_q_points;
                   ++q_point)
                {
                  const double neumann_value =
                    (exact_solution.gradient(
                       fe_face_values.quadrature_point(q_point)) *
                     fe_face_values.normal_vector(q_point));

                  for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    cell_rhs(i) +=
                      (fe_face_values.shape_value(i, q_point) * // phi_i(x_q)
                       neumann_value *                          // g(x_q)
                       fe_face_values.JxW(q_point));            // dx
                }
            }

        // Now that we have the contributions of the present cell, we can
        // transfer it to the global matrix and right hand side vector, as in
        // the examples before:
        cell->get_dof_indices(local_dof_indices);
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            for (unsigned int j = 0; j < dofs_per_cell; ++j)
              system_matrix.add(local_dof_indices[i],
                                local_dof_indices[j],
                                cell_matrix(i, j));

            system_rhs(local_dof_indices[i]) += cell_rhs(i);
          }
      }

    // Likewise, elimination and treatment of boundary values has been shown
    // previously.
    //
    // We note, however that now the boundary indicator for which we
    // interpolate boundary values (denoted by the second parameter to
    // <code>interpolate_boundary_values</code>) does not represent the whole
    // boundary any more. Rather, it is that portion of the boundary which we
    // have not assigned another indicator (see below). The degrees of freedom
    // at the boundary that do not belong to Gamma1 are therefore excluded
    // from the interpolation of boundary values, just as we want.
    hanging_node_constraints.condense(system_matrix);
    hanging_node_constraints.condense(system_rhs);

    std::map<types::global_dof_index, double> boundary_values;
    VectorTools::interpolate_boundary_values(dof_handler,
                                             types::boundary_id(0),
                                             Solution<dim>(),
                                             boundary_values);
    MatrixTools::apply_boundary_values(boundary_values,
                                       system_matrix,
                                       solution,
                                       system_rhs);
  }


  // @sect4{HelmholtzProblem::solve}

  // Solving the system of equations is done in the same way as before:
  template <int dim>
  void HelmholtzProblem<dim>::solve()
  {
    SolverControl            solver_control(1000, 1e-6 * system_rhs.l2_norm());
    SolverCG<Vector<double>> cg(solver_control);

    PreconditionSSOR<SparseMatrix<double>> preconditioner;
    preconditioner.initialize(system_matrix, 1.2);

    cg.solve(system_matrix, solution, system_rhs, preconditioner);

    hanging_node_constraints.distribute(solution);
  }


  // @sect4{HelmholtzProblem::refine_grid}

  // Now for the function doing grid refinement. Depending on the refinement
  // mode passed to the constructor, we do global or adaptive refinement.
  //
  // Global refinement is simple, so there is not much to comment on.  In case
  // of adaptive refinement, we use the same functions and classes as in the
  // previous example program. Note that one could treat Neumann boundaries
  // differently than Dirichlet boundaries, and one should in fact do so here
  // since we have Neumann boundary conditions on part of the boundaries, but
  // since we don't have a function here that describes the Neumann values (we
  // only construct these values from the exact solution when assembling the
  // matrix), we omit this detail even though doing this in a strictly correct
  // way would not be hard to add.
  //
  // At the end of the switch, we have a default case that simply says
  // `DEAL_II_ASSERT_UNREACHABLE()`. This macro raises an error whenever
  // the program reaches this point; the program is then aborted.
  // This is intentional: Right now we have only
  // implemented two refinement strategies (global and adaptive), but someone
  // might want to add a third strategy (for example adaptivity with a
  // different refinement criterion) and add a third member to the enumeration
  // that determines the refinement mode. If it weren't for the default case
  // of the switch statement, this function would simply run to its end
  // without doing anything. This is most likely not what was intended. One of
  // the [defensive programming
  // techniques](https://en.wikipedia.org/wiki/Defensive_programming) that you
  // will find all over the deal.II library is therefore to always have default
  // cases that abort, to make sure that values not considered when listing the
  // cases in the switch statement are eventually caught, and forcing
  // programmers to add code to handle them. The documentation of
  // DEAL_II_ASSERT_UNREACHABLE() shows other examples of how this macro can be
  // used.
  //
  // We will use this same technique in other places further down as well.
  template <int dim>
  void HelmholtzProblem<dim>::refine_grid()
  {
    switch (refinement_mode)
      {
        case global_refinement:
          {
            triangulation.refine_global(1);
            break;
          }

        case adaptive_refinement:
          {
            Vector<float> estimated_error_per_cell(
              triangulation.n_active_cells());

            KellyErrorEstimator<dim>::estimate(
              dof_handler,
              QGauss<dim - 1>(fe->degree + 1),
              std::map<types::boundary_id, const Function<dim> *>(),
              solution,
              estimated_error_per_cell);

            GridRefinement::refine_and_coarsen_fixed_number(
              triangulation, estimated_error_per_cell, 0.3, 0.03);

            triangulation.execute_coarsening_and_refinement();

            break;
          }

        default:
          {
            DEAL_II_ASSERT_UNREACHABLE();
          }
      }
  }


  // @sect4{HelmholtzProblem::process_solution}

  // Finally we want to process the solution after it has been computed. For
  // this, we integrate the error in various (semi-)norms, and we generate
  // tables that will later be used to display the convergence against the
  // continuous solution in a nice format.
  template <int dim>
  void HelmholtzProblem<dim>::process_solution(const unsigned int cycle)
  {
    // Our first task is to compute error norms. In order to integrate the
    // difference between computed numerical solution and the continuous
    // solution (described by the Solution class defined at the top of this
    // file), we first need a vector that will hold the norm of the error on
    // each cell. Since accuracy with 16 digits is not so important for these
    // quantities, we save some memory by using <code>float</code> instead of
    // <code>double</code> values.
    //
    // The next step is to use a function from the library which computes the
    // error in the L2 norm on each cell.  We have to pass it the DoF handler
    // object, the vector holding the nodal values of the numerical solution,
    // the continuous solution as a function object, the vector into which it
    // shall place the norm of the error on each cell, a quadrature rule by
    // which this norm shall be computed, and the type of norm to be
    // used. Here, we use a Gauss formula with three points in each space
    // direction, and compute the L2 norm.
    //
    // Finally, we want to get the global L2 norm. This can of course be
    // obtained by summing the squares of the norms on each cell, and taking
    // the square root of that value. This is equivalent to taking the l2
    // (lower case <code>l</code>) norm of the vector of norms on each cell:
    Vector<float> difference_per_cell(triangulation.n_active_cells());
    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      Solution<dim>(),
                                      difference_per_cell,
                                      QGauss<dim>(fe->degree + 1),
                                      VectorTools::L2_norm);
    const double L2_error =
      VectorTools::compute_global_error(triangulation,
                                        difference_per_cell,
                                        VectorTools::L2_norm);

    // By same procedure we get the H1 semi-norm. We re-use the
    // <code>difference_per_cell</code> vector since it is no longer used
    // after computing the <code>L2_error</code> variable above. The global
    // $H^1$ semi-norm error is then computed by taking the sum of squares
    // of the errors on each individual cell, and then the square root of
    // it -- an operation that is conveniently performed by
    // VectorTools::compute_global_error.
    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      Solution<dim>(),
                                      difference_per_cell,
                                      QGauss<dim>(fe->degree + 1),
                                      VectorTools::H1_seminorm);
    const double H1_error =
      VectorTools::compute_global_error(triangulation,
                                        difference_per_cell,
                                        VectorTools::H1_seminorm);

    // Finally, we compute the maximum norm. Of course, we can't actually
    // compute the true maximum of the error over *all* points in the domain,
    // but only the maximum over a finite set of evaluation points that, for
    // convenience, we will still call "quadrature points" and represent by
    // an object of type Quadrature even though we do not actually perform any
    // integration.
    //
    // There is then the question of what points precisely we want to evaluate
    // at. It turns out that the result we get depends quite sensitively on the
    // "quadrature" points being used. There is also the issue of
    // superconvergence: Finite element solutions are, on some meshes and for
    // polynomial degrees $k\ge 2$, particularly accurate at the node points as
    // well as at Gauss-Lobatto points, much more accurate than at randomly
    // chosen points. (See
    // @cite Li2019 and the discussion and references in Section 1.2 for more
    // information on this.) In other words, if we are interested in finding
    // the largest difference $u(\mathbf x)-u_h(\mathbf x)$, then we ought to
    // look at points $\mathbf x$ that are specifically not of this "special"
    // kind of points and we should specifically not use
    // `QGauss(fe->degree+1)` to define where we evaluate. Rather, we use a
    // special quadrature rule that is obtained by iterating the trapezoidal
    // rule by the degree of the finite element times two plus one in each space
    // direction. Note that the constructor of the QIterated class takes a
    // one-dimensional quadrature rule and a number that tells it how often it
    // shall repeat this rule in each space direction.
    //
    // Using this special quadrature rule, we can then try to find the maximal
    // error on each cell. Finally, we compute the global L infinity error
    // from the L infinity errors on each cell with a call to
    // VectorTools::compute_global_error.
    const QTrapezoid<1>  q_trapez;
    const QIterated<dim> q_iterated(q_trapez, fe->degree * 2 + 1);
    VectorTools::integrate_difference(dof_handler,
                                      solution,
                                      Solution<dim>(),
                                      difference_per_cell,
                                      q_iterated,
                                      VectorTools::Linfty_norm);
    const double Linfty_error =
      VectorTools::compute_global_error(triangulation,
                                        difference_per_cell,
                                        VectorTools::Linfty_norm);

    // After all these errors have been computed, we finally write some
    // output. In addition, we add the important data to the TableHandler by
    // specifying the key of the column and the value.  Note that it is not
    // necessary to define column keys beforehand -- it is sufficient to just
    // add values, and columns will be introduced into the table in the order
    // values are added the first time.
    const unsigned int n_active_cells = triangulation.n_active_cells();
    const unsigned int n_dofs         = dof_handler.n_dofs();

    std::cout << "Cycle " << cycle << ':' << std::endl
              << "   Number of active cells:       " << n_active_cells
              << std::endl
              << "   Number of degrees of freedom: " << n_dofs << std::endl;

    convergence_table.add_value("cycle", cycle);
    convergence_table.add_value("cells", n_active_cells);
    convergence_table.add_value("dofs", n_dofs);
    convergence_table.add_value("L2", L2_error);
    convergence_table.add_value("H1", H1_error);
    convergence_table.add_value("Linfty", Linfty_error);
  }


  // @sect4{HelmholtzProblem::run}

  // As in previous example programs, the <code>run</code> function controls
  // the flow of execution. The basic layout is as in previous examples: an
  // outer loop over successively refined grids, and in this loop first
  // problem setup, assembling the linear system, solution, and
  // post-processing.
  //
  // The first task in the main loop is creation and refinement of grids. This
  // is as in previous examples, with the only difference that we want to have
  // part of the boundary marked as Neumann type, rather than Dirichlet.
  //
  // For this, we will use the following convention: Faces belonging to Gamma1
  // will have the boundary indicator <code>0</code> (which is the default, so
  // we don't have to set it explicitly), and faces belonging to Gamma2 will
  // use <code>1</code> as boundary indicator.  To set these values, we loop
  // over all cells, then over all faces of a given cell, check whether it is
  // part of the boundary that we want to denote by Gamma2, and if so set its
  // boundary indicator to <code>1</code>. For the present program, we
  // consider the left and bottom boundaries as Gamma2. We determine whether a
  // face is part of that boundary by asking whether the x or y coordinates
  // (i.e. vector components 0 and 1) of the midpoint of a face equals -1, up
  // to some small wiggle room that we have to give since it is instable to
  // compare floating point numbers that are subject to round off in
  // intermediate computations.
  //
  // It is worth noting that we have to loop over all cells here, not only the
  // active ones. The reason is that upon refinement, newly created faces
  // inherit the boundary indicator of their parent face. If we now only set
  // the boundary indicator for active faces, coarsen some cells and refine
  // them later on, they will again have the boundary indicator of the parent
  // cell which we have not modified, instead of the one we
  // intended. Consequently, we have to change the boundary indicators of
  // faces of all cells on Gamma2, whether they are active or not.
  // Alternatively, we could of course have done this job on the coarsest mesh
  // (i.e. before the first refinement step) and refined the mesh only after
  // that.
  template <int dim>
  void HelmholtzProblem<dim>::run()
  {
    const unsigned int n_cycles =
      (refinement_mode == global_refinement) ? 5 : 9;
    for (unsigned int cycle = 0; cycle < n_cycles; ++cycle)
      {
        if (cycle == 0)
          {
            GridGenerator::hyper_cube(triangulation, -1., 1.);
            triangulation.refine_global(3);

            for (const auto &cell : triangulation.cell_iterators())
              for (const auto &face : cell->face_iterators())
                {
                  const auto center = face->center();
                  if ((std::fabs(center[0] - (-1.0)) < 1e-12) ||
                      (std::fabs(center[1] - (-1.0)) < 1e-12))
                    face->set_boundary_id(1);
                }
          }
        else
          refine_grid();


        // The next steps are already known from previous examples. This is
        // mostly the basic set-up of every finite element program:
        setup_system();

        assemble_system();
        solve();

        // The last step in this chain of function calls is usually the
        // evaluation of the computed solution for the quantities one is
        // interested in. This is done in the following function. Since the
        // function generates output that indicates the number of the present
        // refinement step, we pass this number as an argument.
        process_solution(cycle);
      }

    // @sect5{Output of graphical data}

    // After the last iteration we output the solution on the finest
    // grid. This is done using the following sequence of statements which we
    // have already discussed in previous examples. The first step is to
    // generate a suitable filename (called <code>vtk_filename</code> here,
    // since we want to output data in VTK format; we add the prefix to
    // distinguish the filename from that used for other output files further
    // down below). Here, we augment the name by the mesh refinement
    // algorithm, and as above we make sure that we abort the program if
    // another refinement method is added and not handled by the following
    // switch statement:
    std::string vtk_filename;
    switch (refinement_mode)
      {
        case global_refinement:
          vtk_filename = "solution-global";
          break;
        case adaptive_refinement:
          vtk_filename = "solution-adaptive";
          break;
        default:
          DEAL_II_ASSERT_UNREACHABLE();
      }

    // We augment the filename by a postfix denoting the finite element which
    // we have used in the computation. To this end, the finite element base
    // class stores the maximal polynomial degree of shape functions in each
    // coordinate variable as a variable <code>degree</code>, which we append
    // as "-q1", "-q2", etc., to the filename.
    vtk_filename += "-q" + std::to_string(fe->degree);

    // Once we have the base name for the output file, we add an extension
    // appropriate for VTK output, open a file, and add the solution vector to
    // the object that will do the actual output:
    vtk_filename += ".vtk";
    std::ofstream output(vtk_filename);

    DataOut<dim> data_out;
    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution, "solution");

    // Now building the intermediate format as before is the next step. We
    // introduce one more feature of deal.II here. The background is the
    // following: in some of the runs of this function, we have used
    // biquadratic finite elements. However, since almost all output formats
    // only support bilinear data, the data is written only bilinear, and
    // information is consequently lost.  Of course, we can't change the
    // format in which graphic programs accept their inputs, but we can write
    // the data differently such that we more closely resemble the information
    // available in the quadratic approximation. We can, for example, write
    // each cell as four sub-cells with bilinear data each, such that we have
    // nine data points for each cell in the triangulation. The graphic
    // programs will, of course, display this data still only bilinear, but at
    // least we have given some more of the information we have.
    //
    // In order to allow writing more than one sub-cell per actual cell, the
    // <code>build_patches</code> function accepts a parameter (the default is
    // <code>1</code>, which is why you haven't seen this parameter in
    // previous examples). This parameter denotes into how many sub-cells per
    // space direction each cell shall be subdivided for output. For example,
    // if you give <code>2</code>, this leads to 4 cells in 2d and 8 cells in
    // 3d. For quadratic elements, two sub-cells per space direction is
    // obviously the right choice, so this is what we choose. In general, for
    // elements of polynomial order <code>q</code>, we use <code>q</code>
    // subdivisions, and the order of the elements is determined in the same
    // way as above.
    //
    // With the intermediate format so generated, we can then actually write
    // the graphical output:
    data_out.build_patches(fe->degree);
    data_out.write_vtk(output);

    // @sect5{Output of convergence tables}

    // After graphical output, we would also like to generate tables from the
    // error computations we have done in
    // <code>process_solution</code>. There, we have filled a table object
    // with the number of cells for each refinement step as well as the errors
    // in different norms.

    // For a nicer textual output of this data, one may want to set the
    // precision with which the values will be written upon output. We use 3
    // digits for this, which is usually sufficient for error norms. By
    // default, data is written in fixed point notation. However, for columns
    // one would like to see in scientific notation another function call sets
    // the <code>scientific_flag</code> to <code>true</code>, leading to
    // floating point representation of numbers.
    convergence_table.set_precision("L2", 3);
    convergence_table.set_precision("H1", 3);
    convergence_table.set_precision("Linfty", 3);

    convergence_table.set_scientific("L2", true);
    convergence_table.set_scientific("H1", true);
    convergence_table.set_scientific("Linfty", true);

    // For the output of a table into a LaTeX file, the default captions of
    // the columns are the keys given as argument to the
    // <code>add_value</code> functions. To have TeX captions that differ from
    // the default ones you can specify them by the following function calls.
    // Note, that `\\' is reduced to `\' by the compiler such that the real
    // TeX caption is, e.g., `$L^\infty$-error'.
    convergence_table.set_tex_caption("cells", "\\# cells");
    convergence_table.set_tex_caption("dofs", "\\# dofs");
    convergence_table.set_tex_caption("L2", "$L^2$-error");
    convergence_table.set_tex_caption("H1", "$H^1$-error");
    convergence_table.set_tex_caption("Linfty", "$L^\\infty$-error");

    // Finally, the default LaTeX format for each column of the table is `c'
    // (centered). To specify a different (e.g. `right') one, the following
    // function may be used:
    convergence_table.set_tex_format("cells", "r");
    convergence_table.set_tex_format("dofs", "r");

    // After this, we can finally write the table to the standard output
    // stream <code>std::cout</code> (after one extra empty line, to make
    // things look prettier). Note, that the output in text format is quite
    // simple and that captions may not be printed directly above the specific
    // columns.
    std::cout << std::endl;
    convergence_table.write_text(std::cout);

    // The table can also be written into a LaTeX file.  The (nicely)
    // formatted table can be viewed after calling `latex filename.tex` and
    // whatever output viewer you prefer, where filename is the name of the file
    // to which we will write output. We construct the file name in the same way
    // as before, but with a different prefix "error":
    std::string error_filename = "error";
    switch (refinement_mode)
      {
        case global_refinement:
          error_filename += "-global";
          break;
        case adaptive_refinement:
          error_filename += "-adaptive";
          break;
        default:
          DEAL_II_ASSERT_UNREACHABLE();
      }

    error_filename += "-q" + std::to_string(fe->degree);
    error_filename += ".tex";
    std::ofstream error_table_file(error_filename);

    convergence_table.write_tex(error_table_file);


    // @sect5{Further table manipulations}

    // In case of global refinement, it might be of interest to also output
    // the convergence rates. This may be done by the functionality the
    // ConvergenceTable offers over the regular TableHandler. However, we do
    // it only for global refinement, since for adaptive refinement the
    // determination of something like an order of convergence is somewhat
    // more involved. While we are at it, we also show a few other things that
    // can be done with tables.
    if (refinement_mode == global_refinement)
      {
        // The first thing is that one can group individual columns together
        // to form so-called super columns. Essentially, the columns remain
        // the same, but the ones that were grouped together will get a
        // caption running across all columns in a group. For example, let's
        // merge the "cycle" and "cells" columns into a super column named "n
        // cells":
        convergence_table.add_column_to_supercolumn("cycle", "n cells");
        convergence_table.add_column_to_supercolumn("cells", "n cells");

        // Next, it isn't necessary to always output all columns, or in the
        // order in which they were originally added during the run.
        // Selecting and re-ordering the columns works as follows (note that
        // this includes super columns):
        std::vector<std::string> new_order;
        new_order.emplace_back("n cells");
        new_order.emplace_back("H1");
        new_order.emplace_back("L2");
        convergence_table.set_column_order(new_order);

        // For everything that happened to the ConvergenceTable until this
        // point, it would have been sufficient to use a simple
        // TableHandler. Indeed, the ConvergenceTable is derived from the
        // TableHandler but it offers the additional functionality of
        // automatically evaluating convergence rates. For example, here is
        // how we can let the table compute reduction and convergence rates
        // (convergence rates are the binary logarithm of the reduction rate):
        convergence_table.evaluate_convergence_rates(
          "L2", ConvergenceTable::reduction_rate);
        convergence_table.evaluate_convergence_rates(
          "L2", ConvergenceTable::reduction_rate_log2);
        convergence_table.evaluate_convergence_rates(
          "H1", ConvergenceTable::reduction_rate);
        convergence_table.evaluate_convergence_rates(
          "H1", ConvergenceTable::reduction_rate_log2);
        // Each of these function calls produces an additional column that is
        // merged with the original column (in our example the `L2' and the
        // `H1' column) to a supercolumn.

        // Finally, we want to write this convergence chart again, first to
        // the screen and then, in LaTeX format, to disk. The filename is
        // again constructed as above.
        std::cout << std::endl;
        convergence_table.write_text(std::cout);

        std::string conv_filename = "convergence";
        switch (refinement_mode)
          {
            case global_refinement:
              conv_filename += "-global";
              break;
            case adaptive_refinement:
              conv_filename += "-adaptive";
              break;
            default:
              DEAL_II_ASSERT_UNREACHABLE();
          }
        conv_filename += "-q" + std::to_string(fe->degree);
        conv_filename += ".tex";

        std::ofstream table_file(conv_filename);
        convergence_table.write_tex(table_file);
      }
  }

  // The final step before going to <code>main()</code> is then to close the
  // namespace <code>Step7</code> into which we have put everything we needed
  // for this program:
} // namespace Step7

// @sect3{Main function}

// The main function is mostly as before. The only difference is that we solve
// three times, once for Q1 and adaptive refinement, once for Q1 elements and
// global refinement, and once for Q2 elements and global refinement.
//
// Since we instantiate several template classes below for two space
// dimensions, we make this more generic by declaring a constant at the
// beginning of the function denoting the number of space dimensions. If you
// want to run the program in 1d or 2d, you will then only have to change this
// one instance, rather than all uses below:
int main()
{
  const unsigned int dim = 2;

  try
    {
      using namespace dealii;
      using namespace Step7;

      // Now for the three calls to the main class. Each call is blocked into
      // curly braces in order to destroy the respective objects (i.e. the
      // finite element and the HelmholtzProblem object) at the end of the
      // block and before we go to the next run. This avoids conflicts with
      // variable names, and also makes sure that memory is released
      // immediately after one of the three runs has finished, and not only at
      // the end of the <code>try</code> block.
      {
        std::cout << "Solving with Q1 elements, adaptive refinement"
                  << std::endl
                  << "============================================="
                  << std::endl
                  << std::endl;

        const FE_Q<dim>       fe(1);
        HelmholtzProblem<dim> helmholtz_problem_2d(
          fe, HelmholtzProblem<dim>::adaptive_refinement);

        helmholtz_problem_2d.run();

        std::cout << std::endl;
      }

      {
        std::cout << "Solving with Q1 elements, global refinement" << std::endl
                  << "===========================================" << std::endl
                  << std::endl;

        const FE_Q<dim>       fe(1);
        HelmholtzProblem<dim> helmholtz_problem_2d(
          fe, HelmholtzProblem<dim>::global_refinement);

        helmholtz_problem_2d.run();

        std::cout << std::endl;
      }

      {
        std::cout << "Solving with Q2 elements, global refinement" << std::endl
                  << "===========================================" << std::endl
                  << std::endl;

        const FE_Q<dim>       fe(2);
        HelmholtzProblem<dim> helmholtz_problem_2d(
          fe, HelmholtzProblem<dim>::global_refinement);

        helmholtz_problem_2d.run();

        std::cout << std::endl;
      }
      {
        std::cout << "Solving with Q2 elements, adaptive refinement"
                  << std::endl
                  << "===========================================" << std::endl
                  << std::endl;

        const FE_Q<dim>       fe(2);
        HelmholtzProblem<dim> helmholtz_problem_2d(
          fe, HelmholtzProblem<dim>::adaptive_refinement);

        helmholtz_problem_2d.run();

        std::cout << std::endl;
      }
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
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
      std::cerr << std::endl
                << std::endl
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
