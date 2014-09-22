// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2014 by the deal.II authors
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

#ifndef __deal2__matrix_tools_h
#define __deal2__matrix_tools_h


#include <deal.II/base/config.h>
#include <deal.II/base/function.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/dofs/function_map.h>

#include <map>

DEAL_II_NAMESPACE_OPEN


// forward declarations
template <int dim> class Quadrature;


template<typename number> class Vector;
template<typename number> class FullMatrix;
template<typename number> class SparseMatrix;

template <typename number> class BlockSparseMatrix;
template <typename Number> class BlockVector;

template <int dim, int spacedim> class Mapping;
template <int dim, int spacedim> class DoFHandler;
template <int dim, int spacedim> class MGDoFHandler;
template <int dim, int spacedim> class FEValues;

namespace hp
{
  template <int> class QCollection;
  template <int, int> class MappingCollection;
  template <int, int> class DoFHandler;
}


#ifdef DEAL_II_WITH_PETSC
namespace PETScWrappers
{
  class SparseMatrix;
  class Vector;
  namespace MPI
  {
    class SparseMatrix;
    class BlockSparseMatrix;
    class Vector;
    class BlockVector;
  }
}
#endif

#ifdef DEAL_II_WITH_TRILINOS
namespace TrilinosWrappers
{
  class SparseMatrix;
  class Vector;
  class BlockSparseMatrix;
  class BlockVector;
  namespace MPI
  {
    class Vector;
    class BlockVector;
  }
}
#endif


/**
 * This namespace provides functions that assemble certain standard matrices for a
 * given triangulation, using a given finite element, a given mapping and a
 * quadrature formula.
 *
 *
 * <h3>Conventions for all functions</h3>
 *
 * There exist two versions of each function. One with a Mapping
 * argument and one without. If a code uses a mapping different from
 * MappingQ1 the functions <em>with</em> mapping argument should
 * be used. Code that uses only MappingQ1 may also use the
 * functions <em>without</em> Mapping argument. Each of these
 * latter functions create a MappingQ1 object and just call the
 * respective functions with that object as mapping argument.
 *
 * All functions take a sparse matrix object to hold the matrix to be
 * created. The functions assume that the matrix is initialized with a
 * sparsity pattern (SparsityPattern) corresponding to the given degree
 * of freedom handler, i.e. the sparsity structure is already as needed.
 * You can do this by calling the DoFTools::make_sparsity_pattern()
 * function.
 *
 * Furthermore it is assumed that no relevant data is in the matrix. Some
 * entries will be overwritten and some others will contain invalid data if
 * the matrix wasn't empty before. Therefore you may want to clear the matrix
 * before assemblage.
 *
 * By default, all created matrices are `raw': they are not condensed,
 * i.e. hanging nodes are not eliminated. The reason is that you may
 * want to add several matrices and could then condense afterwards
 * only once, instead of for every matrix. To actually do computations
 * with these matrices, you have to condense the matrix using the
 * ConstraintMatrix::condense function; you also have to
 * condense the right hand side accordingly and distribute the
 * solution afterwards. Alternatively, you can give an optional argument
 * ConstraintMatrix that writes cell matrix (and vector) entries with
 * distribute_local_to_global into the global matrix and vector. This way,
 * adding several matrices from different sources is more complicated and
 * you should make sure that you do not mix different ways of applying
 * constraints. Particular caution is necessary when the given
 * constraint matrix contains inhomogeneous constraints: In that case, the
 * matrix assembled this way must be the only matrix (or you need to
 * assemble the <b>same</b> right hand side for <b>every</b> matrix
 * you generate and add together).
 *
 * If you want to use boundary conditions with the matrices generated
 * by the functions of this class in addition to the ones in a possible
 * constraint matrix, you have to use a function like
 * <tt>ProblemBase<>::apply_dirichlet_bc</tt> to matrix and right hand
 * side.
 *
 *
 * <h3>Supported matrices</h3>
 *
 * At present there are functions to create the following matrices:
 * <ul>
 * <li> @p create_mass_matrix: create the matrix with entries
 *   $m_{ij} = \int_\Omega \phi_i(x) \phi_j(x) dx$ by numerical
 *   quadrature. Here, the $\phi_i$ are the basis functions of the
 *   finite element space given.
 *
 *   A coefficient may be given to evaluate
 *   $m_{ij} = \int_\Omega a(x) \phi_i(x) \phi_j(x) dx$ instead.
 *
 * <li> @p create_laplace_matrix: create the matrix with entries
 *   $a_{ij} = \int_\Omega \nabla\phi_i(x) \nabla\phi_j(x) dx$ by
 *   numerical quadrature.
 *
 *   Again, a coefficient may be given to evaluate
 *   $a_{ij} = \int_\Omega a(x) \nabla\phi_i(x) \nabla\phi_j(x) dx$ instead.
 * </ul>
 *
 * Make sure that the order of the Quadrature formula given to these
 * functions is sufficiently high to compute the matrices with the
 * required accuracy. For the choice of this quadrature rule you need
 * to take into account the polynomial degree of the FiniteElement
 * basis functions, the roughness of the coefficient @p a, as well as
 * the degree of the given @p Mapping (if any).
 *
 * Note, that for vector-valued elements the mass matrix and the
 * laplace matrix is implemented in such a way that each component
 * couples only with itself, i.e. there is no coupling of shape
 * functions belonging to different components. If the degrees of
 * freedom have been sorted according to their vector component (e.g.,
 * using DoFRenumbering::component_wise()), then the resulting
 * matrices will be block diagonal.
 *
 * If the finite element for which the mass matrix or the Laplace
 * matrix is to be built has more than one component, the functions
 * accept a single coefficient as well as a vector valued coefficient
 * function. For the latter case, the number of components must
 * coincide with the number of components of the system finite
 * element.
 *
 *
 * <h3>Matrices on the boundary</h3>
 *
 * The create_boundary_mass_matrix() creates the matrix with entries $m_{ij} =
 * \int_{\Gamma} \phi_i \phi_j dx$, where $\Gamma$ is the union of boundary
 * parts with indicators contained in a FunctionMap passed to the function
 * (i.e. if you want to set up the mass matrix for the parts of the boundary
 * with indicators zero and 2, you pass the function a map of <tt>unsigned
 * char</tt>s as parameter @p boundary_functions containing the keys zero and
 * 2). The size of the matrix is equal to the number of degrees of freedom
 * that have support on the boundary, i.e. it is <em>not</em> a matrix on all
 * degrees of freedom, but only a subset. (The $\phi_i$ in the formula are
 * this subsect of basis functions which have at least part of their support
 * on $\Gamma$.) In order to determine which shape functions are to be
 * considered, and in order to determine in which order, the function takes a
 * @p dof_to_boundary_mapping; this object maps global DoF numbers to a
 * numbering of the degrees of freedom located on the boundary, and can be
 * obtained using the function DoFTools::map_dof_to_boundary_indices().
 *
 * In order to work, the function needs a matrix of the correct size, built on
 * top of a corresponding sparsity pattern. Since we only work on a subset of
 * the degrees of freedom, we can't use the matrices and sparsity patterns
 * that are created for the entire set of degrees of freedom. Rather, you
 * should use the DoFHandler::make_boundary_sparsity_pattern() function to
 * create the correct sparsity pattern, and build a matrix on top of it.
 *
 * Note that at present there is no function that computes the mass matrix for
 * <em>all</em> shape functions, though such a function would be trivial to
 * implement.
 *
 *
 * <h3>Right hand sides</h3>
 *
 * In many cases, you will not only want to build the matrix, but also
 * a right hand side, which will give a vector with
 * $f_i = \int_\Omega f(x) \phi_i(x) dx$. For this purpose, each function
 * exists in two versions, one only building the matrix and one also
 * building the right hand side vector. If you want to create a right
 * hand side vector without creating a matrix, you can use the
 * VectorTools::create_right_hand_side() function. The use of the
 * latter may be useful if you want to create many right hand side
 * vectors.
 *
 * @ingroup numerics
 * @author Wolfgang Bangerth, 1998, Ralf Hartmann, 2001
 */
namespace MatrixCreator
{
  /**
   * Assemble the mass matrix. If no coefficient is given (i.e., if
   * the pointer to a function object is zero as it is by default),
   * the coefficient is taken as being constant and equal to one.
   *
   * If the library is configured to use multithreading, this function works
   * in parallel.
   *
   * The optional argument @p constraints allows to apply constraints
   * on the resulting matrix directly. Note, however, that this
   * becomes difficult when you have inhomogeneous constraints and
   * later want to add several such matrices, for example in time
   * dependent settings such as the main loop of step-26.
   *
   * See the general doc of this class for more information.
   */
  template <int dim, typename number, int spacedim>
  void create_mass_matrix (const Mapping<dim, spacedim>       &mapping,
                           const DoFHandler<dim,spacedim>    &dof,
                           const Quadrature<dim>    &q,
                           SparseMatrix<number>     &matrix,
                           const Function<spacedim> *const a = 0,
                           const ConstraintMatrix   &constraints = ConstraintMatrix());

  /**
   * Calls the create_mass_matrix() function, see above, with
   * <tt>mapping=MappingQ1@<dim@>()</tt>.
   */
  template <int dim, typename number, int spacedim>
  void create_mass_matrix (const DoFHandler<dim,spacedim>    &dof,
                           const Quadrature<dim>    &q,
                           SparseMatrix<number>     &matrix,
                           const Function<spacedim> *const a = 0,
                           const ConstraintMatrix   &constraints = ConstraintMatrix());

  /**
   * Assemble the mass matrix and a right hand side vector. If no
   * coefficient is given (i.e., if the pointer to a function object
   * is zero as it is by default), the coefficient is taken as being
   * constant and equal to one.
   *
   * If the library is configured to use multithreading, this function works
   * in parallel.
   *
   * The optional argument @p constraints allows to apply constraints on the
   * resulting matrix directly. Note, however, that this
   * becomes difficult when you have inhomogeneous constraints and
   * later want to add several such matrices, for example in time
   * dependent settings such as the main loop of step-26.
   *
   * See the general doc of this class for more information.
   */
  template <int dim, typename number, int spacedim>
  void create_mass_matrix (const Mapping<dim, spacedim>   &mapping,
                           const DoFHandler<dim,spacedim> &dof,
                           const Quadrature<dim>    &q,
                           SparseMatrix<number>     &matrix,
                           const Function<spacedim> &rhs,
                           Vector<double>           &rhs_vector,
                           const Function<spacedim> *const a = 0,
                           const ConstraintMatrix   &constraints = ConstraintMatrix());

  /**
   * Calls the create_mass_matrix() function, see above, with
   * <tt>mapping=MappingQ1@<dim@>()</tt>.
   */
  template <int dim, typename number, int spacedim>
  void create_mass_matrix (const DoFHandler<dim,spacedim> &dof,
                           const Quadrature<dim>    &q,
                           SparseMatrix<number>     &matrix,
                           const Function<spacedim> &rhs,
                           Vector<double>           &rhs_vector,
                           const Function<spacedim> *const a = 0,
                           const ConstraintMatrix   &constraints = ConstraintMatrix());

  /**
   * Same function as above, but for hp objects.
   */
  template <int dim, typename number, int spacedim>
  void create_mass_matrix (const hp::MappingCollection<dim,spacedim>       &mapping,
                           const hp::DoFHandler<dim,spacedim>    &dof,
                           const hp::QCollection<dim>    &q,
                           SparseMatrix<number>     &matrix,
                           const Function<spacedim> *const a = 0,
                           const ConstraintMatrix   &constraints = ConstraintMatrix());

  /**
   * Same function as above, but for hp objects.
   */
  template <int dim, typename number, int spacedim>
  void create_mass_matrix (const hp::DoFHandler<dim,spacedim>    &dof,
                           const hp::QCollection<dim>    &q,
                           SparseMatrix<number>     &matrix,
                           const Function<spacedim> *const a = 0,
                           const ConstraintMatrix   &constraints = ConstraintMatrix());

  /**
   * Same function as above, but for hp objects.
   */
  template <int dim, typename number, int spacedim>
  void create_mass_matrix (const hp::MappingCollection<dim,spacedim> &mapping,
                           const hp::DoFHandler<dim,spacedim> &dof,
                           const hp::QCollection<dim> &q,
                           SparseMatrix<number>     &matrix,
                           const Function<spacedim> &rhs,
                           Vector<double>           &rhs_vector,
                           const Function<spacedim> *const a = 0,
                           const ConstraintMatrix   &constraints = ConstraintMatrix());

  /**
   * Same function as above, but for hp objects.
   */
  template <int dim, typename number, int spacedim>
  void create_mass_matrix (const hp::DoFHandler<dim,spacedim> &dof,
                           const hp::QCollection<dim> &q,
                           SparseMatrix<number>     &matrix,
                           const Function<spacedim> &rhs,
                           Vector<double>           &rhs_vector,
                           const Function<spacedim> *const a = 0,
                           const ConstraintMatrix   &constraints = ConstraintMatrix());


  /**
   * Assemble the mass matrix and a right hand side vector along the boundary.
   *
   * The matrix is assumed to already be initialized with a suiting sparsity
   * pattern (the DoFHandler provides an appropriate function).
   *
   * If the library is configured to use multithreading, this function works
   * in parallel.
   *
   * @arg @p weight: an optional weight for the computation of the mass
   * matrix. If no weight is given, it is set to one.
   *
   * @arg @p component_mapping: if the components in @p boundary_functions and
   * @p dof do not coincide, this vector allows them to be remapped. If the
   * vector is not empty, it has to have one entry for each component in @p
   * dof. This entry is the component number in @p boundary_functions that
   * should be used for this component in @p dof. By default, no remapping is
   * applied.
   *
   * @todo This function does not work for finite elements with cell-dependent
   * shape functions.
   */
  template <int dim, int spacedim>

  void create_boundary_mass_matrix (const Mapping<dim, spacedim>       &mapping,
                                    const DoFHandler<dim,spacedim>    &dof,
                                    const Quadrature<dim-1>  &q,
                                    SparseMatrix<double>     &matrix,
                                    const typename FunctionMap<spacedim>::type &boundary_functions,
                                    Vector<double>           &rhs_vector,
                                    std::vector<types::global_dof_index> &dof_to_boundary_mapping,
                                    const Function<spacedim> *const weight = 0,
                                    std::vector<unsigned int> component_mapping = std::vector<unsigned int>());

//                                    * Same function, but for 1d.
//                                    */
//
//     void create_boundary_mass_matrix (const Mapping<1,1>       &mapping,
//                                    const DoFHandler<1,1>    &dof,
//                                    const Quadrature<0>    &q,
//                                    SparseMatrix<double>   &matrix,
//                                    const FunctionMap<1>::type &boundary_functions,
//                                    Vector<double>         &rhs_vector,
//                                    std::vector<types::global_dof_index>&dof_to_boundary_mapping,
//                                    const Function<1> * const a = 0);
// //codimension 1
//
//     void create_boundary_mass_matrix (const Mapping<1,2>       &mapping,
//                                    const DoFHandler<1,2>    &dof,
//                                    const Quadrature<0>    &q,
//                                    SparseMatrix<double>   &matrix,
//                                    const FunctionMap<2>::type &boundary_functions,
//                                    Vector<double>         &rhs_vector,
//                                    std::vector<types::global_dof_index>&dof_to_boundary_mapping,
//                                    const Function<2> * const a = 0);



  /**
   * Calls the create_boundary_mass_matrix() function, see above, with
   * <tt>mapping=MappingQ1@<dim@>()</tt>.
   */
  template <int dim, int spacedim>

  void create_boundary_mass_matrix (const DoFHandler<dim,spacedim>    &dof,
                                    const Quadrature<dim-1>  &q,
                                    SparseMatrix<double>     &matrix,
                                    const typename FunctionMap<spacedim>::type        &boundary_functions,
                                    Vector<double>           &rhs_vector,
                                    std::vector<types::global_dof_index> &dof_to_boundary_mapping,
                                    const Function<spacedim> *const a = 0,
                                    std::vector<unsigned int> component_mapping = std::vector<unsigned int>());

  /**
   * Same function as above, but for hp objects.
   */
  template <int dim, int spacedim>

  void create_boundary_mass_matrix (const hp::MappingCollection<dim,spacedim>       &mapping,
                                    const hp::DoFHandler<dim,spacedim>    &dof,
                                    const hp::QCollection<dim-1>  &q,
                                    SparseMatrix<double>     &matrix,
                                    const typename FunctionMap<spacedim>::type &boundary_functions,
                                    Vector<double>           &rhs_vector,
                                    std::vector<types::global_dof_index> &dof_to_boundary_mapping,
                                    const Function<spacedim> *const a = 0,
                                    std::vector<unsigned int> component_mapping = std::vector<unsigned int>());

  /**
   * Same function as above, but for hp objects.
   */
//
//     void create_boundary_mass_matrix (const hp::MappingCollection<1,1>       &mapping,
//                                    const hp::DoFHandler<1,1>    &dof,
//                                    const hp::QCollection<0>    &q,
//                                    SparseMatrix<double>   &matrix,
//                                    const FunctionMap<1>::type &boundary_functions,
//                                    Vector<double>         &rhs_vector,
//                                    std::vector<types::global_dof_index>&dof_to_boundary_mapping,
//                                    const Function<1> * const a = 0);

  /**
   * Same function as above, but for hp objects.
   */
  template <int dim, int spacedim>

  void create_boundary_mass_matrix (const hp::DoFHandler<dim,spacedim>    &dof,
                                    const hp::QCollection<dim-1>  &q,
                                    SparseMatrix<double>     &matrix,
                                    const typename FunctionMap<spacedim>::type        &boundary_functions,
                                    Vector<double>           &rhs_vector,
                                    std::vector<types::global_dof_index> &dof_to_boundary_mapping,
                                    const Function<spacedim> *const a = 0,
                                    std::vector<unsigned int> component_mapping = std::vector<unsigned int>());

  /**
   * Assemble the Laplace matrix. If no coefficient is given (i.e., if
   * the pointer to a function object is zero as it is by default),
   * the coefficient is taken as being constant and equal to one.
   *
   * If the library is configured to use multithreading, this function works
   * in parallel.
   *
   * The optional argument @p constraints allows to apply constraints on the
   * resulting matrix directly. Note, however, that this
   * becomes difficult when you have inhomogeneous constraints and
   * later want to add several such matrices, for example in time
   * dependent settings such as the main loop of step-26.
   *
   * See the general doc of this class for more information.
   */
  template <int dim, int spacedim>
  void create_laplace_matrix (const Mapping<dim, spacedim>   &mapping,
                              const DoFHandler<dim,spacedim> &dof,
                              const Quadrature<dim>    &q,
                              SparseMatrix<double>     &matrix,
                              const Function<spacedim> *const a = 0,
                              const ConstraintMatrix   &constraints = ConstraintMatrix());

  /**
   * Calls the create_laplace_matrix() function, see above, with
   * <tt>mapping=MappingQ1@<dim@>()</tt>.
   */
  template <int dim, int spacedim>
  void create_laplace_matrix (const DoFHandler<dim,spacedim> &dof,
                              const Quadrature<dim>    &q,
                              SparseMatrix<double>     &matrix,
                              const Function<spacedim> *const a = 0,
                              const ConstraintMatrix   &constraints = ConstraintMatrix());

  /**
   * Assemble the Laplace matrix and a right hand side vector. If no
   * coefficient is given, it is assumed to be constant one.
   *
   * If the library is configured to use multithreading, this function works
   * in parallel.
   *
   * The optional argument @p constraints allows to apply constraints on the
   * resulting matrix directly. Note, however, that this
   * becomes difficult when you have inhomogeneous constraints and
   * later want to add several such matrices, for example in time
   * dependent settings such as the main loop of step-26.
   *
   * See the general doc of this class for more information.
   */
  template <int dim, int spacedim>
  void create_laplace_matrix (const Mapping<dim, spacedim>   &mapping,
                              const DoFHandler<dim,spacedim> &dof,
                              const Quadrature<dim>    &q,
                              SparseMatrix<double>     &matrix,
                              const Function<spacedim> &rhs,
                              Vector<double>           &rhs_vector,
                              const Function<spacedim> *const a = 0,
                              const ConstraintMatrix   &constraints = ConstraintMatrix());

  /**
   * Calls the create_laplace_matrix() function, see above, with
   * <tt>mapping=MappingQ1@<dim@>()</tt>.
   */
  template <int dim, int spacedim>
  void create_laplace_matrix (const DoFHandler<dim,spacedim> &dof,
                              const Quadrature<dim>    &q,
                              SparseMatrix<double>     &matrix,
                              const Function<spacedim> &rhs,
                              Vector<double>           &rhs_vector,
                              const Function<spacedim> *const a = 0,
                              const ConstraintMatrix   &constraints = ConstraintMatrix());

  /**
   * Like the functions above, but for hp dof handlers, mappings, and
   * quadrature collections.
   */
  template <int dim, int spacedim>
  void create_laplace_matrix (const hp::MappingCollection<dim,spacedim> &mapping,
                              const hp::DoFHandler<dim,spacedim> &dof,
                              const hp::QCollection<dim>    &q,
                              SparseMatrix<double>     &matrix,
                              const Function<spacedim> *const a = 0,
                              const ConstraintMatrix   &constraints = ConstraintMatrix());

  /**
   * Like the functions above, but for hp dof handlers, mappings, and
   * quadrature collections.
   */
  template <int dim, int spacedim>
  void create_laplace_matrix (const hp::DoFHandler<dim,spacedim> &dof,
                              const hp::QCollection<dim>    &q,
                              SparseMatrix<double>     &matrix,
                              const Function<spacedim> *const a = 0,
                              const ConstraintMatrix   &constraints = ConstraintMatrix());

  /**
   * Like the functions above, but for hp dof handlers, mappings, and
   * quadrature collections.
   */
  template <int dim, int spacedim>
  void create_laplace_matrix (const hp::MappingCollection<dim,spacedim> &mapping,
                              const hp::DoFHandler<dim,spacedim> &dof,
                              const hp::QCollection<dim>    &q,
                              SparseMatrix<double>     &matrix,
                              const Function<spacedim>      &rhs,
                              Vector<double>           &rhs_vector,
                              const Function<spacedim> *const a = 0,
                              const ConstraintMatrix   &constraints = ConstraintMatrix());

  /**
   * Like the functions above, but for hp dof handlers, mappings, and
   * quadrature collections.
   */
  template <int dim, int spacedim>
  void create_laplace_matrix (const hp::DoFHandler<dim,spacedim> &dof,
                              const hp::QCollection<dim>    &q,
                              SparseMatrix<double>     &matrix,
                              const Function<spacedim>      &rhs,
                              Vector<double>           &rhs_vector,
                              const Function<spacedim> *const a = 0,
                              const ConstraintMatrix   &constraints = ConstraintMatrix());

  /**
   * Exception
   */
  DeclException0 (ExcComponentMismatch);
}



/**
 * Provide a collection of functions operating on matrices. These include
 * the application of boundary conditions to a linear system of equations
 * and others.
 *
 *
 * <h3>Boundary conditions</h3>
 *
 * The apply_boundary_values() function inserts boundary conditions
 * into a system of equations.  To actually do this you have to
 * specify a list of degree of freedom indices along with the values
 * these degrees of freedom shall assume. To see how to get such a
 * list, see the discussion of the
 * VectorTools::interpolate_boundary_values function.
 *
 * There are two ways to incorporate fixed degrees of freedom such as boundary
 * nodes into a linear system, as discussed below.
 *
 *
 * <h3>Global elimination</h3>
 *
 * In the first method, we first assemble the global linear system without
 * respect for fixed degrees of freedom, and in a second step eliminate them
 * again from the linear system. The inclusion into the assembly process is as
 * follows: when the matrix and vectors are set up, a list of nodes subject to
 * Dirichlet bc is made and matrix and vectors are modified accordingly. This
 * is done by deleting all entries in the matrix in the line of this degree of
 * freedom, setting the main diagonal entry to a suitable positive value and
 * the right hand side element to a value so that the solution of the linear
 * system will have the boundary value at this node. To decouple the remaining
 * linear system of equations and to make the system symmetric again (at least
 * if it was before), one Gauss elimination step is performed with this line,
 * by adding this (now almost empty) line to all other lines which couple with
 * the given degree of freedom and thus eliminating all coupling between this
 * degree of freedom and others. Now the respective column also consists only
 * of zeroes, apart from the main diagonal entry. Alternatively, the functions
 * in this class take a boolean parameter that allows to omit this last step,
 * if symmetry of the resulting linear system is not required. Note that
 * usually even CG can cope with a non-symmetric linear system with this
 * particular structure.
 *
 * Finding which rows contain an entry in the column for which we are
 * presently performing a Gauss elimination step is either difficult
 * or very simple, depending on the circumstances. If the sparsity
 * pattern is symmetric (whether the matrix is symmetric is irrelevant
 * here), then we can infer the rows which have a nonzero entry in the
 * present column by looking at which columns in the present row are
 * nonempty. In this case, we only need to look into a fixed number of
 * rows and need not search all rows. On the other hand, if the
 * sparsity pattern is nonsymmetric, then we need to use an iterative
 * solver which can handle nonsymmetric matrices in any case, so there
 * may be no need to do the Gauss elimination anyway. In fact, this is
 * the way the function works: it takes a parameter
 * (@p eliminate_columns) that specifies whether the sparsity pattern
 * is symmetric; if so, then the column is eliminated and the right
 * hand side is also modified accordingly. If not, then only the row
 * is deleted and the column is not touched at all, and all right hand
 * side values apart from the one corresponding to the present row
 * remain unchanged.
 *
 * If the sparsity pattern for your matrix is non-symmetric, you must
 * set the value of this parameter to @p false in any case, since then
 * we can't eliminate the column without searching all rows, which
 * would be too expensive (if @p N be the number of rows, and @p m the
 * number of nonzero elements per row, then eliminating one column is
 * an <tt>O(N*log(m))</tt> operation, since searching in each row takes
 * <tt>log(m)</tt> operations). If your sparsity pattern is symmetric, but
 * your matrix is not, then you might specify @p false as well. If your
 * sparsity pattern and matrix are both symmetric, you might want to
 * specify @p true (the complexity of eliminating one row is then
 * <tt>O(m*log(m))</tt>, since we only have to search @p m rows for the
 * respective element of the column). Given the fact that @p m is
 * roughly constant, irrespective of the discretization, and that the
 * number of boundary nodes is <tt>sqrt(N)</tt> in 2d, the algorithm for
 * symmetric sparsity patterns is <tt>O(sqrt(N)*m*log(m))</tt>, while it
 * would be <tt>O(N*sqrt(N)*log(m))</tt> for the general case; the latter
 * is too expensive to be performed.
 *
 * It seems as if we had to make clear not to overwrite the lines of
 * other boundary nodes when doing the Gauss elimination
 * step. However, since we reset the right hand side when passing such
 * a node, it is not a problem to change the right hand side values of
 * other boundary nodes not yet processed. It would be a problem to
 * change those entries of nodes already processed, but since the
 * matrix entry of the present column on the row of an already
 * processed node is zero, the Gauss step does not change the right
 * hand side. We need therefore not take special care of other
 * boundary nodes.
 *
 * To make solving faster, we preset the solution vector with the right
 * boundary values (as to why this is necessary, see the discussion below in
 * the description of local elimination). It it not clear whether the deletion
 * of coupling between the boundary degree of freedom and other dofs really
 * forces the corresponding entry in the solution vector to have the right
 * value when using iterative solvers, since their search directions may
 * contain components in the direction of the boundary node. For this reason,
 * we perform a very simple line balancing by not setting the main diagonal
 * entry to unity, but rather to the value it had before deleting this line,
 * or to the first nonzero main diagonal entry if it is zero for some reason.
 * Of course we have to change the right hand side appropriately. This is not
 * a very good strategy, but it at least should give the main diagonal entry a
 * value in the right order of dimension, which makes the solution process a
 * bit more stable. A refined algorithm would set the entry to the mean of the
 * other diagonal entries, but this seems to be too expensive.
 *
 * In some cases, it might be interesting to solve several times with
 * the same matrix, but for different right hand sides or boundary
 * values. However, since the modification for boundary values of the
 * right hand side vector depends on the original matrix, this is not
 * possible without storing the original matrix somewhere and applying
 * the @p apply_boundary_conditions function to a copy of it each
 * time we want to solve. In that case, you can use the
 * FilteredMatrix class in the @p LAC sublibrary. There you can
 * also find a formal (mathematical) description of the process of
 * modifying the matrix and right hand side vectors for boundary
 * values.
 *
 *
 * <h3>Local elimination</h3>
 *
 * The second way of handling boundary values is to modify the local
 * matrix and vector contributions appropriately before transferring
 * them into the global sparse matrix and vector. This is what
 * local_apply_boundary_values() does. The advantage is that we save
 * the call to the apply_boundary_values function (which is expensive
 * because it has to work on sparse data structures). On the other
 * hand, the local_apply_boundary_values() function is called many
 * times, even if we only have a very small number of fixed boundary
 * nodes, and the main drawback is that this function doesn't work as
 * expected if there are hanging nodes that also need to be
 * treated. The reason that this function doesn't work is that it is
 * meant to be run before distribution into the global matrix,
 * i.e. before hanging nodes are distributed; since hanging nodes can
 * be constrained to a boundary node, the treatment of hanging nodes
 * can add entries again to rows and columns corresponding to boundary
 * values and that we have already vacated in the local elimination
 * step. To make things worse, in 3d constrained nodes can even lie on
 * the boundary. Thus, it is imperative that boundary node elimination
 * happens @em after hanging node elimination, but this can't be
 * achieved with local elimination of boundary nodes unless there are
 * no hanging node constraints at all.
 *
 * Local elimination has one additional drawback: we don't have access
 * to the solution vector, only to the local contributions to the
 * matrix and right hand side. The problem with this is subtle, but
 * can lead to very hard to find difficulties: when we eliminate a
 * degree of freedom, we delete the row and column of this unknown,
 * and set the diagonal entry to some positive value. To make the
 * problem more or less well-conditioned, we set this diagonal entry
 * to the absolute value of its prior value if that was non-zero, or
 * to the average magnitude of all other nonzero diagonal
 * elements. Then we set the right hand side value such that the
 * resulting solution entry has the right value as given by the
 * boundary values. Since we add these contributions up over all local
 * contributions, the diagonal entry and the respective value in the
 * right hand side are added up correspondingly, so that the entry in
 * the solution of the linear system is still valid.
 *
 * A problem arises, however, if the diagonal entries so chosen are not
 * appropriate for the linear system. Consider, for example, a mixed Laplace
 * problem with matrix <tt>[[A B][C^T 0]]</tt>, where we only specify boundary
 * values for the second component of the solution. In the mixed formulation,
 * the stress-strain tensor only appears in either the matrix @p B or @p C, so
 * one of them may be significantly larger or smaller than the other one. Now,
 * if we eliminate boundary values, we delete some rows and columns, but we
 * also introduce a few entries on the diagonal of the lower right block, so
 * that we get the system <tt>[[A' B'][C'^T X]]</tt>. The diagonal entries in
 * the matrix @p X will be of the same order of magnitude as those in @p
 * A. Now, if we solve this system in the Schur complement formulation, we
 * have to invert the matrix <tt>X-C'^TA'^{-1}B'</tt>. Deleting rows and
 * columns above makes sure that boundary nodes indeed have empty rows and
 * columns in the Schur complement as well, except for the entries in @p
 * X. However, the entries in @p X may be of significantly different orders of
 * magnitude than those in <tt>C'^TA'^{-1}B'</tt>! If this is the case, we may
 * run into trouble with iterative solvers. For example, assume that we start
 * with zero entries in the solution vector and that the entries in @p X are
 * several orders of magnitude too small; in this case, iterative solvers will
 * compute the residual vector in each step and form correction vectors, but
 * since the entries in @p X are so small, the residual contributions for
 * boundary nodes are really small, despite the fact that the boundary nodes
 * are still at values close to zero and not in accordance with the prescribed
 * boundary values. Since the residual is so small, the corrections the
 * iterative solver computes are very small, and in the end the solver will
 * indicate convergence to a small total residual with the boundary values
 * still being significantly wrong.
 *
 * We avoid this problem in the global elimination process described above by
 * 'priming' the solution vector with the correct values for boundary
 * nodes. However, we can't do this for the local elimination
 * process. Therefore, if you experience a problem like the one above, you
 * need to either increase the diagonal entries in @p X to a size that matches
 * those in the other part of the Schur complement, or, simpler, prime the
 * solution vector before you start the solver.
 *
 * In conclusion, local elimination of boundary nodes only works if
 * there are no hanging nodes and even then doesn't always work fully
 * satisfactorily.
 *
 * @ingroup numerics
 * @author Wolfgang Bangerth, 1998, 2000, 2004, 2005
 */
namespace MatrixTools
{
  /**
   * Import namespace MatrixCreator for
   * backward compatibility with older
   * versions of deal.II in which these
   * namespaces were classes and class
   * MatrixTools was publicly derived from
   * class MatrixCreator.
   */
  using namespace MatrixCreator;

  /**
   * Apply Dirichlet boundary conditions
   * to the system matrix and vectors
   * as described in the general
   * documentation.
   */
  template <typename number>
  void
  apply_boundary_values (const std::map<types::global_dof_index,double> &boundary_values,
                         SparseMatrix<number>  &matrix,
                         Vector<number>        &solution,
                         Vector<number>        &right_hand_side,
                         const bool             eliminate_columns = true);

  /**
   * Apply Dirichlet boundary
   * conditions to the system
   * matrix and vectors as
   * described in the general
   * documentation. This function
   * works for block sparse
   * matrices and block vectors
   */
  template <typename number>
  void
  apply_boundary_values (const std::map<types::global_dof_index,double> &boundary_values,
                         BlockSparseMatrix<number>           &matrix,
                         BlockVector<number>                 &solution,
                         BlockVector<number>                 &right_hand_side,
                         const bool           eliminate_columns = true);

#ifdef DEAL_II_WITH_PETSC
  /**
   * Apply Dirichlet boundary conditions to
   * the system matrix and vectors as
   * described in the general
   * documentation. This function works on
   * the classes that are used to wrap
   * PETSc objects.
   *
   * Note that this function is not very
   * efficient: it needs to alternatingly
   * read and write into the matrix, a
   * situation that PETSc does not handle
   * too well. In addition, we only get rid
   * of rows corresponding to boundary
   * nodes, but the corresponding case of
   * deleting the respective columns
   * (i.e. if @p eliminate_columns is @p
   * true) is not presently implemented,
   * and probably will never because it is
   * too expensive without direct access to
   * the PETSc data structures. (This leads
   * to the situation where the action
   * indicated by the default value of the
   * last argument is actually not
   * implemented; that argument has
   * <code>true</code> as its default value
   * to stay consistent with the other
   * functions of same name in this class.)
   * A third reason against this function
   * is that it doesn't handle the case
   * where the matrix is distributed across
   * an MPI system.
   *
   * This function is used in
   * step-17 and
   * step-18.
   */
  void
  apply_boundary_values (const std::map<types::global_dof_index,double> &boundary_values,
                         PETScWrappers::SparseMatrix  &matrix,
                         PETScWrappers::Vector  &solution,
                         PETScWrappers::Vector  &right_hand_side,
                         const bool             eliminate_columns = true);

  /**
   * Same function, but for parallel PETSc
   * matrices.
   */
  void
  apply_boundary_values (const std::map<types::global_dof_index,double> &boundary_values,
                         PETScWrappers::MPI::SparseMatrix  &matrix,
                         PETScWrappers::MPI::Vector  &solution,
                         PETScWrappers::MPI::Vector  &right_hand_side,
                         const bool             eliminate_columns = true);

  /**
   * Same function, but for
   * parallel PETSc matrices. Note
   * that this function only
   * operates on the local range of
   * the parallel matrix, i.e. it
   * only eliminates rows
   * corresponding to degrees of
   * freedom for which the row is
   * stored on the present
   * processor. All other boundary
   * nodes are ignored, and it
   * doesn't matter whether they
   * are present in the first
   * argument to this function or
   * not. A consequence of this,
   * however, is that this function
   * has to be called from all
   * processors that participate in
   * sharing the contents of the
   * given matrices and vectors. It
   * is also implied that the local
   * range for all objects passed
   * to this function is the same.
   */
  void
  apply_boundary_values (const std::map<types::global_dof_index,double> &boundary_values,
                         PETScWrappers::MPI::SparseMatrix  &matrix,
                         PETScWrappers::Vector       &solution,
                         PETScWrappers::MPI::Vector  &right_hand_side,
                         const bool             eliminate_columns = true);

  /**
   * Same as above but for BlockSparseMatrix.
   */
  void
  apply_boundary_values (const std::map<types::global_dof_index,double>  &boundary_values,
                         PETScWrappers::MPI::BlockSparseMatrix &matrix,
                         PETScWrappers::MPI::BlockVector        &solution,
                         PETScWrappers::MPI::BlockVector        &right_hand_side,
                         const bool       eliminate_columns = true);

#endif

#ifdef DEAL_II_WITH_TRILINOS
  /**
   * Apply Dirichlet boundary
   * conditions to the system matrix
   * and vectors as described in the
   * general documentation. This
   * function works on the classes
   * that are used to wrap Trilinos
   * objects.
   *
   * Note that this function is not
   * very efficient: it needs to
   * alternatingly read and write
   * into the matrix, a situation
   * that Trilinos does not handle
   * too well. In addition, we only
   * get rid of rows corresponding to
   * boundary nodes, but the
   * corresponding case of deleting
   * the respective columns (i.e. if
   * @p eliminate_columns is @p true)
   * is not presently implemented,
   * and probably will never because
   * it is too expensive without
   * direct access to the Trilinos
   * data structures. (This leads to
   * the situation where the action
   * indicated by the default value
   * of the last argument is actually
   * not implemented; that argument
   * has <code>true</code> as its
   * default value to stay consistent
   * with the other functions of same
   * name in this class.)  A third
   * reason against this function is
   * that it doesn't handle the case
   * where the matrix is distributed
   * across an MPI system.
   */
  void
  apply_boundary_values (const std::map<types::global_dof_index,double> &boundary_values,
                         TrilinosWrappers::SparseMatrix  &matrix,
                         TrilinosWrappers::Vector        &solution,
                         TrilinosWrappers::Vector        &right_hand_side,
                         const bool             eliminate_columns = true);

  /**
   * This function does the same as
   * the one above, except now
   * working on block structures.
   */
  void
  apply_boundary_values (const std::map<types::global_dof_index,double> &boundary_values,
                         TrilinosWrappers::BlockSparseMatrix  &matrix,
                         TrilinosWrappers::BlockVector        &solution,
                         TrilinosWrappers::BlockVector        &right_hand_side,
                         const bool                eliminate_columns = true);

  /**
   * Apply Dirichlet boundary
   * conditions to the system matrix
   * and vectors as described in the
   * general documentation. This
   * function works on the classes
   * that are used to wrap Trilinos
   * objects.
   *
   * Note that this function is not
   * very efficient: it needs to
   * alternatingly read and write
   * into the matrix, a situation
   * that Trilinos does not handle
   * too well. In addition, we only
   * get rid of rows corresponding to
   * boundary nodes, but the
   * corresponding case of deleting
   * the respective columns (i.e. if
   * @p eliminate_columns is @p true)
   * is not presently implemented,
   * and probably will never because
   * it is too expensive without
   * direct access to the Trilinos
   * data structures. (This leads to
   * the situation where the action
   * indicated by the default value
   * of the last argument is actually
   * not implemented; that argument
   * has <code>true</code> as its
   * default value to stay consistent
   * with the other functions of same
   * name in this class.) This
   * function does work on MPI vector
   * types.
   */
  void
  apply_boundary_values (const std::map<types::global_dof_index,double> &boundary_values,
                         TrilinosWrappers::SparseMatrix  &matrix,
                         TrilinosWrappers::MPI::Vector   &solution,
                         TrilinosWrappers::MPI::Vector   &right_hand_side,
                         const bool             eliminate_columns = true);

  /**
   * This function does the same as
   * the one above, except now working
   * on block structures.
   */
  void
  apply_boundary_values (const std::map<types::global_dof_index,double> &boundary_values,
                         TrilinosWrappers::BlockSparseMatrix  &matrix,
                         TrilinosWrappers::MPI::BlockVector   &solution,
                         TrilinosWrappers::MPI::BlockVector   &right_hand_side,
                         const bool                eliminate_columns = true);
#endif

  /**
   * Rather than applying boundary
   * values to the global matrix
   * and vector after creating the
   * global matrix, this function
   * does so during assembly, by
   * modifying the local matrix and
   * vector contributions. If you
   * call this function on all
   * local contributions, the
   * resulting matrix will have the
   * same entries, and the final
   * call to
   * apply_boundary_values() on the
   * global system will not be
   * necessary.
   *
   * Since this function does not
   * have to work on the
   * complicated data structures of
   * sparse matrices, it is
   * relatively cheap. It may
   * therefore be a win if you have
   * many fixed degrees of freedom
   * (e.g. boundary nodes), or if
   * access to the sparse matrix is
   * expensive (e.g. for block
   * sparse matrices, or for PETSc
   * or trilinos
   * matrices). However, it doesn't
   * work as expected if there are
   * also hanging nodes to be
   * considered. More caveats are
   * listed in the general
   * documentation of this class.
   */
  void
  local_apply_boundary_values (const std::map<types::global_dof_index,double> &boundary_values,
                               const std::vector<types::global_dof_index> &local_dof_indices,
                               FullMatrix<double> &local_matrix,
                               Vector<double>     &local_rhs,
                               const bool          eliminate_columns);

  /**
   * Exception
   */
  DeclException0 (ExcBlocksDontMatch);
}



DEAL_II_NAMESPACE_CLOSE

#endif
