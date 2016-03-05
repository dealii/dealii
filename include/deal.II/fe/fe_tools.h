// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2016 by the deal.II authors
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

#ifndef dealii__fe_tools_H
#define dealii__fe_tools_H



#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/symmetric_tensor.h>

#include <vector>
#include <string>


DEAL_II_NAMESPACE_OPEN

template <typename number> class FullMatrix;
template <int dim> class Quadrature;
template <int dim, int spacedim> class FiniteElement;
template <int dim, int spacedim> class DoFHandler;
template <int dim> class FiniteElementData;
class ConstraintMatrix;



/*!@addtogroup feall */
/*@{*/


/**
 * This namespace offers interpolations and extrapolations of discrete
 * functions of one @p FiniteElement @p fe1 to another @p FiniteElement @p
 * fe2.
 *
 * It also provides the local interpolation matrices that interpolate on each
 * cell. Furthermore it provides the difference matrix $id-I_h$ that is needed
 * for evaluating $(id-I_h)z$ for e.g. the dual solution $z$.
 *
 * For more information about the <tt>spacedim</tt> template parameter check
 * the documentation of FiniteElement or the one of Triangulation.
 *
 * @author Wolfgang Bangerth, Ralf Hartmann, Guido Kanschat; 2000, 2003, 2004,
 * 2005, 2006
 */
namespace FETools
{
  /**
   * A base class for factory objects creating finite elements of a given
   * degree. Derived classes are called whenever one wants to have a
   * transparent way to create a finite element object.
   *
   * This class is used in the FETools::get_fe_from_name() and
   * FETools::add_fe_name() functions.
   *
   * @author Guido Kanschat, 2006
   */
  template <int dim, int spacedim=dim>
  class FEFactoryBase : public Subscriptor
  {
  public:
    /**
     * Create a FiniteElement and return a pointer to it.
     */
    virtual FiniteElement<dim,spacedim> *
    get (const unsigned int degree) const = 0;

    /**
     * Create a FiniteElement from a quadrature formula (currently only
     * implemented for FE_Q) and return a pointer to it.
     */

    virtual FiniteElement<dim,spacedim> *
    get (const Quadrature<1> &quad) const = 0;
    /**
     * Virtual destructor doing nothing but making the compiler happy.
     */
    virtual ~FEFactoryBase();
  };

  /**
   * A concrete class for factory objects creating finite elements of a given
   * degree.
   *
   * The class's get() function generates a finite element object of the type
   * given as template argument, and with the degree (however the finite
   * element class wishes to interpret this number) given as argument to
   * get().
   *
   * @author Guido Kanschat, 2006
   */
  template <class FE>
  class FEFactory : public FEFactoryBase<FE::dimension,FE::space_dimension>
  {
  public:
    /**
     * Create a FiniteElement and return a pointer to it.
     */
    virtual FiniteElement<FE::dimension,FE::space_dimension> *
    get (const unsigned int degree) const;

    /**
     * Create a FiniteElement from a quadrature formula (currently only
     * implemented for FE_Q) and return a pointer to it.
     */
    virtual FiniteElement<FE::dimension,FE::space_dimension> *
    get (const Quadrature<1> &quad) const;
  };

  /**
   * @warning In most cases, you will probably want to use
   * compute_base_renumbering().
   *
   * Compute the vector required to renumber the dofs of a cell by component.
   * Furthermore, compute the vector storing the start indices of each
   * component in the local block vector.
   *
   * The second vector is organized such that there is a vector for each base
   * element containing the start index for each component served by this base
   * element.
   *
   * While the first vector is checked to have the correct size, the second
   * one is reinitialized for convenience.
   */
  template<int dim, int spacedim>
  void compute_component_wise(
    const FiniteElement<dim,spacedim>                &fe,
    std::vector<unsigned int>               &renumbering,
    std::vector<std::vector<unsigned int> > &start_indices);

  /**
   * Compute the vector required to renumber the dofs of a cell by block.
   * Furthermore, compute the vector storing either the start indices or the
   * size of each local block vector.
   *
   * If the @p bool parameter is true, @p block_data is filled with the start
   * indices of each local block. If it is false, then the block sizes are
   * returned.
   *
   * The vector <tt>renumbering</tt> will be indexed by the standard numbering
   * of local degrees of freedom, namely first first vertex, then second
   * vertex, after vertices lines, quads, and hexes. For each index, the entry
   * indicates the index which this degree of freedom receives in a numbering
   * scheme, where the first block is numbered completely before the second.
   */
  template<int dim, int spacedim>
  void compute_block_renumbering (
    const FiniteElement<dim,spacedim>  &fe,
    std::vector<types::global_dof_index> &renumbering,
    std::vector<types::global_dof_index> &block_data,
    bool return_start_indices = true);

  /**
   * @name Generation of local matrices
   * @{
   */
  /**
   * Gives the interpolation matrix that interpolates a @p fe1- function to a
   * @p fe2-function on each cell. The interpolation_matrix needs to be of
   * size <tt>(fe2.dofs_per_cell, fe1.dofs_per_cell)</tt>.
   *
   * Note, that if the finite element space @p fe1 is a subset of the finite
   * element space @p fe2 then the @p interpolation_matrix is an embedding
   * matrix.
   */
  template <int dim, typename number, int spacedim>
  void
  get_interpolation_matrix(const FiniteElement<dim,spacedim> &fe1,
                           const FiniteElement<dim,spacedim> &fe2,
                           FullMatrix<number> &interpolation_matrix);

  /**
   * Gives the interpolation matrix that interpolates a @p fe1- function to a
   * @p fe2-function, and interpolates this to a second @p fe1-function on
   * each cell. The interpolation_matrix needs to be of size
   * <tt>(fe1.dofs_per_cell, fe1.dofs_per_cell)</tt>.
   *
   * Note, that this function only makes sense if the finite element space due
   * to @p fe1 is not a subset of the finite element space due to @p fe2, as
   * if it were a subset then the @p interpolation_matrix would be only the
   * unit matrix.
   */
  template <int dim, typename number, int spacedim>
  void
  get_back_interpolation_matrix(const FiniteElement<dim,spacedim> &fe1,
                                const FiniteElement<dim,spacedim> &fe2,
                                FullMatrix<number> &interpolation_matrix);

  /**
   * Gives the unit matrix minus the back interpolation matrix.  The @p
   * difference_matrix needs to be of size <tt>(fe1.dofs_per_cell,
   * fe1.dofs_per_cell)</tt>.
   *
   * This function gives the matrix that transforms a @p fe1 function $z$ to
   * $z-I_hz$ where $I_h$ denotes the interpolation operator from the @p fe1
   * space to the @p fe2 space. This matrix hence is useful to evaluate error-
   * representations where $z$ denotes the dual solution.
   */
  template <int dim, typename number, int spacedim>
  void
  get_interpolation_difference_matrix(const FiniteElement<dim,spacedim> &fe1,
                                      const FiniteElement<dim,spacedim> &fe2,
                                      FullMatrix<number> &difference_matrix);

  /**
   * Compute the local $L^2$-projection matrix from fe1 to fe2.
   */
  template <int dim, typename number, int spacedim>
  void get_projection_matrix(const FiniteElement<dim,spacedim> &fe1,
                             const FiniteElement<dim,spacedim> &fe2,
                             FullMatrix<number> &matrix);

  /**
   * Compute the matrix of nodal values of a finite element applied to all its
   * shape functions.
   *
   * This function is supposed to help building finite elements from
   * polynomial spaces and should be called inside the constructor of an
   * element. Applied to a completely initialized finite element, the result
   * should be the unit matrix by definition of the node values.
   *
   * Using this matrix allows the construction of the basis of shape functions
   * in two steps.
   *
   * <ol>
   *
   * <li>Define the space of shape functions using an arbitrary basis
   * <i>w<sub>j</sub></i> and compute the matrix <i>M</i> of node functionals
   * <i>N<sub>i</sub></i> applied to these basis functions.
   *
   * <li>Compute the basis <i>v<sub>j</sub></i> of the finite element shape
   * function space by applying <i>M<sup>-1</sup></i> to the basis
   * <i>w<sub>j</sub></i>.
   *
   * </ol>
   *
   * @note The FiniteElement must provide generalized support points and and
   * interpolation functions.
   */
  template <int dim, int spacedim>
  void compute_node_matrix(FullMatrix<double> &M,
                           const FiniteElement<dim,spacedim> &fe);

  /**
   * For all possible (isotropic and anisotropic) refinement cases compute the
   * embedding matrices from a coarse cell to the child cells. Each column of
   * the resulting matrices contains the representation of a coarse grid basis
   * function by the fine grid basis; the matrices are split such that there
   * is one matrix for every child.
   *
   * This function computes the coarse grid function in a sufficiently large
   * number of quadrature points and fits the fine grid functions using least
   * squares approximation. Therefore, the use of this function is restricted
   * to the case that the finite element spaces are actually nested.
   *
   * Note, that <code>matrices[refinement_case-1][child]</code> includes the
   * embedding (or prolongation) matrix of child <code>child</code> for the
   * RefinementCase <code>refinement_case</code>. Here, we use
   * <code>refinement_case-1</code> instead of <code>refinement_case</code> as
   * for RefinementCase::no_refinement(=0) there are no prolongation matrices
   * available.
   *
   * Typically this function is called by the various implementations of
   * FiniteElement classes in order to fill the respective
   * FiniteElement::prolongation matrices.
   *
   * @param fe The finite element class for which we compute the embedding
   * matrices.
   *
   * @param matrices A reference to RefinementCase<dim>::isotropic_refinement
   * vectors of FullMatrix objects. Each vector corresponds to one
   * RefinementCase @p refinement_case and is of the vector size
   * GeometryInfo<dim>::n_children(refinement_case). This is the format used
   * in FiniteElement, where we want to use this function mostly.
   *
   * @param isotropic_only Set to <code>true</code> if you only want to
   * compute matrices for isotropic refinement.
   *
   * @param threshold is the gap allowed in the least squares algorithm
   * computing the embedding.
   */
  template <int dim, typename number, int spacedim>
  void compute_embedding_matrices(const FiniteElement<dim,spacedim> &fe,
                                  std::vector<std::vector<FullMatrix<number> > > &matrices,
                                  const bool isotropic_only = false,
                                  const double threshold = 1.e-12);

  /**
   * Compute the embedding matrices on faces needed for constraint matrices.
   *
   * @param fe The finite element for which to compute these matrices.
   *
   * @param matrices An array of <i>GeometryInfo<dim>::subfaces_per_face =
   * 2<sup>dim-1</sup></i> FullMatrix objects,holding the embedding matrix for
   * each subface.
   *
   * @param face_coarse The number of the face on the coarse side of the face
   * for which this is computed.
   *
   * @param face_fine The number of the face on the refined side of the face
   * for which this is computed.
   *
   * @param threshold is the gap allowed in the least squares algorithm
   * computing the embedding.
   *
   * @warning This function will be used in computing constraint matrices. It
   * is not sufficiently tested yet.
   */
  template <int dim, typename number, int spacedim>
  void
  compute_face_embedding_matrices(const FiniteElement<dim,spacedim> &fe,
                                  FullMatrix<number> (&matrices)[GeometryInfo<dim>::max_children_per_face],
                                  const unsigned int face_coarse,
                                  const unsigned int face_fine,
                                  const double threshold = 1.e-12);

  /**
   * For all possible (isotropic and anisotropic) refinement cases compute the
   * <i>L<sup>2</sup></i>-projection matrices from the children to a coarse
   * cell.
   *
   * Note, that <code>matrices[refinement_case-1][child]</code> includes the
   * projection (or restriction) matrix of child <code>child</code> for the
   * RefinementCase <code>refinement_case</code>. Here, we use
   * <code>refinement_case-1</code> instead of <code>refinement_case</code> as
   * for RefinementCase::no_refinement(=0) there are no projection matrices
   * available.
   *
   * Typically this function is called by the various implementations of
   * FiniteElement classes in order to fill the respective
   * FiniteElement::restriction matrices.
   *
   * @arg fe The finite element class for which we compute the projection
   * matrices.  @arg matrices A reference to
   * <tt>RefinementCase<dim>::isotropic_refinement</tt> vectors of FullMatrix
   * objects. Each vector corresponds to one RefinementCase @p refinement_case
   * and is of the vector size
   * <tt>GeometryInfo<dim>::n_children(refinement_case)</tt>. This is the
   * format used in FiniteElement, where we want to use this function mostly.
   *
   * @arg isotropic_only Set to <code>true</code> if you only want to compute
   * matrices for isotropic refinement.
   */
  template <int dim, typename number, int spacedim>
  void compute_projection_matrices(
    const FiniteElement<dim,spacedim> &fe,
    std::vector<std::vector<FullMatrix<number> > > &matrices,
    const bool isotropic_only = false);

  /**
   * Projects scalar data defined in quadrature points to a finite element
   * space on a single cell.
   *
   * What this function does is the following: assume that there is scalar
   * data <tt>u<sub>q</sub>, 0 <= q < Q:=quadrature.size()</tt> defined at the
   * quadrature points of a cell, with the points defined by the given
   * <tt>rhs_quadrature</tt> object. We may then want to ask for that finite
   * element function (on a single cell) <tt>v<sub>h</sub></tt> in the finite-
   * dimensional space defined by the given FE object that is the projection
   * of <tt>u</tt> in the following sense:
   *
   * Usually, the projection <tt>v<sub>h</sub></tt> is that function that
   * satisfies <tt>(v<sub>h</sub>,w)=(u,w)</tt> for all discrete test
   * functions <tt>w</tt>. In the present case, we can't evaluate the right
   * hand side, since <tt>u</tt> is only defined in the quadrature points
   * given by <tt>rhs_quadrature</tt>, so we replace it by a quadrature
   * approximation. Likewise, the left hand side is approximated using the
   * <tt>lhs_quadrature</tt> object; if this quadrature object is chosen
   * appropriately, then the integration of the left hand side can be done
   * exactly, without any approximation. The use of different quadrature
   * objects is necessary if the quadrature object for the right hand side has
   * too few quadrature points -- for example, if data <tt>q</tt> is only
   * defined at the cell center, then the corresponding one-point quadrature
   * formula is obviously insufficient to approximate the scalar product on
   * the left hand side by a definite form.
   *
   * After these quadrature approximations, we end up with a nodal
   * representation <tt>V<sub>h</sub></tt> of <tt>v<sub>h</sub></tt> that
   * satisfies the following system of linear equations: <tt>M V<sub>h</sub> =
   * Q U</tt>, where <tt>M<sub>ij</sub>=(phi_i,phi_j)</tt> is the mass matrix
   * approximated by <tt>lhs_quadrature</tt>, and <tt>Q</tt> is the matrix
   * <tt>Q<sub>iq</sub>=phi<sub>i</sub>(x<sub>q</sub>) w<sub>q</sub></tt>
   * where <tt>w<sub>q</sub></tt> are quadrature weights; <tt>U</tt> is the
   * vector of quadrature point data <tt>u<sub>q</sub></tt>.
   *
   * In order to then get the nodal representation <tt>V<sub>h</sub></tt> of
   * the projection of <tt>U</tt>, one computes <tt>V<sub>h</sub> = X U,
   * X=M<sup>-1</sup> Q</tt>. The purpose of this function is to compute the
   * matrix <tt>X</tt> and return it through the last argument of this
   * function.
   *
   * Note that this function presently only supports scalar data. An extension
   * of the mass matrix is of course trivial, but one has to define the order
   * of data in the vector <tt>U</tt> if it contains vector valued data in all
   * quadrature points.
   *
   * A use for this function is described in the introduction to the step-18
   * example program.
   *
   * The opposite of this function, interpolation of a finite element function
   * onto quadrature points is essentially what the
   * <tt>FEValues::get_function_values</tt> functions do; to make things a
   * little simpler, the
   * <tt>FETools::compute_interpolation_to_quadrature_points_matrix</tt>
   * provides the matrix form of this.
   *
   * Note that this function works on a single cell, rather than an entire
   * triangulation. In effect, it therefore doesn't matter if you use a
   * continuous or discontinuous version of the finite element.
   *
   * It is worth noting that there are a few confusing cases of this function.
   * The first one is that it really only makes sense to project onto a finite
   * element that has at most as many degrees of freedom per cell as there are
   * quadrature points; the projection of N quadrature point data into a space
   * with M>N unknowns is well-defined, but often yields funny and non-
   * intuitive results. Secondly, one would think that if the quadrature point
   * data is defined in the support points of the finite element, i.e. the
   * quadrature points of <tt>ths_quadrature</tt> equal
   * <tt>fe.get_unit_support_points()</tt>, then the projection should be the
   * identity, i.e. each degree of freedom of the finite element equals the
   * value of the given data in the support point of the corresponding shape
   * function. However, this is not generally the case: while the matrix
   * <tt>Q</tt> in that case is the identity matrix, the mass matrix
   * <tt>M</tt> is not equal to the identity matrix, except for the special
   * case that the quadrature formula <tt>lhs_quadrature</tt> also has its
   * quadrature points in the support points of the finite element.
   *
   * Finally, this function only defines a cell wise projection, while one
   * frequently wants to apply it to all cells in a triangulation. However, if
   * it is applied to one cell after the other, the results from later cells
   * may overwrite nodal values computed already from previous cells if
   * degrees of freedom live on the interfaces between cells. The function is
   * therefore most useful for discontinuous elements.
   */
  template <int dim, int spacedim>
  void
  compute_projection_from_quadrature_points_matrix (const FiniteElement<dim,spacedim> &fe,
                                                    const Quadrature<dim>    &lhs_quadrature,
                                                    const Quadrature<dim>    &rhs_quadrature,
                                                    FullMatrix<double>       &X);

  /**
   * Given a (scalar) local finite element function, compute the matrix that
   * maps the vector of nodal values onto the vector of values of this
   * function at quadrature points as given by the second argument. In a
   * sense, this function does the opposite of the
   * FETools::compute_projection_from_quadrature_points_matrix function.
   */
  template <int dim, int spacedim>
  void
  compute_interpolation_to_quadrature_points_matrix (const FiniteElement<dim,spacedim> &fe,
                                                     const Quadrature<dim>    &quadrature,
                                                     FullMatrix<double>       &I_q);

  /**
   * Computes the projection of tensorial (first-order tensor) data stored at
   * the quadrature points @p vector_of_tensors_at_qp to data @p
   * vector_of_tensors_at_nodes at the support points of the cell.  The data
   * in @p vector_of_tensors_at_qp is ordered sequentially following the
   * quadrature point numbering.  The size of @p vector_of_tensors_at_qp must
   * correspond to the number of columns of @p projection_matrix.  The size of
   * @p vector_of_tensors_at_nodes must correspond to the number of rows of @p
   * vector_of_tensors_at_nodes .  The projection matrix @p projection_matrix
   * describes the projection of scalar data from the quadrature points and
   * can be obtained from the
   * FETools::compute_projection_from_quadrature_points_matrix function.
   */
  template <int dim>
  void
  compute_projection_from_quadrature_points(
    const FullMatrix<double>    &projection_matrix,
    const std::vector< Tensor<1, dim > >    &vector_of_tensors_at_qp,
    std::vector< Tensor<1, dim > >          &vector_of_tensors_at_nodes);



  /**
   * same as last function but for a @p SymmetricTensor .
   */
  template <int dim>
  void
  compute_projection_from_quadrature_points(
    const FullMatrix<double>    &projection_matrix,
    const std::vector< SymmetricTensor<2, dim > >   &vector_of_tensors_at_qp,
    std::vector< SymmetricTensor<2, dim > >         &vector_of_tensors_at_nodes);




  /**
   * This method implements the
   * FETools::compute_projection_from_quadrature_points_matrix method for
   * faces of a mesh.  The matrix that it returns, X, is face specific and its
   * size is fe.dofs_per_cell by rhs_quadrature.size().  The dimension, dim
   * must be larger than 1 for this class, since Quadrature<dim-1> objects are
   * required. See the documentation on the Quadrature class for more
   * information.
   */
  template <int dim, int spacedim>
  void
  compute_projection_from_face_quadrature_points_matrix (const FiniteElement<dim, spacedim> &fe,
                                                         const Quadrature<dim-1>    &lhs_quadrature,
                                                         const Quadrature<dim-1>    &rhs_quadrature,
                                                         const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
                                                         const unsigned int          face,
                                                         FullMatrix<double>         &X);



  //@}
  /**
   * @name Functions which should be in DoFTools
   */
  //@{
  /**
   * Gives the interpolation of a the @p dof1-function @p u1 to a @p
   * dof2-function @p u2. @p dof1 and @p dof2 need to be DoFHandlers based on
   * the same triangulation.
   *
   * If the elements @p fe1 and @p fe2 are either both continuous or both
   * discontinuous then this interpolation is the usual point interpolation.
   * The same is true if @p fe1 is a continuous and @p fe2 is a discontinuous
   * finite element. For the case that @p fe1 is a discontinuous and @p fe2 is
   * a continuous finite element there is no point interpolation defined at
   * the discontinuities.  Therefore the mean value is taken at the DoF values
   * on the discontinuities.
   *
   * Note that for continuous elements on grids with hanging nodes (i.e.
   * locally refined grids) this function does not give the expected output.
   * Indeed, the resulting output vector does not necessarily respect
   * continuity requirements at hanging nodes: if, for example, you are
   * interpolating a Q2 field to a Q1 field, then at hanging nodes the output
   * field will have the function value of the input field, which however is
   * not usually the mean value of the two adjacent nodes. It is thus not part
   * of the Q1 function space on the whole triangulation, although it is of
   * course Q1 on each cell.
   *
   * For this case (continuous elements on grids with hanging nodes), please
   * use the @p interpolate() function with an additional ConstraintMatrix
   * argument, see below, or make the field conforming yourself by calling the
   * @p distribute function of your hanging node constraints object.
   */
  template <int dim, int spacedim,
            template <int, int> class DoFHandlerType1,
            template <int, int> class DoFHandlerType2,
            class InVector, class OutVector>
  void
  interpolate (const DoFHandlerType1<dim,spacedim> &dof1,
               const InVector                      &u1,
               const DoFHandlerType2<dim,spacedim> &dof2,
               OutVector                           &u2);

  /**
   * Gives the interpolation of a the @p dof1-function @p u1 to a @p
   * dof2-function @p u2. @p dof1 and @p dof2 need to be DoFHandlers (or
   * hp::DoFHandlers) based on the same triangulation.  @p constraints is a
   * hanging node constraints object corresponding to @p dof2. This object is
   * particular important when interpolating onto continuous elements on grids
   * with hanging nodes (locally refined grids).
   *
   * If the elements @p fe1 and @p fe2 are either both continuous or both
   * discontinuous then this interpolation is the usual point interpolation.
   * The same is true if @p fe1 is a continuous and @p fe2 is a discontinuous
   * finite element. For the case that @p fe1 is a discontinuous and @p fe2 is
   * a continuous finite element there is no point interpolation defined at
   * the discontinuities.  Therefore the mean value is taken at the DoF values
   * at the discontinuities.
   */
  template <int dim, int spacedim,
            template <int, int> class DoFHandlerType1,
            template <int, int> class DoFHandlerType2,
            class InVector, class OutVector>
  void interpolate (const DoFHandlerType1<dim,spacedim> &dof1,
                    const InVector                      &u1,
                    const DoFHandlerType2<dim,spacedim> &dof2,
                    const ConstraintMatrix              &constraints,
                    OutVector                           &u2);

  /**
   * Gives the interpolation of the @p fe1-function @p u1 to a @p
   * fe2-function, and interpolates this to a second @p fe1-function named @p
   * u1_interpolated.
   *
   * Note, that this function does not work on continuous elements at hanging
   * nodes. For that case use the @p back_interpolate function, below, that
   * takes an additional @p ConstraintMatrix object.
   *
   * Furthermore note, that for the specific case when the finite element
   * space corresponding to @p fe1 is a subset of the finite element space
   * corresponding to @p fe2, this function is simply an identity mapping.
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  void back_interpolate (const DoFHandler<dim,spacedim>    &dof1,
                         const InVector                    &u1,
                         const FiniteElement<dim,spacedim> &fe2,
                         OutVector                         &u1_interpolated);

  /**
   * Same as last function, except that the dof handler objects might be of
   * type @p hp::DoFHandler.
   */
  template <int dim,
            template <int> class DoFHandlerType,
            class InVector, class OutVector, int spacedim>
  void back_interpolate (const DoFHandlerType<dim>         &dof1,
                         const InVector                    &u1,
                         const FiniteElement<dim,spacedim> &fe2,
                         OutVector                         &u1_interpolated);

  /**
   * Gives the interpolation of the @p dof1-function @p u1 to a @p
   * dof2-function, and interpolates this to a second @p dof1-function named
   * @p u1_interpolated.  @p constraints1 and @p constraints2 are the hanging
   * node constraints corresponding to @p dof1 and @p dof2, respectively.
   * These objects are particular important when continuous elements on grids
   * with hanging nodes (locally refined grids) are involved.
   *
   * Furthermore note, that for the specific case when the finite element
   * space corresponding to @p dof1 is a subset of the finite element space
   * corresponding to @p dof2, this function is simply an identity mapping.
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  void back_interpolate (const DoFHandler<dim,spacedim> &dof1,
                         const ConstraintMatrix         &constraints1,
                         const InVector                 &u1,
                         const DoFHandler<dim,spacedim> &dof2,
                         const ConstraintMatrix         &constraints2,
                         OutVector                      &u1_interpolated);

  /**
   * Gives $(Id-I_h)z_1$ for a given @p dof1-function $z_1$, where $I_h$ is
   * the interpolation from @p fe1 to @p fe2. The result $(Id-I_h)z_1$ is
   * written into @p z1_difference.
   *
   * Note, that this function does not work for continuous elements at hanging
   * nodes. For that case use the @p interpolation_difference function, below,
   * that takes an additional @p ConstraintMatrix object.
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  void interpolation_difference(const DoFHandler<dim,spacedim>    &dof1,
                                const InVector                    &z1,
                                const FiniteElement<dim,spacedim> &fe2,
                                OutVector                         &z1_difference);

  /**
   * Gives $(Id-I_h)z_1$ for a given @p dof1-function $z_1$, where $I_h$ is
   * the interpolation from @p fe1 to @p fe2. The result $(Id-I_h)z_1$ is
   * written into @p z1_difference.  @p constraints1 and @p constraints2 are
   * the hanging node constraints corresponding to @p dof1 and @p dof2,
   * respectively. These objects are particular important when continuous
   * elements on grids with hanging nodes (locally refined grids) are
   * involved.
   *
   * For parallel computations with PETSc, supply @p z1 with ghost elements
   * and @p z1_difference without ghost elements.
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  void interpolation_difference(const DoFHandler<dim,spacedim> &dof1,
                                const ConstraintMatrix         &constraints1,
                                const InVector                 &z1,
                                const DoFHandler<dim,spacedim> &dof2,
                                const ConstraintMatrix         &constraints2,
                                OutVector                      &z1_difference);



  /**
   * $L^2$ projection for discontinuous elements. Operates the same direction
   * as interpolate.
   *
   * The global projection can be computed by local matrices if the finite
   * element spaces are discontinuous. With continuous elements, this is
   * impossible, since a global mass matrix must be inverted.
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  void project_dg (const DoFHandler<dim,spacedim> &dof1,
                   const InVector                 &u1,
                   const DoFHandler<dim,spacedim> &dof2,
                   OutVector                      &u2);

  /**
   * Gives the patchwise extrapolation of a @p dof1 function @p z1 to a @p
   * dof2 function @p z2.  @p dof1 and @p dof2 need to be DoFHandler objects
   * based on the same triangulation. This function is used, for example, for
   * extrapolating patchwise a piecewise linear solution to a piecewise
   * quadratic solution.
   *
   * The function's name is historical and probably not particularly well
   * chosen. The function performs the following operations, one after the
   * other:
   *
   * - It interpolates directly from every cell of @p dof1 to the
   * corresponding cell of `dof2` using the interpolation matrix of the finite
   * element spaces used on these cells and provided by the finite element
   * objects involved. This step is done using the FETools::interpolate()
   * function. - It then performs a loop over all non-active cells of `dof2`.
   * If such a non-active cell has at least one active child, then we call the
   * children of this cell a "patch". We then interpolate from the children of
   * this patch to the patch, using the finite element space associated with
   * `dof2` and immediately interpolate back to the children. In essence, this
   * information throws away all information in the solution vector that lives
   * on a scale smaller than the patch cell. - Since we traverse non-active
   * cells from the coarsest to the finest levels, we may find patches that
   * correspond to child cells of previously treated patches if the mesh had
   * been refined adaptively (this cannot happen if the mesh has been refined
   * globally because there the children of a patch are all active). We also
   * perform the operation described above on these patches, but it is easy to
   * see that on patches that are children of previously treated patches, the
   * operation is now the identity operation (since it interpolates from the
   * children of the current patch a function that had previously been
   * interpolated to these children from an even coarser patch). Consequently,
   * this does not alter the solution vector any more.
   *
   * The name of the function originates from the fact that it can be used to
   * construct a representation of a function of higher polynomial degree on a
   * once coarser mesh. For example, if you imagine that you start with a
   * $Q_1$ function on globally refined mesh, and that @p dof2 is associated
   * with a $Q_2$ element, then this function computes the equivalent of the
   * operator $I_{2h}^{(2)}$ interpolating the original piecewise linear
   * function onto a quadratic function on a once coarser mesh with mesh size
   * $2h$ (but representing this function on the original mesh with size $h$).
   * If the exact solution is sufficiently smooth, then
   * $u^\ast=I_{2h}^{(2)}u_h$ is typically a better approximation to the exact
   * solution $u$ of the PDE than $u_h$ is. In other words, this function
   * provides a postprocessing step that improves the solution in a similar
   * way one often obtains by extrapolating a sequence of solutions,
   * explaining the origin of the function's name.
   *
   * @note The resulting field does not satisfy continuity requirements of the
   * given finite elements if the algorithm outlined above is used. When you
   * use continuous elements on grids with hanging nodes, please use the @p
   * extrapolate function with an additional ConstraintMatrix argument, see
   * below.
   *
   * @note Since this function operates on patches of cells, it requires that
   * the underlying grid is refined at least once for every coarse grid cell.
   * If this is not the case, an exception will be raised.
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  void extrapolate (const DoFHandler<dim,spacedim> &dof1,
                    const InVector                 &z1,
                    const DoFHandler<dim,spacedim> &dof2,
                    OutVector                      &z2);

  /**
   * Gives the patchwise extrapolation of a @p dof1 function @p z1 to a @p
   * dof2 function @p z2.  @p dof1 and @p dof2 need to be DoFHandler objects
   * based on the same triangulation.  @p constraints is a hanging node
   * constraints object corresponding to @p dof2. This object is necessary
   * when interpolating onto continuous elements on grids with hanging nodes
   * (locally refined grids).
   *
   * Otherwise, the function does the same as the other @p extrapolate
   * function above (for which the documentation provides an extensive
   * description of its operation).
   */
  template <int dim, class InVector, class OutVector, int spacedim>
  void extrapolate (const DoFHandler<dim,spacedim> &dof1,
                    const InVector                 &z1,
                    const DoFHandler<dim,spacedim> &dof2,
                    const ConstraintMatrix         &constraints,
                    OutVector                      &z2);
  //@}
  /**
   * The numbering of the degrees of freedom in continuous finite elements is
   * hierarchic, i.e. in such a way that we first number the vertex dofs, in
   * the order of the vertices as defined by the triangulation, then the line
   * dofs in the order and respecting the direction of the lines, then the
   * dofs on quads, etc. However, we could have, as well, numbered them in a
   * lexicographic way, i.e. with indices first running in x-direction, then
   * in y-direction and finally in z-direction. Discontinuous elements of
   * class FE_DGQ() are numbered in this way, for example.
   *
   * This function constructs a table which lexicographic index each degree of
   * freedom in the hierarchic numbering would have. It operates on the
   * continuous finite element given as first argument, and outputs the
   * lexicographic indices in the second.
   *
   * Note that since this function uses specifics of the continuous finite
   * elements, it can only operate on FiniteElementData<dim> objects inherent
   * in FE_Q(). However, this function does not take a FE_Q object as it is
   * also invoked by the FE_Q() constructor.
   *
   * It is assumed that the size of the output argument already matches the
   * correct size, which is equal to the number of degrees of freedom in the
   * finite element.
   */

  template <int dim>
  void
  hierarchic_to_lexicographic_numbering (unsigned int degree,
                                         std::vector<unsigned int> &h2l);

  template <int dim>
  void
  hierarchic_to_lexicographic_numbering (const FiniteElementData<dim> &fe_data,
                                         std::vector<unsigned int>    &h2l);

  /**
   * Like the previous function but instead of returning its result through
   * the last argument return it as a value.
   */
  template <int dim>
  std::vector<unsigned int>
  hierarchic_to_lexicographic_numbering (const FiniteElementData<dim> &fe_data);

  /**
   * This is the reverse function to the above one, generating the map from
   * the lexicographic to the hierarchical numbering. All the remarks made
   * about the above function are also valid here.
   */
  template <int dim>
  void
  lexicographic_to_hierarchic_numbering (const FiniteElementData<dim> &fe_data,
                                         std::vector<unsigned int>    &l2h);

  /**
   * Like the previous function but instead of returning its result through
   * the last argument return it as a value.
   */
  template <int dim>
  std::vector<unsigned int>
  lexicographic_to_hierarchic_numbering (const FiniteElementData<dim> &fe_data);

  /**
   * Parse the name of a finite element and generate a finite element object
   * accordingly. The parser ignores space characters between words (things
   * matching the regular expression [A-Za-z0-9_]).
   *
   * The name must be in the form which is returned by the
   * FiniteElement::get_name function, where dimension template parameters
   * &lt;2&gt; etc. can be omitted. Alternatively, the explicit number can be
   * replaced by <tt>dim</tt> or <tt>d</tt>. If a number is given, it
   * <b>must</b> match the template parameter of this function.
   *
   * The names of FESystem elements follow the pattern
   * <code>FESystem[FE_Base1^p1-FE_Base2^p2]</code> The powers <code>p1</code>
   * etc. may either be numbers or can be replaced by <tt>dim</tt> or
   * <tt>d</tt>.
   *
   *
   * If no finite element can be reconstructed from this string, an exception
   * of type @p FETools::ExcInvalidFEName is thrown.
   *
   * The function returns a pointer to a newly create finite element. It is in
   * the caller's responsibility to destroy the object pointed to at an
   * appropriate later time.
   *
   * Since the value of the template argument can't be deduced from the
   * (string) argument given to this function, you have to explicitly specify
   * it when you call this function.
   *
   * This function knows about all the standard elements defined in the
   * library. However, it doesn't by default know about elements that you may
   * have defined in your program. To make your own elements known to this
   * function, use the add_fe_name() function.  This function does not work if
   * one wants to get a codimension 1 finite element.
   */
  template <int dim, int spacedim>
  FiniteElement<dim, spacedim> *
  get_fe_by_name (const std::string &name);


  /**
   * @deprecated Use get_fe_by_name() with two template parameters instead
   */
  template <int dim>
  FiniteElement<dim,dim> *
  get_fe_from_name (const std::string &name);


  /**
   * Extend the list of finite elements that can be generated by
   * get_fe_from_name() by the one given as @p name. If get_fe_from_name() is
   * later called with this name, it will use the object given as second
   * argument to create a finite element object.
   *
   * The format of the @p name parameter should include the name of a finite
   * element. However, it is safe to use either the class name alone or to use
   * the result of FiniteElement::get_name (which includes the space dimension
   * as well as the polynomial degree), since everything after the first non-
   * name character will be ignored.
   *
   * The FEFactory object should be an object newly created with <tt>new</tt>.
   * FETools will take ownership of this object and delete it once it is not
   * used anymore.
   *
   * In most cases, if you want objects of type <code>MyFE</code> be created
   * whenever the name <code>my_fe</code> is given to get_fe_from_name, you
   * will want the second argument to this function be of type
   * FEFactory@<MyFE@>, but you can of course create your custom finite
   * element factory class.
   *
   * This function takes over ownership of the object given as second
   * argument, i.e. you should never attempt to destroy it later on. The
   * object will be deleted at the end of the program's lifetime.
   *
   * If the name of the element is already in use, an exception is thrown.
   * Thus, functionality of get_fe_from_name() can only be added, not changed.
   *
   * @note This function manipulates a global table (one table for each space
   * dimension). It is thread safe in the sense that every access to this
   * table is secured by a lock. Nevertheless, since each name can be added
   * only once, user code has to make sure that only one thread adds a new
   * element.
   *
   * Note also that this table exists once for each space dimension. If you
   * have a program that works with finite elements in different space
   * dimensions (for example,
   * @ref step_4 "step-4"
   * does something like this), then you should call this function for each
   * space dimension for which you want your finite element added to the map.
   */
  template <int dim, int spacedim>
  void add_fe_name (const std::string &name,
                    const FEFactoryBase<dim,spacedim> *factory);

  /**
   * The string used for get_fe_from_name() cannot be translated to a finite
   * element.
   *
   * Either the string is badly formatted or you are using a custom element
   * that must be added using add_fe_name() first.
   *
   * @ingroup Exceptions
   */
  DeclException1 (ExcInvalidFEName,
                  std::string,
                  << "Can't re-generate a finite element from the string '"
                  << arg1 << "'.");

  /**
   * The string used for get_fe_from_name() cannot be translated to a finite
   * element.
   *
   * Dimension arguments in finite element names should be avoided. If they
   * are there, the dimension should be <tt>dim</tt> or <tt>d</tt>. Here, you
   * gave a numeric dimension argument, which does not match the template
   * dimension of the finite element class.
   *
   * @ingroup Exceptions
   */
  DeclException2 (ExcInvalidFEDimension,
                  char, int,
                  << "The dimension " << arg1
                  << " in the finite element string must match "
                  << "the space dimension "
                  << arg2 << ".");

  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcInvalidFE);

  /**
   * The finite element must be
   * @ref GlossPrimitive "primitive".
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcFENotPrimitive);
  /**
   * Exception
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcTriangulationMismatch);

  /**
   * A continuous element is used on a mesh with hanging nodes, but the
   * constraint matrices are missing.
   *
   * @ingroup Exceptions
   */
  DeclException1 (ExcHangingNodesNotAllowed,
                  int,
                  << "You are using continuous elements on a grid with "
                  << "hanging nodes but without providing hanging node "
                  << "constraints. Use the respective function with "
                  << "additional ConstraintMatrix argument(s), instead."
                  << (arg1?"":""));
  /**
   * You need at least two grid levels.
   *
   * @ingroup Exceptions
   */
  DeclException0 (ExcGridNotRefinedAtLeastOnce);
  /**
   * The dimensions of the matrix used did not match the expected dimensions.
   *
   * @ingroup Exceptions
   */
  DeclException4 (ExcMatrixDimensionMismatch,
                  int, int, int, int,
                  << "This is a " << arg1 << "x" << arg2 << " matrix, "
                  << "but should be a " << arg3 << "x" << arg4 << " matrix.");

  /**
   * Exception thrown if an embedding matrix was computed inaccurately.
   *
   * @ingroup Exceptions
   */
  DeclException1(ExcLeastSquaresError, double,
                 << "Least squares fit leaves a gap of " << arg1);

  /**
   * Exception thrown if one variable may not be greater than another.
   *
   * @ingroup Exceptions
   */
  DeclException2 (ExcNotGreaterThan,
                  int,  int,
                  << arg1 << " must be greater than " << arg2);
}


#ifndef DOXYGEN

namespace FETools
{
  template <class FE>
  FiniteElement<FE::dimension, FE::space_dimension> *
  FEFactory<FE>::get (const unsigned int degree) const
  {
    return new FE(degree);
  }
}

#endif

/*@}*/

DEAL_II_NAMESPACE_CLOSE

/*----------------------------   fe_tools.h     ---------------------------*/
/* end of #ifndef dealii__fe_tools_H */
#endif
/*----------------------------   fe_tools.h     ---------------------------*/
