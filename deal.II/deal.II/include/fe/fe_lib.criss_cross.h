/*----------------------------   fe_lib.criss_cross.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __fe_lib_criss_cross_H
#define __fe_lib_criss_cross_H
/*----------------------------   fe_lib.criss_cross.h     ---------------------------*/


#include <fe/fe.h>
#include <base/quadrature.h>


/**
 * This class implements a rather unusual macro element, the so-called
 * criss-cross element. Its purpose is mostly to demonstrate the absence
 * of superconvergence effects on triangular meshes where at each vertex
 * more or less than six elements meet, but never exactly six.
 *
 * The construction of the element is best explained in 2d. Consider a
 * quadrilateral with basis functions at each vertex and one at the
 * crossing-point of the two diagonals. The element is divided by the
 * diagonals into four triangles and assume that each vertex basis
 * function has support only on the two triangles adjacent to the
 * respective vertex and is constant zero on the other two triangles;
 * they are linear on each of the triangles and globally continuous.
 * The center basis function lives on each of the four triangles, is
 * linear on each triangles and vanishes at the faces of the quadrilateral.
 *
 * Now, on the unit element, these basis functions are the same as for
 * a triangular trial space, namely the class of ${\cal P}_1$ Lagrange
 * elements. Due to the arrangement of the four triangles on the
 * quadrilateral, it is clear that considering the whole triangulation
 * of the domain, always four triangles meet at the points which
 * correspond with the centers of the quadrilaterals and $2*s$ triangles
 * meet at the vertices of the quadrilaterals, if $s$ is the number of
 * quadrilaterals meeting there. Thus, in most cases the number of
 * triangles meeting are four or eight, which effectively destroys
 * superconvergence at nodes.
 *
 * This element is not quite equivalent from beginning to the linear
 * triangular elements. The reason for this is that if we use a bilinear
 * mapping of the unit quadrilateral to the real cell, the diagonals will
 * in general not be straight lines. Therefore, the shape functions will
 * in general not be linear on the real cell, unlike for the linear
 * triangular element, which uses a linear mapping. The missing linearity
 * makes assemblage of matrices a bit more complicated, since the gradient
 * is not constant and we need more than one quadrature point, as well
 * as some other dubtle difficulties. This problem can, however, be cured
 * using a trick: the usual transformation from unit coordinates $\vec\xi$
 * to real coordinates $\vec x(\vec\xi)$ looks like
 * $$
 *   \vec x(\vec\xi) = \sum_{i=0}^3 \phi_i^L(\vec\xi) \vec x_i
 * $$
 * with $\phi_i^L$ being the bilinear basis functions associated with the
 * vertices and $\vec x_i$ being the coordinates of the vertices in real
 * space. Now, we could also choose
 * $$
 *   \vec x(\vec\xi) = \sum_{i=0}^4 \phi_i(\vec\xi) \vec x_i
 * $$
 * with the basis functions $\phi_i$ of this element, the four vertices
 * in real space $\vec x_0..\vec x_3$ and an interior point in real space
 * $\vec x_4$. We can choose the interior point quite arbitrarily and it
 * will become clear in a moment how we have to do so. First let us note
 * that because the vertex basis functions are linear on the faces,
 * because they vanish on the two faces not adjacent to the associated
 * vertex and because the center basis function vanishes at the four
 * faces, the four sides of the unit cell are mapped to straight lines
 * in real space, just like for the bilinear mapping.
 *
 * Now, to ensure that the mapping of each of the four triangles to the
 * real space is linear, we have to require that the two diagonals are
 * mapped to straight lines. One necessary condition for this is, that the
 * center point of the unit cell is mapped to the crossing point of the
 * two diagonals in real space. Therefore, we choose $\vec x_4$ to be
 * this point. Then we note, that because the vertex basis functions vanish
 * on the diagonal not through the vertex and are constant zero beyond that,
 * the mapping of the line from the center to a vertex is described entirely
 * by the basis function of that vertex and the center basis function; but
 * because they both are linear on that line, the line is also a straight
 * one in real space. This proves that by this construction of the mapping
 * between unit and real cell, the mapping of each of the four triangles
 * is linear (note that this does not hold for the whole element; the
 * mapping of the quadrilaterals is only piecewise linear and globally
 * continuous). It also proves that the trial space using this element
 * is equivalent to the trial space using triangles and linear elements.
 *
 * Since in one space dimension, this element equals two linear elements,
 * i.e. a linear ansatz on a mesh once more refined than the present one,
 * this element is not implemented for one dimension. There may be an
 * analogue to the criss-cross element in more than two space dimensions,
 * but it is not implemented at present.
 *
 * As stated above, the element is not really a good one. It may, however,
 * serve to study superconvergence effects and also to satisfy the author's
 * curiosity. At least for the first of these two reasons, it is better
 * suited than using a genuine triangulation of the domain (i.e. using real
 * triangles rather than subdividing quadrilaterals as shown above), since
 * the construction of triangulations with four or eight cells meeting at
 * each vertex is certainly not feasible other than by hand, while the
 * decomposition of a domain using quadrilaterals is easier.
 *
 *
 * \section{Hanging nodes}
 *
 * Hanging nodes are handled exactly like for any other element. It should
 * however be noted that the support of basis functions get quite
 * complicated in the presence of hanging nodes, as the following figure
 * depicts:
 * \begin{verbatim}
 *   *-----------------*--------*----
 *   |                /|\       |
 *   |              /..|.\      |
 *   |            /....|...\    |
 *   |          /......|.....\  |
 *   |         /.......|.......\|
 *   |       /.........*--------*----
 *   |      /..........|......./|
 *   |    /............|....../ |
 *   |   /.............|..../   |
 *   | /...............|.....\  |
 *   |/................|.......\|
 *   *-----------------o--------*-----
 * \end{verbatim}
 * The dotted area is the support of the basis function associated with the
 * bottom middle vertex (denoted by #o#) after the hanging node in the center
 * of the `picture' was eliminated. This strange structure of the support
 * should not pose too many problems in practice, it is only note here for
 * completeness and for curiosity.
 *
 *
 * \section{Experience with the criss-cross element}
 *
 * Experience is that the error for the criss cross element shows
 * the same convergence rate as the linear Lagrange element ($h^2$ for the
 * $L^2$ error, $h$ for the $H^1$ error). The $L^2$ error is about the same
 * size for the same number of elements as for the linear element; since
 * the criss-cross elements has about twice as many degrees of freedom as
 * the linear element for the same triangulation, the $L^2$ error really
 * is about twice as large as a function of the number of degrees of freedom.
 *
 * Converse to that, the $H^1$ error is about a factor of 1.2 smaller for
 * the same number of degrees of freedoms.
 *
 * Apart from all this data, it must not be forgotten that we cannot
 * expect superconvergence neither in the Gauss points not in the vertices.
 * Thus we may improve the accuracy of the solution obtained with the linear
 * element by a postprocess, while we can't do so for the criss-cross element.
 *
 * All given data refer to a Poisson equation with nonhomogeneous boundary
 * values on the unit disk (resp. a triangulation of that) and hanging nodes.
 *
 *
 * \section{Using quadrature formulae for this element}
 *
 * When using one of the usual quadrature formulae, a common problem is
 * that some of the quadrature points lie on the interfaces of the
 * triangles. For this reason, there is a family of quadrature formulae
 * defined below, names \ref{QCrissCross1} and higher order, which
 * resemble the quadrature formulae used on triangular domains, but
 * taken four-fold, i.e. for each of the four subtriangles.
 *
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class FECrissCross : public FiniteElement<dim> {
  public:
    				     /**
				      * Constructor
				      */
    FECrissCross ();

				     /**
				      * Return the value of the #i#th shape
				      * function at point #p# on the unit cell.
				      */
    virtual double shape_value(const unsigned int i,
			       const Point<dim>& p) const;

				     /**
				      * Return the gradient of the #i#th shape
				      * function at point #p# on the unit cell.
				      */
    virtual Tensor<1,dim> shape_grad(const unsigned int i,
				     const Point<dim>& p) const;

				     /**
				      * Return the tensor of second derivatives
				      * of the #i#th shape function at
				      * point #p# on the unit cell.
				      *
				      * The second derivatives are zero almost
				      * everywhere for this element; however,
				      * they are singular at the diagonals, so
				      * when trying to use this tensor, you
				      * should take special care and you may
				      * need to do some evaluation by hand.
				      */
    virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
					   const Point<dim>   &p) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_unit_support_points (vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_support_points (const DoFHandler<dim>::cell_iterator &cell,
				     const Boundary<dim> &boundary,
				     vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_support_points (const DoFHandler<dim>::face_iterator &face,
					  const Boundary<dim> &boundary,
					  vector<Point<dim> > &support_points) const;

    				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_local_mass_matrix (const DoFHandler<dim>::cell_iterator &cell,
					const Boundary<dim> &boundary,
					dFMatrix &local_mass_matrix) const;

    				     /**
				      * Return the value of the #i#th shape
				      * function at point #p# on the unit cell.
				      * Here, the (bi-)linear basis functions
				      * are meant, which are used for the
				      * computation of the transformation from
				      * unit cell to real space cell.
				      */
    virtual double shape_value_transform (const unsigned int i,
					  const Point<dim> &p) const;

				     /**
				      * Return the gradient of the #i#th shape
				      * function at point #p# on the unit cell.
				      * Here, the (bi-)linear basis functions
				      * are meant, which are used for the
				      * computation of the transformation from
				      * unit cell to real space cell.
				      */
    virtual Tensor<1,dim> shape_grad_transform (const unsigned int i,
						const Point<dim> &p) const;

    				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      *
				      * In two spatial dimensions, this function
				      * simply returns the length of the face.
				      */
    virtual void get_face_jacobians (const DoFHandler<dim>::face_iterator &face,
				     const Boundary<dim>         &boundary,
				     const vector<Point<dim-1> > &unit_points,
				     vector<double>      &face_jacobi_determinants) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      *
				      * In two spatial dimensions, this function
				      * simply returns half the length of the
				      * whole face.
				      */
    virtual void get_subface_jacobians (const DoFHandler<dim>::face_iterator &face,
					const unsigned int           subface_no,
					const vector<Point<dim-1> > &unit_points,
					vector<double>      &face_jacobi_determinants) const;

				     /**
				      * Return the normal vectors to the
				      * face with number #face_no# of #cell#.
				      *
				      * For linear finite elements, this function
				      * is particularly simple since all normal
				      * vectors are equal and can easiliy be
				      * computed from the direction of the face
				      * without using the transformation (Jacobi)
				      * matrix, at least for two dimensions.
				      *
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_normal_vectors (const DoFHandler<dim>::cell_iterator &cell,
				     const unsigned int          face_no,
				     const Boundary<dim>         &boundary,
				     const vector<Point<dim-1> > &unit_points,
				     vector<Point<dim> >         &normal_vectors) const;

				     /**
				      * Return the normal vectors to the
				      * subface with number #subface_no# of
				      * the face with number #face_no# of #cell#.
				      *
				      * For linear finite elements, this function
				      * is particularly simple since all normal
				      * vectors are equal and can easiliy be
				      * computed from the direction of the face
				      * without using the transformation (Jacobi)
				      * matrix, at least for two dimensions.
				      *
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_normal_vectors (const DoFHandler<dim>::cell_iterator &cell,
				     const unsigned int           face_no,
				     const unsigned int           subface_no,
				     const vector<Point<dim-1> > &unit_points,
				     vector<Point<dim> >         &normal_vectors) const;    

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      *
				      * For one dimensional elements, this
				      * function simply passes through to
				      * the one implemented in the base class.
				      * For higher dimensional finite elements
				      * we use linear mappings and therefore
				      * the boundary object is ignored since
				      * the boundary is approximated using
				      * piecewise multilinear boundary segments.
				      */
    virtual void fill_fe_values (const DoFHandler<dim>::cell_iterator &cell,
				 const vector<Point<dim> >            &unit_points,
				 vector<Tensor<2,dim> >               &jacobians,
				 const bool              compute_jacobians,
				 vector<Tensor<3,dim> > &jacobians_grad,
				 const bool              compute_jacobians_grad,
				 vector<Point<dim> >    &support_points,
				 const bool              compute_support_points,
				 vector<Point<dim> >    &q_points,
				 const bool              compute_q_points,
				 const dFMatrix         &shape_values_transform,
				 const vector<vector<Tensor<1,dim> > > &shape_grad_transform,
				 const Boundary<dim> &boundary) const;

    DeclException0 (ExcNotUseful);
};




/**
 * Quadrature formula for the criss-cross element. This quadrature
 * formula uses one point at the barycenter of each of the four subtriangles.
 *
 * For the same reason as for the criss-cross element itself, this
 * formula is not implemented for one space dimension.
 */
template <int dim>
class QCrissCross1 : public Quadrature<dim> {
  public:
    QCrissCross1 ();

    DeclException0 (ExcNotUseful);
};



/*----------------------------   fe_lib.criss_cross.h     ---------------------------*/
/* end of #ifndef __fe_lib_criss_cross_H */
#endif
/*----------------------------   fe_lib.criss_cross.h     ---------------------------*/
