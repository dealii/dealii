/*----------------------------   fe_lib.h     ---------------------------*/
/*      $Id$                 */
#ifndef __fe_lib_H
#define __fe_lib_H
/*----------------------------   fe_lib.h     ---------------------------*/


#include <fe/fe.h>



/**
 * Define a (bi-, tri-, etc)linear finite element in #dim# space dimensions,
 * along with (bi-, tri-)linear (therefore isoparametric) transforms from the
 * unit cell to the real cell.
 *
 * The linear, isoparametric mapping from a point $\vec \xi$ on the unit cell
 * to a point $\vec x$ on the real cell is defined as
 * $$ \vec x(\vec \xi)  = \sum_j {\vec p_j} N_j(\xi) $$
 * where $\vec p_j$ is the vector to the $j$th corner point of the cell in
 * real space and $N_j(\vec \xi)$ is the value of the basis function associated
 * with the $j$th corner point, on the unit cell at point $\vec \xi$. The sum
 * over $j$ runs over all corner points.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class FELinear : public FELinearMapping<dim> {
  public:
				     /**
				      * Constructor
				      */
    FELinear ();

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
    virtual Point<dim> shape_grad(const unsigned int i,
				  const Point<dim>& p) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_unit_ansatz_points (vector<Point<dim> > &ansatz_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_ansatz_points (const DoFHandler<dim>::cell_iterator &cell,
				    const Boundary<dim> &boundary,
				    vector<Point<dim> > &ansatz_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_ansatz_points (const DoFHandler<dim>::face_iterator &face,
					 const Boundary<dim> &boundary,
					 vector<Point<dim> > &ansatz_points) const;

    				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_local_mass_matrix (const DoFHandler<dim>::cell_iterator &cell,
					const Boundary<dim> &boundary,
					dFMatrix &local_mass_matrix) const;
};




/**
 * Define a (bi-, tri-, etc)quadratic finite element in #dim# space dimensions.
 * A linear (subparametric) mapping from the unit cell
 * to the real cell is implemented.
 */
template <int dim>
class FEQuadraticSub : public FELinearMapping<dim> {
  public:
				     /**
				      * Constructor
				      */
    FEQuadraticSub ();

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
    virtual Point<dim> shape_grad(const unsigned int i,
				  const Point<dim>& p) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_unit_ansatz_points (vector<Point<dim> > &ansatz_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_ansatz_points (const DoFHandler<dim>::cell_iterator &cell,
				    const Boundary<dim> &boundary,
				    vector<Point<dim> > &ansatz_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_ansatz_points (const DoFHandler<dim>::face_iterator &face,
					 const Boundary<dim> &boundary,
					 vector<Point<dim> > &ansatz_points) const;

    				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_local_mass_matrix (const DoFHandler<dim>::cell_iterator &cell,
					const Boundary<dim> &boundary,
					dFMatrix &local_mass_matrix) const;
};




/**
 * Define a (bi-, tri-, etc)cubic finite element in #dim# space dimensions.
 * A linear (subparametric) mapping from the unit cell
 * to the real cell is implemented.
 *
 * The numbering of degrees of freedom in one spatial dimension is as follows:
 * \begin{verbatim}
 *   0--2--3--1
 * \end{verbatim}
 *
 * The numbering of degrees of freedom in two spatial dimension is as follows:
 * \begin{verbatim}
 *   3--8--9--2
 *   |        |
 *   11 15 14 7
 *   |        |
 *   10 12 13 6
 *   |        |
 *   0--4--5--1
 * \end{verbatim}
 * Note the reverse ordering of degrees of freedom on the left and upper
 * line and the counterclockwise numbering of the interior degrees of
 * freedom.
 */
template <int dim>
class FECubicSub : public FELinearMapping<dim> {
  public:
				     /**
				      * Constructor
				      */
    FECubicSub ();

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
    virtual Point<dim> shape_grad(const unsigned int i,
				  const Point<dim>& p) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_unit_ansatz_points (vector<Point<dim> > &ansatz_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_ansatz_points (const DoFHandler<dim>::cell_iterator &cell,
				    const Boundary<dim> &boundary,
				    vector<Point<dim> > &ansatz_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_ansatz_points (const DoFHandler<dim>::face_iterator &face,
					 const Boundary<dim> &boundary,
					 vector<Point<dim> > &ansatz_points) const;

    				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_local_mass_matrix (const DoFHandler<dim>::cell_iterator &cell,
					const Boundary<dim> &boundary,
					dFMatrix &local_mass_matrix) const;
};



/**
 * Define a (bi-, tri-, etc)quartic finite element in #dim# space dimensions.
 * A linear (subparametric) mapping from the unit cell
 * to the real cell is implemented.
 *
 * The numbering of degrees of freedom in one spatial dimension is as follows:
 * \begin{verbatim}
 *   0--2--3--4--1
 * \end{verbatim}
 *
 * The numbering of degrees of freedom in two spatial dimension is as follows:
 * \begin{verbatim}
 *   3--10-11-12-2
 *   |           |
 *   15 19 22 18 9
 *   |           |
 *   14 23 24 21 8
 *   |           |
 *   13 16 20 17 7
 *   |           |
 *   0--4--5--6--1
 * \end{verbatim}
 * Note the reverse ordering of degrees of freedom on the left and upper
 * line and the numbering of the interior degrees of
 * freedom.
 */
template <int dim>
class FEQuarticSub : public FELinearMapping<dim> {
  public:
				     /**
				      * Constructor
				      */
    FEQuarticSub ();

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
    virtual Point<dim> shape_grad(const unsigned int i,
				  const Point<dim>& p) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_unit_ansatz_points (vector<Point<dim> > &ansatz_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_ansatz_points (const DoFHandler<dim>::cell_iterator &cell,
				    const Boundary<dim> &boundary,
				    vector<Point<dim> > &ansatz_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_ansatz_points (const DoFHandler<dim>::face_iterator &face,
					 const Boundary<dim> &boundary,
					 vector<Point<dim> > &ansatz_points) const;

    				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_local_mass_matrix (const DoFHandler<dim>::cell_iterator &cell,
					const Boundary<dim> &boundary,
					dFMatrix &local_mass_matrix) const;
};








/*----------------------------   fe_lib.h     ---------------------------*/
/* end of #ifndef __fe_lib_H */
#endif
/*----------------------------   fe_lib.h     ---------------------------*/


