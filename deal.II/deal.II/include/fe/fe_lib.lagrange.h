/*----------------------------   fe_lib.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __fe_lib_lagrange_H
#define __fe_lib_lagrange_H
/*----------------------------   fe_lib.h     ---------------------------*/


#include <fe/q1_mapping.h>

/**
 * Isoparametric Q1 finite element in #dim# space dimensions.
 *
 * The linear, isoparametric mapping from a point $\vec \xi$ on the unit cell
 * to a point $\vec x$ on the real cell is defined as
 * $$ \vec x(\vec \xi)  = \sum_j {\vec p_j} N_j(\xi) $$
 * where $\vec p_j$ is the vector to the $j$th corner point of the cell in
 * real space and $N_j(\vec \xi)$ is the value of the basis function associated
 * with the $j$th corner point, on the unit cell at point $\vec \xi$. The sum
 * over $j$ runs over all corner points.
 *
 * The number of degrees of freedom equal the number of the respective vertex
 * of the cell
 *
 * @author Wolfgang Bangerth, 1998, 1999
 */
template <int dim>
class FEQ1 : public FEQ1Mapping<dim>
{
  public:
				     /**
				      * Constructor
				      */
    FEQ1 ();
  protected:
				     /**
				      * Constructor that is called by the
				      * constructor of the derived
				      * #FEDG_Q1# class.
				      * It uses  no dofs in the vertices and 
				      * $2^d$ dofs per cell. No constraint
				      * matrices are build.
				      * For more detail see class #FEDGLinear#.
				      */
    FEQ1 (const int);
    
  public:
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
				      * For linear elements, all second
				      * derivatives on the unit cell are zero.
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
				     vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_support_points (const DoFHandler<dim>::face_iterator &face,
					  vector<Point<dim> > &support_points) const;

    				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      *
				      * Please note that as allowed in the
				      * documentation of the base class,
				      * this function does not implement
				      * the setting up of the local mass
				      * matrix in three space dimensions
				      * because of too high computational
				      * costs. The specified exception
				      * is thrown instead.
				      */
    virtual void get_local_mass_matrix (const DoFHandler<dim>::cell_iterator &cell,
					FullMatrix<double> &local_mass_matrix) const;

  private:

				     /**
				      * This function is simply singled out of
				      * the constructor; it sets up the
				      * #restriction# and #prolongation#
				      * matrices. Since we have two constructors
				      * which need this functionality, we
				      * provide a single function for this.
				      */
    void initialize_matrices ();
};




/**
 * Subparametric Q2 finite element in #dim# space dimensions.
 * A Q1 mapping from the unit cell
 * to the real cell is implemented.
 *
 * The numbering of the degrees of freedom is as follows:
 * \begin{itemize}
 * \item 1D case:
 *   \begin{verbatim}
 *      0---2---1
 *   \end{verbatim}
 *
 * \item 2D case:
 *   \begin{verbatim}
 *      3---6---2
 *      |       |
 *      7   8   5
 *      |       |
 *      0---4---1
 *   \end{verbatim}
 *
 * \item 3D case:
 *   \begin{verbatim}
 *         7--14---6        7--14---6
 *        /|       |       /       /|
 *      19 |       13     19      1813
 *      /  15      |     /       /  |
 *     3   |       |    3---10--2   |
 *     |   4--12---5    |       |   5
 *     |  /       /     |       9  /
 *    11 16      17     11      | 17
 *     |/       /       |       |/
 *     0---8---1        0---8---1
 *
 *         *-------*        *-------*
 *        /|       |       /       /|
 *       / |  21   |      /  24   / |
 *      /  |       |     /       /  |
 *     *   |       |    *-------*   |
 *     |25 *-------*    |       |23 *
 *     |  /       /     |   20  |  /
 *     | /  22   /      |       | /
 *     |/       /       |       |/
 *     *-------*        *-------* 
 *   \end{verbatim}
 *   The center vertex has number 26.
 *
 *   The respective coordinate values of the support points of the degrees
 *   of freedom are as follows:
 *   \begin{itemize}
 *   \item Index 0: #[0, 0, 0]#;
 *   \item Index 1: #[1, 0, 0]#;
 *   \item Index 2: #[1, 0, 1]#;
 *   \item Index 3: #[0, 0, 1]#;
 *   \item Index 4: #[0, 1, 0]#;
 *   \item Index 5: #[1, 1, 0]#;
 *   \item Index 6: #[1, 1, 1]#;
 *   \item Index 7: #[0, 1, 1]#;
 *   \item Index 8: #[1/2, 0, 0]#;
 *   \item Index 9: #[1, 0, 1/2]#;
 *   \item Index 10: # [1/2, 0, 1]#;
 *   \item Index 11: # [0, 0, 1/2]#;
 *   \item Index 12: # [1/2, 1, 0]#;
 *   \item Index 13: # [1, 1, 1/2]#;
 *   \item Index 14: # [1/2, 1, 1]#;
 *   \item Index 15: # [0, 1, 1/2]#;
 *   \item Index 16: # [0, 1/2, 0]#;
 *   \item Index 17: # [1, 1/2, 0]#;
 *   \item Index 18: # [1, 1/2, 1]#;
 *   \item Index 19: # [0, 1/2, 1]#;
 *   \item Index 20: # [1/2, 0, 1/2]#;
 *   \item Index 21: # [1/2, 1, 1/2]#;
 *   \item Index 22: # [1/2, 1/2, 0]#;
 *   \item Index 23: # [1, 1/2, 1/2]#;
 *   \item Index 24: # [1/2, 1/2, 1]#;
 *   \item Index 25: # [0, 1/2, 1/2]#;
 *   \item Index 26: # [1/2, 1/2, 1/2]#; 
 *   \end{itemize}
 * \end{itemize}
 *
 * @author Wolfgang Bangerth, 1998, 1999
 */
template <int dim>
class FEQ2 : public FEQ1Mapping<dim>
{
  public:
				     /**
				      * Constructor
				      */
    FEQ2 ();
  protected:
				     /**
				      * Constructor that is called by the
				      * constructor of the derived
				      * #FEDGQuadraticSub# class.
				      * It uses no dofs in the vertices, no
				      * dofs in the lines and 
				      * $3^d$ dofs per cell. No constraint
				      * matrices are build.
				      * For more detail see class 
				      * #FEDGQuadraticSub#.
				      */
    FEQ2 (const int);
    
  public:
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
				     vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_support_points (const DoFHandler<dim>::face_iterator &face,
					  vector<Point<dim> > &support_points) const;

    				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      *
				      * Please note that as allowed in the
				      * documentation of the base class,
				      * this function does not implement
				      * the setting up of the local mass
				      * matrix in three space dimensions
				      * because of too high computational
				      * costs. The specified exception
				      * is thrown instead.
				      */
    virtual void get_local_mass_matrix (const DoFHandler<dim>::cell_iterator &cell,
					FullMatrix<double> &local_mass_matrix) const;

  private:

				     /**
				      * This function is simply singled out of
				      * the constructor; it sets up the
				      * #restriction# and #prolongation#
				      * matrices. Since we have two constructors
				      * which need this functionality, we
				      * provide a single function for this.
				      */
    void initialize_matrices ();
};




/**
 * Subparametric Q3 finite element in #dim# space dimensions.
 * A Q1 mapping from the unit cell
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
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class FEQ3 : public FEQ1Mapping<dim>
{
  public:
				     /**
				      * Constructor
				      */
    FEQ3 ();

  protected:
				     /**
				      * Constructor that is called by the
				      * constructor of the derived
				      * #FEDGCubicSub# class.
				      * It uses  no dofs in the vertices and 
				      * $4^d$ dofs per cell. No constraint
				      * matrices are build.
				      * For more detail see class
				      * #FEDGCubicSub#.
				      */
    FEQ3 (const int);

  public:
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
				     vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_support_points (const DoFHandler<dim>::face_iterator &face,
					  vector<Point<dim> > &support_points) const;

    				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      *
				      * Please note that as allowed in the
				      * documentation of the base class,
				      * this function does not implement
				      * the setting up of the local mass
				      * matrix in three space dimensions
				      * because of too high computational
				      * costs. The specified exception
				      * is thrown instead.
				      */
    virtual void get_local_mass_matrix (const DoFHandler<dim>::cell_iterator &cell,
					FullMatrix<double> &local_mass_matrix) const;

  private:

				     /**
				      * This function is simply singled out of
				      * the constructor; it sets up the
				      * #restriction# and #prolongation#
				      * matrices. Since we have two constructors
				      * which need this functionality, we
				      * provide a single function for this.
				      */
    void initialize_matrices ();
};



/**
 * Subparametric Q4 finite element in #dim# space dimensions.
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
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class FEQ4 : public FEQ1Mapping<dim>
{
  public:
				     /**
				      * Constructor
				      */
    FEQ4 ();

  protected:
				     /**
				      * Constructor that is called by the
				      * constructor of the derived
				      * #FEDGQuarticSub# class.
				      * It uses  no dofs in the vertices and 
				      * $5^d$ dofs per cell. No constraint
				      * matrices are build.
				      * For more detail see class
				      * #FEDGQuarticSub#.
				      */
    FEQ4 (const int);

  public:
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
				     vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_support_points (const DoFHandler<dim>::face_iterator &face,
					  vector<Point<dim> > &support_points) const;

    				     /**
				      * Refer to the base class for detailed
				      * information on what this function does.
				      *
				      * Please note that as allowed in the
				      * documentation of the base class,
				      * this function does not implement
				      * the setting up of the local mass
				      * matrix in three space dimensions
				      * because of too high computational
				      * costs. The specified exception
				      * is thrown instead.
				      */
    virtual void get_local_mass_matrix (const DoFHandler<dim>::cell_iterator &cell,
					FullMatrix<double> &local_mass_matrix) const;

  private:

				     /**
				      * This function is simply singled out of
				      * the constructor; it sets up the
				      * #restriction# and #prolongation#
				      * matrices. Since we have two constructors
				      * which need this functionality, we
				      * provide a single function for this.
				      */
    void initialize_matrices ();
};








/*----------------------------   fe_lib.h     ---------------------------*/
/* end of #ifndef __fe_lib_lagrange_H */
#endif
/*----------------------------   fe_lib.h     ---------------------------*/


