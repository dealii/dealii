/*----------------------------   fe_lib.system.h     ---------------------------*/
/*      $Id$                 */
#ifndef __deal_fe_system_H
#define __deal_fe_system_H
/*----------------------------   fe_lib.system.h     ---------------------------*/


#include <fe/fe.h>
#include <vector>
#include <pair.h>



/**
 * This class provides an interface to group several equal elements together
 * into one. To the outside world, the resulting object looks just like
 * a usual finite element object, which is composed of several other finite
 * elements of the same class each.
 *
 * Basically, this composed finite element has #N# times as many degrees of
 * freedom (and therefore also #N# times as many shape functions) as a single
 * object of the underlying finite element would have had. Among these,
 * always #N# have the same properties, i.e. are represented by the same
 * shape functions. These #N# shape functions for each degree of freedom
 * of the basic finite element are numbered consecutively, i.e. for
 * the common case of a velocity #(u,v,w)#, the sequence of basis functions
 * will be #u1, v1, w1, u2, v2, w2, ..., uN, vN, wN# compared to
 *  #u1, ..., uN, v1, ..., vN, w1, ...wN#.
 *
 * Using this scheme, the overall numbering of degrees of freedom is as
 * follows: for each subobject (vertex, line, quad, or hex), the degrees
 * of freedom are numbered such that we run over all subelements first,
 * before turning for the next dof on this subobject or for the next subobject.
 * For example, for the bicubic element in one space dimension, and for
 * two subobjects grouped together by this class, the ordering for
 * the system #s=(u,v)# is:
 * \begin{itemize}
 * \item First vertex: #u0, v0 = s0, s1#
 * \item Second vertex: #u1, v1 = s2, s3#
 * \item First degree of freedom on the line (=cell):
 *   #u2, v2 = s3, s4#
 * \item Second degree of freedom on the line:
 *   #u3, v3 = s5, s6#.
 * \end{itemize}
 *
 * In the most cases, the composed element behaves as if it were a usual element
 * with more degrees of freedom. However the underlying structure is visible in
 * the restriction, prolongation and interface constraint matrices, which do not
 * couple the degrees of freedom of the subobject. E.g. the continuity requirement
 * is imposed for the shape functions of the subobjects separately; no requirement
 * exist between shape functions of different subobjects, i.e. in the above
 * example: on a hanging node, the respective value of the #u# velocity is only
 * coupled to #u# at the vertices and the line on the larger cell next to this
 * vertex, there is no interaction with #v# and #w# of this or the other cell.
 *
 * Likewise, the matrix computed by the #get_local_mass_matrix# function, which
 * originally is defined to be $m_{ij} = \int_K \phi_i \phi_j dx$ contains
 * only those $m_{ij}$ for which the respective shape functions belong to the
 * same subobject, all other entries are set to zero. The matrix therefore is
 * a block matrix, where each block is a diagonal matrix with entries equal to
 * the entry at this block's position in the local mass matrix of a single
 * finite element object. This behaviour is consistent with one common use
 * of the mass matrix, which is in projecting functions onto the grid; in this
 * case, one wants to project each component of the function (here it is a vector
 * function) to the respective component of the finite element, without interaction
 * of the different components.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1999
 */
template <int dim>
class FESystem //<dim>
  : public FiniteElement<dim>
{
				     /**
				      * Copy constructor prohibited.
				      */
    FESystem(const FESystem<dim>&);

  public:

				     /**
				      * Constructor. Take a finite element type
				      * and the number of elements you want to
				      * group together using this class.
				      *
				      * In fact, the object #fe# is not used,
				      * apart from getting the number of dofs
				      * per vertex, line, etc for that finite
				      * element class. The objects creates its
				      * own copy of the finite element object
				      * at construction time (but after
				      * the initialization of the base class
				      * #FiniteElement#, which is why we need
				      * a valid finite element object passed
				      * to the constructor).
				      *
				      * Obviously, the template finite element
				      * class needs to be of the same dimension
				      * as is this object.
				      */
    template <class FE>
    FESystem (const FE &fe, const unsigned int n_elements);

				     /** 
				      * Constructor for mixed
				      * discretizations with two
				      * base elements.
				      *
				      * See the other constructor.
				      */
    template <class FE1, class FE2>
    FESystem (const FE1 &fe1, const unsigned int n1,
	      const FE2 &fe2, const unsigned int n2);

    				     /** 
				      * Constructor for mixed
				      * discretizations with three
				      * base elements.
				      *
				      * See the other constructor.
				      */
    template <class FE1, class FE2, class FE3>
    FESystem (const FE1 &fe1, const unsigned int n1,
	      const FE2 &fe2, const unsigned int n2,
	      const FE3 &fe3, const unsigned int n3);

				     /**
				      * Destructor.
				      */
    virtual ~FESystem ();

    				     /**
				      * Return the value of the #i#th shape
				      * function at point #p# on the unit cell.
				      *
				      * For an element composed of #N#
				      * subelements, the first #N# shape
				      * functions refer to the zeroth shape
				      * function of the underlying object,
				      * the shape functions #N..2N-1# refer
				      * to the base shape function with
				      * number #1#, and so on. The #i# shape
				      * function therefore equals the
				      * #i/N# the shape function of the
				      * base object.
				      */
    virtual double shape_value(const unsigned int i,
			       const Point<dim>  &p) const;

				     /**
				      * Return the gradient of the #i#th shape
				      * function at point #p# on the unit cell.
				      *
				      * For the ordering of shape functions
				      * refer to the #shape_value# function.
				      */
    virtual Tensor<1,dim> shape_grad(const unsigned int i,
				     const Point<dim>& p) const;

				     /**
				      * Return the tensor of second derivatives
				      * of the #i#th shape function at
				      * point #p# on the unit cell.
				      *
				      * For the ordering of shape functions
				      * refer to the #shape_value# function.
				      */
    virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
					   const Point<dim>   &p) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      *
				      * For the ordering of shape functions
				      * refer to the #shape_value# function.
				      */
    virtual void get_unit_support_points (vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      *
				      * For the ordering of shape functions
				      * refer to the #shape_value# function.
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
				      * Fill the local mass matrix. The elements
				      * of this matrix are the integrals
				      * $\int_K \phi_i \phi_j dx$ over a given
				      * cell $K$. However, here only those
				      * elements of the matrix are set for which
				      * the shape functions $\phi_i$ and
				      * $\phi_j$ belong to the same subelement,
				      * i.e. the resulting matrix is a block
				      * matrix where each block is a diagonal
				      * matrix with diagonal values equal to
				      * the respective entry of the local mass
				      * matrix for the underlying finite element
				      * class. This definition of the mass
				      * matrix for systems of finite elements
				      * is consistent with the use of the matrix
				      * for the projection of initial values and
				      * the like, where the subelements are not
				      * coupled to each other. Also in most
				      * other cases you will not want the
				      * coupling terms to appear in the mass
				      * matrix.
				      *
				      * If the shape functions of this element
				      * were numbered such that the first
				      * numbers are for the shape functions of
				      * the first subelement, then those for
				      * the second subelement, and so on, then
				      * the mass matrix generated by this
				      * function would be a block diagonal
				      * matrix with each block being the mass
				      * matrix of the base finite element.
				      *
				      * Refer to the base class for more
				      * information on this function.
				      */
    virtual void get_local_mass_matrix (const DoFHandler<dim>::cell_iterator &cell,
					FullMatrix<double> &local_mass_matrix) const;

				     /**
				      * Return the value of the #i#th shape
				      * function of the transformation mapping
				      * from unit cell to real cell. Since
				      * the transform functions are not
				      * touched when clustering several finite
				      * element objects together using this
				      * class, this function simply passes down
				      * the call to the respective function of
				      * the underlying element.
				      */
    virtual double shape_value_transform (const unsigned int i,
					  const Point<dim> &p) const;

				     /**
				      * Same as above: return gradient of the
				      * #i#th shape function for the mapping
				      * from unit to real cell.
				      */
    virtual Tensor<1,dim> shape_grad_transform (const unsigned int i,
						const Point<dim> &p) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      *
				      * Since this function is only about the
				      * mapping from unit to real cell, it
				      * is not affected by putting several
				      * equal elements together, so this
				      * function simply passes down to the
				      * underlying object.
				      */
    virtual void get_face_jacobians (const DoFHandler<dim>::face_iterator &face,
				     const vector<Point<dim-1> > &unit_points,
				     vector<double>      &face_jacobi_determinants) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      *
				      * Since this function is only about the
				      * mapping from unit to real cell, it
				      * is not affected by putting several
				      * equal elements together, so this
				      * function simply passes down to the
				      * underlying object.
				      */
    virtual void get_subface_jacobians (const DoFHandler<dim>::face_iterator &face,
					const unsigned int           subface_no,
					const vector<Point<dim-1> > &unit_points,
					vector<double>      &face_jacobi_determinants) const;

				     /**
				      * Return the normal vectors to the
				      * face with number #face_no# of #cell#.
				      *
				      * Since this function is only about the
				      * mapping from unit to real cell, it
				      * is not affected by putting several
				      * equal elements together, so this
				      * function simply passes down to the
				      * underlying object.
				      *
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_normal_vectors (const DoFHandler<dim>::cell_iterator &cell,
				     const unsigned int          face_no,
				     const vector<Point<dim-1> > &unit_points,
				     vector<Point<dim> >         &normal_vectors) const;

				     /**
				      * Return the normal vectors to the
				      * subface with number #subface_no# of
				      * the face with number #face_no# of #cell#.
				      *
				      * Since this function is only about the
				      * mapping from unit to real cell, it
				      * is not affected by putting several
				      * equal elements together, so this
				      * function simply passes down to the
				      * underlying object.
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
				      * Implementation of the corresponding function of #FiniteElement#.
				      */
    virtual void fill_fe_values (const DoFHandler<dim>::cell_iterator &cell,
				 const vector<Point<dim> >            &unit_points,
				 vector<Tensor<2,dim> >               &jacobians,
				 const bool              compute_jacobians,
				 vector<Tensor<3,dim> > &jacobians_grad,
				 const bool              compute_jacobians_grad,
				 vector<Point<dim> > &support_points,
				 const bool           compute_support_points,
				 vector<Point<dim> > &q_points,
				 const bool           compute_q_points,
				 const FullMatrix<double>  &shape_values_transform,
				 const vector<vector<Tensor<1,dim> > > &shape_grad_transform) const;
    
				     /** 
				      * Number of different base
				      * elements of this object.
				      *
				      * Since these objects can have
				      * multiplicity and subobjects
				      * themselves, this may be
				      * smaller than the total number
				      * of finite elements composed
				      * into this structure.
				      */
    virtual unsigned n_base_elements() const;

				     /**
				      * How often is a composing element used.
				      *
				      */
    unsigned element_multiplicity(unsigned index) const;

				     /**
				      * Access to a composing element.
				      *
				      * If you assemble your system
				      * matrix, you usually will not
				      * want to have an FEValues object
				      * with a lot of equal entries. Ok,
				      * so initialize your FEValues with
				      * the #base_element# you get by
				      * this function. In a mixed
				      * discretization, you can choose
				      * the different base element types
				      * by index.
				      *
				      */
    virtual const FiniteElement<dim>& base_element(unsigned index) const;

  private:

				     /**
				      * Pairs of multiplicity and element type.
				      */
    typedef pair<const FiniteElement<dim> *, unsigned > ElementPair;
    
				     /**
				      * Pointer to underlying finite
				      * element classes.
				      *
				      * This object contains a pointer
				      * to each contributing element
				      * of a mixed discretization and
				      * its multiplicity. It is
				      * created by the constructor and
				      * constant afterwards.
				      */
    vector< ElementPair > base_elements;

				     /**
				      * The base element establishing a
				      * component.
				      *
				      * This table converts a
				      * component number to the
				      * #base_element# number. While
				      * component information contains
				      * multiplicity of base elements,
				      * the result allows access to
				      * shape functions of the base
				      * element.
				      */
    vector<unsigned> component_to_base_table;
    
				     /**
				      * Helper function used in the constructor:
				      * take a #FiniteElementData# object
				      * and return an object of the same type
				      * with the number of degrees of
				      * freedom per vertex, line, etc.
				      * multiplied by #n#. Don't touch the
				      * number of functions for the
				      * transformation from unit to real
				      * cell.
				      */
    static FiniteElementData<dim>
    multiply_dof_numbers (const FiniteElementData<dim> &fe_data,
			  const unsigned int            N);
    
				     /**
				      * Same as above for mixed elements
				      * with two different sub-elements.
				      */
    static FiniteElementData<dim>
    multiply_dof_numbers (const FiniteElementData<dim> &fe1,
			  const unsigned int            N1,
			  const FiniteElementData<dim> &fe2,
			  const unsigned int            N2);

    				     /**
				      * Same as above for mixed elements
				      * with three different sub-elements.
				      */
    static FiniteElementData<dim>
    multiply_dof_numbers (const FiniteElementData<dim> &fe1,
			  const unsigned int            N1,
			  const FiniteElementData<dim> &fe2,
			  const unsigned int            N2,
			  const FiniteElementData<dim> &fe3,
			  const unsigned int            N3);

				     /**
				      * This function is simply singled out of
				      * the constructor. It sets up the
				      * index table for the system as well as
				      * #restriction# and #prolongation#
				      * matrices. Since the operation of this
				      * function can be done without explicit
				      * knowledge of the data type of the
				      * underlying finite element class, we
				      * don't want to have this function in
				      * the general template definition in
				      * the #.h# file.
				      */
    void initialize();

				     /**
				      * Used by #initialize#.
				      */
    void build_cell_table();
    
				     /**
				      * Used by #initialize#.
				      */
    void build_face_table();

				     /**
				      * Used by #initialize#.
				      */
    void build_interface_constraints ();
    
				     /**
				      *Exception.
				      */
    DeclException0(ExcElementTransformNotEqual);
};





/* ------------------------- template functions ------------------------- */

template<int dim>
inline unsigned
FESystem<dim>::n_base_elements() const
{
  return base_elements.size();
}



template<int dim>
inline unsigned
FESystem<dim>::element_multiplicity(unsigned index) const
{
  Assert (index < base_elements.size(), 
	  ExcIndexRange(index, 0, base_elements.size()));
  return base_elements[index].second;
}



template <int dim>
inline const FiniteElement<dim>&
FESystem<dim>::base_element(unsigned index) const
{
  Assert (index < base_elements.size(), 
	  ExcIndexRange(index, 0, base_elements.size()));
  return *base_elements[index].first;
}



template <int dim>
template <class FE>
FESystem<dim>::FESystem (const FE &fe, const unsigned int n_elements) :
		FiniteElement (multiply_dof_numbers(fe, n_elements)),
		base_elements(1),
		component_to_base_table(n_components)
{
  base_elements[0] = ElementPair(new FE, n_elements);
  base_elements[0].first -> subscribe ();
  initialize ();
};



template <int dim>
template <class FE1, class FE2>
FESystem<dim>::FESystem (const FE1 &fe1, const unsigned int n1,
			 const FE2 &fe2, const unsigned int n2)
		:
		FiniteElement (multiply_dof_numbers(fe1, n1, fe2, n2)),
		base_elements(2),
		component_to_base_table(n_components)
{
  Assert(fe1.n_transform_functions == fe2.n_transform_functions,
	 ExcElementTransformNotEqual());
  
  base_elements[0] = ElementPair(new FE1, n1);
  base_elements[0].first -> subscribe ();
  base_elements[1] = ElementPair(new FE2, n2);
  base_elements[1].first -> subscribe ();
  initialize ();
};



template <int dim>
template <class FE1, class FE2, class FE3>
FESystem<dim>::FESystem (const FE1 &fe1, const unsigned int n1,
			 const FE2 &fe2, const unsigned int n2,
			 const FE3 &fe3, const unsigned int n3)
		:
		FiniteElement (multiply_dof_numbers(fe1, n1,
						    fe2, n2,
						    fe3, n3)),
		base_elements(3),
		component_to_base_table(n_components)
{
  Assert(fe1.n_transform_functions == fe2.n_transform_functions,
	 ExcElementTransformNotEqual());
  Assert(fe1.n_transform_functions == fe3.n_transform_functions,
	 ExcElementTransformNotEqual());
  
  base_elements[0] = ElementPair(new FE1, n1);
  base_elements[0].first -> subscribe ();
  base_elements[1] = ElementPair(new FE2, n2);
  base_elements[1].first -> subscribe ();
  base_elements[2] = ElementPair(new FE3, n3);
  base_elements[2].first -> subscribe ();
  initialize ();
};




/*----------------------------  fe_lib.system.h  ---------------------------*/
/* end of #ifndef __fe_system_H */
#endif
/*----------------------------  fe_lib.system.h  ---------------------------*/
