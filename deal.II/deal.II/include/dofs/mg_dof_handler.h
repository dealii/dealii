/*----------------------------   mg_dof.h     ---------------------------*/
/*      $Id$                 */
#ifndef __mg_dof_H
#define __mg_dof_H
/*----------------------------   mg_dof.h     ---------------------------*/


#include <grid/dof.h>



// forward declarations
template <int dim, typename BaseClass> class MGDoFLineAccessor;
template <int dim, typename BaseClass> class MGDoFLineAccessor;
template <int dim, typename BaseClass> class MGDoFQuadAccessor;
template <int dim, typename BaseClass> class MGDoFQuadAccessor;





template <int dim>
class MGDoFHandler : public DoFHandler<dim> {
  public:
    				     /**
				      * Constructor. Take #tria# as the
				      * triangulation to work on.
				      */
    MGDoFHandler (Triangulation<dim> *tria);

				     /**
				      * Destructor
				      */
    virtual ~MGDoFHandler ();
    
				     /**
				      * Go through the triangulation and
				      * distribute the degrees of freedoms
				      * needed for the given finite element
				      * according to the given distribution
				      * method. We first call the #DoFHandler#'s
				      * function and then distribute the
				      * levelwise numbers.
				      *
				      * A copy of the transferred finite
				      * element is stored.
				      *
				      * This function uses the user flags of the
				      * triangulation object, so make sure you
				      * don't need them after calling this
				      * function, or if so store them.
				      */
    virtual void distribute_dofs (const FiniteElementBase<dim> &);

				     /**
				      * Renumber the degrees of freedom according
				      * to the given scheme, eventually using
				      * constraint information and the given
				      * starting points. The starting points
				      * default to an empty list, the use of
				      * constraint information defaults to
				      * false.
				      *
				      * See the general documentation of the
				      * #DoFHandler# class for more details.
				      */
    virtual void renumber_dofs (const RenumberingMethod method,
				const bool use_constraints         = false,
				const vector<int> &starting_points = vector<int>());


  private:

				     /** We need for each vertex a list of the
				      * degree of freedom indices on each of the
				      * levels this vertex lives on. Since most
				      * vertices live only on a few levels, it
				      * is not economical to reserve indices for
				      * all the levels there are; rather, we
				      * create an object which holds the indices
				      * on those levels only where the vertex
				      * lives. To construct such an array, it
				      * is necessary to know beforehand which
				      * is the coarsest level the vertex lives
				      * on, how many levels it lives on and
				      * how many dofs there are on each vertex.
				      * If we have this information, we can
				      * allocate exactly the amount of memory
				      * which is needed and need not handle
				      * growing arrays and the like.
				      */
    class MGVertexDoFs {
      public:
	
					 /**
					  * Constructor. Allocates memory and
					  * sets all indices to #-1#.
					  */
	MGVertexDoFs (const unsigned int coarsest_level,
		      const unsigned int n_levels,
		      const unsigned int dofs_per_vertex);

					 /**
					  * Destructor
					  */
	~MGVertexDoFs ();

					 /**
					  * Set the index with number
					  * #dof_number# of this vertex on
					  * #level# to the given index. To
					  * compute the position in the array,
					  * one has to specify how many
					  * dofs per vertex there are. It is
					  * not checked that the level number
					  * is below the number of the finest
					  * level this vertex lives on.
					  *
					  * The function is inline, so should
					  * be reasonably fast.
					  */
	void set_index (const unsigned int level,
			const unsigned int dof_number,
			const unsigned int dofs_per_vertex,
			const unsigned int index);
	
					 /**
					  * Return the index with number
					  * #dof_number# of this vertex on
					  * #level#. To
					  * compute the position in the array,
					  * one has to specify how many
					  * dofs per vertex there are. It is
					  * not checked that the level number
					  * is below the number of the finest
					  * level this vertex lives on.
					  *
					  * The function is inline, so should
					  * be reasonably fast.
					  */
	int get_index (const unsigned int level,
		       const unsigned int dof_number,
		       const unsigned int dofs_per_vertex) const;
	
					 /**
					  * Exception.
					  */
	DeclException0 (ExcNoMemory);
					 /**
					  * Exception.
					  */
	DeclException0 (ExcInvalidIndex);
					 /**
					  * Exception.
					  */
	DeclException1 (ExcInvalidLevel,
			int,
			<< "The given level number " << arg1 << " is below the "
			<< "coarsest level this vertex lives on.");
	
      private:
					 /**
					  * Store the coarsest level this
					  * vertex lives on. This is used
					  * as an offset when accessing the
					  * dofs of a specific level.
					  */
	const unsigned int coarsest_level;

					 /**
					  * Array holding the indices.
					  */
	int *const         indices;
    };

    
				     /**
				      * Reserve enough space for the MG dof
				      * indices for a given triangulation.
				      */
    void reserve_space ();

    				     /**
				      * Space to store the DoF numbers for the
				      * different levels. Unlike the #levels#
				      * object in the #DoFHandler#, these are
				      * not global numbers but rather are
				      * numbers which start from zero on each
				      * level.
				      */
    vector<DoFLevel<dim>*>    mg_levels;

				     /**
				      * For each vertex there is a list of
				      * indices of the degrees of freedom indices
				      * on the different levels it lives on and
				      * which are these levels.
				      */
    vector<MGVertexDoFs>      mg_vertex_dofs;
    
				     /**
				      * Vectors storing the number of degrees of
				      * freedom on each level.
				      */
    vector<unsigned int>      mg_used_dofs;

    friend class MGDoFLineAccessor<dim, LineAccessor<dim> >;
    friend class MGDoFLineAccessor<dim, CellAccessor<dim> >;
    friend class MGDoFQuadAccessor<dim, QuadAccessor<dim> >;
    friend class MGDoFQuadAccessor<dim, CellAccessor<dim> >;
};

    





/* ----------------------- Inline functions of MGVertexDoFs -------------------*/


template <int dim>
inline
void MGDoFHandler<dim>::MGVertexDoFs::set_index  (const unsigned int level,
						  const unsigned int dof_number,
						  const unsigned int dofs_per_vertex,
						  const unsigned int index) {
  Assert (level >= coarsest_level, ExcInvalidLevel(level));
  Assert (dof_number < dofs_per_vertex, ExcInvalidIndex ());
  
  indices[(level-coarsest_level)*dofs_per_vertex + dof_number] = index;
};




template <int dim>
inline
int MGDoFHandler<dim>::MGVertexDoFs::get_index  (const unsigned int level,
						 const unsigned int dof_number,
						 const unsigned int dofs_per_vertex) const {
  Assert (level >= coarsest_level, ExcInvalidLevel(level));
  Assert (dof_number < dofs_per_vertex, ExcInvalidIndex ());
  
  return indices[(level-coarsest_level)*dofs_per_vertex + dof_number];
};



/*----------------------------   mg_dof.h     ---------------------------*/
/* end of #ifndef __mg_dof_H */
#endif
/*----------------------------   mg_dof.h     ---------------------------*/
