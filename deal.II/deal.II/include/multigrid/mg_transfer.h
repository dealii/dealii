//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------
#ifndef __deal2__mg_transfer_h
#define __deal2__mg_transfer_h


#include <lac/block_vector.h>
#ifdef DEAL_PREFER_MATRIX_EZ
#  include <lac/sparse_matrix_ez.h>
#  include <lac/block_sparse_matrix_ez.h>
#else
#  include <lac/sparsity_pattern.h>
#  include <lac/block_sparsity_pattern.h>
#endif
#include <multigrid/mg_base.h>
#include <multigrid/mg_level_object.h>

#include <boost_local/shared_ptr.hpp>


template <int dim> class MGDoFHandler;

/*
 * MGTransferBase is defined in mg_base.h
 */

/**
 * Implementation of the @p{MGTransferBase} interface for which the transfer
 * operations are prebuilt upon construction of the object of this class as
 * matrices. This is the fast way, since it only needs to build the operation
 * once by looping over all cells and storing the result in a matrix for
 * each level, but requires additional memory.
 *
 * See @ref{MGTransferBase} to find out which of the transfer classes
 * is best for your needs.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1999-2003
 */
template <typename number>
class MGTransferPrebuilt : public MGTransferBase<Vector<number> > 
{
  public:
				     /**
				      * Destructor.
				      */
    virtual ~MGTransferPrebuilt ();
    
				     /**
				      * Actually build the prolongation
				      * matrices for each level.
				      */
    template <int dim>
    void build_matrices (const MGDoFHandler<dim> &mg_dof);

				     /**
				      * Prolongate a vector from level
				      * @p{to_level-1} to level
				      * @p{to_level}. The previous
				      * content of @p{dst} is
				      * overwritten.
				      *
				      * @p{src} is assumed to be a vector with
				      * as many elements as there are degrees
				      * of freedom on the coarser level of
				      * the two involved levels, while @p{src}
				      * shall have as many elements as there
				      * are degrees of freedom on the finer
				      * level.
				      */
    virtual void prolongate (const unsigned int    to_level,
			     Vector<number>       &dst,
			     const Vector<number> &src) const;

				     /**
				      * Restrict a vector from level
				      * @p{from_level} to level
				      * @p{from_level-1} and add this
				      * restriction to
				      * @p{dst}. Obviously, if the
				      * refined region on level
				      * @p{from_level} is smaller than
				      * that on level @p{from_level-1},
				      * some degrees of freedom in
				      * @p{dst} are not covered and will
				      * not be altered. For the other
				      * degress of freedom, the result
				      * of the restriction is added.
				      *
				      * @p{src} is assumed to be a vector with
				      * as many elements as there are degrees
				      * of freedom on the finer level of
				      * the two involved levels, while @p{src}
				      * shall have as many elements as there
				      * are degrees of freedom on the coarser
				      * level.
				      */
    virtual void restrict_and_add (const unsigned int    from_level,
				   Vector<number>       &dst,
				   const Vector<number> &src) const;

    				     /**
				      * Transfer from a vector on the
				      * global grid to vectors defined
				      * on each of the levels
				      * separately, i.a. an @p{MGVector}.
				      */
    template <int dim, class InVector>
    void
    copy_to_mg (const MGDoFHandler<dim>        &mg_dof,
		MGLevelObject<Vector<number> > &dst,
		const InVector                 &src) const;

				     /**
				      * Transfer from multi-level vector to
				      * normal vector.
				      *
				      * Copies data from active
				      * portions of an MGVector into
				      * the respective positions of a
				      * @p{Vector<number>}. In order to
				      * keep the result consistent,
				      * constrained degrees of freedom
				      * are set to zero.
				      */
    template <int dim, class OutVector>
    void
    copy_from_mg (const MGDoFHandler<dim>              &mg_dof,
		  OutVector                            &dst,
		  const MGLevelObject<Vector<number> > &src) const;

				     /**
				      * Add a multi-level vector to a
				      * normal vector.
				      *
				      * Works as the previous
				      * function, but probably not for
				      * continuous elements.
				      */
    template <int dim, class OutVector>
    void
    copy_from_mg_add (const MGDoFHandler<dim>              &mg_dof,
		      OutVector                            &dst,
		      const MGLevelObject<Vector<number> > &src) const;

				     /**
				      * Finite element does not
				      * provide prolongation matrices.
				      */
    DeclException0(ExcNoProlongation);
    
				     /**
				      * Call @p{build_matrices}
				      * function first.
				      */
    DeclException0(ExcMatricesNotBuilt);

    				     /**
				      * Memory used by this object.
				      */
    unsigned int memory_consumption () const;
    

  private:

				   /**
				    * Sizes of the multi-level vectors.
				    */
    std::vector<unsigned int> sizes;

#ifdef DEAL_PREFER_MATRIX_EZ
				     /**
				      * The actual prolongation matrix.
				      * column indices belong to the
				      * dof indices of the mother cell,
				      * i.e. the coarse level.
				      * while row indices belong to the
				      * child cell, i.e. the fine level.
				      */
    std::vector<boost::shared_ptr<SparseMatrixEZ<double> > > prolongation_matrices;
#else
				     /**
				      * Sparsity patterns for transfer
				      * matrices.
				      */
    std::vector<boost::shared_ptr<SparsityPattern> >   prolongation_sparsities;

				     /**
				      * The actual prolongation matrix.
				      * column indices belong to the
				      * dof indices of the mother cell,
				      * i.e. the coarse level.
				      * while row indices belong to the
				      * child cell, i.e. the fine level.
				      */
    std::vector<boost::shared_ptr<SparseMatrix<double> > > prolongation_matrices;
#endif

				     /**
				      * Structure that is used to
				      * disambiguate calls to
				      * @p{copy_to_mg} for 1d and
				      * non-1d. We provide two
				      * functions of @p{copy_to_mg},
				      * where the 1d function takes an
				      * argument of type
				      * @p{is_1d<true>} and the other
				      * one of type @p{is_1d<false>}.
				      */
    template <bool> struct is_1d {};
    
				     /**
				      * Implementation of the
				      * @p{copy_to_mg} function for
				      * 1d. We have to resort to some
				      * template trickery because we
				      * can't specialize template
				      * functions on the (outer)
				      * template of the class, without
				      * also fully specializing the
				      * inner (member function)
				      * template parameters. However,
				      * it can be done by adding the
				      * additional argument that
				      * converts template
				      * specialization into function
				      * overloading.
				      */
    template <int dim, class InVector>
    void
    copy_to_mg (const MGDoFHandler<dim>        &mg_dof,
		MGLevelObject<Vector<number> > &dst,
		const InVector                 &src,
		const is_1d<true>              &) const;

				     /**
				      * Same for all other space
				      * dimensions.
				      */
    template <int dim, class InVector>
    void
    copy_to_mg (const MGDoFHandler<dim>        &mg_dof,
		MGLevelObject<Vector<number> > &dst,
		const InVector                 &src,
		const is_1d<false>             &) const;
};


/**
 * Implementation of matrix generation for @ref{MGTransferBlock} and
 * @p{MGTransferSelect}.
 *
 * @author Guido Kanschat, 2001-2003
 */
class MGTransferBlockBase
{
  public:
    				     /**
				      * Memory used by this object.
				      */
    unsigned int memory_consumption () const;
    

  protected:  
				     /**
				      * Actually build the prolongation
				      * matrices for each level.
				      *
				      * If a field selected is given,
				      * only matrices for these
				      * components are to be built. By
				      * default, all matrices are
				      * built.
				      *
				      * This function is only called
				      * by derived classes. These can
				      * also set the member variable
				      * @p{target_component} for
				      * re-ordering and grouping of
				      * components.
				      */
    template <int dim>
    void build_matrices (const MGDoFHandler<dim>& mg_dof,
			 const std::vector<bool>& selected);

				   /**
				    * Flag of selected components.
				    */
    std::vector<bool> selected;

				     /**
				      * Target component if renumbering is required.
				      */
    std::vector<unsigned int> target_component;
    
				     /**
				      * Target component if
				      * renumbering of level vectors
				      * is required.
				      */
    std::vector<unsigned int> mg_target_component;
    
				   /**
				    * Sizes of the multi-level vectors.
				    */
    std::vector<std::vector<unsigned int> > sizes;
  
				   /**
				    * Start index of each component.
				    */
    std::vector<unsigned int> component_start;
  
				   /**
				    * Start index of each component on
				    * all levels.
				    */
    std::vector<std::vector<unsigned int> > mg_component_start;

				     /**
				      * Call @p{build_matrices}
				      * function first.
				      */
    DeclException0(ExcMatricesNotBuilt);

#ifdef DEAL_PREFER_MATRIX_EZ
  protected:
  
				     /**
				      * The actual prolongation matrix.
				      * column indices belong to the
				      * dof indices of the mother cell,
				      * i.e. the coarse level.
				      * while row indices belong to the
				      * child cell, i.e. the fine level.
				      */
    std::vector<boost::shared_ptr<BlockSparseMatrixEZ<double> > > prolongation_matrices;
#else
  private:
    std::vector<boost::shared_ptr<BlockSparsityPattern> >   prolongation_sparsities;

  protected:
  
				     /**
				      * The actual prolongation matrix.
				      * column indices belong to the
				      * dof indices of the mother cell,
				      * i.e. the coarse level.
				      * while row indices belong to the
				      * child cell, i.e. the fine level.
				      */
    std::vector<boost::shared_ptr<BlockSparseMatrix<double> > > prolongation_matrices;
#endif
};

//TODO:[GK] Update this class

/**
 * Implementation of the @p{MGTransferBase} interface for block
 * matrices and block vectors.
 *
 * Warning! Due to additional requirements on MGTransferSelect, the
 * implementation in the base class has changed. Therefore, this class
 * is left in an untested state. If you use it and you encounter
 * problems, please contact Guido Kanschat.
 *
 * In addition to the functionality of
 * @ref{MGTransferPrebuilt}, the operation may be restricted to
 * certain blocks of the vector.
 *
 * If the restricted mode is chosen, block vectors used in the
 * transfer routines may only have as many components as there are
 * @p{true}s in the selected-field.
 *
 * See @ref{MGTransferBase} to find out which of the transfer classes
 * is best for your needs.
 *
 * @author Guido Kanschat, 2001, 2002
 */
template <typename number>
class MGTransferBlock : public MGTransferBase<BlockVector<number> >,
			private MGTransferBlockBase
{
  public:
				     /**
				      * Destructor.
				      */
    virtual ~MGTransferBlock ();
    
				     /**
				      * Build the prolongation
				      * matrices for each level.
				      *
				      * If a field selected is given,
				      * only matrices for these
				      * components are to be built. By
				      * default, all matrices are
				      * built.
				      *
				      * The last argument
				      * @p{target_component} allows
				      * grouping of components the
				      * same way as in
				      * @p{DoFRenumbering::component_wise}.
				      * 
				      * This function is a front-end
				      * for the same function in
				      * @ref{MGTransferBlockBase}.
				      */
    template <int dim>
    void build_matrices (const MGDoFHandler<dim> &mg_dof,
			 std::vector<bool> selected = std::vector<bool>(),
			 const std::vector<unsigned int>& target_component
			 = std::vector<unsigned int>(),
			 const std::vector<unsigned int>& mg_target_component
			 =std::vector<unsigned int>());

				     /**
				      * Prolongate a vector from level
				      * @p{to_level-1} to level
				      * @p{to_level}. The previous
				      * content of @p{dst} is
				      * overwritten.
				      *
				      * @p{src} is assumed to be a vector with
				      * as many elements as there are degrees
				      * of freedom on the coarser level of
				      * the two involved levels, while @p{src}
				      * shall have as many elements as there
				      * are degrees of freedom on the finer
				      * level.
				      */
    virtual void prolongate (const unsigned int    to_level,
			     BlockVector<number>       &dst,
			     const BlockVector<number> &src) const;

				     /**
				      * Restrict a vector from level
				      * @p{from_level} to level
				      * @p{from_level-1} and add this
				      * restriction to
				      * @p{dst}. Obviously, if the
				      * refined region on level
				      * @p{from_level} is smaller than
				      * that on level @p{from_level-1},
				      * some degrees of freedom in
				      * @p{dst} are not covered and will
				      * not be altered. For the other
				      * degress of freedom, the result
				      * of the restriction is added.
				      *
				      * @p{src} is assumed to be a vector with
				      * as many elements as there are degrees
				      * of freedom on the finer level of
				      * the two involved levels, while @p{src}
				      * shall have as many elements as there
				      * are degrees of freedom on the coarser
				      * level.
				      */
    virtual void restrict_and_add (const unsigned int    from_level,
				   BlockVector<number>       &dst,
				   const BlockVector<number> &src) const;

    				     /**
				      * Transfer from a vector on the
				      * global grid to vectors defined
				      * on each of the levels
				      * separately, i.a. an @p{MGVector}.
				      */
    template <int dim, class InVector>
    void
    copy_to_mg (const MGDoFHandler<dim>             &mg_dof,
		MGLevelObject<BlockVector<number> > &dst,
		const InVector                      &src) const;

				     /**
				      * Transfer from multi-level vector to
				      * normal vector.
				      *
				      * Copies data from active
				      * portions of an MGVector into
				      * the respective positions of a
				      * @p{Vector<number>}. In order to
				      * keep the result consistent,
				      * constrained degrees of freedom
				      * are set to zero.
				      */
    template <int dim, class OutVector>
    void
    copy_from_mg (const MGDoFHandler<dim>                   &mg_dof,
		  OutVector                                 &dst,
		  const MGLevelObject<BlockVector<number> > &src) const;

				     /**
				      * Add a multi-level vector to a
				      * normal vector.
				      *
				      * Works as the previous
				      * function, but probably not for
				      * continuous elements.
				      */
    template <int dim, class OutVector>
    void
    copy_from_mg_add (const MGDoFHandler<dim>                   &mg_dof,
		      OutVector                                 &dst,
		      const MGLevelObject<BlockVector<number> > &src) const;

    MGTransferBlockBase::memory_consumption;
    
  private:
				     /**
				      * Structure that is used to
				      * disambiguate calls to
				      * @p{copy_to_mg} for 1d and
				      * non-1d. We provide two
				      * functions of @p{copy_to_mg},
				      * where the 1d function takes an
				      * argument of type
				      * @p{is_1d<true>} and the other
				      * one of type @p{is_1d<false>}.
				      */
    template <bool> struct is_1d {};

				     /**
				      * Implementation of the
				      * @p{copy_to_mg} function for
				      * 1d. We have to resort to some
				      * template trickery because we
				      * can't specialize template
				      * functions on the (outer)
				      * template of the class, without
				      * also fully specializing the
				      * inner (member function)
				      * template parameters. However,
				      * it can be done by adding the
				      * additional argument that
				      * converts template
				      * specialization into function
				      * overloading.
				      */
    template <int dim, class InVector>
    void
    copy_to_mg (const MGDoFHandler<dim>             &mg_dof,
		MGLevelObject<BlockVector<number> > &dst,
		const InVector                      &src,
		const is_1d<true>                   &) const;

				     /**
				      * Same for all other space
				      * dimensions.
				      */
    template <int dim, class InVector>
    void
    copy_to_mg (const MGDoFHandler<dim>             &mg_dof,
		MGLevelObject<BlockVector<number> > &dst,
		const InVector                      &src,
		const is_1d<false>                  &) const;    
};


//TODO:[GK] Update documentation for copy_* functions

/**
 * Implementation of the @p{MGTransferBase} interface for block
 * matrices and simple vectors. This class uses @ref{MGTransferBlockBase}
 * selecting a single component or grouping several components into a
 * single block. The transfer operators themselves are implemented for
 * Vector and BlockVector objects.
 *
 * See @ref{MGTransferBase} to find out which of the transfer classes
 * is best for your needs.
 *
 * @author Guido Kanschat, 2001, 2002, 2003
 */
template <typename number>
class MGTransferSelect : public MGTransferBase<Vector<number> >,
			 private MGTransferBlockBase
{
  public:
				     /**
				      * Destructor.
				      */
    virtual ~MGTransferSelect ();
    
				     /**
				      * Actually build the prolongation
				      * matrices for grouped components.
				      *
				      * @p{selected} tells the copy
				      * functions operating on single
				      * vectors, which component this
				      * vector belongs to.
				      *
				      * @p{mg_selected} is the number
				      * of the component for which the
				      * transfer matrices should be
				      * built.
				      *
				      * The argument
				      * @p{target_component}
				      * corresponds to the grouping
				      * mechanism of
				      * @p{DoFRenumbering::component_wise(...)},
				      * which should be used to create
				      * the corresponding block
				      * structure in matrices and
				      * vectors.
				      *
				      * This function is a front-end
				      * for the same function in
				      * @ref{MGTransferBlockBase}.
				      */
    template <int dim>
    void build_matrices (const MGDoFHandler<dim> &mg_dof,
			 unsigned int selected,
			 unsigned int mg_selected,
			 const std::vector<unsigned int>& target_component
			 = std::vector<unsigned int>(),
			 const std::vector<unsigned int>& mg_target_component
			 = std::vector<unsigned int>());

				     /**
				      * Change selected
				      * component. Handle with care!
				      */
    void select (const unsigned int component);
    
				     /**
				      * Prolongate a vector from level
				      * @p{to_level-1} to level
				      * @p{to_level}. The previous
				      * content of @p{dst} is
				      * overwritten.
				      *
				      * @p{src} is assumed to be a vector with
				      * as many elements as there are degrees
				      * of freedom on the coarser level of
				      * the two involved levels, while @p{src}
				      * shall have as many elements as there
				      * are degrees of freedom on the finer
				      * level.
				      */
    virtual void prolongate (const unsigned int    to_level,
			     Vector<number>       &dst,
			     const Vector<number> &src) const;

				     /**
				      * Restrict a vector from level
				      * @p{from_level} to level
				      * @p{from_level-1} and add this
				      * restriction to
				      * @p{dst}. Obviously, if the
				      * refined region on level
				      * @p{from_level} is smaller than
				      * that on level @p{from_level-1},
				      * some degrees of freedom in
				      * @p{dst} are not covered and will
				      * not be altered. For the other
				      * degress of freedom, the result
				      * of the restriction is added.
				      *
				      * @p{src} is assumed to be a vector with
				      * as many elements as there are degrees
				      * of freedom on the finer level of
				      * the two involved levels, while @p{src}
				      * shall have as many elements as there
				      * are degrees of freedom on the coarser
				      * level.
				      */
    virtual void restrict_and_add (const unsigned int    from_level,
				   Vector<number>       &dst,
				   const Vector<number> &src) const;

				     /**
				      * Structure that is used to
				      * disambiguate calls to
				      * @p{copy_to_mg} for 1d and
				      * non-1d. We provide two
				      * functions of @p{copy_to_mg},
				      * where the 1d function takes an
				      * argument of type
				      * @p{is_1d<true>} and the other
				      * one of type @p{is_1d<false>}.
				      */
    template <bool> struct is_1d {};
    
    				     /**
				      * Transfer from a vector on the
				      * global grid to vectors defined
				      * on each of the levels
				      * separately, i.a. an @p{MGVector}.
				      */
    template <int dim, typename number2>
    void
    copy_to_mg (const MGDoFHandler<dim>        &mg_dof,
		MGLevelObject<Vector<number> > &dst,
		const Vector<number2>          &src) const;

				     /**
				      * Transfer from multi-level vector to
				      * normal vector.
				      *
				      * Copies data from active
				      * portions of an MGVector into
				      * the respective positions of a
				      * @p{Vector<number>}. In order to
				      * keep the result consistent,
				      * constrained degrees of freedom
				      * are set to zero.
				      */
    template <int dim, typename number2>
    void
    copy_from_mg (const MGDoFHandler<dim>              &mg_dof,
		  Vector<number2>                      &dst,
		  const MGLevelObject<Vector<number> > &src) const;

				     /**
				      * Add a multi-level vector to a
				      * normal vector.
				      *
				      * Works as the previous
				      * function, but probably not for
				      * continuous elements.
				      */
    template <int dim, typename number2>
    void
    copy_from_mg_add (const MGDoFHandler<dim>              &mg_dof,
		      Vector<number2>                      &dst,
		      const MGLevelObject<Vector<number> > &src) const;

    				     /**
				      * Transfer from a vector on the
				      * global grid to vectors defined
				      * on each of the levels
				      * separately, i.a. an @p{MGVector}.
				      */
    template <int dim, typename number2>
    void
    copy_to_mg (const MGDoFHandler<dim>        &mg_dof,
		MGLevelObject<Vector<number> > &dst,
		const BlockVector<number2>     &src) const;

				     /**
				      * Transfer from multi-level vector to
				      * normal vector.
				      *
				      * Copies data from active
				      * portions of an MGVector into
				      * the respective positions of a
				      * @p{Vector<number>}. In order to
				      * keep the result consistent,
				      * constrained degrees of freedom
				      * are set to zero.
				      */
    template <int dim, typename number2>
    void
    copy_from_mg (const MGDoFHandler<dim>              &mg_dof,
		  BlockVector<number2>                 &dst,
		  const MGLevelObject<Vector<number> > &src) const;

				     /**
				      * Add a multi-level vector to a
				      * normal vector.
				      *
				      * Works as the previous
				      * function, but probably not for
				      * continuous elements.
				      */
    template <int dim, typename number2>
    void
    copy_from_mg_add (const MGDoFHandler<dim>              &mg_dof,
		      BlockVector<number2>                 &dst,
		      const MGLevelObject<Vector<number> > &src) const;

				     /**
				      * Memory used by this object.
				      */
    unsigned int memory_consumption () const;
    
  private:
				     /**
				      * Transfer from multi-level vector to
				      * normal vector.
				      *
				      * Copies data from active
				      * portions of an MGVector into
				      * the respective positions of a
				      * @p{Vector<number>}. In order to
				      * keep the result consistent,
				      * constrained degrees of freedom
				      * are set to zero.
				      */
    template <int dim, class OutVector>
    void
    do_copy_from_mg (const MGDoFHandler<dim>              &mg_dof,
		     OutVector                            &dst,
		     const MGLevelObject<Vector<number> > &src,
		     const unsigned int offset) const;

				     /**
				      * Add a multi-level vector to a
				      * normal vector.
				      *
				      * Works as the previous
				      * function, but probably not for
				      * continuous elements.
				      */
    template <int dim, class OutVector>
    void
    do_copy_from_mg_add (const MGDoFHandler<dim>              &mg_dof,
			 OutVector                            &dst,
			 const MGLevelObject<Vector<number> > &src,
			 const unsigned int offset) const;

				     /**
				      * Implementation of the
				      * @p{copy_to_mg} function for
				      * 1d. We have to resort to some
				      * template trickery because we
				      * can't specialize template
				      * functions on the (outer)
				      * template of the class, without
				      * also fully specializing the
				      * inner (member function)
				      * template parameters. However,
				      * it can be done by adding the
				      * additional argument that
				      * converts template
				      * specialization into function
				      * overloading.
				      */
    template <int dim, class InVector>
    void
    do_copy_to_mg (const MGDoFHandler<dim>        &mg_dof,
		   MGLevelObject<Vector<number> > &dst,
		   const InVector                 &src,
		   const unsigned int              offset,
		   const is_1d<true>              &) const;

				     /**
				      * Same for all other space
				      * dimensions.
				      */
    template <int dim, class InVector>
    void
    do_copy_to_mg (const MGDoFHandler<dim>        &mg_dof,
		   MGLevelObject<Vector<number> > &dst,
		   const InVector                 &src,
		   const unsigned int              offset,
		   const is_1d<false>             &) const;

                                     /**
                                      * Selected component.
                                      */
    unsigned int selected_component;
};

//----------------------------------------------------------------------//
template <typename number>
inline void
MGTransferSelect<number>::select(const unsigned int component)
{
  selected_component = component;
}


#endif
