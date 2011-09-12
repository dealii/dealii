//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2008, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__distributed_tria_h
#define __deal2__distributed_tria_h


#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/grid/tria.h>

#include <deal.II/base/std_cxx1x/function.h>

#include <vector>
#include <list>
#include <utility>

#ifdef DEAL_II_USE_P4EST
#include <p4est_connectivity.h>
#include <p4est.h>
#include <p4est_ghost.h>

#include <p8est_connectivity.h>
#include <p8est.h>
#include <p8est_ghost.h>
#endif


DEAL_II_NAMESPACE_OPEN

template <int, int> class Triangulation;

#ifdef DEAL_II_USE_P4EST

namespace internal
{
  namespace DoFHandler
  {
    namespace Policy
    {
      template <int, int> class ParallelDistributed;
    }
  }
}


namespace internal
{
  namespace p4est
  {
				     /**
				      * A structure whose explicit
				      * specializations contain
				      * typedefs to the relevant
				      * p4est_* and p8est_*
				      * types. Using this
				      * structure, for example by
				      * saying
				      * <code>types@<dim@>::connectivity</code>
				      * we can write code in a
				      * dimension independent way,
				      * either referring to
				      * p4est_connectivity_t or
				      * p8est_connectivity_t,
				      * depending on template
				      * argument.
				      */
    template <int> struct types;

    template <>
    struct types<2>
    {
	typedef p4est_connectivity_t connectivity;
	typedef p4est_t              forest;
	typedef p4est_tree_t         tree;
	typedef p4est_quadrant_t     quadrant;
	typedef p4est_topidx_t       topidx;
	typedef p4est_locidx_t       locidx;
	typedef p4est_balance_type_t balance_type;
	typedef p4est_ghost_t        ghost;
    };

    template <>
    struct types<3>
    {
	typedef p8est_connectivity_t connectivity;
	typedef p8est_t              forest;
	typedef p8est_tree_t         tree;
	typedef p8est_quadrant_t     quadrant;
	typedef p4est_topidx_t       topidx;
	typedef p4est_locidx_t       locidx;
	typedef p8est_balance_type_t balance_type;
	typedef p8est_ghost_t        ghost;
    };


				     /**
				      * Initialize the
				      * GeometryInfo<dim>::max_children_per_cell
				      * children of the cell
				      * p4est_cell.
				      */
    template <int dim>
    void
    init_quadrant_children
    (const typename types<dim>::quadrant & p4est_cell,
     typename types<dim>::quadrant (&p4est_children)[GeometryInfo<dim>::max_children_per_cell]);


				     /**
				      * Initialize quadrant to represent a coarse cell.
				      */
    template <int dim>
    void
    init_coarse_quadrant(typename types<dim>::quadrant & quad);



				     /**
				      * Returns whether q1 and q2 are equal
				      */
    template <int dim>
    bool
    quadrant_is_equal (const typename types<dim>::quadrant & q1,
		       const typename types<dim>::quadrant & q2);

				     //TODO: remove these functions from
				     //public interface somehow? [TH]

				     /**
				      * returns whether q1 is an ancestor of q2
				      */
    template <int dim>
    bool
    quadrant_is_ancestor (const typename types<dim>::quadrant & q1,
			  const typename types<dim>::quadrant & q2);
  }
}



namespace parallel
{
  namespace distributed
  {


/**
 * This class acts like the dealii::Triangulation class, but it
 * distributes the mesh across a number of different processors when
 * using MPI. The class's interface does not add a lot to the
 * dealii::Triangulation class but there are a number of difficult
 * algorithms under the hood that ensure we always have a
 * load-balanced, fully distributed mesh. Use of this class is
 * explained in step-40, step-32, the @ref distributed documentation
 * module, as well as the @ref distributed_paper . See there for more
 * information.
 *
 * @note This class does not support anisotropic refinement, because
 * it relies on the p4est library that does not support this. Attempts
 * to refine cells anisotropically will result in errors.
 *
 * @author Wolfgang Bangerth, Timo Heister 2008, 2009, 2010
 * @ingroup distributed
 */
    template <int dim, int spacedim = dim>
    class Triangulation : public dealii::Triangulation<dim,spacedim>
    {
      public:
					 /**
					  * Import the various
					  * iterator typedefs from the
					  * base class.
					  */
	typedef typename dealii::Triangulation<dim,spacedim>::active_cell_iterator active_cell_iterator;
	typedef typename dealii::Triangulation<dim,spacedim>::cell_iterator        cell_iterator;
	typedef typename dealii::Triangulation<dim,spacedim>::raw_cell_iterator    raw_cell_iterator;

					 /**
					  * Constructor.
					  *
					  * @param mpi_communicator denotes
					  * the MPI communicator to be used for
					  * the triangulation.
					  *
					  * @param smooth_grid Degree
					  * and kind of mesh smoothing
					  * to be applied to the
					  * mesh. See the
					  * dealii::Triangulation
					  * class for a description of
					  * the kinds of smoothing
					  * operations that can be
					  * applied.
					  *
					  * @note This class does not
					  * currently support the
					  * <code>check_for_distorted_cells</code>
					  * argument provided by the
					  * base class.
					  *
					  * @note While it is possible to pass
					  * all of the mesh smoothing flags
					  * listed in the base class to
					  * objects of this type, it is not
					  * always possible to honor all of
					  * these smoothing options if they
					  * would require knowledge of
					  * refinement/coarsening flags on
					  * cells not locally owned by this
					  * processor. As a consequence, for
					  * some of these flags, the ultimate
					  * number of cells of the parallel
					  * triangulation may depend on the
					  * number of processors into which it
					  * is partitioned. On the other hand,
					  * if no smoothing flags are passed,
					  * if you always mark the same cells
					  * of the mesh, you will always get
					  * the exact same refined mesh
					  * independent of the number of
					  * processors into which the
					  * triangulation is partitioned.
					  */
	Triangulation (MPI_Comm mpi_communicator,
		       const typename dealii::Triangulation<dim,spacedim>::MeshSmoothing
		       smooth_grid = (dealii::Triangulation<dim,spacedim>::none));

					 /**
					  * Destructor.
					  */
	virtual ~Triangulation ();

					 /**
					  * Reset this triangulation into a
					  * virgin state by deleting all data.
					  *
					  * Note that this operation is only
					  * allowed if no subscriptions to this
					  * object exist any more, such as
					  * DoFHandler objects using it.
					  */
	virtual void clear ();

					 /**
					  * Implementation of the same
					  * function as in the base
					  * class.
					  */
	virtual void copy_triangulation (const dealii::Triangulation<dim, spacedim> &old_tria);

					 /**
					  * Create a triangulation as
					  * documented in the base
					  * class.
					  *
					  * This function also sets up
					  * the various data
					  * structures necessary to
					  * distribute a mesh across a
					  * number of processors. This
					  * will be necessary once the
					  * mesh is being refined,
					  * though we will always keep
					  * the entire coarse mesh
					  * that is generated by this
					  * function on all
					  * processors.
					  */
	virtual void create_triangulation (const std::vector<Point<spacedim> >    &vertices,
					   const std::vector<CellData<dim> > &cells,
					   const SubCellData                 &subcelldata);

					 /**
					  * Coarsen and refine the
					  * mesh according to
					  * refinement and coarsening
					  * flags set.
					  *
					  * Since the current
					  * processor only has control
					  * over those cells it owns
					  * (i.e. the ones for which
					  * <code>cell-@>subdomain_id()
					  * ==
					  * this-@>locally_owned_subdomain()</code>),
					  * refinement and coarsening
					  * flags are only respected
					  * for those locally owned
					  * cells. Flags may be set on
					  * other cells as well (and
					  * may often, in fact, if you
					  * call
					  * Triangulation::prepare_coarsening_and_refinement)
					  * but will be largely
					  * ignored: the decision to
					  * refine the global mesh
					  * will only be affected by
					  * flags set on locally owned
					  * cells.
					  */
	virtual void execute_coarsening_and_refinement ();

					 /**
					  * Return the subdomain id of
					  * those cells that are owned
					  * by the current
					  * processor. All cells in
					  * the triangulation that do
					  * not have this subdomain id
					  * are either owned by
					  * another processor or have
					  * children that only exist
					  * on other processors.
					  */
	types::subdomain_id_t locally_owned_subdomain () const;

					 /**
					  * Return the number of
					  * active cells in the
					  * triangulation that are
					  * locally owned, i.e. that
					  * have a subdomain_id equal
					  * to
					  * locally_owned_subdomain(). Note
					  * that there may be more
					  * active cells in the
					  * triangulation stored on
					  * the present processor,
					  * such as for example ghost
					  * cells, or cells further
					  * away from the locally
					  * owned block of cells but
					  * that are needed to ensure
					  * that the triangulation
					  * that stores this
					  * processor's set of active
					  * cells still remains
					  * balanced with respect to
					  * the 2:1 size ratio of
					  * adjacent cells.
					  *
					  * As a consequence of the remark
					  * above, the result of this function
					  * is always smaller or equal to the
					  * result of the function with the
					  * same name in the ::Triangulation
					  * base class, which includes the
					  * active ghost and artificial cells
					  * (see also @ref GlossArtificialCell
					  * and @ref GlossGhostCell).
					  */
	unsigned int n_locally_owned_active_cells () const;

					 /**
					  * Return the sum over all
					  * processors of the number
					  * of active cells owned by
					  * each processor. This
					  * equals the overall number
					  * of active cells in the
					  * distributed triangulation.
					  */
	unsigned int n_global_active_cells () const;

					 /**
					  * Return the number of
					  * active cells owned by each
					  * of the MPI processes that
					  * contribute to this
					  * triangulation. The element
					  * of this vector indexed by
					  * locally_owned_subdomain()
					  * equals the result of
					  * n_locally_owned_active_cells().
					  */
	const std::vector<unsigned int> &
	n_locally_owned_active_cells_per_processor () const;

					 /**
					  * Return the MPI
					  * communicator used by this
					  * triangulation.
					  */
	MPI_Comm get_communicator () const;

					 /**
					  * Return the local memory
					  * consumption in bytes.
					  */
	virtual std::size_t memory_consumption () const;

					 /**
					  * Return the local memory
					  * consumption contained in the p4est
					  * data structures alone. This is
					  * already contained in
					  * memory_consumption() but made
					  * available separately for debugging
					  * purposes.
					  */
	virtual std::size_t memory_consumption_p4est () const;

					 /**
					  * A collective operation that produces
					  * a sequence of output files with the
					  * given file base name that contain
					  * the mesh in VTK format.
					  *
					  * More than anything else, this
					  * function is useful for debugging the
					  * interface between deal.II and p4est.
					  */
	void write_mesh_vtk (const char *file_basename) const;

					 /**
					  * Produce a check sum of the
					  * triangulation.  This is a
					  * collective operation and
					  * is mostly useful for
					  * debugging purposes.
					  */
	unsigned int get_checksum () const;

					 /**
					  * Used to inform in the callbacks of
					  * register_data_attach() and
					  * notify_ready_to_unpack() how the
					  * cell with the given cell_iterator
					  * is going to change.  Note that
					  * this may me different then the
					  * refine_flag() and coarsen_flag()
					  * in the cell_iterator because of
					  * refinement constraints that this
					  * machine does not see.
					  */
	enum CellStatus
	{
	      CELL_PERSIST, CELL_REFINE, CELL_COARSEN, CELL_INVALID
	};

					 /**
					  * Register a function with
					  * the current Triangulation
					  * object that will be used
					  * to attach data to active
					  * cells before
					  * execute_coarsening_and_refinement(). In
					  * execute_coarsening_and_refinement()
					  * the Triangulation will
					  * call the given function
					  * pointer and provide
					  * @p size bytes to store
					  * data. If necessary, this data will be
					  * transferred to the new
					  * owner of that cell during repartitioning
					  * the tree. See
					  * notify_ready_to_unpack()
					  * on how to retrieve the
					  * data.
					  *
					  * Callers need to store the
					  * return value.  It
					  * specifies an offset of the
					  * position at which data can
					  * later be retrieved during
					  * a call to
					  * notify_ready_to_unpack().
					  */
	unsigned int
	register_data_attach (const std::size_t size,
			      const std_cxx1x::function<void (const cell_iterator &,
							      const CellStatus,
							      void*)> & pack_callback);

					 /**
					  * The given function is called for
					  * each new active cell and supplies
					  * a pointer to the data saved with
					  * register_data_attach().
					  */
	void
	notify_ready_to_unpack (const unsigned int offset,
				const std_cxx1x::function<void (const cell_iterator &,
								const CellStatus,
								const void*)> & unpack_callback);

					 /**
					  * Returns a permutation vector for the order the coarse
					  * cells are handed of to p4est. For example the first
					  * element i in this vector denotes that the first cell
					  * in hierarchical ordering is the ith deal cell starting
					  * from begin(0).
					  */
	const std::vector<unsigned int> &
	get_p4est_tree_to_coarse_cell_permutation() const;

  private:
					 /**
					  * MPI communicator to be
					  * used for the
					  * triangulation. We create a
					  * unique communicator for
					  * this class, which is a
					  * duplicate of the one
					  * passed to the constructor.
					  */
	MPI_Comm mpi_communicator;

					 /**
					  * The subdomain id to be
					  * used for the current
					  * processor.
					  */
	types::subdomain_id_t my_subdomain;

					 /**
					  * A flag that indicates whether the
					  * triangulation has actual content.
					  */
	bool triangulation_has_content;

					 /**
					  * A structure that contains
					  * some numbers about the
					  * distributed triangulation.
					  */
	struct NumberCache
	{
	    std::vector<unsigned int> n_locally_owned_active_cells;
	    unsigned int              n_global_active_cells;
	};

	NumberCache number_cache;

					 /**
					  * A data structure that holds the
					  * connectivity between trees. Since
					  * each tree is rooted in a coarse grid
					  * cell, this data structure holds the
					  * connectivity between the cells of
					  * the coarse grid.
					  */
	typename dealii::internal::p4est::types<dim>::connectivity *connectivity;

					 /**
					  * A data structure that holds the
					  * local part of the global
					  * triangulation.
					  */
	typename dealii::internal::p4est::types<dim>::forest *parallel_forest;

					 /**
					  * A flag that indicates
					  * whether refinement of a
					  * triangulation is currently
					  * in progress. This flag is
					  * used to disambiguate whether
					  * a call to execute_coarsening_and_triangulation
					  * came from the outside or
					  * through a recursive call. While the
					  * first time we want to take
					  * over work to copy things
					  * from a refined p4est, the
					  * other times we don't want to
					  * get in the way as these
					  * latter calls to
					  * Triangulation::execute_coarsening_and_refinement()
					  * are simply there in order to
					  * re-create a triangulation
					  * that matches the p4est.
					  */
	bool refinement_in_progress;


					 /**
					  *
					  */
	unsigned int attached_data_size;
	unsigned int n_attached_datas;
	typedef  std_cxx1x::function<
	  void(typename Triangulation<dim,spacedim>::cell_iterator, CellStatus, void*)
	  > pack_callback_t;

	typedef std::pair<unsigned int, pack_callback_t> callback_pair_t;

	typedef std::list<callback_pair_t> callback_list_t;
	callback_list_t attached_data_pack_callbacks;





					 /**
					  * Two arrays that store which p4est
					  * tree corresponds to which coarse
					  * grid cell and vice versa. We need
					  * these arrays because p4est goes with
					  * the original order of coarse cells
					  * when it sets up its forest, and then
					  * applies the Morton ordering within
					  * each tree. But if coarse grid cells
					  * are badly ordered this may mean that
					  * individual parts of the forest
					  * stored on a local machine may be
					  * split across coarse grid cells that
					  * are not geometrically
					  * close. Consequently, we apply a
					  * Cuthill-McKee preordering to ensure
					  * that the part of the forest stored
					  * by p4est is located on geometrically
					  * close coarse grid cells.
					  */
	std::vector<unsigned int> coarse_cell_to_p4est_tree_permutation;
	std::vector<unsigned int> p4est_tree_to_coarse_cell_permutation;

					 /**
					  * Return a pointer to the p4est
					  * tree that belongs to the given
					  * dealii_coarse_cell_index()
					  */
	typename dealii::internal::p4est::types<dim>::tree *
	init_tree(const int dealii_coarse_cell_index) const;

					 /**
					  * The function that computes the
					  * permutation between the two data
					  * storage schemes.
					  */
	void setup_coarse_cell_to_p4est_tree_permutation ();

					 /**
					  * Take the contents of a newly created
					  * triangulation we are attached to and
					  * copy it to p4est data structures.
					  *
					  * This function exists in 2d
					  * and 3d variants.
					  */
	void copy_new_triangulation_to_p4est (dealii::internal::int2type<2>);
	void copy_new_triangulation_to_p4est (dealii::internal::int2type<3>);

					 /**
					  * Copy the local part of the refined
					  * forest from p4est into the attached
					  * triangulation.
					  */
	void copy_local_forest_to_triangulation ();


					 /**
					  * Update the number_cache
					  * variable after mesh
					  * creation or refinement.
					  */
	void update_number_cache ();

					 /**
					  * Internal function notifying all
					  * registered classes to attach their
					  * data before repartitioning
					  * occurs. Called from
					  * execute_coarsening_and_refinement().
					  */
	void attach_mesh_data();


	template <int, int> friend class dealii::internal::DoFHandler::Policy::ParallelDistributed;
    };


				     /**
				      * Specialization of the general template
				      * for the 1d case. There is currently no
				      * support for distributing 1d
				      * triangulations. Consequently, all this
				      * class does is throw an exception.
				      */
    template <int spacedim>
    class Triangulation<1,spacedim> : public dealii::Triangulation<1,spacedim>
    {
      public:
					 /**
					  * Constructor. The argument denotes
					  * the MPI communicator to be used for
					  * the triangulation.
					  */
	Triangulation (MPI_Comm mpi_communicator);

					 /**
					  * Destructor.
					  */
	virtual ~Triangulation ();

					 /**
					  * Return the MPI
					  * communicator used by this
					  * triangulation.
					  */
	MPI_Comm get_communicator () const;

					 /**
					  * Returns a permutation vector for the order the coarse
					  * cells are handed of to p4est. For example the first
					  * element i in this vector denotes that the first cell
					  * in hierarchical ordering is the ith deal cell starting
					  * from begin(0).
					  */
	const std::vector<unsigned int> &
	get_p4est_tree_to_coarse_cell_permutation() const;

					 /**
					  * Return the subdomain id of
					  * those cells that are owned
					  * by the current
					  * processor. All cells in
					  * the triangulation that do
					  * not have this subdomain id
					  * are either owned by
					  * another processor or have
					  * children that only exist
					  * on other processors.
					  */
	types::subdomain_id_t locally_owned_subdomain () const;

					 /**
					  * Dummy arrays. This class
					  * isn't usable but the
					  * compiler wants to see
					  * these variables at a
					  * couple places anyway.
					  */
	std::vector<unsigned int> coarse_cell_to_p4est_tree_permutation;
	std::vector<unsigned int> p4est_tree_to_coarse_cell_permutation;
    };
  }
}


#else // DEAL_II_USE_P4EST

namespace parallel
{
  namespace distributed
  {
				     /**
				      * Dummy class the compiler chooses for
				      * parallel distributed triangulations if
				      * we didn't actually configure deal.II
				      * with the p4est library. The existence
				      * of this class allows us to refer to
				      * parallel::distributed::Triangulation
				      * objects throughout the library even if
				      * it is disabled.
				      *
				      * Since the constructor of this class is
				      * private, no such objects can actually
				      * be created if we don't have p4est
				      * available.
				      */
    template <int dim, int spacedim = dim>
    class Triangulation : public dealii::Triangulation<dim,spacedim>
    {
      private:
					 /**
					  * Constructor.
					  */
	Triangulation ();

      public:

					 /**
					  * Destructor.
					  */
	virtual ~Triangulation ();

					 /**
					  * Return the subdomain id of
					  * those cells that are owned
					  * by the current
					  * processor. All cells in
					  * the triangulation that do
					  * not have this subdomain id
					  * are either owned by
					  * another processor or have
					  * children that only exist
					  * on other processors.
					  */
	types::subdomain_id_t locally_owned_subdomain () const;
    };
  }
}


#endif


DEAL_II_NAMESPACE_CLOSE

#endif
