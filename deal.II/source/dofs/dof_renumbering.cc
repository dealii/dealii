//---------------------------------------------------------------------------
//    $Id$
//    Version: $name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/fe/fe.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>

#include <deal.II/multigrid/mg_dof_handler.h>
#include <deal.II/multigrid/mg_dof_accessor.h>
#include <deal.II/multigrid/mg_tools.h>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/king_ordering.hpp>
#include <boost/graph/minimum_degree_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <functional>


// for whatever reason, the random_shuffle function used below needs
// lrand48 to be declared when using -ansi as compiler flag (rather
// than do so itself). however, inclusion of <cstdlib> or <stdlib.h>
// does not help, so we declare that function ourselves. Since this
// holds only for some compiler versions, do so conditionally on a
// ./configure-time test
#ifdef DEAL_II_DECLARE_LRAND48
extern "C" long int lrand48 (void);
#endif

DEAL_II_NAMESPACE_OPEN


namespace DoFRenumbering
{
  namespace internal
  {
// The following two classes are defined to be used in the compute_*
// functions. Using them allows to use the same function to compute
// level numbering for multigrid as well as numbering for the global
// system matrix.
    template <class T>
    class WrapDoFIterator : private T
    {
      public:
	typedef typename T::AccessorType AccessorType;

	WrapDoFIterator (const T& t) : T(t) {}

	void get_dof_indices (std::vector<unsigned int>& v) const
	  {
	    (*this)->get_dof_indices(v);
	  }

	template <class T2>
	bool operator != (const T2& i) const
	  {
	    return (! (T::operator==(i)));
	  }
					 // Allow access to these private operators of T
	using T::operator->;
	using T::operator++;
	using T::operator==;
    };



    template <class T>
    class WrapMGDoFIterator : private T
    {
      public:
	typedef typename T::AccessorType AccessorType;

	WrapMGDoFIterator (const T& t) : T(t) {}

	void get_dof_indices (std::vector<unsigned int>& v) const
	  {
	    (*this)->get_mg_dof_indices(v);
	  }

	bool operator != (const WrapMGDoFIterator<T>& i) const
	  {
	    return (! (T::operator==(i)));
	  }
					 // Allow access to these
					 // private operators of T
	using T::operator->;
	using T::operator++;
	using T::operator==;
    };
  }

  namespace boost
  {
#ifndef DEAL_II_BOOST_GRAPH_COMPILER_BUG
    namespace types
    {
      using namespace ::boost;
      using namespace std;

      typedef adjacency_list<vecS, vecS, undirectedS,
			     property<vertex_color_t, default_color_type,
				      property<vertex_degree_t,int> > > Graph;
      typedef graph_traits<Graph>::vertex_descriptor Vertex;
      typedef graph_traits<Graph>::vertices_size_type size_type;

      typedef std::pair<size_type, size_type> Pair;
    }


    namespace internal
    {
      template <class DH>
      void create_graph (const DH                                                       &dof_handler,
			 const bool                                                      use_constraints,
			 types::Graph                                                   &graph,
			 types::property_map<types::Graph,types::vertex_degree_t>::type &graph_degree)
      {
	{
				   // create intermediate sparsity pattern
				   // (faster than directly submitting
				   // indices)
	  ConstraintMatrix constraints;
	  if (use_constraints)
	    DoFTools::make_hanging_node_constraints (dof_handler, constraints);
	  constraints.close ();
	  CompressedSimpleSparsityPattern csp (dof_handler.n_dofs(),
					       dof_handler.n_dofs());
	  DoFTools::make_sparsity_pattern (dof_handler, csp, constraints);

				   // submit the entries to the boost graph
	  for (unsigned int row=0;row<csp.n_rows(); ++row)
	    for (unsigned int col=0; col < csp.row_length(row); ++col)
	      add_edge (row, csp.column_number (row, col), graph);
	}

	types::graph_traits<types::Graph>::vertex_iterator ui, ui_end;

	graph_degree = get(::boost::vertex_degree, graph);
	for (::boost::tie(ui, ui_end) = vertices(graph); ui != ui_end; ++ui)
	  graph_degree[*ui] = degree(*ui, graph);
       }
    }
#endif


    template <class DH>
    void
    Cuthill_McKee (DH&              dof_handler,
		   const bool       reversed_numbering,
		   const bool       use_constraints)
    {
      std::vector<unsigned int> renumbering(dof_handler.n_dofs(),
					    DH::invalid_dof_index);
      compute_Cuthill_McKee(renumbering, dof_handler, reversed_numbering,
			    use_constraints);

				       // actually perform renumbering;
				       // this is dimension specific and
				       // thus needs an own function
      dof_handler.renumber_dofs (renumbering);
    }


    template <class DH>
    void
    compute_Cuthill_McKee (std::vector<unsigned int>& new_dof_indices,
			   const DH        &dof_handler,
			   const bool       reversed_numbering,
			   const bool       use_constraints)
    {
#ifdef DEAL_II_BOOST_GRAPH_COMPILER_BUG
      (void)new_dof_indices;
      (void)dof_handler;
      (void)reversed_numbering;
      (void)use_constraints;
      Assert (false,
	      ExcMessage ("Due to a bug in your compiler, this function triggers an internal "
			  "compiler error and has been disabled. If you need to use the "
			  "function, you need to upgrade to a newer version of the compiler."));
#else
      types::Graph
	graph(dof_handler.n_dofs());
      types::property_map<types::Graph,types::vertex_degree_t>::type
	graph_degree;

      internal::create_graph (dof_handler, use_constraints, graph, graph_degree);

      types::property_map<types::Graph, types::vertex_index_t>::type
	index_map = get(::boost::vertex_index, graph);


      std::vector<types::Vertex> inv_perm(num_vertices(graph));

      if (reversed_numbering == false)
	::boost::cuthill_mckee_ordering(graph, inv_perm.rbegin(),
					get(::boost::vertex_color, graph),
					make_degree_map(graph));
      else
	::boost::cuthill_mckee_ordering(graph, inv_perm.begin(),
					get(::boost::vertex_color, graph),
					make_degree_map(graph));

      for (types::size_type c = 0; c != inv_perm.size(); ++c)
	new_dof_indices[index_map[inv_perm[c]]] = c;

      Assert (std::find (new_dof_indices.begin(), new_dof_indices.end(),
			 DH::invalid_dof_index) == new_dof_indices.end(),
	      ExcInternalError());
#endif
    }



    template <class DH>
    void
    king_ordering (DH&              dof_handler,
		   const bool       reversed_numbering,
		   const bool       use_constraints)
    {
      std::vector<unsigned int> renumbering(dof_handler.n_dofs(),
					    DH::invalid_dof_index);
      compute_king_ordering(renumbering, dof_handler, reversed_numbering,
			    use_constraints);

				       // actually perform renumbering;
				       // this is dimension specific and
				       // thus needs an own function
      dof_handler.renumber_dofs (renumbering);
    }


    template <class DH>
    void
    compute_king_ordering (std::vector<unsigned int>& new_dof_indices,
			   const DH        &dof_handler,
			   const bool       reversed_numbering,
			   const bool       use_constraints)
    {
#ifdef DEAL_II_BOOST_GRAPH_COMPILER_BUG
      (void)new_dof_indices;
      (void)dof_handler;
      (void)reversed_numbering;
      (void)use_constraints;
      Assert (false,
	      ExcMessage ("Due to a bug in your compiler, this function triggers an internal "
			  "compiler error and has been disabled. If you need to use the "
			  "function, you need to upgrade to a newer version of the compiler."));
#else
      types::Graph
	graph(dof_handler.n_dofs());
      types::property_map<types::Graph,types::vertex_degree_t>::type
	graph_degree;

      internal::create_graph (dof_handler, use_constraints, graph, graph_degree);

      types::property_map<types::Graph, types::vertex_index_t>::type
	index_map = get(::boost::vertex_index, graph);


      std::vector<types::Vertex> inv_perm(num_vertices(graph));

      if (reversed_numbering == false)
	::boost::king_ordering(graph, inv_perm.rbegin());
      else
	::boost::king_ordering(graph, inv_perm.begin());

      for (types::size_type c = 0; c != inv_perm.size(); ++c)
	new_dof_indices[index_map[inv_perm[c]]] = c;

      Assert (std::find (new_dof_indices.begin(), new_dof_indices.end(),
			 DH::invalid_dof_index) == new_dof_indices.end(),
	      ExcInternalError());
#endif
    }



    template <class DH>
    void
    minimum_degree (DH&              dof_handler,
		    const bool       reversed_numbering,
		    const bool       use_constraints)
    {
      std::vector<unsigned int> renumbering(dof_handler.n_dofs(),
					    DH::invalid_dof_index);
      compute_minimum_degree(renumbering, dof_handler, reversed_numbering,
			     use_constraints);

				       // actually perform renumbering;
				       // this is dimension specific and
				       // thus needs an own function
      dof_handler.renumber_dofs (renumbering);
    }


    template <class DH>
    void
    compute_minimum_degree (std::vector<unsigned int>& new_dof_indices,
			    const DH        &dof_handler,
			    const bool       reversed_numbering,
			    const bool       use_constraints)
    {
      Assert (use_constraints == false, ExcNotImplemented());

				       // the following code is pretty
				       // much a verbatim copy of the
				       // sample code for the
				       // minimum_degree_ordering manual
				       // page from the BOOST Graph
				       // Library
      using namespace ::boost;

      int delta = 0;

      typedef double Type;

				       // must be BGL directed graph now
      typedef adjacency_list<vecS, vecS, directedS>  Graph;
      typedef graph_traits<Graph>::vertex_descriptor Vertex;

      int n = dof_handler.n_dofs();

      Graph G(n);

      std::vector<unsigned int> dofs_on_this_cell;

      typename DH::active_cell_iterator cell = dof_handler.begin_active(),
					endc = dof_handler.end();

      for (; cell!=endc; ++cell)
	{

	  const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

	  dofs_on_this_cell.resize (dofs_per_cell);

	  cell->get_dof_indices (dofs_on_this_cell);
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      if (dofs_on_this_cell[i] > dofs_on_this_cell[j])
		{
		  add_edge (dofs_on_this_cell[i], dofs_on_this_cell[j], G);
		  add_edge (dofs_on_this_cell[j], dofs_on_this_cell[i], G);
		}
	}


      typedef std::vector<int> Vector;


      Vector inverse_perm(n, 0);

      Vector perm(n, 0);


      Vector supernode_sizes(n, 1);
				       // init has to be 1

      ::boost::property_map<Graph, vertex_index_t>::type
	  id = get(vertex_index, G);


      Vector degree(n, 0);


      minimum_degree_ordering
	(G,
	 make_iterator_property_map(&degree[0], id, degree[0]),
	 &inverse_perm[0],
	 &perm[0],
	 make_iterator_property_map(&supernode_sizes[0], id, supernode_sizes[0]),
	 delta, id);


      for (int i=0; i<n; ++i)
	{
	  Assert (std::find (perm.begin(), perm.end(), i)
		  != perm.end(),
		  ExcInternalError());
	  Assert (std::find (inverse_perm.begin(), inverse_perm.end(), i)
		  != inverse_perm.end(),
		  ExcInternalError());
	  Assert (inverse_perm[perm[i]] == i, ExcInternalError());
	}

      if (reversed_numbering == true)
	std::copy (perm.begin(), perm.end(),
		   new_dof_indices.begin());
      else
	std::copy (inverse_perm.begin(), inverse_perm.end(),
		   new_dof_indices.begin());
    }

  }  // namespace boost



  template <class DH>
  void
  Cuthill_McKee (DH&              dof_handler,
		 const bool       reversed_numbering,
		 const bool       use_constraints,
		 const std::vector<unsigned int> &starting_indices)
  {
    std::vector<unsigned int> renumbering(dof_handler.n_dofs(),
					  DH::invalid_dof_index);
    compute_Cuthill_McKee(renumbering, dof_handler, reversed_numbering,
			  use_constraints, starting_indices);

				     // actually perform renumbering;
				     // this is dimension specific and
				     // thus needs an own function
    dof_handler.renumber_dofs (renumbering);
  }



  template <class DH>
  void
  compute_Cuthill_McKee (std::vector<unsigned int>& new_indices,
			 const DH&                  dof_handler,
			 const bool                 reversed_numbering,
			 const bool                 use_constraints,
			 const std::vector<unsigned int>& starting_indices)
  {
				     // make the connection graph. in
				     // more than 2d use an intermediate
				     // compressed sparsity pattern
				     // since the we don't have very
				     // good estimates for
				     // max_couplings_between_dofs() in
				     // 3d and this then leads to
				     // excessive memory consumption
				     //
				     // note that if constraints are not
				     // requested, then the
				     // 'constraints' object will be
				     // empty, and calling condense with
				     // it is a no-op
    ConstraintMatrix constraints;
    if (use_constraints)
      DoFTools::make_hanging_node_constraints (dof_handler, constraints);
    constraints.close ();

    SparsityPattern sparsity;
    if (DH::dimension < 2)
      {
	sparsity.reinit (dof_handler.n_dofs(),
			 dof_handler.n_dofs(),
			 dof_handler.max_couplings_between_dofs());
	DoFTools::make_sparsity_pattern (dof_handler, sparsity, constraints);
	sparsity.compress();
      }
    else
      {
	CompressedSimpleSparsityPattern csp (dof_handler.n_dofs(),
					     dof_handler.n_dofs());
	DoFTools::make_sparsity_pattern (dof_handler, csp, constraints);
	sparsity.copy_from (csp);
      }

				     // constraints are not needed anymore
    constraints.clear ();

    Assert(new_indices.size() == sparsity.n_rows(),
	   ExcDimensionMismatch(new_indices.size(),
				sparsity.n_rows()));

    SparsityTools::reorder_Cuthill_McKee (sparsity, new_indices,
					  starting_indices);

    if (reversed_numbering)
      new_indices = Utilities::reverse_permutation (new_indices);
  }



  template <int dim>
  void Cuthill_McKee (MGDoFHandler<dim>               &dof_handler,
		      const unsigned int               level,
		      const bool                       reversed_numbering,
		      const std::vector<unsigned int> &starting_indices)
  {
//TODO: we should be doing the same here as in the other compute_CMK function to preserve some memory

				     // make the connection graph
    SparsityPattern sparsity (dof_handler.n_dofs(level),
			      dof_handler.max_couplings_between_dofs());
    MGTools::make_sparsity_pattern (dof_handler, sparsity, level);

    std::vector<unsigned int> new_indices(sparsity.n_rows());
    SparsityTools::reorder_Cuthill_McKee (sparsity, new_indices,
					  starting_indices);

    if (reversed_numbering)
      new_indices = Utilities::reverse_permutation (new_indices);

				     // actually perform renumbering;
				     // this is dimension specific and
				     // thus needs an own function
    dof_handler.renumber_dofs (level, new_indices);
  }



  template <int dim, int spacedim>
  void
  component_wise (DoFHandler<dim,spacedim>        &dof_handler,
		  const std::vector<unsigned int> &component_order_arg)
  {
    std::vector<unsigned int> renumbering (dof_handler.n_locally_owned_dofs(),
					   DoFHandler<dim>::invalid_dof_index);

    typedef
      internal::WrapDoFIterator<typename DoFHandler<dim,spacedim>
                                ::active_cell_iterator>
      ITERATOR;

    typename DoFHandler<dim,spacedim>::active_cell_iterator
      istart = dof_handler.begin_active();
    ITERATOR start = istart;
    const typename DoFHandler<dim,spacedim>::cell_iterator
      end = dof_handler.end();

    const unsigned int result =
      compute_component_wise<dim, spacedim, ITERATOR,
      typename DoFHandler<dim,spacedim>::cell_iterator>
      (renumbering, start, end, component_order_arg);
    if (result == 0)
      return;

				     // verify that the last numbered
				     // degree of freedom is either
				     // equal to the number of degrees
				     // of freedom in total (the
				     // sequential case) or in the
				     // distributed case at least
				     // makes sense
    Assert ((result == dof_handler.n_locally_owned_dofs())
	    ||
	    ((dof_handler.n_locally_owned_dofs() < dof_handler.n_dofs())
	     &&
	     (result <= dof_handler.n_dofs())),
	    ExcRenumberingIncomplete());

    dof_handler.renumber_dofs (renumbering);
  }



  template <int dim>
  void
  component_wise (hp::DoFHandler<dim>             &dof_handler,
		  const std::vector<unsigned int> &component_order_arg)
  {
    std::vector<unsigned int> renumbering (dof_handler.n_dofs(),
					   hp::DoFHandler<dim>::invalid_dof_index);

    typedef
      internal::WrapDoFIterator<typename hp::DoFHandler<dim>::active_cell_iterator> ITERATOR;

    typename hp::DoFHandler<dim>::active_cell_iterator
      istart = dof_handler.begin_active();
    ITERATOR start = istart;
    const typename hp::DoFHandler<dim>::cell_iterator
      end = dof_handler.end();

    const unsigned int result =
      compute_component_wise<dim, dim, ITERATOR,
      typename hp::DoFHandler<dim>::cell_iterator>(renumbering,
						   start, end,
						   component_order_arg);

    if (result == 0) return;

    Assert (result == dof_handler.n_dofs(),
	    ExcRenumberingIncomplete());

    dof_handler.renumber_dofs (renumbering);
  }



  template <int dim>
  void
  component_wise (MGDoFHandler<dim> &dof_handler,
		  const unsigned int level,
		  const std::vector<unsigned int> &component_order_arg)
  {
    std::vector<unsigned int> renumbering (dof_handler.n_dofs(level),
					   DoFHandler<dim>::invalid_dof_index);

    typedef
      internal::WrapMGDoFIterator<typename MGDoFHandler<dim>::cell_iterator> ITERATOR;

    typename MGDoFHandler<dim>::cell_iterator
      istart =dof_handler.begin(level);
    ITERATOR start = istart;
    typename MGDoFHandler<dim>::cell_iterator
      iend = dof_handler.end(level);
    const ITERATOR end = iend;

    const unsigned int result =
      compute_component_wise<dim, dim, ITERATOR, ITERATOR>(
	renumbering, start, end, component_order_arg);

    if (result == 0) return;

    Assert (result == dof_handler.n_dofs(level),
	    ExcRenumberingIncomplete());

    if (renumbering.size()!=0)
      dof_handler.renumber_dofs (level, renumbering);
  }



  template <int dim>
  void
  component_wise (MGDoFHandler<dim> &dof_handler,
		  const std::vector<unsigned int> &component_order_arg)
  {
				     // renumber the non-MG part of
				     // the DoFHandler in parallel to
				     // the MG part. Because
				     // MGDoFHandler::renumber_dofs
				     // uses the user flags we can't
				     // run renumbering on individual
				     // levels in parallel to the
				     // other levels
    void (*non_mg_part) (DoFHandler<dim> &, const std::vector<unsigned int> &)
      = &component_wise<dim>;
    Threads::Task<>
      task = Threads::new_task (non_mg_part, dof_handler, component_order_arg);

    for (unsigned int level=0; level<dof_handler.get_tria().n_levels(); ++level)
      component_wise (dof_handler, level, component_order_arg);

    task.join();
  }



  template <int dim, int spacedim, class ITERATOR, class ENDITERATOR>
  unsigned int
  compute_component_wise (std::vector<unsigned int>& new_indices,
			  const ITERATOR   & start,
			  const ENDITERATOR& end,
			  const std::vector<unsigned int> &component_order_arg)
  {
    const hp::FECollection<dim,spacedim>
      fe_collection (start->get_dof_handler().get_fe ());

				     // do nothing if the FE has only
				     // one component
    if (fe_collection.n_components() == 1)
      {
	new_indices.resize(0);
	return 0;
      }

				     // Copy last argument into a
				     // writable vector.
    std::vector<unsigned int> component_order (component_order_arg);
				     // If the last argument was an
				     // empty vector, set up things to
				     // store components in the order
				     // found in the system.
    if (component_order.size() == 0)
      for (unsigned int i=0; i<fe_collection.n_components(); ++i)
	component_order.push_back (i);

    Assert (component_order.size() == fe_collection.n_components(),
	    ExcDimensionMismatch(component_order.size(), fe_collection.n_components()));

    for (unsigned int i=0; i<component_order.size(); ++i)
      Assert(component_order[i] < fe_collection.n_components(),
	     ExcIndexRange(component_order[i], 0, fe_collection.n_components()));

				     // vector to hold the dof indices on
				     // the cell we visit at a time
    std::vector<unsigned int> local_dof_indices;

				     // prebuilt list to which component
				     // a given dof on a cell
				     // should go. note that we get into
				     // trouble here if the shape
				     // function is not primitive, since
				     // then there is no single vector
				     // component to which it
				     // belongs. in this case, assign it
				     // to the first vector component to
				     // which it belongs
    std::vector<std::vector<unsigned int> > component_list (fe_collection.size());
    for (unsigned int f=0; f<fe_collection.size(); ++f)
      {
	const FiniteElement<dim,spacedim> & fe = fe_collection[f];
	const unsigned int dofs_per_cell = fe.dofs_per_cell;
	component_list[f].resize(dofs_per_cell);
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  if (fe.is_primitive(i))
	    component_list[f][i]
	      = component_order[fe.system_to_component_index(i).first];
	  else
	    {
	      const unsigned int comp = (std::find(fe.get_nonzero_components(i).begin(),
						   fe.get_nonzero_components(i).end(),
						   true) -
					 fe.get_nonzero_components(i).begin());

					   // then associate this degree
					   // of freedom with this
					   // component
	      component_list[f][i] = component_order[comp];
	    }
      }

				     // set up a map where for each
				     // component the respective degrees
				     // of freedom are collected.
				     //
				     // note that this map is sorted by
				     // component but that within each
				     // component it is NOT sorted by
				     // dof index. note also that some
				     // dof indices are entered
				     // multiply, so we will have to
				     // take care of that
    std::vector<std::vector<unsigned int> >
      component_to_dof_map (fe_collection.n_components());
    for (ITERATOR cell=start; cell!=end; ++cell)
      if (cell->is_locally_owned())
	{
					 // on each cell: get dof indices
					 // and insert them into the global
					 // list using their component
	  const unsigned int fe_index = cell->active_fe_index();
	  const unsigned int dofs_per_cell =fe_collection[fe_index].dofs_per_cell;
	  local_dof_indices.resize (dofs_per_cell);
	  cell.get_dof_indices (local_dof_indices);
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    if (start->get_dof_handler().locally_owned_dofs().is_element(local_dof_indices[i]))
	      component_to_dof_map[component_list[fe_index][i]].
		push_back (local_dof_indices[i]);
	}

				     // now we've got all indices sorted
				     // into buckets labelled with their
				     // target component number. we've
				     // only got to traverse this list
				     // and assign the new indices
				     //
				     // however, we first want to sort
				     // the indices entered into the
				     // buckets to preserve the order
				     // within each component and during
				     // this also remove duplicate
				     // entries
				     //
				     // note that we no longer have to
				     // care about non-primitive shape
				     // functions since the buckets
				     // corresponding to the second and
				     // following vector components of a
				     // non-primitive FE will simply be
				     // empty, everything being shoved
				     // into the first one. The same
				     // holds if several components were
				     // joined into a single target.
    for (unsigned int component=0; component<fe_collection.n_components();
	 ++component)
      {
	std::sort (component_to_dof_map[component].begin(),
		   component_to_dof_map[component].end());
	component_to_dof_map[component]
	  .erase (std::unique (component_to_dof_map[component].begin(),
			       component_to_dof_map[component].end()),
		  component_to_dof_map[component].end());
      }

				     // calculate the number of locally owned
				     // DoFs per bucket
    const unsigned int n_buckets = fe_collection.n_components();
    std::vector<unsigned int> shifts(n_buckets);

    if (const parallel::distributed::Triangulation<dim,spacedim> * tria
	= (dynamic_cast<const parallel::distributed::Triangulation<dim,spacedim>*>
	   (&start->get_dof_handler().get_tria())))
      {
#ifdef DEAL_II_USE_P4EST
	std::vector<unsigned int> local_dof_count(n_buckets);

	for (unsigned int c=0; c<n_buckets; ++c)
	  local_dof_count[c] = component_to_dof_map[c].size();


					 // gather information from all CPUs
	std::vector<unsigned int>
	  all_dof_counts(fe_collection.n_components() *
			 Utilities::System::get_n_mpi_processes (tria->get_communicator()));

	MPI_Allgather ( &local_dof_count[0], n_buckets, MPI_UNSIGNED, &all_dof_counts[0],
			n_buckets, MPI_UNSIGNED, tria->get_communicator());

	for (unsigned int i=0; i<n_buckets; ++i)
	  Assert (all_dof_counts[n_buckets*tria->locally_owned_subdomain()+i]
		  ==
		  local_dof_count[i],
		  ExcInternalError());

					 //calculate shifts
	unsigned int cumulated = 0;
	for (unsigned int c=0; c<n_buckets; ++c)
	  {
	    shifts[c]=cumulated;
	    for (types::subdomain_id_t i=0; i<tria->locally_owned_subdomain(); ++i)
	      shifts[c] += all_dof_counts[c+n_buckets*i];
	    for (unsigned int i=0; i<Utilities::System::get_n_mpi_processes (tria->get_communicator()); ++i)
	      cumulated += all_dof_counts[c+n_buckets*i];
	  }
#else
	(void)tria;
	Assert (false, ExcInternalError());
#endif
      }
    else
      {
	shifts[0] = 0;
	for (unsigned int c=1; c<fe_collection.n_components(); ++c)
	  shifts[c] = shifts[c-1] + component_to_dof_map[c-1].size();
      }




				     // now concatenate all the
				     // components in the order the user
				     // desired to see
    unsigned int next_free_index = 0;
    for (unsigned int component=0; component<fe_collection.n_components(); ++component)
      {
	const typename std::vector<unsigned int>::const_iterator
	  begin_of_component = component_to_dof_map[component].begin(),
	  end_of_component   = component_to_dof_map[component].end();

	next_free_index = shifts[component];

	for (typename std::vector<unsigned int>::const_iterator
	       dof_index = begin_of_component;
	     dof_index != end_of_component; ++dof_index)
	  {
	    Assert (start->get_dof_handler().locally_owned_dofs()
		    .index_within_set(*dof_index)
		    <
		    new_indices.size(),
		    ExcInternalError());
	    new_indices[start->get_dof_handler().locally_owned_dofs()
			.index_within_set(*dof_index)]
	      = next_free_index++;
	  }
      }

    return next_free_index;
  }

  namespace
  {
	// helper function for hierarchical()
    template <int dim, class iterator>
    unsigned int
    compute_hierarchical_recursive (
      unsigned int next_free,
      std::vector<unsigned int>& new_indices,
      const iterator & cell,
      const IndexSet & locally_owned)
    {
      if (cell->has_children())
        {
		  //recursion
          for (unsigned int c = 0;c < GeometryInfo<dim>::max_children_per_cell; ++c)
            next_free = compute_hierarchical_recursive<dim> (
                          next_free,
                          new_indices,
                          cell->child (c),
                          locally_owned);
        }
      else
        {
          if (cell->is_locally_owned())
            {
              const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
              std::vector<unsigned int> local_dof_indices (dofs_per_cell);
              cell->get_dof_indices (local_dof_indices);

              for (unsigned int i = 0;i < dofs_per_cell;++i)
                {
                  if (locally_owned.is_element (local_dof_indices[i]))
                    {
                      // this is a locally owned DoF, assign new number if not assigned a number yet
                      unsigned int idx = locally_owned.index_within_set (local_dof_indices[i]);
                      if (new_indices[idx] == DoFHandler<dim>::invalid_dof_index)
                        {
                          new_indices[idx] = locally_owned.nth_index_in_set (next_free);
                          next_free++;
                        }
                    }
                }
            }
        }
      return next_free;
    }
  }



  template <int dim>
  void
  hierarchical (DoFHandler<dim> &dof_handler)
  {
    std::vector<unsigned int> renumbering (dof_handler.n_locally_owned_dofs(),
					   DoFHandler<dim>::invalid_dof_index);

	typename DoFHandler<dim>::cell_iterator cell;

	unsigned int next_free = 0;
	const IndexSet locally_owned = dof_handler.locally_owned_dofs();

    const parallel::distributed::Triangulation<dim> * tria
	= dynamic_cast<const parallel::distributed::Triangulation<dim>*>
	   (&dof_handler.get_tria());

	if (tria)
      {
#ifdef DEAL_II_USE_P4EST
        //this is a distributed Triangulation. We need to traverse the coarse
        //cells in the order p4est does
        for (unsigned int c = 0; c < tria->n_cells (0); ++c)
          {
            unsigned int coarse_cell_index =
              tria->get_p4est_tree_to_coarse_cell_permutation() [c];

            const typename DoFHandler<dim>::cell_iterator
            cell (tria, 0, coarse_cell_index, &dof_handler);

            next_free = compute_hierarchical_recursive<dim> (next_free,
                        renumbering,
                        cell,
                        locally_owned);
          }
#else
        Assert (false, ExcNotImplemented());
#endif
      }
    else
      {
        //this is not a distributed Triangulation. Traverse coarse cells in the
        //normal order
        for (cell = dof_handler.begin (0); cell != dof_handler.end (0); ++cell)
          next_free = compute_hierarchical_recursive<dim> (next_free,
                      renumbering,
                      cell,
                      locally_owned);
      }

				     // verify that the last numbered
				     // degree of freedom is either
				     // equal to the number of degrees
				     // of freedom in total (the
				     // sequential case) or in the
				     // distributed case at least
				     // makes sense
    Assert ((next_free == dof_handler.n_locally_owned_dofs())
	    ||
	    ((dof_handler.n_locally_owned_dofs() < dof_handler.n_dofs())
	     &&
	     (next_free <= dof_handler.n_dofs())),
	    ExcRenumberingIncomplete());

	// make sure that all local DoFs got new numbers assigned
	Assert (std::find (renumbering.begin(), renumbering.end(),
		       numbers::invalid_unsigned_int)
	    == renumbering.end(),
	    ExcInternalError());

	dof_handler.renumber_dofs(renumbering);
  }



  template <class DH>
  void
  sort_selected_dofs_back (DH&                      dof_handler,
			   const std::vector<bool>& selected_dofs)
  {
    std::vector<unsigned int> renumbering(dof_handler.n_dofs(),
					  DH::invalid_dof_index);
    compute_sort_selected_dofs_back(renumbering, dof_handler, selected_dofs);

    dof_handler.renumber_dofs(renumbering);
  }



  template <class DH>
  void
  compute_sort_selected_dofs_back (std::vector<unsigned int>& new_indices,
				   const DH&                  dof_handler,
				   const std::vector<bool>&   selected_dofs)
  {
    const unsigned int n_dofs = dof_handler.n_dofs();
    Assert (selected_dofs.size() == n_dofs,
	    ExcDimensionMismatch (selected_dofs.size(), n_dofs));

				     // re-sort the dofs according to
				     // their selection state
    Assert (new_indices.size() == n_dofs,
	    ExcDimensionMismatch(new_indices.size(), n_dofs));

    const unsigned int   n_selected_dofs = std::count (selected_dofs.begin(),
						       selected_dofs.end(),
						       false);

    unsigned int next_unselected = 0;
    unsigned int next_selected   = n_selected_dofs;
    for (unsigned int i=0; i<n_dofs; ++i)
      if (selected_dofs[i] == false)
	{
	  new_indices[i] = next_unselected;
	  ++next_unselected;
	}
      else
	{
	  new_indices[i] = next_selected;
	  ++next_selected;
	};
    Assert (next_unselected == n_selected_dofs, ExcInternalError());
    Assert (next_selected == n_dofs, ExcInternalError());
  }



  template <class DH>
  void
  cell_wise_dg (DH& dof,
		const std::vector<typename DH::cell_iterator>& cells)
  {
    std::vector<unsigned int> renumbering(dof.n_dofs());
    std::vector<unsigned int> reverse(dof.n_dofs());
    compute_cell_wise_dg(renumbering, reverse, dof, cells);

    dof.renumber_dofs(renumbering);
  }



  template <class DH>
  void
  cell_wise (
    DH& dof,
    const std::vector<typename DH::cell_iterator>& cells)
  {
    std::vector<unsigned int> renumbering(dof.n_dofs());
    std::vector<unsigned int> reverse(dof.n_dofs());
    compute_cell_wise(renumbering, reverse, dof, cells);

    dof.renumber_dofs(renumbering);
  }


//TODO: Discuss if we cannot replace this function by the next
  template <class DH>
  void
  compute_cell_wise_dg (
    std::vector<unsigned int>& new_indices,
    std::vector<unsigned int>& reverse,
    const DH& dof,
    const typename std::vector<typename DH::cell_iterator>& cells)
  {
    Assert(cells.size() == dof.get_tria().n_active_cells(),
	   ExcDimensionMismatch(cells.size(),
				dof.get_tria().n_active_cells()));

    unsigned int n_global_dofs = dof.n_dofs();

				     // Actually, we compute the
				     // inverse of the reordering
				     // vector, called reverse here.
				     // Later, its inverse is computed
				     // into new_indices, which is the
				     // return argument.

    Assert(new_indices.size() == n_global_dofs,
	   ExcDimensionMismatch(new_indices.size(), n_global_dofs));
    Assert(reverse.size() == n_global_dofs,
	   ExcDimensionMismatch(reverse.size(), n_global_dofs));

    std::vector<unsigned int> cell_dofs;

    unsigned int global_index = 0;

    typename std::vector<typename DH::cell_iterator>::const_iterator cell;

    for(cell = cells.begin(); cell != cells.end(); ++cell)
      {
	Assert((*cell)->get_fe().n_dofs_per_face()==0, ExcNotDGFEM());
					 // Determine the number of dofs
					 // on this cell and reinit the
					 // vector storing these
					 // numbers.
	unsigned int n_cell_dofs = (*cell)->get_fe().n_dofs_per_cell();
	cell_dofs.resize(n_cell_dofs);

	(*cell)->get_dof_indices(cell_dofs);

					 // Sort here to make sure that
					 // degrees of freedom inside a
					 // single cell are in the same
					 // order after renumbering.
	std::sort(cell_dofs.begin(), cell_dofs.end());

	for (unsigned int i=0;i<n_cell_dofs;++i)
	  {
	    reverse[global_index++] = cell_dofs[i];
	  }
      }
    Assert(global_index == n_global_dofs, ExcRenumberingIncomplete());

    for (unsigned int i=0;i<reverse.size(); ++i)
      new_indices[reverse[i]] = i;
  }


  template <class DH>
  void
  compute_cell_wise (
    std::vector<unsigned int>& new_indices,
    std::vector<unsigned int>& reverse,
    const DH& dof,
    const typename std::vector<typename DH::cell_iterator>& cells)
  {
    Assert(cells.size() == dof.get_tria().n_active_cells(),
	   ExcDimensionMismatch(cells.size(),
				dof.get_tria().n_active_cells()));

    unsigned int n_global_dofs = dof.n_dofs();

				     // Actually, we compute the
				     // inverse of the reordering
				     // vector, called reverse here.
				     // Later, irs inverse is computed
				     // into new_indices, which is the
				     // return argument.

    Assert(new_indices.size() == n_global_dofs,
	   ExcDimensionMismatch(new_indices.size(), n_global_dofs));
    Assert(reverse.size() == n_global_dofs,
	   ExcDimensionMismatch(reverse.size(), n_global_dofs));

				     // For continuous elements, we must
				     // make sure, that each dof is
				     // reordered only once.
    std::vector<bool> already_sorted(n_global_dofs, false);
    std::vector<unsigned int> cell_dofs;

    unsigned int global_index = 0;

    typename std::vector<typename DH::cell_iterator>::const_iterator cell;

    for(cell = cells.begin(); cell != cells.end(); ++cell)
      {
					 // Determine the number of dofs
					 // on this cell and reinit the
					 // vector storing these
					 // numbers.
	unsigned int n_cell_dofs = (*cell)->get_fe().n_dofs_per_cell();
	cell_dofs.resize(n_cell_dofs);

	(*cell)->get_dof_indices(cell_dofs);

					 // Sort here to make sure that
					 // degrees of freedom inside a
					 // single cell are in the same
					 // order after renumbering.
	std::sort(cell_dofs.begin(), cell_dofs.end());

	for (unsigned int i=0;i<n_cell_dofs;++i)
	  {
	    if (!already_sorted[cell_dofs[i]])
	      {
		already_sorted[cell_dofs[i]] = true;
		reverse[global_index++] = cell_dofs[i];
	      }
	  }
      }
    Assert(global_index == n_global_dofs, ExcRenumberingIncomplete());

    for (unsigned int i=0;i<reverse.size(); ++i)
      new_indices[reverse[i]] = i;
  }



  template <int dim>
  void cell_wise_dg (
    MGDoFHandler<dim>& dof,
    const unsigned int level,
    const typename std::vector<typename MGDoFHandler<dim>::cell_iterator>& cells)
  {
    std::vector<unsigned int> renumbering(dof.n_dofs(level));
    std::vector<unsigned int> reverse(dof.n_dofs(level));

    compute_cell_wise_dg(renumbering, reverse, dof, level, cells);
    dof.renumber_dofs(level, renumbering);
  }



  template <int dim>
  void cell_wise (
    MGDoFHandler<dim>& dof,
    const unsigned int level,
    const typename std::vector<typename MGDoFHandler<dim>::cell_iterator>& cells)
  {
    std::vector<unsigned int> renumbering(dof.n_dofs(level));
    std::vector<unsigned int> reverse(dof.n_dofs(level));

    compute_cell_wise(renumbering, reverse, dof, level, cells);
    dof.renumber_dofs(level, renumbering);
  }



  template <int dim>
  void compute_cell_wise_dg (
    std::vector<unsigned int>& new_order,
    std::vector<unsigned int>& reverse,
    const MGDoFHandler<dim>& dof,
    const unsigned int level,
    const typename std::vector<typename MGDoFHandler<dim>::cell_iterator>& cells)
  {
    Assert(cells.size() == dof.get_tria().n_cells(level),
	   ExcDimensionMismatch(cells.size(),
				dof.get_tria().n_cells(level)));
    switch (dim)
      {
	case 3:
	      Assert(dof.get_fe().n_dofs_per_quad()==0,
		     ExcNotDGFEM());
	case 2:
	      Assert(dof.get_fe().n_dofs_per_line()==0,
		     ExcNotDGFEM());
	default:
	      Assert(dof.get_fe().n_dofs_per_vertex()==0,
		     ExcNotDGFEM());
      }

    Assert (new_order.size() == dof.n_dofs(level),
	    ExcDimensionMismatch(new_order.size(), dof.n_dofs(level)));
    Assert (reverse.size() == dof.n_dofs(level),
	    ExcDimensionMismatch(reverse.size(), dof.n_dofs(level)));

    unsigned int n_global_dofs = dof.n_dofs(level);
    unsigned int n_cell_dofs = dof.get_fe().n_dofs_per_cell();

    std::vector<unsigned int> cell_dofs(n_cell_dofs);

    unsigned int global_index = 0;

    typename std::vector<typename MGDoFHandler<dim>::cell_iterator>::const_iterator cell;

    for(cell = cells.begin(); cell != cells.end(); ++cell)
      {
	Assert ((*cell)->level() == (int) level, ExcInternalError());

	(*cell)->get_mg_dof_indices(cell_dofs);
	std::sort(cell_dofs.begin(), cell_dofs.end());

	for (unsigned int i=0;i<n_cell_dofs;++i)
	  {
	    reverse[global_index++] = cell_dofs[i];
	  }
      }
    Assert(global_index == n_global_dofs, ExcRenumberingIncomplete());

    for (unsigned int i=0;i<new_order.size(); ++i)
      new_order[reverse[i]] = i;
  }



  template <int dim>
  void compute_cell_wise (
    std::vector<unsigned int>& new_order,
    std::vector<unsigned int>& reverse,
    const MGDoFHandler<dim>& dof,
    const unsigned int level,
    const typename std::vector<typename MGDoFHandler<dim>::cell_iterator>& cells)
  {
    Assert(cells.size() == dof.get_tria().n_cells(level),
	   ExcDimensionMismatch(cells.size(),
				dof.get_tria().n_cells(level)));
    Assert (new_order.size() == dof.n_dofs(level),
	    ExcDimensionMismatch(new_order.size(), dof.n_dofs(level)));
    Assert (reverse.size() == dof.n_dofs(level),
	    ExcDimensionMismatch(reverse.size(), dof.n_dofs(level)));

    unsigned int n_global_dofs = dof.n_dofs(level);
    unsigned int n_cell_dofs = dof.get_fe().n_dofs_per_cell();

    std::vector<bool> already_sorted(n_global_dofs, false);
    std::vector<unsigned int> cell_dofs(n_cell_dofs);

    unsigned int global_index = 0;

    typename std::vector<typename MGDoFHandler<dim>::cell_iterator>::const_iterator cell;

    for(cell = cells.begin(); cell != cells.end(); ++cell)
      {
	Assert ((*cell)->level() == (int) level, ExcInternalError());

	(*cell)->get_mg_dof_indices(cell_dofs);
	std::sort(cell_dofs.begin(), cell_dofs.end());

	for (unsigned int i=0;i<n_cell_dofs;++i)
	  {
	    if (!already_sorted[cell_dofs[i]])
	      {
		already_sorted[cell_dofs[i]] = true;
		reverse[global_index++] = cell_dofs[i];
	      }
	  }
      }
    Assert(global_index == n_global_dofs, ExcRenumberingIncomplete());

    for (unsigned int i=0;i<new_order.size(); ++i)
      new_order[reverse[i]] = i;
  }



  template <class DH, int dim>
  void
  downstream_dg (DH& dof, const Point<dim>& direction)
  {
    std::vector<unsigned int> renumbering(dof.n_dofs());
    compute_downstream_dg(renumbering, dof, direction);

    dof.renumber_dofs(renumbering);
  }



  template <class DH, int dim>
  void
  downstream (DH& dof, const Point<dim>& direction,
	      const bool dof_wise_renumbering)
  {
    std::vector<unsigned int> renumbering(dof.n_dofs());
    std::vector<unsigned int> reverse(dof.n_dofs());
    compute_downstream(renumbering, reverse, dof, direction,
		       dof_wise_renumbering);

    dof.renumber_dofs(renumbering);
  }



  template <class DH, int dim>
  void
  compute_downstream_dg (
    std::vector<unsigned int>& new_indices,
    const DH& dof,
    const Point<dim>& direction)
  {
    std::vector<typename DH::cell_iterator>
      ordered_cells(dof.get_tria().n_active_cells());
    const CompareDownstream<typename DH::cell_iterator, dim> comparator(direction);

    typename DH::active_cell_iterator begin = dof.begin_active();
    typename DH::active_cell_iterator end = dof.end();

    copy (begin, end, ordered_cells.begin());
    std::sort (ordered_cells.begin(), ordered_cells.end(), comparator);

    std::vector<unsigned int> reverse(new_indices.size());
    compute_cell_wise_dg(new_indices, reverse, dof, ordered_cells);
  }



  template <class DH, int dim>
  void
  compute_downstream_dg (
    std::vector<unsigned int>& new_indices,
    std::vector<unsigned int>& reverse,
    const DH& dof,
    const Point<dim>& direction)
  {
    std::vector<typename DH::cell_iterator>
      ordered_cells(dof.get_tria().n_active_cells());
    const CompareDownstream<typename DH::cell_iterator, dim> comparator(direction);

    typename DH::active_cell_iterator begin = dof.begin_active();
    typename DH::active_cell_iterator end = dof.end();

    copy (begin, end, ordered_cells.begin());
    std::sort (ordered_cells.begin(), ordered_cells.end(), comparator);

    compute_cell_wise_dg(new_indices, reverse, dof, ordered_cells);
  }


  template <class DH, int dim>
  void
  compute_downstream (
    std::vector<unsigned int>& new_indices,
    std::vector<unsigned int>& reverse,
    const DH& dof,
    const Point<dim>& direction,
    const bool dof_wise_renumbering)
  {
    if (dof_wise_renumbering == false)
      {
	std::vector<typename DH::cell_iterator>
	  ordered_cells(dof.get_tria().n_active_cells());
	const CompareDownstream<typename DH::cell_iterator, dim> comparator(direction);

	typename DH::active_cell_iterator begin = dof.begin_active();
	typename DH::active_cell_iterator end = dof.end();

	copy (begin, end, ordered_cells.begin());
	std::sort (ordered_cells.begin(), ordered_cells.end(), comparator);

	compute_cell_wise(new_indices, reverse, dof, ordered_cells);
      }
    else
      {
				// similar code as for
				// DoFTools::map_dofs_to_support_points, but
				// need to do this for general DH classes and
				// want to be able to sort the result
				// (otherwise, could use something like
				// DoFTools::map_support_points_to_dofs)
	const unsigned int n_dofs = dof.n_dofs();
	std::vector<std::pair<Point<dim>,unsigned int> > support_point_list
	  (n_dofs);

	const hp::FECollection<dim> fe_collection (dof.get_fe ());
	Assert (fe_collection[0].has_support_points(),
		DoFTools::ExcFEHasNoSupportPoints());
	hp::QCollection<dim> quadrature_collection;
	for (unsigned int comp=0; comp<fe_collection.size(); ++comp)
	  {
	    Assert (fe_collection[comp].has_support_points(),
		    DoFTools::ExcFEHasNoSupportPoints());
	    quadrature_collection.push_back
	      (Quadrature<DH::dimension> (fe_collection[comp].
					  get_unit_support_points()));
	  }
	hp::FEValues<DH::dimension,DH::space_dimension>
	  hp_fe_values (fe_collection, quadrature_collection,
			update_quadrature_points);

	std::vector<bool> already_touched (n_dofs, false);

	std::vector<unsigned int> local_dof_indices;
	typename DH::active_cell_iterator begin = dof.begin_active();
	typename DH::active_cell_iterator end = dof.end();
	for ( ; begin != end; ++begin)
	  {
	    const unsigned int dofs_per_cell = begin->get_fe().dofs_per_cell;
	    local_dof_indices.resize (dofs_per_cell);
	    hp_fe_values.reinit (begin);
	    const FEValues<dim> &fe_values =
	      hp_fe_values.get_present_fe_values ();
	    begin->get_dof_indices(local_dof_indices);
	    const std::vector<Point<DH::space_dimension> > & points
	      = fe_values.get_quadrature_points ();
	    for (unsigned int i=0; i<dofs_per_cell; ++i)
	      if (!already_touched[local_dof_indices[i]])
		{
		  support_point_list[local_dof_indices[i]].first = points[i];
		  support_point_list[local_dof_indices[i]].second =
		    local_dof_indices[i];
		  already_touched[local_dof_indices[i]] = true;
		}
	  }

	ComparePointwiseDownstream<dim> comparator (direction);
	std::sort (support_point_list.begin(), support_point_list.end(),
		   comparator);
	for (unsigned int i=0; i<n_dofs; ++i)
	  new_indices[support_point_list[i].second] = i;
      }
  }



  template <int dim>
  void downstream_dg (MGDoFHandler<dim>& dof,
		      const unsigned int level,
		      const Point<dim>&  direction)
  {
    std::vector<unsigned int> renumbering(dof.n_dofs(level));
    std::vector<unsigned int> reverse(dof.n_dofs(level));
    compute_downstream_dg(renumbering, reverse, dof, level, direction);

    dof.renumber_dofs(level, renumbering);
  }



  template <int dim>
  void downstream (MGDoFHandler<dim>& dof,
		   const unsigned int level,
		   const Point<dim>&  direction,
		   const bool         dof_wise_renumbering)
  {
    std::vector<unsigned int> renumbering(dof.n_dofs(level));
    std::vector<unsigned int> reverse(dof.n_dofs(level));
    compute_downstream(renumbering, reverse, dof, level, direction,
		       dof_wise_renumbering);

    dof.renumber_dofs(level, renumbering);
  }



  template <int dim>
  void
  compute_downstream_dg (
    std::vector<unsigned int>& new_indices,
    std::vector<unsigned int>& reverse,
    const MGDoFHandler<dim>& dof,
    const unsigned int level,
    const Point<dim>& direction)
  {
    std::vector<typename MGDoFHandler<dim>::cell_iterator>
      ordered_cells(dof.get_tria().n_cells(level));
    const CompareDownstream<typename MGDoFHandler<dim>::cell_iterator, dim>
      comparator(direction);

    typename MGDoFHandler<dim>::cell_iterator begin = dof.begin(level);
    typename MGDoFHandler<dim>::cell_iterator end = dof.end(level);

    std::copy (begin, end, ordered_cells.begin());
    std::sort (ordered_cells.begin(), ordered_cells.end(), comparator);

    compute_cell_wise_dg(new_indices, reverse, dof, level, ordered_cells);
  }



  template <int dim>
  void
  compute_downstream (
    std::vector<unsigned int>& new_indices,
    std::vector<unsigned int>& reverse,
    const MGDoFHandler<dim>& dof,
    const unsigned int level,
    const Point<dim>& direction,
    const bool dof_wise_renumbering)
  {
    if (dof_wise_renumbering == false)
      {
	std::vector<typename MGDoFHandler<dim>::cell_iterator>
	  ordered_cells(dof.get_tria().n_cells(level));
	const CompareDownstream<typename MGDoFHandler<dim>::cell_iterator, dim> comparator(direction);

	typename MGDoFHandler<dim>::cell_iterator begin = dof.begin(level);
	typename MGDoFHandler<dim>::cell_iterator end = dof.end(level);

	std::copy (begin, end, ordered_cells.begin());
	std::sort (ordered_cells.begin(), ordered_cells.end(), comparator);

	compute_cell_wise(new_indices, reverse, dof, level, ordered_cells);
      }
    else
      {
	Assert (dof.get_fe().has_support_points(),
		DoFTools::ExcFEHasNoSupportPoints());
	const unsigned int n_dofs = dof.n_dofs(level);
	std::vector<std::pair<Point<dim>,unsigned int> > support_point_list
	  (n_dofs);

	Quadrature<dim>   q_dummy(dof.get_fe().get_unit_support_points());
	FEValues<dim,dim> fe_values (dof.get_fe(), q_dummy,
				     update_quadrature_points);

	std::vector<bool> already_touched (dof.n_dofs(), false);

	const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
	std::vector<unsigned int> local_dof_indices (dofs_per_cell);
	typename MGDoFHandler<dim>::cell_iterator begin = dof.begin(level);
	typename MGDoFHandler<dim>::cell_iterator end = dof.end(level);
	for ( ; begin != end; ++begin)
	  {
	    begin->get_mg_dof_indices(local_dof_indices);
	    fe_values.reinit (begin);
	    const std::vector<Point<dim> > & points
	      = fe_values.get_quadrature_points ();
	    for (unsigned int i=0; i<dofs_per_cell; ++i)
	      if (!already_touched[local_dof_indices[i]])
		{
		  support_point_list[local_dof_indices[i]].first = points[i];
		  support_point_list[local_dof_indices[i]].second =
		    local_dof_indices[i];
		  already_touched[local_dof_indices[i]] = true;
		}
	  }

	ComparePointwiseDownstream<dim> comparator (direction);
	std::sort (support_point_list.begin(), support_point_list.end(),
		   comparator);
	for (unsigned int i=0; i<n_dofs; ++i)
	  new_indices[support_point_list[i].second] = i;
      }
  }



/**
 * Provide comparator for DoFCellAccessors
 */
  namespace internal
  {
    template <int dim>
    struct ClockCells
    {
					 /**
					  * Center of rotation.
					  */
	const Point<dim>& center;
					 /**
					  * Revert sorting order.
					  */
	bool counter;

					 /**
					  * Constructor.
					  */
	ClockCells (const Point<dim>& center, bool counter) :
			center(center),
			counter(counter)
	  {}
					 /**
					  * Comparison operator
					  */
	template <class DHCellIterator>
	bool operator () (const DHCellIterator& c1,
			  const DHCellIterator& c2) const
	  {

	    const Point<dim> v1 = c1->center() - center;
	    const Point<dim> v2 = c2->center() - center;
	    const double s1 = std::atan2(v1(0), v1(1));
	    const double s2 = std::atan2(v2(0), v2(1));
	    return ( counter ? (s1>s2) : (s2>s1));
	  }
    };
  }



  template <class DH, int dim>
  void
  clockwise_dg (
    DH& dof,
    const Point<dim>& center,
    const bool counter)
  {
    std::vector<unsigned int> renumbering(dof.n_dofs());
    compute_clockwise_dg(renumbering, dof, center, counter);

    dof.renumber_dofs(renumbering);
  }



  template <class DH, int dim>
  void
  compute_clockwise_dg (
    std::vector<unsigned int>& new_indices,
    const DH& dof,
    const Point<dim>& center,
    const bool counter)
  {
    std::vector<typename DH::cell_iterator>
      ordered_cells(dof.get_tria().n_active_cells());
    internal::ClockCells<dim> comparator(center, counter);

    typename DH::active_cell_iterator begin = dof.begin_active();
    typename DH::active_cell_iterator end = dof.end();

    std::copy (begin, end, ordered_cells.begin());
    std::sort (ordered_cells.begin(), ordered_cells.end(), comparator);

    std::vector<unsigned int> reverse(new_indices.size());
    compute_cell_wise_dg(new_indices, reverse, dof, ordered_cells);
  }



  template <int dim>
  void clockwise_dg (MGDoFHandler<dim>& dof,
		     const unsigned int level,
		     const Point<dim>& center,
		     const bool counter)
  {
    std::vector<typename MGDoFHandler<dim>::cell_iterator>
      ordered_cells(dof.get_tria().n_cells(level));
    internal::ClockCells<dim> comparator(center, counter);

    typename MGDoFHandler<dim>::cell_iterator begin = dof.begin(level);
    typename MGDoFHandler<dim>::cell_iterator end = dof.end(level);

    std::copy (begin, end, ordered_cells.begin());
    std::sort (ordered_cells.begin(), ordered_cells.end(), comparator);

    cell_wise_dg(dof, level, ordered_cells);
  }



  template <class DH>
  void
  random (DH& dof_handler)
  {
    std::vector<unsigned int> renumbering(dof_handler.n_dofs(),
					  DH::invalid_dof_index);
    compute_random(renumbering, dof_handler);

    dof_handler.renumber_dofs(renumbering);
  }



  template <class DH>
  void
  compute_random (
    std::vector<unsigned int>& new_indices,
    const DH&                  dof_handler)
  {
    const unsigned int n_dofs = dof_handler.n_dofs();
    Assert(new_indices.size() == n_dofs,
	   ExcDimensionMismatch(new_indices.size(), n_dofs));

    for (unsigned i=0; i<n_dofs; ++i)
      new_indices[i] = i;

    std::random_shuffle (new_indices.begin(), new_indices.end());
  }



  template <class DH>
  void
  subdomain_wise (DH &dof_handler)
  {
    std::vector<unsigned int> renumbering(dof_handler.n_dofs(),
					  DH::invalid_dof_index);
    compute_subdomain_wise(renumbering, dof_handler);

    dof_handler.renumber_dofs(renumbering);
  }



  template <class DH>
  void
  compute_subdomain_wise (std::vector<unsigned int> &new_dof_indices,
			  const DH                  &dof_handler)
  {
    const unsigned int n_dofs = dof_handler.n_dofs();
    Assert (new_dof_indices.size() == n_dofs,
	    ExcDimensionMismatch (new_dof_indices.size(), n_dofs));

				     // first get the association of each dof
				     // with a subdomain and determine the total
				     // number of subdomain ids used
    std::vector<types::subdomain_id_t> subdomain_association (n_dofs);
    DoFTools::get_subdomain_association (dof_handler,
					 subdomain_association);
    const unsigned int n_subdomains
      = *std::max_element (subdomain_association.begin(),
			   subdomain_association.end()) + 1;

				     // then renumber the subdomains by first
				     // looking at those belonging to subdomain
				     // 0, then those of subdomain 1, etc. note
				     // that the algorithm is stable, i.e. if
				     // two dofs i,j have i<j and belong to the
				     // same subdomain, then they will be in
				     // this order also after reordering
    std::fill (new_dof_indices.begin(), new_dof_indices.end(),
	       numbers::invalid_unsigned_int);
    unsigned int next_free_index = 0;
    for (unsigned int subdomain=0; subdomain<n_subdomains; ++subdomain)
      for (unsigned int i=0; i<n_dofs; ++i)
	if (subdomain_association[i] == subdomain)
	  {
	    Assert (new_dof_indices[i] == numbers::invalid_unsigned_int,
		    ExcInternalError());
	    new_dof_indices[i] = next_free_index;
	    ++next_free_index;
	  }

				     // we should have numbered all dofs
    Assert (next_free_index == n_dofs, ExcInternalError());
    Assert (std::find (new_dof_indices.begin(), new_dof_indices.end(),
		       numbers::invalid_unsigned_int)
	    == new_dof_indices.end(),
	    ExcInternalError());
  }

} // namespace DoFRenumbering



/*-------------- Explicit Instantiations -------------------------------*/
#include "dof_renumbering.inst"


DEAL_II_NAMESPACE_CLOSE
