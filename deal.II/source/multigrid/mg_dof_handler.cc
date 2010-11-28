//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <base/memory_consumption.h>
#include <lac/sparse_matrix.h>
#include <lac/constraint_matrix.h>
#include <dofs/dof_levels.h>
#include <dofs/dof_faces.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>
#include <grid/tria_levels.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria.h>
#include <base/geometry_info.h>
#include <fe/fe.h>
#include <base/exceptions.h>

#include <algorithm>

DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace distributed
  {
    template <int, int> class Triangulation;
  }
}


//TODO[WB]: this class is currently only implemented for dim==spacedim. to
//make the general case happen it should undergo a similar transformation as
//the ::DoFHandler class, i.e. much of the functions here that are dimension
//dependent should be moved into a local Implementation class.
// Some of the following already goes in the right direction
namespace internal
{
  namespace MGDoFHandler
  {
				     // access class
				     // dealii::MGDoFHandler instead
				     // of namespace
				     // internal::MGDoFHandler. same
				     // for dealii::DoFHandler
    using dealii::MGDoFHandler;
    using dealii::DoFHandler;

/**
 * A class with the same purpose as the similarly named class of the
 * Triangulation class. See there for more information.
 */
    struct Implementation
    {

					 /**
					  * Distribute dofs on the given
					  * cell, with new dofs starting
					  * with index
					  * @p next_free_dof. Return the
					  * next unused index number.
					  *
					  * This function is excluded from
					  * the @p distribute_dofs
					  * function since it can not be
					  * implemented dimension
					  * independent.
					  *
					  * Note that unlike for the usual
					  * dofs, here all cells and not
					  * only active ones are allowed.
					  */
	template <int spacedim>
	static
	unsigned int
	distribute_dofs_on_cell (typename MGDoFHandler<1,spacedim>::cell_iterator &cell,
				 unsigned int   next_free_dof)
	  {
	    const unsigned int dim = 1;

					     // distribute dofs of vertices
	    if (cell->get_fe().dofs_per_vertex > 0)
	      for (unsigned int v=0; v<GeometryInfo<1>::vertices_per_cell; ++v)
		{
		  typename MGDoFHandler<dim,spacedim>::cell_iterator neighbor = cell->neighbor(v);

		  if (neighbor.state() == IteratorState::valid)
		    {
						       // has neighbor already been processed?
		      if (neighbor->user_flag_set() &&
			  (neighbor->level() == cell->level()))
							 // copy dofs if the neighbor is on
							 // the same level (only then are
							 // mg dofs the same)
			{
			  if (v==0)
			    for (unsigned int d=0; d<cell->get_fe().dofs_per_vertex; ++d)
			      cell->set_mg_vertex_dof_index (cell->level(), 0, d,
							     neighbor->mg_vertex_dof_index (cell->level(), 1, d));
			  else
			    for (unsigned int d=0; d<cell->get_fe().dofs_per_vertex; ++d)
			      cell->set_mg_vertex_dof_index (cell->level(), 1, d,
							     neighbor->mg_vertex_dof_index (cell->level(), 0, d));

							   // next neighbor
			  continue;
			};
		    };

						   // otherwise: create dofs newly
		  for (unsigned int d=0; d<cell->get_fe().dofs_per_vertex; ++d)
		    cell->set_mg_vertex_dof_index (cell->level(), v, d, next_free_dof++);
		};

					     // dofs of line
	    if (cell->get_fe().dofs_per_line > 0)
	      for (unsigned int d=0; d<cell->get_fe().dofs_per_line; ++d)
		cell->set_mg_dof_index (cell->level(), d, next_free_dof++);

					     // note that this cell has been processed
	    cell->set_user_flag ();

	    return next_free_dof;
	  }


	template <int spacedim>
	static
	unsigned int
	distribute_dofs_on_cell (typename MGDoFHandler<2,spacedim>::cell_iterator &cell,
				 unsigned int   next_free_dof)
	  {
	    const unsigned int dim = 2;
	  if (cell->get_fe().dofs_per_vertex > 0)
					     // number dofs on vertices
	    for (unsigned int vertex=0; vertex<GeometryInfo<2>::vertices_per_cell; ++vertex)
					       // check whether dofs for this
					       // vertex have been distributed
					       // (only check the first dof)
	      if (cell->mg_vertex_dof_index(cell->level(), vertex, 0) == DoFHandler<2>::invalid_dof_index)
		for (unsigned int d=0; d<cell->get_fe().dofs_per_vertex; ++d)
		  cell->set_mg_vertex_dof_index (cell->level(), vertex, d, next_free_dof++);

					   // for the four sides
	  if (cell->get_fe().dofs_per_line > 0)
	    for (unsigned int side=0; side<GeometryInfo<2>::faces_per_cell; ++side)
	      {
		typename MGDoFHandler<dim,spacedim>::line_iterator line = cell->line(side);

						 // distribute dofs if necessary:
						 // check whether line dof is already
						 // numbered (check only first dof)
		if (line->mg_dof_index(cell->level(), 0) == DoFHandler<2>::invalid_dof_index)
						   // if not: distribute dofs
		  for (unsigned int d=0; d<cell->get_fe().dofs_per_line; ++d)
		    line->set_mg_dof_index (cell->level(), d, next_free_dof++);
	      };


					   // dofs of quad
	  if (cell->get_fe().dofs_per_quad > 0)
	    for (unsigned int d=0; d<cell->get_fe().dofs_per_quad; ++d)
	      cell->set_mg_dof_index (cell->level(), d, next_free_dof++);


					   // note that this cell has been processed
	  cell->set_user_flag ();

	  return next_free_dof;
	}


	template <int spacedim>
	static
	unsigned int
	distribute_dofs_on_cell (typename MGDoFHandler<3,spacedim>::cell_iterator &cell,
				 unsigned int   next_free_dof)
	  {
	    const unsigned int dim = 3;
	  if (cell->get_fe().dofs_per_vertex > 0)
					     // number dofs on vertices
	    for (unsigned int vertex=0; vertex<GeometryInfo<3>::vertices_per_cell; ++vertex)
					       // check whether dofs for this
					       // vertex have been distributed
					       // (only check the first dof)
	      if (cell->mg_vertex_dof_index(cell->level(), vertex, 0) == DoFHandler<3>::invalid_dof_index)
		for (unsigned int d=0; d<cell->get_fe().dofs_per_vertex; ++d)
		  cell->set_mg_vertex_dof_index (cell->level(), vertex, d, next_free_dof++);

					   // for the lines
	  if (cell->get_fe().dofs_per_line > 0)
	    for (unsigned int l=0; l<GeometryInfo<3>::lines_per_cell; ++l)
	      {
		typename MGDoFHandler<dim,spacedim>::line_iterator line = cell->line(l);

						 // distribute dofs if necessary:
						 // check whether line dof is already
						 // numbered (check only first dof)
		if (line->mg_dof_index(cell->level(), 0) == DoFHandler<3>::invalid_dof_index)
						   // if not: distribute dofs
		  for (unsigned int d=0; d<cell->get_fe().dofs_per_line; ++d)
		    line->set_mg_dof_index (cell->level(), d, next_free_dof++);
	      };

					   // for the quads
	  if (cell->get_fe().dofs_per_quad > 0)
	    for (unsigned int q=0; q<GeometryInfo<3>::quads_per_cell; ++q)
	      {
		typename MGDoFHandler<dim,spacedim>::quad_iterator quad = cell->quad(q);

						 // distribute dofs if necessary:
						 // check whether line dof is already
						 // numbered (check only first dof)
		if (quad->mg_dof_index(cell->level(), 0) == DoFHandler<3>::invalid_dof_index)
						   // if not: distribute dofs
		  for (unsigned int d=0; d<cell->get_fe().dofs_per_quad; ++d)
		    quad->set_mg_dof_index (cell->level(), d, next_free_dof++);
	      };


					   // dofs of cell
	  if (cell->get_fe().dofs_per_hex > 0)
	    for (unsigned int d=0; d<cell->get_fe().dofs_per_hex; ++d)
	      cell->set_mg_dof_index (cell->level(), d, next_free_dof++);


					   // note that this cell has
					   // been processed
	  cell->set_user_flag ();

	  return next_free_dof;
	}



	template <int spacedim>
	static
	unsigned int
	get_dof_index (const MGDoFHandler<1,spacedim> &mg_dof_handler,
		       const internal::DoFHandler::DoFLevel<1> &mg_level,
		       const internal::DoFHandler::DoFFaces<1> &,
		       const unsigned int       obj_index,
		       const unsigned int       fe_index,
		       const unsigned int       local_index,
		       const int2type<1>)
	  {
	    return mg_level.lines.
	      get_dof_index (mg_dof_handler,
			     obj_index,
			     fe_index,
			     local_index);
	  }


	template <int spacedim>
	static
	void
	set_dof_index (const MGDoFHandler<1,spacedim> &mg_dof_handler,
		       internal::DoFHandler::DoFLevel<1> &mg_level,
		       internal::DoFHandler::DoFFaces<1> &,
		       const unsigned int       obj_index,
		       const unsigned int       fe_index,
		       const unsigned int       local_index,
		       const unsigned int       global_index,
		       const int2type<1>)
	  {
	    mg_level.lines.
	      set_dof_index (mg_dof_handler,
			     obj_index,
			     fe_index,
			     local_index,
			     global_index);

	  }


	template <int spacedim>
	static
	unsigned int
	get_dof_index(const MGDoFHandler<2,spacedim> &mg_dof_handler,
		      const internal::DoFHandler::DoFLevel<2> &,
		      const internal::DoFHandler::DoFFaces<2> &mg_faces,
		      const unsigned int       obj_index,
		      const unsigned int       fe_index,
		      const unsigned int       local_index,
		      const int2type<1>)
	  {
	    return mg_faces.lines.
	      get_dof_index (mg_dof_handler,
			     obj_index,
			     fe_index,
			     local_index);
	  }


	template <int spacedim>
	static
	void
	set_dof_index (const MGDoFHandler<2,spacedim> &mg_dof_handler,
		       internal::DoFHandler::DoFLevel<2> &,
		       internal::DoFHandler::DoFFaces<2> &mg_faces,
		       const unsigned int       obj_index,
		       const unsigned int       fe_index,
		       const unsigned int       local_index,
		       const unsigned int       global_index,
		       const int2type<1>)
	  {
	    mg_faces.lines.
	      set_dof_index (mg_dof_handler,
			     obj_index,
			     fe_index,
			     local_index,
			     global_index);

	  }


	template <int spacedim>
	static
	unsigned int
	get_dof_index (const MGDoFHandler<2,spacedim> &mg_dof_handler,
		       const internal::DoFHandler::DoFLevel<2> &mg_level,
		       const internal::DoFHandler::DoFFaces<2> &,
		       const unsigned int       obj_index,
		       const unsigned int       fe_index,
		       const unsigned int       local_index,
		       const int2type<2>)
	  {
	    return mg_level.quads.
	      get_dof_index (mg_dof_handler,
			     obj_index,
			     fe_index,
			     local_index);
	  }


	template <int spacedim>
	static
	void
	set_dof_index (const MGDoFHandler<2,spacedim> &mg_dof_handler,
		       internal::DoFHandler::DoFLevel<2> &mg_level,
		       internal::DoFHandler::DoFFaces<2> &,
		       const unsigned int       obj_index,
		       const unsigned int       fe_index,
		       const unsigned int       local_index,
		       const unsigned int       global_index,
		       const int2type<2>)
	  {
	    mg_level.quads.
	      set_dof_index (mg_dof_handler,
			     obj_index,
			     fe_index,
			     local_index,
			     global_index);

	  }


	template <int spacedim>
	static
	unsigned int
	get_dof_index (const MGDoFHandler<3,spacedim> &mg_dof_handler,
		       const internal::DoFHandler::DoFLevel<3> &,
		       const internal::DoFHandler::DoFFaces<3> &mg_faces,
		       const unsigned int       obj_index,
		       const unsigned int       fe_index,
		       const unsigned int       local_index,
		       const int2type<1>)
	  {
	    return mg_faces.lines.
	      get_dof_index (mg_dof_handler,
			     obj_index,
			     fe_index,
			     local_index);
	  }


	template <int spacedim>
	static
	void
	set_dof_index (const MGDoFHandler<3,spacedim> &mg_dof_handler,
		       internal::DoFHandler::DoFLevel<3> &,
		       internal::DoFHandler::DoFFaces<3> &mg_faces,
		       const unsigned int       obj_index,
		       const unsigned int       fe_index,
		       const unsigned int       local_index,
		       const unsigned int       global_index,
		       const int2type<1>)
	  {
	    mg_faces.lines.
	      set_dof_index (mg_dof_handler,
			     obj_index,
			     fe_index,
			     local_index,
			     global_index);

	  }


	template <int spacedim>
	static
	unsigned int
	get_dof_index (const MGDoFHandler<3,spacedim> &mg_dof_handler,
		       const internal::DoFHandler::DoFLevel<3> &,
		       const internal::DoFHandler::DoFFaces<3> &mg_faces,
		       const unsigned int       obj_index,
		       const unsigned int       fe_index,
		       const unsigned int       local_index,
		       const int2type<2>)
	  {
	    return mg_faces.quads.
	      get_dof_index (mg_dof_handler,
			     obj_index,
			     fe_index,
			     local_index);
	  }


	template <int spacedim>
	static
	unsigned int
	get_dof_index (const MGDoFHandler<3,spacedim> &mg_dof_handler,
		       const internal::DoFHandler::DoFLevel<3> &mg_level,
		       const internal::DoFHandler::DoFFaces<3> &,
		       const unsigned int       obj_index,
		       const unsigned int       fe_index,
		       const unsigned int       local_index,
		       const int2type<3>)
	  {
	    return mg_level.hexes.
	      get_dof_index (mg_dof_handler,
			     obj_index,
			     fe_index,
			     local_index);
	  }


	template <int spacedim>
	static
	void
	set_dof_index (const MGDoFHandler<3,spacedim> &mg_dof_handler,
		       internal::DoFHandler::DoFLevel<3> &,
		       internal::DoFHandler::DoFFaces<3> &mg_faces,
		       const unsigned int       obj_index,
		       const unsigned int       fe_index,
		       const unsigned int       local_index,
		       const unsigned int       global_index,
		       const int2type<2>)
	  {
	    mg_faces.quads.
	      set_dof_index (mg_dof_handler,
			     obj_index,
			     fe_index,
			     local_index,
			     global_index);

	  }


	template <int spacedim>
	static
	void
	set_dof_index (const MGDoFHandler<3,spacedim> &mg_dof_handler,
		       internal::DoFHandler::DoFLevel<3> &mg_level,
		       internal::DoFHandler::DoFFaces<3> &,
		       const unsigned int       obj_index,
		       const unsigned int       fe_index,
		       const unsigned int       local_index,
		       const unsigned int       global_index,
		       const int2type<3>)
	  {
	    mg_level.hexes.
	      set_dof_index (mg_dof_handler,
			     obj_index,
			     fe_index,
			     local_index,
			     global_index);

	  }




	template <int spacedim>
	static
	void reserve_space (MGDoFHandler<1,spacedim> &mg_dof_handler)
	  {
	    const unsigned int dim = 1;

	    Assert (mg_dof_handler.get_tria().n_levels() > 0,
		    ExcMessage("Invalid triangulation"));

					     //////////////////////////
					     // DESTRUCTION
	    mg_dof_handler.clear_space ();

					     ////////////////////////////
					     // CONSTRUCTION

					     // first allocate space for the
					     // lines on each level
	    for (unsigned int i=0; i<mg_dof_handler.get_tria().n_levels(); ++i)
	      {
		mg_dof_handler.mg_levels.push_back (new internal::DoFHandler::DoFLevel<1>);

		mg_dof_handler.mg_levels.back()->lines.dofs = std::vector<unsigned int>(mg_dof_handler.get_tria().n_raw_lines(i) *
											mg_dof_handler.get_fe().dofs_per_line,
											DoFHandler<1>::invalid_dof_index);
	      }

					     // now allocate space for the
					     // vertices. To this end, we need
					     // to construct as many objects as
					     // there are vertices and let them
					     // allocate enough space for their
					     // vertex indices on the levels they
					     // live on. We need therefore to
					     // count to how many levels a cell
					     // belongs to, which we do by looping
					     // over all cells and storing the
					     // maximum and minimum level each
					     // vertex we pass by  belongs to
	    mg_dof_handler.mg_vertex_dofs.resize (mg_dof_handler.get_tria().n_vertices());

	    std::vector<unsigned int> min_level (mg_dof_handler.get_tria().n_vertices(), mg_dof_handler.get_tria().n_levels());
	    std::vector<unsigned int> max_level (mg_dof_handler.get_tria().n_vertices(), 0);

	    typename dealii::Triangulation<dim,spacedim>::cell_iterator cell = mg_dof_handler.get_tria().begin(),
						       endc = mg_dof_handler.get_tria().end();
	    for (; cell!=endc; ++cell)
	      for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell;
		   ++vertex)
		{
		  const unsigned int vertex_index = cell->vertex_index(vertex);
		  if (min_level[vertex_index] > static_cast<unsigned int>(cell->level()))
		    min_level[vertex_index] = cell->level();
		  if (max_level[vertex_index] < static_cast<unsigned int>(cell->level()))
		    max_level[vertex_index] = cell->level();
		}

					     // now allocate the needed space
	    for (unsigned int vertex=0; vertex<mg_dof_handler.get_tria().n_vertices(); ++vertex)
	      if (mg_dof_handler.get_tria().vertex_used(vertex))
		{
		  Assert (min_level[vertex] < mg_dof_handler.get_tria().n_levels(),   ExcInternalError());
		  Assert (max_level[vertex] >= min_level[vertex], ExcInternalError());

		  mg_dof_handler.mg_vertex_dofs[vertex].init (min_level[vertex],
					       max_level[vertex],
					       mg_dof_handler.get_fe().dofs_per_vertex);
		}
	      else
		{
		  Assert (min_level[vertex] == mg_dof_handler.get_tria().n_levels(),   ExcInternalError());
		  Assert (max_level[vertex] == 0, ExcInternalError());

						   // reset to original state
		  mg_dof_handler.mg_vertex_dofs[vertex].init (1, 0, 0);
		}
	  }



	template <int spacedim>
	static
	void reserve_space (MGDoFHandler<2,spacedim> &mg_dof_handler)
	  {
	    const unsigned int dim = 2;
	    Assert (mg_dof_handler.get_tria().n_levels() > 0,
		    ExcMessage("Invalid triangulation"));

					     ////////////////////////////
					     // DESTRUCTION
	    mg_dof_handler.clear_space ();

					     ////////////////////////////
					     // CONSTRUCTION

					     // first allocate space for the
					     // lines and quads on each level
	    for (unsigned int i=0; i<mg_dof_handler.get_tria().n_levels(); ++i)
	      {
		mg_dof_handler.mg_levels.push_back (new internal::DoFHandler::DoFLevel<2>);

		mg_dof_handler.mg_levels.back()->quads.dofs = std::vector<unsigned int> (mg_dof_handler.get_tria().n_raw_quads(i) *
											 mg_dof_handler.get_fe().dofs_per_quad,
											 DoFHandler<2>::invalid_dof_index);
	      }

	    mg_dof_handler.mg_faces = new internal::DoFHandler::DoFFaces<2>;
	    mg_dof_handler.mg_faces->lines.dofs = std::vector<unsigned int> (mg_dof_handler.get_tria().n_raw_lines() *
									     mg_dof_handler.get_fe().dofs_per_line,
									     DoFHandler<2>::invalid_dof_index);


					     // now allocate space for the
					     // vertices. To this end, we need
					     // to construct as many objects as
					     // there are vertices and let them
					     // allocate enough space for their
					     // vertex indices on the levels they
					     // live on. We need therefore to
					     // count to how many levels a cell
					     // belongs to, which we do by looping
					     // over all cells and storing the
					     // maximum and minimum level each
					     // vertex we pass by  belongs to
	    mg_dof_handler.mg_vertex_dofs.resize (mg_dof_handler.get_tria().n_vertices());

					     // initialize these arrays with
					     // invalid values (min>max)
	    std::vector<unsigned int> min_level (mg_dof_handler.get_tria().n_vertices(),
						 mg_dof_handler.get_tria().n_levels());
	    std::vector<unsigned int> max_level (mg_dof_handler.get_tria().n_vertices(), 0);

	    typename dealii::Triangulation<dim,spacedim>::cell_iterator cell = mg_dof_handler.get_tria().begin(),
						       endc = mg_dof_handler.get_tria().end();
	    for (; cell!=endc; ++cell)
	      for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell;
		   ++vertex)
		{
		  const unsigned int vertex_index = cell->vertex_index(vertex);
		  if (min_level[vertex_index] > static_cast<unsigned int>(cell->level()))
		    min_level[vertex_index] = cell->level();
		  if (max_level[vertex_index] < static_cast<unsigned int>(cell->level()))
		    max_level[vertex_index] = cell->level();
		}


					     // now allocate the needed space
	    for (unsigned int vertex=0; vertex<mg_dof_handler.get_tria().n_vertices(); ++vertex)
	      if (mg_dof_handler.get_tria().vertex_used(vertex))
		{
		  Assert (min_level[vertex] < mg_dof_handler.get_tria().n_levels(), ExcInternalError());
		  Assert (max_level[vertex] >= min_level[vertex], ExcInternalError());

		  mg_dof_handler.mg_vertex_dofs[vertex].init (min_level[vertex],
					       max_level[vertex],
					       mg_dof_handler.get_fe().dofs_per_vertex);
		}
	      else
		{
		  Assert (min_level[vertex] == mg_dof_handler.get_tria().n_levels(),   ExcInternalError());
		  Assert (max_level[vertex] == 0, ExcInternalError());

						   // reset to original state
		  mg_dof_handler.mg_vertex_dofs[vertex].init (1, 0, 0);
		}
	  }



	template <int spacedim>
	static
	void reserve_space (MGDoFHandler<3,spacedim> &mg_dof_handler)
	  {
	    const unsigned int dim = 3;

	    Assert (mg_dof_handler.get_tria().n_levels() > 0,
		    ExcMessage("Invalid triangulation"));

					     ////////////////////////////
					     // DESTRUCTION
	    mg_dof_handler.clear_space ();

					     ////////////////////////////
					     // CONSTRUCTION

					     // first allocate space for the
					     // lines and quads on each level
	    for (unsigned int i=0; i<mg_dof_handler.get_tria().n_levels(); ++i)
	      {
		mg_dof_handler.mg_levels.push_back (new internal::DoFHandler::DoFLevel<3>);

		mg_dof_handler.mg_levels.back()->hexes.dofs
		  = std::vector<unsigned int> (mg_dof_handler.get_tria().n_raw_hexs(i) *
					       mg_dof_handler.get_fe().dofs_per_hex,
					       DoFHandler<3>::invalid_dof_index);
	      }
	    mg_dof_handler.mg_faces = new internal::DoFHandler::DoFFaces<3>;
	    mg_dof_handler.mg_faces->lines.dofs = std::vector<unsigned int> (mg_dof_handler.get_tria().n_raw_lines() *
									     mg_dof_handler.get_fe().dofs_per_line,
									     DoFHandler<3>::invalid_dof_index);
	    mg_dof_handler.mg_faces->quads.dofs = std::vector<unsigned int> (mg_dof_handler.get_tria().n_raw_quads() *
									     mg_dof_handler.get_fe().dofs_per_quad,
									     DoFHandler<3>::invalid_dof_index);


					     // now allocate space for the
					     // vertices. To this end, we need
					     // to construct as many objects as
					     // there are vertices and let them
					     // allocate enough space for their
					     // vertex indices on the levels they
					     // live on. We need therefore to
					     // count to how many levels a cell
					     // belongs to, which we do by looping
					     // over all cells and storing the
					     // maximum and minimum level each
					     // vertex we pass by  belongs to
	    mg_dof_handler.mg_vertex_dofs.resize (mg_dof_handler.get_tria().n_vertices());

	    std::vector<unsigned int> min_level (mg_dof_handler.get_tria().n_vertices(), mg_dof_handler.get_tria().n_levels());
	    std::vector<unsigned int> max_level (mg_dof_handler.get_tria().n_vertices(), 0);

	    typename dealii::Triangulation<dim,spacedim>::cell_iterator cell = mg_dof_handler.get_tria().begin(),
						       endc = mg_dof_handler.get_tria().end();
	    for (; cell!=endc; ++cell)
	      for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell;
		   ++vertex)
		{
		  const unsigned int vertex_index = cell->vertex_index(vertex);
		  if (min_level[vertex_index] > static_cast<unsigned int>(cell->level()))
		    min_level[vertex_index] = cell->level();
		  if (max_level[vertex_index] < static_cast<unsigned int>(cell->level()))
		    max_level[vertex_index] = cell->level();
		}


					     // now allocate the needed space
	    for (unsigned int vertex=0; vertex<mg_dof_handler.get_tria().n_vertices(); ++vertex)
	      if (mg_dof_handler.get_tria().vertex_used(vertex))
		{
		  Assert (min_level[vertex] < mg_dof_handler.get_tria().n_levels(), ExcInternalError());
		  Assert (max_level[vertex] >= min_level[vertex], ExcInternalError());

		  mg_dof_handler.mg_vertex_dofs[vertex].init (min_level[vertex],
					       max_level[vertex],
					       mg_dof_handler.get_fe().dofs_per_vertex);
		}
	      else
		{
		  Assert (min_level[vertex] == mg_dof_handler.get_tria().n_levels(), ExcInternalError());
		  Assert (max_level[vertex] == 0, ExcInternalError());

						   // reset to original state
		  mg_dof_handler.mg_vertex_dofs[vertex].init (1, 0, 0);
		}
	  }


    };
  }
}




/* ------------------------ MGVertexDoFs ----------------------------------- */
//TODO: This seems horrible memory fragmentation!

template <int dim, int spacedim>
MGDoFHandler<dim,spacedim>::MGVertexDoFs::MGVertexDoFs ()
		:
		coarsest_level (numbers::invalid_unsigned_int),
		finest_level (0),
		indices (0)
{}


template <int dim, int spacedim>
void MGDoFHandler<dim,spacedim>::MGVertexDoFs::init (const unsigned int cl,
					    const unsigned int fl,
					    const unsigned int dofs_per_vertex)
{
  if (indices != 0)
    {
      delete[] indices;
      indices = 0;
    }

  coarsest_level = cl;
  finest_level   = fl;

				   // in the case where an invalid
				   // entry is requested, just leave
				   // everything else alone
  if (cl > fl)
    return;

  const unsigned int n_levels = finest_level-coarsest_level + 1;

  indices = new unsigned int[n_levels * dofs_per_vertex];
  Assert (indices != 0, ExcNoMemory ());

  for (unsigned int i=0; i<n_levels*dofs_per_vertex; ++i)
    indices[i] = DoFHandler<dim, spacedim>::invalid_dof_index;
}


template <int dim, int spacedim>
MGDoFHandler<dim,spacedim>::MGVertexDoFs::~MGVertexDoFs ()
{
  delete[] indices;
}


template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::MGVertexDoFs &
MGDoFHandler<dim,spacedim>::MGVertexDoFs::operator = (const MGVertexDoFs &)
{
  Assert (false,
	  ExcMessage ("We don't know how many dofs per vertex there are, so there's "
		      "no way we can copy the 'indices' array"));
  return *this;
}



template <int dim, int spacedim>
unsigned int MGDoFHandler<dim,spacedim>::MGVertexDoFs::get_coarsest_level () const
{
  return coarsest_level;
}


template <int dim, int spacedim>
unsigned int MGDoFHandler<dim,spacedim>::MGVertexDoFs::get_finest_level () const
{
  return finest_level;
}


/* ------------------------ MGDoFHandler ------------------------------------- */

template <int dim, int spacedim>
const unsigned int MGDoFHandler<dim,spacedim>::dimension;


template <int dim, int spacedim>
MGDoFHandler<dim,spacedim>::MGDoFHandler ()
		:
		mg_faces (NULL)
{}



template <int dim, int spacedim>
MGDoFHandler<dim,spacedim>::MGDoFHandler (const Triangulation<dim,spacedim> &tria)
		:
		DoFHandler<dim,spacedim> (tria),
		mg_faces (NULL)
{
  Assert ((dynamic_cast<const parallel::distributed::Triangulation< dim, spacedim >*>
	   (&tria)
	   == 0),
	  ExcMessage ("The given triangulation is parallel distributed but "
		      "this class does not currently support this."));
}


template <int dim, int spacedim>
MGDoFHandler<dim,spacedim>::~MGDoFHandler ()
{
  clear ();
}



template <int dim, int spacedim>
unsigned int
MGDoFHandler<dim,spacedim>::memory_consumption() const
{
  unsigned int mem = DoFHandler<dim,spacedim>::memory_consumption();
  for (unsigned int l=0;l<mg_levels.size();++l)
    mem += mg_levels[l]->memory_consumption();

  mem += MemoryConsumption::memory_consumption(*mg_faces);

  for (unsigned int l=0;l<mg_vertex_dofs.size();++l)
    mem += sizeof(MGVertexDoFs)
	   + (1+mg_vertex_dofs[l].get_finest_level()-mg_vertex_dofs[l].get_coarsest_level())
	   * sizeof(unsigned int);
  mem += MemoryConsumption::memory_consumption(mg_used_dofs);
  return mem;
}


/*------------------------ Cell iterator functions ------------------------*/


template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::raw_cell_iterator
MGDoFHandler<dim,spacedim>::begin_raw (const unsigned int level) const
{
  switch (dim)
    {
      case 1:
	    return begin_raw_line (level);
      case 2:
	    return begin_raw_quad (level);
      case 3:
	    return begin_raw_hex (level);
      default:
	    Assert (false, ExcNotImplemented());
	    return raw_cell_iterator();
    }
}



template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::cell_iterator
MGDoFHandler<dim,spacedim>::begin (const unsigned int level) const
{
  switch (dim)
    {
      case 1:
	    return begin_line (level);
      case 2:
	    return begin_quad (level);
      case 3:
	    return begin_hex (level);
      default:
	    Assert (false, ExcImpossibleInDim(dim));
	    return cell_iterator();
    }
}



template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::active_cell_iterator
MGDoFHandler<dim,spacedim>::begin_active (const unsigned int level) const
{
  switch (dim)
    {
      case 1:
	    return begin_active_line (level);
      case 2:
	    return begin_active_quad (level);
      case 3:
	    return begin_active_hex (level);
      default:
	    Assert (false, ExcNotImplemented());
	    return active_cell_iterator();
    }
}



template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::raw_cell_iterator
MGDoFHandler<dim,spacedim>::last_raw () const
{
  switch (dim)
    {
      case 1:
	    return last_raw_line ();
      case 2:
	    return last_raw_quad ();
      case 3:
	    return last_raw_hex ();
      default:
	    Assert (false, ExcNotImplemented());
	    return raw_cell_iterator();
    }
}



template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::raw_cell_iterator
MGDoFHandler<dim,spacedim>::last_raw (const unsigned int level) const
{
  switch (dim)
    {
      case 1:
	    return last_raw_line (level);
      case 2:
	    return last_raw_quad (level);
      case 3:
	    return last_raw_hex (level);
      default:
	    Assert (false, ExcNotImplemented());
	    return raw_cell_iterator();
    }
}



template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::cell_iterator
MGDoFHandler<dim,spacedim>::last () const
{
  switch (dim)
    {
      case 1:
	    return last_line ();
      case 2:
	    return last_quad ();
      case 3:
	    return last_hex ();
      default:
	    Assert (false, ExcNotImplemented());
	    return cell_iterator();
    }
}



template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::cell_iterator
MGDoFHandler<dim,spacedim>::last (const unsigned int level) const
{
  switch (dim)
    {
      case 1:
	    return last_line (level);
      case 2:
	    return last_quad (level);
      case 3:
	    return last_hex (level);
      default:
	    Assert (false, ExcNotImplemented());
	    return cell_iterator();
    }
}



template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::active_cell_iterator
MGDoFHandler<dim,spacedim>::last_active () const
{
  switch (dim)
    {
      case 1:
	    return last_active_line ();
      case 2:
	    return last_active_quad ();
      case 3:
	    return last_active_hex ();
      default:
	    Assert (false, ExcNotImplemented());
	    return active_cell_iterator();
    }
}



template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::active_cell_iterator
MGDoFHandler<dim,spacedim>::last_active (const unsigned int level) const
{
  switch (dim)
    {
      case 1:
	    return last_active_line (level);
      case 2:
	    return last_active_quad (level);
      case 3:
	    return last_active_hex (level);
      default:
	    Assert (false, ExcNotImplemented());
	    return active_cell_iterator();
    }
}


template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::raw_cell_iterator
MGDoFHandler<dim,spacedim>::end () const
{
  switch (dim)
    {
      case 1:
	    return end_line();
      case 2:
	    return end_quad();
      case 3:
	    return end_hex();
      default:
	    Assert (false, ExcImpossibleInDim(dim));
	    return raw_cell_iterator();
    }
}



template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::raw_cell_iterator
MGDoFHandler<dim, spacedim>::end_raw (const unsigned int level) const
{
  return (level == this->get_tria().n_levels()-1 ?
	  end() :
	  begin_raw (level+1));
}


template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::cell_iterator
MGDoFHandler<dim, spacedim>::end (const unsigned int level) const
{
  return (level == this->get_tria().n_levels()-1 ?
	  cell_iterator(end()) :
	  begin (level+1));
}


template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::active_cell_iterator
MGDoFHandler<dim, spacedim>::end_active (const unsigned int level) const
{
  return (level == this->get_tria().n_levels()-1 ?
	  active_cell_iterator(end()) :
	  begin_active (level+1));
}


/*------------------------ Face iterator functions ------------------------*/


template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::raw_face_iterator
MGDoFHandler<dim,spacedim>::begin_raw_face () const
{
  switch (dim)
    {
      case 1:
	    Assert (false, ExcImpossibleInDim(1));
	    return raw_face_iterator();
      case 2:
	    return begin_raw_line ();
      case 3:
	    return begin_raw_quad ();
      default:
	    Assert (false, ExcNotImplemented());
	    return raw_face_iterator ();
    }
}



template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::face_iterator
MGDoFHandler<dim,spacedim>::begin_face () const
{
  switch (dim)
    {
      case 1:
	    Assert (false, ExcImpossibleInDim(1));
	    return raw_face_iterator();
      case 2:
	    return begin_line ();
      case 3:
	    return begin_quad ();
      default:
	    Assert (false, ExcNotImplemented());
	    return face_iterator ();
    }
}



template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::active_face_iterator
MGDoFHandler<dim,spacedim>::begin_active_face () const
{
  switch (dim)
    {
      case 1:
	    Assert (false, ExcImpossibleInDim(1));
	    return raw_face_iterator();
      case 2:
	    return begin_active_line ();
      case 3:
	    return begin_active_quad ();
      default:
	    Assert (false, ExcNotImplemented());
	    return active_face_iterator ();
    }
}



template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::raw_face_iterator
MGDoFHandler<dim, spacedim>::end_raw_face () const
{
  return end_face();
}



template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::raw_face_iterator
MGDoFHandler<dim,spacedim>::end_face () const
{
  switch (dim)
    {
      case 1:
	    Assert (false, ExcImpossibleInDim(1));
	    return raw_face_iterator();
      case 2:
	    return end_line ();
      case 3:
	    return end_quad ();
      default:
	    Assert (false, ExcNotImplemented());
	    return raw_face_iterator ();
    }
}



template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::active_face_iterator
MGDoFHandler<dim, spacedim>::end_active_face () const
{
  return active_face_iterator(end_face());
}





template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::raw_face_iterator
MGDoFHandler<dim,spacedim>::last_raw_face () const
{
  switch (dim)
    {
      case 1:
	    Assert (false, ExcImpossibleInDim(1));
	    return raw_face_iterator();
      case 2:
	    return last_raw_line ();
      case 3:
	    return last_raw_quad ();
      default:
	    Assert (false, ExcNotImplemented());
	    return raw_face_iterator ();
    }
}



template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::face_iterator
MGDoFHandler<dim,spacedim>::last_face () const
{
  switch (dim)
    {
      case 1:
	    Assert (false, ExcImpossibleInDim(1));
	    return raw_face_iterator();
      case 2:
	    return last_line ();
      case 3:
	    return last_quad ();
      default:
	    Assert (false, ExcNotImplemented());
	    return raw_face_iterator ();
    }
}



template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::active_face_iterator
MGDoFHandler<dim,spacedim>::last_active_face () const
{
  switch (dim)
    {
      case 1:
	    Assert (false, ExcImpossibleInDim(1));
	    return raw_face_iterator();
      case 2:
	    return last_active_line ();
      case 3:
	    return last_active_quad ();
      default:
	    Assert (false, ExcNotImplemented());
	    return raw_face_iterator ();
    }
}


/*------------------------ Line iterator functions ------------------------*/



template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::raw_line_iterator
MGDoFHandler<dim, spacedim>::begin_raw_line (const unsigned int level) const
{
  switch (dim)
    {
      case 1:
	    Assert (level<this->get_tria().n_levels(), ExcInvalidLevel(level));

	    if (this->get_tria().n_raw_lines(level) == 0)
	      return end_line ();

	    return raw_line_iterator (&this->get_tria(),
				      level,
				      0,
				      this);

      default:
	    Assert (level == 0, ExcFacesHaveNoLevel());
	    return raw_line_iterator (&this->get_tria(),
				      0,
				      0,
				      this);
    }
}


template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::line_iterator
MGDoFHandler<dim, spacedim>::begin_line (const unsigned int level) const
{
  				   // level is checked in begin_raw
  raw_line_iterator ri = begin_raw_line (level);
  if (ri.state() != IteratorState::valid)
    return ri;
  while (ri->used() == false)
    if ((++ri).state() != IteratorState::valid)
      return ri;
  return ri;
}



template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::active_line_iterator
MGDoFHandler<dim, spacedim>::begin_active_line (const unsigned int level) const
{
  				   // level is checked in begin_raw
  line_iterator i = begin_line (level);
  if (i.state() != IteratorState::valid)
    return i;
  while (i->has_children())
    if ((++i).state() != IteratorState::valid)
      return i;
  return i;
}



template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::raw_line_iterator
MGDoFHandler<dim, spacedim>::end_line () const
{
  return raw_line_iterator (&this->get_tria(),
			    -1,
			    -1,
			    this);
}



template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::raw_line_iterator
MGDoFHandler<dim, spacedim>::last_raw_line (const unsigned int level) const
{
  switch (dim)
    {
      case 1:
	    Assert (level<this->get_tria().n_levels(), ExcInvalidLevel(level));
	    Assert (this->get_tria().n_raw_lines(level) != 0,
		    ExcEmptyLevel (level));

	    return raw_line_iterator (&this->get_tria(),
				      level,
				      this->get_tria().n_raw_lines(level)-1,
				      this);

      default:
	    Assert (level == 0, ExcFacesHaveNoLevel());
	    return raw_line_iterator (&this->get_tria(),
				      0,
				      this->get_tria().n_raw_lines()-1,
				      this);
    }
}



template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::raw_line_iterator
MGDoFHandler<dim, spacedim>::last_raw_line () const
{
  if (dim == 1)
    return last_raw_line (this->get_tria().n_levels()-1);
  else
    return last_raw_line (0);
}


template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::line_iterator
MGDoFHandler<dim, spacedim>::last_line (const unsigned int level) const
{
  				   // level is checked in last_raw
  raw_line_iterator ri = last_raw_line(level);
  if (ri->used()==true)
    return ri;
  while ((--ri).state() == IteratorState::valid)
    if (ri->used()==true)
      return ri;
  return ri;
}


template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::line_iterator
MGDoFHandler<dim, spacedim>::last_line () const
{
  if (dim == 1)
    return last_line (this->get_tria().n_levels()-1);
  else
    return last_line (0);
}


template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::active_line_iterator
MGDoFHandler<dim, spacedim>::last_active_line (const unsigned int level) const
{
				   // level is checked in last_raw
  line_iterator i = last_line(level);
  if (i->has_children()==false)
    return i;
  while ((--i).state() == IteratorState::valid)
    if (i->has_children()==false)
      return i;
  return i;
}


template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::active_line_iterator
MGDoFHandler<dim, spacedim>::last_active_line () const
{
  if (dim == 1)
    return last_active_line (this->get_tria().n_levels()-1);
  else
    return last_active_line (0);
}


template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::raw_line_iterator
MGDoFHandler<dim, spacedim>::end_raw_line (const unsigned int level) const
{
  Assert (dim == 1 || level == 0, ExcFacesHaveNoLevel());
  if (dim == 1)
    return (level == this->get_tria().n_levels()-1 ?
	    end_line() :
	    begin_raw_line (level+1));
  else
    return end_line();
}


template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::line_iterator
MGDoFHandler<dim, spacedim>::end_line (const unsigned int level) const
{
  Assert (dim == 1 || level == 0, ExcFacesHaveNoLevel());
  if (dim == 1)
    return (level == this->get_tria().n_levels()-1 ?
	    line_iterator(end_line()) :
	    begin_line (level+1));
  else
    return line_iterator(end_line());
}


template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::active_line_iterator
MGDoFHandler<dim, spacedim>::end_active_line (const unsigned int level) const
{
  Assert (dim == 1 || level == 0, ExcFacesHaveNoLevel());
  if (dim == 1)
    return (level == this->get_tria().n_levels()-1 ?
	    active_line_iterator(end_line()) :
	    begin_active_line (level+1));
  else
    return active_line_iterator(end_line());
}



/*------------------------ Quad iterator functions ------------------------*/


template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::raw_quad_iterator
MGDoFHandler<dim,spacedim>::begin_raw_quad (const unsigned int level) const
{
  switch (dim)
    {
      case 1:
	    Assert (false, ExcImpossibleInDim(1));
	    return raw_hex_iterator();
      case 2:
      {
	Assert (level<this->get_tria().n_levels(), ExcInvalidLevel(level));

	if (this->get_tria().n_raw_quads(level) == 0)
	  return end_quad();

	return raw_quad_iterator (&this->get_tria(),
				  level,
				  0,
				  this);
      }

      case 3:
      {
	Assert (level == 0, ExcFacesHaveNoLevel());

	return raw_quad_iterator (&this->get_tria(),
				  0,
				  0,
				  this);
      }


      default:
	    Assert (false, ExcNotImplemented());
	    return raw_hex_iterator();
    }
}



template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::quad_iterator
MGDoFHandler<dim,spacedim>::begin_quad (const unsigned int level) const
{
				   // level is checked in begin_raw
  raw_quad_iterator ri = begin_raw_quad (level);
  if (ri.state() != IteratorState::valid)
    return ri;
  while (ri->used() == false)
    if ((++ri).state() != IteratorState::valid)
      return ri;
  return ri;
}



template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::active_quad_iterator
MGDoFHandler<dim,spacedim>::begin_active_quad (const unsigned int level) const
{
				   // level is checked in begin_raw
  quad_iterator i = begin_quad (level);
  if (i.state() != IteratorState::valid)
    return i;
  while (i->has_children())
    if ((++i).state() != IteratorState::valid)
      return i;
  return i;
}



template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::raw_quad_iterator
MGDoFHandler<dim, spacedim>::end_raw_quad (const unsigned int level) const
{
  Assert (dim == 2 || level == 0, ExcFacesHaveNoLevel());
  if (dim == 2)
    return (level == this->get_tria().n_levels()-1 ?
	    end_quad() :
	    begin_raw_quad (level+1));
  else
    return end_quad();
}



template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::quad_iterator
MGDoFHandler<dim, spacedim>::end_quad (const unsigned int level) const
{
  Assert (dim == 2 || level == 0, ExcFacesHaveNoLevel());
  if (dim == 2)
    return (level == this->get_tria().n_levels()-1 ?
	    quad_iterator(end_quad()) :
	    begin_quad (level+1));
  else
    return quad_iterator(end_quad());
}


template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::active_quad_iterator
MGDoFHandler<dim, spacedim>::end_active_quad (const unsigned int level) const
{
  Assert(dim == 2 || level == 0, ExcFacesHaveNoLevel());
  if (dim == 2)
    return (level == this->get_tria().n_levels()-1 ?
	    active_quad_iterator(end_quad()) :
	    begin_active_quad (level+1));
  else
    return active_quad_iterator(end_quad());
}



template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::raw_quad_iterator
MGDoFHandler<dim,spacedim>::end_quad () const
{
  return raw_quad_iterator (&this->get_tria(),
			    -1,
			    -1,
			    this);
}



template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::raw_quad_iterator
MGDoFHandler<dim,spacedim>::last_raw_quad (const unsigned int level) const
{
  switch (dim)
    {
      case 1:
	    Assert (false, ExcImpossibleInDim(1));
	    return raw_quad_iterator();
      case 2:
	    Assert (level<this->get_tria().n_levels(),
		    ExcInvalidLevel(level));
	    Assert (this->get_tria().n_raw_quads(level) != 0,
		    ExcEmptyLevel (level));
	    return raw_quad_iterator (&this->get_tria(),
				      level,
				      this->get_tria().n_raw_quads(level)-1,
				      this);
      case 3:
	    Assert (level == 0, ExcFacesHaveNoLevel());
	    return raw_quad_iterator (&this->get_tria(),
				      0,
				      this->get_tria().n_raw_quads()-1,
				      this);
      default:
	    Assert (false, ExcNotImplemented());
	    return raw_quad_iterator();
    }
}



template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::raw_quad_iterator
MGDoFHandler<dim,spacedim>::last_raw_quad () const
{
  if (dim == 2)
    return last_raw_quad (this->get_tria().n_levels()-1);
  else
    return last_raw_quad (0);
}



template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::quad_iterator
MGDoFHandler<dim,spacedim>::last_quad (const unsigned int level) const
{
				   // level is checked in last_raw
  raw_quad_iterator ri = last_raw_quad(level);
  if (ri->used()==true)
    return ri;
  while ((--ri).state() == IteratorState::valid)
    if (ri->used()==true)
      return ri;
  return ri;
}



template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::quad_iterator
MGDoFHandler<dim,spacedim>::last_quad () const
{
  if (dim == 2)
    return last_quad (this->get_tria().n_levels()-1);
  else
    return last_quad (0);
}



template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::active_quad_iterator
MGDoFHandler<dim,spacedim>::last_active_quad (const unsigned int level) const
{
				   // level is checked in last_raw
  quad_iterator i = last_quad(level);
  if (i->has_children()==false)
    return i;
  while ((--i).state() == IteratorState::valid)
    if (i->has_children()==false)
      return i;
  return i;
}



template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::active_quad_iterator
MGDoFHandler<dim,spacedim>::last_active_quad () const
{
  if (dim == 2)
    return last_active_quad (this->get_tria().n_levels()-1);
  else
    return last_active_quad (0);
}


/*------------------------ Hex iterator functions ------------------------*/


template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::raw_hex_iterator
MGDoFHandler<dim,spacedim>::begin_raw_hex (const unsigned int level) const
{
  switch (dim)
    {
      case 1:
      case 2:
	    Assert (false, ExcImpossibleInDim(1));
	    return raw_hex_iterator();
      case 3:
      {
	Assert (level<this->get_tria().n_levels(), ExcInvalidLevel(level));

	if (this->get_tria().n_raw_hexs(level) == 0)
	  return end_hex();

	return raw_hex_iterator (&this->get_tria(),
				 level,
				 0,
				 this);
      }

      default:
	    Assert (false, ExcNotImplemented());
	    return raw_hex_iterator();
    }
}



template <int dim, int spacedim>
typename MGDoFHandler<dim,spacedim>::hex_iterator
MGDoFHandler<dim,spacedim>::begin_hex (const unsigned int level) const
{
				   // level is checked in begin_raw
  raw_hex_iterator ri = begin_raw_hex (level);
  if (ri.state() != IteratorState::valid)
    return ri;
  while (ri->used() == false)
    if ((++ri).state() != IteratorState::valid)
      return ri;
  return ri;
}



template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::active_hex_iterator
MGDoFHandler<dim, spacedim>::begin_active_hex (const unsigned int level) const
{
				   // level is checked in begin_raw
  hex_iterator i = begin_hex (level);
  if (i.state() != IteratorState::valid)
    return i;
  while (i->has_children())
    if ((++i).state() != IteratorState::valid)
      return i;
  return i;
}



template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::raw_hex_iterator
MGDoFHandler<dim, spacedim>::end_raw_hex (const unsigned int level) const
{
  return (level == this->get_tria().n_levels()-1 ?
	  end_hex() :
	  begin_raw_hex (level+1));
}


template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::hex_iterator
MGDoFHandler<dim, spacedim>::end_hex (const unsigned int level) const
{
  return (level == this->get_tria().n_levels()-1 ?
	  hex_iterator(end_hex()) :
	  begin_hex (level+1));
}


template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::active_hex_iterator
MGDoFHandler<dim, spacedim>::end_active_hex (const unsigned int level) const
{
  return (level == this->get_tria().n_levels()-1 ?
	  active_hex_iterator(end_hex()) :
	  begin_active_hex (level+1));
}



template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::raw_hex_iterator
MGDoFHandler<dim, spacedim>::end_hex () const
{
  return raw_hex_iterator (&this->get_tria(),
			   -1,
			   -1,
			   this);
}



template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::raw_hex_iterator
MGDoFHandler<dim, spacedim>::last_raw_hex (const unsigned int level) const
{
  switch (dim)
    {
      case 1:
      case 2:
	    Assert (false, ExcImpossibleInDim(dim));
	    return raw_hex_iterator();

      case 3:
	    Assert (level<this->get_tria().n_levels(),
		    ExcInvalidLevel(level));
	    Assert (this->get_tria().n_raw_hexs(level) != 0,
		    ExcEmptyLevel (level));

	    return raw_hex_iterator (&this->get_tria(),
				     level,
				     this->get_tria().n_raw_hexs(level)-1,
				     this);
      default:
	    Assert (false, ExcNotImplemented());
	    return raw_hex_iterator();
    }
}



template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::raw_hex_iterator
MGDoFHandler<dim, spacedim>::last_raw_hex () const
{
  return last_raw_hex (this->get_tria().n_levels()-1);
}



template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::hex_iterator
MGDoFHandler<dim, spacedim>::last_hex (const unsigned int level) const
{
				   // level is checked in last_raw
  raw_hex_iterator ri = last_raw_hex(level);
  if (ri->used()==true)
    return ri;
  while ((--ri).state() == IteratorState::valid)
    if (ri->used()==true)
      return ri;
  return ri;
}



template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::hex_iterator
MGDoFHandler<dim, spacedim>::last_hex () const
{
  return last_hex (this->get_tria().n_levels()-1);
}



template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::active_hex_iterator
MGDoFHandler<dim, spacedim>::last_active_hex (const unsigned int level) const
{
				   // level is checked in last_raw
  hex_iterator i = last_hex(level);
  if (i->has_children()==false)
    return i;
  while ((--i).state() == IteratorState::valid)
    if (i->has_children()==false)
      return i;
  return i;
}



template <int dim, int spacedim>
typename MGDoFHandler<dim, spacedim>::active_hex_iterator
MGDoFHandler<dim, spacedim>::last_active_hex () const
{
  return last_active_hex (this->get_tria().n_levels()-1);
}




//---------------------------------------------------------------------------


template <int dim, int spacedim>
template <int structdim>
unsigned int
MGDoFHandler<dim,spacedim>::
get_dof_index (const unsigned int       obj_level,
	       const unsigned int       obj_index,
	       const unsigned int       fe_index,
	       const unsigned int       local_index) const
{
  return internal::MGDoFHandler::Implementation::
    get_dof_index (*this,
		   *this->mg_levels[obj_level],
		   *this->mg_faces,
		   obj_index,
		   fe_index,
		   local_index,
		   internal::int2type<structdim>());
}



template <int dim, int spacedim>
template <int structdim>
void
MGDoFHandler<dim,spacedim>::
set_dof_index (const unsigned int obj_level,
	       const unsigned int obj_index,
	       const unsigned int fe_index,
	       const unsigned int local_index,
	       const unsigned int global_index) const
{
  internal::MGDoFHandler::Implementation::
    set_dof_index (*this,
		   *this->mg_levels[obj_level],
		   *this->mg_faces,
		   obj_index,
		   fe_index,
		   local_index,
		   global_index,
		   internal::int2type<structdim>());
}



template <int dim, int spacedim>
void MGDoFHandler<dim,spacedim>::distribute_dofs (const FiniteElement<dim,spacedim> &fe,
					 const unsigned int        offset)
{  
				   // first distribute global dofs
  DoFHandler<dim,spacedim>::distribute_dofs (fe);


				   // reserve space for the MG dof numbers
  reserve_space ();
  mg_used_dofs.resize (this->get_tria().n_levels(), 0);

				   // Clear user flags because we will
				   // need them. But first we save
				   // them and make sure that we
				   // restore them later such that at
				   // the end of this function the
				   // Triangulation will be in the
				   // same state as it was at the
				   // beginning of this function.
  std::vector<bool> user_flags;
  this->get_tria().save_user_flags(user_flags);
  const_cast<Triangulation<dim,spacedim> &>(this->get_tria()).clear_user_flags ();

				   // now distribute indices on each level
				   // separately
  for (unsigned int level=0; level<this->get_tria().n_levels(); ++level)
    {
      unsigned int next_free_dof = offset;
      cell_iterator cell = begin(level),
		    endc = end(level);

      for (; cell != endc; ++cell)
	next_free_dof = internal::MGDoFHandler::Implementation::distribute_dofs_on_cell<spacedim> (cell, next_free_dof);

      mg_used_dofs[level] = next_free_dof;
    }

				   // finally restore the user flags
  const_cast<Triangulation<dim,spacedim> &>(this->get_tria())
    .load_user_flags(user_flags);

  this->block_info_object.initialize(*this, true);
}



template <int dim, int spacedim>
void
MGDoFHandler<dim,spacedim>::clear ()
{
				   // release own memory
  clear_space ();

				   // let base class release its mem
				   // as well
  DoFHandler<dim,spacedim>::clear ();
}



template <int dim, int spacedim>
unsigned int MGDoFHandler<dim,spacedim>::n_dofs (const unsigned int level) const {
  Assert (level < mg_used_dofs.size(), ExcInvalidLevel(level));

  return mg_used_dofs[level];
}



template <>
void MGDoFHandler<1>::renumber_dofs (const unsigned int level,
				     const std::vector<unsigned int> &new_numbers) {
  Assert (new_numbers.size() == n_dofs(level), DoFHandler<1>::ExcRenumberingIncomplete());

				   // note that we can not use cell iterators
				   // in this function since then we would
				   // renumber the dofs on the interface of
				   // two cells more than once. Anyway, this
				   // ways it's not only more correct but also
				   // faster
  for (std::vector<MGVertexDoFs>::iterator i=mg_vertex_dofs.begin();
       i!=mg_vertex_dofs.end(); ++i)
				     // if the present vertex lives on
				     // the present level
    if ((i->get_coarsest_level() <= level) &&
	(i->get_finest_level() >= level))
      for (unsigned int d=0; d<this->get_fe().dofs_per_vertex; ++d)
	i->set_index (level, d, this->get_fe().dofs_per_vertex,
		      new_numbers[i->get_index (level, d,
						this->get_fe().dofs_per_vertex)]);

  for (std::vector<unsigned int>::iterator i=mg_levels[level]->lines.dofs.begin();
       i!=mg_levels[level]->lines.dofs.end(); ++i)
    {
      if (*i != DoFHandler<1>::invalid_dof_index)
	{
	  Assert(*i<new_numbers.size(), ExcInternalError());
	  *i = new_numbers[*i];
	}
    }
}



template <>
void MGDoFHandler<2>::renumber_dofs (const unsigned int  level,
				     const std::vector<unsigned int>  &new_numbers) {
  Assert (new_numbers.size() == n_dofs(level),
	  DoFHandler<2>::ExcRenumberingIncomplete());

  for (std::vector<MGVertexDoFs>::iterator i=mg_vertex_dofs.begin();
       i!=mg_vertex_dofs.end(); ++i)
				     // if the present vertex lives on
				     // the present level
    if ((i->get_coarsest_level() <= level) &&
	(i->get_finest_level() >= level))
      for (unsigned int d=0; d<this->get_fe().dofs_per_vertex; ++d)
	i->set_index (level, d, this->get_fe().dofs_per_vertex,
		      new_numbers[i->get_index (level, d,
						this->get_fe().dofs_per_vertex)]);

  if (this->get_fe().dofs_per_line > 0)
    {
				       // save user flags as they will be modified
      std::vector<bool> user_flags;
      this->get_tria().save_user_flags(user_flags);
      const_cast<Triangulation<2> &>(this->get_tria()).clear_user_flags ();

				       // flag all lines adjacent to cells of the current
				       // level, as those lines logically belong to the same
				       // level as the cell, at least for for isotropic
				       // refinement
      cell_iterator cell = begin(level),
		endcell  = end(level);
      for ( ; cell != endcell; ++cell)
	for (unsigned int line=0; line < GeometryInfo<2>::faces_per_cell; ++line)
	  cell->face(line)->set_user_flag();

      line_iterator line = begin_line(),
		 endline = end_line();

      for( ; line != endline; ++line)
	if (line->user_flag_set())
	  {
	    for (unsigned int d=0; d<this->get_fe().dofs_per_line; ++d)
	      line->set_mg_dof_index (level, d, new_numbers[line->mg_dof_index(level, d)]);
	    line->clear_user_flag();
	  }
				       // finally, restore user flags
      const_cast<Triangulation<2> &>(this->get_tria()).load_user_flags (user_flags);
    }

  for (std::vector<unsigned int>::iterator i=mg_levels[level]->quads.dofs.begin();
       i!=mg_levels[level]->quads.dofs.end(); ++i)
    {
      if (*i != DoFHandler<2>::invalid_dof_index)
	{
	  Assert(*i<new_numbers.size(), ExcInternalError());
	  *i = new_numbers[*i];
	}
    }
}



template <>
void MGDoFHandler<3>::renumber_dofs (const unsigned int  level,
				     const std::vector<unsigned int>  &new_numbers) {
  Assert (new_numbers.size() == n_dofs(level),
	  DoFHandler<3>::ExcRenumberingIncomplete());

  for (std::vector<MGVertexDoFs>::iterator i=mg_vertex_dofs.begin();
       i!=mg_vertex_dofs.end(); ++i)
				     // if the present vertex lives on
				     // the present level
    if ((i->get_coarsest_level() <= level) &&
	(i->get_finest_level() >= level))
      for (unsigned int d=0; d<this->get_fe().dofs_per_vertex; ++d)
	i->set_index (level, d, this->get_fe().dofs_per_vertex,
		      new_numbers[i->get_index (level, d,
						this->get_fe().dofs_per_vertex)]);

				   // LINE DoFs
  if (this->get_fe().dofs_per_line > 0)
    {
				       // save user flags as they will be modified
      std::vector<bool> user_flags;
      this->get_tria().save_user_flags(user_flags);
      const_cast<Triangulation<3> &>(this->get_tria()).clear_user_flags ();

				       // flag all lines adjacent to cells of the current
				       // level, as those lines logically belong to the same
				       // level as the cell, at least for for isotropic
				       // refinement

      cell_iterator cell = begin(level),
		endcell  = end(level);
      for ( ; cell != endcell ; ++cell)
	for (unsigned int line=0; line < GeometryInfo<3>::lines_per_cell; ++line)
	  cell->line(line)->set_user_flag();


      line_iterator line = begin_line(),
		 endline = end_line();

      for( ; line != endline; ++line)
	if (line->user_flag_set())
	  {
	    for (unsigned int d=0; d<this->get_fe().dofs_per_line; ++d)
	      line->set_mg_dof_index (level, d, new_numbers[line->mg_dof_index(level, d)]);
	    line->clear_user_flag();
	  }
				       // finally, restore user flags
      const_cast<Triangulation<3> &>(this->get_tria()).load_user_flags (user_flags);
    }

				   // QUAD DoFs
  if (this->get_fe().dofs_per_quad > 0)
    {
				       // save user flags as they will be modified
      std::vector<bool> user_flags;
      this->get_tria().save_user_flags(user_flags);
      const_cast<Triangulation<3> &>(this->get_tria()).clear_user_flags ();

				       // flag all quads adjacent to cells of the current
				       // level, as those lines logically belong to the same
				       // level as the cell, at least for for isotropic
				       // refinement
      cell_iterator cell = begin(level),
		endcell  = end(level);
      for ( ; cell != endcell ; ++cell)
	for (unsigned int quad=0; quad < GeometryInfo<3>::faces_per_cell; ++quad)
	  cell->face(quad)->set_user_flag();

      quad_iterator quad = begin_quad(),
		 endline = end_quad();

      for( ; quad != endline; ++quad)
	if (quad->user_flag_set())
	  {
	    for (unsigned int d=0; d<this->get_fe().dofs_per_quad; ++d)
	      quad->set_mg_dof_index (level, d, new_numbers[quad->mg_dof_index(level, d)]);
	    quad->clear_user_flag();
	  }
				       // finally, restore user flags
      const_cast<Triangulation<3> &>(this->get_tria()).load_user_flags (user_flags);
    }

				   //HEX DoFs
  for (std::vector<unsigned int>::iterator i=mg_levels[level]->hexes.dofs.begin();
       i!=mg_levels[level]->hexes.dofs.end(); ++i)
    {
      if (*i != DoFHandler<3>::invalid_dof_index)
	{
	  Assert(*i<new_numbers.size(), ExcInternalError());
	  *i = new_numbers[*i];
	}
    }
}



template <int dim, int spacedim>
void MGDoFHandler<dim,spacedim>::reserve_space ()
{
  internal::MGDoFHandler::Implementation::reserve_space (*this);
}



template <int dim, int spacedim>
void MGDoFHandler<dim,spacedim>::clear_space ()
{
                                   // delete all levels and set them up
                                   // newly, since vectors are
                                   // troublesome if you want to change
                                   // their size
  for (unsigned int i=0; i<mg_levels.size(); ++i)
    delete mg_levels[i];
  mg_levels.clear ();
  delete mg_faces;
  mg_faces=NULL;

				   // also delete vector of vertex
				   // indices
  std::vector<MGVertexDoFs> tmp;
  std::swap (mg_vertex_dofs, tmp);
}



// explicit instantiations
#include "mg_dof_handler.inst"


DEAL_II_NAMESPACE_CLOSE
