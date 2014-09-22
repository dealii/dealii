// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2014 by the deal.II authors
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

#ifndef __deal2__data_out_h
#define __deal2__data_out_h



#include <deal.II/base/config.h>
#include <deal.II/numerics/data_out_dof_data.h>

#include <deal.II/base/std_cxx11/shared_ptr.h>

DEAL_II_NAMESPACE_OPEN

template <int, int> class FEValuesBase;

namespace internal
{
  namespace DataOut
  {
    /**
     * A derived class for use in the DataOut class. This is a class for the
     * AdditionalData kind of data structure discussed in the documentation of
     * the WorkStream context.
     */
    template <int dim, int spacedim>
    struct ParallelData : public ParallelDataBase<dim,spacedim>
    {
      ParallelData (const unsigned int n_datasets,
                    const unsigned int n_subdivisions,
                    const std::vector<unsigned int> &n_postprocessor_outputs,
                    const Mapping<dim,spacedim> &mapping,
                    const std::vector<std_cxx11::shared_ptr<dealii::hp::FECollection<dim,spacedim> > > &finite_elements,
                    const UpdateFlags update_flags,
                    const std::vector<std::vector<unsigned int> > &cell_to_patch_index_map);

      std::vector<Point<spacedim> > patch_evaluation_points;

      const std::vector<std::vector<unsigned int> > *cell_to_patch_index_map;
    };
  }
}



/**
 * This class is the main class to provide output of data described by finite
 * element fields defined on a collection of cells.
 *
 * This class is an actual implementation of the functionality proposed by
 * the DataOut_DoFData class. It offers a function build_patches() that
 * generates the patches to be written in some graphics format from the data
 * stored in the base class. Most of the interface and an example of its
 * use is described in the documentation of the base class.
 *
 * The only thing this class offers is the function build_patches() which
 * loops over all cells of the triangulation stored by the
 * attach_dof_handler() function of the base class (with the exception of
 * cells of parallel::distributed::Triangulation objects that are not owned by
 * the current processor) and converts the data on these to actual patches
 * which are the objects that are later output by the functions of the base
 * classes. You can give a parameter to the function which determines how many
 * subdivisions in each coordinate direction are to be performed, i.e. of how
 * many subcells each patch shall consist. Default is one, but you may want to
 * choose a higher number for higher order elements, for example two for
 * quadratic elements, three for cubic elements three, and so on. The purpose
 * of this parameter is because most graphics programs do not allow to specify
 * higher order polynomial functions in the file formats: only data at
 * vertices can be plotted and is then shown as a bilinear interpolation
 * within the interior of cells. This may be insufficient if you have higher
 * order finite elements, and the only way to achieve better output is to
 * subdivide each cell of the mesh into several cells for graphical output. Of
 * course, what you get to see is still a bilinear interpolation on each cell
 * of the output (where these cells are not subdivisions of the cells of the
 * triangulation in use) due to the same limitations in output formats, but at
 * least a bilinear interpolation of a higher order polynomial on a finer
 * mesh.
 *
 * Note that after having called build_patches() once, you can call one or
 * more of the write() functions of DataOutInterface. You can therefore
 * output the same data in more than one format without having to rebuild
 * the patches.
 *
 *
 * <h3>User interface information</h3>
 *
 * The base classes of this class, DataOutBase, DataOutInterface and
 * DataOut_DoFData() offer several interfaces of their own. Refer to the
 * DataOutBase class's documentation for a discussion of the different
 * output formats presently supported, DataOutInterface for ways of
 * selecting which format to use upon output at run-time and without
 * the need to adapt your program when new formats become available, as
 * well as for flags to determine aspects of output. The DataOut_DoFData()
 * class's documentation has an example of using nodal data to generate
 * output.
 *
 *
 * <h3>Extensions</h3>
 *
 * By default, this class produces patches for all active cells. Sometimes,
 * this is not what you want, maybe because they are simply too many (and too
 * small to be seen individually) or because you only want to see a certain
 * region of the domain (for example in parallel programs such as the step-18
 * example program), or for some other reason.
 *
 * For this, internally build_patches() does not generate
 * the sequence of cells to be converted into patches itself, but relies
 * on the two functions first_cell() and next_cell(). By default, they
 * return the first active cell, and the next active cell, respectively.
 * Since they are @p virtual functions, you can write your own class
 * derived from DataOut in which you overload these two functions to select other
 * cells for output. This may, for example, include only cells that are
 * in parts of a domain (e.g., if you don't care about the solution elsewhere,
 * think for example a buffer region in which you attenuate outgoing waves
 * in the Perfectly Matched Layer method) or if you don't want output to
 * be generated at all levels of an adaptively refined mesh because this
 * creates too much data (in this case, the set of cells returned by
 * your implementations of first_cell() and next_cell() will include
 * non-active cells, and DataOut::build_patches() will simply take
 * interpolated values of the solution instead of the exact values on these
 * cells children for output). Once you derive your own class, you would
 * just create an object of this type instead of an object of type DataOut,
 * and everything else will remain the same.
 *
 * The two functions are not constant, so you may store information within
 * your derived class about the last accessed cell. This is useful if the
 * information of the last cell which was accessed is not sufficient to
 * determine the next one.
 *
 * There is one caveat, however: if you have cell data (in contrast to nodal,
 * or dof, data) such as error indicators, then you must make sure that
 * first_cell() and next_cell() only walk over active cells, since cell data
 * cannot be interpolated to a coarser cell. If you do have cell data and use
 * this pair of functions and they return a non-active cell, then an exception
 * will be thrown.
 *
 * @pre This class only makes sense if the first template
 * argument, <code>dim</code> equals the dimension of the
 * DoFHandler type given as the second template argument, i.e., if
 * <code>dim == DH::dimension</code>. This redundancy is a historical
 * relic from the time where the library had only a single DoFHandler
 * class and this class consequently only a single template argument.
 *
 * @ingroup output
 * @author Wolfgang Bangerth, 1999
 */
template <int dim, class DH=DoFHandler<dim> >
class DataOut : public DataOut_DoFData<DH, DH::dimension, DH::space_dimension>
{
public:
  /**
   * Typedef to the iterator type of the dof handler class under
   * consideration.
   */
  typedef typename DataOut_DoFData<DH,DH::dimension,DH::space_dimension>::cell_iterator cell_iterator;
  typedef typename DataOut_DoFData<DH,DH::dimension,DH::space_dimension>::active_cell_iterator active_cell_iterator;

  /**
   * Enumeration describing the region of the domain in which curved cells
   * shall be created.
   */
  enum CurvedCellRegion
  {
    no_curved_cells,
    curved_boundary,
    curved_inner_cells
  };

  /**
   * This is the central function of this class since it builds the list of
   * patches to be written by the low-level functions of the base class. See
   * the general documentation of this class for further information.
   *
   * The default value <tt>0</tt> of <tt>n_subdivisions</tt> indicates that
   * the value stored in DataOutInterface::default_subdivisions is to be used.
   */
  virtual void build_patches (const unsigned int n_subdivisions = 0);

  /**
   * Same as above, except that the additional first parameter defines a
   * mapping that is to be used in the generation of output. If
   * <tt>n_subdivisions>1</tt>, the points interior of subdivided patches
   * which originate from cells at the boundary of the domain can be computed
   * using the mapping, i.e. a higher order mapping leads to a representation
   * of a curved boundary by using more subdivisions. Some mappings like
   * MappingQEulerian result in curved cells in the interior of the
   * domain. However, there is nor easy way to get this information from the
   * Mapping. Thus the last argument @p curved_region take one of three values
   * resulting in no curved cells at all, curved cells at the boundary
   * (default) or curved cells in the whole domain.
   *
   * Even for non-curved cells the mapping argument can be used for the
   * Eulerian mappings (see class MappingQ1Eulerian) where a mapping is used
   * not only to determine the position of points interior to a cell, but also
   * of the vertices.  It offers an opportunity to watch the solution on a
   * deformed triangulation on which the computation was actually carried out,
   * even if the mesh is internally stored in its undeformed configuration and
   * the deformation is only tracked by an additional vector that holds the
   * deformation of each vertex.
   *
   * @todo The @p mapping argument should be replaced by a
   * hp::MappingCollection in case of a hp::DoFHandler.
   */
  virtual void build_patches (const Mapping<DH::dimension,DH::space_dimension> &mapping,
                              const unsigned int n_subdivisions = 0,
                              const CurvedCellRegion curved_region = curved_boundary);

  /**
   * Return the first cell which we want output for. The default
   * implementation returns the first active cell, but you might want to
   * return other cells in a derived class.
   */
  virtual cell_iterator first_cell ();

  /**
   * Return the next cell after @p cell which we want output for.  If there
   * are no more cells, <tt>#dofs->end()</tt> shall be returned.
   *
   * The default implementation returns the next active cell, but you might
   * want to return other cells in a derived class. Note that the default
   * implementation assumes that the given @p cell is active, which is
   * guaranteed as long as first_cell() is also used from the default
   * implementation. Overloading only one of the two functions might not be a
   * good idea.
   */
  virtual cell_iterator next_cell (const cell_iterator &cell);

  /**
   * Exception
   */
  DeclException1 (ExcInvalidNumberOfSubdivisions,
                  int,
                  << "The number of subdivisions per patch, " << arg1
                  << ", is not valid.");

private:

  /**
   * Return the first cell produced by the first_cell()/next_cell() function
   * pair that is locally owned. If this object operates on a non-distributed
   * triangulation, the result equals what first_cell() returns.
   */
  cell_iterator first_locally_owned_cell ();

  /**
   * Return the next cell produced by the next_cell() function that is locally
   * owned. If this object operates on a non-distributed triangulation, the
   * result equals what first_cell() returns.
   */
  cell_iterator next_locally_owned_cell (const cell_iterator &cell);

  /**
   * Build one patch. This function is called in a WorkStream context.
   *
   * The result is written into the patch variable.
   */
  void build_one_patch (const std::pair<cell_iterator, unsigned int> *cell_and_index,
                        internal::DataOut::ParallelData<DH::dimension, DH::space_dimension> &data,
                        ::dealii::DataOutBase::Patch<DH::dimension, DH::space_dimension> &patch,
                        const CurvedCellRegion curved_cell_region,
                        std::vector<dealii::DataOutBase::Patch<DH::dimension, DH::space_dimension> > &patches);
};



DEAL_II_NAMESPACE_CLOSE

#endif
