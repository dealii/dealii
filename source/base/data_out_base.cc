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


//TODO: Do neighbors for dx and povray smooth triangles

//////////////////////////////////////////////////////////////////////
// Remarks on the implementations
//
// Variable names: in most functions, variable names have been
// standardized in the following way:
//
// n1, n2, ni Number of points in coordinate direction 1, 2, i
//    will be 1 if i>=dim
//
// i1, i2, ii Loop variable running up to ni
//
// d1, d2, di Multiplicators for ii to find positions in the
//    array of nodes.
//////////////////////////////////////////////////////////////////////

#include <deal.II/base/data_out_base.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/std_cxx11/shared_ptr.h>
#include <deal.II/base/mpi.h>

#include <cstring>
#include <algorithm>
#include <iomanip>
#include <ctime>
#include <cmath>
#include <set>
#include <sstream>
#include <fstream>

// we use uint32_t and uint8_t below, which are declared here:
#include <stdint.h>

#ifdef DEAL_II_WITH_ZLIB
#  include <zlib.h>
#endif

#ifdef DEAL_II_WITH_HDF5
#include <hdf5.h>
#endif

DEAL_II_NAMESPACE_OPEN


// we need the following exception from a global function, so can't declare it
// in the usual way inside a class
namespace
{
  DeclException2 (ExcUnexpectedInput,
                  std::string, std::string,
                  << "Unexpected input: expected line\n  <"
                  << arg1
                  << ">\nbut got\n  <"
                  << arg2 << ">");
}


namespace
{
#ifdef DEAL_II_WITH_ZLIB
  // the functions in this namespace are
  // taken from the libb64 project, see
  // http://sourceforge.net/projects/libb64
  //
  // libb64 has been placed in the public
  // domain
  namespace base64
  {
    typedef enum
    {
      step_A, step_B, step_C
    } base64_encodestep;

    typedef struct
    {
      base64_encodestep step;
      char result;
    } base64_encodestate;

    void base64_init_encodestate(base64_encodestate *state_in)
    {
      state_in->step = step_A;
      state_in->result = 0;
    }

    inline
    char base64_encode_value(char value_in)
    {
      static const char *encoding
        = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
      if (value_in > 63) return '=';
      return encoding[(int)value_in];
    }

    int base64_encode_block(const char *plaintext_in,
                            int length_in,
                            char *code_out,
                            base64_encodestate *state_in)
    {
      const char *plainchar = plaintext_in;
      const char *const plaintextend = plaintext_in + length_in;
      char *codechar = code_out;
      char result;
      char fragment;

      result = state_in->result;

      switch (state_in->step)
        {
          while (1)
            {
            case step_A:
              if (plainchar == plaintextend)
                {
                  state_in->result = result;
                  state_in->step = step_A;
                  return codechar - code_out;
                }
              fragment = *plainchar++;
              result = (fragment & 0x0fc) >> 2;
              *codechar++ = base64_encode_value(result);
              result = (fragment & 0x003) << 4;
            case step_B:
              if (plainchar == plaintextend)
                {
                  state_in->result = result;
                  state_in->step = step_B;
                  return codechar - code_out;
                }
              fragment = *plainchar++;
              result |= (fragment & 0x0f0) >> 4;
              *codechar++ = base64_encode_value(result);
              result = (fragment & 0x00f) << 2;
            case step_C:
              if (plainchar == plaintextend)
                {
                  state_in->result = result;
                  state_in->step = step_C;
                  return codechar - code_out;
                }
              fragment = *plainchar++;
              result |= (fragment & 0x0c0) >> 6;
              *codechar++ = base64_encode_value(result);
              result  = (fragment & 0x03f) >> 0;
              *codechar++ = base64_encode_value(result);
            }
        }
      /* control should not reach here */
      return codechar - code_out;
    }

    int base64_encode_blockend(char *code_out, base64_encodestate *state_in)
    {
      char *codechar = code_out;

      switch (state_in->step)
        {
        case step_B:
          *codechar++ = base64_encode_value(state_in->result);
          *codechar++ = '=';
          *codechar++ = '=';
          break;
        case step_C:
          *codechar++ = base64_encode_value(state_in->result);
          *codechar++ = '=';
          break;
        case step_A:
          break;
        }
      *codechar++ = '\0';

      return codechar - code_out;
    }
  }


  /**
   * Do a base64 encoding of the given data.
   *
   * The function allocates memory as
   * necessary and returns a pointer to
   * it. The calling function must release
   * this memory again.
   */
  char *
  encode_block (const char *data,
                const int   data_size)
  {
    base64::base64_encodestate state;
    base64::base64_init_encodestate(&state);

    char *encoded_data = new char[2*data_size+1];

    const int encoded_length_data
      = base64::base64_encode_block (data, data_size,
                                     encoded_data, &state);
    base64::base64_encode_blockend (encoded_data + encoded_length_data,
                                    &state);

    return encoded_data;
  }
#endif


#ifdef DEAL_II_WITH_ZLIB
  /**
   * Do a zlib compression followed
   * by a base64 encoding of the
   * given data. The result is then
   * written to the given stream.
   */
  template <typename T>
  void write_compressed_block (const std::vector<T> &data,
                               std::ostream         &output_stream)
  {
    if (data.size() != 0)
      {
        // allocate a buffer for compressing
        // data and do so
        uLongf compressed_data_length
          = compressBound (data.size() * sizeof(T));
        char *compressed_data = new char[compressed_data_length];
        int err = compress2 ((Bytef *) compressed_data,
                             &compressed_data_length,
                             (const Bytef *) &data[0],
                             data.size() * sizeof(T),
                             Z_BEST_COMPRESSION);
        Assert (err == Z_OK, ExcInternalError());

        // now encode the compression header
        const uint32_t compression_header[4]
          = { 1,                                   /* number of blocks */
              (uint32_t)(data.size() * sizeof(T)), /* size of block */
              (uint32_t)(data.size() * sizeof(T)), /* size of last block */
              (uint32_t)compressed_data_length
            }; /* list of compressed sizes of blocks */

        char *encoded_header = encode_block ((char *)&compression_header[0],
                                             4 * sizeof(compression_header[0]));
        output_stream << encoded_header;
        delete[] encoded_header;

        // next do the compressed
        // data encoding in base64
        char *encoded_data = encode_block (compressed_data,
                                           compressed_data_length);
        delete[] compressed_data;

        output_stream << encoded_data;
        delete[] encoded_data;
      }
  }
#endif
}


// some declarations of functions and locally used classes
namespace DataOutBase
{
  namespace
  {
    /**
     * Class holding the data of one cell of a patch in two space
     * dimensions for output. It is the projection of a cell in
     * three-dimensional space (two coordinates, one height value) to
     * the direction of sight.
     */
    class SvgCell
    {
    public:

      // Center of the cell (three-dimensional)
      Point<3> center;

      /**
       * Vector of vertices of this cell (three-dimensional)
       */
      Point<3> vertices[4];

      /**
       * Depth into the picture, which is defined as the distance from
       * an observer at an the origin in direction of the line of sight.
       */
      float depth;

      /**
       * Vector of vertices of this cell (projected, two-dimensional).
       */
      Point<2> projected_vertices[4];

      // Center of the cell (projected, two-dimensional)
      Point<2> projected_center;

      /**
       * Comparison operator for sorting.
       */
      bool operator < (const SvgCell &) const;
    };

    bool SvgCell::operator < (const SvgCell &e) const
    {
      // note the "wrong" order in
      // which we sort the elements
      return depth > e.depth;
    }



    /**
     * Class holding the data of one cell of a patch in two space
     * dimensions for output. It is the projection of a cell in
     * three-dimensional space (two coordinates, one height value) to
     * the direction of sight.
     */
    class EpsCell2d
    {
    public:

      /**
       * Vector of vertices of this cell.
       */
      Point<2> vertices[4];

      /**
       * Data value from which the actual colors will be computed by the
       * colorization function stated in the <tt>EpsFlags</tt> class.
       */
      float color_value;

      /**
       * Depth into the picture, which is defined as the distance from
       * an observer at an the origin in direction of the line of sight.
       */
      float depth;

      /**
       * Comparison operator for sorting.
       */
      bool operator < (const EpsCell2d &) const;
    };


    /**
     * This is a helper function for the write_gmv() function. There,
     * the data in the patches needs to be copied around as output is
     * one variable globally at a time, rather than all data on each
     * vertex at a time. This copying around can be done detached from
     * the main thread, and is thus moved into this separate function.
     *
     * Note that because of the similarity of the formats, this function
     * is also used by the Vtk and Tecplot output functions.
     */
    template <int dim, int spacedim>
    void
    write_gmv_reorder_data_vectors (const std::vector<Patch<dim,spacedim> > &patches,
                                    Table<2,double>                         &data_vectors)
    {
      // unlike in the main function, we
      // don't have here the data_names
      // field, so we initialize it with
      // the number of data sets in the
      // first patch. the equivalence of
      // these two definitions is checked
      // in the main function.

      // we have to take care, however, whether the
      // points are appended to the end of the
      // patch->data table
      const unsigned int n_data_sets
        =patches[0].points_are_available ? (patches[0].data.n_rows() - spacedim) : patches[0].data.n_rows();

      Assert (data_vectors.size()[0] == n_data_sets,
              ExcInternalError());

      // loop over all patches
      unsigned int next_value = 0;
      for (typename std::vector<Patch<dim,spacedim> >::const_iterator patch=patches.begin();
           patch != patches.end(); ++patch)
        {
          const unsigned int n_subdivisions = patch->n_subdivisions;

          Assert ((patch->data.n_rows() == n_data_sets && !patch->points_are_available) ||
                  (patch->data.n_rows() == n_data_sets+spacedim && patch->points_are_available),
                  ExcDimensionMismatch (patch->points_are_available
                                        ?
                                        (n_data_sets + spacedim)
                                        :
                                        n_data_sets,
                                        patch->data.n_rows()));
          Assert ((n_data_sets == 0)
                  ||
                  (patch->data.n_cols() == Utilities::fixed_power<dim>(n_subdivisions+1)),
                  ExcInvalidDatasetSize (patch->data.n_cols(), n_subdivisions+1));

          for (unsigned int i=0; i<patch->data.n_cols(); ++i, ++next_value)
            for (unsigned int data_set=0; data_set<n_data_sets; ++data_set)
              data_vectors[data_set][next_value] = patch->data(data_set,i);
        }

      for (unsigned int data_set=0; data_set<n_data_sets; ++data_set)
        Assert (data_vectors[data_set].size() == next_value,
                ExcInternalError());
    }
  }
}

//----------------------------------------------------------------------//
// DataOutFilter class member functions
//----------------------------------------------------------------------//

template<int dim>
void DataOutBase::DataOutFilter::write_point(const unsigned int &index, const Point<dim> &p)
{
  Map3DPoint::const_iterator  it;
  unsigned int      internal_ind;
  Point<3>      int_pt;

  for (int d=0; d<3; ++d) int_pt(d) = (d < dim ? p(d) : 0);
  node_dim = dim;
  it = existing_points.find(int_pt);

  // If the point isn't in the set, or we're not filtering duplicate points, add it
  if (it == existing_points.end() || !flags.filter_duplicate_vertices)
    {
      internal_ind = existing_points.size();
      existing_points.insert(std::make_pair(int_pt, internal_ind));
    }
  else
    {
      internal_ind = it->second;
    }
  // Now add the index to the list of filtered points
  filtered_points[index] = internal_ind;
}

void DataOutBase::DataOutFilter::internal_add_cell(const unsigned int &cell_index, const unsigned int &pt_index)
{
  filtered_cells[cell_index] = filtered_points[pt_index];
}

void DataOutBase::DataOutFilter::fill_node_data(std::vector<double> &node_data) const
{
  Map3DPoint::const_iterator  it;

  node_data.resize(existing_points.size()*node_dim);

  for (it=existing_points.begin(); it!=existing_points.end(); ++it)
    {
      for (int d=0; d<node_dim; ++d) node_data[node_dim*it->second+d] = it->first(d);
    }
}

void DataOutBase::DataOutFilter::fill_cell_data(const unsigned int &local_node_offset, std::vector<unsigned int> &cell_data) const
{
  std::map<unsigned int, unsigned int>::const_iterator  it;

  cell_data.resize(filtered_cells.size());

  for (it=filtered_cells.begin(); it!=filtered_cells.end(); ++it)
    {
      cell_data[it->first] = it->second+local_node_offset;
    }
}

template<int dim>
void
DataOutBase::DataOutFilter::write_cell(
  unsigned int index,
  unsigned int start,
  unsigned int d1,
  unsigned int d2,
  unsigned int d3)
{
  unsigned int base_entry = index * GeometryInfo<dim>::vertices_per_cell;
  n_cell_verts = GeometryInfo<dim>::vertices_per_cell;
  internal_add_cell(base_entry+0, start);
  internal_add_cell(base_entry+1, start+d1);
  if (dim>=2)
    {
      internal_add_cell(base_entry+2, start+d2+d1);
      internal_add_cell(base_entry+3, start+d2);
      if (dim>=3)
        {
          internal_add_cell(base_entry+4, start+d3);
          internal_add_cell(base_entry+5, start+d3+d1);
          internal_add_cell(base_entry+6, start+d3+d2+d1);
          internal_add_cell(base_entry+7, start+d3+d2);
        }
    }
}

void DataOutBase::DataOutFilter::write_data_set(const std::string &name, const unsigned int &dimension, const unsigned int &set_num, const Table<2,double> &data_vectors)
{
  unsigned int    num_verts = existing_points.size();
  unsigned int    i, r, d, new_dim;

  // HDF5/XDMF output only supports 1D or 3D output, so force rearrangement if needed
  if (flags.xdmf_hdf5_output && dimension != 1) new_dim = 3;
  else new_dim = dimension;

  // Record the data set name, dimension, and allocate space for it
  data_set_names.push_back(name);
  data_set_dims.push_back(new_dim);
  data_sets.push_back(std::vector<double>(new_dim*num_verts));

  // TODO: averaging, min/max, etc for merged vertices
  for (i=0; i<filtered_points.size(); ++i)
    {
      for (d=0; d<new_dim; ++d)
        {
          r = filtered_points[i];
          if (d < dimension) data_sets.back()[r*new_dim+d] = data_vectors(set_num+d, i);
          else data_sets.back()[r*new_dim+d] = 0;
        }
    }
}


//----------------------------------------------------------------------//
//Auxiliary data
//----------------------------------------------------------------------//

namespace
{
  const char *gmv_cell_type[4] =
  {
    "", "line 2", "quad 4", "hex 8"
  };

  const char *ucd_cell_type[4] =
  {
    "", "line", "quad", "hex"
  };

  const char *tecplot_cell_type[4] =
  {
    "", "lineseg", "quadrilateral", "brick"
  };

#ifdef DEAL_II_HAVE_TECPLOT
  const unsigned int tecplot_binary_cell_type[4] =
  {
    0, 0, 1, 3
  };
#endif

  // NOTE: The dimension of the array is choosen to 5 to allow the choice
  // DataOutBase<deal_II_dimension,deal_II_dimension+1> in general
  // Wolfgang supposed that we don't need it in general, but however this
  // choice avoids a -Warray-bounds check warning
  const unsigned int vtk_cell_type[5] =
  {
    0, 3, 9, 12, static_cast<unsigned int>(-1)
  };

//----------------------------------------------------------------------//
//Auxiliary functions
//----------------------------------------------------------------------//
//For a given patch, compute the node interpolating the corner nodes
//linearly at the point (xstep, ystep, zstep)*1./n_subdivisions.
//If the points are saved in the patch->data member, return the
//saved point instead

//TODO: Make this function return its value, rather than using a reference
// as first argument; take a reference for 'patch', not a pointer
  template <int dim, int spacedim>
  inline
  void
  compute_node(
    Point<spacedim> &node,
    const DataOutBase::Patch<dim,spacedim> *patch,
    const unsigned int xstep,
    const unsigned int ystep,
    const unsigned int zstep,
    const unsigned int n_subdivisions)
  {
    if (patch->points_are_available)
      {
        unsigned int point_no=0;
        // note: switch without break !
        switch (dim)
          {
          case 3:
            Assert (zstep<n_subdivisions+1, ExcIndexRange(zstep,0,n_subdivisions+1));
            point_no+=(n_subdivisions+1)*(n_subdivisions+1)*zstep;
          case 2:
            Assert (ystep<n_subdivisions+1, ExcIndexRange(ystep,0,n_subdivisions+1));
            point_no+=(n_subdivisions+1)*ystep;
          case 1:
            Assert (xstep<n_subdivisions+1, ExcIndexRange(xstep,0,n_subdivisions+1));
            point_no+=xstep;

            // break here for dim<=3
            break;

          default:
            Assert (false, ExcNotImplemented());
          }
        for (unsigned int d=0; d<spacedim; ++d)
          node[d]=patch->data(patch->data.size(0)-spacedim+d,point_no);
      }
    else
      {
        // perform a dim-linear interpolation
        const double stepsize=1./n_subdivisions,
                     xfrac=xstep*stepsize;

        node = (patch->vertices[1] * xfrac) + (patch->vertices[0] * (1-xfrac));
        if (dim>1)
          {
            const double yfrac=ystep*stepsize;
            node*= 1-yfrac;
            node += ((patch->vertices[3] * xfrac) + (patch->vertices[2] * (1-xfrac))) * yfrac;
            if (dim>2)
              {
                const double zfrac=zstep*stepsize;
                node *= (1-zfrac);
                node += (((patch->vertices[5] * xfrac) + (patch->vertices[4] * (1-xfrac)))
                         * (1-yfrac) +
                         ((patch->vertices[7] * xfrac) + (patch->vertices[6] * (1-xfrac)))
                         * yfrac) * zfrac;
              }
          }
      }
  }



  template<int dim, int spacedim>
  static
  void
  compute_sizes(const std::vector<DataOutBase::Patch<dim, spacedim> > &patches,
                unsigned int &n_nodes,
                unsigned int &n_cells)
  {
    n_nodes = 0;
    n_cells = 0;
    for (typename std::vector<DataOutBase::Patch<dim,spacedim> >::const_iterator patch=patches.begin();
         patch!=patches.end(); ++patch)
      {
        n_nodes += Utilities::fixed_power<dim>(patch->n_subdivisions+1);
        n_cells += Utilities::fixed_power<dim>(patch->n_subdivisions);
      }
  }

  /**
   * Class for writing basic
   * entities in @ref
   * SoftwareOpenDX format,
   * depending on the flags.
   */
  class DXStream
  {
  public:
    /**
     * Constructor, storing
     * persistent values for
     * later use.
     */
    DXStream (std::ostream &stream,
              const DataOutBase::DXFlags flags);

    /**
     * Output operator for points.
     */
    template <int dim>
    void write_point (const unsigned int index,
                      const Point<dim> &);

    /**
     * Do whatever is necessary to
     * terminate the list of points.
     */
    void flush_points ();

    /**
     * Write dim-dimensional cell
     * with first vertex at
     * number start and further
     * vertices offset by the
     * specified values. Values
     * not needed are ignored.
     *
     * The order of vertices for
     * these cells in different
     * dimensions is
     * <ol>
     * <li> [0,1]
     * <li> [0,2,1,3]
     * <li> [0,4,2,6,1,5,3,7]
     * </ol>
     */
    template <int dim>
    void write_cell (const unsigned int index,
                     const unsigned int start,
                     const unsigned int x_offset,
                     const unsigned int y_offset,
                     const unsigned int z_offset);

    /**
     * Do whatever is necessary to
     * terminate the list of cells.
     */
    void flush_cells ();

    /**
     * Write a complete set of
     * data for a single node.
     *
     * The index given as first
     * argument indicates the
     * number of a data set, as
     * some output formats require
     * this number to be printed.
     */
    template<typename data>
    void write_dataset (const unsigned int       index,
                        const std::vector<data> &values);

    /**
     * Forwarding of output stream
     */
    template <typename T>
    std::ostream &operator<< (const T &);

  private:
    /**
     * The ostream to use. Since
     * the life span of these
     * objects is small, we use a
     * very simple storage
     * technique.
     */
    std::ostream &stream;

    /**
     * The flags controlling the output
     */
    const DataOutBase::DXFlags flags;
  };

  /**
   * Class for writing basic
   * entities in @ref SoftwareGMV
   * format, depending on the
   * flags.
   */
  class GmvStream
  {
  public:
    /**
     * Constructor, storing
     * persistent values for
     * later use.
     */
    GmvStream (std::ostream &stream,
               const DataOutBase::GmvFlags flags);

    /**
     * Output operator for points.
     */
    template <int dim>
    void write_point (const unsigned int index,
                      const Point<dim> &);

    /**
     * Do whatever is necessary to
     * terminate the list of points.
     */
    void flush_points ();

    /**
     * Write dim-dimensional cell
     * with first vertex at
     * number start and further
     * vertices offset by the
     * specified values. Values
     * not needed are ignored.
     *
     * The order of vertices for
     * these cells in different
     * dimensions is
     * <ol>
     * <li> [0,1]
     * <li> [0,1,3,2]
     * <li> [0,1,3,2,4,5,7,6]
     * </ol>
     */
    template <int dim>
    void write_cell(const unsigned int index,
                    const unsigned int start,
                    const unsigned int x_offset,
                    const unsigned int y_offset,
                    const unsigned int z_offset);

    /**
     * Do whatever is necessary to
     * terminate the list of cells.
     */
    void flush_cells ();

    /**
     * Forwarding of output stream
     */
    template <typename T>
    std::ostream &operator<< (const T &);

    /**
     * Since GMV reads the x, y
     * and z coordinates in
     * separate fields, we enable
     * write() to output only a
     * single selected component
     * at once and do this dim
     * times for the whole set of
     * nodes. This integer can be
     * used to select the
     * component written.
     */
    unsigned int selected_component;

  private:
    /**
     * The ostream to use. Since
     * the life span of these
     * objects is small, we use a
     * very simple storage
     * technique.
     */
    std::ostream &stream;

    /**
     * The flags controlling the output
     */
    const DataOutBase::GmvFlags flags;
  };

  /**
   * Class for writing basic
   * entities in @ref
   * SoftwareTecplot format,
   * depending on the flags.
   */
  class TecplotStream
  {
  public:
    /**
     * Constructor, storing
     * persistent values for
     * later use.
     */
    TecplotStream (std::ostream &stream, const DataOutBase::TecplotFlags flags);

    /**
     * Output operator for points.
     */
    template <int dim>
    void write_point (const unsigned int index,
                      const Point<dim> &);

    /**
     * Do whatever is necessary to
     * terminate the list of points.
     */
    void flush_points ();

    /**
     * Write dim-dimensional cell
     * with first vertex at
     * number start and further
     * vertices offset by the
     * specified values. Values
     * not needed are ignored.
     *
     * The order of vertices for
     * these cells in different
     * dimensions is
     * <ol>
     * <li> [0,1]
     * <li> [0,1,3,2]
     * <li> [0,1,3,2,4,5,7,6]
     * </ol>
     */
    template <int dim>
    void write_cell(const unsigned int index,
                    const unsigned int start,
                    const unsigned int x_offset,
                    const unsigned int y_offset,
                    const unsigned int z_offset);

    /**
     * Do whatever is necessary to
     * terminate the list of cells.
     */
    void flush_cells ();

    /**
     * Forwarding of output stream
     */
    template <typename T>
    std::ostream &operator<< (const T &);

    /**
     * Since TECPLOT reads the x, y
     * and z coordinates in
     * separate fields, we enable
     * write() to output only a
     * single selected component
     * at once and do this dim
     * times for the whole set of
     * nodes. This integer can be
     * used to select the
     * component written.
     */
    unsigned int selected_component;

  private:
    /**
     * The ostream to use. Since
     * the life span of these
     * objects is small, we use a
     * very simple storage
     * technique.
     */
    std::ostream &stream;

    /**
     * The flags controlling the output
     */
    const DataOutBase::TecplotFlags flags;
  };

  /**
   * Class for writing basic
   * entities in UCD format for
   * @ref SoftwareAVS, depending on
   * the flags.
   */
  class UcdStream
  {
  public:
    /**
     * Constructor, storing
     * persistent values for
     * later use.
     */
    UcdStream (std::ostream &stream,
               const DataOutBase::UcdFlags flags);

    /**
     * Output operator for points.
     */
    template <int dim>
    void write_point (const unsigned int index,
                      const Point<dim> &);

    /**
     * Do whatever is necessary to
     * terminate the list of points.
     */
    void flush_points ();

    /**
     * Write dim-dimensional cell
     * with first vertex at
     * number start and further
     * vertices offset by the
     * specified values. Values
     * not needed are ignored.
     *
     * The additional offset 1 is
     * added inside this
     * function.
     *
     * The order of vertices for
     * these cells in different
     * dimensions is
     * <ol>
     * <li> [0,1]
     * <li> [0,1,3,2]
     * <li> [0,1,5,4,2,3,7,6]
     * </ol>
     */
    template <int dim>
    void write_cell(const unsigned int index,
                    const unsigned int start,
                    const unsigned int x_offset,
                    const unsigned int y_offset,
                    const unsigned int z_offset);

    /**
     * Do whatever is necessary to
     * terminate the list of cells.
     */
    void flush_cells ();

    /**
     * Write a complete set of
     * data for a single node.
     *
     * The index given as first
     * argument indicates the
     * number of a data set, as
     * some output formats require
     * this number to be printed.
     */
    template<typename data>
    void write_dataset (const unsigned int       index,
                        const std::vector<data> &values);

    /**
     * Forwarding of output stream
     */
    template <typename T>
    std::ostream &operator<< (const T &);
  private:
    /**
     * The ostream to use. Since
     * the life span of these
     * objects is small, we use a
     * very simple storage
     * technique.
     */
    std::ostream &stream;

    /**
     * The flags controlling the output
     */
    const DataOutBase::UcdFlags flags;
  };

  /**
   * Class for writing basic
   * entities in @ref SoftwareVTK
   * format, depending on the
   * flags.
   */
  class VtkStream
  {
  public:
    /**
     * Constructor, storing
     * persistent values for
     * later use.
     */
    VtkStream (std::ostream &stream,
               const DataOutBase::VtkFlags flags);

    /**
     * Output operator for points.
     */
    template <int dim>
    void write_point (const unsigned int index,
                      const Point<dim> &);

    /**
     * Do whatever is necessary to
     * terminate the list of points.
     */
    void flush_points ();

    /**
     * Write dim-dimensional cell
     * with first vertex at
     * number start and further
     * vertices offset by the
     * specified values. Values
     * not needed are ignored.
     *
     * The order of vertices for
     * these cells in different
     * dimensions is
     * <ol>
     * <li> [0,1]
     * <li> []
     * <li> []
     * </ol>
     */
    template <int dim>
    void write_cell(const unsigned int index,
                    const unsigned int start,
                    const unsigned int x_offset,
                    const unsigned int y_offset,
                    const unsigned int z_offset);

    /**
     * Do whatever is necessary to
     * terminate the list of cells.
     */
    void flush_cells ();

    /**
     * Forwarding of output stream
     */
    template <typename T>
    std::ostream &operator<< (const T &);

  private:
    /**
     * The ostream to use. Since
     * the life span of these
     * objects is small, we use a
     * very simple storage
     * technique.
     */
    std::ostream &stream;

    /**
     * The flags controlling the output
     */
    const DataOutBase::VtkFlags flags;
  };


  class VtuStream
  {
  public:
    /**
     * Constructor, storing
     * persistent values for
     * later use.
     */
    VtuStream (std::ostream &stream,
               const DataOutBase::VtkFlags flags);

    /**
     * Output operator for points.
     */
    template <int dim>
    void write_point (const unsigned int index,
                      const Point<dim> &);

    /**
     * Do whatever is necessary to
     * terminate the list of points.
     */
    void flush_points ();

    /**
     * Write dim-dimensional cell
     * with first vertex at
     * number start and further
     * vertices offset by the
     * specified values. Values
     * not needed are ignored.
     *
     * The order of vertices for
     * these cells in different
     * dimensions is
     * <ol>
     * <li> [0,1]
     * <li> []
     * <li> []
     * </ol>
     */
    template <int dim>
    void write_cell(const unsigned int index,
                    const unsigned int start,
                    const unsigned int x_offset,
                    const unsigned int y_offset,
                    const unsigned int z_offset);

    /**
     * Do whatever is necessary to
     * terminate the list of cells.
     */
    void flush_cells ();

    /**
     * Forwarding of output stream
     */
    template <typename T>
    std::ostream &operator<< (const T &);

    /**
     * Forwarding of output stream.
     *
     * If libz was found during
     * configuration, this operator
     * compresses and encodes the
     * entire data
     * block. Otherwise, it simply
     * writes it element by
     * element.
     */
    template <typename T>
    std::ostream &operator<< (const std::vector<T> &);

  private:
    /**
     * The ostream to use. Since
     * the life span of these
     * objects is small, we use a
     * very simple storage
     * technique.
     */
    std::ostream &stream;

    /**
     * The flags controlling the output
     */
    const DataOutBase::VtkFlags flags;

    /**
     * A list of vertices and
     * cells, to be used in case we
     * want to compress the data.
     *
     * The data types of these
     * arrays needs to match what
     * we print in the XML-preamble
     * to the respective parts of
     * VTU files (e.g. Float64 and
     * Int32)
     */
    std::vector<double>  vertices;
    std::vector<int32_t> cells;
  };


//----------------------------------------------------------------------//

  DXStream::DXStream(std::ostream &out,
                     const DataOutBase::DXFlags f)
    :
    stream(out), flags(f)
  {}


  template<int dim>
  void
  DXStream::write_point (const unsigned int,
                         const Point<dim> &p)
  {
    if (flags.coordinates_binary)
      {
        float data[dim];
        for (unsigned int d=0; d<dim; ++d)
          data[d] = p(d);
        stream.write(reinterpret_cast<const char *>(data),
                     dim * sizeof(*data));
      }
    else
      {
        for (unsigned int d=0; d<dim; ++d)
          stream << p(d) << '\t';
        stream << '\n';
      }
  }


  void
  DXStream::flush_points ()
  {}


  template<int dim>
  void
  DXStream::write_cell(
    unsigned int,
    unsigned int start,
    unsigned int d1,
    unsigned int d2,
    unsigned int d3)
  {
    int nodes[1<<dim];
    nodes[GeometryInfo<dim>::dx_to_deal[0]] = start;
    nodes[GeometryInfo<dim>::dx_to_deal[1]] = start+d1;
    if (dim>=2)
      {
        // Add shifted line in y direction
        nodes[GeometryInfo<dim>::dx_to_deal[2]] = start+d2;
        nodes[GeometryInfo<dim>::dx_to_deal[3]] = start+d2+d1;
        if (dim>=3)
          {
            // Add shifted quad in z direction
            nodes[GeometryInfo<dim>::dx_to_deal[4]] = start+d3;
            nodes[GeometryInfo<dim>::dx_to_deal[5]] = start+d3+d1;
            nodes[GeometryInfo<dim>::dx_to_deal[6]] = start+d3+d2;
            nodes[GeometryInfo<dim>::dx_to_deal[7]] = start+d3+d2+d1;
          }
      }

    if (flags.int_binary)
      stream.write(reinterpret_cast<const char *>(nodes),
                   (1<<dim) * sizeof(*nodes));
    else
      {
        const unsigned int final = (1<<dim) - 1;
        for (unsigned int i=0; i<final ; ++i)
          stream << nodes[i] << '\t';
        stream << nodes[final] << '\n';
      }
  }


  void
  DXStream::flush_cells ()
  {}


  template<typename data>
  inline
  void
  DXStream::write_dataset(const unsigned int     /*index*/,
                          const std::vector<data> &values)
  {
    if (flags.data_binary)
      {
        stream.write(reinterpret_cast<const char *>(&values[0]),
                     values.size()*sizeof(data));
      }
    else
      {
        for (unsigned int i=0; i<values.size(); ++i)
          stream << '\t' << values[i];
        stream << '\n';
      }
  }



//----------------------------------------------------------------------//

  GmvStream::GmvStream (std::ostream &out,
                        const DataOutBase::GmvFlags f)
    :
    selected_component(numbers::invalid_unsigned_int),
    stream(out), flags(f)
  {}


  template<int dim>
  void
  GmvStream::write_point (const unsigned int,
                          const Point<dim> &p)
  {
    Assert(selected_component != numbers::invalid_unsigned_int,
           ExcNotInitialized());
    stream << p(selected_component) << ' ';
  }


  void
  GmvStream::flush_points ()
  {}


  template<int dim>
  void
  GmvStream::write_cell(
    unsigned int,
    unsigned int s,
    unsigned int d1,
    unsigned int d2,
    unsigned int d3)
  {
    // Vertices are numbered starting
    // with one.
    const unsigned int start=s+1;
    stream << gmv_cell_type[dim] << '\n';

    stream << start << '\t'
           << start+d1;
    if (dim>=2)
      {
        stream << '\t' << start+d2+d1
               << '\t' << start+d2;
        if (dim>=3)
          {
            stream << '\t' << start+d3
                   << '\t' << start+d3+d1
                   << '\t' << start+d3+d2+d1
                   << '\t' << start+d3+d2;
          }
      }
    stream << '\n';
  }



  void
  GmvStream::flush_cells ()
  {}


//----------------------------------------------------------------------//

  TecplotStream::TecplotStream(std::ostream &out, const DataOutBase::TecplotFlags f)
    :
    selected_component(numbers::invalid_unsigned_int),
    stream(out), flags(f)
  {}


  template<int dim>
  void
  TecplotStream::write_point (const unsigned int,
                              const Point<dim> &p)
  {
    Assert(selected_component != numbers::invalid_unsigned_int,
           ExcNotInitialized());
    stream << p(selected_component) << '\n';
  }


  void
  TecplotStream::flush_points ()
  {}


  template<int dim>
  void
  TecplotStream::write_cell(
    unsigned int,
    unsigned int s,
    unsigned int d1,
    unsigned int d2,
    unsigned int d3)
  {
    const unsigned int start = s+1;

    stream << start << '\t'
           << start+d1;
    if (dim>=2)
      {
        stream << '\t' << start+d2+d1
               << '\t' << start+d2;
        if (dim>=3)
          {
            stream << '\t' << start+d3
                   << '\t' << start+d3+d1
                   << '\t' << start+d3+d2+d1
                   << '\t' << start+d3+d2;
          }
      }
    stream << '\n';
  }



  void
  TecplotStream::flush_cells ()
  {}



//----------------------------------------------------------------------//

  UcdStream::UcdStream(std::ostream &out, const DataOutBase::UcdFlags f)
    :
    stream(out), flags(f)
  {}


  template<int dim>
  void
  UcdStream::write_point (const unsigned int index,
                          const Point<dim> &p)
  {
    stream << index+1
           << "   ";
    // write out coordinates
    for (unsigned int i=0; i<dim; ++i)
      stream << p(i) << ' ';
    // fill with zeroes
    for (unsigned int i=dim; i<3; ++i)
      stream << "0 ";
    stream << '\n';
  }



  void
  UcdStream::flush_points ()
  {}


  template<int dim>
  void
  UcdStream::write_cell(
    unsigned int index,
    unsigned int start,
    unsigned int d1,
    unsigned int d2,
    unsigned int d3)
  {
    int nodes[1<<dim];
    nodes[GeometryInfo<dim>::ucd_to_deal[0]] = start;
    nodes[GeometryInfo<dim>::ucd_to_deal[1]] = start+d1;
    if (dim>=2)
      {
        // Add shifted line in y direction
        nodes[GeometryInfo<dim>::ucd_to_deal[2]] = start+d2;
        nodes[GeometryInfo<dim>::ucd_to_deal[3]] = start+d2+d1;
        if (dim>=3)
          {
            // Add shifted quad in z direction
            nodes[GeometryInfo<dim>::ucd_to_deal[4]] = start+d3;
            nodes[GeometryInfo<dim>::ucd_to_deal[5]] = start+d3+d1;
            nodes[GeometryInfo<dim>::ucd_to_deal[6]] = start+d3+d2;
            nodes[GeometryInfo<dim>::ucd_to_deal[7]] = start+d3+d2+d1;
          }
      }

    // Write out all cells and remember
    // that all indices must be shifted
    // by one.
    stream << index+1 << "\t0 " << ucd_cell_type[dim];
    const unsigned int final = (1<<dim);
    for (unsigned int i=0; i<final ; ++i)
      stream << '\t' << nodes[i]+1;
    stream << '\n';
  }



  void
  UcdStream::flush_cells ()
  {}


  template<typename data>
  inline
  void
  UcdStream::write_dataset(const unsigned int       index,
                           const std::vector<data> &values)
  {
    stream << index+1;
    for (unsigned int i=0; i<values.size(); ++i)
      stream << '\t' << values[i];
    stream << '\n';
  }



//----------------------------------------------------------------------//

  VtkStream::VtkStream(std::ostream &out, const DataOutBase::VtkFlags f)
    :
    stream(out), flags(f)
  {}


  template<int dim>
  void
  VtkStream::write_point (const unsigned int,
                          const Point<dim> &p)
  {
    // write out coordinates
    stream << p;
    // fill with zeroes
    for (unsigned int i=dim; i<3; ++i)
      stream << " 0";
    stream << '\n';
  }



  void
  VtkStream::flush_points ()
  {}


  template<int dim>
  void
  VtkStream::write_cell(
    unsigned int,
    unsigned int start,
    unsigned int d1,
    unsigned int d2,
    unsigned int d3)
  {
    stream << GeometryInfo<dim>::vertices_per_cell << '\t'
           << start << '\t'
           << start+d1;
    if (dim>=2)
      {
        stream << '\t' << start+d2+d1
               << '\t' << start+d2;
        if (dim>=3)
          {
            stream << '\t' << start+d3
                   << '\t' << start+d3+d1
                   << '\t' << start+d3+d2+d1
                   << '\t' << start+d3+d2;
          }
      }
    stream << '\n';
  }


  void
  VtkStream::flush_cells ()
  {}



  VtuStream::VtuStream(std::ostream &out, const DataOutBase::VtkFlags f)
    :
    stream(out), flags(f)
  {}


  template<int dim>
  void
  VtuStream::write_point (const unsigned int,
                          const Point<dim> &p)
  {
#if !defined(DEAL_II_WITH_ZLIB)
    // write out coordinates
    stream << p;
    // fill with zeroes
    for (unsigned int i=dim; i<3; ++i)
      stream << " 0";
    stream << '\n';
#else
    // if we want to compress, then
    // first collect all the data in
    // an array
    for (unsigned int i=0; i<dim; ++i)
      vertices.push_back(p[i]);
    for (unsigned int i=dim; i<3; ++i)
      vertices.push_back(0);
#endif
  }


  void
  VtuStream::flush_points ()
  {
#ifdef DEAL_II_WITH_ZLIB
    // compress the data we have in
    // memory and write them to the
    // stream. then release the data
    *this << vertices << '\n';
    vertices.clear ();
#endif
  }


  template<int dim>
  void
  VtuStream::write_cell(
    unsigned int,
    unsigned int start,
    unsigned int d1,
    unsigned int d2,
    unsigned int d3)
  {
#if !defined(DEAL_II_WITH_ZLIB)
    stream << start << '\t'
           << start+d1;
    if (dim>=2)
      {
        stream << '\t' << start+d2+d1
               << '\t' << start+d2;
        if (dim>=3)
          {
            stream << '\t' << start+d3
                   << '\t' << start+d3+d1
                   << '\t' << start+d3+d2+d1
                   << '\t' << start+d3+d2;
          }
      }
    stream << '\n';
#else
    cells.push_back (start);
    cells.push_back (start+d1);
    if (dim>=2)
      {
        cells.push_back (start+d2+d1);
        cells.push_back (start+d2);
        if (dim>=3)
          {
            cells.push_back (start+d3);
            cells.push_back (start+d3+d1);
            cells.push_back (start+d3+d2+d1);
            cells.push_back (start+d3+d2);
          }
      }
#endif
  }



  void
  VtuStream::flush_cells ()
  {
#ifdef DEAL_II_WITH_ZLIB
    // compress the data we have in
    // memory and write them to the
    // stream. then release the data
    *this << cells << '\n';
    cells.clear ();
#endif
  }


  template <typename T>
  std::ostream &
  VtuStream::operator<< (const std::vector<T> &data)
  {
#ifdef DEAL_II_WITH_ZLIB
    // compress the data we have in
    // memory and write them to the
    // stream. then release the data
    write_compressed_block (data, stream);
#else
    for (unsigned int i=0; i<data.size(); ++i)
      stream << data[i] << ' ';
#endif

    return stream;
  }

  template <typename T>
  std::ostream &
  DXStream::operator<< (const T &t)
  {
    stream << t;
    return stream;
  }


  template <typename T>
  std::ostream &
  GmvStream::operator<< (const T &t)
  {
    stream << t;
    return stream;
  }


  template <typename T>
  std::ostream &
  UcdStream::operator<< (const T &t)
  {
    stream << t;
    return stream;
  }


  template <typename T>
  std::ostream &
  VtkStream::operator<< (const T &t)
  {
    stream << t;
    return stream;
  }
}


//----------------------------------------------------------------------//

namespace DataOutBase
{
  template <int dim, int spacedim>
  const unsigned int Patch<dim,spacedim>::space_dim;

  const unsigned int Deal_II_IntermediateFlags::format_version;



  template <int dim, int spacedim>
  const unsigned int Patch<dim,spacedim>::no_neighbor;


  template <int dim, int spacedim>
  Patch<dim,spacedim>::Patch ()
    :
    patch_index(no_neighbor),
    n_subdivisions (1),
    points_are_available(false)
    // all the other data has a
    // constructor of its own, except
    // for the "neighbors" field, which
    // we set to invalid values.
  {
    for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
      neighbors[i] = no_neighbor;

    Assert (dim<=spacedim, ExcIndexRange(dim,0,spacedim));
    Assert (spacedim<=3, ExcNotImplemented());
  }



  template <int dim, int spacedim>
  bool
  Patch<dim,spacedim>::operator == (const Patch &patch) const
  {
//TODO: make tolerance relative
    const double epsilon=3e-16;
    for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
      if (vertices[i].distance(patch.vertices[i]) > epsilon)
        return false;

    for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
      if (neighbors[i] != patch.neighbors[i])
        return false;

    if (patch_index != patch.patch_index)
      return false;

    if (n_subdivisions != patch.n_subdivisions)
      return false;

    if (points_are_available != patch.points_are_available)
      return false;

    if (data.n_rows() != patch.data.n_rows())
      return false;

    if (data.n_cols() != patch.data.n_cols())
      return false;

    for (unsigned int i=0; i<data.n_rows(); ++i)
      for (unsigned int j=0; j<data.n_cols(); ++j)
        if (data[i][j] != patch.data[i][j])
          return false;

    return true;
  }



  template <int dim, int spacedim>
  std::size_t
  Patch<dim,spacedim>::memory_consumption () const
  {
    return (sizeof(vertices) / sizeof(vertices[0]) *
            MemoryConsumption::memory_consumption(vertices[0])
            +
            MemoryConsumption::memory_consumption(n_subdivisions)
            +
            MemoryConsumption::memory_consumption(data)
            +
            MemoryConsumption::memory_consumption(points_are_available));
  }



  UcdFlags::UcdFlags (const bool write_preamble)
    :
    write_preamble (write_preamble)
  {}



  PovrayFlags::PovrayFlags (const bool smooth,
                            const bool bicubic_patch,
                            const bool external_data)
    :
    smooth (smooth),
    bicubic_patch(bicubic_patch),
    external_data(external_data)
  {}


  DataOutFilterFlags::DataOutFilterFlags (const bool filter_duplicate_vertices,
                                          const bool xdmf_hdf5_output) :
    filter_duplicate_vertices(filter_duplicate_vertices),
    xdmf_hdf5_output(xdmf_hdf5_output)
  {}


  void DataOutFilterFlags::declare_parameters (ParameterHandler &prm)
  {
    prm.declare_entry ("Filter duplicate vertices", "false",
                       Patterns::Bool(),
                       "Whether to remove duplicate vertex values.");
    prm.declare_entry ("XDMF HDF5 output", "false",
                       Patterns::Bool(),
                       "Whether the data will be used in an XDMF/HDF5 combination.");
  }



  void DataOutFilterFlags::parse_parameters (const ParameterHandler &prm)
  {
    filter_duplicate_vertices = prm.get_bool ("Filter duplicate vertices");
    xdmf_hdf5_output = prm.get_bool ("XDMF HDF5 output");
  }



  std::size_t
  DataOutFilterFlags::memory_consumption () const
  {
    // only simple data elements, so
    // use sizeof operator
    return sizeof (*this);
  }



  DXFlags::DXFlags (const bool write_neighbors,
                    const bool int_binary,
                    const bool coordinates_binary,
                    const bool data_binary)
    :
    write_neighbors(write_neighbors),
    int_binary(int_binary),
    coordinates_binary(coordinates_binary),
    data_binary(data_binary),
    data_double(false)
  {}


  void DXFlags::declare_parameters (ParameterHandler &prm)
  {
    prm.declare_entry ("Write neighbors", "true",
                       Patterns::Bool(),
                       "A boolean field indicating whether neighborship "
                       "information between cells is to be written to the "
                       "OpenDX output file");
    prm.declare_entry ("Integer format", "ascii",
                       Patterns::Selection("ascii|32|64"),
                       "Output format of integer numbers, which is "
                       "either a text representation (ascii) or binary integer "
                       "values of 32 or 64 bits length");
    prm.declare_entry ("Coordinates format", "ascii",
                       Patterns::Selection("ascii|32|64"),
                       "Output format of vertex coordinates, which is "
                       "either a text representation (ascii) or binary "
                       "floating point values of 32 or 64 bits length");
    prm.declare_entry ("Data format", "ascii",
                       Patterns::Selection("ascii|32|64"),
                       "Output format of data values, which is "
                       "either a text representation (ascii) or binary "
                       "floating point values of 32 or 64 bits length");
  }



  void DXFlags::parse_parameters (const ParameterHandler &prm)
  {
    write_neighbors = prm.get_bool ("Write neighbors");
//TODO:[GK] Read the new  parameters
  }



  std::size_t
  DXFlags::memory_consumption () const
  {
    // only simple data elements, so
    // use sizeof operator
    return sizeof (*this);
  }




  void UcdFlags::declare_parameters (ParameterHandler &prm)
  {
    prm.declare_entry ("Write preamble", "true",
                       Patterns::Bool(),
                       "A flag indicating whether a comment should be "
                       "written to the beginning of the output file "
                       "indicating date and time of creation as well "
                       "as the creating program");
  }



  void UcdFlags::parse_parameters (const ParameterHandler &prm)
  {
    write_preamble = prm.get_bool ("Write preamble");
  }


  std::size_t
  UcdFlags::memory_consumption () const
  {
    // only simple data elements, so
    // use sizeof operator
    return sizeof (*this);
  }



  GnuplotFlags::GnuplotFlags ()
    :
    dummy (0)
  {}



  void GnuplotFlags::declare_parameters (ParameterHandler &/*prm*/)
  {}



  void GnuplotFlags::parse_parameters (const ParameterHandler &/*prm*/) const
  {}



  size_t
  GnuplotFlags::memory_consumption () const
  {
    // only simple data elements, so
    // use sizeof operator
    return sizeof (*this);
  }


  SvgFlags::SvgFlags (const unsigned int height_vector,
                      const int azimuth_angle,
                      const int polar_angle,
                      const unsigned int line_thickness,
                      const bool margin,
                      const bool draw_colorbar)
    :
    height(4000),
    width(0),
    height_vector(height_vector),
    azimuth_angle(azimuth_angle),
    polar_angle(polar_angle),
    line_thickness(line_thickness),
    margin(margin),
    draw_colorbar(draw_colorbar)
  {}


  std::size_t
  SvgFlags::memory_consumption () const
  {
    // only simple data elements, so
    // use sizeof operator
    return sizeof (*this);
  }




  void PovrayFlags::declare_parameters (ParameterHandler &prm)
  {
    prm.declare_entry ("Use smooth triangles", "false",
                       Patterns::Bool(),
                       "A flag indicating whether POVRAY should use smoothed "
                       "triangles instead of the usual ones");
    prm.declare_entry ("Use bicubic patches", "false",
                       Patterns::Bool(),
                       "Whether POVRAY should use bicubic patches");
    prm.declare_entry ("Include external file", "true",
                       Patterns::Bool (),
                       "Whether camera and lightling information should "
                       "be put into an external file \"data.inc\" or into "
                       "the POVRAY input file");
  }



  void PovrayFlags::parse_parameters (const ParameterHandler &prm)
  {
    smooth        = prm.get_bool ("Use smooth triangles");
    bicubic_patch = prm.get_bool ("Use bicubic patches");
    external_data = prm.get_bool ("Include external file");
  }



  std::size_t
  PovrayFlags::memory_consumption () const
  {
    // only simple data elements, so
    // use sizeof operator
    return sizeof (*this);
  }



  EpsFlags::EpsFlags (const unsigned int  height_vector,
                      const unsigned int  color_vector,
                      const SizeType      size_type,
                      const unsigned int  size,
                      const double        line_width,
                      const double        azimut_angle,
                      const double        turn_angle,
                      const double        z_scaling,
                      const bool          draw_mesh,
                      const bool          draw_cells,
                      const bool          shade_cells,
                      const ColorFunction color_function)
    :
    height_vector(height_vector),
    color_vector(color_vector),
    size_type(size_type),
    size(size),
    line_width(line_width),
    azimut_angle(azimut_angle),
    turn_angle(turn_angle),
    z_scaling(z_scaling),
    draw_mesh(draw_mesh),
    draw_cells(draw_cells),
    shade_cells(shade_cells),
    color_function(color_function)
  {}



  EpsFlags::RgbValues
  EpsFlags::default_color_function (const double x,
                                    const double xmin,
                                    const double xmax)
  {
    RgbValues rgb_values = { 0,0,0 };

// A difficult color scale:
//     xmin          = black  (1)
// 3/4*xmin+1/4*xmax = blue   (2)
// 1/2*xmin+1/2*xmax = green  (3)
// 1/4*xmin+3/4*xmax = red    (4)
//              xmax = white  (5)
// Makes the following color functions:
//
// red      green    blue
//       __
//      /      /\  /  /\    /
// ____/    __/  \/  /  \__/

//     { 0                                (1) - (3)
// r = { ( 4*x-2*xmin+2*xmax)/(xmax-xmin) (3) - (4)
//     { 1                                (4) - (5)
//
//     { 0                                (1) - (2)
// g = { ( 4*x-3*xmin-  xmax)/(xmax-xmin) (2) - (3)
//     { (-4*x+  xmin+3*xmax)/(xmax-xmin) (3) - (4)
//     { ( 4*x-  xmin-3*xmax)/(xmax-xmin) (4) - (5)
//
//     { ( 4*x-4*xmin       )/(xmax-xmin) (1) - (2)
// b = { (-4*x+2*xmin+2*xmax)/(xmax-xmin) (2) - (3)
//     { 0                                (3) - (4)
//     { ( 4*x-  xmin-3*xmax)/(xmax-xmin) (4) - (5)

    double sum   =   xmax+  xmin;
    double sum13 =   xmin+3*xmax;
    double sum22 = 2*xmin+2*xmax;
    double sum31 = 3*xmin+  xmax;
    double dif = xmax-xmin;
    double rezdif = 1.0/dif;

    int where;

    if (x<(sum31)/4)
      where = 0;
    else if (x<(sum22)/4)
      where = 1;
    else if (x<(sum13)/4)
      where = 2;
    else
      where = 3;

    if (dif!=0)
      {
        switch (where)
          {
          case 0:
            rgb_values.red   = 0;
            rgb_values.green = 0;
            rgb_values.blue  = (x-xmin)*4.*rezdif;
            break;
          case 1:
            rgb_values.red   = 0;
            rgb_values.green = (4*x-3*xmin-xmax)*rezdif;
            rgb_values.blue  = (sum22-4.*x)*rezdif;
            break;
          case 2:
            rgb_values.red   = (4*x-2*sum)*rezdif;
            rgb_values.green = (xmin+3*xmax-4*x)*rezdif;
            rgb_values.blue  = 0;
            break;
          case 3:
            rgb_values.red   = 1;
            rgb_values.green = (4*x-xmin-3*xmax)*rezdif;
            rgb_values.blue  = (4.*x-sum13)*rezdif;
          default:
            break;
          }
      }
    else // White
      rgb_values.red = rgb_values.green = rgb_values.blue = 1;

    return rgb_values;
  }



  EpsFlags::RgbValues
  EpsFlags::grey_scale_color_function (const double x,
                                       const double xmin,
                                       const double xmax)
  {
    EpsFlags::RgbValues rgb_values;
    rgb_values.red = rgb_values.blue = rgb_values.green
                                       = (x-xmin)/(xmax-xmin);
    return rgb_values;
  }



  EpsFlags::RgbValues
  EpsFlags::reverse_grey_scale_color_function (const double x,
                                               const double xmin,
                                               const double xmax)
  {
    EpsFlags::RgbValues rgb_values;
    rgb_values.red = rgb_values.blue = rgb_values.green
                                       = 1-(x-xmin)/(xmax-xmin);
    return rgb_values;
  }



  bool EpsCell2d::operator < (const EpsCell2d &e) const
  {
    // note the "wrong" order in
    // which we sort the elements
    return depth > e.depth;
  }



  void EpsFlags::declare_parameters (ParameterHandler &prm)
  {
    prm.declare_entry ("Index of vector for height", "0",
                       Patterns::Integer(),
                       "Number of the input vector that is to be used to "
                       "generate height information");
    prm.declare_entry ("Index of vector for color", "0",
                       Patterns::Integer(),
                       "Number of the input vector that is to be used to "
                       "generate color information");
    prm.declare_entry ("Scale to width or height", "width",
                       Patterns::Selection ("width|height"),
                       "Whether width or height should be scaled to match "
                       "the given size");
    prm.declare_entry ("Size (width or height) in eps units", "300",
                       Patterns::Integer(),
                       "The size (width or height) to which the eps output "
                       "file is to be scaled");
    prm.declare_entry ("Line widths in eps units", "0.5",
                       Patterns::Double(),
                       "The width in which the postscript renderer is to "
                       "plot lines");
    prm.declare_entry ("Azimut angle", "60",
                       Patterns::Double(0,180),
                       "Angle of the viewing position against the vertical "
                       "axis");
    prm.declare_entry ("Turn angle", "30",
                       Patterns::Double(0,360),
                       "Angle of the viewing direction against the y-axis");
    prm.declare_entry ("Scaling for z-axis", "1",
                       Patterns::Double (),
                       "Scaling for the z-direction relative to the scaling "
                       "used in x- and y-directions");
    prm.declare_entry ("Draw mesh lines", "true",
                       Patterns::Bool(),
                       "Whether the mesh lines, or only the surface should be "
                       "drawn");
    prm.declare_entry ("Fill interior of cells", "true",
                       Patterns::Bool(),
                       "Whether only the mesh lines, or also the interior of "
                       "cells should be plotted. If this flag is false, then "
                       "one can see through the mesh");
    prm.declare_entry ("Color shading of interior of cells", "true",
                       Patterns::Bool(),
                       "Whether the interior of cells shall be shaded");
    prm.declare_entry ("Color function", "default",
                       Patterns::Selection ("default|grey scale|reverse grey scale"),
                       "Name of a color function used to colorize mesh lines "
                       "and/or cell interiors");
  }



  void EpsFlags::parse_parameters (const ParameterHandler &prm)
  {
    height_vector = prm.get_integer ("Index of vector for height");
    color_vector  = prm.get_integer ("Index of vector for color");
    if (prm.get ("Scale to width or height") == "width")
      size_type   = width;
    else
      size_type   = height;
    size          = prm.get_integer ("Size (width or height) in eps units");
    line_width    = prm.get_double ("Line widths in eps units");
    azimut_angle  = prm.get_double ("Azimut angle");
    turn_angle    = prm.get_double ("Turn angle");
    z_scaling     = prm.get_double ("Scaling for z-axis");
    draw_mesh     = prm.get_bool ("Draw mesh lines");
    draw_cells    = prm.get_bool ("Fill interior of cells");
    shade_cells   = prm.get_bool ("Color shading of interior of cells");
    if (prm.get("Color function") == "default")
      color_function = &default_color_function;
    else if (prm.get("Color function") == "grey scale")
      color_function = &grey_scale_color_function;
    else if (prm.get("Color function") == "reverse grey scale")
      color_function = &reverse_grey_scale_color_function;
    else
      // we shouldn't get here, since
      // the parameter object should
      // already have checked that the
      // given value is valid
      Assert (false, ExcInternalError());
  }



  std::size_t
  EpsFlags::memory_consumption () const
  {
    // only simple data elements, so
    // use sizeof operator
    return sizeof (*this);
  }



  GmvFlags::GmvFlags ()
  {}



  void GmvFlags::declare_parameters (ParameterHandler &/*prm*/)
  {}



  void GmvFlags::parse_parameters (const ParameterHandler &/*prm*/) const
  {}


  std::size_t
  GmvFlags::memory_consumption () const
  {
    // only simple data elements, so
    // use sizeof operator
    return sizeof (*this);
  }



  TecplotFlags::
  TecplotFlags (const char *tecplot_binary_file_name,
                const char *zone_name)
    :
    tecplot_binary_file_name(tecplot_binary_file_name),
    zone_name(zone_name)
  {}



  void TecplotFlags::declare_parameters (ParameterHandler &/*prm*/)
  {}



  void TecplotFlags::parse_parameters (const ParameterHandler &/*prm*/) const
  {}


  std::size_t
  TecplotFlags::memory_consumption () const
  {
    // only simple data elements, so
    // use sizeof operator
    return sizeof (*this);
  }



  VtkFlags::VtkFlags (const double time,
                      const unsigned int cycle,
                      const bool print_date_and_time)
    :
    time (time),
    cycle (cycle),
    print_date_and_time (print_date_and_time)
  {}



  void VtkFlags::declare_parameters (ParameterHandler &/*prm*/)
  {}



  void VtkFlags::parse_parameters (const ParameterHandler &/*prm*/) const
  {}



  std::size_t
  VtkFlags::memory_consumption () const
  {
    // only simple data elements, so
    // use sizeof operator
    return sizeof (*this);
  }


  Deal_II_IntermediateFlags::Deal_II_IntermediateFlags ()
    :
    dummy (0)
  {}



  void Deal_II_IntermediateFlags::declare_parameters (ParameterHandler &/*prm*/)
  {}



  void Deal_II_IntermediateFlags::parse_parameters (const ParameterHandler &/*prm*/) const
  {}


  std::size_t
  Deal_II_IntermediateFlags::memory_consumption () const
  {
    // only simple data elements, so
    // use sizeof operator
    return sizeof (*this);
  }



  OutputFormat
  parse_output_format (const std::string &format_name)
  {
    if (format_name == "none")
      return none;

    if (format_name == "dx")
      return dx;

    if (format_name == "ucd")
      return ucd;

    if (format_name == "gnuplot")
      return gnuplot;

    if (format_name == "povray")
      return povray;

    if (format_name == "eps")
      return eps;

    if (format_name == "gmv")
      return gmv;

    if (format_name == "tecplot")
      return tecplot;

    if (format_name == "tecplot_binary")
      return tecplot_binary;

    if (format_name == "vtk")
      return vtk;

    if (format_name == "vtu")
      return vtu;

    if (format_name == "deal.II intermediate")
      return deal_II_intermediate;

    if (format_name == "hdf5")
      return hdf5;

    AssertThrow (false,
                 ExcMessage ("The given file format name is not recognized: <"
                             + format_name + ">"));

    // return something invalid
    return OutputFormat(-1);
  }



  std::string
  get_output_format_names ()
  {
    return "none|dx|ucd|gnuplot|povray|eps|gmv|tecplot|tecplot_binary|vtk|vtu|hdf5|svg|deal.II intermediate";
  }



  std::string
  default_suffix (const OutputFormat output_format)
  {
    switch (output_format)
      {
      case none:
        return "";
      case dx:
        return ".dx";
      case ucd:
        return ".inp";
      case gnuplot:
        return ".gnuplot";
      case povray:
        return ".pov";
      case eps:
        return ".eps";
      case gmv:
        return ".gmv";
      case tecplot:
        return ".dat";
      case tecplot_binary:
        return ".plt";
      case vtk:
        return ".vtk";
      case vtu:
        return ".vtu";
      case deal_II_intermediate:
        return ".d2";
      case hdf5:
        return ".h5";
      case svg:
        return ".svg";
      default:
        Assert (false, ExcNotImplemented());
        return "";
      }
  }


//----------------------------------------------------------------------//

  template <int dim, int spacedim, typename STREAM>
  void
  write_nodes (const std::vector<Patch<dim,spacedim> > &patches,
               STREAM &out)
  {
    Assert (dim<=3, ExcNotImplemented());
    unsigned int count = 0;
    // We only need this point below,
    // but it does not harm to declare
    // it here.
    Point<spacedim> node;

    for (typename std::vector<Patch<dim,spacedim> >::const_iterator
         patch=patches.begin();
         patch!=patches.end(); ++patch)
      {
        const unsigned int n_subdivisions = patch->n_subdivisions;
        const unsigned int n = n_subdivisions+1;
        // Length of loops in all
        // dimensions. If a dimension
        // is not used, a loop of
        // length one will do the job.
        const unsigned int n1 = (dim>0) ? n : 1;
        const unsigned int n2 = (dim>1) ? n : 1;
        const unsigned int n3 = (dim>2) ? n : 1;

        for (unsigned int i3=0; i3<n3; ++i3)
          for (unsigned int i2=0; i2<n2; ++i2)
            for (unsigned int i1=0; i1<n1; ++i1)
              {
                compute_node(node, &*patch,
                             i1,
                             i2,
                             i3,
                             n_subdivisions);
                out.write_point(count++, node);
              }
      }
    out.flush_points ();
  }

  template <int dim, int spacedim, typename STREAM>
  void
  write_cells (const std::vector<Patch<dim,spacedim> > &patches,
               STREAM &out)
  {
    Assert (dim<=3, ExcNotImplemented());
    unsigned int count = 0;
    unsigned int first_vertex_of_patch = 0;
    // Array to hold all the node
    // numbers of a cell. 8 is
    // sufficient for 3D
    for (typename std::vector<Patch<dim,spacedim> >::const_iterator
         patch=patches.begin();
         patch!=patches.end(); ++patch)
      {
        const unsigned int n_subdivisions = patch->n_subdivisions;
        const unsigned int n = n_subdivisions+1;
        // Length of loops in all dimensons
        const unsigned int n1 = (dim>0) ? n_subdivisions : 1;
        const unsigned int n2 = (dim>1) ? n_subdivisions : 1;
        const unsigned int n3 = (dim>2) ? n_subdivisions : 1;
        // Offsets of outer loops
        unsigned int d1 = 1;
        unsigned int d2 = n;
        unsigned int d3 = n*n;
        for (unsigned int i3=0; i3<n3; ++i3)
          for (unsigned int i2=0; i2<n2; ++i2)
            for (unsigned int i1=0; i1<n1; ++i1)
              {
                const unsigned int offset = first_vertex_of_patch+i3*d3+i2*d2+i1*d1;
                // First write line in x direction
                out.template write_cell<dim>(count++, offset, d1, d2, d3);
              }
        // finally update the number
        // of the first vertex of this patch
        first_vertex_of_patch += Utilities::fixed_power<dim>(n_subdivisions+1);
      }

    out.flush_cells ();
  }


  template <int dim, int spacedim, class STREAM>
  void
  write_data (
    const std::vector<Patch<dim,spacedim> > &patches,
    unsigned int n_data_sets,
    const bool double_precision,
    STREAM &out)
  {
    Assert (dim<=3, ExcNotImplemented());
    unsigned int count = 0;

    for (typename std::vector<Patch<dim,spacedim> >::const_iterator patch
         = patches.begin();
         patch != patches.end(); ++patch)
      {
        const unsigned int n_subdivisions = patch->n_subdivisions;
        const unsigned int n = n_subdivisions+1;
        // Length of loops in all dimensions
        Assert ((patch->data.n_rows() == n_data_sets && !patch->points_are_available) ||
                (patch->data.n_rows() == n_data_sets+spacedim && patch->points_are_available),
                ExcDimensionMismatch (patch->points_are_available
                                      ?
                                      (n_data_sets + spacedim)
                                      :
                                      n_data_sets,
                                      patch->data.n_rows()));
        Assert (patch->data.n_cols() == Utilities::fixed_power<dim>(n),
                ExcInvalidDatasetSize (patch->data.n_cols(), n));

        std::vector<float>  floats(n_data_sets);
        std::vector<double> doubles(n_data_sets);

        // Data is already in
        // lexicographic ordering
        for (unsigned int i=0; i<Utilities::fixed_power<dim>(n); ++i, ++count)
          if (double_precision)
            {
              for (unsigned int data_set=0; data_set<n_data_sets; ++data_set)
                doubles[data_set] = patch->data(data_set, i);
              out.write_dataset(count, doubles);
            }
          else
            {
              for (unsigned int data_set=0; data_set<n_data_sets; ++data_set)
                floats[data_set] = patch->data(data_set, i);
              out.write_dataset(count, floats);
            }
      }
  }



  namespace
  {
    /**
     * This function projects a three-dimensional point (Point<3> point)
     * onto a two-dimensional image plane, specified by the position of
     * the camera viewing system (Point<3> camera_position), camera
     * direction (Point<3> camera_position), camera horizontal (Point<3>
     * camera_horizontal, necessary for the correct alignment of the
     * later images), and the focus of the camera (float camera_focus).
     */
    Point<2> svg_project_point(Point<3> point, Point<3> camera_position, Point<3> camera_direction, Point<3> camera_horizontal, float camera_focus)
    {
      Point<3> camera_vertical;
      camera_vertical[0] = camera_horizontal[1] * camera_direction[2] - camera_horizontal[2] * camera_direction[1];
      camera_vertical[1] = camera_horizontal[2] * camera_direction[0] - camera_horizontal[0] * camera_direction[2];
      camera_vertical[2] = camera_horizontal[0] * camera_direction[1] - camera_horizontal[1] * camera_direction[0];

      float phi;
      phi  = camera_focus;
      phi /= (point[0] - camera_position[0]) * camera_direction[0] + (point[1] - camera_position[1]) * camera_direction[1] + (point[2] - camera_position[2]) * camera_direction[2];

      Point<3> projection;
      projection[0] = camera_position[0] + phi * (point[0] - camera_position[0]);
      projection[1] = camera_position[1] + phi * (point[1] - camera_position[1]);
      projection[2] = camera_position[2] + phi * (point[2] - camera_position[2]);

      Point<2> projection_decomposition;
      projection_decomposition[0]  = (projection[0] - camera_position[0] - camera_focus * camera_direction[0]) * camera_horizontal[0];
      projection_decomposition[0] += (projection[1] - camera_position[1] - camera_focus * camera_direction[1]) * camera_horizontal[1];
      projection_decomposition[0] += (projection[2] - camera_position[2] - camera_focus * camera_direction[2]) * camera_horizontal[2];

      projection_decomposition[1]  = (projection[0] - camera_position[0] - camera_focus * camera_direction[0]) * camera_vertical[0];
      projection_decomposition[1] += (projection[1] - camera_position[1] - camera_focus * camera_direction[1]) * camera_vertical[1];
      projection_decomposition[1] += (projection[2] - camera_position[2] - camera_focus * camera_direction[2]) * camera_vertical[2];

      return projection_decomposition;
    }


    /**
     * Function to compute the gradient parameters for a triangle with
     * given values for the vertices.
     */
    Point<6> svg_get_gradient_parameters(Point<3> points[])
    {
      Point<3> v_min, v_max, v_inter;

      // Use the Bubblesort algorithm to sort the points with respect to the third coordinate
      for (int i = 0; i < 2; ++i)
        {
          for (int j = 0; j < 2-i; ++j)
            {
              if (points[j][2] > points[j + 1][2])
                {
                  Point<3> temp = points[j];
                  points[j] = points[j+1];
                  points[j+1] = temp;
                }
            }
        }

      // save the related three-dimensional vectors v_min, v_inter, and v_max
      v_min = points[0];
      v_inter = points[1];
      v_max = points[2];

      Point<2> A[2];
      Point<2> b, gradient;

      // determine the plane offset c
      A[0][0] = v_max[0] - v_min[0];
      A[0][1] = v_inter[0] - v_min[0];
      A[1][0] = v_max[1] - v_min[1];
      A[1][1] = v_inter[1] - v_min[1];

      b[0] = - v_min[0];
      b[1] = - v_min[1];

      double x, sum;
      bool col_change = false;

      if (A[0][0] == 0)
        {
          col_change = true;

          A[0][0] = A[0][1];
          A[0][1] = 0;

          double temp = A[1][0];
          A[1][0] = A[1][1];
          A[1][1] = temp;
        }

      for (unsigned int k = 0; k < 1; k++)
        {
          for (unsigned int i = k+1; i < 2; i++)
            {
              x = A[i][k] / A[k][k];

              for (unsigned int j = k+1; j < 2; j++) A[i][j] = A[i][j] - A[k][j] * x;

              b[i] = b[i] - b[k]*x;

            }
        }

      b[1] = b[1] / A[1][1];

      for (int i = 0; i >= 0; i--)
        {
          sum = b[i];

          for (unsigned int j = i+1; j < 2; j++) sum = sum - A[i][j] * b[j];

          b[i] = sum / A[i][i];
        }

      if (col_change)
        {
          double temp = b[0];
          b[0] = b[1];
          b[1] = temp;
        }

      double c = b[0] * (v_max[2] - v_min[2]) + b[1] * (v_inter[2] - v_min[2]) + v_min[2];

      // Determine the first entry of the gradient (phi, cf. documentation)
      A[0][0] = v_max[0] - v_min[0];
      A[0][1] = v_inter[0] - v_min[0];
      A[1][0] = v_max[1] - v_min[1];
      A[1][1] = v_inter[1] - v_min[1];

      b[0] = 1.0 - v_min[0];
      b[1] = - v_min[1];

      col_change = false;

      if (A[0][0] == 0)
        {
          col_change = true;

          A[0][0] = A[0][1];
          A[0][1] = 0;

          double temp = A[1][0];
          A[1][0] = A[1][1];
          A[1][1] = temp;
        }

      for (unsigned int k = 0; k < 1; k++)
        {
          for (unsigned int i = k+1; i < 2; i++)
            {
              x = A[i][k] / A[k][k];

              for (unsigned int j = k+1; j < 2; j++) A[i][j] = A[i][j] - A[k][j] * x;

              b[i] = b[i] - b[k] * x;

            }
        }

      b[1]=b[1] / A[1][1];

      for (int i = 0; i >= 0; i--)
        {
          sum = b[i];

          for (unsigned int j = i+1; j < 2; j++) sum = sum - A[i][j]*b[j];

          b[i] = sum / A[i][i];
        }

      if (col_change)
        {
          double temp = b[0];
          b[0] = b[1];
          b[1] = temp;
        }

      gradient[0] = b[0] * (v_max[2] - v_min[2]) + b[1] * (v_inter[2] - v_min[2]) - c + v_min[2];

      // determine the second entry of the gradient
      A[0][0] = v_max[0] - v_min[0];
      A[0][1] = v_inter[0] - v_min[0];
      A[1][0] = v_max[1] - v_min[1];
      A[1][1] = v_inter[1] - v_min[1];

      b[0] = - v_min[0];
      b[1] = 1.0 - v_min[1];

      col_change = false;

      if (A[0][0] == 0)
        {
          col_change = true;

          A[0][0] = A[0][1];
          A[0][1] = 0;

          double temp = A[1][0];
          A[1][0] = A[1][1];
          A[1][1] = temp;
        }

      for (unsigned int k = 0; k < 1; k++)
        {
          for (unsigned int i = k+1; i < 2; i++)
            {
              x = A[i][k] / A[k][k];

              for (unsigned int j = k+1; j < 2; j++) A[i][j] = A[i][j] - A[k][j] * x;

              b[i] = b[i] - b[k] * x;
            }
        }

      b[1] = b[1] / A[1][1];

      for (int i = 0; i >= 0; i--)
        {
          sum = b[i];

          for (unsigned int j = i+1; j < 2; j++) sum = sum - A[i][j] * b[j];

          b[i] = sum / A[i][i];
        }

      if (col_change)
        {
          double temp = b[0];
          b[0] = b[1];
          b[1] = temp;
        }

      gradient[1] = b[0] * (v_max[2] - v_min[2]) + b[1] * (v_inter[2] - v_min[2]) - c + v_min[2];

      // normalize the gradient
      double gradient_norm = sqrt(pow(gradient[0], 2.0) + pow(gradient[1], 2.0));
      gradient[0] /= gradient_norm;
      gradient[1] /= gradient_norm;

      double lambda = - gradient[0] * (v_min[0] - v_max[0]) - gradient[1] * (v_min[1] - v_max[1]);

      Point<6> gradient_parameters(true);

      gradient_parameters[0] = v_min[0];
      gradient_parameters[1] = v_min[1];

      gradient_parameters[2] = v_min[0] + lambda * gradient[0];
      gradient_parameters[3] = v_min[1] + lambda * gradient[1];

      gradient_parameters[4] = v_min[2];
      gradient_parameters[5] = v_max[2];

      return gradient_parameters;
    }
  }



  template <int dim, int spacedim>
  void write_ucd (const std::vector<Patch<dim,spacedim> > &patches,
                  const std::vector<std::string>          &data_names,
                  const std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> > &,
                  const UcdFlags                          &flags,
                  std::ostream                            &out)
  {
    AssertThrow (out, ExcIO());

#ifndef DEAL_II_WITH_MPI
    // verify that there are indeed
    // patches to be written out. most
    // of the times, people just forget
    // to call build_patches when there
    // are no patches, so a warning is
    // in order. that said, the
    // assertion is disabled if we
    // support MPI since then it can
    // happen that on the coarsest
    // mesh, a processor simply has no
    // cells it actually owns, and in
    // that case it is legit if there
    // are no patches
    Assert (patches.size() > 0, ExcNoPatches());
#else
    if (patches.size() == 0)
      return;
#endif

    const unsigned int n_data_sets = data_names.size();

    UcdStream ucd_out(out, flags);

    // first count the number of cells
    // and cells for later use
    unsigned int n_nodes;
    unsigned int n_cells;
    compute_sizes<dim,spacedim> (patches, n_nodes, n_cells);
    ///////////////////////
    // preamble
    if (flags.write_preamble)
      {
        std::time_t  time1= std::time (0);
        std::tm     *time = std::localtime(&time1);
        out << "# This file was generated by the deal.II library." << '\n'
            << "# Date =  "
            << time->tm_year+1900 << "/"
            << time->tm_mon+1 << "/"
            << time->tm_mday << '\n'
            << "# Time =  "
            << time->tm_hour << ":"
            << std::setw(2) << time->tm_min << ":"
            << std::setw(2) << time->tm_sec << '\n'
            << "#" << '\n'
            << "# For a description of the UCD format see the AVS Developer's guide."
            << '\n'
            << "#" << '\n';
      }

    // start with ucd data
    out << n_nodes << ' '
        << n_cells << ' '
        << n_data_sets << ' '
        << 0 << ' '                  // no cell data at present
        << 0                         // no model data
        << '\n';

    write_nodes(patches, ucd_out);
    out << '\n';

    write_cells(patches, ucd_out);
    out << '\n';

    /////////////////////////////
    // now write data
    if (n_data_sets != 0)
      {
        out << n_data_sets << "    ";    // number of vectors
        for (unsigned int i=0; i<n_data_sets; ++i)
          out << 1 << ' ';               // number of components;
        // only 1 supported presently
        out << '\n';

        for (unsigned int data_set=0; data_set<n_data_sets; ++data_set)
          out << data_names[data_set]
              << ",dimensionless"      // no units supported at present
              << '\n';

        write_data(patches, n_data_sets, true, ucd_out);
      }
    // make sure everything now gets to
    // disk
    out.flush ();

    // assert the stream is still ok
    AssertThrow (out, ExcIO());
  }


  template <int dim, int spacedim>
  void write_dx (const std::vector<Patch<dim,spacedim> > &patches,
                 const std::vector<std::string>          &data_names,
                 const std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> > &,
                 const DXFlags                           &flags,
                 std::ostream                            &out)
  {
    AssertThrow (out, ExcIO());

#ifndef DEAL_II_WITH_MPI
    // verify that there are indeed
    // patches to be written out. most
    // of the times, people just forget
    // to call build_patches when there
    // are no patches, so a warning is
    // in order. that said, the
    // assertion is disabled if we
    // support MPI since then it can
    // happen that on the coarsest
    // mesh, a processor simply has no
    // cells it actually owns, and in
    // that case it is legit if there
    // are no patches
    Assert (patches.size() > 0, ExcNoPatches());
#else
    if (patches.size() == 0)
      return;
#endif
    // Stream with special features for dx output
    DXStream dx_out(out, flags);

    // Variable counting the offset of
    // binary data.
    unsigned int offset = 0;

    const unsigned int n_data_sets = data_names.size();

    // first count the number of cells
    // and cells for later use
    unsigned int n_nodes;
    unsigned int n_cells;
    compute_sizes<dim,spacedim>(patches, n_nodes, n_cells);
    // start with vertices order is
    // lexicographical, x varying
    // fastest
    out << "object \"vertices\" class array type float rank 1 shape " << spacedim
        << " items " << n_nodes;

    if (flags.coordinates_binary)
      {
        out << " lsb ieee data 0" << '\n';
        offset += n_nodes * spacedim *  sizeof(float);
      }
    else
      {
        out << " data follows" << '\n';
        write_nodes(patches, dx_out);
      }

    ///////////////////////////////
    // first write the coordinates of all vertices

    /////////////////////////////////////////
    // write cells
    out << "object \"cells\" class array type int rank 1 shape "
        << GeometryInfo<dim>::vertices_per_cell
        << " items " << n_cells;

    if (flags.int_binary)
      {
        out << " lsb binary data " << offset << '\n';
        offset += n_cells * sizeof (int);
      }
    else
      {
        out << " data follows" << '\n';
        write_cells(patches, dx_out);
        out << '\n';
      }


    out << "attribute \"element type\" string \"";
    if (dim==1) out << "lines";
    if (dim==2) out << "quads";
    if (dim==3) out << "cubes";
    out << "\"" << '\n'
        << "attribute \"ref\" string \"positions\"" << '\n';

//TODO:[GK] Patches must be of same size!
    /////////////////////////////
    // write neighbor information
    if (flags.write_neighbors)
      {
        out << "object \"neighbors\" class array type int rank 1 shape "
            << GeometryInfo<dim>::faces_per_cell
            << " items " << n_cells
            << " data follows";

        for (typename std::vector<Patch<dim,spacedim> >::const_iterator
             patch=patches.begin();
             patch!=patches.end(); ++patch)
          {
            const unsigned int n = patch->n_subdivisions;
            const unsigned int n1 = (dim>0) ? n : 1;
            const unsigned int n2 = (dim>1) ? n : 1;
            const unsigned int n3 = (dim>2) ? n : 1;
            unsigned int cells_per_patch = Utilities::fixed_power<dim>(n);
            unsigned int dx = 1;
            unsigned int dy = n;
            unsigned int dz = n*n;

            const unsigned int patch_start = patch->patch_index * cells_per_patch;

            for (unsigned int i3=0; i3<n3; ++i3)
              for (unsigned int i2=0; i2<n2; ++i2)
                for (unsigned int i1=0; i1<n1; ++i1)
                  {
                    const unsigned int nx = i1*dx;
                    const unsigned int ny = i2*dy;
                    const unsigned int nz = i3*dz;

                    out << '\n';
                    // Direction -x
                    // Last cell in row
                    // of other patch
                    if (i1==0)
                      {
                        const unsigned int nn = patch->neighbors[0];
                        out << '\t';
                        if (nn != patch->no_neighbor)
                          out << (nn*cells_per_patch+ny+nz+dx*(n-1));
                        else
                          out << "-1";
                      }
                    else
                      {
                        out << '\t'
                            << patch_start+nx-dx+ny+nz;
                      }
                    // Direction +x
                    // First cell in row
                    // of other patch
                    if (i1 == n-1)
                      {
                        const unsigned int nn = patch->neighbors[1];
                        out << '\t';
                        if (nn != patch->no_neighbor)
                          out << (nn*cells_per_patch+ny+nz);
                        else
                          out << "-1";
                      }
                    else
                      {
                        out << '\t'
                            << patch_start+nx+dx+ny+nz;
                      }
                    if (dim<2)
                      continue;
                    // Direction -y
                    if (i2==0)
                      {
                        const unsigned int nn = patch->neighbors[2];
                        out << '\t';
                        if (nn != patch->no_neighbor)
                          out << (nn*cells_per_patch+nx+nz+dy*(n-1));
                        else
                          out << "-1";
                      }
                    else
                      {
                        out << '\t'
                            << patch_start+nx+ny-dy+nz;
                      }
                    // Direction +y
                    if (i2 == n-1)
                      {
                        const unsigned int nn = patch->neighbors[3];
                        out << '\t';
                        if (nn != patch->no_neighbor)
                          out << (nn*cells_per_patch+nx+nz);
                        else
                          out << "-1";
                      }
                    else
                      {
                        out << '\t'
                            << patch_start+nx+ny+dy+nz;
                      }
                    if (dim<3)
                      continue;

                    // Direction -z
                    if (i3==0)
                      {
                        const unsigned int nn = patch->neighbors[4];
                        out << '\t';
                        if (nn != patch->no_neighbor)
                          out << (nn*cells_per_patch+nx+ny+dz*(n-1));
                        else
                          out << "-1";
                      }
                    else
                      {
                        out << '\t'
                            << patch_start+nx+ny+nz-dz;
                      }
                    // Direction +z
                    if (i3 == n-1)
                      {
                        const unsigned int nn = patch->neighbors[5];
                        out << '\t';
                        if (nn != patch->no_neighbor)
                          out << (nn*cells_per_patch+nx+ny);
                        else
                          out << "-1";
                      }
                    else
                      {
                        out << '\t'
                            << patch_start+nx+ny+nz+dz;
                      }
                  }
            out << '\n';
          }
      }
    /////////////////////////////
    // now write data
    if (n_data_sets != 0)
      {
        out << "object \"data\" class array type float rank 1 shape "
            << n_data_sets
            << " items " << n_nodes;

        if (flags.data_binary)
          {
            out << " lsb ieee data " << offset << '\n';
            offset += n_data_sets * n_nodes * ((flags.data_double)
                                               ? sizeof(double)
                                               : sizeof(float));
          }
        else
          {
            out << " data follows" << '\n';
            write_data(patches, n_data_sets, flags.data_double, dx_out);
          }

        // loop over all patches
        out << "attribute \"dep\" string \"positions\"" << '\n';
      }
    else
      {
        out << "object \"data\" class constantarray type float rank 0 items " << n_nodes << " data follows"
            << '\n' << '0' << '\n';
      }

    // no model data

    out << "object \"deal data\" class field" << '\n'
        << "component \"positions\" value \"vertices\"" << '\n'
        << "component \"connections\" value \"cells\"" << '\n'
        << "component \"data\" value \"data\"" << '\n';

    if (flags.write_neighbors)
      out << "component \"neighbors\" value \"neighbors\"" << '\n';

    if (true)
      {
        std::time_t  time1= std::time (0);
        std::tm     *time = std::localtime(&time1);
        out << "attribute \"created\" string \""
            << time->tm_year+1900 << "/"
            << time->tm_mon+1 << "/"
            << time->tm_mday
            << ' '
            << time->tm_hour << ":"
            << std::setw(2) << time->tm_min << ":"
            << std::setw(2) << time->tm_sec << '"' << '\n';
      }

    out << "end" << '\n';
    // Write all binary data now
    if (flags.coordinates_binary)
      write_nodes(patches, dx_out);
    if (flags.int_binary)
      write_cells(patches, dx_out);
    if (flags.data_binary)
      write_data(patches, n_data_sets, flags.data_double, dx_out);

    // make sure everything now gets to
    // disk
    out.flush ();

    // assert the stream is still ok
    AssertThrow (out, ExcIO());
  }



  template <int dim, int spacedim>
  void write_gnuplot (const std::vector<Patch<dim,spacedim> > &patches,
                      const std::vector<std::string>          &data_names,
                      const std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> > &,
                      const GnuplotFlags                      &/*flags*/,
                      std::ostream                            &out)
  {
    AssertThrow (out, ExcIO());

#ifndef DEAL_II_WITH_MPI
    // verify that there are indeed
    // patches to be written out. most
    // of the times, people just forget
    // to call build_patches when there
    // are no patches, so a warning is
    // in order. that said, the
    // assertion is disabled if we
    // support MPI since then it can
    // happen that on the coarsest
    // mesh, a processor simply has no
    // cells it actually owns, and in
    // that case it is legit if there
    // are no patches
    Assert (patches.size() > 0, ExcNoPatches());
#else
    if (patches.size() == 0)
      return;
#endif

    const unsigned int n_data_sets = data_names.size();

    // write preamble
    if (true)
      {
        // block this to have local
        // variables destroyed after
        // use
        const std::time_t  time1= std::time (0);
        const std::tm     *time = std::localtime(&time1);
        out << "# This file was generated by the deal.II library." << '\n'
            << "# Date =  "
            << time->tm_year+1900 << "/"
            << time->tm_mon+1 << "/"
            << time->tm_mday << '\n'
            << "# Time =  "
            << time->tm_hour << ":"
            << std::setw(2) << time->tm_min << ":"
            << std::setw(2) << time->tm_sec << '\n'
            << "#" << '\n'
            << "# For a description of the GNUPLOT format see the GNUPLOT manual."
            << '\n'
            << "#" << '\n'
            << "# ";

        switch (spacedim)
          {
          case 1:
            out << "<x> ";
            break;
          case 2:
            out << "<x> <y> ";
            break;
          case 3:
            out << "<x> <y> <z> ";
            break;

          default:
            Assert (false, ExcNotImplemented());
          }

        for (unsigned int i=0; i<data_names.size(); ++i)
          out << '<' << data_names[i] << "> ";
        out << '\n';
      }


    // loop over all patches
    for (typename std::vector<Patch<dim,spacedim> >::const_iterator patch=patches.begin();
         patch != patches.end(); ++patch)
      {
        const unsigned int n_subdivisions = patch->n_subdivisions;
        const unsigned int n = n_subdivisions+1;
        // Length of loops in all dimensions
        const unsigned int n1 = (dim>0) ? n : 1;
        const unsigned int n2 = (dim>1) ? n : 1;
        const unsigned int n3 = (dim>2) ? n : 1;
        unsigned int d1 = 1;
        unsigned int d2 = n;
        unsigned int d3 = n*n;

        Assert ((patch->data.n_rows() == n_data_sets && !patch->points_are_available) ||
                (patch->data.n_rows() == n_data_sets+spacedim && patch->points_are_available),
                ExcDimensionMismatch (patch->points_are_available
                                      ?
                                      (n_data_sets + spacedim)
                                      :
                                      n_data_sets,
                                      patch->data.n_rows()));
        Assert (patch->data.n_cols() == Utilities::fixed_power<dim>(n),
                ExcInvalidDatasetSize (patch->data.n_cols(), n_subdivisions+1));

        Point<spacedim> this_point;
        Point<spacedim> node;
        if (dim<3)
          {
            for (unsigned int i2=0; i2<n2; ++i2)
              {
                for (unsigned int i1=0; i1<n1; ++i1)
                  {
                    // compute coordinates for
                    // this patch point
                    compute_node(node, &*patch, i1, i2, 0, n_subdivisions);
                    out << node << ' ';

                    for (unsigned int data_set=0; data_set<n_data_sets; ++data_set)
                      out << patch->data(data_set,i1*d1+i2*d2) << ' ';
                    out << '\n';
                  }
                // end of row in patch
                if (dim>1)
                  out << '\n';
              }
            // end of patch
            if (dim==1)
              out << '\n';
            out << '\n';
          }
        else if (dim==3)
          {
            // for all grid points: draw
            // lines into all positive
            // coordinate directions if
            // there is another grid point
            // there
            for (unsigned int i3=0; i3<n3; ++i3)
              for (unsigned int i2=0; i2<n2; ++i2)
                for (unsigned int i1=0; i1<n1; ++i1)
                  {
                    // compute coordinates for
                    // this patch point
                    compute_node(this_point, &*patch, i1, i2, i3, n_subdivisions);
                    // line into positive x-direction
                    // if possible
                    if (i1 < n_subdivisions)
                      {
                        // write point here
                        // and its data
                        out << this_point;
                        for (unsigned int data_set=0; data_set<n_data_sets; ++data_set)
                          out  << ' '
                               << patch->data(data_set,i1*d1+i2*d2+i3*d3);
                        out << '\n';

                        // write point there
                        // and its data
                        compute_node(node, &*patch, i1+1, i2, i3, n_subdivisions);
                        out << node;

                        for (unsigned int data_set=0; data_set<n_data_sets; ++data_set)
                          out  << ' '
                               << patch->data(data_set,(i1+1)*d1+i2*d2+i3*d3);
                        out << '\n';

                        // end of line
                        out << '\n'
                            << '\n';
                      }

                    // line into positive y-direction
                    // if possible
                    if (i2 < n_subdivisions)
                      {
                        // write point here
                        // and its data
                        out << this_point;
                        for (unsigned int data_set=0; data_set<n_data_sets; ++data_set)
                          out  << ' '
                               << patch->data(data_set, i1*d1+i2*d2+i3*d3);
                        out << '\n';

                        // write point there
                        // and its data
                        compute_node(node, &*patch, i1, i2+1, i3, n_subdivisions);
                        out << node;

                        for (unsigned int data_set=0; data_set<n_data_sets; ++data_set)
                          out  << ' '
                               << patch->data(data_set,i1*d1+(i2+1)*d2+i3*d3);
                        out << '\n';

                        // end of line
                        out << '\n'
                            << '\n';
                      }

                    // line into positive z-direction
                    // if possible
                    if (i3 < n_subdivisions)
                      {
                        // write point here
                        // and its data
                        out << this_point;
                        for (unsigned int data_set=0; data_set<n_data_sets; ++data_set)
                          out  << ' '
                               << patch->data(data_set,i1*d1+i2*d2+i3*d3);
                        out << '\n';

                        // write point there
                        // and its data
                        compute_node(node, &*patch, i1, i2, i3+1, n_subdivisions);
                        out << node;

                        for (unsigned int data_set=0; data_set<n_data_sets; ++data_set)
                          out  << ' '
                               << patch->data(data_set,i1*d1+i2*d2+(i3+1)*d3);
                        out << '\n';
                        // end of line
                        out << '\n'
                            << '\n';
                      }

                  }
          }
        else
          Assert (false, ExcNotImplemented());
      }
    // make sure everything now gets to
    // disk
    out.flush ();

    AssertThrow (out, ExcIO());
  }



  template <int dim, int spacedim>
  void write_povray (const std::vector<Patch<dim,spacedim> > &patches,
                     const std::vector<std::string>          &data_names,
                     const std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> > &,
                     const PovrayFlags                       &flags,
                     std::ostream                            &out)
  {
    AssertThrow (out, ExcIO());

#ifndef DEAL_II_WITH_MPI
    // verify that there are indeed
    // patches to be written out. most
    // of the times, people just forget
    // to call build_patches when there
    // are no patches, so a warning is
    // in order. that said, the
    // assertion is disabled if we
    // support MPI since then it can
    // happen that on the coarsest
    // mesh, a processor simply has no
    // cells it actually owns, and in
    // that case it is legit if there
    // are no patches
    Assert (patches.size() > 0, ExcNoPatches());
#else
    if (patches.size() == 0)
      return;
#endif
    Assert (dim==2, ExcNotImplemented());        // only for 2-D surfaces on a 2-D plane
    Assert (spacedim==2, ExcNotImplemented());

    const unsigned int n_data_sets = data_names.size();

    // write preamble
    if (true)
      {
        // block this to have local
        // variables destroyed after use
        const std::time_t  time1= std::time (0);
        const std::tm     *time = std::localtime(&time1);
        out << "/* This file was generated by the deal.II library." << '\n'
            << "   Date =  "
            << time->tm_year+1900 << "/"
            << time->tm_mon+1 << "/"
            << time->tm_mday << '\n'
            << "   Time =  "
            << time->tm_hour << ":"
            << std::setw(2) << time->tm_min << ":"
            << std::setw(2) << time->tm_sec << '\n'
            << '\n'
            << "   For a description of the POVRAY format see the POVRAY manual."
            << '\n'
            << "*/ " << '\n';

        // include files
        out << "#include \"colors.inc\" " << '\n'
            << "#include \"textures.inc\" " << '\n';


        // use external include file for textures,
        // camera and light
        if (flags.external_data)
          out << "#include \"data.inc\" " << '\n';
        else                          // all definitions in data file
          {
            // camera
            out << '\n' << '\n'
                << "camera {"            << '\n'
                << "  location <1,4,-7>" << '\n'
                << "  look_at <0,0,0>"   << '\n'
                << "  angle 30"          << '\n'
                << "}"                   << '\n';

            // light
            out << '\n'
                << "light_source {"      << '\n'
                << "  <1,4,-7>"      << '\n'
                << "  color Grey"        << '\n'
                << "}"                   << '\n';
            out << '\n'
                << "light_source {"      << '\n'
                << "  <0,20,0>"      << '\n'
                << "  color White"       << '\n'
                << "}"                   << '\n';
          }
      }

    // max. and min. heigth of solution
    Assert(patches.size()>0, ExcInternalError());
    double hmin=patches[0].data(0,0);
    double hmax=patches[0].data(0,0);

    for (typename std::vector<Patch<dim,spacedim> >::const_iterator patch=patches.begin();
         patch != patches.end(); ++patch)
      {
        const unsigned int n_subdivisions = patch->n_subdivisions;

        Assert ((patch->data.n_rows() == n_data_sets && !patch->points_are_available) ||
                (patch->data.n_rows() == n_data_sets+spacedim && patch->points_are_available),
                ExcDimensionMismatch (patch->points_are_available
                                      ?
                                      (n_data_sets + spacedim)
                                      :
                                      n_data_sets,
                                      patch->data.n_rows()));
        Assert (patch->data.n_cols() == Utilities::fixed_power<dim>(n_subdivisions+1),
                ExcInvalidDatasetSize (patch->data.n_cols(), n_subdivisions+1));

        for (unsigned int i=0; i<n_subdivisions+1; ++i)
          for (unsigned int j=0; j<n_subdivisions+1; ++j)
            {
              const int dl = i*(n_subdivisions+1)+j;
              if (patch->data(0,dl)<hmin)
                hmin=patch->data(0,dl);
              if (patch->data(0,dl)>hmax)
                hmax=patch->data(0,dl);
            }
      }

    out << "#declare HMIN=" << hmin << ";" << '\n'
        << "#declare HMAX=" << hmax << ";" << '\n' << '\n';

    if (!flags.external_data)
      {
        // texture with scaled niveau lines
        // 10 lines in the surface
        out << "#declare Tex=texture{" << '\n'
            << "  pigment {" << '\n'
            << "    gradient y" << '\n'
            << "    scale y*(HMAX-HMIN)*" << 0.1 << '\n'
            << "    color_map {" << '\n'
            << "      [0.00 color Light_Purple] " << '\n'
            << "      [0.95 color Light_Purple] " << '\n'
            << "      [1.00 color White]    " << '\n'
            << "} } }" << '\n' << '\n';
      }

    if (!flags.bicubic_patch)
      {
        // start of mesh header
        out << '\n'
            << "mesh {" << '\n';
      }

    // loop over all patches
    for (typename std::vector<Patch<dim,spacedim> >::const_iterator patch=patches.begin();
         patch != patches.end(); ++patch)
      {
        const unsigned int n_subdivisions = patch->n_subdivisions;
        const unsigned int n = n_subdivisions+1;
        const unsigned int d1=1;
        const unsigned int d2=n;

        Assert ((patch->data.n_rows() == n_data_sets && !patch->points_are_available) ||
                (patch->data.n_rows() == n_data_sets+spacedim && patch->points_are_available),
                ExcDimensionMismatch (patch->points_are_available
                                      ?
                                      (n_data_sets + spacedim)
                                      :
                                      n_data_sets,
                                      patch->data.n_rows()));
        Assert (patch->data.n_cols() == Utilities::fixed_power<dim>(n),
                ExcInvalidDatasetSize (patch->data.n_cols(), n_subdivisions+1));


        std::vector<Point<spacedim> > ver(n*n);

        for (unsigned int i2=0; i2<n; ++i2)
          for (unsigned int i1=0; i1<n; ++i1)
            {
              // compute coordinates for
              // this patch point, storing in ver
              compute_node(ver[i1*d1+i2*d2], &*patch, i1, i2, 0, n_subdivisions);
            }


        if (!flags.bicubic_patch)
          {
            // approximate normal
            // vectors in patch
            std::vector<Point<3> > nrml;
            // only if smooth triangles are used
            if (flags.smooth)
              {
                nrml.resize(n*n);
                // These are
                // difference
                // quotients of
                // the surface
                // mapping. We
                // take them
                // symmetric
                // inside the
                // patch and
                // one-sided at
                // the edges
                Point<3> h1,h2;
                // Now compute normals in every point
                for (unsigned int i=0; i<n; ++i)
                  for (unsigned int j=0; j<n; ++j)
                    {
                      const unsigned int il = (i==0) ? i : (i-1);
                      const unsigned int ir = (i==n_subdivisions) ? i : (i+1);
                      const unsigned int jl = (j==0) ? j : (j-1);
                      const unsigned int jr = (j==n_subdivisions) ? j : (j+1);

                      h1(0)=ver[ir*d1+j*d2](0) - ver[il*d1+j*d2](0);
                      h1(1)=patch->data(0,ir*d1+j*d2)-
                            patch->data(0,il*d1+j*d2);
                      h1(2)=ver[ir*d1+j*d2](1) - ver[il*d1+j*d2](1);

                      h2(0)=ver[i*d1+jr*d2](0) - ver[i*d1+jl*d2](0);
                      h2(1)=patch->data(0,i*d1+jr*d2)-
                            patch->data(0,i*d1+jl*d2);
                      h2(2)=ver[i*d1+jr*d2](1) - ver[i*d1+jl*d2](1);

                      nrml[i*d1+j*d2](0)=h1(1)*h2(2)-h1(2)*h2(1);
                      nrml[i*d1+j*d2](1)=h1(2)*h2(0)-h1(0)*h2(2);
                      nrml[i*d1+j*d2](2)=h1(0)*h2(1)-h1(1)*h2(0);

                      // normalize Vector
                      double norm=std::sqrt(
                                    std::pow(nrml[i*d1+j*d2](0),2.)+
                                    std::pow(nrml[i*d1+j*d2](1),2.)+
                                    std::pow(nrml[i*d1+j*d2](2),2.));

                      if (nrml[i*d1+j*d2](1)<0)
                        norm*=-1.;

                      for (unsigned int k=0; k<3; ++k)
                        nrml[i*d1+j*d2](k)/=norm;
                    }
              }

            // setting up triangles
            for (unsigned int i=0; i<n_subdivisions; ++i)
              for (unsigned int j=0; j<n_subdivisions; ++j)
                {
                  // down/left vertex of triangle
                  const int dl = i*d1+j*d2;
                  if (flags.smooth)
                    {
                      // writing smooth_triangles

                      // down/right triangle
                      out << "smooth_triangle {" << '\n' << "\t<"
                          << ver[dl](0) << ","
                          << patch->data(0,dl) << ","
                          << ver[dl](1) << ">, <"
                          << nrml[dl](0) << ", "
                          << nrml[dl](1) << ", "
                          << nrml[dl](2)
                          << ">," << '\n';
                      out << " \t<"
                          << ver[dl+d1](0) << ","
                          << patch->data(0,dl+d1)  << ","
                          << ver[dl+d1](1) << ">, <"
                          << nrml[dl+d1](0) << ", "
                          << nrml[dl+d1](1) << ", "
                          << nrml[dl+d1](2)
                          << ">," << '\n';
                      out << "\t<"
                          << ver[dl+d1+d2](0) << ","
                          << patch->data(0,dl+d1+d2)  << ","
                          << ver[dl+d1+d2](1) << ">, <"
                          << nrml[dl+d1+d2](0) << ", "
                          << nrml[dl+d1+d2](1) << ", "
                          << nrml[dl+d1+d2](2)
                          << ">}" << '\n';

                      // upper/left triangle
                      out << "smooth_triangle {" << '\n' << "\t<"
                          << ver[dl](0) << ","
                          << patch->data(0,dl) << ","
                          << ver[dl](1) << ">, <"
                          << nrml[dl](0) << ", "
                          << nrml[dl](1) << ", "
                          << nrml[dl](2)
                          << ">," << '\n';
                      out << "\t<"
                          << ver[dl+d1+d2](0) << ","
                          << patch->data(0,dl+d1+d2)  << ","
                          << ver[dl+d1+d2](1) << ">, <"
                          << nrml[dl+d1+d2](0) << ", "
                          << nrml[dl+d1+d2](1) << ", "
                          << nrml[dl+d1+d2](2)
                          << ">," << '\n';
                      out << "\t<"
                          << ver[dl+d2](0) << ","
                          << patch->data(0,dl+d2)  << ","
                          << ver[dl+d2](1) << ">, <"
                          << nrml[dl+d2](0) << ", "
                          << nrml[dl+d2](1) << ", "
                          << nrml[dl+d2](2)
                          << ">}" << '\n';
                    }
                  else
                    {
                      // writing standard triangles
                      // down/right triangle
                      out << "triangle {" << '\n' << "\t<"
                          << ver[dl](0) << ","
                          << patch->data(0,dl) << ","
                          << ver[dl](1) << ">," << '\n';
                      out << "\t<"
                          << ver[dl+d1](0) << ","
                          << patch->data(0,dl+d1)  << ","
                          << ver[dl+d1](1) << ">," << '\n';
                      out << "\t<"
                          << ver[dl+d1+d2](0) << ","
                          << patch->data(0,dl+d1+d2)  << ","
                          << ver[dl+d1+d2](1) << ">}" << '\n';

                      // upper/left triangle
                      out << "triangle {" << '\n' << "\t<"
                          << ver[dl](0) << ","
                          << patch->data(0,dl) << ","
                          << ver[dl](1) << ">," << '\n';
                      out << "\t<"
                          << ver[dl+d1+d2](0) << ","
                          << patch->data(0,dl+d1+d2)  << ","
                          << ver[dl+d1+d2](1) << ">," << '\n';
                      out << "\t<"
                          << ver[dl+d2](0) << ","
                          << patch->data(0,dl+d2)  << ","
                          << ver[dl+d2](1) << ">}" << '\n';
                    }
                }
          }
        else
          {
            // writing bicubic_patch
            Assert (n_subdivisions==3, ExcDimensionMismatch(n_subdivisions,3));
            out << '\n'
                << "bicubic_patch {" << '\n'
                << "  type 0" << '\n'
                << "  flatness 0" << '\n'
                << "  u_steps 0" << '\n'
                << "  v_steps 0" << '\n';
            for (int i=0; i<16; ++i)
              {
                out << "\t<" << ver[i](0) << "," << patch->data(0,i) << "," << ver[i](1) << ">";
                if (i!=15) out << ",";
                out << '\n';
              }
            out << "  texture {Tex}" <<  '\n'
                << "}" << '\n';
          }
      }

    if (!flags.bicubic_patch)
      {
        // the end of the mesh
        out << "  texture {Tex}" << '\n'
            << "}" << '\n'
            << '\n';
      }

    // make sure everything now gets to
    // disk
    out.flush ();

    AssertThrow (out, ExcIO());
  }



  template <int dim, int spacedim>
  void write_eps (const std::vector<Patch<dim,spacedim> > &patches,
                  const std::vector<std::string>          &/*data_names*/,
                  const std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> > &,
                  const EpsFlags                          &flags,
                  std::ostream                            &out)
  {
    // not implemented, see the documentation of the function
    AssertThrow (dim==2, ExcNotImplemented());
  }


  template <int spacedim>
  void write_eps (const std::vector<Patch<2,spacedim> > &patches,
                  const std::vector<std::string>          &/*data_names*/,
                  const std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> > &,
                  const EpsFlags                          &flags,
                  std::ostream                            &out)
  {
    Assert (out, ExcIO());

#ifndef DEAL_II_WITH_MPI
    // verify that there are indeed
    // patches to be written out. most
    // of the times, people just forget
    // to call build_patches when there
    // are no patches, so a warning is
    // in order. that said, the
    // assertion is disabled if we
    // support MPI since then it can
    // happen that on the coarsest
    // mesh, a processor simply has no
    // cells it actually owns, and in
    // that case it is legit if there
    // are no patches
    Assert (patches.size() > 0, ExcNoPatches());
#else
    if (patches.size() == 0)
      return;
#endif

    const unsigned int old_precision = out.precision();

    // set up an array of cells to be
    // written later. this array holds the
    // cells of all the patches as
    // projected to the plane perpendicular
    // to the line of sight.
    //
    // note that they are kept sorted by
    // the set, where we chose the value
    // of the center point of the cell
    // along the line of sight as value
    // for sorting
    std::multiset<EpsCell2d> cells;

    // two variables in which we
    // will store the minimum and
    // maximum values of the field
    // to be used for colorization
    //
    // preset them by 0 to calm down the
    // compiler; they are initialized later
    double min_color_value=0, max_color_value=0;

    // Array for z-coordinates of points.
    // The elevation determined by a function if spacedim=2
    // or the z-cooridate of the grid point if spacedim=3
    double heights[4] = { 0, 0, 0, 0 };

    // compute the cells for output and
    // enter them into the set above
    // note that since dim==2, we
    // have exactly four vertices per
    // patch and per cell
    for (typename std::vector<Patch<2,spacedim> >::const_iterator patch=patches.begin();
         patch!=patches.end(); ++patch)
      {
        const unsigned int n_subdivisions = patch->n_subdivisions;
        const unsigned int n = n_subdivisions+1;
        const unsigned int d1 = 1;
        const unsigned int d2 = n;

        for (unsigned int i2=0; i2<n_subdivisions; ++i2)
          for (unsigned int i1=0; i1<n_subdivisions; ++i1)
            {
              Point<spacedim> points[4];
              compute_node(points[0], &*patch, i1, i2, 0, n_subdivisions);
              compute_node(points[1], &*patch, i1+1, i2, 0, n_subdivisions);
              compute_node(points[2], &*patch, i1, i2+1, 0, n_subdivisions);
              compute_node(points[3], &*patch, i1+1, i2+1, 0, n_subdivisions);

              switch (spacedim)
                {
                case 2:
                  Assert ((flags.height_vector < patch->data.n_rows()) ||
                          patch->data.n_rows() == 0,
                          ExcIndexRange (flags.height_vector, 0,
                                         patch->data.n_rows()));
                  heights[0] = patch->data.n_rows() != 0 ?
                               patch->data(flags.height_vector,i1*d1 + i2*d2) * flags.z_scaling
                               : 0;
                  heights[1] = patch->data.n_rows() != 0 ?
                               patch->data(flags.height_vector,(i1+1)*d1 + i2*d2) * flags.z_scaling
                               : 0;
                  heights[2] = patch->data.n_rows() != 0 ?
                               patch->data(flags.height_vector,i1*d1 + (i2+1)*d2) * flags.z_scaling
                               : 0;
                  heights[3] = patch->data.n_rows() != 0 ?
                               patch->data(flags.height_vector,(i1+1)*d1 + (i2+1)*d2) * flags.z_scaling
                               : 0;

                  break;
                case 3:
                  // Copy z-coordinates into the height vector
                  for (unsigned int i=0; i<4; ++i)
                    heights[i] = points[i](2);
                  break;
                default:
                  Assert(false, ExcNotImplemented());
                }


              // now compute the projection of
              // the bilinear cell given by the
              // four vertices and their heights
              // and write them to a proper
              // cell object. note that we only
              // need the first two components
              // of the projected position for
              // output, but we need the value
              // along the line of sight for
              // sorting the cells for back-to-
              // front-output
              //
              // this computation was first written
              // by Stefan Nauber. please no-one
              // ask me why it works that way (or
              // may be not), especially not about
              // the angles and the sign of
              // the height field, I don't know
              // it.
              EpsCell2d eps_cell;
              const double pi = numbers::PI;
              const double cx = -std::cos(pi-flags.azimut_angle * 2*pi / 360.),
                           cz = -std::cos(flags.turn_angle * 2*pi / 360.),
                           sx = std::sin(pi-flags.azimut_angle * 2*pi / 360.),
                           sz = std::sin(flags.turn_angle * 2*pi / 360.);
              for (unsigned int vertex=0; vertex<4; ++vertex)
                {
                  const double x = points[vertex](0),
                               y = points[vertex](1),
                               z = -heights[vertex];

                  eps_cell.vertices[vertex](0) = -   cz*x+   sz*y;
                  eps_cell.vertices[vertex](1) = -cx*sz*x-cx*cz*y-sx*z;

                  //      ( 1 0    0 )
                  // D1 = ( 0 cx -sx )
                  //      ( 0 sx  cx )

                  //      ( cy 0 sy )
                  // Dy = (  0 1  0 )
                  //      (-sy 0 cy )

                  //      ( cz -sz 0 )
                  // Dz = ( sz  cz 0 )
                  //      (  0   0 1 )

//       ( cz -sz 0 )( 1 0    0 )(x)   ( cz*x-sz*(cx*y-sx*z)+0*(sx*y+cx*z) )
// Dxz = ( sz  cz 0 )( 0 cx -sx )(y) = ( sz*x+cz*(cx*y-sx*z)+0*(sx*y+cx*z) )
//       (  0   0 1 )( 0 sx  cx )(z)   (  0*x+  *(cx*y-sx*z)+1*(sx*y+cx*z) )
                }

              // compute coordinates of
              // center of cell
              const Point<spacedim> center_point
                = (points[0] + points[1] + points[2] + points[3]) / 4;
              const double center_height
                = -(heights[0] + heights[1] + heights[2] + heights[3]) / 4;

              // compute the depth into
              // the picture
              eps_cell.depth = -sx*sz*center_point(0)
                               -sx*cz*center_point(1)
                               +cx*center_height;

              if (flags.draw_cells && flags.shade_cells)
                {
                  Assert ((flags.color_vector < patch->data.n_rows()) ||
                          patch->data.n_rows() == 0,
                          ExcIndexRange (flags.color_vector, 0,
                                         patch->data.n_rows()));
                  const double color_values[4]
                    = { patch->data.n_rows() != 0 ?
                        patch->data(flags.color_vector,i1 *d1 + i2 *d2)       : 1,

                        patch->data.n_rows() != 0 ?
                        patch->data(flags.color_vector,(i1+1)*d1 + i2 *d2)   : 1,

                        patch->data.n_rows() != 0 ?
                        patch->data(flags.color_vector,i1 *d1 + (i2+1)*d2)     : 1,

                        patch->data.n_rows() != 0 ?
                        patch->data(flags.color_vector,(i1+1)*d1 + (i2+1)*d2) : 1
                      };

                  // set color value to average of the value
                  // at the vertices
                  eps_cell.color_value = (color_values[0] +
                                          color_values[1] +
                                          color_values[3] +
                                          color_values[2]) / 4;

                  // update bounds of color
                  // field
                  if (patch == patches.begin())
                    min_color_value = max_color_value = eps_cell.color_value;
                  else
                    {
                      min_color_value = (min_color_value < eps_cell.color_value ?
                                         min_color_value : eps_cell.color_value);
                      max_color_value = (max_color_value > eps_cell.color_value ?
                                         max_color_value : eps_cell.color_value);
                    }
                }

              // finally add this cell
              cells.insert (eps_cell);
            }
      }

    // find out minimum and maximum x and
    // y coordinates to compute offsets
    // and scaling factors
    double x_min = cells.begin()->vertices[0](0);
    double x_max = x_min;
    double y_min = cells.begin()->vertices[0](1);
    double y_max = y_min;

    for (typename std::multiset<EpsCell2d>::const_iterator
         cell=cells.begin();
         cell!=cells.end(); ++cell)
      for (unsigned int vertex=0; vertex<4; ++vertex)
        {
          x_min = std::min (x_min, cell->vertices[vertex](0));
          x_max = std::max (x_max, cell->vertices[vertex](0));
          y_min = std::min (y_min, cell->vertices[vertex](1));
          y_max = std::max (y_max, cell->vertices[vertex](1));
        }

    // scale in x-direction such that
    // in the output 0 <= x <= 300.
    // don't scale in y-direction to
    // preserve the shape of the
    // triangulation
    const double scale = (flags.size /
                          (flags.size_type==EpsFlags::width ?
                           x_max - x_min :
                           y_min - y_max));

    const Point<2> offset(x_min, y_min);


    // now write preamble
    if (true)
      {
        // block this to have local
        // variables destroyed after
        // use
        std::time_t  time1= std::time (0);
        std::tm     *time = std::localtime(&time1);
        out << "%!PS-Adobe-2.0 EPSF-1.2" << '\n'
            << "%%Title: deal.II Output" << '\n'
            << "%%Creator: the deal.II library" << '\n'
            << "%%Creation Date: "
            << time->tm_year+1900 << "/"
            << time->tm_mon+1 << "/"
            << time->tm_mday << " - "
            << time->tm_hour << ":"
            << std::setw(2) << time->tm_min << ":"
            << std::setw(2) << time->tm_sec << '\n'
            << "%%BoundingBox: "
            // lower left corner
            << "0 0 "
            // upper right corner
            << static_cast<unsigned int>( (x_max-x_min) * scale + 0.5)
            << ' '
            << static_cast<unsigned int>( (y_max-y_min) * scale + 0.5)
            << '\n';

        // define some abbreviations to keep
        // the output small:
        // m=move turtle to
        // l=define a line
        // s=set rgb color
        // sg=set gray value
        // lx=close the line and plot the line
        // lf=close the line and fill the interior
        out << "/m {moveto} bind def"      << '\n'
            << "/l {lineto} bind def"      << '\n'
            << "/s {setrgbcolor} bind def" << '\n'
            << "/sg {setgray} bind def"    << '\n'
            << "/lx {lineto closepath stroke} bind def" << '\n'
            << "/lf {lineto closepath fill} bind def"   << '\n';

        out << "%%EndProlog" << '\n'
            << '\n';
        // set fine lines
        out << flags.line_width << " setlinewidth" << '\n';
        // allow only five digits
        // for output (instead of the
        // default six); this should suffice
        // even for fine grids, but reduces
        // the file size significantly
        out << std::setprecision (5);
      }

    // check if min and max
    // values for the color are
    // actually different. If
    // that is not the case (such
    // things happen, for
    // example, in the very first
    // time step of a time
    // dependent problem, if the
    // initial values are zero),
    // all values are equal, and
    // then we can draw
    // everything in an arbitrary
    // color. Thus, change one of
    // the two values arbitrarily
    if (max_color_value == min_color_value)
      max_color_value = min_color_value+1;

    // now we've got all the information
    // we need. write the cells.
    // note: due to the ordering, we
    // traverse the list of cells
    // back-to-front
    for (typename std::multiset<EpsCell2d>::const_iterator
         cell=cells.begin();
         cell!=cells.end(); ++cell)
      {
        if (flags.draw_cells)
          {
            if (flags.shade_cells)
              {
                const EpsFlags::RgbValues rgb_values
                  = (*flags.color_function) (cell->color_value,
                                             min_color_value,
                                             max_color_value);

                // write out color
                if (rgb_values.is_grey())
                  out << rgb_values.red << " sg ";
                else
                  out << rgb_values.red   << ' '
                      << rgb_values.green << ' '
                      << rgb_values.blue  << " s ";
              }
            else
              out << "1 sg ";

            out << (cell->vertices[0]-offset) * scale << " m "
                << (cell->vertices[1]-offset) * scale << " l "
                << (cell->vertices[3]-offset) * scale << " l "
                << (cell->vertices[2]-offset) * scale << " lf"
                << '\n';
          }

        if (flags.draw_mesh)
          out << "0 sg "      // draw lines in black
              << (cell->vertices[0]-offset) * scale << " m "
              << (cell->vertices[1]-offset) * scale << " l "
              << (cell->vertices[3]-offset) * scale << " l "
              << (cell->vertices[2]-offset) * scale << " lx"
              << '\n';
      }
    out << "showpage" << '\n';
    // make sure everything now gets to
    // disk
    out << std::setprecision(old_precision);
    out.flush ();

    AssertThrow (out, ExcIO());
  }



  template <int dim, int spacedim>
  void write_gmv (const std::vector<Patch<dim,spacedim> > &patches,
                  const std::vector<std::string>          &data_names,
                  const std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> > &,
                  const GmvFlags                          &flags,
                  std::ostream                            &out)
  {
    Assert(dim<=3, ExcNotImplemented());
    AssertThrow (out, ExcIO());

#ifndef DEAL_II_WITH_MPI
    // verify that there are indeed
    // patches to be written out. most
    // of the times, people just forget
    // to call build_patches when there
    // are no patches, so a warning is
    // in order. that said, the
    // assertion is disabled if we
    // support MPI since then it can
    // happen that on the coarsest
    // mesh, a processor simply has no
    // cells it actually owns, and in
    // that case it is legit if there
    // are no patches
    Assert (patches.size() > 0, ExcNoPatches());
#else
    if (patches.size() == 0)
      return;
#endif

    GmvStream gmv_out(out, flags);
    const unsigned int n_data_sets = data_names.size();
    // check against # of data sets in
    // first patch. checks against all
    // other patches are made in
    // write_gmv_reorder_data_vectors
    Assert ((patches[0].data.n_rows() == n_data_sets && !patches[0].points_are_available) ||
            (patches[0].data.n_rows() == n_data_sets+spacedim && patches[0].points_are_available),
            ExcDimensionMismatch (patches[0].points_are_available
                                  ?
                                  (n_data_sets + spacedim)
                                  :
                                  n_data_sets,
                                  patches[0].data.n_rows()));

    ///////////////////////
    // preamble
    out << "gmvinput ascii"
        << '\n'
        << '\n';

    // first count the number of cells
    // and cells for later use
    unsigned int n_nodes;
    unsigned int n_cells;
    compute_sizes<dim,spacedim>(patches, n_nodes, n_cells);

    // in gmv format the vertex
    // coordinates and the data have an
    // order that is a bit unpleasant
    // (first all x coordinates, then
    // all y coordinate, ...; first all
    // data of variable 1, then
    // variable 2, etc), so we have to
    // copy the data vectors a bit around
    //
    // note that we copy vectors when
    // looping over the patches since we
    // have to write them one variable
    // at a time and don't want to use
    // more than one loop
    //
    // this copying of data vectors can
    // be done while we already output
    // the vertices, so do this on a
    // separate task and when wanting
    // to write out the data, we wait
    // for that task to finish
    Table<2,double> data_vectors (n_data_sets, n_nodes);
    void (*fun_ptr) (const std::vector<Patch<dim,spacedim> > &,
                     Table<2,double> &)
      = &write_gmv_reorder_data_vectors<dim,spacedim>;
    Threads::Task<> reorder_task = Threads::new_task (fun_ptr, patches, data_vectors);

    ///////////////////////////////
    // first make up a list of used
    // vertices along with their
    // coordinates
    //
    // note that we have to print
    // 3 dimensions
    out << "nodes " << n_nodes << '\n';
    for (unsigned int d=0; d<spacedim; ++d)
      {
        gmv_out.selected_component = d;
        write_nodes(patches, gmv_out);
        out << '\n';
      }
    gmv_out.selected_component = numbers::invalid_unsigned_int;

    for (unsigned int d=spacedim; d<3; ++d)
      {
        for (unsigned int i=0; i<n_nodes; ++i)
          out << "0 ";
        out << '\n';
      }

    /////////////////////////////////
    // now for the cells. note that
    // vertices are counted from 1 onwards
    out << "cells " << n_cells << '\n';
    write_cells(patches, gmv_out);

    ///////////////////////////////////////
    // data output.
    out << "variable" << '\n';

    // now write the data vectors to
    // @p{out} first make sure that all
    // data is in place
    reorder_task.join ();

    // then write data.
    // the '1' means: node data (as opposed
    // to cell data, which we do not
    // support explicitly here)
    for (unsigned int data_set=0; data_set<n_data_sets; ++data_set)
      {
        out << data_names[data_set] << " 1" << '\n';
        std::copy (data_vectors[data_set].begin(),
                   data_vectors[data_set].end(),
                   std::ostream_iterator<double>(out, " "));
        out << '\n'
            << '\n';
      }



    // end of variable section
    out << "endvars" << '\n';

    // end of output
    out << "endgmv"
        << '\n';

    // make sure everything now gets to
    // disk
    out.flush ();

    // assert the stream is still ok
    AssertThrow (out, ExcIO());
  }



  template <int dim, int spacedim>
  void write_tecplot (const std::vector<Patch<dim,spacedim> > &patches,
                      const std::vector<std::string>          &data_names,
                      const std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> > &,
                      const TecplotFlags                      &flags,
                      std::ostream                            &out)
  {
    AssertThrow (out, ExcIO());

#ifndef DEAL_II_WITH_MPI
    // verify that there are indeed
    // patches to be written out. most
    // of the times, people just forget
    // to call build_patches when there
    // are no patches, so a warning is
    // in order. that said, the
    // assertion is disabled if we
    // support MPI since then it can
    // happen that on the coarsest
    // mesh, a processor simply has no
    // cells it actually owns, and in
    // that case it is legit if there
    // are no patches
    Assert (patches.size() > 0, ExcNoPatches());
#else
    if (patches.size() == 0)
      return;
#endif

    TecplotStream tecplot_out(out, flags);

    const unsigned int n_data_sets = data_names.size();
    // check against # of data sets in
    // first patch. checks against all
    // other patches are made in
    // write_gmv_reorder_data_vectors
    Assert ((patches[0].data.n_rows() == n_data_sets && !patches[0].points_are_available) ||
            (patches[0].data.n_rows() == n_data_sets+spacedim && patches[0].points_are_available),
            ExcDimensionMismatch (patches[0].points_are_available
                                  ?
                                  (n_data_sets + spacedim)
                                  :
                                  n_data_sets,
                                  patches[0].data.n_rows()));

    // first count the number of cells
    // and cells for later use
    unsigned int n_nodes;
    unsigned int n_cells;
    compute_sizes<dim,spacedim>(patches, n_nodes, n_cells);

    ///////////
    // preamble
    {
      std::time_t  time1= std::time (0);
      std::tm     *time = std::localtime(&time1);
      out << "# This file was generated by the deal.II library." << '\n'
          << "# Date =  "
          << time->tm_year+1900 << "/"
          << time->tm_mon+1 << "/"
          << time->tm_mday << '\n'
          << "# Time =  "
          << time->tm_hour << ":"
          << std::setw(2) << time->tm_min << ":"
          << std::setw(2) << time->tm_sec << '\n'
          << "#" << '\n'
          << "# For a description of the Tecplot format see the Tecplot documentation."
          << '\n'
          << "#" << '\n';


      out << "Variables=";

      switch (spacedim)
        {
        case 1:
          out << "\"x\"";
          break;
        case 2:
          out << "\"x\", \"y\"";
          break;
        case 3:
          out << "\"x\", \"y\", \"z\"";
          break;
        default:
          Assert (false, ExcNotImplemented());
        }

      for (unsigned int data_set=0; data_set<n_data_sets; ++data_set)
        out << ", \"" << data_names[data_set] << "\"";

      out << '\n';

      out << "zone ";
      if (flags.zone_name)
        out << "t=\"" << flags.zone_name << "\" ";

      out << "f=feblock, n=" << n_nodes << ", e=" << n_cells
          << ", et=" << tecplot_cell_type[dim] << '\n';
    }


    // in Tecplot FEBLOCK format the vertex
    // coordinates and the data have an
    // order that is a bit unpleasant
    // (first all x coordinates, then
    // all y coordinate, ...; first all
    // data of variable 1, then
    // variable 2, etc), so we have to
    // copy the data vectors a bit around
    //
    // note that we copy vectors when
    // looping over the patches since we
    // have to write them one variable
    // at a time and don't want to use
    // more than one loop
    //
    // this copying of data vectors can
    // be done while we already output
    // the vertices, so do this on a
    // separate task and when wanting
    // to write out the data, we wait
    // for that task to finish

    Table<2,double> data_vectors (n_data_sets, n_nodes);

    void (*fun_ptr) (const std::vector<Patch<dim,spacedim> > &,
                     Table<2,double> &)
      = &write_gmv_reorder_data_vectors<dim,spacedim>;
    Threads::Task<> reorder_task = Threads::new_task (fun_ptr, patches, data_vectors);

    ///////////////////////////////
    // first make up a list of used
    // vertices along with their
    // coordinates


    for (unsigned int d=0; d<spacedim; ++d)
      {
        tecplot_out.selected_component = d;
        write_nodes(patches, tecplot_out);
        out << '\n';
      }


    ///////////////////////////////////////
    // data output.
    //
    // now write the data vectors to
    // @p{out} first make sure that all
    // data is in place
    reorder_task.join ();

    // then write data.
    for (unsigned int data_set=0; data_set<n_data_sets; ++data_set)
      {
        std::copy (data_vectors[data_set].begin(),
                   data_vectors[data_set].end(),
                   std::ostream_iterator<double>(out, "\n"));
        out << '\n';
      }

    write_cells(patches, tecplot_out);

    // make sure everything now gets to
    // disk
    out.flush ();

    // assert the stream is still ok
    AssertThrow (out, ExcIO());
  }



//---------------------------------------------------------------------------
// Macros for handling Tecplot API data

#ifdef DEAL_II_HAVE_TECPLOT

  namespace
  {
    class TecplotMacros
    {
    public:
      TecplotMacros(const unsigned int n_nodes = 0,
                    const unsigned int n_vars = 0,
                    const unsigned int n_cells = 0,
                    const unsigned int n_vert = 0);
      ~TecplotMacros();
      float &nd(const unsigned int i, const unsigned int j);
      int    &cd(const unsigned int i, const unsigned int j);
      std::vector<float> nodalData;
      std::vector<int>   connData;
    private:
      unsigned int n_nodes;
      unsigned int n_vars;
      unsigned int n_cells;
      unsigned int n_vert;
    };


    inline
    TecplotMacros::TecplotMacros(const unsigned int n_nodes,
                                 const unsigned int n_vars,
                                 const unsigned int n_cells,
                                 const unsigned int n_vert)
      :
      n_nodes(n_nodes),
      n_vars(n_vars),
      n_cells(n_cells),
      n_vert(n_vert)
    {
      nodalData.resize(n_nodes*n_vars);
      connData.resize(n_cells*n_vert);
    }



    inline
    TecplotMacros::~TecplotMacros()
    {}



    inline
    float &TecplotMacros::nd (const unsigned int i,
                              const unsigned int j)
    {
      return nodalData[i*n_nodes+j];
    }



    inline
    int &TecplotMacros::cd (const unsigned int i,
                            const unsigned int j)
    {
      return connData[i+j*n_vert];
    }

  }


#endif
//---------------------------------------------------------------------------



  template <int dim, int spacedim>
  void write_tecplot_binary (const std::vector<Patch<dim,spacedim> > &patches,
                             const std::vector<std::string>          &data_names,
                             const std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> > &vector_data_ranges,
                             const TecplotFlags                      &flags,
                             std::ostream                            &out)
  {

#ifndef DEAL_II_HAVE_TECPLOT

    // simply call the ASCII output
    // function if the Tecplot API
    // isn't present
    write_tecplot (patches, data_names, vector_data_ranges, flags, out);
    return;

#else

    // Tecplot binary output only good
    // for 2D & 3D
    if (dim == 1)
      {
        write_tecplot (patches, data_names, vector_data_ranges, flags, out);
        return;
      }

    // if the user hasn't specified a
    // file name we should call the
    // ASCII function and use the
    // ostream @p{out} instead of doing
    // something silly later
    char *file_name = (char *) flags.tecplot_binary_file_name;

    if (file_name == NULL)
      {
        // At least in debug mode we
        // should tell users why they
        // don't get tecplot binary
        // output
        Assert(false, ExcMessage("Specify the name of the tecplot_binary"
                                 " file through the TecplotFlags interface."));
        write_tecplot (patches, data_names, vector_data_ranges, flags, out);
        return;
      }


    AssertThrow (out, ExcIO());

#ifndef DEAL_II_WITH_MPI
    // verify that there are indeed
    // patches to be written out. most
    // of the times, people just forget
    // to call build_patches when there
    // are no patches, so a warning is
    // in order. that said, the
    // assertion is disabled if we
    // support MPI since then it can
    // happen that on the coarsest
    // mesh, a processor simply has no
    // cells it actually owns, and in
    // that case it is legit if there
    // are no patches
    Assert (patches.size() > 0, ExcNoPatches());
#else
    if (patches.size() == 0)
      return;
#endif

    const unsigned int n_data_sets = data_names.size();
    // check against # of data sets in
    // first patch. checks against all
    // other patches are made in
    // write_gmv_reorder_data_vectors
    Assert ((patches[0].data.n_rows() == n_data_sets && !patches[0].points_are_available) ||
            (patches[0].data.n_rows() == n_data_sets+spacedim && patches[0].points_are_available),
            ExcDimensionMismatch (patches[0].points_are_available
                                  ?
                                  (n_data_sets + spacedim)
                                  :
                                  n_data_sets,
                                  patches[0].data.n_rows()));

    // first count the number of cells
    // and cells for later use
    unsigned int n_nodes;
    unsigned int n_cells;
    compute_sizes<dim,spacedim>(patches, n_nodes, n_cells);
    // local variables only needed to write Tecplot
    // binary output files
    const unsigned int vars_per_node  = (spacedim+n_data_sets),
                       nodes_per_cell = GeometryInfo<dim>::vertices_per_cell;

    TecplotMacros tm(n_nodes, vars_per_node, n_cells, nodes_per_cell);

    int is_double = 0,
        tec_debug = 0,
        cell_type = tecplot_binary_cell_type[dim];

    std::string tec_var_names;
    switch (spacedim)
      {
      case 2:
        tec_var_names  = "x y";
        break;
      case 3:
        tec_var_names  = "x y z";
        break;
      default:
        Assert(false, ExcNotImplemented());
      }

    for (unsigned int data_set=0; data_set<n_data_sets; ++data_set)
      {
        tec_var_names += " ";
        tec_var_names += data_names[data_set];
      }
    // in Tecplot FEBLOCK format the vertex
    // coordinates and the data have an
    // order that is a bit unpleasant
    // (first all x coordinates, then
    // all y coordinate, ...; first all
    // data of variable 1, then
    // variable 2, etc), so we have to
    // copy the data vectors a bit around
    //
    // note that we copy vectors when
    // looping over the patches since we
    // have to write them one variable
    // at a time and don't want to use
    // more than one loop
    //
    // this copying of data vectors can
    // be done while we already output
    // the vertices, so do this on a
    // separate task and when wanting
    // to write out the data, we wait
    // for that task to finish
    Table<2,double> data_vectors (n_data_sets, n_nodes);

    void (*fun_ptr) (const std::vector<Patch<dim,spacedim> > &,
                     Table<2,double> &)
      = &write_gmv_reorder_data_vectors<dim,spacedim>;
    Threads::Task<> reorder_task = Threads::new_task (fun_ptr, patches, data_vectors);

    ///////////////////////////////
    // first make up a list of used
    // vertices along with their
    // coordinates
    for (unsigned int d=1; d<=spacedim; ++d)
      {
        unsigned int entry=0;

        for (typename std::vector<Patch<dim,spacedim> >::const_iterator patch=patches.begin();
             patch!=patches.end(); ++patch)
          {
            const unsigned int n_subdivisions = patch->n_subdivisions;

            switch (dim)
              {
              case 2:
              {
                for (unsigned int j=0; j<n_subdivisions+1; ++j)
                  for (unsigned int i=0; i<n_subdivisions+1; ++i)
                    {
                      const double x_frac = i * 1./n_subdivisions,
                                   y_frac = j * 1./n_subdivisions;

                      tm.nd((d-1),entry) = static_cast<float>(
                                             (((patch->vertices[1](d-1) * x_frac) +
                                               (patch->vertices[0](d-1) * (1-x_frac))) * (1-y_frac) +
                                              ((patch->vertices[3](d-1) * x_frac) +
                                               (patch->vertices[2](d-1) * (1-x_frac))) * y_frac)
                                           );
                      entry++;
                    }
                break;
              }

              case 3:
              {
                for (unsigned int j=0; j<n_subdivisions+1; ++j)
                  for (unsigned int k=0; k<n_subdivisions+1; ++k)
                    for (unsigned int i=0; i<n_subdivisions+1; ++i)
                      {
                        const double x_frac = i * 1./n_subdivisions,
                                     y_frac = k * 1./n_subdivisions,
                                     z_frac = j * 1./n_subdivisions;

                        // compute coordinates for
                        // this patch point
                        tm.nd((d-1),entry) = static_cast<float>(
                                               ((((patch->vertices[1](d-1) * x_frac) +
                                                  (patch->vertices[0](d-1) * (1-x_frac))) * (1-y_frac) +
                                                 ((patch->vertices[3](d-1) * x_frac) +
                                                  (patch->vertices[2](d-1) * (1-x_frac))) * y_frac)   * (1-z_frac) +
                                                (((patch->vertices[5](d-1) * x_frac) +
                                                  (patch->vertices[4](d-1) * (1-x_frac))) * (1-y_frac) +
                                                 ((patch->vertices[7](d-1) * x_frac) +
                                                  (patch->vertices[6](d-1) * (1-x_frac))) * y_frac)   * z_frac)
                                             );
                        entry++;
                      }
                break;
              }

              default:
                Assert (false, ExcNotImplemented());
              }
          }
      }


    ///////////////////////////////////////
    // data output.
    //
    reorder_task.join ();

    // then write data.
    for (unsigned int data_set=0; data_set<n_data_sets; ++data_set)
      for (unsigned int entry=0; entry<data_vectors[data_set].size(); entry++)
        tm.nd((spacedim+data_set),entry) = static_cast<float>(data_vectors[data_set][entry]);




    /////////////////////////////////
    // now for the cells. note that
    // vertices are counted from 1 onwards
    unsigned int first_vertex_of_patch = 0;
    unsigned int elem=0;

    for (typename std::vector<Patch<dim,spacedim> >::const_iterator patch=patches.begin();
         patch!=patches.end(); ++patch)
      {
        const unsigned int n_subdivisions = patch->n_subdivisions;
        const unsigned int n = n_subdivisions+1;
        const unsigned int d1=1;
        const unsigned int d2=n;
        const unsigned int d3=n*n;
        // write out the cells making
        // up this patch
        switch (dim)
          {
          case 2:
          {
            for (unsigned int i2=0; i2<n_subdivisions; ++i2)
              for (unsigned int i1=0; i1<n_subdivisions; ++i1)
                {
                  tm.cd(0,elem) = first_vertex_of_patch+(i1  )*d1+(i2  )*d2+1;
                  tm.cd(1,elem) = first_vertex_of_patch+(i1+1)*d1+(i2  )*d2+1;
                  tm.cd(2,elem) = first_vertex_of_patch+(i1+1)*d1+(i2+1)*d2+1;
                  tm.cd(3,elem) = first_vertex_of_patch+(i1  )*d1+(i2+1)*d2+1;

                  elem++;
                }
            break;
          }

          case 3:
          {
            for (unsigned int i3=0; i3<n_subdivisions; ++i3)
              for (unsigned int i2=0; i2<n_subdivisions; ++i2)
                for (unsigned int i1=0; i1<n_subdivisions; ++i1)
                  {
                    // note: vertex indices start with 1!


                    tm.cd(0,elem) = first_vertex_of_patch+(i1  )*d1+(i2  )*d2+(i3  )*d3+1;
                    tm.cd(1,elem) = first_vertex_of_patch+(i1+1)*d1+(i2  )*d2+(i3  )*d3+1;
                    tm.cd(2,elem) = first_vertex_of_patch+(i1+1)*d1+(i2+1)*d2+(i3  )*d3+1;
                    tm.cd(3,elem) = first_vertex_of_patch+(i1  )*d1+(i2+1)*d2+(i3  )*d3+1;
                    tm.cd(4,elem) = first_vertex_of_patch+(i1  )*d1+(i2  )*d2+(i3+1)*d3+1;
                    tm.cd(5,elem) = first_vertex_of_patch+(i1+1)*d1+(i2  )*d2+(i3+1)*d3+1;
                    tm.cd(6,elem) = first_vertex_of_patch+(i1+1)*d1+(i2+1)*d2+(i3+1)*d3+1;
                    tm.cd(7,elem) = first_vertex_of_patch+(i1  )*d1+(i2+1)*d2+(i3+1)*d3+1;

                    elem++;
                  }
            break;
          }

          default:
            Assert (false, ExcNotImplemented());
          }


        // finally update the number
        // of the first vertex of this patch
        first_vertex_of_patch += Utilities::fixed_power<dim>(n);
      }


    {
      int ierr      = 0,
          num_nodes = static_cast<int>(n_nodes),
          num_cells = static_cast<int>(n_cells);

      char dot[2] = {'.', 0};
      // Unfortunately, TECINI takes a
      // char *, but c_str() gives a
      // const char *.  As we don't do
      // anything else with
      // tec_var_names following
      // const_cast is ok
      char *var_names=const_cast<char *> (tec_var_names.c_str());
      ierr = TECINI (NULL,
                     var_names,
                     file_name,
                     dot,
                     &tec_debug,
                     &is_double);

      Assert (ierr == 0, ExcErrorOpeningTecplotFile(file_name));

      char FEBLOCK[] = {'F','E','B','L','O','C','K',0};
      ierr = TECZNE (NULL,
                     &num_nodes,
                     &num_cells,
                     &cell_type,
                     FEBLOCK,
                     NULL);

      Assert (ierr == 0, ExcTecplotAPIError());

      int total = (vars_per_node*num_nodes);

      ierr = TECDAT (&total,
                     &tm.nodalData[0],
                     &is_double);

      Assert (ierr == 0, ExcTecplotAPIError());

      ierr = TECNOD (&tm.connData[0]);

      Assert (ierr == 0, ExcTecplotAPIError());

      ierr = TECEND ();

      Assert (ierr == 0, ExcTecplotAPIError());
    }
#endif
  }



  template <int dim, int spacedim>
  void
  write_vtk (const std::vector<Patch<dim,spacedim> > &patches,
             const std::vector<std::string>          &data_names,
             const std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> > &vector_data_ranges,
             const VtkFlags                          &flags,
             std::ostream                            &out)
  {
    AssertThrow (out, ExcIO());

#ifndef DEAL_II_WITH_MPI
    // verify that there are indeed
    // patches to be written out. most
    // of the times, people just forget
    // to call build_patches when there
    // are no patches, so a warning is
    // in order. that said, the
    // assertion is disabled if we
    // support MPI since then it can
    // happen that on the coarsest
    // mesh, a processor simply has no
    // cells it actually owns, and in
    // that case it is legit if there
    // are no patches
    Assert (patches.size() > 0, ExcNoPatches());
#else
    if (patches.size() == 0)
      return;
#endif

    VtkStream vtk_out(out, flags);

    const unsigned int n_data_sets = data_names.size();
    // check against # of data sets in
    // first patch. checks against all
    // other patches are made in
    // write_gmv_reorder_data_vectors
    Assert ((patches[0].data.n_rows() == n_data_sets && !patches[0].points_are_available) ||
            (patches[0].data.n_rows() == n_data_sets+spacedim && patches[0].points_are_available),
            ExcDimensionMismatch (patches[0].points_are_available
                                  ?
                                  (n_data_sets + spacedim)
                                  :
                                  n_data_sets,
                                  patches[0].data.n_rows()));

    ///////////////////////
    // preamble
    {
      std::time_t  time1= std::time (0);
      std::tm     *time = std::localtime(&time1);
      out << "# vtk DataFile Version 3.0"
          << '\n'
          << "#This file was generated by the deal.II library";
      if (flags.print_date_and_time)
        out << " on "
            << time->tm_year+1900 << "/"
            << time->tm_mon+1 << "/"
            << time->tm_mday << " at "
            << time->tm_hour << ":"
            << std::setw(2) << time->tm_min << ":"
            << std::setw(2) << time->tm_sec;
      else
        out << ".";
      out << '\n'
          << "ASCII"
          << '\n';
      // now output the data header
      out << "DATASET UNSTRUCTURED_GRID\n"
          << '\n';
    }

    // if desired, output time and cycle of the simulation, following
    // the instructions at
    // http://www.visitusers.org/index.php?title=Time_and_Cycle_in_VTK_files
    {
      const unsigned int
      n_metadata = ((flags.cycle != std::numeric_limits<unsigned int>::min() ? 1 : 0)
                    +
                    (flags.time != std::numeric_limits<double>::min() ? 1 : 0));
      if (n_metadata > 0)
        out << "FIELD FieldData " << n_metadata << "\n";

      if (flags.cycle != std::numeric_limits<unsigned int>::min())
        {
          out << "CYCLE 1 1 int\n"
              << flags.cycle << "\n";
        }
      if (flags.time != std::numeric_limits<double>::min())
        {
          out << "TIME 1 1 double\n"
              << flags.time << "\n";
        }
    }

    // first count the number of cells
    // and cells for later use
    unsigned int n_nodes;
    unsigned int n_cells;
    compute_sizes<dim,spacedim>(patches, n_nodes, n_cells);
    // in gmv format the vertex
    // coordinates and the data have an
    // order that is a bit unpleasant
    // (first all x coordinates, then
    // all y coordinate, ...; first all
    // data of variable 1, then
    // variable 2, etc), so we have to
    // copy the data vectors a bit around
    //
    // note that we copy vectors when
    // looping over the patches since we
    // have to write them one variable
    // at a time and don't want to use
    // more than one loop
    //
    // this copying of data vectors can
    // be done while we already output
    // the vertices, so do this on a
    // separate task and when wanting
    // to write out the data, we wait
    // for that task to finish
    Table<2,double> data_vectors (n_data_sets, n_nodes);

    void (*fun_ptr) (const std::vector<Patch<dim,spacedim> > &,
                     Table<2,double> &)
      = &write_gmv_reorder_data_vectors<dim,spacedim>;
    Threads::Task<> reorder_task = Threads::new_task (fun_ptr, patches, data_vectors);

    ///////////////////////////////
    // first make up a list of used
    // vertices along with their
    // coordinates
    //
    // note that we have to print
    // d=1..3 dimensions
    out << "POINTS " << n_nodes << " double" << '\n';
    write_nodes(patches, vtk_out);
    out << '\n';
    /////////////////////////////////
    // now for the cells
    out << "CELLS " << n_cells << ' '
        << n_cells*(GeometryInfo<dim>::vertices_per_cell+1)
        << '\n';
    write_cells(patches, vtk_out);
    out << '\n';
    // next output the types of the
    // cells. since all cells are
    // the same, this is simple
    out << "CELL_TYPES " << n_cells << '\n';
    for (unsigned int i=0; i<n_cells; ++i)
      out << ' ' << vtk_cell_type[dim];
    out << '\n';
    ///////////////////////////////////////
    // data output.

    // now write the data vectors to
    // @p{out} first make sure that all
    // data is in place
    reorder_task.join ();

    // then write data.  the
    // 'POINT_DATA' means: node data
    // (as opposed to cell data, which
    // we do not support explicitly
    // here). all following data sets
    // are point data
    out << "POINT_DATA " << n_nodes
        << '\n';

    // when writing, first write out
    // all vector data, then handle the
    // scalar data sets that have been
    // left over
    std::vector<bool> data_set_written (n_data_sets, false);
    for (unsigned int n_th_vector=0; n_th_vector<vector_data_ranges.size(); ++n_th_vector)
      {
        AssertThrow (std_cxx11::get<1>(vector_data_ranges[n_th_vector]) >=
                     std_cxx11::get<0>(vector_data_ranges[n_th_vector]),
                     ExcLowerRange (std_cxx11::get<1>(vector_data_ranges[n_th_vector]),
                                    std_cxx11::get<0>(vector_data_ranges[n_th_vector])));
        AssertThrow (std_cxx11::get<1>(vector_data_ranges[n_th_vector]) < n_data_sets,
                     ExcIndexRange (std_cxx11::get<1>(vector_data_ranges[n_th_vector]),
                                    0, n_data_sets));
        AssertThrow (std_cxx11::get<1>(vector_data_ranges[n_th_vector]) + 1
                     - std_cxx11::get<0>(vector_data_ranges[n_th_vector]) <= 3,
                     ExcMessage ("Can't declare a vector with more than 3 components "
                                 "in VTK"));

        // mark these components as already written:
        for (unsigned int i=std_cxx11::get<0>(vector_data_ranges[n_th_vector]);
             i<=std_cxx11::get<1>(vector_data_ranges[n_th_vector]);
             ++i)
          data_set_written[i] = true;

        // write the
        // header. concatenate all the
        // component names with double
        // underscores unless a vector
        // name has been specified
        out << "VECTORS ";

        if (std_cxx11::get<2>(vector_data_ranges[n_th_vector]) != "")
          out << std_cxx11::get<2>(vector_data_ranges[n_th_vector]);
        else
          {
            for (unsigned int i=std_cxx11::get<0>(vector_data_ranges[n_th_vector]);
                 i<std_cxx11::get<1>(vector_data_ranges[n_th_vector]);
                 ++i)
              out << data_names[i] << "__";
            out << data_names[std_cxx11::get<1>(vector_data_ranges[n_th_vector])];
          }

        out << " double"
            << '\n';

        // now write data. pad all
        // vectors to have three
        // components
        for (unsigned int n=0; n<n_nodes; ++n)
          {
            switch (std_cxx11::get<1>(vector_data_ranges[n_th_vector]) -
                    std_cxx11::get<0>(vector_data_ranges[n_th_vector]))
              {
              case 0:
                out << data_vectors(std_cxx11::get<0>(vector_data_ranges[n_th_vector]), n) << " 0 0"
                    << '\n';
                break;

              case 1:
                out << data_vectors(std_cxx11::get<0>(vector_data_ranges[n_th_vector]),   n) << ' '<< data_vectors(std_cxx11::get<0>(vector_data_ranges[n_th_vector])+1, n) << " 0"
                    << '\n';
                break;
              case 2:
                out << data_vectors(std_cxx11::get<0>(vector_data_ranges[n_th_vector]),   n) << ' '<< data_vectors(std_cxx11::get<0>(vector_data_ranges[n_th_vector])+1, n) << ' '<< data_vectors(std_cxx11::get<0>(vector_data_ranges[n_th_vector])+2, n)
                    << '\n';
                break;

              default:
                // VTK doesn't
                // support
                // anything else
                // than vectors
                // with 1, 2, or
                // 3 components
                Assert (false, ExcInternalError());
              }
          }
      }

    // now do the left over scalar data sets
    for (unsigned int data_set=0; data_set<n_data_sets; ++data_set)
      if (data_set_written[data_set] == false)
        {
          out << "SCALARS "
              << data_names[data_set]
              << " double 1"
              << '\n'
              << "LOOKUP_TABLE default"
              << '\n';
          std::copy (data_vectors[data_set].begin(),
                     data_vectors[data_set].end(),
                     std::ostream_iterator<double>(out, " "));
          out << '\n';
        }

    // make sure everything now gets to
    // disk
    out.flush ();

    // assert the stream is still ok
    AssertThrow (out, ExcIO());
  }


  void write_vtu_header (std::ostream &out,
                         const VtkFlags &flags)
  {
    AssertThrow (out, ExcIO());
    std::time_t  time1= std::time (0);
    std::tm     *time = std::localtime(&time1);
    out << "<?xml version=\"1.0\" ?> \n";
    out << "<!-- \n";
    out << "# vtk DataFile Version 3.0"
        << '\n'
        << "#This file was generated by the deal.II library";
    if (flags.print_date_and_time)
      out << " on "
          << time->tm_year+1900 << "/"
          << time->tm_mon+1 << "/"
          << time->tm_mday << " at "
          << time->tm_hour << ":"
          << std::setw(2) << time->tm_min << ":"
          << std::setw(2) << time->tm_sec;
    else
      out << ".";
    out << "\n-->\n";
    out << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\"";
#ifdef DEAL_II_WITH_ZLIB
    out << " compressor=\"vtkZLibDataCompressor\"";
#endif
#ifdef DEAL_II_WORDS_BIGENDIAN
    out << " byte_order=\"BigEndian\"";
#else
    out << " byte_order=\"LittleEndian\"";
#endif
    out << ">";
    out << '\n';
    out << "<UnstructuredGrid>";
    out << '\n';
  }



  void write_vtu_footer (std::ostream &out)
  {
    AssertThrow (out, ExcIO());
    out << " </UnstructuredGrid>\n";
    out << "</VTKFile>\n";
  }



  template <int dim, int spacedim>
  void
  write_vtu (const std::vector<Patch<dim,spacedim> > &patches,
             const std::vector<std::string>          &data_names,
             const std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> > &vector_data_ranges,
             const VtkFlags                          &flags,
             std::ostream                            &out)
  {
    write_vtu_header(out, flags);
    write_vtu_main (patches, data_names, vector_data_ranges, flags, out);
    write_vtu_footer(out);

    out << std::flush;
  }


  template <int dim, int spacedim>
  void write_vtu_main (const std::vector<Patch<dim,spacedim> > &patches,
                       const std::vector<std::string>          &data_names,
                       const std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> > &vector_data_ranges,
                       const VtkFlags                          &flags,
                       std::ostream                            &out)
  {
    AssertThrow (out, ExcIO());

#ifndef DEAL_II_WITH_MPI
    // verify that there are indeed
    // patches to be written out. most
    // of the times, people just forget
    // to call build_patches when there
    // are no patches, so a warning is
    // in order. that said, the
    // assertion is disabled if we
    // support MPI since then it can
    // happen that on the coarsest
    // mesh, a processor simply has no
    // cells it actually owns, and in
    // that case it is legit if there
    // are no patches
    Assert (patches.size() > 0, ExcNoPatches());
#else
    if (patches.size() == 0)
      {
        // we still need to output a valid vtu file, because other CPUs
        // might output data. This is the minimal file that is accepted by paraview and visit.
        // if we remove the field definitions, visit is complaining.
        out << "<Piece NumberOfPoints=\"0\" NumberOfCells=\"0\" >\n"
            << "<Cells>\n"
            << "<DataArray type=\"UInt8\" Name=\"types\"></DataArray>\n"
            << "</Cells>\n"
            << "  <PointData Scalars=\"scalars\">\n";
        std::vector<bool> data_set_written (data_names.size(), false);
        for (unsigned int n_th_vector=0; n_th_vector<vector_data_ranges.size(); ++n_th_vector)
          {
            // mark these components as already
            // written:
            for (unsigned int i=std_cxx11::get<0>(vector_data_ranges[n_th_vector]);
                 i<=std_cxx11::get<1>(vector_data_ranges[n_th_vector]);
                 ++i)
              data_set_written[i] = true;

            // write the
            // header. concatenate all the
            // component names with double
            // underscores unless a vector
            // name has been specified
            out << "    <DataArray type=\"Float64\" Name=\"";

            if (std_cxx11::get<2>(vector_data_ranges[n_th_vector]) != "")
              out << std_cxx11::get<2>(vector_data_ranges[n_th_vector]);
            else
              {
                for (unsigned int i=std_cxx11::get<0>(vector_data_ranges[n_th_vector]);
                     i<std_cxx11::get<1>(vector_data_ranges[n_th_vector]);
                     ++i)
                  out << data_names[i] << "__";
                out << data_names[std_cxx11::get<1>(vector_data_ranges[n_th_vector])];
              }

            out << "\" NumberOfComponents=\"3\"></DataArray>\n";
          }

        for (unsigned int data_set=0; data_set<data_names.size(); ++data_set)
          if (data_set_written[data_set] == false)
            {
              out << "    <DataArray type=\"Float64\" Name=\""
                  << data_names[data_set]
                  << "\"></DataArray>\n";
            }

        out << "  </PointData>\n";
        out << "</Piece>\n";

        out << std::flush;

        return;
      }
#endif

    // first up: metadata
    //
    // if desired, output time and cycle of the simulation, following
    // the instructions at
    // http://www.visitusers.org/index.php?title=Time_and_Cycle_in_VTK_files
    {
      const unsigned int
      n_metadata = ((flags.cycle != std::numeric_limits<unsigned int>::min() ? 1 : 0)
                    +
                    (flags.time != std::numeric_limits<double>::min() ? 1 : 0));
      if (n_metadata > 0)
        out << "<FieldData>\n";

      if (flags.cycle != std::numeric_limits<unsigned int>::min())
        {
          out << "<DataArray type=\"Float32\" Name=\"CYCLE\" NumberOfTuples=\"1\" format=\"ascii\">"
              << flags.cycle
              << "</DataArray>\n";
        }
      if (flags.time != std::numeric_limits<double>::min())
        {
          out << "<DataArray type=\"Float32\" Name=\"TIME\" NumberOfTuples=\"1\" format=\"ascii\">"
              << flags.time
              << "</DataArray>\n";
        }

      if (n_metadata > 0)
        out << "</FieldData>\n";
    }


    VtuStream vtu_out(out, flags);

    const unsigned int n_data_sets = data_names.size();
    // check against # of data sets in
    // first patch. checks against all
    // other patches are made in
    // write_gmv_reorder_data_vectors
    Assert ((patches[0].data.n_rows() == n_data_sets && !patches[0].points_are_available) ||
            (patches[0].data.n_rows() == n_data_sets+spacedim && patches[0].points_are_available),
            ExcDimensionMismatch (patches[0].points_are_available
                                  ?
                                  (n_data_sets + spacedim)
                                  :
                                  n_data_sets,
                                  patches[0].data.n_rows()));


#ifdef DEAL_II_WITH_ZLIB
    const char *ascii_or_binary = "binary";
#else
    const char *ascii_or_binary = "ascii";
#endif


    // first count the number of cells
    // and cells for later use
    unsigned int n_nodes;
    unsigned int n_cells;
    compute_sizes<dim,spacedim>(patches, n_nodes, n_cells);
    // in gmv format the vertex
    // coordinates and the data have an
    // order that is a bit unpleasant
    // (first all x coordinates, then
    // all y coordinate, ...; first all
    // data of variable 1, then
    // variable 2, etc), so we have to
    // copy the data vectors a bit around
    //
    // note that we copy vectors when
    // looping over the patches since we
    // have to write them one variable
    // at a time and don't want to use
    // more than one loop
    //
    // this copying of data vectors can
    // be done while we already output
    // the vertices, so do this on a
    // separate task and when wanting
    // to write out the data, we wait
    // for that task to finish
    Table<2,double> data_vectors (n_data_sets, n_nodes);

    void (*fun_ptr) (const std::vector<Patch<dim,spacedim> > &,
                     Table<2,double> &)
      = &write_gmv_reorder_data_vectors<dim,spacedim>;
    Threads::Task<> reorder_task = Threads::new_task (fun_ptr, patches,
                                                      data_vectors);

    ///////////////////////////////
    // first make up a list of used
    // vertices along with their
    // coordinates
    //
    // note that according to the standard, we
    // have to print d=1..3 dimensions, even if
    // we are in reality in 2d, for example
    out << "<Piece NumberOfPoints=\"" << n_nodes
        <<"\" NumberOfCells=\"" << n_cells << "\" >\n";
    out << "  <Points>\n";
    out << "    <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\""
        << ascii_or_binary << "\">\n";
    write_nodes(patches, vtu_out);
    out << "    </DataArray>\n";
    out << "  </Points>\n\n";
    /////////////////////////////////
    // now for the cells
    out << "  <Cells>\n";
    out << "    <DataArray type=\"Int32\" Name=\"connectivity\" format=\""
        << ascii_or_binary << "\">\n";
    write_cells(patches, vtu_out);
    out << "    </DataArray>\n";

    // XML VTU format uses offsets; this is
    // different than the VTK format, which
    // puts the number of nodes per cell in
    // front of the connectivity list.
    out << "    <DataArray type=\"Int32\" Name=\"offsets\" format=\""
        << ascii_or_binary << "\">\n";

    std::vector<int32_t> offsets (n_cells);
    for (unsigned int i=0; i<n_cells; ++i)
      offsets[i] = (i+1)*GeometryInfo<dim>::vertices_per_cell;
    vtu_out << offsets;
    out << "\n";
    out << "    </DataArray>\n";

    // next output the types of the
    // cells. since all cells are
    // the same, this is simple
    out << "    <DataArray type=\"UInt8\" Name=\"types\" format=\""
        << ascii_or_binary << "\">\n";

    {
      // uint8_t might be a typedef to unsigned
      // char which is then not printed as
      // ascii integers
#ifdef DEAL_II_WITH_ZLIB
      std::vector<uint8_t> cell_types (n_cells,
                                       static_cast<uint8_t>(vtk_cell_type[dim]));
#else
      std::vector<unsigned int> cell_types (n_cells,
                                            vtk_cell_type[dim]);
#endif
      // this should compress well :-)
      vtu_out << cell_types;
    }
    out << "\n";
    out << "    </DataArray>\n";
    out << "  </Cells>\n";


    ///////////////////////////////////////
    // data output.

    // now write the data vectors to
    // @p{out} first make sure that all
    // data is in place
    reorder_task.join ();

    // then write data.  the
    // 'POINT_DATA' means: node data
    // (as opposed to cell data, which
    // we do not support explicitly
    // here). all following data sets
    // are point data
    out << "  <PointData Scalars=\"scalars\">\n";

    // when writing, first write out
    // all vector data, then handle the
    // scalar data sets that have been
    // left over
    std::vector<bool> data_set_written (n_data_sets, false);
    for (unsigned int n_th_vector=0; n_th_vector<vector_data_ranges.size(); ++n_th_vector)
      {
        AssertThrow (std_cxx11::get<1>(vector_data_ranges[n_th_vector]) >=
                     std_cxx11::get<0>(vector_data_ranges[n_th_vector]),
                     ExcLowerRange (std_cxx11::get<1>(vector_data_ranges[n_th_vector]),
                                    std_cxx11::get<0>(vector_data_ranges[n_th_vector])));
        AssertThrow (std_cxx11::get<1>(vector_data_ranges[n_th_vector]) < n_data_sets,
                     ExcIndexRange (std_cxx11::get<1>(vector_data_ranges[n_th_vector]),
                                    0, n_data_sets));
        AssertThrow (std_cxx11::get<1>(vector_data_ranges[n_th_vector]) + 1
                     - std_cxx11::get<0>(vector_data_ranges[n_th_vector]) <= 3,
                     ExcMessage ("Can't declare a vector with more than 3 components "
                                 "in VTK"));

        // mark these components as already
        // written:
        for (unsigned int i=std_cxx11::get<0>(vector_data_ranges[n_th_vector]);
             i<=std_cxx11::get<1>(vector_data_ranges[n_th_vector]);
             ++i)
          data_set_written[i] = true;

        // write the
        // header. concatenate all the
        // component names with double
        // underscores unless a vector
        // name has been specified
        out << "    <DataArray type=\"Float64\" Name=\"";

        if (std_cxx11::get<2>(vector_data_ranges[n_th_vector]) != "")
          out << std_cxx11::get<2>(vector_data_ranges[n_th_vector]);
        else
          {
            for (unsigned int i=std_cxx11::get<0>(vector_data_ranges[n_th_vector]);
                 i<std_cxx11::get<1>(vector_data_ranges[n_th_vector]);
                 ++i)
              out << data_names[i] << "__";
            out << data_names[std_cxx11::get<1>(vector_data_ranges[n_th_vector])];
          }

        out << "\" NumberOfComponents=\"3\" format=\""
            << ascii_or_binary << "\">\n";

        // now write data. pad all
        // vectors to have three
        // components
        std::vector<double> data;
        data.reserve (n_nodes*dim);

        for (unsigned int n=0; n<n_nodes; ++n)
          {
            switch (std_cxx11::get<1>(vector_data_ranges[n_th_vector]) -
                    std_cxx11::get<0>(vector_data_ranges[n_th_vector]))
              {
              case 0:
                data.push_back (data_vectors(std_cxx11::get<0>(vector_data_ranges[n_th_vector]), n));
                data.push_back (0);
                data.push_back (0);
                break;

              case 1:
                data.push_back (data_vectors(std_cxx11::get<0>(vector_data_ranges[n_th_vector]),   n));
                data.push_back (data_vectors(std_cxx11::get<0>(vector_data_ranges[n_th_vector])+1, n));
                data.push_back (0);
                break;
              case 2:
                data.push_back (data_vectors(std_cxx11::get<0>(vector_data_ranges[n_th_vector]),   n));
                data.push_back (data_vectors(std_cxx11::get<0>(vector_data_ranges[n_th_vector])+1, n));
                data.push_back (data_vectors(std_cxx11::get<0>(vector_data_ranges[n_th_vector])+2, n));
                break;

              default:
                // VTK doesn't
                // support
                // anything else
                // than vectors
                // with 1, 2, or
                // 3 components
                Assert (false, ExcInternalError());
              }
          }
        vtu_out << data;
        out << "    </DataArray>\n";
      }

    // now do the left over scalar data sets
    for (unsigned int data_set=0; data_set<n_data_sets; ++data_set)
      if (data_set_written[data_set] == false)
        {
          out << "    <DataArray type=\"Float64\" Name=\""
              << data_names[data_set]
              << "\" format=\""
              << ascii_or_binary << "\">\n";

          std::vector<double> data (data_vectors[data_set].begin(),
                                    data_vectors[data_set].end());
          vtu_out << data;
          out << "    </DataArray>\n";
        }

    out << "  </PointData>\n";

    // Finish up writing a valid XML file
    out << " </Piece>\n";

    // make sure everything now gets to
    // disk
    out.flush ();

    // assert the stream is still ok
    AssertThrow (out, ExcIO());
  }


  template <int dim, int spacedim>
  void write_svg (const std::vector<Patch<dim,spacedim> > &,
                  const std::vector<std::string> &,
                  const std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> > &,
                  const SvgFlags &,
                  std::ostream &)
  {
    Assert (false, ExcNotImplemented());
  }

  template <int spacedim>
  void write_svg (const std::vector<Patch<2,spacedim> > &patches,
                  const std::vector<std::string> &data_names,
                  const std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> > &vector_data_ranges,
                  const SvgFlags &flags,
                  std::ostream &out)
  {
    const int dim = 2;
    const unsigned int height = flags.height;
    unsigned int width = flags.width;

    // margin around the plotted area
    unsigned int margin_in_percent = 0;
    if (flags.margin) margin_in_percent = 5;


    // determine the bounding box in the model space
    double x_dimension, y_dimension, z_dimension;

    typename std::vector<Patch<dim,spacedim> >::const_iterator patch = patches.begin();

    unsigned int n_subdivisions = patch->n_subdivisions;
    unsigned int n = n_subdivisions + 1;
    const unsigned int d1 = 1;
    const unsigned int d2 = n;

    Point<spacedim> projected_point;
    Point<spacedim> projected_points[4];

    Point<2> projection_decomposition;
    Point<2> projection_decompositions[4];

    compute_node(projected_point, &*patch, 0, 0, 0, n_subdivisions);

    Assert ((flags.height_vector < patch->data.n_rows()) ||
            patch->data.n_rows() == 0,
            ExcIndexRange (flags.height_vector, 0, patch->data.n_rows()));

    double x_min = projected_point[0];
    double x_max = x_min;
    double y_min = projected_point[1];
    double y_max = y_min;
    double z_min = patch->data.n_rows() != 0 ? patch->data(flags.height_vector,0) : 0;
    double z_max = z_min;

    // iterate over the patches
    for (; patch != patches.end(); ++patch)
      {
        n_subdivisions = patch->n_subdivisions;
        n = n_subdivisions + 1;

        for (unsigned int i2 = 0; i2 < n_subdivisions; ++i2)
          {
            for (unsigned int i1 = 0; i1 < n_subdivisions; ++i1)
              {
                compute_node(projected_points[0], &*patch, i1, i2, 0, n_subdivisions);
                compute_node(projected_points[1], &*patch, i1+1, i2, 0, n_subdivisions);
                compute_node(projected_points[2], &*patch, i1, i2+1, 0, n_subdivisions);
                compute_node(projected_points[3], &*patch, i1+1, i2+1, 0, n_subdivisions);

                x_min = std::min(x_min, (double)projected_points[0][0]);
                x_min = std::min(x_min, (double)projected_points[1][0]);
                x_min = std::min(x_min, (double)projected_points[2][0]);
                x_min = std::min(x_min, (double)projected_points[3][0]);

                x_max = std::max(x_max, (double)projected_points[0][0]);
                x_max = std::max(x_max, (double)projected_points[1][0]);
                x_max = std::max(x_max, (double)projected_points[2][0]);
                x_max = std::max(x_max, (double)projected_points[3][0]);

                y_min = std::min(y_min, (double)projected_points[0][1]);
                y_min = std::min(y_min, (double)projected_points[1][1]);
                y_min = std::min(y_min, (double)projected_points[2][1]);
                y_min = std::min(y_min, (double)projected_points[3][1]);

                y_max = std::max(y_max, (double)projected_points[0][1]);
                y_max = std::max(y_max, (double)projected_points[1][1]);
                y_max = std::max(y_max, (double)projected_points[2][1]);
                y_max = std::max(y_max, (double)projected_points[3][1]);

                Assert ((flags.height_vector < patch->data.n_rows()) ||
                        patch->data.n_rows() == 0,
                        ExcIndexRange (flags.height_vector, 0, patch->data.n_rows()));

                z_min = std::min(z_min, (double)patch->data(flags.height_vector, i1*d1 + i2*d2));
                z_min = std::min(z_min, (double)patch->data(flags.height_vector, (i1+1)*d1 + i2*d2));
                z_min = std::min(z_min, (double)patch->data(flags.height_vector, i1*d1 + (i2+1)*d2));
                z_min = std::min(z_min, (double)patch->data(flags.height_vector, (i1+1)*d1 + (i2+1)*d2));

                z_max = std::max(z_max, (double)patch->data(flags.height_vector, i1*d1 + i2*d2));
                z_max = std::max(z_max, (double)patch->data(flags.height_vector, (i1+1)*d1 + i2*d2));
                z_max = std::max(z_max, (double)patch->data(flags.height_vector, i1*d1 + (i2+1)*d2));
                z_max = std::max(z_max, (double)patch->data(flags.height_vector, (i1+1)*d1 + (i2+1)*d2));
              }
          }
      }

    x_dimension = x_max - x_min;
    y_dimension = y_max - y_min;
    z_dimension = z_max - z_min;


// set initial camera position
    Point<3> camera_position(true);
    Point<3> camera_direction(true);
    Point<3> camera_horizontal(true);
    float camera_focus = 0;

    // translate camera from the origin to the initial position
    camera_position[0] = 0.;
    camera_position[1] = 0.;
    camera_position[2] = z_min + 2. * z_dimension;

    camera_direction[0] = 0.;
    camera_direction[1] = 0.;
    camera_direction[2] = - 1.;

    camera_horizontal[0] = 1.;
    camera_horizontal[1] = 0.;
    camera_horizontal[2] = 0.;

    camera_focus = .5 * z_dimension;

    Point<3> camera_position_temp;
    Point<3> camera_direction_temp;
    Point<3> camera_horizontal_temp;

    const float angle_factor = 3.14159265 / 180.;

    // (I) rotate the camera to the chosen polar angle
    camera_position_temp[1] = cos(angle_factor * flags.polar_angle) * camera_position[1] - sin(angle_factor * flags.polar_angle) * camera_position[2];
    camera_position_temp[2] = sin(angle_factor * flags.polar_angle) * camera_position[1] + cos(angle_factor * flags.polar_angle) * camera_position[2];

    camera_direction_temp[1] = cos(angle_factor * flags.polar_angle) * camera_direction[1] - sin(angle_factor * flags.polar_angle) * camera_direction[2];
    camera_direction_temp[2] = sin(angle_factor * flags.polar_angle) * camera_direction[1] + cos(angle_factor * flags.polar_angle) * camera_direction[2];

    camera_horizontal_temp[1] = cos(angle_factor * flags.polar_angle) * camera_horizontal[1] - sin(angle_factor * flags.polar_angle) * camera_horizontal[2];
    camera_horizontal_temp[2] = sin(angle_factor * flags.polar_angle) * camera_horizontal[1] + cos(angle_factor * flags.polar_angle) * camera_horizontal[2];

    camera_position[1] = camera_position_temp[1];
    camera_position[2] = camera_position_temp[2];

    camera_direction[1] = camera_direction_temp[1];
    camera_direction[2] = camera_direction_temp[2];

    camera_horizontal[1] = camera_horizontal_temp[1];
    camera_horizontal[2] = camera_horizontal_temp[2];

    // (II) rotate the camera to the chosen azimuth angle
    camera_position_temp[0] = cos(angle_factor * flags.azimuth_angle) * camera_position[0] - sin(angle_factor * flags.azimuth_angle) * camera_position[1];
    camera_position_temp[1] = sin(angle_factor * flags.azimuth_angle) * camera_position[0] + cos(angle_factor * flags.azimuth_angle) * camera_position[1];

    camera_direction_temp[0] = cos(angle_factor * flags.azimuth_angle) * camera_direction[0] - sin(angle_factor * flags.azimuth_angle) * camera_direction[1];
    camera_direction_temp[1] = sin(angle_factor * flags.azimuth_angle) * camera_direction[0] + cos(angle_factor * flags.azimuth_angle) * camera_direction[1];

    camera_horizontal_temp[0] = cos(angle_factor * flags.azimuth_angle) * camera_horizontal[0] - sin(angle_factor * flags.azimuth_angle) * camera_horizontal[1];
    camera_horizontal_temp[1] = sin(angle_factor * flags.azimuth_angle) * camera_horizontal[0] + cos(angle_factor * flags.azimuth_angle) * camera_horizontal[1];

    camera_position[0] = camera_position_temp[0];
    camera_position[1] = camera_position_temp[1];

    camera_direction[0] = camera_direction_temp[0];
    camera_direction[1] = camera_direction_temp[1];

    camera_horizontal[0] = camera_horizontal_temp[0];
    camera_horizontal[1] = camera_horizontal_temp[1];

    // (III) translate the camera
    camera_position[0] = x_min + .5 * x_dimension;
    camera_position[1] = y_min + .5 * y_dimension;

    camera_position[0] += (z_min + 2. * z_dimension) * sin(angle_factor * flags.polar_angle) * sin(angle_factor * flags.azimuth_angle);
    camera_position[1] -= (z_min + 2. * z_dimension) * sin(angle_factor * flags.polar_angle) * cos(angle_factor * flags.azimuth_angle);


// determine the bounding box on the projection plane
    double x_min_perspective, y_min_perspective;
    double x_max_perspective, y_max_perspective;
    double x_dimension_perspective, y_dimension_perspective;

    patch = patches.begin();

    n_subdivisions = patch->n_subdivisions;
    n = n_subdivisions + 1;

    Point<3> point(true);

    compute_node(projected_point, &*patch, 0, 0, 0, n_subdivisions);

    Assert ((flags.height_vector < patch->data.n_rows()) ||
            patch->data.n_rows() == 0,
            ExcIndexRange (flags.height_vector, 0, patch->data.n_rows()));

    point[0] = projected_point[0];
    point[1] = projected_point[1];
    point[2] = patch->data.n_rows() != 0 ? patch->data(flags.height_vector, 0) : 0;

    projection_decomposition = svg_project_point(point, camera_position, camera_direction, camera_horizontal, camera_focus);

    x_min_perspective = projection_decomposition[0];
    x_max_perspective = projection_decomposition[0];
    y_min_perspective = projection_decomposition[1];
    y_max_perspective = projection_decomposition[1];

    // iterate over the patches
    for (; patch != patches.end(); ++patch)
      {
        n_subdivisions = patch->n_subdivisions;
        n = n_subdivisions + 1;

        for (unsigned int i2 = 0; i2 < n_subdivisions; ++i2)
          {
            for (unsigned int i1 = 0; i1 < n_subdivisions; ++i1)
              {
                Point<spacedim> projected_vertices[4];
                Point<3> vertices[4];

                compute_node(projected_vertices[0], &*patch, i1, i2, 0, n_subdivisions);
                compute_node(projected_vertices[1], &*patch, i1+1, i2, 0, n_subdivisions);
                compute_node(projected_vertices[2], &*patch, i1, i2+1, 0, n_subdivisions);
                compute_node(projected_vertices[3], &*patch, i1+1, i2+1, 0, n_subdivisions);

                Assert ((flags.height_vector < patch->data.n_rows()) ||
                        patch->data.n_rows() == 0,
                        ExcIndexRange (flags.height_vector, 0, patch->data.n_rows()));

                vertices[0][0] = projected_vertices[0][0];
                vertices[0][1] = projected_vertices[0][1];
                vertices[0][2] = patch->data.n_rows() != 0 ? patch->data(0,i1*d1 + i2*d2) : 0;

                vertices[1][0] = projected_vertices[1][0];
                vertices[1][1] = projected_vertices[1][1];
                vertices[1][2] = patch->data.n_rows() != 0 ? patch->data(0,(i1+1)*d1 + i2*d2) : 0;

                vertices[2][0] = projected_vertices[2][0];
                vertices[2][1] = projected_vertices[2][1];
                vertices[2][2] = patch->data.n_rows() != 0 ? patch->data(0,i1*d1 + (i2+1)*d2) : 0;

                vertices[3][0] = projected_vertices[3][0];
                vertices[3][1] = projected_vertices[3][1];
                vertices[3][2] = patch->data.n_rows() != 0 ? patch->data(0,(i1+1)*d1 + (i2+1)*d2) : 0;

                projection_decompositions[0] = svg_project_point(vertices[0], camera_position, camera_direction, camera_horizontal, camera_focus);
                projection_decompositions[1] = svg_project_point(vertices[1], camera_position, camera_direction, camera_horizontal, camera_focus);
                projection_decompositions[2] = svg_project_point(vertices[2], camera_position, camera_direction, camera_horizontal, camera_focus);
                projection_decompositions[3] = svg_project_point(vertices[3], camera_position, camera_direction, camera_horizontal, camera_focus);

                x_min_perspective = std::min(x_min_perspective, (double)projection_decompositions[0][0]);
                x_min_perspective = std::min(x_min_perspective, (double)projection_decompositions[1][0]);
                x_min_perspective = std::min(x_min_perspective, (double)projection_decompositions[2][0]);
                x_min_perspective = std::min(x_min_perspective, (double)projection_decompositions[3][0]);

                x_max_perspective = std::max(x_max_perspective, (double)projection_decompositions[0][0]);
                x_max_perspective = std::max(x_max_perspective, (double)projection_decompositions[1][0]);
                x_max_perspective = std::max(x_max_perspective, (double)projection_decompositions[2][0]);
                x_max_perspective = std::max(x_max_perspective, (double)projection_decompositions[3][0]);

                y_min_perspective = std::min(y_min_perspective, (double)projection_decompositions[0][1]);
                y_min_perspective = std::min(y_min_perspective, (double)projection_decompositions[1][1]);
                y_min_perspective = std::min(y_min_perspective, (double)projection_decompositions[2][1]);
                y_min_perspective = std::min(y_min_perspective, (double)projection_decompositions[3][1]);

                y_max_perspective = std::max(y_max_perspective, (double)projection_decompositions[0][1]);
                y_max_perspective = std::max(y_max_perspective, (double)projection_decompositions[1][1]);
                y_max_perspective = std::max(y_max_perspective, (double)projection_decompositions[2][1]);
                y_max_perspective = std::max(y_max_perspective, (double)projection_decompositions[3][1]);
              }
          }
      }

    x_dimension_perspective = x_max_perspective - x_min_perspective;
    y_dimension_perspective = y_max_perspective - y_min_perspective;

    std::multiset<SvgCell> cells;

    // iterate over the patches
    for (patch = patches.begin(); patch != patches.end(); ++patch)
      {
        n_subdivisions = patch->n_subdivisions;
        n = n_subdivisions + 1;

        for (unsigned int i2 = 0; i2 < n_subdivisions; ++i2)
          {
            for (unsigned int i1 = 0; i1 < n_subdivisions; ++i1)
              {
                Point<spacedim> projected_vertices[4];
                SvgCell cell;

                compute_node(projected_vertices[0], &*patch, i1, i2, 0, n_subdivisions);
                compute_node(projected_vertices[1], &*patch, i1+1, i2, 0, n_subdivisions);
                compute_node(projected_vertices[2], &*patch, i1, i2+1, 0, n_subdivisions);
                compute_node(projected_vertices[3], &*patch, i1+1, i2+1, 0, n_subdivisions);

                Assert ((flags.height_vector < patch->data.n_rows()) ||
                        patch->data.n_rows() == 0,
                        ExcIndexRange (flags.height_vector, 0, patch->data.n_rows()));

                cell.vertices[0][0] = projected_vertices[0][0];
                cell.vertices[0][1] = projected_vertices[0][1];
                cell.vertices[0][2] = patch->data.n_rows() != 0 ? patch->data(0,i1*d1 + i2*d2) : 0;

                cell.vertices[1][0] = projected_vertices[1][0];
                cell.vertices[1][1] = projected_vertices[1][1];
                cell.vertices[1][2] = patch->data.n_rows() != 0 ? patch->data(0,(i1+1)*d1 + i2*d2) : 0;

                cell.vertices[2][0] = projected_vertices[2][0];
                cell.vertices[2][1] = projected_vertices[2][1];
                cell.vertices[2][2] = patch->data.n_rows() != 0 ? patch->data(0,i1*d1 + (i2+1)*d2) : 0;

                cell.vertices[3][0] = projected_vertices[3][0];
                cell.vertices[3][1] = projected_vertices[3][1];
                cell.vertices[3][2] = patch->data.n_rows() != 0 ? patch->data(0,(i1+1)*d1 + (i2+1)*d2) : 0;

                cell.projected_vertices[0] = svg_project_point(cell.vertices[0], camera_position, camera_direction, camera_horizontal, camera_focus);
                cell.projected_vertices[1] = svg_project_point(cell.vertices[1], camera_position, camera_direction, camera_horizontal, camera_focus);
                cell.projected_vertices[2] = svg_project_point(cell.vertices[2], camera_position, camera_direction, camera_horizontal, camera_focus);
                cell.projected_vertices[3] = svg_project_point(cell.vertices[3], camera_position, camera_direction, camera_horizontal, camera_focus);

                cell.center = .25 * (cell.vertices[0] + cell.vertices[1] + cell.vertices[2] + cell.vertices[3]);
                cell.projected_center = svg_project_point(cell.center, camera_position, camera_direction, camera_horizontal, camera_focus);

                cell.depth = cell.center.distance(camera_position);

                cells.insert(cell);
              }
          }
      }


    // write the svg file
    if (width==0)
      width = static_cast<unsigned int>(.5 + height * (x_dimension_perspective / y_dimension_perspective));
    unsigned int additional_width = 0;

    if (flags.draw_colorbar) additional_width = static_cast<unsigned int>(.5 + height * .3); // additional width for colorbar

    // basic svg header and background rectangle
    out << "<svg width=\"" << width + additional_width << "\" height=\"" << height << "\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">" << '\n'
        << " <rect width=\"" << width + additional_width << "\" height=\"" << height << "\" style=\"fill:white\"/>" << '\n' << '\n';

    unsigned int triangle_counter = 0;

    // write the cells in the correct order
    for (typename std::multiset<SvgCell>::const_iterator cell = cells.begin(); cell != cells.end(); ++cell)
      {
        Point<3> points3d_triangle[3];

        for (unsigned int triangle_index = 0; triangle_index < 4; triangle_index++)
          {
            switch (triangle_index)
              {
              case 0:
                points3d_triangle[0] = cell->vertices[0], points3d_triangle[1] = cell->vertices[1], points3d_triangle[2] = cell->center;
                break;
              case 1:
                points3d_triangle[0] = cell->vertices[1], points3d_triangle[1] = cell->vertices[3], points3d_triangle[2] = cell->center;
                break;
              case 2:
                points3d_triangle[0] = cell->vertices[3], points3d_triangle[1] = cell->vertices[2], points3d_triangle[2] = cell->center;
                break;
              case 3:
                points3d_triangle[0] = cell->vertices[2], points3d_triangle[1] = cell->vertices[0], points3d_triangle[2] = cell->center;
                break;
              default:
                break;
              }

            Point<6> gradient_param = svg_get_gradient_parameters(points3d_triangle);

            double start_h = .667 - ((gradient_param[4] - z_min) / z_dimension) * .667;
            double stop_h = .667 - ((gradient_param[5] - z_min) / z_dimension) * .667;

            unsigned int start_r = 0;
            unsigned int start_g = 0;
            unsigned int start_b = 0;

            unsigned int stop_r = 0;
            unsigned int stop_g = 0;
            unsigned int stop_b = 0;

            unsigned int start_i = static_cast<unsigned int>(start_h * 6.);
            unsigned int stop_i = static_cast<unsigned int>(stop_h * 6.);

            double start_f = start_h * 6. - start_i;
            double start_q = 1. - start_f;

            double stop_f = stop_h * 6. - stop_i;
            double stop_q = 1. - stop_f;

            switch (start_i % 6)
              {
              case 0:
                start_r = 255, start_g = static_cast<unsigned int>(.5 + 255. * start_f);
                break;
              case 1:
                start_r = static_cast<unsigned int>(.5 + 255. * start_q), start_g = 255;
                break;
              case 2:
                start_g = 255, start_b = static_cast<unsigned int>(.5 + 255. * start_f);
                break;
              case 3:
                start_g = static_cast<unsigned int>(.5 + 255. * start_q), start_b = 255;
                break;
              case 4:
                start_r = static_cast<unsigned int>(.5 + 255. * start_f), start_b = 255;
                break;
              case 5:
                start_r = 255, start_b = static_cast<unsigned int>(.5 + 255. * start_q);
                break;
              default:
                break;
              }

            switch (stop_i % 6)
              {
              case 0:
                stop_r = 255, stop_g = static_cast<unsigned int>(.5 + 255. * stop_f);
                break;
              case 1:
                stop_r = static_cast<unsigned int>(.5 + 255. * stop_q), stop_g = 255;
                break;
              case 2:
                stop_g = 255, stop_b = static_cast<unsigned int>(.5 + 255. * stop_f);
                break;
              case 3:
                stop_g = static_cast<unsigned int>(.5 + 255. * stop_q), stop_b = 255;
                break;
              case 4:
                stop_r = static_cast<unsigned int>(.5 + 255. * stop_f), stop_b = 255;
                break;
              case 5:
                stop_r = 255, stop_b = static_cast<unsigned int>(.5 + 255. * stop_q);
                break;
              default:
                break;
              }

            Point<3> gradient_start_point_3d, gradient_stop_point_3d;

            gradient_start_point_3d[0] = gradient_param[0];
            gradient_start_point_3d[1] = gradient_param[1];
            gradient_start_point_3d[2] = gradient_param[4];

            gradient_stop_point_3d[0] = gradient_param[2];
            gradient_stop_point_3d[1] = gradient_param[3];
            gradient_stop_point_3d[2] = gradient_param[5];

            Point<2> gradient_start_point = svg_project_point(gradient_start_point_3d, camera_position, camera_direction, camera_horizontal, camera_focus);
            Point<2> gradient_stop_point = svg_project_point(gradient_stop_point_3d, camera_position, camera_direction, camera_horizontal, camera_focus);

            // define linear gradient
            out << "  <linearGradient id=\"" << triangle_counter << "\" gradientUnits=\"userSpaceOnUse\" "
                << "x1=\""
                << static_cast<unsigned int>(.5 + ((gradient_start_point[0] - x_min_perspective) / x_dimension_perspective) * (width - (width/100.) * 2. * margin_in_percent) + ((width/100.) * margin_in_percent))
                << "\" "
                << "y1=\""
                << static_cast<unsigned int>(.5 + height - (height/100.) * margin_in_percent - ((gradient_start_point[1] - y_min_perspective) / y_dimension_perspective) * (height - (height/100.) * 2. * margin_in_percent))
                << "\" "
                << "x2=\""
                << static_cast<unsigned int>(.5 + ((gradient_stop_point[0] - x_min_perspective) / x_dimension_perspective) * (width - (width/100.) * 2. * margin_in_percent) + ((width/100.) * margin_in_percent))
                << "\" "
                << "y2=\""
                << static_cast<unsigned int>(.5 + height - (height/100.) * margin_in_percent - ((gradient_stop_point[1] - y_min_perspective) / y_dimension_perspective) * (height - (height/100.) * 2. * margin_in_percent))
                << "\""
                << ">" << '\n'
                << "   <stop offset=\"0\" style=\"stop-color:rgb(" << start_r << "," << start_g << "," << start_b << ")\"/>" << '\n'
                << "   <stop offset=\"1\" style=\"stop-color:rgb(" << stop_r << "," << stop_g << "," << stop_b << ")\"/>" << '\n'
                << "  </linearGradient>" << '\n';

            // draw current triangle
            double x1 = 0, y1 = 0, x2 = 0, y2 = 0;
            double x3 = cell->projected_center[0];
            double y3 = cell->projected_center[1];

            switch (triangle_index)
              {
              case 0:
                x1 = cell->projected_vertices[0][0], y1 = cell->projected_vertices[0][1], x2 = cell->projected_vertices[1][0], y2 = cell->projected_vertices[1][1];
                break;
              case 1:
                x1 = cell->projected_vertices[1][0], y1 = cell->projected_vertices[1][1], x2 = cell->projected_vertices[3][0], y2 = cell->projected_vertices[3][1];
                break;
              case 2:
                x1 = cell->projected_vertices[3][0], y1 = cell->projected_vertices[3][1], x2 = cell->projected_vertices[2][0], y2 = cell->projected_vertices[2][1];
                break;
              case 3:
                x1 = cell->projected_vertices[2][0], y1 = cell->projected_vertices[2][1], x2 = cell->projected_vertices[0][0], y2 = cell->projected_vertices[0][1];
                break;
              default:
                break;
              }

            out << "  <path d=\"M "
                << static_cast<unsigned int>(.5 + ((x1 - x_min_perspective) / x_dimension_perspective) * (width - (width/100.) * 2. * margin_in_percent) + ((width/100.) * margin_in_percent))
                << ' '
                << static_cast<unsigned int>(.5 + height - (height/100.) * margin_in_percent - ((y1 - y_min_perspective) / y_dimension_perspective) * (height - (height/100.) * 2. * margin_in_percent))
                << " L "
                << static_cast<unsigned int>(.5 + ((x2 - x_min_perspective) / x_dimension_perspective) * (width - (width/100.) * 2. * margin_in_percent) + ((width/100.) * margin_in_percent))
                << ' '
                << static_cast<unsigned int>(.5 + height - (height/100.) * margin_in_percent - ((y2 - y_min_perspective) / y_dimension_perspective) * (height - (height/100.) * 2. * margin_in_percent))
                << " L "
                << static_cast<unsigned int>(.5 + ((x3 - x_min_perspective) / x_dimension_perspective) * (width - (width/100.) * 2. * margin_in_percent) + ((width/100.) * margin_in_percent))
                << ' '
                << static_cast<unsigned int>(.5 + height - (height/100.) * margin_in_percent - ((y3 - y_min_perspective) / y_dimension_perspective) * (height - (height/100.) * 2. * margin_in_percent))
                << " L "
                << static_cast<unsigned int>(.5 + ((x1 - x_min_perspective) / x_dimension_perspective) * (width - (width/100.) * 2. * margin_in_percent) + ((width/100.) * margin_in_percent))
                << ' '
                << static_cast<unsigned int>(.5 + height - (height/100.) * margin_in_percent - ((y1 - y_min_perspective) / y_dimension_perspective) * (height - (height/100.) * 2. * margin_in_percent))
                << "\" style=\"stroke:black; fill:url(#" << triangle_counter << "); stroke-width:" << flags.line_thickness << "\"/>" << '\n';

            triangle_counter++;
          }
      }


// draw the colorbar
    if (flags.draw_colorbar)
      {
        out << '\n' << " <!-- colorbar -->" << '\n';

        unsigned int element_height = static_cast<unsigned int>(((height/100.) * (71. - 2.*margin_in_percent)) / 4);
        unsigned int element_width = static_cast<unsigned int>(.5 + (height/100.) * 2.5);

        additional_width = 0;
        if (!flags.margin) additional_width = static_cast<unsigned int>(.5 + (height/100.) * 2.5);

        for (unsigned int index = 0; index < 4; index++)
          {
            double start_h = .667 - ((index+1) / 4.) * .667;
            double stop_h = .667 - (index / 4.) * .667;

            unsigned int start_r = 0;
            unsigned int start_g = 0;
            unsigned int start_b = 0;

            unsigned int stop_r = 0;
            unsigned int stop_g = 0;
            unsigned int stop_b = 0;

            unsigned int start_i = static_cast<unsigned int>(start_h * 6.);
            unsigned int stop_i = static_cast<unsigned int>(stop_h * 6.);

            double start_f = start_h * 6. - start_i;
            double start_q = 1. - start_f;

            double stop_f = stop_h * 6. - stop_i;
            double stop_q = 1. - stop_f;

            switch (start_i % 6)
              {
              case 0:
                start_r = 255, start_g = static_cast<unsigned int>(.5 + 255. * start_f);
                break;
              case 1:
                start_r = static_cast<unsigned int>(.5 + 255. * start_q), start_g = 255;
                break;
              case 2:
                start_g = 255, start_b = static_cast<unsigned int>(.5 + 255. * start_f);
                break;
              case 3:
                start_g = static_cast<unsigned int>(.5 + 255. * start_q), start_b = 255;
                break;
              case 4:
                start_r = static_cast<unsigned int>(.5 + 255. * start_f), start_b = 255;
                break;
              case 5:
                start_r = 255, start_b = static_cast<unsigned int>(.5 + 255. * start_q);
                break;
              default:
                break;
              }

            switch (stop_i % 6)
              {
              case 0:
                stop_r = 255, stop_g = static_cast<unsigned int>(.5 + 255. * stop_f);
                break;
              case 1:
                stop_r = static_cast<unsigned int>(.5 + 255. * stop_q), stop_g = 255;
                break;
              case 2:
                stop_g = 255, stop_b = static_cast<unsigned int>(.5 + 255. * stop_f);
                break;
              case 3:
                stop_g = static_cast<unsigned int>(.5 + 255. * stop_q), stop_b = 255;
                break;
              case 4:
                stop_r = static_cast<unsigned int>(.5 + 255. * stop_f), stop_b = 255;
                break;
              case 5:
                stop_r = 255, stop_b = static_cast<unsigned int>(.5 + 255. * stop_q);
                break;
              default:
                break;
              }

            // define gradient
            out << "  <linearGradient id=\"colorbar_" << index << "\" gradientUnits=\"userSpaceOnUse\" "
                << "x1=\"" << width + additional_width << "\" "
                << "y1=\"" << static_cast<unsigned int>(.5 + (height/100.) * (margin_in_percent + 29)) + (3-index) * element_height << "\" "
                << "x2=\"" << width + additional_width << "\" "
                << "y2=\"" << static_cast<unsigned int>(.5 + (height/100.) * (margin_in_percent + 29)) + (4-index) * element_height << "\""
                << ">" << '\n'
                << "   <stop offset=\"0\" style=\"stop-color:rgb(" << start_r << "," << start_g << "," << start_b << ")\"/>" << '\n'
                << "   <stop offset=\"1\" style=\"stop-color:rgb(" << stop_r << "," << stop_g << "," << stop_b << ")\"/>" << '\n'
                << "  </linearGradient>" << '\n';

            // draw box corresponding to the gradient above
            out << "  <rect"
                << " x=\"" << width + additional_width
                << "\" y=\"" << static_cast<unsigned int>(.5 + (height/100.) * (margin_in_percent + 29)) + (3-index) * element_height
                << "\" width=\"" << element_width
                << "\" height=\"" << element_height
                << "\" style=\"stroke:black; stroke-width:2; fill:url(#colorbar_" << index << ")\"/>" << '\n';
          }

        for (unsigned int index = 0; index < 5; index++)
          {
            out << "  <text x=\"" << width + additional_width + static_cast<unsigned int>(1.5 * element_width)
                << "\" y=\"" << static_cast<unsigned int>(.5 + (height/100.) * (margin_in_percent + 29) + (4.-index) * element_height + 30.) << "\""
                << " style=\"text-anchor:start; font-size:80; font-family:Helvetica";

            if (index == 0 || index == 4) out << "; font-weight:bold";

            out << "\">" << (float)(((int)((z_min + index * (z_dimension / 4.))*10000))/10000.);

            if (index == 4) out << " max";
            if (index == 0) out << " min";

            out << "</text>" << '\n';
          }
      }

    // finalize the svg file
    out << '\n' << "</svg>";
    out.flush();

  }



  template <int dim, int spacedim>
  void
  write_deal_II_intermediate (const std::vector<Patch<dim,spacedim> > &patches,
                              const std::vector<std::string>          &data_names,
                              const std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> > &vector_data_ranges,
                              const Deal_II_IntermediateFlags         &/*flags*/,
                              std::ostream                            &out)
  {
    AssertThrow (out, ExcIO());

    // first write tokens indicating the
    // template parameters. we need this in
    // here because we may want to read in data
    // again even if we don't know in advance
    // the template parameters, see step-19
    out << dim << ' ' << spacedim << '\n';

    // then write a header
    out << "[deal.II intermediate format graphics data]" << '\n'
        << "[written by " << DEAL_II_PACKAGE_NAME << " " << DEAL_II_PACKAGE_VERSION << "]" << '\n'
        << "[Version: " << Deal_II_IntermediateFlags::format_version << "]" << '\n';

    out << data_names.size() << '\n';
    for (unsigned int i=0; i<data_names.size(); ++i)
      out << data_names[i] << '\n';

    out << patches.size() << '\n';
    for (unsigned int i=0; i<patches.size(); ++i)
      out << patches[i] << '\n';

    out << vector_data_ranges.size() << '\n';
    for (unsigned int i=0; i<vector_data_ranges.size(); ++i)
      out << std_cxx11::get<0>(vector_data_ranges[i]) << ' '
          << std_cxx11::get<1>(vector_data_ranges[i]) << '\n'
          << std_cxx11::get<2>(vector_data_ranges[i]) << '\n';

    out << '\n';
    // make sure everything now gets to
    // disk
    out.flush ();
  }



  std::pair<unsigned int, unsigned int>
  determine_intermediate_format_dimensions (std::istream &input)
  {
    Assert (input, ExcIO());

    unsigned int dim, spacedim;
    input >> dim >> spacedim;

    return std::make_pair (dim, spacedim);
  }
} // namespace DataOutBase




/* --------------------------- class DataOutInterface ---------------------- */


template <int dim, int spacedim>
DataOutInterface<dim,spacedim>::DataOutInterface ()
  : default_subdivisions(1)
{}


template <int dim, int spacedim>
DataOutInterface<dim,spacedim>::~DataOutInterface ()
{}




template <int dim, int spacedim>
void DataOutInterface<dim,spacedim>::write_dx (std::ostream &out) const
{
  DataOutBase::write_dx (get_patches(), get_dataset_names(),
                         get_vector_data_ranges(),
                         dx_flags, out);
}



template <int dim, int spacedim>
void DataOutInterface<dim,spacedim>::write_ucd (std::ostream &out) const
{
  DataOutBase::write_ucd (get_patches(), get_dataset_names(),
                          get_vector_data_ranges(),
                          ucd_flags, out);
}



template <int dim, int spacedim>
void DataOutInterface<dim,spacedim>::write_gnuplot (std::ostream &out) const
{
  DataOutBase::write_gnuplot (get_patches(), get_dataset_names(),
                              get_vector_data_ranges(),
                              gnuplot_flags, out);
}



template <int dim, int spacedim>
void DataOutInterface<dim,spacedim>::write_povray (std::ostream &out) const
{
  DataOutBase::write_povray (get_patches(), get_dataset_names(),
                             get_vector_data_ranges(),
                             povray_flags, out);
}



template <int dim, int spacedim>
void DataOutInterface<dim,spacedim>::write_eps (std::ostream &out) const
{
  DataOutBase::write_eps (get_patches(), get_dataset_names(),
                          get_vector_data_ranges(),
                          eps_flags, out);
}



template <int dim, int spacedim>
void DataOutInterface<dim,spacedim>::write_gmv (std::ostream &out) const
{
  DataOutBase::write_gmv (get_patches(), get_dataset_names(),
                          get_vector_data_ranges(),
                          gmv_flags, out);
}



template <int dim, int spacedim>
void DataOutInterface<dim,spacedim>::write_tecplot (std::ostream &out) const
{
  DataOutBase::write_tecplot (get_patches(), get_dataset_names(),
                              get_vector_data_ranges(),
                              tecplot_flags, out);
}



template <int dim, int spacedim>
void DataOutInterface<dim,spacedim>::write_tecplot_binary (std::ostream &out) const
{
  DataOutBase::write_tecplot_binary (get_patches(), get_dataset_names(),
                                     get_vector_data_ranges(),
                                     tecplot_flags, out);
}



template <int dim, int spacedim>
void DataOutInterface<dim,spacedim>::write_vtk (std::ostream &out) const
{
  DataOutBase::write_vtk (get_patches(), get_dataset_names(),
                          get_vector_data_ranges(),
                          vtk_flags, out);
}

template <int dim, int spacedim>
void DataOutInterface<dim,spacedim>::write_vtu (std::ostream &out) const
{
  DataOutBase::write_vtu (get_patches(), get_dataset_names(),
                          get_vector_data_ranges(),
                          vtk_flags, out);
}

template <int dim, int spacedim>
void DataOutInterface<dim,spacedim>::write_svg (std::ostream &out) const
{
  DataOutBase::write_svg (get_patches(), get_dataset_names(),
                          get_vector_data_ranges(),
                          svg_flags, out);
}

template <int dim, int spacedim>
void DataOutInterface<dim,spacedim>::write_vtu_in_parallel (const char *filename, MPI_Comm comm) const
{
#ifndef DEAL_II_WITH_MPI
  //without MPI fall back to the normal way to write a vtu file:
  (void)comm;

  std::ofstream f(filename);
  write_vtu (f);
#else

  int myrank, nproc, err;
  MPI_Comm_rank(comm, &myrank);
  MPI_Comm_size(comm, &nproc);

  MPI_Info info;
  MPI_Info_create(&info);
  MPI_File fh;
  err = MPI_File_open(comm, const_cast<char *>(filename),
                      MPI_MODE_CREATE | MPI_MODE_WRONLY, info, &fh);
  AssertThrow(err==0, ExcMessage("Unable to open file with MPI_File_open!"));


  MPI_File_set_size(fh, 0); // delete the file contents
  // this barrier is necessary, because otherwise others might already
  // write while one core is still setting the size to zero.
  MPI_Barrier(comm);
  MPI_Info_free(&info);

  unsigned int header_size;

  //write header
  if (myrank==0)
    {
      std::stringstream ss;
      DataOutBase::write_vtu_header(ss, vtk_flags);
      header_size = ss.str().size();
      MPI_File_write(fh, const_cast<char *>(ss.str().c_str()), header_size, MPI_CHAR, NULL);
    }

  MPI_Bcast(&header_size, 1, MPI_INT, 0, comm);

  MPI_File_seek_shared( fh, header_size, MPI_SEEK_SET );
  {
    std::stringstream ss;
    DataOutBase::write_vtu_main (get_patches(), get_dataset_names(),
                                 get_vector_data_ranges(),
                                 vtk_flags, ss);
    MPI_File_write_ordered(fh, const_cast<char *>(ss.str().c_str()), ss.str().size(), MPI_CHAR, NULL);
  }

  //write footer
  if (myrank==0)
    {
      std::stringstream ss;
      DataOutBase::write_vtu_footer(ss);
      unsigned int footer_size = ss.str().size();
      MPI_File_write_shared(fh, const_cast<char *>(ss.str().c_str()), footer_size, MPI_CHAR, NULL);
    }
  MPI_File_close( &fh );
#endif
}


template <int dim, int spacedim>
void
DataOutInterface<dim,spacedim>::
write_pvd_record (std::ostream &out,
                  const std::vector<std::pair<double,std::string> >  &times_and_names) const
{
  AssertThrow (out, ExcIO());

  out << "<?xml version=\"1.0\"?>\n";

  std::time_t  time1= std::time (0);
  std::tm     *time = std::localtime(&time1);
  out << "<!--\n";
  out << "#This file was generated by the deal.II library on "
      << time->tm_year+1900 << "/"
      << time->tm_mon+1 << "/"
      << time->tm_mday << " at "
      << time->tm_hour << ":"
      << std::setw(2) << time->tm_min << ":"
      << std::setw(2) << time->tm_sec
      << "\n-->\n";

  out << "<VTKFile type=\"Collection\" version=\"0.1\" ByteOrder=\"LittleEndian\">\n";
  out << "  <Collection>\n";

  for (unsigned int i=0; i<times_and_names.size(); ++i)
    out << "    <DataSet timestep=\"" << times_and_names[i].first
        << "\" group=\"\" part=\"0\" file=\"" << times_and_names[i].second
        << "\"/>\n";

  out << "  </Collection>\n";
  out << "</VTKFile>\n";

  out.flush();

  AssertThrow (out, ExcIO());
}


template <int dim, int spacedim>
void
DataOutInterface<dim,spacedim>::write_pvtu_record (std::ostream &out,
                                                   const std::vector<std::string> &piece_names) const
{
  AssertThrow (out, ExcIO());

  const std::vector<std::string> data_names = get_dataset_names();
  const std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> > vector_data_ranges
    = get_vector_data_ranges();

  const unsigned int n_data_sets = data_names.size();

  out << "<?xml version=\"1.0\"?>\n";

  std::time_t  time1= std::time (0);
  std::tm     *time = std::localtime(&time1);
  out << "<!--\n";
  out << "#This file was generated by the deal.II library on "
      << time->tm_year+1900 << "/"
      << time->tm_mon+1 << "/"
      << time->tm_mday << " at "
      << time->tm_hour << ":"
      << std::setw(2) << time->tm_min << ":"
      << std::setw(2) << time->tm_sec
      << "\n-->\n";

  out << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
  out << "  <PUnstructuredGrid GhostLevel=\"0\">\n";
  out << "    <PPointData Scalars=\"scalars\">\n";

  // We need to output in the same order as
  // the write_vtu function does:
  std::vector<bool> data_set_written (n_data_sets, false);
  for (unsigned int n_th_vector=0; n_th_vector<vector_data_ranges.size(); ++n_th_vector)
    {
      AssertThrow (std_cxx11::get<1>(vector_data_ranges[n_th_vector]) >=
                   std_cxx11::get<0>(vector_data_ranges[n_th_vector]),
                   ExcLowerRange (std_cxx11::get<1>(vector_data_ranges[n_th_vector]),
                                  std_cxx11::get<0>(vector_data_ranges[n_th_vector])));
      AssertThrow (std_cxx11::get<1>(vector_data_ranges[n_th_vector]) < n_data_sets,
                   ExcIndexRange (std_cxx11::get<1>(vector_data_ranges[n_th_vector]),
                                  0, n_data_sets));
      AssertThrow (std_cxx11::get<1>(vector_data_ranges[n_th_vector]) + 1
                   - std_cxx11::get<0>(vector_data_ranges[n_th_vector]) <= 3,
                   ExcMessage ("Can't declare a vector with more than 3 components "
                               "in VTK"));

      // mark these components as already
      // written:
      for (unsigned int i=std_cxx11::get<0>(vector_data_ranges[n_th_vector]);
           i<=std_cxx11::get<1>(vector_data_ranges[n_th_vector]);
           ++i)
        data_set_written[i] = true;

      // write the
      // header. concatenate all the
      // component names with double
      // underscores unless a vector
      // name has been specified
      out << "    <PDataArray type=\"Float64\" Name=\"";

      if (std_cxx11::get<2>(vector_data_ranges[n_th_vector]) != "")
        out << std_cxx11::get<2>(vector_data_ranges[n_th_vector]);
      else
        {
          for (unsigned int i=std_cxx11::get<0>(vector_data_ranges[n_th_vector]);
               i<std_cxx11::get<1>(vector_data_ranges[n_th_vector]);
               ++i)
            out << data_names[i] << "__";
          out << data_names[std_cxx11::get<1>(vector_data_ranges[n_th_vector])];
        }

      out << "\" NumberOfComponents=\"3\" format=\"ascii\"/>\n";
    }

  for (unsigned int data_set=0; data_set<n_data_sets; ++data_set)
    if (data_set_written[data_set] == false)
      {
        out << "    <PDataArray type=\"Float64\" Name=\""
            << data_names[data_set]
            << "\" format=\"ascii\"/>\n";
      }

  out << "    </PPointData>\n";

  out << "    <PPoints>\n";
  out << "      <PDataArray type=\"Float64\" NumberOfComponents=\"3\"/>\n";
  out << "    </PPoints>\n";

  for (unsigned int i=0; i<piece_names.size(); ++i)
    out << "    <Piece Source=\"" << piece_names[i] << "\"/>\n";

  out << "  </PUnstructuredGrid>\n";
  out << "</VTKFile>\n";

  out.flush();

  // assert the stream is still ok
  AssertThrow (out, ExcIO());
}



template <int dim, int spacedim>
void
DataOutInterface<dim,spacedim>::write_visit_record (std::ostream &out,
                                                    const std::vector<std::string> &piece_names) const
{
  out << "!NBLOCKS " << piece_names.size() << '\n';
  for (unsigned int i=0; i<piece_names.size(); ++i)
    out << piece_names[i] << '\n';

  out << std::flush;
}



template <int dim, int spacedim>
void
DataOutInterface<dim,spacedim>::write_visit_record (std::ostream &out,
                                                    const std::vector<std::vector<std::string> > &piece_names) const
{
  AssertThrow (out, ExcIO());

  if (piece_names.size() == 0)
    return;

  const double nblocks = piece_names[0].size();
  Assert(nblocks > 0, ExcMessage("piece_names should be a vector of nonempty vectors.") )

  out << "!NBLOCKS " << nblocks << '\n';
  for (std::vector<std::vector<std::string> >::const_iterator domain = piece_names.begin(); domain != piece_names.end(); ++domain)
    {
      Assert(domain->size() == nblocks, ExcMessage("piece_names should be a vector of equal sized vectors.") )
      for (std::vector<std::string>::const_iterator subdomain = domain->begin(); subdomain != domain->end(); ++subdomain)
        out << *subdomain << '\n';
    }

  out << std::flush;
}



template <int dim, int spacedim>
void DataOutInterface<dim,spacedim>::
write_deal_II_intermediate (std::ostream &out) const
{
  DataOutBase::write_deal_II_intermediate (get_patches(), get_dataset_names(),
                                           get_vector_data_ranges(),
                                           deal_II_intermediate_flags, out);
}


template <int dim, int spacedim>
XDMFEntry DataOutInterface<dim,spacedim>::
create_xdmf_entry (const std::string &h5_filename, const double cur_time, MPI_Comm comm) const
{
  DataOutBase::DataOutFilter  data_filter(DataOutBase::DataOutFilterFlags(false, true));
  write_filtered_data(data_filter);
  return create_xdmf_entry(data_filter, h5_filename, cur_time, comm);
}



template <int dim, int spacedim>
XDMFEntry DataOutInterface<dim,spacedim>::
create_xdmf_entry (const DataOutBase::DataOutFilter &data_filter,
                   const std::string &h5_filename, const double cur_time, MPI_Comm comm) const
{
  return create_xdmf_entry(data_filter, h5_filename, h5_filename, cur_time, comm);
}



template <int dim, int spacedim>
XDMFEntry DataOutInterface<dim,spacedim>::
create_xdmf_entry (const DataOutBase::DataOutFilter &data_filter,
                   const std::string &h5_mesh_filename,
                   const std::string &h5_solution_filename,
                   const double cur_time,
                   MPI_Comm comm) const
{
  unsigned int    local_node_cell_count[2], global_node_cell_count[2];
  int             myrank;

#ifndef DEAL_II_WITH_HDF5
  // throw an exception, but first make
  // sure the compiler does not warn about
  // the now unused function arguments
  (void)data_filter;
  (void)h5_mesh_filename;
  (void)h5_solution_filename;
  (void)cur_time;
  (void)comm;
  AssertThrow(false, ExcMessage ("XDMF support requires HDF5 to be turned on."));
#endif
  AssertThrow(dim == 2 || dim == 3, ExcMessage ("XDMF only supports 2 or 3 dimensions."));

  local_node_cell_count[0] = data_filter.n_nodes();
  local_node_cell_count[1] = data_filter.n_cells();

  // And compute the global total
#ifdef DEAL_II_WITH_MPI
  MPI_Comm_rank(comm, &myrank);
  MPI_Allreduce(local_node_cell_count, global_node_cell_count, 2, MPI_UNSIGNED, MPI_SUM, comm);
#else
  myrank = 0;
  global_node_cell_count[0] = local_node_cell_count[0];
  global_node_cell_count[1] = local_node_cell_count[1];
#endif

  // Output the XDMF file only on the root process
  if (myrank == 0)
    {
      XDMFEntry       entry(h5_mesh_filename, h5_solution_filename, cur_time, global_node_cell_count[0], global_node_cell_count[1], dim);
      unsigned int  n_data_sets = data_filter.n_data_sets();

      // The vector names generated here must match those generated in the HDF5 file
      unsigned int    i;
      for (i=0; i<n_data_sets; ++i)
        {
          entry.add_attribute(data_filter.get_data_set_name(i), data_filter.get_data_set_dim(i));
        }

      return entry;
    }
  else
    {
      return XDMFEntry();
    }
}

template <int dim, int spacedim>
void DataOutInterface<dim,spacedim>::
write_xdmf_file (const std::vector<XDMFEntry> &entries,
                 const std::string &filename,
                 MPI_Comm comm) const
{
  int             myrank;

#ifdef DEAL_II_WITH_MPI
  MPI_Comm_rank(comm, &myrank);
#else
  (void)comm;
  myrank = 0;
#endif

  // Only rank 0 process writes the XDMF file
  if (myrank == 0)
    {
      std::ofstream                               xdmf_file(filename.c_str());
      std::vector<XDMFEntry>::const_iterator      it;

      xdmf_file << "<?xml version=\"1.0\" ?>\n";
      xdmf_file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
      xdmf_file << "<Xdmf Version=\"2.0\">\n";
      xdmf_file << "  <Domain>\n";
      xdmf_file << "    <Grid Name=\"CellTime\" GridType=\"Collection\" CollectionType=\"Temporal\">\n";

      // Write out all the entries indented
      for (it=entries.begin(); it!=entries.end(); ++it)
        xdmf_file << it->get_xdmf_content(3);

      xdmf_file << "    </Grid>\n";
      xdmf_file << "  </Domain>\n";
      xdmf_file << "</Xdmf>\n";

      xdmf_file.close();
    }
}


/*
 * Get the XDMF content associated with this entry.
 * If the entry is not valid, this returns an empty string.
 */
std::string XDMFEntry::get_xdmf_content(const unsigned int indent_level) const
{
  std::stringstream   ss;
  std::map<std::string, unsigned int>::const_iterator     it;

  if (!valid) return "";

  ss << indent(indent_level+0) << "<Grid Name=\"mesh\" GridType=\"Uniform\">\n";
  ss << indent(indent_level+1) << "<Time Value=\"" << entry_time << "\"/>\n";
  ss << indent(indent_level+1) << "<Geometry GeometryType=\"" << (dimension == 2 ? "XY" : "XYZ" ) << "\">\n";
  ss << indent(indent_level+2) << "<DataItem Dimensions=\"" << num_nodes << " " << dimension << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n";
  ss << indent(indent_level+3) << h5_mesh_filename << ":/nodes\n";
  ss << indent(indent_level+2) << "</DataItem>\n";
  ss << indent(indent_level+1) << "</Geometry>\n";
  // If we have cells defined, use a quadrilateral (2D) or hexahedron (3D) topology
  if (num_cells > 0)
    {
      ss << indent(indent_level+1) << "<Topology TopologyType=\"" << (dimension == 2 ? "Quadrilateral" : "Hexahedron") << "\" NumberOfElements=\"" << num_cells << "\">\n";
      ss << indent(indent_level+2) << "<DataItem Dimensions=\"" << num_cells << " " << (2 << (dimension-1)) << "\" NumberType=\"UInt\" Format=\"HDF\">\n";
      ss << indent(indent_level+3) << h5_mesh_filename << ":/cells\n";
      ss << indent(indent_level+2) << "</DataItem>\n";
      ss << indent(indent_level+1) << "</Topology>\n";
    }
  else
    {
      // Otherwise, we assume the points are isolated in space and use a Polyvertex topology
      ss << indent(indent_level+1) << "<Topology TopologyType=\"Polyvertex\" NumberOfElements=\"" << num_nodes << "\">\n";
      ss << indent(indent_level+1) << "</Topology>\n";
    }

  for (it=attribute_dims.begin(); it!=attribute_dims.end(); ++it)
    {
      ss << indent(indent_level+1) << "<Attribute Name=\"" << it->first << "\" AttributeType=\"" << (it->second > 1 ? "Vector" : "Scalar") << "\" Center=\"Node\">\n";
      // Vectors must have 3 elements even for 2D models
      ss << indent(indent_level+2) << "<DataItem Dimensions=\"" << num_nodes << " " << (it->second > 1 ? 3 : 1) << "\" NumberType=\"Float\" Precision=\"8\" Format=\"HDF\">\n";
      ss << indent(indent_level+3) << h5_sol_filename << ":/" << it->first << "\n";
      ss << indent(indent_level+2) << "</DataItem>\n";
      ss << indent(indent_level+1) << "</Attribute>\n";
    }

  ss << indent(indent_level+0) << "</Grid>\n";

  return ss.str();
}

/*
 * Write the data in this DataOutInterface to a DataOutFilter object.
 * Filtering is performed based on the DataOutFilter flags.
 */
template <int dim, int spacedim>
void DataOutInterface<dim,spacedim>::
write_filtered_data (DataOutBase::DataOutFilter &filtered_data) const
{
  DataOutBase::write_filtered_data(get_patches(), get_dataset_names(),
                                   get_vector_data_ranges(),
                                   filtered_data);
}

template <int dim, int spacedim>
void DataOutBase::write_filtered_data (const std::vector<Patch<dim,spacedim> > &patches,
                                       const std::vector<std::string>          &data_names,
                                       const std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> > &vector_data_ranges,
                                       DataOutBase::DataOutFilter &filtered_data)
{
  const unsigned int n_data_sets = data_names.size();
  unsigned int    n_node, n_cell;
  Table<2,double> data_vectors;
  Threads::Task<> reorder_task;

#ifndef DEAL_II_WITH_MPI
  // verify that there are indeed
  // patches to be written out. most
  // of the times, people just forget
  // to call build_patches when there
  // are no patches, so a warning is
  // in order. that said, the
  // assertion is disabled if we
  // support MPI since then it can
  // happen that on the coarsest
  // mesh, a processor simply has no
  // cells it actually owns, and in
  // that case it is legit if there
  // are no patches
  Assert (patches.size() > 0, ExcNoPatches());
#endif

  compute_sizes<dim,spacedim>(patches, n_node, n_cell);

  data_vectors = Table<2,double> (n_data_sets, n_node);
  void (*fun_ptr) (const std::vector<Patch<dim,spacedim> > &, Table<2,double> &) = &DataOutBase::template write_gmv_reorder_data_vectors<dim,spacedim>;
  reorder_task = Threads::new_task (fun_ptr, patches, data_vectors);

  // Write the nodes/cells to the DataOutFilter object.
  write_nodes(patches, filtered_data);
  write_cells(patches, filtered_data);

  // Ensure reordering is done before we output data set values
  reorder_task.join ();

  // when writing, first write out
  // all vector data, then handle the
  // scalar data sets that have been
  // left over
  unsigned int    i, n_th_vector, data_set, pt_data_vector_dim;
  std::string     vector_name;
  for (n_th_vector=0,data_set=0; data_set<n_data_sets;)
    {
      // Advance n_th_vector to at least the current data set we are on
      while (n_th_vector < vector_data_ranges.size() && std_cxx11::get<0>(vector_data_ranges[n_th_vector]) < data_set) n_th_vector++;

      // Determine the dimension of this data
      if (n_th_vector < vector_data_ranges.size() && std_cxx11::get<0>(vector_data_ranges[n_th_vector]) == data_set)
        {
          // Multiple dimensions
          pt_data_vector_dim = std_cxx11::get<1>(vector_data_ranges[n_th_vector]) - std_cxx11::get<0>(vector_data_ranges[n_th_vector])+1;

          // Ensure the dimensionality of the data is correct
          AssertThrow (std_cxx11::get<1>(vector_data_ranges[n_th_vector]) >= std_cxx11::get<0>(vector_data_ranges[n_th_vector]),
                       ExcLowerRange (std_cxx11::get<1>(vector_data_ranges[n_th_vector]), std_cxx11::get<0>(vector_data_ranges[n_th_vector])));
          AssertThrow (std_cxx11::get<1>(vector_data_ranges[n_th_vector]) < n_data_sets,
                       ExcIndexRange (std_cxx11::get<1>(vector_data_ranges[n_th_vector]), 0, n_data_sets));

          // Determine the vector name
          // Concatenate all the
          // component names with double
          // underscores unless a vector
          // name has been specified
          if (std_cxx11::get<2>(vector_data_ranges[n_th_vector]) != "")
            {
              vector_name = std_cxx11::get<2>(vector_data_ranges[n_th_vector]);
            }
          else
            {
              vector_name = "";
              for (i=std_cxx11::get<0>(vector_data_ranges[n_th_vector]); i<std_cxx11::get<1>(vector_data_ranges[n_th_vector]); ++i)
                vector_name += data_names[i] + "__";
              vector_name += data_names[std_cxx11::get<1>(vector_data_ranges[n_th_vector])];
            }
        }
      else
        {
          // One dimension
          pt_data_vector_dim = 1;
          vector_name = data_names[data_set];
        }

      // Write data to the filter object
      filtered_data.write_data_set(vector_name, pt_data_vector_dim, data_set, data_vectors);

      // Advance the current data set
      data_set += pt_data_vector_dim;
    }
}

template <int dim, int spacedim>
void DataOutInterface<dim,spacedim>::
write_hdf5_parallel (const std::string &filename, MPI_Comm comm) const
{
  DataOutBase::DataOutFilter  data_filter(DataOutBase::DataOutFilterFlags(false, true));
  write_filtered_data(data_filter);
  write_hdf5_parallel(data_filter, filename, comm);
}

template <int dim, int spacedim>
void DataOutInterface<dim,spacedim>::
write_hdf5_parallel (const DataOutBase::DataOutFilter &data_filter,
                     const std::string &filename, MPI_Comm comm) const
{
  DataOutBase::write_hdf5_parallel(get_patches(), data_filter, filename, comm);
}

template <int dim, int spacedim>
void DataOutInterface<dim,spacedim>::
write_hdf5_parallel (const DataOutBase::DataOutFilter &data_filter,
                     const bool write_mesh_file, const std::string &mesh_filename, const std::string &solution_filename, MPI_Comm comm) const
{
  DataOutBase::write_hdf5_parallel(get_patches(), data_filter, write_mesh_file, mesh_filename, solution_filename, comm);
}

template <int dim, int spacedim>
void DataOutBase::write_hdf5_parallel (const std::vector<Patch<dim,spacedim> > &patches,
                                       const DataOutBase::DataOutFilter &data_filter,
                                       const std::string &filename,
                                       MPI_Comm comm)
{
  write_hdf5_parallel(patches, data_filter, true, filename, filename, comm);
}

template <int dim, int spacedim>
void DataOutBase::write_hdf5_parallel (const std::vector<Patch<dim,spacedim> > &patches,
                                       const DataOutBase::DataOutFilter &data_filter,
                                       const bool write_mesh_file,
                                       const std::string &mesh_filename,
                                       const std::string &solution_filename,
                                       MPI_Comm comm)
{
#ifndef DEAL_II_WITH_HDF5
  // throw an exception, but first make
  // sure the compiler does not warn about
  // the now unused function arguments
  (void)patches;
  (void)data_filter;
  (void)write_mesh_file;
  (void)mesh_filename;
  (void)solution_filename;
  (void)comm;
  AssertThrow(false, ExcMessage ("HDF5 support is disabled."));
#else
#ifndef DEAL_II_COMPILER_SUPPORTS_MPI
  // verify that there are indeed
  // patches to be written out. most
  // of the times, people just forget
  // to call build_patches when there
  // are no patches, so a warning is
  // in order. that said, the
  // assertion is disabled if we
  // support MPI since then it can
  // happen that on the coarsest
  // mesh, a processor simply has no
  // cells it actually owns, and in
  // that case it is legit if there
  // are no patches
  Assert (data_filter.n_nodes() > 0, ExcNoPatches());
#else
  hid_t           h5_mesh_file_id=-1, h5_solution_file_id, file_plist_id, plist_id;
  hid_t           node_dataspace, node_dataset, node_file_dataspace, node_memory_dataspace;
  hid_t           cell_dataspace, cell_dataset, cell_file_dataspace, cell_memory_dataspace;
  hid_t           pt_data_dataspace, pt_data_dataset, pt_data_file_dataspace, pt_data_memory_dataspace;
  herr_t          status;
  unsigned int    local_node_cell_count[2], global_node_cell_count[2], global_node_cell_offsets[2];
  hsize_t         count[2], offset[2], node_ds_dim[2], cell_ds_dim[2];
  std::vector<double>          node_data_vec;
  std::vector<unsigned int>    cell_data_vec;

  // If HDF5 is not parallel and we're using multiple processes, abort
#ifndef H5_HAVE_PARALLEL
#  ifdef DEAL_II_WITH_MPI
  int world_size;
  MPI_Comm_size(comm, &world_size);
  AssertThrow (world_size <= 1,
               ExcMessage ("Serial HDF5 output on multiple processes is not yet supported."));
#  endif
#endif

  local_node_cell_count[0] = data_filter.n_nodes();
  local_node_cell_count[1] = data_filter.n_cells();

  // Create file access properties
  file_plist_id = H5Pcreate(H5P_FILE_ACCESS);
  AssertThrow(file_plist_id != -1, ExcIO());
  // If MPI is enabled *and* HDF5 is parallel, we can do parallel output
#ifdef DEAL_II_WITH_MPI
#ifdef H5_HAVE_PARALLEL
  // Set the access to use the specified MPI_Comm object
  status = H5Pset_fapl_mpio(file_plist_id, comm, MPI_INFO_NULL);
  AssertThrow(status >= 0, ExcIO());
#endif
#endif

  // Compute the global total number of nodes/cells
  // And determine the offset of the data for this process
#ifdef DEAL_II_WITH_MPI
  MPI_Allreduce(local_node_cell_count, global_node_cell_count, 2, MPI_UNSIGNED, MPI_SUM, comm);
  MPI_Scan(local_node_cell_count, global_node_cell_offsets, 2, MPI_UNSIGNED, MPI_SUM, comm);
  global_node_cell_offsets[0] -= local_node_cell_count[0];
  global_node_cell_offsets[1] -= local_node_cell_count[1];
#else
  global_node_cell_offsets[0] = global_node_cell_offsets[1] = 0;
#endif

  // Create the property list for a collective write
  plist_id = H5Pcreate(H5P_DATASET_XFER);
  AssertThrow(plist_id >= 0, ExcIO());
#ifdef DEAL_II_WITH_MPI
#ifdef H5_HAVE_PARALLEL
  status = H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
  AssertThrow(status >= 0, ExcIO());
#endif
#endif

  if (write_mesh_file)
    {
      // Overwrite any existing files (change this to an option?)
      h5_mesh_file_id = H5Fcreate(mesh_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, file_plist_id);
      AssertThrow(h5_mesh_file_id >= 0, ExcIO());

      // Create the dataspace for the nodes and cells
      node_ds_dim[0] = global_node_cell_count[0];
      node_ds_dim[1] = dim;
      node_dataspace = H5Screate_simple(2, node_ds_dim, NULL);
      AssertThrow(node_dataspace >= 0, ExcIO());

      cell_ds_dim[0] = global_node_cell_count[1];
      cell_ds_dim[1] = GeometryInfo<dim>::vertices_per_cell;
      cell_dataspace = H5Screate_simple(2, cell_ds_dim, NULL);
      AssertThrow(cell_dataspace >= 0, ExcIO());

      // Create the dataset for the nodes and cells
#if H5Gcreate_vers == 1
      node_dataset = H5Dcreate(h5_mesh_file_id, "nodes", H5T_NATIVE_DOUBLE, node_dataspace, H5P_DEFAULT);
#else
      node_dataset = H5Dcreate(h5_mesh_file_id, "nodes", H5T_NATIVE_DOUBLE, node_dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
      AssertThrow(node_dataset >= 0, ExcIO());
#if H5Gcreate_vers == 1
      cell_dataset = H5Dcreate(h5_mesh_file_id, "cells", H5T_NATIVE_UINT, cell_dataspace, H5P_DEFAULT);
#else
      cell_dataset = H5Dcreate(h5_mesh_file_id, "cells", H5T_NATIVE_UINT, cell_dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
      AssertThrow(cell_dataset >= 0, ExcIO());

      // Close the node and cell dataspaces since we're done with them
      status = H5Sclose(node_dataspace);
      AssertThrow(status >= 0, ExcIO());
      status = H5Sclose(cell_dataspace);
      AssertThrow(status >= 0, ExcIO());

      // Create the data subset we'll use to read from memory
      count[0] = local_node_cell_count[0];
      count[1] = dim;
      offset[0] = global_node_cell_offsets[0];
      offset[1] = 0;
      node_memory_dataspace = H5Screate_simple(2, count, NULL);
      AssertThrow(node_memory_dataspace >= 0, ExcIO());

      // Select the hyperslab in the file
      node_file_dataspace = H5Dget_space(node_dataset);
      AssertThrow(node_file_dataspace >= 0, ExcIO());
      status = H5Sselect_hyperslab(node_file_dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
      AssertThrow(status >= 0, ExcIO());

      // And repeat for cells
      count[0] = local_node_cell_count[1];
      count[1] = GeometryInfo<dim>::vertices_per_cell;
      offset[0] = global_node_cell_offsets[1];
      offset[1] = 0;
      cell_memory_dataspace = H5Screate_simple(2, count, NULL);
      AssertThrow(cell_memory_dataspace >= 0, ExcIO());

      cell_file_dataspace = H5Dget_space(cell_dataset);
      AssertThrow(cell_file_dataspace >= 0, ExcIO());
      status = H5Sselect_hyperslab(cell_file_dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
      AssertThrow(status >= 0, ExcIO());

      // And finally, write the node data
      data_filter.fill_node_data(node_data_vec);
      status = H5Dwrite(node_dataset, H5T_NATIVE_DOUBLE, node_memory_dataspace, node_file_dataspace, plist_id, &node_data_vec[0]);
      AssertThrow(status >= 0, ExcIO());
      node_data_vec.clear();

      // And the cell data
      data_filter.fill_cell_data(global_node_cell_offsets[0], cell_data_vec);
      status = H5Dwrite(cell_dataset, H5T_NATIVE_UINT, cell_memory_dataspace, cell_file_dataspace, plist_id, &cell_data_vec[0]);
      AssertThrow(status >= 0, ExcIO());
      cell_data_vec.clear();

      // Close the file dataspaces
      status = H5Sclose(node_file_dataspace);
      AssertThrow(status >= 0, ExcIO());
      status = H5Sclose(cell_file_dataspace);
      AssertThrow(status >= 0, ExcIO());

      // Close the memory dataspaces
      status = H5Sclose(node_memory_dataspace);
      AssertThrow(status >= 0, ExcIO());
      status = H5Sclose(cell_memory_dataspace);
      AssertThrow(status >= 0, ExcIO());

      // Close the datasets
      status = H5Dclose(node_dataset);
      AssertThrow(status >= 0, ExcIO());
      status = H5Dclose(cell_dataset);
      AssertThrow(status >= 0, ExcIO());

      // If the filenames are different, we need to close the mesh file
      if (mesh_filename != solution_filename)
        {
          status = H5Fclose(h5_mesh_file_id);
          AssertThrow(status >= 0, ExcIO());
        }
    }

  // If the filenames are identical, continue with the same file
  if (mesh_filename == solution_filename && write_mesh_file)
    {
      h5_solution_file_id = h5_mesh_file_id;
    }
  else
    {
      // Otherwise we need to open a new file
      h5_solution_file_id = H5Fcreate(solution_filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, file_plist_id);
      AssertThrow(h5_solution_file_id >= 0, ExcIO());
    }

  // when writing, first write out
  // all vector data, then handle the
  // scalar data sets that have been
  // left over
  unsigned int    i, pt_data_vector_dim;
  std::string     vector_name;
  for (i=0; i<data_filter.n_data_sets(); ++i)
    {
      // Allocate space for the point data
      // Must be either 1D or 3D
      pt_data_vector_dim = data_filter.get_data_set_dim(i);
      vector_name = data_filter.get_data_set_name(i);

      // Create the dataspace for the point data
      node_ds_dim[0] = global_node_cell_count[0];
      node_ds_dim[1] = pt_data_vector_dim;
      pt_data_dataspace = H5Screate_simple(2, node_ds_dim, NULL);
      AssertThrow(pt_data_dataspace >= 0, ExcIO());

#if H5Gcreate_vers == 1
      pt_data_dataset = H5Dcreate(h5_solution_file_id, vector_name.c_str(), H5T_NATIVE_DOUBLE, pt_data_dataspace, H5P_DEFAULT);
#else
      pt_data_dataset = H5Dcreate(h5_solution_file_id, vector_name.c_str(), H5T_NATIVE_DOUBLE, pt_data_dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif
      AssertThrow(pt_data_dataset >= 0, ExcIO());

      // Create the data subset we'll use to read from memory
      count[0] = local_node_cell_count[0];
      count[1] = pt_data_vector_dim;
      offset[0] = global_node_cell_offsets[0];
      offset[1] = 0;
      pt_data_memory_dataspace = H5Screate_simple(2, count, NULL);
      AssertThrow(pt_data_memory_dataspace >= 0, ExcIO());

      // Select the hyperslab in the file
      pt_data_file_dataspace = H5Dget_space(pt_data_dataset);
      AssertThrow(pt_data_file_dataspace >= 0, ExcIO());
      status = H5Sselect_hyperslab(pt_data_file_dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
      AssertThrow(status >= 0, ExcIO());

      // And finally, write the data
      status = H5Dwrite(pt_data_dataset, H5T_NATIVE_DOUBLE, pt_data_memory_dataspace, pt_data_file_dataspace, plist_id, data_filter.get_data_set(i));
      AssertThrow(status >= 0, ExcIO());

      // Close the dataspaces
      status = H5Sclose(pt_data_dataspace);
      AssertThrow(status >= 0, ExcIO());
      status = H5Sclose(pt_data_memory_dataspace);
      AssertThrow(status >= 0, ExcIO());
      status = H5Sclose(pt_data_file_dataspace);
      AssertThrow(status >= 0, ExcIO());
      // Close the dataset
      status = H5Dclose(pt_data_dataset);
      AssertThrow(status >= 0, ExcIO());
    }

  // Close the file property list
  status = H5Pclose(file_plist_id);
  AssertThrow(status >= 0, ExcIO());

  // Close the parallel access
  status = H5Pclose(plist_id);
  AssertThrow(status >= 0, ExcIO());

  // Close the file
  status = H5Fclose(h5_solution_file_id);
  AssertThrow(status >= 0, ExcIO());
#endif
#endif
}



template <int dim, int spacedim>
void
DataOutInterface<dim,spacedim>::write (std::ostream &out,
                                       const DataOutBase::OutputFormat output_format_) const
{
  DataOutBase::OutputFormat output_format = output_format_;
  if (output_format == DataOutBase::default_format)
    output_format = default_fmt;

  switch (output_format)
    {
    case DataOutBase::none:
      break;

    case DataOutBase::dx:
      write_dx (out);
      break;

    case DataOutBase::ucd:
      write_ucd (out);
      break;

    case DataOutBase::gnuplot:
      write_gnuplot (out);
      break;

    case DataOutBase::povray:
      write_povray (out);
      break;

    case DataOutBase::eps:
      write_eps (out);
      break;

    case DataOutBase::gmv:
      write_gmv (out);
      break;

    case DataOutBase::tecplot:
      write_tecplot (out);
      break;

    case DataOutBase::tecplot_binary:
      write_tecplot_binary (out);
      break;

    case DataOutBase::vtk:
      write_vtk (out);
      break;

    case DataOutBase::vtu:
      write_vtu (out);
      break;

    case DataOutBase::svg:
      write_svg (out);
      break;

    case DataOutBase::deal_II_intermediate:
      write_deal_II_intermediate (out);
      break;

    default:
      Assert (false, ExcNotImplemented());
    }
}



template <int dim, int spacedim>
void
DataOutInterface<dim,spacedim>::set_default_format(const DataOutBase::OutputFormat fmt)
{
  Assert (fmt != DataOutBase::default_format, ExcNotImplemented());
  default_fmt = fmt;
}



template <int dim, int spacedim>
void
DataOutInterface<dim,spacedim>::set_flags (const DataOutBase::DXFlags &flags)
{
  dx_flags = flags;
}



template <int dim, int spacedim>
void
DataOutInterface<dim,spacedim>::set_flags (const DataOutBase::UcdFlags &flags)
{
  ucd_flags = flags;
}



template <int dim, int spacedim>
void
DataOutInterface<dim,spacedim>::set_flags (const DataOutBase::GnuplotFlags &flags)
{
  gnuplot_flags = flags;
}



template <int dim, int spacedim>
void
DataOutInterface<dim,spacedim>::set_flags (const DataOutBase::PovrayFlags &flags)
{
  povray_flags = flags;
}



template <int dim, int spacedim>
void
DataOutInterface<dim,spacedim>::set_flags (const DataOutBase::EpsFlags &flags)
{
  eps_flags = flags;
}



template <int dim, int spacedim>
void
DataOutInterface<dim,spacedim>::set_flags (const DataOutBase::GmvFlags &flags)
{
  gmv_flags = flags;
}



template <int dim, int spacedim>
void
DataOutInterface<dim,spacedim>::set_flags (const DataOutBase::TecplotFlags &flags)
{
  tecplot_flags = flags;
}



template <int dim, int spacedim>
void
DataOutInterface<dim,spacedim>::set_flags (const DataOutBase::VtkFlags &flags)
{
  vtk_flags = flags;
}



template <int dim, int spacedim>
void
DataOutInterface<dim,spacedim>::set_flags (const DataOutBase::SvgFlags &flags)
{
  svg_flags = flags;
}



template <int dim, int spacedim>
void
DataOutInterface<dim,spacedim>::set_flags (const DataOutBase::Deal_II_IntermediateFlags &flags)
{
  deal_II_intermediate_flags = flags;
}



template <int dim, int spacedim>
std::string
DataOutInterface<dim,spacedim>::
default_suffix (const DataOutBase::OutputFormat output_format) const
{
  if (output_format == DataOutBase::default_format)
    return DataOutBase::default_suffix (default_fmt);
  else
    return DataOutBase::default_suffix (output_format);
}



template <int dim, int spacedim>
void
DataOutInterface<dim,spacedim>::declare_parameters (ParameterHandler &prm)
{
  prm.declare_entry ("Output format", "gnuplot",
                     Patterns::Selection (DataOutBase::get_output_format_names ()),
                     "A name for the output format to be used");
  prm.declare_entry("Subdivisions", "1", Patterns::Integer(),
                    "Number of subdivisions of each mesh cell");

  prm.enter_subsection ("DX output parameters");
  DataOutBase::DXFlags::declare_parameters (prm);
  prm.leave_subsection ();

  prm.enter_subsection ("UCD output parameters");
  DataOutBase::UcdFlags::declare_parameters (prm);
  prm.leave_subsection ();

  prm.enter_subsection ("Gnuplot output parameters");
  DataOutBase::GnuplotFlags::declare_parameters (prm);
  prm.leave_subsection ();

  prm.enter_subsection ("Povray output parameters");
  DataOutBase::PovrayFlags::declare_parameters (prm);
  prm.leave_subsection ();

  prm.enter_subsection ("Eps output parameters");
  DataOutBase::EpsFlags::declare_parameters (prm);
  prm.leave_subsection ();

  prm.enter_subsection ("Gmv output parameters");
  DataOutBase::GmvFlags::declare_parameters (prm);
  prm.leave_subsection ();

  prm.enter_subsection ("Tecplot output parameters");
  DataOutBase::TecplotFlags::declare_parameters (prm);
  prm.leave_subsection ();

  prm.enter_subsection ("Vtk output parameters");
  DataOutBase::VtkFlags::declare_parameters (prm);
  prm.leave_subsection ();


  prm.enter_subsection ("deal.II intermediate output parameters");
  DataOutBase::Deal_II_IntermediateFlags::declare_parameters (prm);
  prm.leave_subsection ();
}



template <int dim, int spacedim>
void
DataOutInterface<dim,spacedim>::parse_parameters (ParameterHandler &prm)
{
  const std::string &output_name = prm.get ("Output format");
  default_fmt = DataOutBase::parse_output_format (output_name);
  default_subdivisions = prm.get_integer ("Subdivisions");

  prm.enter_subsection ("DX output parameters");
  dx_flags.parse_parameters (prm);
  prm.leave_subsection ();

  prm.enter_subsection ("UCD output parameters");
  ucd_flags.parse_parameters (prm);
  prm.leave_subsection ();

  prm.enter_subsection ("Gnuplot output parameters");
  gnuplot_flags.parse_parameters (prm);
  prm.leave_subsection ();

  prm.enter_subsection ("Povray output parameters");
  povray_flags.parse_parameters (prm);
  prm.leave_subsection ();

  prm.enter_subsection ("Eps output parameters");
  eps_flags.parse_parameters (prm);
  prm.leave_subsection ();

  prm.enter_subsection ("Gmv output parameters");
  gmv_flags.parse_parameters (prm);
  prm.leave_subsection ();

  prm.enter_subsection ("Tecplot output parameters");
  tecplot_flags.parse_parameters (prm);
  prm.leave_subsection ();

  prm.enter_subsection ("Vtk output parameters");
  vtk_flags.parse_parameters (prm);
  prm.leave_subsection ();

  prm.enter_subsection ("deal.II intermediate output parameters");
  deal_II_intermediate_flags.parse_parameters (prm);
  prm.leave_subsection ();
}



template <int dim, int spacedim>
std::size_t
DataOutInterface<dim,spacedim>::memory_consumption () const
{
  return (sizeof (default_fmt) +
          MemoryConsumption::memory_consumption (dx_flags) +
          MemoryConsumption::memory_consumption (ucd_flags) +
          MemoryConsumption::memory_consumption (gnuplot_flags) +
          MemoryConsumption::memory_consumption (povray_flags) +
          MemoryConsumption::memory_consumption (eps_flags) +
          MemoryConsumption::memory_consumption (gmv_flags) +
          MemoryConsumption::memory_consumption (tecplot_flags) +
          MemoryConsumption::memory_consumption (vtk_flags) +
          MemoryConsumption::memory_consumption (svg_flags) +
          MemoryConsumption::memory_consumption (deal_II_intermediate_flags));
}



template <int dim, int spacedim>
std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> >
DataOutInterface<dim,spacedim>::get_vector_data_ranges () const
{
  return std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> >();
}




// ---------------------------------------------- DataOutReader ----------

template <int dim, int spacedim>
void
DataOutReader<dim,spacedim>::read (std::istream &in)
{
  Assert (in, ExcIO());

  // first empty previous content
  {
    std::vector<typename dealii::DataOutBase::Patch<dim,spacedim> >
    tmp;
    tmp.swap (patches);
  }
  {
    std::vector<std::string> tmp;
    tmp.swap (dataset_names);
  }
  {
    std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> > tmp;
    tmp.swap (vector_data_ranges);
  }

  // then check that we have the
  // correct header of this
  // file. both the first and second
  // real lines have to match, as
  // well as the dimension
  // information written before that
  // and the Version information
  // written in the third line
  {
    std::pair<unsigned int, unsigned int>
    dimension_info
      = DataOutBase::determine_intermediate_format_dimensions (in);
    AssertThrow ((dimension_info.first  == dim) &&
                 (dimension_info.second == spacedim),
                 ExcIncompatibleDimensions (dimension_info.first, dim,
                                            dimension_info.second, spacedim));

    // read to the end of the line
    std::string tmp;
    getline (in, tmp);
  }

  {
    std::string header;
    getline (in, header);

    std::ostringstream s;
    s << "[deal.II intermediate format graphics data]";

    Assert (header == s.str(), ExcUnexpectedInput(s.str(),header));
  }
  {
    std::string header;
    getline (in, header);

    std::ostringstream s;
    s << "[written by " << DEAL_II_PACKAGE_NAME << " " << DEAL_II_PACKAGE_VERSION << "]";

    Assert (header == s.str(), ExcUnexpectedInput(s.str(),header));
  }
  {
    std::string header;
    getline (in, header);

    std::ostringstream s;
    s << "[Version: " << dealii::DataOutBase::Deal_II_IntermediateFlags::format_version << "]";

    Assert (header == s.str(),
            ExcMessage("Invalid or incompatible file format. Intermediate format "
                       "files can only be read by the same deal.II version as they "
                       "are written by."));
  }

  // then read the rest of the data
  unsigned int n_datasets;
  in >> n_datasets;
  dataset_names.resize (n_datasets);
  for (unsigned int i=0; i<n_datasets; ++i)
    in >> dataset_names[i];

  unsigned int n_patches;
  in >> n_patches;
  patches.resize (n_patches);
  for (unsigned int i=0; i<n_patches; ++i)
    in >> patches[i];

  unsigned int n_vector_data_ranges;
  in >> n_vector_data_ranges;
  vector_data_ranges.resize (n_vector_data_ranges);
  for (unsigned int i=0; i<n_vector_data_ranges; ++i)
    {
      in >> std_cxx11::get<0>(vector_data_ranges[i])
         >> std_cxx11::get<1>(vector_data_ranges[i]);

      // read in the name of that vector
      // range. because it is on a separate
      // line, we first need to read to the
      // end of the previous line (nothing
      // should be there any more after we've
      // read the previous two integers) and
      // then read the entire next line for
      // the name
      std::string name;
      getline(in, name);
      getline(in, name);
      std_cxx11::get<2>(vector_data_ranges[i]) = name;
    }

  Assert (in, ExcIO());
}



template <int dim, int spacedim>
void
DataOutReader<dim,spacedim>::
merge (const DataOutReader<dim,spacedim> &source)
{
  typedef typename dealii::DataOutBase::Patch<dim,spacedim> Patch;


  const std::vector<Patch> source_patches = source.get_patches ();
  Assert (patches.size () != 0,        ExcNoPatches ());
  Assert (source_patches.size () != 0, ExcNoPatches ());
  // check equality of component
  // names
  Assert (get_dataset_names() == source.get_dataset_names(),
          ExcIncompatibleDatasetNames());

  // check equality of the vector data
  // specifications
  Assert (get_vector_data_ranges().size() ==
          source.get_vector_data_ranges().size(),
          ExcMessage ("Both sources need to declare the same components "
                      "as vectors."));
  for (unsigned int i=0; i<get_vector_data_ranges().size(); ++i)
    {
      Assert (std_cxx11::get<0>(get_vector_data_ranges()[i]) ==
              std_cxx11::get<0>(source.get_vector_data_ranges()[i]),
              ExcMessage ("Both sources need to declare the same components "
                          "as vectors."));
      Assert (std_cxx11::get<1>(get_vector_data_ranges()[i]) ==
              std_cxx11::get<1>(source.get_vector_data_ranges()[i]),
              ExcMessage ("Both sources need to declare the same components "
                          "as vectors."));
      Assert (std_cxx11::get<2>(get_vector_data_ranges()[i]) ==
              std_cxx11::get<2>(source.get_vector_data_ranges()[i]),
              ExcMessage ("Both sources need to declare the same components "
                          "as vectors."));
    }

  // make sure patches are compatible
  Assert (patches[0].n_subdivisions == source_patches[0].n_subdivisions,
          ExcIncompatiblePatchLists());
  Assert (patches[0].data.n_rows() == source_patches[0].data.n_rows(),
          ExcIncompatiblePatchLists());
  Assert (patches[0].data.n_cols() == source_patches[0].data.n_cols(),
          ExcIncompatiblePatchLists());

  // merge patches. store old number
  // of elements, since we need to
  // adjust patch numbers, etc
  // afterwards
  const unsigned int old_n_patches = patches.size();
  patches.insert (patches.end(),
                  source_patches.begin(),
                  source_patches.end());

  // adjust patch numbers
  for (unsigned int i=old_n_patches; i<patches.size(); ++i)
    patches[i].patch_index += old_n_patches;

  // adjust patch neighbors
  for (unsigned int i=old_n_patches; i<patches.size(); ++i)
    for (unsigned int n=0; n<GeometryInfo<dim>::faces_per_cell; ++n)
      if (patches[i].neighbors[n] != Patch::no_neighbor)
        patches[i].neighbors[n] += old_n_patches;
}



template <int dim, int spacedim>
const std::vector<typename dealii::DataOutBase::Patch<dim,spacedim> > &
DataOutReader<dim,spacedim>::get_patches () const
{
  return patches;
}



template <int dim, int spacedim>
std::vector<std::string>
DataOutReader<dim,spacedim>::get_dataset_names () const
{
  return dataset_names;
}



template <int dim, int spacedim>
std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> >
DataOutReader<dim,spacedim>::get_vector_data_ranges () const
{
  return vector_data_ranges;
}



namespace DataOutBase
{
  template <int dim, int spacedim>
  std::ostream &
  operator << (std::ostream                           &out,
               const Patch<dim,spacedim> &patch)
  {
    // write a header line
    out << "[deal.II intermediate Patch<" << dim << ',' << spacedim << ">]"
        << '\n';

    // then write all the data that is
    // in this patch
    for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
      out << patch.vertices[GeometryInfo<dim>::ucd_to_deal[i]] << ' ';
    out << '\n';

    for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
      out << patch.neighbors[i] << ' ';
    out << '\n';

    out << patch.patch_index << ' ' << patch.n_subdivisions
        << '\n';

    out << patch.points_are_available<<'\n';

    out << patch.data.n_rows() << ' ' << patch.data.n_cols() << '\n';
    for (unsigned int i=0; i<patch.data.n_rows(); ++i)
      for (unsigned int j=0; j<patch.data.n_cols(); ++j)
        out << patch.data[i][j] << ' ';
    out << '\n';
    out << '\n';

    return out;
  }


  template <int dim, int spacedim>
  std::istream &
  operator >> (std::istream                     &in,
               Patch<dim,spacedim> &patch)
  {
    Assert (in, ExcIO());

    // read a header line and compare
    // it to what we usually
    // write. skip all lines that
    // contain only blanks at the start
    {
      std::string header;
      do
        {
          getline (in, header);
          while ((header.size() != 0) &&
                 (header[header.size()-1] == ' '))
            header.erase(header.size()-1);
        }
      while ((header == "") && in);

      std::ostringstream s;
      s << "[deal.II intermediate Patch<" << dim << ',' << spacedim << ">]";

      Assert (header == s.str(), ExcUnexpectedInput(s.str(),header));
    }


    // then read all the data that is
    // in this patch
    for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
      in >> patch.vertices[GeometryInfo<dim>::ucd_to_deal[i]];

    for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
      in >> patch.neighbors[i];

    in >> patch.patch_index >> patch.n_subdivisions;

    in >> patch.points_are_available;

    unsigned int n_rows, n_cols;
    in >> n_rows >> n_cols;
    patch.data.reinit (n_rows, n_cols);
    for (unsigned int i=0; i<patch.data.n_rows(); ++i)
      for (unsigned int j=0; j<patch.data.n_cols(); ++j)
        in >> patch.data[i][j];

    Assert (in, ExcIO());

    return in;
  }
}



// explicit instantiations
#include "data_out_base.inst"

DEAL_II_NAMESPACE_CLOSE
