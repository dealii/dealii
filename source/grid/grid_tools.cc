// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2016 by the deal.II authors
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

#include <deal.II/base/std_cxx11/array.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/filtered_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/distributed/tria_base.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/numerics/matrix_tools.h>

#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <cmath>
#include <numeric>
#include <list>
#include <set>


DEAL_II_NAMESPACE_OPEN


namespace GridTools
{

  template <int dim, int spacedim>
  double
  diameter (const Triangulation<dim, spacedim> &tria)
  {
    // we can't deal with distributed meshes since we don't have all
    // vertices locally. there is one exception, however: if the mesh has
    // never been refined. the way to test this is not to ask
    // tria.n_levels()==1, since this is something that can happen on one
    // processor without being true on all. however, we can ask for the
    // global number of active cells and use that
#if defined(DEAL_II_WITH_P4EST) && defined(DEBUG)
    if (const parallel::distributed::Triangulation<dim,spacedim> *p_tria
        = dynamic_cast<const parallel::distributed::Triangulation<dim,spacedim>*>(&tria))
      Assert (p_tria->n_global_active_cells() == tria.n_cells(0),
              ExcNotImplemented());
#endif

    // the algorithm used simply traverses all cells and picks out the
    // boundary vertices. it may or may not be faster to simply get all
    // vectors, don't mark boundary vertices, and compute the distances
    // thereof, but at least as the mesh is refined, it seems better to
    // first mark boundary nodes, as marking is O(N) in the number of
    // cells/vertices, while computing the maximal distance is O(N*N)
    const std::vector<Point<spacedim> > &vertices = tria.get_vertices ();
    std::vector<bool> boundary_vertices (vertices.size(), false);

    typename Triangulation<dim,spacedim>::active_cell_iterator
    cell = tria.begin_active();
    const typename Triangulation<dim,spacedim>::active_cell_iterator
    endc = tria.end();
    for (; cell!=endc; ++cell)
      for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
        if (cell->face(face)->at_boundary ())
          for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_face; ++i)
            boundary_vertices[cell->face(face)->vertex_index(i)] = true;

    // now traverse the list of boundary vertices and check distances.
    // since distances are symmetric, we only have to check one half
    double max_distance_sqr = 0;
    std::vector<bool>::const_iterator pi = boundary_vertices.begin();
    const unsigned int N = boundary_vertices.size();
    for (unsigned int i=0; i<N; ++i, ++pi)
      {
        std::vector<bool>::const_iterator pj = pi+1;
        for (unsigned int j=i+1; j<N; ++j, ++pj)
          if ((*pi==true) && (*pj==true) &&
              ((vertices[i]-vertices[j]).norm_square() > max_distance_sqr))
            max_distance_sqr = (vertices[i]-vertices[j]).norm_square();
      };

    return std::sqrt(max_distance_sqr);
  }



  template <int dim, int spacedim>
  double
  volume (const Triangulation<dim, spacedim> &triangulation,
          const Mapping<dim,spacedim> &mapping)
  {
    // get the degree of the mapping if possible. if not, just assume 1
    const unsigned int mapping_degree
      = (dynamic_cast<const MappingQ<dim,spacedim>*>(&mapping) != 0 ?
         dynamic_cast<const MappingQ<dim,spacedim>*>(&mapping)->get_degree() :
         1);

    // then initialize an appropriate quadrature formula
    const QGauss<dim> quadrature_formula (mapping_degree + 1);
    const unsigned int n_q_points = quadrature_formula.size();

    // we really want the JxW values from the FEValues object, but it
    // wants a finite element. create a cheap element as a dummy
    // element
    FE_Nothing<dim,spacedim> dummy_fe;
    FEValues<dim,spacedim> fe_values (mapping, dummy_fe, quadrature_formula,
                                      update_JxW_values);

    typename Triangulation<dim,spacedim>::active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();

    double local_volume = 0;

    // compute the integral quantities by quadrature
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          for (unsigned int q=0; q<n_q_points; ++q)
            local_volume += fe_values.JxW(q);
        }

    double global_volume = 0;

#ifdef DEAL_II_WITH_MPI
    if (const parallel::Triangulation<dim,spacedim> *p_tria
        = dynamic_cast<const parallel::Triangulation<dim,spacedim>*>(&triangulation))
      global_volume = Utilities::MPI::sum (local_volume, p_tria->get_communicator());
    else
#endif
      global_volume = local_volume;

    return global_volume;
  }


  template <>
  double
  cell_measure<3>(const std::vector<Point<3> > &all_vertices,
                  const unsigned int (&vertex_indices)[GeometryInfo<3>::vertices_per_cell])
  {
    // note that this is the
    // cell_measure based on the new
    // deal.II numbering. When called
    // from inside GridReordering make
    // sure that you reorder the
    // vertex_indices before
    const double x[8] = { all_vertices[vertex_indices[0]](0),
                          all_vertices[vertex_indices[1]](0),
                          all_vertices[vertex_indices[2]](0),
                          all_vertices[vertex_indices[3]](0),
                          all_vertices[vertex_indices[4]](0),
                          all_vertices[vertex_indices[5]](0),
                          all_vertices[vertex_indices[6]](0),
                          all_vertices[vertex_indices[7]](0)
                        };
    const double y[8] = { all_vertices[vertex_indices[0]](1),
                          all_vertices[vertex_indices[1]](1),
                          all_vertices[vertex_indices[2]](1),
                          all_vertices[vertex_indices[3]](1),
                          all_vertices[vertex_indices[4]](1),
                          all_vertices[vertex_indices[5]](1),
                          all_vertices[vertex_indices[6]](1),
                          all_vertices[vertex_indices[7]](1)
                        };
    const double z[8] = { all_vertices[vertex_indices[0]](2),
                          all_vertices[vertex_indices[1]](2),
                          all_vertices[vertex_indices[2]](2),
                          all_vertices[vertex_indices[3]](2),
                          all_vertices[vertex_indices[4]](2),
                          all_vertices[vertex_indices[5]](2),
                          all_vertices[vertex_indices[6]](2),
                          all_vertices[vertex_indices[7]](2)
                        };

    /*
      This is the same Maple script as in the barycenter method above
      except of that here the shape functions tphi[0]-tphi[7] are ordered
      according to the lexicographic numbering.

      x := array(0..7):
      y := array(0..7):
      z := array(0..7):
      tphi[0] := (1-xi)*(1-eta)*(1-zeta):
      tphi[1] :=     xi*(1-eta)*(1-zeta):
      tphi[2] := (1-xi)*    eta*(1-zeta):
      tphi[3] :=     xi*    eta*(1-zeta):
      tphi[4] := (1-xi)*(1-eta)*zeta:
      tphi[5] :=     xi*(1-eta)*zeta:
      tphi[6] := (1-xi)*    eta*zeta:
      tphi[7] :=     xi*    eta*zeta:
      x_real := sum(x[s]*tphi[s], s=0..7):
      y_real := sum(y[s]*tphi[s], s=0..7):
      z_real := sum(z[s]*tphi[s], s=0..7):
      with (linalg):
      J := matrix(3,3, [[diff(x_real, xi), diff(x_real, eta), diff(x_real, zeta)],
      [diff(y_real, xi), diff(y_real, eta), diff(y_real, zeta)],
      [diff(z_real, xi), diff(z_real, eta), diff(z_real, zeta)]]):
      detJ := det (J):

      measure := simplify ( int ( int ( int (detJ, xi=0..1), eta=0..1), zeta=0..1)):

      readlib(C):

      C(measure, optimized);

      The C code produced by this maple script is further optimized by
      hand. In particular, division by 12 is performed only once, not
      hundred of times.
    */

    const double t3 = y[3]*x[2];
    const double t5 = z[1]*x[5];
    const double t9 = z[3]*x[2];
    const double t11 = x[1]*y[0];
    const double t14 = x[4]*y[0];
    const double t18 = x[5]*y[7];
    const double t20 = y[1]*x[3];
    const double t22 = y[5]*x[4];
    const double t26 = z[7]*x[6];
    const double t28 = x[0]*y[4];
    const double t34 = z[3]*x[1]*y[2]+t3*z[1]-t5*y[7]+y[7]*x[4]*z[6]+t9*y[6]-t11*z[4]-t5*y[3]-t14*z[2]+z[1]*x[4]*y[0]-t18*z[3]+t20*z[0]-t22*z[0]-y[0]*x[5]*z[4]-t26*y[3]+t28*z[2]-t9*y[1]-y[1]*x[4]*z[0]-t11*z[5];
    const double t37 = y[1]*x[0];
    const double t44 = x[1]*y[5];
    const double t46 = z[1]*x[0];
    const double t49 = x[0]*y[2];
    const double t52 = y[5]*x[7];
    const double t54 = x[3]*y[7];
    const double t56 = x[2]*z[0];
    const double t58 = x[3]*y[2];
    const double t64 = -x[6]*y[4]*z[2]-t37*z[2]+t18*z[6]-x[3]*y[6]*z[2]+t11*z[2]+t5*y[0]+t44*z[4]-t46*y[4]-t20*z[7]-t49*z[6]-t22*z[1]+t52*z[3]-t54*z[2]-t56*y[4]-t58*z[0]+y[1]*x[2]*z[0]+t9*y[7]+t37*z[4];
    const double t66 = x[1]*y[7];
    const double t68 = y[0]*x[6];
    const double t70 = x[7]*y[6];
    const double t73 = z[5]*x[4];
    const double t76 = x[6]*y[7];
    const double t90 = x[4]*z[0];
    const double t92 = x[1]*y[3];
    const double t95 = -t66*z[3]-t68*z[2]-t70*z[2]+t26*y[5]-t73*y[6]-t14*z[6]+t76*z[2]-t3*z[6]+x[6]*y[2]*z[4]-z[3]*x[6]*y[2]+t26*y[4]-t44*z[3]-x[1]*y[2]*z[0]+x[5]*y[6]*z[4]+t54*z[5]+t90*y[2]-t92*z[2]+t46*y[2];
    const double t102 = x[2]*y[0];
    const double t107 = y[3]*x[7];
    const double t114 = x[0]*y[6];
    const double t125 = y[0]*x[3]*z[2]-z[7]*x[5]*y[6]-x[2]*y[6]*z[4]+t102*z[6]-t52*z[6]+x[2]*y[4]*z[6]-t107*z[5]-t54*z[6]+t58*z[6]-x[7]*y[4]*z[6]+t37*z[5]-t114*z[4]+t102*z[4]-z[1]*x[2]*y[0]+t28*z[6]-y[5]*x[6]*z[4]-z[5]*x[1]*y[4]-t73*y[7];
    const double t129 = z[0]*x[6];
    const double t133 = y[1]*x[7];
    const double t145 = y[1]*x[5];
    const double t156 = t90*y[6]-t129*y[4]+z[7]*x[2]*y[6]-t133*z[5]+x[5]*y[3]*z[7]-t26*y[2]-t70*z[3]+t46*y[3]+z[5]*x[7]*y[4]+z[7]*x[3]*y[6]-t49*z[4]+t145*z[7]-x[2]*y[7]*z[6]+t70*z[5]+t66*z[5]-z[7]*x[4]*y[6]+t18*z[4]+x[1]*y[4]*z[0];
    const double t160 = x[5]*y[4];
    const double t165 = z[1]*x[7];
    const double t178 = z[1]*x[3];
    const double t181 = t107*z[6]+t22*z[7]+t76*z[3]+t160*z[1]-x[4]*y[2]*z[6]+t70*z[4]+t165*y[5]+x[7]*y[2]*z[6]-t76*z[5]-t76*z[4]+t133*z[3]-t58*z[1]+y[5]*x[0]*z[4]+t114*z[2]-t3*z[7]+t20*z[2]+t178*y[7]+t129*y[2];
    const double t207 = t92*z[7]+t22*z[6]+z[3]*x[0]*y[2]-x[0]*y[3]*z[2]-z[3]*x[7]*y[2]-t165*y[3]-t9*y[0]+t58*z[7]+y[3]*x[6]*z[2]+t107*z[2]+t73*y[0]-x[3]*y[5]*z[7]+t3*z[0]-t56*y[6]-z[5]*x[0]*y[4]+t73*y[1]-t160*z[6]+t160*z[0];
    const double t228 = -t44*z[7]+z[5]*x[6]*y[4]-t52*z[4]-t145*z[4]+t68*z[4]+t92*z[5]-t92*z[0]+t11*z[3]+t44*z[0]+t178*y[5]-t46*y[5]-t178*y[0]-t145*z[0]-t20*z[5]-t37*z[3]-t160*z[7]+t145*z[3]+x[4]*y[6]*z[2];

    return (t34+t64+t95+t125+t156+t181+t207+t228)/12.;
  }



  template <>
  double
  cell_measure(const std::vector<Point<2> > &all_vertices,
               const unsigned int (&vertex_indices) [GeometryInfo<2>::vertices_per_cell])
  {
    /*
      Get the computation of the measure by this little Maple script. We
      use the blinear mapping of the unit quad to the real quad. However,
      every transformation mapping the unit faces to straight lines should
      do.

      Remember that the area of the quad is given by
      \int_K 1 dx dy  = \int_{\hat K} |det J| d(xi) d(eta)

      # x and y are arrays holding the x- and y-values of the four vertices
      # of this cell in real space.
      x := array(0..3);
      y := array(0..3);
      z := array(0..3);
      tphi[0] := (1-xi)*(1-eta):
      tphi[1] :=     xi*(1-eta):
      tphi[2] := (1-xi)*eta:
      tphi[3] :=     xi*eta:
      x_real := sum(x[s]*tphi[s], s=0..3):
      y_real := sum(y[s]*tphi[s], s=0..3):
      z_real := sum(z[s]*tphi[s], s=0..3):

      Jxi := <diff(x_real,xi)  | diff(y_real,xi) | diff(z_real,xi)>;
      Jeta := <diff(x_real,eta)| diff(y_real,eta)| diff(z_real,eta)>;
      with(VectorCalculus):
      J := CrossProduct(Jxi, Jeta);
      detJ := sqrt(J[1]^2 + J[2]^2 +J[3]^2);

      # measure := evalf (Int (Int (detJ, xi=0..1, method = _NCrule ) , eta=0..1, method = _NCrule  ) ):
      # readlib(C):

      # C(measure, optimized);

      additional optimizaton: divide by 2 only one time
    */

    const double x[4] = { all_vertices[vertex_indices[0]](0),
                          all_vertices[vertex_indices[1]](0),
                          all_vertices[vertex_indices[2]](0),
                          all_vertices[vertex_indices[3]](0)
                        };

    const double y[4] = { all_vertices[vertex_indices[0]](1),
                          all_vertices[vertex_indices[1]](1),
                          all_vertices[vertex_indices[2]](1),
                          all_vertices[vertex_indices[3]](1)
                        };

    return (-x[1]*y[0]+x[1]*y[3]+y[0]*x[2]+x[0]*y[1]-x[0]*y[2]-y[1]*x[3]-x[2]*y[3]+x[3]*y[2])/2;

  }




  template <int dim>
  double
  cell_measure(const std::vector<Point<dim> > &,
               const unsigned int ( &) [GeometryInfo<dim>::vertices_per_cell])
  {
    Assert(false, ExcNotImplemented());
    return 0.;
  }



  template <int dim, int spacedim>
  void
  delete_unused_vertices (std::vector<Point<spacedim> >    &vertices,
                          std::vector<CellData<dim> > &cells,
                          SubCellData                 &subcelldata)
  {
    Assert(subcelldata.check_consistency(dim),
           ExcMessage("Invalid SubCellData supplied according to ::check_consistency(). "
                      "This is caused by data containing objects for the wrong dimension."));

    // first check which vertices are actually used
    std::vector<bool> vertex_used (vertices.size(), false);
    for (unsigned int c=0; c<cells.size(); ++c)
      for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
        {
          Assert(cells[c].vertices[v] < vertices.size(),
                 ExcMessage("Invalid vertex index encountered! cells["
                            + Utilities::int_to_string(c)
                            + "].vertices["
                            + Utilities::int_to_string(v)
                            + "]="
                            + Utilities::int_to_string(cells[c].vertices[v])
                            + " is invalid, because only "
                            + Utilities::int_to_string(vertices.size())
                            + " vertices were supplied."));
          vertex_used[cells[c].vertices[v]] = true;
        }


    // then renumber the vertices that are actually used in the same order as
    // they were beforehand
    const unsigned int invalid_vertex = numbers::invalid_unsigned_int;
    std::vector<unsigned int> new_vertex_numbers (vertices.size(), invalid_vertex);
    unsigned int next_free_number = 0;
    for (unsigned int i=0; i<vertices.size(); ++i)
      if (vertex_used[i] == true)
        {
          new_vertex_numbers[i] = next_free_number;
          ++next_free_number;
        }

    // next replace old vertex numbers by the new ones
    for (unsigned int c=0; c<cells.size(); ++c)
      for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
        cells[c].vertices[v] = new_vertex_numbers[cells[c].vertices[v]];

    // same for boundary data
    for (unsigned int c=0; c<subcelldata.boundary_lines.size(); ++c)
      for (unsigned int v=0; v<GeometryInfo<1>::vertices_per_cell; ++v)
        {
          Assert(subcelldata.boundary_lines[c].vertices[v] < new_vertex_numbers.size(),
                 ExcMessage("Invalid vertex index in subcelldata.boundary_lines. "
                            "subcelldata.boundary_lines["
                            + Utilities::int_to_string(c)
                            + "].vertices["
                            + Utilities::int_to_string(v)
                            + "]="
                            + Utilities::int_to_string(subcelldata.boundary_lines[c].vertices[v])
                            + " is invalid, because only "
                            + Utilities::int_to_string(vertices.size())
                            + " vertices were supplied."));
          subcelldata.boundary_lines[c].vertices[v]
            = new_vertex_numbers[subcelldata.boundary_lines[c].vertices[v]];
        }

    for (unsigned int c=0; c<subcelldata.boundary_quads.size(); ++c)
      for (unsigned int v=0; v<GeometryInfo<2>::vertices_per_cell; ++v)
        {
          Assert(subcelldata.boundary_quads[c].vertices[v] < new_vertex_numbers.size(),
                 ExcMessage("Invalid vertex index in subcelldata.boundary_quads. "
                            "subcelldata.boundary_quads["
                            + Utilities::int_to_string(c)
                            + "].vertices["
                            + Utilities::int_to_string(v)
                            + "]="
                            + Utilities::int_to_string(subcelldata.boundary_quads[c].vertices[v])
                            + " is invalid, because only "
                            + Utilities::int_to_string(vertices.size())
                            + " vertices were supplied."));

          subcelldata.boundary_quads[c].vertices[v]
            = new_vertex_numbers[subcelldata.boundary_quads[c].vertices[v]];
        }

    // finally copy over the vertices which we really need to a new array and
    // replace the old one by the new one
    std::vector<Point<spacedim> > tmp;
    tmp.reserve (std::count(vertex_used.begin(), vertex_used.end(), true));
    for (unsigned int v=0; v<vertices.size(); ++v)
      if (vertex_used[v] == true)
        tmp.push_back (vertices[v]);
    swap (vertices, tmp);
  }



  template <int dim, int spacedim>
  void
  delete_duplicated_vertices (std::vector<Point<spacedim> >    &vertices,
                              std::vector<CellData<dim> > &cells,
                              SubCellData                 &subcelldata,
                              std::vector<unsigned int>   &considered_vertices,
                              double                       tol)
  {
    // create a vector of vertex
    // indices. initialize it to the identity,
    // later on change that if necessary.
    std::vector<unsigned int> new_vertex_numbers(vertices.size());
    for (unsigned int i=0; i<vertices.size(); ++i)
      new_vertex_numbers[i] = i;

    // if the considered_vertices vector is
    // empty, consider all vertices
    if (considered_vertices.size()==0)
      considered_vertices = new_vertex_numbers;

    Assert(considered_vertices.size() <= vertices.size(),
           ExcInternalError());


    // now loop over all vertices to be
    // considered and try to find an identical
    // one
    for (unsigned int i=0; i<considered_vertices.size(); ++i)
      {
        Assert(considered_vertices[i]<vertices.size(),
               ExcInternalError());
        if (new_vertex_numbers[considered_vertices[i]]!=considered_vertices[i])
          // this vertex has been identified with
          // another one already, skip it in the
          // test
          continue;
        // this vertex is not identified with
        // another one so far. search in the list
        // of remaining vertices. if a duplicate
        // vertex is found, set the new vertex
        // index for that vertex to this vertex'
        // index.
        for (unsigned int j=i+1; j<considered_vertices.size(); ++j)
          {
            bool equal=true;
            for (unsigned int d=0; d<spacedim; ++d)
              equal &= (fabs(vertices[considered_vertices[j]](d)-vertices[considered_vertices[i]](d))<tol);
            if (equal)
              {
                new_vertex_numbers[considered_vertices[j]]=considered_vertices[i];
                // we do not suppose, that there might be another duplicate
                // vertex, so break here
                break;
              }
          }
      }

    // now we got a renumbering list. simply
    // renumber all vertices (non-duplicate
    // vertices get renumbered to themselves, so
    // nothing bad happens). after that, the
    // duplicate vertices will be unused, so call
    // delete_unused_vertices() to do that part
    // of the job.
    for (unsigned int c=0; c<cells.size(); ++c)
      for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
        cells[c].vertices[v]=new_vertex_numbers[cells[c].vertices[v]];

    delete_unused_vertices(vertices, cells, subcelldata);
  }



// define some transformations in an anonymous namespace
  namespace
  {
    template <int spacedim>
    class Shift
    {
    public:
      Shift (const Tensor<1,spacedim> &shift)
        :
        shift(shift)
      {}
      Point<spacedim> operator() (const Point<spacedim> p) const
      {
        return p+shift;
      }
    private:
      const Tensor<1,spacedim> shift;
    };


    // the following class is only
    // needed in 2d, so avoid trouble
    // with compilers warning otherwise
    class Rotate2d
    {
    public:
      Rotate2d (const double angle)
        :
        angle(angle)
      {}
      Point<2> operator() (const Point<2> &p) const
      {
        return Point<2> (std::cos(angle)*p(0) - std::sin(angle) * p(1),
                         std::sin(angle)*p(0) + std::cos(angle) * p(1));
      }
    private:
      const double angle;
    };

    // Transformation to rotate around one of the cartesian axes.
    class Rotate3d
    {
    public:
      Rotate3d (const double angle,
                const unsigned int axis)
        :
        angle(angle),
        axis(axis)
      {}

      Point<3> operator() (const Point<3> &p) const
      {
        if (axis==0)
          return Point<3> (p(0),
                           std::cos(angle)*p(1) - std::sin(angle) * p(2),
                           std::sin(angle)*p(1) + std::cos(angle) * p(2));
        else if (axis==1)
          return Point<3> (std::cos(angle)*p(0) + std::sin(angle) * p(2),
                           p(1),
                           -std::sin(angle)*p(0) + std::cos(angle) * p(2));
        else
          return Point<3> (std::cos(angle)*p(0) - std::sin(angle) * p(1),
                           std::sin(angle)*p(0) + std::cos(angle) * p(1),
                           p(2));
      }
    private:
      const double angle;
      const unsigned int axis;
    };

    template <int spacedim>
    class Scale
    {
    public:
      Scale (const double factor)
        :
        factor(factor)
      {}
      Point<spacedim> operator() (const Point<spacedim> p) const
      {
        return p*factor;
      }
    private:
      const double factor;
    };
  }


  template <int dim, int spacedim>
  void
  shift (const Tensor<1,spacedim>   &shift_vector,
         Triangulation<dim, spacedim> &triangulation)
  {
    transform (Shift<spacedim>(shift_vector), triangulation);
  }



  void
  rotate (const double      angle,
          Triangulation<2> &triangulation)
  {
    transform (Rotate2d(angle), triangulation);
  }

  template<int dim>
  void
  rotate (const double      angle,
          const unsigned int axis,
          Triangulation<dim,3> &triangulation)
  {
    Assert(axis<3, ExcMessage("Invalid axis given!"));

    transform (Rotate3d(angle, axis), triangulation);
  }

  template <int dim, int spacedim>
  void
  scale (const double        scaling_factor,
         Triangulation<dim, spacedim> &triangulation)
  {
    Assert (scaling_factor>0, ExcScalingFactorNotPositive (scaling_factor));
    transform (Scale<spacedim>(scaling_factor), triangulation);
  }


  namespace
  {
    /**
     * Solve the Laplace equation for the @p laplace_transform function for one
     * of the @p dim space dimensions. Factorized into a function of its own
     * in order to allow parallel execution.
     */
    void laplace_solve (const SparseMatrix<double> &S,
                        const std::map<unsigned int,double> &m,
                        Vector<double> &u)
    {
      const unsigned int n_dofs=S.n();
      FilteredMatrix<Vector<double> > SF (S);
      PreconditionJacobi<SparseMatrix<double> > prec;
      prec.initialize(S, 1.2);
      FilteredMatrix<Vector<double> > PF (prec);

      SolverControl control (n_dofs, 1.e-10, false, false);
      GrowingVectorMemory<Vector<double> > mem;
      SolverCG<Vector<double> > solver (control, mem);

      Vector<double> f(n_dofs);

      SF.add_constraints(m);
      SF.apply_constraints (f, true);
      solver.solve(SF, u, f, PF);
    }
  }



  // Implementation for 1D only
  template <>
  void laplace_transform (const std::map<unsigned int,Point<1> > &,
                          Triangulation<1> &,
                          const Function<1> *)
  {
    Assert(false, ExcNotImplemented());
  }


  // Implementation for dimensions except 1
  template <int dim>
  void
  laplace_transform (const std::map<unsigned int,Point<dim> > &new_points,
                     Triangulation<dim> &triangulation,
                     const Function<dim> *coefficient)
  {
    // first provide everything that is needed for solving a Laplace
    // equation.
    FE_Q<dim> q1(1);

    DoFHandler<dim> dof_handler(triangulation);
    dof_handler.distribute_dofs(q1);

    DynamicSparsityPattern dsp (dof_handler.n_dofs (),
                                dof_handler.n_dofs ());
    DoFTools::make_sparsity_pattern (dof_handler, dsp);
    dsp.compress ();

    SparsityPattern sparsity_pattern;
    sparsity_pattern.copy_from (dsp);
    sparsity_pattern.compress ();

    SparseMatrix<double> S(sparsity_pattern);

    QGauss<dim> quadrature(4);

    MatrixCreator::create_laplace_matrix(StaticMappingQ1<dim>::mapping, dof_handler, quadrature, S, coefficient);

    // set up the boundary values for
    // the laplace problem
    std::vector<std::map<unsigned int,double> > m(dim);
    typename std::map<unsigned int,Point<dim> >::const_iterator map_end=new_points.end();

    // fill these maps using the data
    // given by new_points
    typename DoFHandler<dim>::cell_iterator cell=dof_handler.begin_active(),
                                            endc=dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
          {
            const typename DoFHandler<dim>::face_iterator face=cell->face(face_no);

            // loop over all vertices of the cell and see if it is listed in the map
            // given as first argument of the function
            for (unsigned int vertex_no=0;
                 vertex_no<GeometryInfo<dim>::vertices_per_face; ++vertex_no)
              {
                const unsigned int vertex_index=face->vertex_index(vertex_no);

                const typename std::map<unsigned int,Point<dim> >::const_iterator map_iter
                  = new_points.find(vertex_index);

                if (map_iter!=map_end)
                  for (unsigned int i=0; i<dim; ++i)
                    m[i].insert(std::pair<unsigned int,double> (
                                  face->vertex_dof_index(vertex_no, 0), map_iter->second(i)));
              }
          }
      }

    // solve the dim problems with
    // different right hand sides.
    Vector<double> us[dim];
    for (unsigned int i=0; i<dim; ++i)
      us[i].reinit (dof_handler.n_dofs());

    // solve linear systems in parallel
    Threads::TaskGroup<> tasks;
    for (unsigned int i=0; i<dim; ++i)
      tasks += Threads::new_task (&laplace_solve,
                                  S, m[i], us[i]);
    tasks.join_all ();

    // change the coordinates of the
    // points of the triangulation
    // according to the computed values
    for (cell=dof_handler.begin_active(); cell!=endc; ++cell)
      for (unsigned int vertex_no=0;
           vertex_no<GeometryInfo<dim>::vertices_per_cell; ++vertex_no)
        {
          Point<dim> &v=cell->vertex(vertex_no);
          const unsigned int dof_index=cell->vertex_dof_index(vertex_no, 0);
          for (unsigned int i=0; i<dim; ++i)
            v(i)=us[i](dof_index);
        }
  }

  template <int dim, int spacedim>
  std::map<unsigned int, Point<spacedim> >
  get_all_vertices_at_boundary (const Triangulation<dim, spacedim> &tria)
  {
    std::map<unsigned int, Point<spacedim> > vertex_map;
    typename Triangulation<dim,spacedim>::active_cell_iterator
    cell = tria.begin_active(),
    endc = tria.end();
    for (; cell!=endc; ++cell)
      {
        for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
          {
            const typename Triangulation<dim, spacedim>::face_iterator &face
              = cell->face(i);
            if (face->at_boundary())
              {
                for (unsigned j = 0; j < GeometryInfo<dim>::vertices_per_face; ++j)
                  {
                    const Point<spacedim> &vertex = face->vertex(j);
                    const unsigned int vertex_index = face->vertex_index(j);
                    vertex_map[vertex_index] = vertex;
                  }
              }
          }
      }
    return vertex_map;
  }

  /**
    * Distort a triangulation in
    * some random way.
    */
  template <int dim, int spacedim>
  void
  distort_random (const double                 factor,
                  Triangulation<dim,spacedim> &triangulation,
                  const bool                   keep_boundary)
  {
    // if spacedim>dim we need to make sure that we perturb
    // points but keep them on
    // the manifold. however, this isn't implemented right now
    Assert (spacedim == dim, ExcNotImplemented());


    // find the smallest length of the
    // lines adjacent to the
    // vertex. take the initial value
    // to be larger than anything that
    // might be found: the diameter of
    // the triangulation, here
    // estimated by adding up the
    // diameters of the coarse grid
    // cells.
    double almost_infinite_length = 0;
    for (typename Triangulation<dim,spacedim>::cell_iterator
         cell=triangulation.begin(0); cell!=triangulation.end(0); ++cell)
      almost_infinite_length += cell->diameter();

    std::vector<double> minimal_length (triangulation.n_vertices(),
                                        almost_infinite_length);

    // also note if a vertex is at the boundary
    std::vector<bool>   at_boundary (triangulation.n_vertices(), false);
    for (typename Triangulation<dim,spacedim>::active_cell_iterator
         cell=triangulation.begin_active(); cell!=triangulation.end(); ++cell)
      if (cell->is_locally_owned())
        {
          if (dim>1)
            {
              for (unsigned int i=0; i<GeometryInfo<dim>::lines_per_cell; ++i)
                {
                  const typename Triangulation<dim,spacedim>::line_iterator line
                    = cell->line(i);

                  if (keep_boundary && line->at_boundary())
                    {
                      at_boundary[line->vertex_index(0)] = true;
                      at_boundary[line->vertex_index(1)] = true;
                    }

                  minimal_length[line->vertex_index(0)]
                    = std::min(line->diameter(),
                               minimal_length[line->vertex_index(0)]);
                  minimal_length[line->vertex_index(1)]
                    = std::min(line->diameter(),
                               minimal_length[line->vertex_index(1)]);
                }
            }
          else //dim==1
            {
              if (keep_boundary)
                for (unsigned int vertex=0; vertex<2; ++vertex)
                  if (cell->at_boundary(vertex) == true)
                    at_boundary[cell->vertex_index(vertex)] = true;

              minimal_length[cell->vertex_index(0)]
                = std::min(cell->diameter(),
                           minimal_length[cell->vertex_index(0)]);
              minimal_length[cell->vertex_index(1)]
                = std::min(cell->diameter(),
                           minimal_length[cell->vertex_index(1)]);
            }
        }

    // create a random number generator for the interval [-1,1]. we use
    // this to make sure the distribution we get is repeatable, i.e.,
    // if you call the function twice on the same mesh, then you will
    // get the same mesh. this would not be the case if you used
    // the rand() function, which carries around some internal state
    boost::random::mt19937 rng;
    boost::random::uniform_real_distribution<> uniform_distribution(-1,1);

    // If the triangulation is distributed, we need to
    // exchange the moved vertices across mpi processes
    if (parallel::distributed::Triangulation< dim, spacedim > *distributed_triangulation
        = dynamic_cast<parallel::distributed::Triangulation<dim,spacedim>*> (&triangulation))
      {
        const std::vector<bool> locally_owned_vertices = get_locally_owned_vertices(triangulation);
        std::vector<bool>       vertex_moved (triangulation.n_vertices(), false);

        // Next move vertices on locally owned cells
        for (typename Triangulation<dim,spacedim>::active_cell_iterator
             cell=triangulation.begin_active(); cell!=triangulation.end(); ++cell)
          if (cell->is_locally_owned())
            {
              for (unsigned int vertex_no=0; vertex_no<GeometryInfo<dim>::vertices_per_cell;
                   ++vertex_no)
                {
                  const unsigned global_vertex_no = cell->vertex_index(vertex_no);

                  // ignore this vertex if we shall keep the boundary and
                  // this vertex *is* at the boundary, if it is already moved
                  // or if another process moves this vertex
                  if ((keep_boundary && at_boundary[global_vertex_no])
                      || vertex_moved[global_vertex_no]
                      || !locally_owned_vertices[global_vertex_no])
                    continue;

                  // first compute a random shift vector
                  Point<spacedim> shift_vector;
                  for (unsigned int d=0; d<spacedim; ++d)
                    shift_vector(d) = uniform_distribution(rng);

                  shift_vector *= factor * minimal_length[global_vertex_no] /
                                  std::sqrt(shift_vector.square());

                  // finally move the vertex
                  cell->vertex(vertex_no) += shift_vector;
                  vertex_moved[global_vertex_no] = true;
                }
            }

#ifdef DEAL_II_WITH_P4EST
        distributed_triangulation
        ->communicate_locally_moved_vertices(locally_owned_vertices);
#else
        (void)distributed_triangulation;
        Assert (false, ExcInternalError());
#endif
      }
    else
      // if this is a sequential triangulation, we could in principle
      // use the algorithm above, but we'll use an algorithm that we used
      // before the parallel::distributed::Triangulation was introduced
      // in order to preserve backward compatibility
      {
        // loop over all vertices and compute their new locations
        const unsigned int n_vertices = triangulation.n_vertices();
        std::vector<Point<spacedim> > new_vertex_locations (n_vertices);
        const std::vector<Point<spacedim> > &old_vertex_locations
          = triangulation.get_vertices();

        for (unsigned int vertex=0; vertex<n_vertices; ++vertex)
          {
            // ignore this vertex if we will keep the boundary and
            // this vertex *is* at the boundary
            if (keep_boundary && at_boundary[vertex])
              new_vertex_locations[vertex] = old_vertex_locations[vertex];
            else
              {
                // compute a random shift vector
                Point<spacedim> shift_vector;
                for (unsigned int d=0; d<spacedim; ++d)
                  shift_vector(d) = uniform_distribution(rng);

                shift_vector *= factor * minimal_length[vertex] /
                                std::sqrt(shift_vector.square());

                // record new vertex location
                new_vertex_locations[vertex] = old_vertex_locations[vertex] + shift_vector;
              }
          }

        // now do the actual move of the vertices
        for (typename Triangulation<dim,spacedim>::active_cell_iterator
             cell=triangulation.begin_active(); cell!=triangulation.end(); ++cell)
          for (unsigned int vertex_no=0;
               vertex_no<GeometryInfo<dim>::vertices_per_cell; ++vertex_no)
            cell->vertex(vertex_no) = new_vertex_locations[cell->vertex_index(vertex_no)];
      }

    // Correct hanging nodes if necessary
    if (dim>=2)
      {
        // We do the same as in GridTools::transform
        //
        // exclude hanging nodes at the boundaries of artificial cells:
        // these may belong to ghost cells for which we know the exact
        // location of vertices, whereas the artificial cell may or may
        // not be further refined, and so we cannot know whether
        // the location of the hanging node is correct or not
        typename Triangulation<dim,spacedim>::active_cell_iterator
        cell = triangulation.begin_active(),
        endc = triangulation.end();
        for (; cell!=endc; ++cell)
          if (!cell->is_artificial())
            for (unsigned int face=0;
                 face<GeometryInfo<dim>::faces_per_cell; ++face)
              if (cell->face(face)->has_children() &&
                  !cell->face(face)->at_boundary())
                {
                  // this face has hanging nodes
                  if (dim==2)
                    cell->face(face)->child(0)->vertex(1)
                      = (cell->face(face)->vertex(0) +
                         cell->face(face)->vertex(1)) / 2;
                  else if (dim==3)
                    {
                      cell->face(face)->child(0)->vertex(1)
                        = .5*(cell->face(face)->vertex(0)
                              +cell->face(face)->vertex(1));
                      cell->face(face)->child(0)->vertex(2)
                        = .5*(cell->face(face)->vertex(0)
                              +cell->face(face)->vertex(2));
                      cell->face(face)->child(1)->vertex(3)
                        = .5*(cell->face(face)->vertex(1)
                              +cell->face(face)->vertex(3));
                      cell->face(face)->child(2)->vertex(3)
                        = .5*(cell->face(face)->vertex(2)
                              +cell->face(face)->vertex(3));

                      // center of the face
                      cell->face(face)->child(0)->vertex(3)
                        = .25*(cell->face(face)->vertex(0)
                               +cell->face(face)->vertex(1)
                               +cell->face(face)->vertex(2)
                               +cell->face(face)->vertex(3));
                    }
                }
      }
  }



  template <int dim, template <int, int> class MeshType, int spacedim>
  unsigned int
  find_closest_vertex (const MeshType<dim,spacedim> &mesh,
                       const Point<spacedim>        &p)
  {
    // first get the underlying
    // triangulation from the
    // mesh and determine vertices
    // and used vertices
    const Triangulation<dim, spacedim> &tria = mesh.get_triangulation();

    const std::vector< Point<spacedim> > &vertices = tria.get_vertices();
    const std::vector< bool       > &used     = tria.get_used_vertices();

    // At the beginning, the first
    // used vertex is the closest one
    std::vector<bool>::const_iterator first =
      std::find(used.begin(), used.end(), true);

    // Assert that at least one vertex
    // is actually used
    Assert(first != used.end(), ExcInternalError());

    unsigned int best_vertex = std::distance(used.begin(), first);
    double       best_dist   = (p - vertices[best_vertex]).norm_square();

    // For all remaining vertices, test
    // whether they are any closer
    for (unsigned int j = best_vertex+1; j < vertices.size(); j++)
      if (used[j])
        {
          double dist = (p - vertices[j]).norm_square();
          if (dist < best_dist)
            {
              best_vertex = j;
              best_dist   = dist;
            }
        }

    return best_vertex;
  }


  template<int dim, template<int, int> class MeshType, int spacedim>
#ifndef _MSC_VER
  std::vector<typename MeshType<dim, spacedim>::active_cell_iterator>
#else
  std::vector<typename dealii::internal::ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim> >::type>
#endif
  find_cells_adjacent_to_vertex(const MeshType<dim,spacedim> &mesh,
                                const unsigned int            vertex)
  {
    // make sure that the given vertex is
    // an active vertex of the underlying
    // triangulation
    Assert(vertex < mesh.get_triangulation().n_vertices(),
           ExcIndexRange(0,mesh.get_triangulation().n_vertices(),vertex));
    Assert(mesh.get_triangulation().get_used_vertices()[vertex],
           ExcVertexNotUsed(vertex));

    // use a set instead of a vector
    // to ensure that cells are inserted only
    // once
    std::set<typename dealii::internal::ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim> >::type> adjacent_cells;

    typename dealii::internal::ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim> >::type
    cell = mesh.begin_active(),
    endc = mesh.end();

    // go through all active cells and look if the vertex is part of that cell
    //
    // in 1d, this is all we need to care about. in 2d/3d we also need to worry
    // that the vertex might be a hanging node on a face or edge of a cell; in
    // this case, we would want to add those cells as well on whose faces the
    // vertex is located but for which it is not a vertex itself.
    //
    // getting this right is a lot simpler in 2d than in 3d. in 2d, a hanging
    // node can only be in the middle of a face and we can query the neighboring
    // cell from the current cell. on the other hand, in 3d a hanging node
    // vertex can also be on an edge but there can be many other cells on
    // this edge and we can not access them from the cell we are currently
    // on.
    //
    // so, in the 3d case, if we run the algorithm as in 2d, we catch all
    // those cells for which the vertex we seek is on a *subface*, but we
    // miss the case of cells for which the vertex we seek is on a
    // sub-edge for which there is no corresponding sub-face (because the
    // immediate neighbor behind this face is not refined), see for example
    // the bits/find_cells_adjacent_to_vertex_6 testcase. thus, if we
    // haven't yet found the vertex for the current cell we also need to
    // look at the mid-points of edges
    //
    // as a final note, deciding whether a neighbor is actually coarser is
    // simple in the case of isotropic refinement (we just need to look at
    // the level of the current and the neighboring cell). however, this
    // isn't so simple if we have used anisotropic refinement since then
    // the level of a cell is not indicative of whether it is coarser or
    // not than the current cell. ultimately, we want to add all cells on
    // which the vertex is, independent of whether they are coarser or
    // finer and so in the 2d case below we simply add *any* *active* neighbor.
    // in the worst case, we add cells multiple times to the adjacent_cells
    // list, but std::set throws out those cells already entered
    for (; cell != endc; ++cell)
      {
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; v++)
          if (cell->vertex_index(v) == vertex)
            {
              // OK, we found a cell that contains
              // the given vertex. We add it
              // to the list.
              adjacent_cells.insert(cell);

              // as explained above, in 2+d we need to check whether
              // this vertex is on a face behind which there is a
              // (possibly) coarser neighbor. if this is the case,
              // then we need to also add this neighbor
              if (dim >= 2)
                for (unsigned int vface = 0; vface < dim; vface++)
                  {
                    const unsigned int face =
                      GeometryInfo<dim>::vertex_to_face[v][vface];

                    if (!cell->at_boundary(face)
                        &&
                        cell->neighbor(face)->active())
                      {
                        // there is a (possibly) coarser cell behind a
                        // face to which the vertex belongs. the
                        // vertex we are looking at is then either a
                        // vertex of that coarser neighbor, or it is a
                        // hanging node on one of the faces of that
                        // cell. in either case, it is adjacent to the
                        // vertex, so add it to the list as well (if
                        // the cell was already in the list then the
                        // std::set makes sure that we get it only
                        // once)
                        adjacent_cells.insert (cell->neighbor(face));
                      }
                  }

              // in any case, we have found a cell, so go to the next cell
              goto next_cell;
            }

        // in 3d also loop over the edges
        if (dim >= 3)
          {
            for (unsigned int e=0; e<GeometryInfo<dim>::lines_per_cell; ++e)
              if (cell->line(e)->has_children())
                // the only place where this vertex could have been
                // hiding is on the mid-edge point of the edge we
                // are looking at
                if (cell->line(e)->child(0)->vertex_index(1) == vertex)
                  {
                    adjacent_cells.insert(cell);

                    // jump out of this tangle of nested loops
                    goto next_cell;
                  }
          }

        // in more than 3d we would probably have to do the same as
        // above also for even lower-dimensional objects
        Assert (dim <= 3, ExcNotImplemented());

        // move on to the next cell if we have found the
        // vertex on the current one
next_cell:
        ;
      }

    // if this was an active vertex then there needs to have been
    // at least one cell to which it is adjacent!
    Assert (adjacent_cells.size() > 0, ExcInternalError());

    // return the result as a vector, rather than the set we built above
    return
      std::vector<typename dealii::internal::ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim> >::type>
      (adjacent_cells.begin(), adjacent_cells.end());
  }



  namespace
  {
    template <int dim, template<int, int> class MeshType, int spacedim>
    void find_active_cell_around_point_internal
    (const MeshType<dim,spacedim> &mesh,
#ifndef _MSC_VER
     std::set<typename MeshType<dim, spacedim>::active_cell_iterator> &searched_cells,
     std::set<typename MeshType<dim, spacedim>::active_cell_iterator> &adjacent_cells)
#else
     std::set<typename dealii::internal::ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim> >::type> &searched_cells,
     std::set<typename dealii::internal::ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim> >::type> &adjacent_cells)
#endif
    {
#ifndef _MSC_VER
      typedef typename MeshType<dim, spacedim>::active_cell_iterator cell_iterator;
#else
      typedef typename dealii::internal::ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim> >::type cell_iterator;
#endif

      // update the searched cells
      searched_cells.insert(adjacent_cells.begin(), adjacent_cells.end());
      // now we to collect all neighbors
      // of the cells in adjacent_cells we
      // have not yet searched.
      std::set<cell_iterator> adjacent_cells_new;

      typename std::set<cell_iterator>::const_iterator
      cell = adjacent_cells.begin(),
      endc = adjacent_cells.end();
      for (; cell != endc; ++cell)
        {
          std::vector<cell_iterator> active_neighbors;
          get_active_neighbors<MeshType<dim, spacedim> >(*cell, active_neighbors);
          for (unsigned int i=0; i<active_neighbors.size(); ++i)
            if (searched_cells.find(active_neighbors[i]) == searched_cells.end())
              adjacent_cells_new.insert(active_neighbors[i]);
        }
      adjacent_cells.clear();
      adjacent_cells.insert(adjacent_cells_new.begin(), adjacent_cells_new.end());
      if (adjacent_cells.size() == 0)
        {
          // we haven't found any other cell that would be a
          // neighbor of a previously found cell, but we know
          // that we haven't checked all cells yet. that means
          // that the domain is disconnected. in that case,
          // choose the first previously untouched cell we
          // can find
          cell_iterator it = mesh.begin_active();
          for ( ; it!=mesh.end(); ++it)
            if (searched_cells.find(it) == searched_cells.end())
              {
                adjacent_cells.insert(it);
                break;
              }
        }
    }
  }

  template <int dim, template<int, int> class MeshType, int spacedim>
#ifndef _MSC_VER
  typename MeshType<dim, spacedim>::active_cell_iterator
#else
  typename dealii::internal::ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim> >::type
#endif
  find_active_cell_around_point (const MeshType<dim,spacedim> &mesh,
                                 const Point<spacedim>        &p)
  {
    return
      find_active_cell_around_point<dim,MeshType,spacedim>
      (StaticMappingQ1<dim,spacedim>::mapping,
       mesh, p).first;
  }


  template <int dim, template <int, int> class MeshType, int spacedim>
#ifndef _MSC_VER
  std::pair<typename MeshType<dim, spacedim>::active_cell_iterator, Point<dim> >
#else
  std::pair<typename dealii::internal::ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim> >::type, Point<dim> >
#endif
  find_active_cell_around_point (const Mapping<dim,spacedim>  &mapping,
                                 const MeshType<dim,spacedim> &mesh,
                                 const Point<spacedim>        &p)
  {
    typedef typename dealii::internal::ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim> >::type active_cell_iterator;

    // The best distance is set to the
    // maximum allowable distance from
    // the unit cell; we assume a
    // max. deviation of 1e-10
    double best_distance = 1e-10;
    int    best_level = -1;
    std::pair<active_cell_iterator, Point<dim> > best_cell;

    // Find closest vertex and determine
    // all adjacent cells
    std::vector<active_cell_iterator> adjacent_cells_tmp
      = find_cells_adjacent_to_vertex(mesh,
                                      find_closest_vertex(mesh, p));

    // Make sure that we have found
    // at least one cell adjacent to vertex.
    Assert(adjacent_cells_tmp.size()>0, ExcInternalError());

    // Copy all the cells into a std::set
    std::set<active_cell_iterator> adjacent_cells (adjacent_cells_tmp.begin(),
                                                   adjacent_cells_tmp.end());
    std::set<active_cell_iterator> searched_cells;

    // Determine the maximal number of cells
    // in the grid.
    // As long as we have not found
    // the cell and have not searched
    // every cell in the triangulation,
    // we keep on looking.
    const unsigned int n_active_cells = mesh.get_triangulation().n_active_cells();
    bool found = false;
    unsigned int cells_searched = 0;
    while (!found && cells_searched < n_active_cells)
      {
        typename std::set<active_cell_iterator>::const_iterator
        cell = adjacent_cells.begin(),
        endc = adjacent_cells.end();
        for (; cell != endc; ++cell)
          {
            try
              {
                const Point<dim> p_cell = mapping.transform_real_to_unit_cell(*cell, p);

                // calculate the infinity norm of
                // the distance vector to the unit cell.
                const double dist = GeometryInfo<dim>::distance_to_unit_cell(p_cell);

                // We compare if the point is inside the
                // unit cell (or at least not too far
                // outside). If it is, it is also checked
                // that the cell has a more refined state
                if ((dist < best_distance)
                    ||
                    ((dist == best_distance)
                     &&
                     ((*cell)->level() > best_level)))
                  {
                    found         = true;
                    best_distance = dist;
                    best_level    = (*cell)->level();
                    best_cell     = std::make_pair(*cell, p_cell);
                  }
              }
            catch (typename MappingQGeneric<dim,spacedim>::ExcTransformationFailed &)
              {
                // ok, the transformation
                // failed presumably
                // because the point we
                // are looking for lies
                // outside the current
                // cell. this means that
                // the current cell can't
                // be the cell around the
                // point, so just ignore
                // this cell and move on
                // to the next
              }
          }

        // update the number of cells searched
        cells_searched += adjacent_cells.size();

        // if we have not found the cell in
        // question and have not yet searched every
        // cell, we expand our search to
        // all the not already searched neighbors of
        // the cells in adjacent_cells. This is
        // what find_active_cell_around_point_internal
        // is for.
        if (!found && cells_searched < n_active_cells)
          {
            find_active_cell_around_point_internal<dim,MeshType,spacedim>
            (mesh, searched_cells, adjacent_cells);
          }
      }

    AssertThrow (best_cell.first.state() == IteratorState::valid,
                 ExcPointNotFound<spacedim>(p));

    return best_cell;
  }



  template <int dim, int spacedim>
  std::pair<typename hp::DoFHandler<dim,spacedim>::active_cell_iterator, Point<dim> >
  find_active_cell_around_point (const hp::MappingCollection<dim,spacedim>   &mapping,
                                 const hp::DoFHandler<dim,spacedim> &mesh,
                                 const Point<spacedim>     &p)
  {
    Assert ((mapping.size() == 1) ||
            (mapping.size() == mesh.get_fe().size()),
            ExcMessage ("Mapping collection needs to have either size 1 "
                        "or size equal to the number of elements in "
                        "the FECollection."));

    typedef typename hp::DoFHandler<dim,spacedim>::active_cell_iterator cell_iterator;

    std::pair<cell_iterator, Point<dim> > best_cell;
    //If we have only one element in the MappingCollection,
    //we use find_active_cell_around_point using only one
    //mapping.
    if (mapping.size() == 1)
      best_cell = find_active_cell_around_point(mapping[0], mesh, p);
    else
      {


        // The best distance is set to the
        // maximum allowable distance from
        // the unit cell; we assume a
        // max. deviation of 1e-10
        double best_distance = 1e-10;
        int    best_level = -1;


        // Find closest vertex and determine
        // all adjacent cells
        unsigned int vertex = find_closest_vertex(mesh, p);

        std::vector<cell_iterator> adjacent_cells_tmp =
          find_cells_adjacent_to_vertex(mesh, vertex);

        // Make sure that we have found
        // at least one cell adjacent to vertex.
        Assert(adjacent_cells_tmp.size()>0, ExcInternalError());

        // Copy all the cells into a std::set
        std::set<cell_iterator> adjacent_cells(adjacent_cells_tmp.begin(), adjacent_cells_tmp.end());
        std::set<cell_iterator> searched_cells;

        // Determine the maximal number of cells
        // in the grid.
        // As long as we have not found
        // the cell and have not searched
        // every cell in the triangulation,
        // we keep on looking.
        const unsigned int n_cells = mesh.get_triangulation().n_cells();
        bool found = false;
        unsigned int cells_searched = 0;
        while (!found && cells_searched < n_cells)
          {
            typename std::set<cell_iterator>::const_iterator
            cell = adjacent_cells.begin(),
            endc = adjacent_cells.end();
            for (; cell != endc; ++cell)
              {
                try
                  {
                    const Point<dim> p_cell = mapping[(*cell)->active_fe_index()].transform_real_to_unit_cell(*cell, p);


                    // calculate the infinity norm of
                    // the distance vector to the unit cell.
                    const double dist = GeometryInfo<dim>::distance_to_unit_cell(p_cell);

                    // We compare if the point is inside the
                    // unit cell (or at least not too far
                    // outside). If it is, it is also checked
                    // that the cell has a more refined state
                    if (dist < best_distance ||
                        (dist == best_distance && (*cell)->level() > best_level))
                      {
                        found       = true;
                        best_distance = dist;
                        best_level    = (*cell)->level();
                        best_cell     = std::make_pair(*cell, p_cell);
                      }
                  }
                catch (typename MappingQGeneric<dim,spacedim>::ExcTransformationFailed &)
                  {
                    // ok, the transformation
                    // failed presumably
                    // because the point we
                    // are looking for lies
                    // outside the current
                    // cell. this means that
                    // the current cell can't
                    // be the cell around the
                    // point, so just ignore
                    // this cell and move on
                    // to the next
                  }
              }
            //udpate the number of cells searched
            cells_searched += adjacent_cells.size();
            // if we have not found the cell in
            // question and have not yet searched every
            // cell, we expand our search to
            // all the not already searched neighbors of
            // the cells in adjacent_cells.
            if (!found && cells_searched < n_cells)
              {
                find_active_cell_around_point_internal<dim,hp::DoFHandler,spacedim>
                (mesh, searched_cells, adjacent_cells);
              }

          }
      }

    AssertThrow (best_cell.first.state() == IteratorState::valid,
                 ExcPointNotFound<spacedim>(p));

    return best_cell;
  }


  namespace
  {

    template<class MeshType>
    bool
    contains_locally_owned_cells (const std::vector<typename MeshType::active_cell_iterator> &cells)
    {
      for (typename std::vector<typename MeshType::active_cell_iterator>::const_iterator
           it = cells.begin(); it != cells.end(); ++it)
        {
          if ((*it)->is_locally_owned())
            return true;
        }
      return false;
    }

    template<class MeshType>
    bool
    contains_artificial_cells (const std::vector<typename MeshType::active_cell_iterator> &cells)
    {
      for (typename std::vector<typename MeshType::active_cell_iterator>::const_iterator
           it = cells.begin(); it != cells.end(); ++it)
        {
          if ((*it)->is_artificial())
            return true;
        }
      return false;
    }

  }



  template <class MeshType>
  std::vector<typename MeshType::active_cell_iterator>
  compute_active_cell_halo_layer
  (const MeshType                                                                    &mesh,
   const std_cxx11::function<bool (const typename MeshType::active_cell_iterator &)> &predicate)
  {
    std::vector<typename MeshType::active_cell_iterator> active_halo_layer;
    std::vector<bool> locally_active_vertices_on_subdomain (mesh.get_triangulation().n_vertices(),
                                                            false);

    // Find the cells for which the predicate is true
    // These are the cells around which we wish to construct
    // the halo layer
    for (typename MeshType::active_cell_iterator
         cell = mesh.begin_active();
         cell != mesh.end(); ++cell)
      if (predicate(cell)) // True predicate --> Part of subdomain
        for (unsigned int v=0; v<GeometryInfo<MeshType::dimension>::vertices_per_cell; ++v)
          locally_active_vertices_on_subdomain[cell->vertex_index(v)] = true;

    // Find the cells that do not conform to the predicate
    // but share a vertex with the selected subdomain
    // These comprise the halo layer
    for (typename MeshType::active_cell_iterator
         cell = mesh.begin_active();
         cell != mesh.end(); ++cell)
      if (!predicate(cell)) // False predicate --> Potential halo cell
        for (unsigned int v=0; v<GeometryInfo<MeshType::dimension>::vertices_per_cell; ++v)
          if (locally_active_vertices_on_subdomain[cell->vertex_index(v)] == true)
            {
              active_halo_layer.push_back(cell);
              break;
            }

    return active_halo_layer;
  }



  template <class MeshType>
  std::vector<typename MeshType::active_cell_iterator>
  compute_ghost_cell_halo_layer (const MeshType &mesh)
  {
    std_cxx11::function<bool (const typename MeshType::active_cell_iterator &)> predicate
      = IteratorFilters::LocallyOwnedCell();

    const std::vector<typename MeshType::active_cell_iterator>
    active_halo_layer = compute_active_cell_halo_layer (mesh, predicate);

    // Check that we never return locally owned or artificial cells
    // What is left should only be the ghost cells
    Assert(contains_locally_owned_cells<MeshType>(active_halo_layer) == false,
           ExcMessage("Halo layer contains locally owned cells"));
    Assert(contains_artificial_cells<MeshType>(active_halo_layer) == false,
           ExcMessage("Halo layer contains artificial cells"));

    return active_halo_layer;
  }



  template <int dim, int spacedim>
  std::vector<std::set<typename Triangulation<dim,spacedim>::active_cell_iterator> >
  vertex_to_cell_map(const Triangulation<dim,spacedim> &triangulation)
  {
    std::vector<std::set<typename Triangulation<dim,spacedim>::active_cell_iterator> >
    vertex_to_cell_map(triangulation.n_vertices());
    typename Triangulation<dim,spacedim>::active_cell_iterator cell = triangulation.begin_active(),
                                                               endc = triangulation.end();
    for (; cell!=endc; ++cell)
      for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
        vertex_to_cell_map[cell->vertex_index(i)].insert(cell);

    // Take care of hanging nodes
    cell = triangulation.begin_active();
    for (; cell!=endc; ++cell)
      {
        for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
          {
            if ((cell->at_boundary(i)==false) && (cell->neighbor(i)->active()))
              {
                typename Triangulation<dim,spacedim>::active_cell_iterator adjacent_cell =
                  cell->neighbor(i);
                for (unsigned int j=0; j<GeometryInfo<dim>::vertices_per_face; ++j)
                  vertex_to_cell_map[cell->face(i)->vertex_index(j)].insert(adjacent_cell);
              }
          }

        // in 3d also loop over the edges
        if (dim==3)
          {
            for (unsigned int i=0; i<GeometryInfo<dim>::lines_per_cell; ++i)
              if (cell->line(i)->has_children())
                // the only place where this vertex could have been
                // hiding is on the mid-edge point of the edge we
                // are looking at
                vertex_to_cell_map[cell->line(i)->child(0)->vertex_index(1)].insert(cell);
          }
      }

    return vertex_to_cell_map;
  }



  template <int dim, int spacedim>
  std::map<unsigned int,types::global_vertex_index>
  compute_local_to_global_vertex_index_map(
    const parallel::distributed::Triangulation<dim,spacedim> &triangulation)
  {
    std::map<unsigned int,types::global_vertex_index> local_to_global_vertex_index;

#ifndef DEAL_II_WITH_MPI

    // without MPI, this function doesn't make sense because on cannot
    // use parallel::distributed::Triangulation in any meaninful
    // way
    (void)triangulation;
    Assert (false, ExcMessage ("This function does not make any sense "
                               "for parallel::distributed::Triangulation "
                               "objects if you do not have MPI enabled."));

#else

    typedef typename Triangulation<dim,spacedim>::active_cell_iterator active_cell_iterator;
    const std::vector<std::set<active_cell_iterator> > vertex_to_cell =
      vertex_to_cell_map(triangulation);

    // Create a local index for the locally "owned" vertices
    types::global_vertex_index next_index = 0;
    unsigned int max_cellid_size = 0;
    std::set<std::pair<types::subdomain_id,types::global_vertex_index> > vertices_added;
    std::map<types::subdomain_id,std::set<unsigned int> > vertices_to_recv;
    std::map<types::subdomain_id,std::vector<std_cxx11::tuple<types::global_vertex_index,
        types::global_vertex_index,std::string> > > vertices_to_send;
    active_cell_iterator cell = triangulation.begin_active(),
                         endc = triangulation.end();
    std::set<active_cell_iterator> missing_vert_cells;
    std::set<unsigned int> used_vertex_index;
    for (; cell!=endc; ++cell)
      {
        if (cell->is_locally_owned())
          {
            for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
              {
                types::subdomain_id lowest_subdomain_id = cell->subdomain_id();
                typename std::set<active_cell_iterator>::iterator
                adjacent_cell = vertex_to_cell[cell->vertex_index(i)].begin(),
                end_adj_cell = vertex_to_cell[cell->vertex_index(i)].end();
                for (; adjacent_cell!=end_adj_cell; ++adjacent_cell)
                  lowest_subdomain_id = std::min(lowest_subdomain_id,
                                                 (*adjacent_cell)->subdomain_id());

                // See if I "own" this vertex
                if (lowest_subdomain_id==cell->subdomain_id())
                  {
                    // Check that the vertex we are working on a vertex that has not be
                    // dealt with yet
                    if (used_vertex_index.find(cell->vertex_index(i))==used_vertex_index.end())
                      {
                        // Set the local index
                        local_to_global_vertex_index[cell->vertex_index(i)] = next_index++;

                        // Store the information that will be sent to the adjacent cells
                        // on other subdomains
                        adjacent_cell = vertex_to_cell[cell->vertex_index(i)].begin();
                        std::vector<types::subdomain_id> subdomain_ids;
                        for (; adjacent_cell!=end_adj_cell; ++adjacent_cell)
                          if ((*adjacent_cell)->subdomain_id()!=cell->subdomain_id())
                            {
                              std::pair<types::subdomain_id,types::global_vertex_index>
                              tmp((*adjacent_cell)->subdomain_id(), cell->vertex_index(i));
                              if (vertices_added.find(tmp)==vertices_added.end())
                                {
                                  vertices_to_send[(*adjacent_cell)->subdomain_id()].push_back(
                                    std_cxx11::tuple<types::global_vertex_index,types::global_vertex_index,
                                    std::string> (i,cell->vertex_index(i),
                                                  cell->id().to_string()));
                                  if (cell->id().to_string().size() > max_cellid_size)
                                    max_cellid_size = cell->id().to_string().size();
                                  vertices_added.insert(tmp);
                                }
                            }
                        used_vertex_index.insert(cell->vertex_index(i));
                      }
                  }
                else
                  {
                    // We don't own the vertex so we will receive its global index
                    vertices_to_recv[lowest_subdomain_id].insert(cell->vertex_index(i));
                    missing_vert_cells.insert(cell);
                  }
              }
          }

        // Some hanging nodes are vertices of ghost cells. They need to be
        // received.
        if (cell->is_ghost())
          {
            for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
              {
                if (cell->at_boundary(i)==false)
                  {
                    if (cell->neighbor(i)->active())
                      {
                        typename Triangulation<dim,spacedim>::active_cell_iterator adjacent_cell =
                          cell->neighbor(i);
                        if ((adjacent_cell->is_locally_owned()))
                          {
                            types::subdomain_id adj_subdomain_id = adjacent_cell->subdomain_id();
                            if (cell->subdomain_id()<adj_subdomain_id)
                              for (unsigned int j=0; j<GeometryInfo<dim>::vertices_per_face; ++j)
                                {
                                  vertices_to_recv[cell->subdomain_id()].insert(cell->face(i)->vertex_index(j));
                                  missing_vert_cells.insert(cell);
                                }
                          }
                      }
                  }
              }
          }
      }

    // Get the size of the largest CellID string
    max_cellid_size = Utilities::MPI::max(max_cellid_size, triangulation.get_communicator());

    // Make indices global by getting the number of vertices owned by each
    // processors and shifting the indices accordingly
    const unsigned int n_cpu = Utilities::MPI::n_mpi_processes(triangulation.get_communicator());
    std::vector<types::global_vertex_index> indices(n_cpu);
    MPI_Allgather(&next_index, 1, DEAL_II_DOF_INDEX_MPI_TYPE, &indices[0],
                  indices.size(), DEAL_II_DOF_INDEX_MPI_TYPE, triangulation.get_communicator());
    const types::global_vertex_index shift = std::accumulate(&indices[0],
                                                             &indices[0]+triangulation.locally_owned_subdomain(),0);

    std::map<unsigned int,types::global_vertex_index>::iterator
    global_index_it = local_to_global_vertex_index.begin(),
    global_index_end = local_to_global_vertex_index.end();
    for (; global_index_it!=global_index_end; ++global_index_it)
      global_index_it->second += shift;

    // In a first message, send the global ID of the vertices and the local
    // positions in the cells. In a second messages, send the cell ID as a
    // resize string. This is done in two messages so that types are not mixed

    // Send the first message
    std::vector<std::vector<types::global_vertex_index> > vertices_send_buffers(
      vertices_to_send.size());
    std::vector<MPI_Request> first_requests(vertices_to_send.size());
    typename std::map<types::subdomain_id,
             std::vector<std_cxx11::tuple<types::global_vertex_index,
             types::global_vertex_index,std::string> > >::iterator
             vert_to_send_it = vertices_to_send.begin(),
             vert_to_send_end = vertices_to_send.end();
    for (unsigned int i=0; vert_to_send_it!=vert_to_send_end;
         ++vert_to_send_it, ++i)
      {
        int destination = vert_to_send_it->first;
        const unsigned int n_vertices = vert_to_send_it->second.size();
        const int buffer_size = 2*n_vertices;
        vertices_send_buffers[i].resize(buffer_size);

        // fill the buffer
        for (unsigned int j=0; j<n_vertices; ++j)
          {
            vertices_send_buffers[i][2*j] = std_cxx11::get<0>(vert_to_send_it->second[j]);
            vertices_send_buffers[i][2*j+1] =
              local_to_global_vertex_index[std_cxx11::get<1>(vert_to_send_it->second[j])];
          }

        // Send the message
        MPI_Isend(&vertices_send_buffers[i][0],buffer_size,DEAL_II_VERTEX_INDEX_MPI_TYPE,
                  destination, 0, triangulation.get_communicator(), &first_requests[i]);
      }

    // Receive the first message
    std::vector<std::vector<types::global_vertex_index> > vertices_recv_buffers(
      vertices_to_recv.size());
    typename std::map<types::subdomain_id,std::set<unsigned int> >::iterator
    vert_to_recv_it = vertices_to_recv.begin(),
    vert_to_recv_end = vertices_to_recv.end();
    for (unsigned int i=0; vert_to_recv_it!=vert_to_recv_end; ++vert_to_recv_it, ++i)
      {
        int source = vert_to_recv_it->first;
        const unsigned int n_vertices = vert_to_recv_it->second.size();
        const int buffer_size = 2*n_vertices;
        vertices_recv_buffers[i].resize(buffer_size);

        // Receive the message
        MPI_Recv(&vertices_recv_buffers[i][0],buffer_size,DEAL_II_VERTEX_INDEX_MPI_TYPE,
                 source, 0, triangulation.get_communicator(), MPI_STATUS_IGNORE);
      }


    // Send second message
    std::vector<std::vector<char> > cellids_send_buffers(vertices_to_send.size());
    std::vector<MPI_Request> second_requests(vertices_to_send.size());
    vert_to_send_it = vertices_to_send.begin();
    for (unsigned int i=0; vert_to_send_it!=vert_to_send_end;
         ++vert_to_send_it, ++i)
      {
        int destination = vert_to_send_it->first;
        const unsigned int n_vertices = vert_to_send_it->second.size();
        const int buffer_size = max_cellid_size*n_vertices;
        cellids_send_buffers[i].resize(buffer_size);

        // fill the buffer
        unsigned int pos = 0;
        for (unsigned int j=0; j<n_vertices; ++j)
          {
            std::string cell_id = std_cxx11::get<2>(vert_to_send_it->second[j]);
            for (unsigned int k=0; k<max_cellid_size; ++k, ++pos)
              {
                if (k<cell_id.size())
                  cellids_send_buffers[i][pos] = cell_id[k];
                // if necessary fill up the reserved part of the buffer with an
                // invalid value
                else
                  cellids_send_buffers[i][pos] = '-';
              }
          }

        // Send the message
        MPI_Isend(&cellids_send_buffers[i][0], buffer_size, MPI_CHAR,
                  destination, 0, triangulation.get_communicator(), &second_requests[i]);
      }

    // Receive the second message
    std::vector<std::vector<char> > cellids_recv_buffers(vertices_to_recv.size());
    vert_to_recv_it = vertices_to_recv.begin();
    for (unsigned int i=0; vert_to_recv_it!=vert_to_recv_end; ++vert_to_recv_it, ++i)
      {
        int source = vert_to_recv_it->first;
        const unsigned int n_vertices = vert_to_recv_it->second.size();
        const int buffer_size = max_cellid_size*n_vertices;
        cellids_recv_buffers[i].resize(buffer_size);

        // Receive the message
        MPI_Recv(&cellids_recv_buffers[i][0],buffer_size, MPI_CHAR,
                 source, 0, triangulation.get_communicator(), MPI_STATUS_IGNORE);
      }


    // Match the data received with the required vertices
    vert_to_recv_it = vertices_to_recv.begin();
    for (unsigned int i=0; vert_to_recv_it!=vert_to_recv_end; ++i, ++vert_to_recv_it)
      {
        for (unsigned int j=0; j<vert_to_recv_it->second.size(); ++j)
          {
            const unsigned int local_pos_recv = vertices_recv_buffers[i][2*j];
            const types::global_vertex_index global_id_recv = vertices_recv_buffers[i][2*j+1];
            const std::string cellid_recv(&cellids_recv_buffers[i][max_cellid_size*j],
                                          &cellids_recv_buffers[i][max_cellid_size*(j+1)]);
            bool found = false;
            typename std::set<active_cell_iterator>::iterator
            cell_set_it = missing_vert_cells.begin(),
            end_cell_set = missing_vert_cells.end();
            for (; (found==false) && (cell_set_it!=end_cell_set); ++cell_set_it)
              {
                typename std::set<active_cell_iterator>::iterator
                candidate_cell = vertex_to_cell[(*cell_set_it)->vertex_index(i)].begin(),
                end_cell = vertex_to_cell[(*cell_set_it)->vertex_index(i)].end();
                for (; candidate_cell!=end_cell; ++candidate_cell)
                  {
                    std::string current_cellid = (*candidate_cell)->id().to_string();
                    current_cellid.resize(max_cellid_size,'-');
                    if (current_cellid.compare(cellid_recv)==0)
                      {
                        local_to_global_vertex_index[(*candidate_cell)->vertex_index(local_pos_recv)] =
                          global_id_recv;
                        found = true;

                        break;
                      }
                  }
              }
          }
      }
#endif

    return local_to_global_vertex_index;
  }



  template <int dim, int spacedim>
  void
  get_face_connectivity_of_cells (const Triangulation<dim,spacedim> &triangulation,
                                  DynamicSparsityPattern            &cell_connectivity)
  {
    cell_connectivity.reinit (triangulation.n_active_cells(),
                              triangulation.n_active_cells());

    // create a map pair<lvl,idx> -> SparsityPattern index
    // TODO: we are no longer using user_indices for this because we can get
    // pointer/index clashes when saving/restoring them. The following approach
    // works, but this map can get quite big. Not sure about more efficient solutions.
    std::map< std::pair<unsigned int,unsigned int>, unsigned int >
    indexmap;
    for (typename dealii::internal::ActiveCellIterator<dim, spacedim, Triangulation<dim, spacedim> >::type
         cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell)
      indexmap[std::pair<unsigned int,unsigned int>(cell->level(),cell->index())] = cell->active_cell_index();

    // next loop over all cells and their neighbors to build the sparsity
    // pattern. note that it's a bit hard to enter all the connections when a
    // neighbor has children since we would need to find out which of its
    // children is adjacent to the current cell. this problem can be omitted
    // if we only do something if the neighbor has no children -- in that case
    // it is either on the same or a coarser level than we are. in return, we
    // have to add entries in both directions for both cells
    for (typename dealii::internal::ActiveCellIterator<dim, spacedim, Triangulation<dim, spacedim> >::type
         cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell)
      {
        const unsigned int index = cell->active_cell_index();
        cell_connectivity.add (index, index);
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          if ((cell->at_boundary(f) == false)
              &&
              (cell->neighbor(f)->has_children() == false))
            {
              unsigned int other_index = indexmap.find(
                                           std::pair<unsigned int,unsigned int>(cell->neighbor(f)->level(),cell->neighbor(f)->index()))->second;
              cell_connectivity.add (index, other_index);
              cell_connectivity.add (other_index, index);
            }
      }
  }



  template <int dim, int spacedim>
  void
  get_face_connectivity_of_cells (const Triangulation<dim,spacedim> &triangulation,
                                  SparsityPattern                   &cell_connectivity)
  {
    DynamicSparsityPattern dsp;
    get_face_connectivity_of_cells(triangulation, dsp);
    cell_connectivity.copy_from(dsp);
  }



  template <int dim, int spacedim>
  void
  get_vertex_connectivity_of_cells (const Triangulation<dim,spacedim> &triangulation,
                                    DynamicSparsityPattern            &cell_connectivity)
  {
    std::vector<std::vector<unsigned int> > vertex_to_cell(triangulation.n_vertices());
    for (typename Triangulation<dim,spacedim>::active_cell_iterator cell=
           triangulation.begin_active(); cell != triangulation.end(); ++cell)
      {
        for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
          vertex_to_cell[cell->vertex_index(v)].push_back(cell->active_cell_index());
      }

    cell_connectivity.reinit (triangulation.n_active_cells(),
                              triangulation.n_active_cells());
    for (typename Triangulation<dim,spacedim>::active_cell_iterator cell=
           triangulation.begin_active(); cell != triangulation.end(); ++cell)
      {
        for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
          for (unsigned int n=0; n<vertex_to_cell[cell->vertex_index(v)].size(); ++n)
            cell_connectivity.add(cell->active_cell_index(), vertex_to_cell[cell->vertex_index(v)][n]);
      }
  }



  template <int dim, int spacedim>
  void
  partition_triangulation (const unsigned int           n_partitions,
                           Triangulation<dim,spacedim> &triangulation)
  {
    Assert ((dynamic_cast<parallel::distributed::Triangulation<dim,spacedim>*>
             (&triangulation)
             == 0),
            ExcMessage ("Objects of type parallel::distributed::Triangulation "
                        "are already partitioned implicitly and can not be "
                        "partitioned again explicitly."));
    Assert (n_partitions > 0, ExcInvalidNumberOfPartitions(n_partitions));

    // check for an easy return
    if (n_partitions == 1)
      {
        for (typename dealii::internal::ActiveCellIterator<dim, spacedim, Triangulation<dim, spacedim> >::type
             cell = triangulation.begin_active();
             cell != triangulation.end(); ++cell)
          cell->set_subdomain_id (0);
        return;
      }

    // we decompose the domain by first
    // generating the connection graph of all
    // cells with their neighbors, and then
    // passing this graph off to METIS.
    // finally defer to the other function for
    // partitioning and assigning subdomain ids
    SparsityPattern cell_connectivity;
    get_face_connectivity_of_cells (triangulation, cell_connectivity);

    partition_triangulation (n_partitions,
                             cell_connectivity,
                             triangulation);
  }



  template <int dim, int spacedim>
  void
  partition_triangulation (const unsigned int           n_partitions,
                           const SparsityPattern        &cell_connection_graph,
                           Triangulation<dim,spacedim>  &triangulation)
  {
    Assert ((dynamic_cast<parallel::distributed::Triangulation<dim,spacedim>*>
             (&triangulation)
             == 0),
            ExcMessage ("Objects of type parallel::distributed::Triangulation "
                        "are already partitioned implicitly and can not be "
                        "partitioned again explicitly."));
    Assert (n_partitions > 0, ExcInvalidNumberOfPartitions(n_partitions));
    Assert (cell_connection_graph.n_rows() == triangulation.n_active_cells(),
            ExcMessage ("Connectivity graph has wrong size"));
    Assert (cell_connection_graph.n_cols() == triangulation.n_active_cells(),
            ExcMessage ("Connectivity graph has wrong size"));

    // check for an easy return
    if (n_partitions == 1)
      {
        for (typename dealii::internal::ActiveCellIterator<dim, spacedim, Triangulation<dim, spacedim> >::type
             cell = triangulation.begin_active();
             cell != triangulation.end(); ++cell)
          cell->set_subdomain_id (0);
        return;
      }

    // partition this connection graph and get
    // back a vector of indices, one per degree
    // of freedom (which is associated with a
    // cell)
    std::vector<unsigned int> partition_indices (triangulation.n_active_cells());
    SparsityTools::partition (cell_connection_graph, n_partitions,  partition_indices);

    // finally loop over all cells and set the
    // subdomain ids
    for (typename dealii::internal::ActiveCellIterator<dim, spacedim, Triangulation<dim, spacedim> >::type
         cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell)
      cell->set_subdomain_id (partition_indices[cell->active_cell_index()]);
  }



  template <int dim, int spacedim>
  void
  get_subdomain_association (const Triangulation<dim, spacedim>  &triangulation,
                             std::vector<types::subdomain_id> &subdomain)
  {
    Assert (subdomain.size() == triangulation.n_active_cells(),
            ExcDimensionMismatch (subdomain.size(),
                                  triangulation.n_active_cells()));
    for (typename Triangulation<dim, spacedim>::active_cell_iterator
         cell = triangulation.begin_active(); cell!=triangulation.end(); ++cell)
      subdomain[cell->active_cell_index()] = cell->subdomain_id();
  }



  template <int dim, int spacedim>
  unsigned int
  count_cells_with_subdomain_association (const Triangulation<dim, spacedim> &triangulation,
                                          const types::subdomain_id       subdomain)
  {
    unsigned int count = 0;
    for (typename Triangulation<dim, spacedim>::active_cell_iterator
         cell = triangulation.begin_active();
         cell!=triangulation.end(); ++cell)
      if (cell->subdomain_id() == subdomain)
        ++count;

    return count;
  }



  template <int dim, int spacedim>
  std::vector<bool>
  get_locally_owned_vertices (const Triangulation<dim,spacedim> &triangulation)
  {
    // start with all vertices
    std::vector<bool> locally_owned_vertices = triangulation.get_used_vertices();

    // if the triangulation is distributed, eliminate those that
    // are owned by other processors -- either because the vertex is
    // on an artificial cell, or because it is on a ghost cell with
    // a smaller subdomain
    if (const parallel::distributed::Triangulation<dim,spacedim> *tr
        = dynamic_cast<const parallel::distributed::Triangulation<dim,spacedim> *>
        (&triangulation))
      for (typename dealii::internal::ActiveCellIterator<dim, spacedim, Triangulation<dim, spacedim> >::type
           cell = triangulation.begin_active();
           cell != triangulation.end(); ++cell)
        if (cell->is_artificial()
            ||
            (cell->is_ghost() &&
             (cell->subdomain_id() < tr->locally_owned_subdomain())))
          for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
            locally_owned_vertices[cell->vertex_index(v)] = false;

    return locally_owned_vertices;
  }



  template <typename MeshType>
  std::list<std::pair<typename MeshType::cell_iterator,
      typename MeshType::cell_iterator> >
      get_finest_common_cells (const MeshType &mesh_1,
                               const MeshType &mesh_2)
  {
    Assert (have_same_coarse_mesh (mesh_1, mesh_2),
            ExcMessage ("The two meshes must be represent triangulations that "
                        "have the same coarse meshes"));

    // the algorithm goes as follows:
    // first, we fill a list with pairs
    // of iterators common to the two
    // meshes on the coarsest
    // level. then we traverse the
    // list; each time, we find a pair
    // of iterators for which both
    // correspond to non-active cells,
    // we delete this item and push the
    // pairs of iterators to their
    // children to the back. if these
    // again both correspond to
    // non-active cells, we will get to
    // the later on for further
    // consideration
    typedef
    std::list<std::pair<typename MeshType::cell_iterator,
        typename MeshType::cell_iterator> >
        CellList;

    CellList cell_list;

    // first push the coarse level cells
    typename MeshType::cell_iterator
    cell_1 = mesh_1.begin(0),
    cell_2 = mesh_2.begin(0);
    for (; cell_1 != mesh_1.end(0); ++cell_1, ++cell_2)
      cell_list.push_back (std::make_pair (cell_1, cell_2));

    // then traverse list as described
    // above
    typename CellList::iterator cell_pair = cell_list.begin();
    while (cell_pair != cell_list.end())
      {
        // if both cells in this pair
        // have children, then erase
        // this element and push their
        // children instead
        if (cell_pair->first->has_children()
            &&
            cell_pair->second->has_children())
          {
            Assert(cell_pair->first->refinement_case()==
                   cell_pair->second->refinement_case(), ExcNotImplemented());
            for (unsigned int c=0; c<cell_pair->first->n_children(); ++c)
              cell_list.push_back (std::make_pair (cell_pair->first->child(c),
                                                   cell_pair->second->child(c)));

            // erasing an iterator
            // keeps other iterators
            // valid, so already
            // advance the present
            // iterator by one and then
            // delete the element we've
            // visited before
            const typename CellList::iterator previous_cell_pair = cell_pair;
            ++cell_pair;

            cell_list.erase (previous_cell_pair);
          }
        else
          // both cells are active, do
          // nothing
          ++cell_pair;
      }

    // just to make sure everything is ok,
    // validate that all pairs have at least one
    // active iterator or have different
    // refinement_cases
    for (cell_pair = cell_list.begin(); cell_pair != cell_list.end(); ++cell_pair)
      Assert (cell_pair->first->active()
              ||
              cell_pair->second->active()
              ||
              (cell_pair->first->refinement_case()
               != cell_pair->second->refinement_case()),
              ExcInternalError());

    return cell_list;
  }

  template <int dim, int spacedim>
  bool
  have_same_coarse_mesh (const Triangulation<dim, spacedim> &mesh_1,
                         const Triangulation<dim, spacedim> &mesh_2)
  {
    // make sure the two meshes have
    // the same number of coarse cells
    if (mesh_1.n_cells (0) != mesh_2.n_cells (0))
      return false;

    // if so, also make sure they have
    // the same vertices on the cells
    // of the coarse mesh
    typename Triangulation<dim, spacedim>::cell_iterator
    cell_1 = mesh_1.begin(0),
    cell_2 = mesh_2.begin(0),
    endc   = mesh_1.end(0);
    for (; cell_1!=endc; ++cell_1, ++cell_2)
      for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
        if (cell_1->vertex(v) != cell_2->vertex(v))
          return false;

    // if we've gotten through all
    // this, then the meshes really
    // seem to have a common coarse
    // mesh
    return true;
  }



  template <typename MeshType>
  bool
  have_same_coarse_mesh (const MeshType &mesh_1,
                         const MeshType &mesh_2)
  {
    return have_same_coarse_mesh (mesh_1.get_triangulation(),
                                  mesh_2.get_triangulation());
  }



  template <int dim, int spacedim>
  double
  minimal_cell_diameter (const Triangulation<dim, spacedim> &triangulation)
  {
    double min_diameter = triangulation.begin_active()->diameter();
    for (typename Triangulation<dim, spacedim>::active_cell_iterator
         cell = triangulation.begin_active(); cell != triangulation.end();
         ++cell)
      min_diameter = std::min (min_diameter,
                               cell->diameter());
    return min_diameter;
  }



  template <int dim, int spacedim>
  double
  maximal_cell_diameter (const Triangulation<dim, spacedim> &triangulation)
  {
    double max_diameter = triangulation.begin_active()->diameter();
    for (typename Triangulation<dim, spacedim>::active_cell_iterator
         cell = triangulation.begin_active(); cell != triangulation.end();
         ++cell)
      max_diameter = std::max (max_diameter,
                               cell->diameter());
    return max_diameter;
  }



  namespace internal
  {
    namespace FixUpDistortedChildCells
    {
      // compute the mean square
      // deviation of the alternating
      // forms of the children of the
      // given object from that of
      // the object itself. for
      // objects with
      // structdim==spacedim, the
      // alternating form is the
      // determinant of the jacobian,
      // whereas for faces with
      // structdim==spacedim-1, the
      // alternating form is the
      // (signed and scaled) normal
      // vector
      //
      // this average square
      // deviation is computed for an
      // object where the center node
      // has been replaced by the
      // second argument to this
      // function
      template <typename Iterator, int spacedim>
      double
      objective_function (const Iterator &object,
                          const Point<spacedim> &object_mid_point)
      {
        const unsigned int structdim = Iterator::AccessorType::structure_dimension;
        Assert (spacedim == Iterator::AccessorType::dimension,
                ExcInternalError());

        // everything below is wrong
        // if not for the following
        // condition
        Assert (object->refinement_case() == RefinementCase<structdim>::isotropic_refinement,
                ExcNotImplemented());
        // first calculate the
        // average alternating form
        // for the parent cell/face
        Point<spacedim> parent_vertices
        [GeometryInfo<structdim>::vertices_per_cell];
        Tensor<spacedim-structdim,spacedim> parent_alternating_forms
        [GeometryInfo<structdim>::vertices_per_cell];

        for (unsigned int i=0; i<GeometryInfo<structdim>::vertices_per_cell; ++i)
          parent_vertices[i] = object->vertex(i);

        GeometryInfo<structdim>::alternating_form_at_vertices (parent_vertices,
                                                               parent_alternating_forms);

        const Tensor<spacedim-structdim,spacedim>
        average_parent_alternating_form
          = std::accumulate (&parent_alternating_forms[0],
                             &parent_alternating_forms[GeometryInfo<structdim>::vertices_per_cell],
                             Tensor<spacedim-structdim,spacedim>());

        // now do the same
        // computation for the
        // children where we use the
        // given location for the
        // object mid point instead of
        // the one the triangulation
        // currently reports
        Point<spacedim> child_vertices
        [GeometryInfo<structdim>::max_children_per_cell]
        [GeometryInfo<structdim>::vertices_per_cell];
        Tensor<spacedim-structdim,spacedim> child_alternating_forms
        [GeometryInfo<structdim>::max_children_per_cell]
        [GeometryInfo<structdim>::vertices_per_cell];

        for (unsigned int c=0; c<object->n_children(); ++c)
          for (unsigned int i=0; i<GeometryInfo<structdim>::vertices_per_cell; ++i)
            child_vertices[c][i] = object->child(c)->vertex(i);

        // replace mid-object
        // vertex. note that for
        // child i, the mid-object
        // vertex happens to have the
        // number
        // max_children_per_cell-i
        for (unsigned int c=0; c<object->n_children(); ++c)
          child_vertices[c][GeometryInfo<structdim>::max_children_per_cell-c-1]
            = object_mid_point;

        for (unsigned int c=0; c<object->n_children(); ++c)
          GeometryInfo<structdim>::alternating_form_at_vertices (child_vertices[c],
                                                                 child_alternating_forms[c]);

        // on a uniformly refined
        // hypercube object, the child
        // alternating forms should
        // all be smaller by a factor
        // of 2^structdim than the
        // ones of the parent. as a
        // consequence, we'll use the
        // squared deviation from
        // this ideal value as an
        // objective function
        double objective = 0;
        for (unsigned int c=0; c<object->n_children(); ++c)
          for (unsigned int i=0; i<GeometryInfo<structdim>::vertices_per_cell; ++i)
            objective += (child_alternating_forms[c][i] -
                          average_parent_alternating_form/std::pow(2.,1.*structdim))
                         .norm_square();

        return objective;
      }


      /**
       * Return the location of the midpoint
       * of the 'f'th face (vertex) of this 1d
       * object.
       */
      template <typename Iterator>
      Point<Iterator::AccessorType::space_dimension>
      get_face_midpoint (const Iterator &object,
                         const unsigned int f,
                         dealii::internal::int2type<1>)
      {
        return object->vertex(f);
      }



      /**
       * Return the location of the midpoint
       * of the 'f'th face (line) of this 2d
       * object.
       */
      template <typename Iterator>
      Point<Iterator::AccessorType::space_dimension>
      get_face_midpoint (const Iterator &object,
                         const unsigned int f,
                         dealii::internal::int2type<2>)
      {
        return object->line(f)->center();
      }



      /**
       * Return the location of the midpoint
       * of the 'f'th face (quad) of this 3d
       * object.
       */
      template <typename Iterator>
      Point<Iterator::AccessorType::space_dimension>
      get_face_midpoint (const Iterator &object,
                         const unsigned int f,
                         dealii::internal::int2type<3>)
      {
        return object->face(f)->center();
      }




      /**
       * Compute the minimal diameter of an
       * object by looking for the minimal
       * distance between the mid-points of
       * its faces. This minimal diameter is
       * used to determine the step length
       * for our grid cell improvement
       * algorithm, and it should be small
       * enough that the point moves around
       * within the cell even if it is highly
       * elongated -- thus, the diameter of
       * the object is not a good measure,
       * while the minimal diameter is. Note
       * that the algorithm below works for
       * both cells that are long rectangles
       * with parallel sides where the
       * nearest distance is between opposite
       * edges as well as highly slanted
       * parallelograms where the shortest
       * distance is between neighboring
       * edges.
       */
      template <typename Iterator>
      double
      minimal_diameter (const Iterator &object)
      {
        const unsigned int
        structdim = Iterator::AccessorType::structure_dimension;

        double diameter = object->diameter();
        for (unsigned int f=0;
             f<GeometryInfo<structdim>::faces_per_cell;
             ++f)
          for (unsigned int e=f+1;
               e<GeometryInfo<structdim>::faces_per_cell;
               ++e)
            diameter = std::min (diameter,
                                 get_face_midpoint
                                 (object, f,
                                  dealii::internal::int2type<structdim>())
                                 .distance (get_face_midpoint
                                            (object,
                                             e,
                                             dealii::internal::int2type<structdim>())));

        return diameter;
      }



      /**
       * Try to fix up a single cell. Return
       * whether we succeeded with this.
       *
       * The second argument indicates
       * whether we need to respect the
       * manifold/boundary on which this
       * object lies when moving around its
       * mid-point.
       */
      template <typename Iterator>
      bool
      fix_up_object (const Iterator &object,
                     const bool respect_manifold)
      {
        const Boundary<Iterator::AccessorType::dimension,
              Iterator::AccessorType::space_dimension>
              *manifold = (respect_manifold ?
                           &object->get_boundary() :
                           0);

        const unsigned int structdim = Iterator::AccessorType::structure_dimension;
        const unsigned int spacedim  = Iterator::AccessorType::space_dimension;

        // right now we can only deal
        // with cells that have been
        // refined isotropically
        // because that is the only
        // case where we have a cell
        // mid-point that can be moved
        // around without having to
        // consider boundary
        // information
        Assert (object->has_children(), ExcInternalError());
        Assert (object->refinement_case() == RefinementCase<structdim>::isotropic_refinement,
                ExcNotImplemented());

        // get the current location of
        // the object mid-vertex:
        Point<spacedim> object_mid_point
          = object->child(0)->vertex (GeometryInfo<structdim>::max_children_per_cell-1);

        // now do a few steepest descent
        // steps to reduce the objective
        // function. compute the diameter in
        // the helper function above
        unsigned int iteration = 0;
        const double diameter = minimal_diameter (object);

        // current value of objective
        // function and initial delta
        double current_value = objective_function (object, object_mid_point);
        double initial_delta = 0;

        do
          {
            // choose a step length
            // that is initially 1/4
            // of the child objects'
            // diameter, and a sequence
            // whose sum does not
            // converge (to avoid
            // premature termination of
            // the iteration)
            const double step_length = diameter / 4 / (iteration + 1);

            // compute the objective
            // function's derivative using a
            // two-sided difference formula
            // with eps=step_length/10
            Tensor<1,spacedim> gradient;
            for (unsigned int d=0; d<spacedim; ++d)
              {
                const double eps = step_length/10;

                Tensor<1,spacedim> h;
                h[d] = eps/2;

                if (respect_manifold == false)
                  gradient[d]
                    = ((objective_function (object, object_mid_point + h)
                        -
                        objective_function (object, object_mid_point - h))
                       /
                       eps);
                else
                  gradient[d]
                    = ((objective_function (object,
                                            manifold->project_to_surface(object,
                                                                         object_mid_point + h))
                        -
                        objective_function (object,
                                            manifold->project_to_surface(object,
                                                                         object_mid_point - h)))
                       /
                       eps);
              }

            // sometimes, the
            // (unprojected) gradient
            // is perpendicular to
            // the manifold, but we
            // can't go there if
            // respect_manifold==true. in
            // that case, gradient=0,
            // and we simply need to
            // quite the loop here
            if (gradient.norm() == 0)
              break;

            // so we need to go in
            // direction -gradient. the
            // optimal value of the
            // objective function is
            // zero, so assuming that
            // the model is quadratic
            // we would have to go
            // -2*val/||gradient|| in
            // this direction, make
            // sure we go at most
            // step_length into this
            // direction
            object_mid_point -= std::min(2 * current_value / (gradient*gradient),
                                         step_length / gradient.norm()) *
                                gradient;

            if (respect_manifold == true)
              object_mid_point = manifold->project_to_surface(object,
                                                              object_mid_point);

            // compute current value of the
            // objective function
            const double previous_value = current_value;
            current_value = objective_function (object, object_mid_point);

            if (iteration == 0)
              initial_delta = (previous_value - current_value);

            // stop if we aren't moving much
            // any more
            if ((iteration >= 1) &&
                ((previous_value - current_value < 0)
                 ||
                 (std::fabs (previous_value - current_value)
                  <
                  0.001 * initial_delta)))
              break;

            ++iteration;
          }
        while (iteration < 20);

        // verify that the new
        // location is indeed better
        // than the one before. check
        // this by comparing whether
        // the minimum value of the
        // products of parent and
        // child alternating forms is
        // positive. for cells this
        // means that the
        // determinants have the same
        // sign, for faces that the
        // face normals of parent and
        // children point in the same
        // general direction
        double old_min_product, new_min_product;

        Point<spacedim> parent_vertices
        [GeometryInfo<structdim>::vertices_per_cell];
        for (unsigned int i=0; i<GeometryInfo<structdim>::vertices_per_cell; ++i)
          parent_vertices[i] = object->vertex(i);

        Tensor<spacedim-structdim,spacedim> parent_alternating_forms
        [GeometryInfo<structdim>::vertices_per_cell];
        GeometryInfo<structdim>::alternating_form_at_vertices (parent_vertices,
                                                               parent_alternating_forms);

        Point<spacedim> child_vertices
        [GeometryInfo<structdim>::max_children_per_cell]
        [GeometryInfo<structdim>::vertices_per_cell];

        for (unsigned int c=0; c<object->n_children(); ++c)
          for (unsigned int i=0; i<GeometryInfo<structdim>::vertices_per_cell; ++i)
            child_vertices[c][i] = object->child(c)->vertex(i);

        Tensor<spacedim-structdim,spacedim> child_alternating_forms
        [GeometryInfo<structdim>::max_children_per_cell]
        [GeometryInfo<structdim>::vertices_per_cell];

        for (unsigned int c=0; c<object->n_children(); ++c)
          GeometryInfo<structdim>::alternating_form_at_vertices (child_vertices[c],
                                                                 child_alternating_forms[c]);

        old_min_product = child_alternating_forms[0][0] * parent_alternating_forms[0];
        for (unsigned int c=0; c<object->n_children(); ++c)
          for (unsigned int i=0; i<GeometryInfo<structdim>::vertices_per_cell; ++i)
            for (unsigned int j=0; j<GeometryInfo<structdim>::vertices_per_cell; ++j)
              old_min_product =
                std::min<double> (old_min_product,
                                  child_alternating_forms[c][i] *
                                  parent_alternating_forms[j]);

        // for the new minimum value,
        // replace mid-object
        // vertex. note that for child
        // i, the mid-object vertex
        // happens to have the number
        // max_children_per_cell-i
        for (unsigned int c=0; c<object->n_children(); ++c)
          child_vertices[c][GeometryInfo<structdim>::max_children_per_cell-c-1]
            = object_mid_point;

        for (unsigned int c=0; c<object->n_children(); ++c)
          GeometryInfo<structdim>::alternating_form_at_vertices (child_vertices[c],
                                                                 child_alternating_forms[c]);

        new_min_product = child_alternating_forms[0][0] * parent_alternating_forms[0];
        for (unsigned int c=0; c<object->n_children(); ++c)
          for (unsigned int i=0; i<GeometryInfo<structdim>::vertices_per_cell; ++i)
            for (unsigned int j=0; j<GeometryInfo<structdim>::vertices_per_cell; ++j)
              new_min_product =
                std::min<double> (new_min_product,
                                  child_alternating_forms[c][i] *
                                  parent_alternating_forms[j]);

        // if new minimum value is
        // better than before, then set the
        // new mid point. otherwise
        // return this object as one of
        // those that can't apparently
        // be fixed
        if (new_min_product >= old_min_product)
          object->child(0)->vertex (GeometryInfo<structdim>::max_children_per_cell-1)
            = object_mid_point;

        // return whether after this
        // operation we have an object that
        // is well oriented
        return (std::max (new_min_product, old_min_product) > 0);
      }



      void fix_up_faces (const dealii::Triangulation<1,1>::cell_iterator &,
                         dealii::internal::int2type<1>,
                         dealii::internal::int2type<1>)
      {
        // nothing to do for the faces of
        // cells in 1d
      }



      // possibly fix up the faces of
      // a cell by moving around its
      // mid-points
      template <int structdim, int spacedim>
      void fix_up_faces (const typename dealii::Triangulation<structdim,spacedim>::cell_iterator &cell,
                         dealii::internal::int2type<structdim>,
                         dealii::internal::int2type<spacedim>)
      {
        // see if we first can fix up
        // some of the faces of this
        // object. we can mess with
        // faces if and only if it is
        // not at the boundary (since
        // otherwise the location of
        // the face mid-point has been
        // determined by the boundary
        // object) and if the
        // neighboring cell is not even
        // more refined than we are
        // (since in that case the
        // sub-faces have themselves
        // children that we can't move
        // around any more). however,
        // the latter case shouldn't
        // happen anyway: if the
        // current face is distorted
        // but the neighbor is even
        // more refined, then the face
        // had been deformed before
        // already, and had been
        // ignored at the time; we
        // should then also be able to
        // ignore it this time as well
        for (unsigned int f=0; f<GeometryInfo<structdim>::faces_per_cell; ++f)
          {
            Assert (cell->face(f)->has_children(), ExcInternalError());
            Assert (cell->face(f)->refinement_case() ==
                    RefinementCase<structdim-1>::isotropic_refinement,
                    ExcInternalError());

            bool subface_is_more_refined = false;
            for (unsigned int g=0; g<GeometryInfo<structdim>::max_children_per_face; ++g)
              if (cell->face(f)->child(g)->has_children())
                {
                  subface_is_more_refined = true;
                  break;
                }

            if (subface_is_more_refined == true)
              continue;

            // so, now we finally know
            // that we can do something
            // about this face
            fix_up_object (cell->face(f), cell->at_boundary(f));
          }
      }


    } /* namespace FixUpDistortedChildCells */
  } /* namespace internal */


  template <int dim, int spacedim>
  typename Triangulation<dim,spacedim>::DistortedCellList

  fix_up_distorted_child_cells (const typename Triangulation<dim,spacedim>::DistortedCellList &distorted_cells,
                                Triangulation<dim,spacedim> &/*triangulation*/)
  {
    typename Triangulation<dim,spacedim>::DistortedCellList unfixable_subset;

    // loop over all cells that we have
    // to fix up
    for (typename std::list<typename Triangulation<dim,spacedim>::cell_iterator>::const_iterator
         cell_ptr = distorted_cells.distorted_cells.begin();
         cell_ptr != distorted_cells.distorted_cells.end(); ++cell_ptr)
      {
        const typename Triangulation<dim,spacedim>::cell_iterator
        cell = *cell_ptr;

        internal::FixUpDistortedChildCells
        ::fix_up_faces (cell,
                        dealii::internal::int2type<dim>(),
                        dealii::internal::int2type<spacedim>());

        // fix up the object. we need to
        // respect the manifold if the cell is
        // embedded in a higher dimensional
        // space; otherwise, like a hex in 3d,
        // every point within the cell interior
        // is fair game
        if (! internal::FixUpDistortedChildCells::fix_up_object (cell,
                                                                 (dim < spacedim)))
          unfixable_subset.distorted_cells.push_back (cell);
      }

    return unfixable_subset;
  }



  template <class MeshType>
  std::vector<typename MeshType::active_cell_iterator>
  get_patch_around_cell(const typename MeshType::active_cell_iterator &cell)
  {
    Assert (cell->is_locally_owned(),
            ExcMessage ("This function only makes sense if the cell for "
                        "which you are asking for a patch, is locally "
                        "owned."));

    std::vector<typename MeshType::active_cell_iterator> patch;
    patch.push_back (cell);
    for (unsigned int face_number=0; face_number<GeometryInfo<MeshType::dimension>::faces_per_cell; ++face_number)
      if (cell->face(face_number)->at_boundary()==false)
        {
          if (cell->neighbor(face_number)->has_children() == false)
            patch.push_back (cell->neighbor(face_number));
          else
            // the neighbor is refined. in 2d/3d, we can simply ask for the children
            // of the neighbor because they can not be further refined and,
            // consequently, the children is active
            if (MeshType::dimension > 1)
              {
                for (unsigned int subface=0; subface<cell->face(face_number)->n_children(); ++subface)
                  patch.push_back (cell->neighbor_child_on_subface (face_number, subface));
              }
            else
              {
                // in 1d, we need to work a bit harder: iterate until we find
                // the child by going from cell to child to child etc
                typename MeshType::cell_iterator neighbor
                  = cell->neighbor (face_number);
                while (neighbor->has_children())
                  neighbor = neighbor->child(1-face_number);

                Assert (neighbor->neighbor(1-face_number) == cell, ExcInternalError());
                patch.push_back (neighbor);
              }
        }
    return patch;
  }



  template <class Container>
  std::vector<typename Container::cell_iterator>
  get_cells_at_coarsest_common_level (
    const std::vector<typename Container::active_cell_iterator> &patch)
  {
    Assert (patch.size() > 0, ExcMessage("Vector containing patch cells should not be an empty vector!"));
    // In order to extract the set of cells with the coarsest common level from the give vector of cells:
    // First it finds the number associated with the minimum level of refinmenet, namely "min_level"
    int min_level = patch[0]->level();

    for (unsigned int i=0; i<patch.size(); ++i)
      min_level = std::min (min_level, patch[i]->level() );
    std::set<typename Container::cell_iterator>  uniform_cells;
    typename std::vector<typename Container::active_cell_iterator>::const_iterator  patch_cell;
    // it loops through all cells of the input vector
    for (patch_cell=patch.begin(); patch_cell!=patch.end () ; ++patch_cell)
      {
        // If the refinement level of each cell i the loop be equal to the min_level, so that
        // that cell inserted into the set of uniform_cells, as the set of cells with the coarsest common refinement level
        if ((*patch_cell)->level() == min_level)
          uniform_cells.insert (*patch_cell);
        else
          // If not, it asks for the parent of the cell, until it finds the parent cell
          // with the refinement level equal to the min_level and inserts that parent cell into the
          // the set of uniform_cells, as the set of cells with the coarsest common refinement level.
          {
            typename Container::cell_iterator parent = *patch_cell;

            while (parent->level() > min_level)
              parent = parent-> parent();
            uniform_cells.insert (parent);
          }
      }

    return std::vector<typename Container::cell_iterator> (uniform_cells.begin(),
                                                           uniform_cells.end());
  }



  template <class Container>
  void build_triangulation_from_patch(const std::vector<typename Container::active_cell_iterator> &patch,
                                      Triangulation<Container::dimension,Container::space_dimension> &local_triangulation,
                                      std::map<typename Triangulation<Container::dimension,Container::space_dimension>::active_cell_iterator,
                                      typename Container::active_cell_iterator> &patch_to_global_tria_map)

  {
    const std::vector<typename Container::cell_iterator> uniform_cells =
      get_cells_at_coarsest_common_level <Container> (patch);
    // First it creates triangulation from the vector of "uniform_cells"
    local_triangulation.clear();
    std::vector<Point<Container::space_dimension> > vertices;
    const unsigned int n_uniform_cells=uniform_cells.size();
    std::vector<CellData<Container::dimension> > cells(n_uniform_cells);
    unsigned int k=0;// for enumerating cells
    unsigned int i=0;// for enumerating vertices
    typename std::vector<typename Container::cell_iterator>::const_iterator uniform_cell;
    for (uniform_cell=uniform_cells.begin(); uniform_cell!=uniform_cells.end(); ++uniform_cell)
      {
        bool repeat_vertex;
        for (unsigned int v=0; v<GeometryInfo<Container::dimension>::vertices_per_cell; ++v)
          {
            Point<Container::space_dimension> position=(*uniform_cell)->vertex (v);
            repeat_vertex=false;

            for (unsigned int m=0; m<i; ++m)
              {
                if (position == vertices[m])
                  {
                    repeat_vertex=true;
                    cells[k].vertices[v]=m;
                    break;
                  }
              }
            if (repeat_vertex==false)
              {
                vertices.push_back(position);
                cells[k].vertices[v]=i;
                i=i+1;
              }

          }//for vertices_per_cell
        k=k+1;
      }
    local_triangulation.create_triangulation(vertices,cells,SubCellData());
    Assert (local_triangulation.n_active_cells() == uniform_cells.size(), ExcInternalError());
    local_triangulation.clear_user_flags ();
    unsigned int index=0;
    // Create a map between cells of class DofHandler into class Triangulation
    std::map<typename Triangulation<Container::dimension,Container::space_dimension>::cell_iterator,
        typename Container::cell_iterator> patch_to_global_tria_map_tmp;
    for (typename Triangulation<Container::dimension,Container::space_dimension>::cell_iterator coarse_cell = local_triangulation.begin();
         coarse_cell != local_triangulation.end(); ++coarse_cell, ++index)
      {
        patch_to_global_tria_map_tmp.insert (std::make_pair(coarse_cell, uniform_cells[index]));
        // To ensure that the cells with the same coordinates (here, we compare their centers) are mapped into each other.

        Assert(coarse_cell->center().distance( uniform_cells[index]->center())<=1e-15*coarse_cell->diameter(),
               ExcInternalError());
      }
    bool refinement_necessary;
    // In this loop we start to do refinement on the above coarse triangulation to reach
    // to the same level of refinement as the patch cells are really on
    do
      {
        refinement_necessary = false;
        for (typename Triangulation<Container::dimension,Container::space_dimension>::active_cell_iterator
             active_tria_cell = local_triangulation.begin_active();
             active_tria_cell != local_triangulation.end(); ++active_tria_cell)
          {
            if (patch_to_global_tria_map_tmp[active_tria_cell]->has_children())
              {
                active_tria_cell -> set_refine_flag();
                refinement_necessary = true;
              }
            else for (unsigned int i=0; i<patch.size(); ++i)
                {
                  // Even though vertices may not be exactly the same, the
                  // appropriate cells will match since == for TriAccessors
                  // checks only cell level and index.
                  if (patch_to_global_tria_map_tmp[active_tria_cell]==patch[i])
                    {
                      // adjust the cell vertices of the local_triangulation to
                      // match cell vertices of the global triangulation
                      for (unsigned int v=0; v<GeometryInfo<Container::dimension>::vertices_per_cell; ++v)
                        active_tria_cell->vertex(v) = patch[i]->vertex(v);

                      Assert(active_tria_cell->center().distance(patch_to_global_tria_map_tmp[active_tria_cell]->center())
                             <=1e-15*active_tria_cell->diameter(), ExcInternalError());

                      active_tria_cell->set_user_flag();
                      break;
                    }
                }
          }

        if (refinement_necessary)
          {
            local_triangulation.execute_coarsening_and_refinement ();

            for (typename Triangulation<Container::dimension,Container::space_dimension>::cell_iterator
                 cell = local_triangulation.begin();
                 cell != local_triangulation.end(); ++cell)
              {

                if (patch_to_global_tria_map_tmp.find(cell)!=patch_to_global_tria_map_tmp.end())
                  {
                    if (cell-> has_children())
                      {
                        // Note: Since the cell got children, then it should not be in the map anymore
                        // children may be added into the map, instead

                        // these children may not yet be in the map
                        for (unsigned int c=0; c<cell->n_children(); ++c)
                          {
                            if (patch_to_global_tria_map_tmp.find(cell->child(c)) ==
                                patch_to_global_tria_map_tmp.end())
                              {
                                patch_to_global_tria_map_tmp.insert (std::make_pair(cell->child(c),
                                                                                    patch_to_global_tria_map_tmp[cell]->child(c)));

                                // One might be tempted to assert that the cell
                                // being added here has the same center as the
                                // equivalent cell in the global triangulation,
                                // but it may not be the case.  For triangulations
                                // that have been perturbed or smoothed, the cell
                                // indices and levels may be the same, but the
                                // vertex locations may not.  We adjust
                                // the vertices of the cells that have no
                                // children (ie the active cells) to be
                                // consistent with the global triangulation
                                // later on and add assertions at that time
                                // to guarantee the cells in the
                                // local_triangulation are physically at the same
                                // locations of the cells in the patch of the
                                // global triangulation.

                              }
                          }
                        // The parent cell whose children were added
                        // into the map should be deleted from the map
                        patch_to_global_tria_map_tmp.erase(cell);
                      }
                  }
              }
          }

      }
    while (refinement_necessary);


    // Last assertion check to make sure we have the right cells and centers
    // in the map, and hence the correct vertices of the triangulation
    for (typename Triangulation<Container::dimension,Container::space_dimension>::cell_iterator
         cell = local_triangulation.begin();
         cell != local_triangulation.end(); ++cell)
      {
        if (cell->user_flag_set() )
          {
            Assert(patch_to_global_tria_map_tmp.find(cell) != patch_to_global_tria_map_tmp.end(),
                   ExcInternalError() );

            Assert(cell->center().distance( patch_to_global_tria_map_tmp[cell]->center())<=1e-15*cell->diameter(),
                   ExcInternalError());
          }
      }


    typename std::map<typename Triangulation<Container::dimension,Container::space_dimension>::cell_iterator,
             typename Container::cell_iterator>::iterator map_tmp_it =
               patch_to_global_tria_map_tmp.begin(),map_tmp_end = patch_to_global_tria_map_tmp.end();
    // Now we just need to take the temporary map of pairs of type cell_iterator "patch_to_global_tria_map_tmp"
    // making pair of active_cell_iterators so that filling out the final map "patch_to_global_tria_map"
    for (; map_tmp_it!=map_tmp_end; ++map_tmp_it)
      patch_to_global_tria_map[map_tmp_it->first] = map_tmp_it->second;
  }




  template <class DoFHandlerType>
  std::map< types::global_dof_index,std::vector<typename DoFHandlerType::active_cell_iterator> >
  get_dof_to_support_patch_map(DoFHandlerType &dof_handler)
  {

    // This is the map from global_dof_index to
    // a set of cells on patch.  We first map into
    // a set because it is very likely that we
    // will attempt to add a cell more than once
    // to a particular patch and we want to preserve
    // uniqueness of cell iterators. std::set does this
    // automatically for us.  Later after it is all
    // constructed, we will copy to a map of vectors
    // since that is the prefered output for other
    // functions.
    std::map< types::global_dof_index,std::set<typename DoFHandlerType::active_cell_iterator> > dof_to_set_of_cells_map;

    std::vector<types::global_dof_index> local_dof_indices;
    std::vector<types::global_dof_index> local_face_dof_indices;
    std::vector<types::global_dof_index> local_line_dof_indices;

    // a place to save the dof_handler user flags and restore them later
    // to maintain const of dof_handler.
    std::vector<bool> user_flags;


    // in 3d, we need pointers from active lines to the
    // active parent lines, so we construct it as needed.
    std::map<typename DoFHandlerType::active_line_iterator, typename DoFHandlerType::line_iterator > lines_to_parent_lines_map;
    if (DoFHandlerType::dimension == 3)
      {

        // save user flags as they will be modified and then later restored
        dof_handler.get_triangulation().save_user_flags(user_flags);
        const_cast<dealii::Triangulation<DoFHandlerType::dimension,DoFHandlerType::space_dimension> &>(dof_handler.get_triangulation()).clear_user_flags ();


        typename DoFHandlerType::active_cell_iterator cell = dof_handler.begin_active(),
                                                      endc = dof_handler.end();
        for (; cell!=endc; ++cell)
          {
            // We only want lines that are locally_relevant
            // although it doesn't hurt to have lines that
            // are children of ghost cells since there are
            // few and we don't have to use them.
            if (cell->is_artificial() == false)
              {
                for (unsigned int l=0; l<GeometryInfo<DoFHandlerType::dimension>::lines_per_cell; ++l)
                  if (cell->line(l)->has_children())
                    for (unsigned int c=0; c<cell->line(l)->n_children(); ++c)
                      {
                        lines_to_parent_lines_map[cell->line(l)->child(c)] = cell->line(l);
                        // set flags to know that child
                        // line has an active parent.
                        cell->line(l)->child(c)->set_user_flag();
                      }
              }
          }
      }


    // We loop through all cells and add cell to the
    // map for the dofs that it immediately touches
    // and then account for all the other dofs of
    // which it is a part, mainly the ones that must
    // be added on account of adaptivity hanging node
    // constraints.
    typename DoFHandlerType::active_cell_iterator cell = dof_handler.begin_active(),
                                                  endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        // Need to loop through all cells that could
        // be in the patch of dofs on locally_owned
        // cells including ghost cells
        if (cell->is_artificial() == false)
          {
            const unsigned int n_dofs_per_cell = cell->get_fe().dofs_per_cell;
            local_dof_indices.resize(n_dofs_per_cell);

            // Take care of adding cell pointer to each
            // dofs that exists on cell.
            cell->get_dof_indices(local_dof_indices);
            for (unsigned int i=0; i< n_dofs_per_cell; ++i )
              dof_to_set_of_cells_map[local_dof_indices[i]].insert(cell);

            // In the case of the adjacent cell (over
            // faces or edges) being more refined, we
            // want to add all of the children to the
            // patch since the support function at that
            // dof could be non-zero along that entire
            // face (or line).

            // Take care of dofs on neighbor faces
            for (unsigned int f=0; f<GeometryInfo<DoFHandlerType::dimension>::faces_per_cell; ++f)
              {
                if (cell->face(f)->has_children())
                  {
                    for (unsigned int c=0; c<cell->face(f)->n_children(); ++c)
                      {
                        //  Add cell to dofs of all subfaces
                        //
                        //   *-------------------*----------*---------*
                        //   |                   | add cell |         |
                        //   |                   |<- to dofs|         |
                        //   |                   |of subface|         |
                        //   |        cell       *----------*---------*
                        //   |                   | add cell |         |
                        //   |                   |<- to dofs|         |
                        //   |                   |of subface|         |
                        //   *-------------------*----------*---------*
                        //
                        Assert (cell->face(f)->child(c)->has_children() == false, ExcInternalError());

                        const unsigned int n_dofs_per_face = cell->get_fe().dofs_per_face;
                        local_face_dof_indices.resize(n_dofs_per_face);

                        cell->face(f)->child(c)->get_dof_indices(local_face_dof_indices);
                        for (unsigned int i=0; i< n_dofs_per_face; ++i )
                          dof_to_set_of_cells_map[local_face_dof_indices[i]].insert(cell);
                      }
                  }
                else if ((cell->face(f)->at_boundary() == false) && (cell->neighbor_is_coarser(f)))
                  {

                    // Add cell to dofs of parent face and all
                    // child faces of parent face
                    //
                    //   *-------------------*----------*---------*
                    //   |                   |          |         |
                    //   |                   |   cell   |         |
                    //   |      add cell     |          |         |
                    //   |      to dofs   -> *----------*---------*
                    //   |      of parent    | add cell |         |
                    //   |       face        |<- to dofs|         |
                    //   |                   |of subface|         |
                    //   *-------------------*----------*---------*
                    //

                    // Add cell to all dofs of parent face
                    std::pair<unsigned int, unsigned int> neighbor_face_no_subface_no = cell->neighbor_of_coarser_neighbor(f);
                    unsigned int face_no = neighbor_face_no_subface_no.first;
                    unsigned int subface = neighbor_face_no_subface_no.second;

                    const unsigned int n_dofs_per_face = cell->get_fe().dofs_per_face;
                    local_face_dof_indices.resize(n_dofs_per_face);

                    cell->neighbor(f)->face(face_no)->get_dof_indices(local_face_dof_indices);
                    for (unsigned int i=0; i< n_dofs_per_face; ++i )
                      dof_to_set_of_cells_map[local_face_dof_indices[i]].insert(cell);

                    // Add cell to all dofs of children of
                    // parent face
                    for (unsigned int c=0; c<cell->neighbor(f)->face(face_no)->n_children(); ++c)
                      {
                        if (c != subface) // don't repeat work on dofs of original cell
                          {
                            const unsigned int n_dofs_per_face = cell->get_fe().dofs_per_face;
                            local_face_dof_indices.resize(n_dofs_per_face);

                            Assert (cell->neighbor(f)->face(face_no)->child(c)->has_children() == false, ExcInternalError());
                            cell->neighbor(f)->face(face_no)->child(c)->get_dof_indices(local_face_dof_indices);
                            for (unsigned int i=0; i<n_dofs_per_face; ++i )
                              dof_to_set_of_cells_map[local_face_dof_indices[i]].insert(cell);
                          }
                      }
                  }
              }


            // If 3d, take care of dofs on lines in the
            // same pattern as faces above. That is, if
            // a cell's line has children, distribute
            // cell to dofs of children of line,  and
            // if cell's line has an active parent, then
            // distribute cell to dofs on parent line
            // and dofs on all children of parent line.
            if (DoFHandlerType::dimension == 3)
              {
                for (unsigned int l=0; l<GeometryInfo<DoFHandlerType::dimension>::lines_per_cell; ++l)
                  {
                    if (cell->line(l)->has_children())
                      {
                        for (unsigned int c=0; c<cell->line(l)->n_children(); ++c)
                          {
                            Assert (cell->line(l)->child(c)->has_children() == false, ExcInternalError());

                            // dofs_per_line returns number of dofs
                            // on line not including the vertices of the line.
                            const unsigned int n_dofs_per_line = 2*cell->get_fe().dofs_per_vertex
                                                                 + cell->get_fe().dofs_per_line;
                            local_line_dof_indices.resize(n_dofs_per_line);

                            cell->line(l)->child(c)->get_dof_indices(local_line_dof_indices);
                            for (unsigned int i=0; i<n_dofs_per_line; ++i )
                              dof_to_set_of_cells_map[local_line_dof_indices[i]].insert(cell);
                          }
                      }
                    // user flag was set above to denote that
                    // an active parent line exists so add
                    // cell to dofs of parent and all it's
                    // children
                    else if (cell->line(l)->user_flag_set() == true)
                      {
                        typename DoFHandlerType::line_iterator parent_line = lines_to_parent_lines_map[cell->line(l)];
                        Assert (parent_line->has_children(), ExcInternalError() );

                        // dofs_per_line returns number of dofs
                        // on line not including the vertices of the line.
                        const unsigned int n_dofs_per_line = 2*cell->get_fe().dofs_per_vertex
                                                             + cell->get_fe().dofs_per_line;
                        local_line_dof_indices.resize(n_dofs_per_line);

                        parent_line->get_dof_indices(local_line_dof_indices);
                        for (unsigned int i=0; i<n_dofs_per_line; ++i )
                          dof_to_set_of_cells_map[local_line_dof_indices[i]].insert(cell);

                        for (unsigned int c=0; c<parent_line->n_children(); ++c)
                          {
                            Assert (parent_line->child(c)->has_children() == false, ExcInternalError());

                            const unsigned int n_dofs_per_line = 2*cell->get_fe().dofs_per_vertex
                                                                 + cell->get_fe().dofs_per_line;
                            local_line_dof_indices.resize(n_dofs_per_line);

                            parent_line->child(c)->get_dof_indices(local_line_dof_indices);
                            for (unsigned int i=0; i<n_dofs_per_line; ++i )
                              dof_to_set_of_cells_map[local_line_dof_indices[i]].insert(cell);
                          }


                      }
                  } // for lines l
              }// if DoFHandlerType::dimension == 3
          }// if cell->is_locally_owned()
      }// for cells


    if (DoFHandlerType::dimension == 3)
      {
        // finally, restore user flags that were changed above
        // to when we constructed the pointers to parent of lines
        // Since dof_handler is const, we must leave it unchanged.
        const_cast<dealii::Triangulation<DoFHandlerType::dimension,DoFHandlerType::space_dimension> &>(dof_handler.get_triangulation()).load_user_flags (user_flags);
      }

    // Finally, we copy map of sets to
    // map of vectors using the std::vector::assign() function
    std::map< types::global_dof_index, std::vector<typename DoFHandlerType::active_cell_iterator> > dof_to_cell_patches;

    typename std::map<types::global_dof_index, std::set< typename DoFHandlerType::active_cell_iterator> >::iterator
    it = dof_to_set_of_cells_map.begin(),
    it_end = dof_to_set_of_cells_map.end();
    for ( ; it!=it_end; ++it)
      dof_to_cell_patches[it->first].assign( it->second.begin(), it->second.end() );

    return dof_to_cell_patches;
  }



  /*
   * Internally used in orthogonal_equality
   *
   * An orthogonal equality test for points:
   *
   * point1 and point2 are considered equal, if
   *   matrix.point1 + offset - point2
   * is parallel to the unit vector in <direction>
   */
  template<int spacedim>
  inline bool orthogonal_equality (const Point<spacedim>    &point1,
                                   const Point<spacedim>    &point2,
                                   const int                 direction,
                                   const Tensor<1,spacedim> &offset,
                                   const FullMatrix<double> &matrix)
  {
    Assert (0<=direction && direction<spacedim,
            ExcIndexRange (direction, 0, spacedim));

    Assert(matrix.m() == matrix.n(), ExcInternalError());

    Point<spacedim> distance;

    if (matrix.m() == spacedim)
      for (int i = 0; i < spacedim; ++i)
        for (int j = 0; j < spacedim; ++j)
          distance(i) += matrix(i,j) * point1(j);
    else
      distance = point1;

    distance += offset - point2;

    for (int i = 0; i < spacedim; ++i)
      {
        // Only compare coordinate-components != direction:
        if (i == direction)
          continue;

        if (fabs(distance(i)) > 1.e-10)
          return false;
      }

    return true;
  }


  /*
   * Internally used in orthogonal_equality
   *
   * A lookup table to transform vertex matchings to orientation flags of
   * the form (face_orientation, face_flip, face_rotation)
   *
   * See the comment on the next function as well as the detailed
   * documentation of make_periodicity_constraints and
   * collect_periodic_faces for details
   */
  template<int dim> struct OrientationLookupTable {};

  template<> struct OrientationLookupTable<1>
  {
    typedef std_cxx11::array<unsigned int, GeometryInfo<1>::vertices_per_face> MATCH_T;
    static inline std::bitset<3> lookup (const MATCH_T &)
    {
      // The 1D case is trivial
      return 1; // [true ,false,false]
    }
  };

  template<> struct OrientationLookupTable<2>
  {
    typedef std_cxx11::array<unsigned int, GeometryInfo<2>::vertices_per_face> MATCH_T;
    static inline std::bitset<3> lookup (const MATCH_T &matching)
    {
      // In 2D matching faces (=lines) results in two cases: Either
      // they are aligned or flipped. We store this "line_flip"
      // property somewhat sloppy as "face_flip"
      // (always: face_orientation = true, face_rotation = false)

      static const MATCH_T m_tff = {{ 0 , 1 }};
      if (matching == m_tff) return 1;           // [true ,false,false]
      static const MATCH_T m_ttf = {{ 1 , 0 }};
      if (matching == m_ttf) return 3;           // [true ,true ,false]
      AssertThrow(false, ExcInternalError());
      // what follows is dead code, but it avoids warnings about the lack
      // of a return value
      return 0;
    }
  };

  template<> struct OrientationLookupTable<3>
  {
    typedef std_cxx11::array<unsigned int, GeometryInfo<3>::vertices_per_face> MATCH_T;
    static inline std::bitset<3> lookup (const MATCH_T &matching)
    {
      // The full fledged 3D case. *Yay*
      // See the documentation in include/deal.II/base/geometry_info.h
      // as well as the actual implementation in source/grid/tria.cc
      // for more details...

      static const MATCH_T m_tff = {{ 0 , 1 , 2 , 3 }};
      if (matching == m_tff) return 1;                   // [true ,false,false]
      static const MATCH_T m_tft = {{ 1 , 3 , 0 , 2 }};
      if (matching == m_tft) return 5;                   // [true ,false,true ]
      static const MATCH_T m_ttf = {{ 3 , 2 , 1 , 0 }};
      if (matching == m_ttf) return 3;                   // [true ,true ,false]
      static const MATCH_T m_ttt = {{ 2 , 0 , 3 , 1 }};
      if (matching == m_ttt) return 7;                   // [true ,true ,true ]
      static const MATCH_T m_fff = {{ 0 , 2 , 1 , 3 }};
      if (matching == m_fff) return 0;                   // [false,false,false]
      static const MATCH_T m_fft = {{ 2 , 3 , 0 , 1 }};
      if (matching == m_fft) return 4;                   // [false,false,true ]
      static const MATCH_T m_ftf = {{ 3 , 1 , 2 , 0 }};
      if (matching == m_ftf) return 2;                   // [false,true ,false]
      static const MATCH_T m_ftt = {{ 1 , 0 , 3 , 2 }};
      if (matching == m_ftt) return 6;                   // [false,true ,true ]
      AssertThrow(false, ExcInternalError());
      // what follows is dead code, but it avoids warnings about the lack
      // of a return value
      return 0;
    }
  };



  template<typename FaceIterator>
  inline bool
  orthogonal_equality (std::bitset<3>     &orientation,
                       const FaceIterator &face1,
                       const FaceIterator &face2,
                       const int          direction,
                       const Tensor<1,FaceIterator::AccessorType::space_dimension> &offset,
                       const FullMatrix<double> &matrix)
  {
    Assert(matrix.m() == matrix.n(),
           ExcMessage("The supplied matrix must be a square matrix"));

    static const int dim = FaceIterator::AccessorType::dimension;

    // Do a full matching of the face vertices:

    std_cxx11::
    array<unsigned int, GeometryInfo<dim>::vertices_per_face> matching;

    std::set<unsigned int> face2_vertices;
    for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_face; ++i)
      face2_vertices.insert(i);

    for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_face; ++i)
      {
        for (std::set<unsigned int>::iterator it = face2_vertices.begin();
             it != face2_vertices.end();
             it++)
          {
            if (orthogonal_equality(face1->vertex(i),face2->vertex(*it),
                                    direction, offset, matrix))
              {
                matching[i] = *it;
                face2_vertices.erase(it);
                break; // jump out of the innermost loop
              }
          }
      }

    // And finally, a lookup to determine the ordering bitmask:
    if (face2_vertices.empty())
      orientation = OrientationLookupTable<dim>::lookup(matching);

    return face2_vertices.empty();
  }



  template<typename FaceIterator>
  inline bool
  orthogonal_equality (const FaceIterator &face1,
                       const FaceIterator &face2,
                       const int          direction,
                       const Tensor<1,FaceIterator::AccessorType::space_dimension> &offset,
                       const FullMatrix<double> &matrix)
  {
    // Call the function above with a dummy orientation array
    std::bitset<3> dummy;
    return orthogonal_equality (dummy, face1, face2, direction, offset, matrix);
  }



  /*
   * Internally used in collect_periodic_faces
   */
  template<typename CellIterator>
  void
  match_periodic_face_pairs
  (std::set<std::pair<CellIterator, unsigned int> > &pairs1,
   std::set<std::pair<typename identity<CellIterator>::type, unsigned int> > &pairs2,
   const int                                        direction,
   std::vector<PeriodicFacePair<CellIterator> >     &matched_pairs,
   const dealii::Tensor<1,CellIterator::AccessorType::space_dimension> &offset,
   const FullMatrix<double>                         &matrix)
  {
    static const int space_dim = CellIterator::AccessorType::space_dimension;
    (void)space_dim;
    Assert (0<=direction && direction<space_dim,
            ExcIndexRange (direction, 0, space_dim));

    Assert (pairs1.size() == pairs2.size(),
            ExcMessage ("Unmatched faces on periodic boundaries"));

    unsigned int n_matches = 0;

    // Match with a complexity of O(n^2). This could be improved...
    std::bitset<3> orientation;
    typedef typename std::set
    <std::pair<CellIterator, unsigned int> >::const_iterator PairIterator;
    for (PairIterator it1 = pairs1.begin(); it1 != pairs1.end(); ++it1)
      {
        for (PairIterator it2 = pairs2.begin(); it2 != pairs2.end(); ++it2)
          {
            const CellIterator cell1 = it1->first;
            const CellIterator cell2 = it2->first;
            const unsigned int face_idx1 = it1->second;
            const unsigned int face_idx2 = it2->second;
            if (GridTools::orthogonal_equality(orientation,
                                               cell1->face(face_idx1),
                                               cell2->face(face_idx2),
                                               direction, offset,
                                               matrix))
              {
                // We have a match, so insert the matching pairs and
                // remove the matched cell in pairs2 to speed up the
                // matching:
                const PeriodicFacePair<CellIterator> matched_face =
                {
                  {cell1, cell2},
                  {face_idx1, face_idx2},
                  orientation,
                  matrix
                };
                matched_pairs.push_back(matched_face);
                pairs2.erase(it2);
                ++n_matches;
                break;
              }
          }
      }

    //Assure that all faces are matched
    AssertThrow (n_matches == pairs1.size() && pairs2.size() == 0,
                 ExcMessage ("Unmatched faces on periodic boundaries"));
  }



  template<typename MeshType>
  void
  collect_periodic_faces
  (const MeshType                        &mesh,
   const types::boundary_id               b_id1,
   const types::boundary_id               b_id2,
   const int                              direction,
   std::vector<PeriodicFacePair<typename MeshType::cell_iterator> > &matched_pairs,
   const Tensor<1,MeshType::space_dimension> &offset,
   const FullMatrix<double>              &matrix)
  {
    static const int dim = MeshType::dimension;
    static const int space_dim = MeshType::space_dimension;
    (void)dim;
    (void)space_dim;
    Assert (0<=direction && direction<space_dim,
            ExcIndexRange (direction, 0, space_dim));

    // Loop over all cells on the highest level and collect all boundary
    // faces belonging to b_id1 and b_id2:

    std::set<std::pair<typename MeshType::cell_iterator, unsigned int> > pairs1;
    std::set<std::pair<typename MeshType::cell_iterator, unsigned int> > pairs2;

    for (typename MeshType::cell_iterator cell = mesh.begin(0);
         cell != mesh.end(0); ++cell)
      {
        for (unsigned int i = 0; i < GeometryInfo<dim>::faces_per_cell; ++i)
          {
            const typename MeshType::face_iterator face = cell->face(i);
            if (face->at_boundary() && face->boundary_id() == b_id1)
              {
                const std::pair<typename MeshType::cell_iterator, unsigned int> pair1
                  = std::make_pair(cell, i);
                pairs1.insert(pair1);
              }

            if (face->at_boundary() && face->boundary_id() == b_id2)
              {
                const std::pair<typename MeshType::cell_iterator, unsigned int> pair2
                  = std::make_pair(cell, i);
                pairs2.insert(pair2);
              }
          }
      }

    Assert (pairs1.size() == pairs2.size(),
            ExcMessage ("Unmatched faces on periodic boundaries"));

    // and call match_periodic_face_pairs that does the actual matching:
    match_periodic_face_pairs(pairs1, pairs2, direction, matched_pairs, offset,
                              matrix);
  }



  template<typename MeshType>
  void
  collect_periodic_faces
  (const MeshType                        &mesh,
   const types::boundary_id               b_id,
   const int                              direction,
   std::vector<PeriodicFacePair<typename MeshType::cell_iterator> > &matched_pairs,
   const Tensor<1,MeshType::space_dimension> &offset,
   const FullMatrix<double>              &matrix)
  {
    static const int dim = MeshType::dimension;
    static const int space_dim = MeshType::space_dimension;
    (void)dim;
    (void)space_dim;
    Assert (0<=direction && direction<space_dim,
            ExcIndexRange (direction, 0, space_dim));

    Assert(dim == space_dim,
           ExcNotImplemented());

    // Loop over all cells on the highest level and collect all boundary
    // faces 2*direction and 2*direction*1:

    std::set<std::pair<typename MeshType::cell_iterator, unsigned int> > pairs1;
    std::set<std::pair<typename MeshType::cell_iterator, unsigned int> > pairs2;

    for (typename MeshType::cell_iterator cell = mesh.begin(0);
         cell != mesh.end(0); ++cell)
      {
        const typename MeshType::face_iterator face_1 = cell->face(2*direction);
        const typename MeshType::face_iterator face_2 = cell->face(2*direction+1);

        if (face_1->at_boundary() && face_1->boundary_id() == b_id)
          {
            const std::pair<typename MeshType::cell_iterator, unsigned int> pair1
              = std::make_pair(cell, 2*direction);
            pairs1.insert(pair1);
          }

        if (face_2->at_boundary() && face_2->boundary_id() == b_id)
          {
            const std::pair<typename MeshType::cell_iterator, unsigned int> pair2
              = std::make_pair(cell, 2*direction+1);
            pairs2.insert(pair2);
          }
      }

    Assert (pairs1.size() == pairs2.size(),
            ExcMessage ("Unmatched faces on periodic boundaries"));


#ifdef DEBUG
    const unsigned int size_old = matched_pairs.size();
#endif

    // and call match_periodic_face_pairs that does the actual matching:
    match_periodic_face_pairs(pairs1, pairs2, direction, matched_pairs, offset,
                              matrix);

#ifdef DEBUG
    //check for standard orientation
    const unsigned int size_new = matched_pairs.size();
    for (unsigned int i = size_old; i < size_new; ++i)
      {
        Assert(matched_pairs[i].orientation == 1,
               ExcMessage("Found a face match with non standard orientation. "
                          "This function is only suitable for meshes with cells "
                          "in default orientation"));
      }
#endif
  }



  template <int dim, int spacedim>
  void copy_boundary_to_manifold_id(Triangulation<dim, spacedim> &tria,
                                    const bool reset_boundary_ids)
  {
    // in 3d, we not only have to copy boundary ids of faces, but also of edges
    // because we see them twice (once from each adjacent boundary face),
    // we cannot immediately reset their boundary ids. thus, copy first
    // and reset later
    if (dim >= 3)
      for (typename Triangulation<dim,spacedim>::active_cell_iterator
           cell=tria.begin_active();
           cell != tria.end(); ++cell)
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          if (cell->face(f)->at_boundary())
            for (signed int e=0; e<static_cast<signed int>(GeometryInfo<dim>::lines_per_face); ++e)
              cell->face(f)->line(e)->set_manifold_id
              (static_cast<types::manifold_id>(cell->face(f)->line(e)->boundary_id()));

    // now do cells
    for (typename Triangulation<dim,spacedim>::active_cell_iterator
         cell=tria.begin_active();
         cell != tria.end(); ++cell)
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        if (cell->face(f)->at_boundary())
          {
            // copy boundary to manifold ids
            cell->face(f)->set_manifold_id
            (static_cast<types::manifold_id>(cell->face(f)->boundary_id()));

            // then reset boundary ids if so desired, and in 3d also that
            // of edges
            if (reset_boundary_ids == true)
              {
                cell->face(f)->set_boundary_id(0);
                if (dim >= 3)
                  for (signed int e=0; e<static_cast<signed int>(GeometryInfo<dim>::lines_per_face); ++e)
                    cell->face(f)->line(e)->set_boundary_id(0);
              }
          }
  }



  template <int dim, int spacedim>
  void copy_material_to_manifold_id(Triangulation<dim, spacedim> &tria,
                                    const bool compute_face_ids)
  {
    typename Triangulation<dim,spacedim>::active_cell_iterator
    cell=tria.begin_active(), endc=tria.end();

    for (; cell != endc; ++cell)
      {
        cell->set_manifold_id(cell->material_id());
        if (compute_face_ids == true)
          {
            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
              {
                if (cell->at_boundary(f) == false)
                  cell->face(f)->set_manifold_id
                  (std::min(cell->material_id(),
                            cell->neighbor(f)->material_id()));
                else
                  cell->face(f)->set_manifold_id(cell->material_id());
              }
          }
      }
  }

  template<int dim, int spacedim>
  std::pair<unsigned int, double>
  get_longest_direction(typename Triangulation<dim, spacedim>::active_cell_iterator cell)
  {
    double max_ratio = 1;
    unsigned int index = 0;

    for (unsigned int i = 0; i < dim; ++i)
      for (unsigned int j = i+1; j < dim; ++j)
        {
          unsigned int ax = i % dim;
          unsigned int next_ax = j % dim;

          double ratio =  cell->extent_in_direction(ax) / cell->extent_in_direction(next_ax);

          if ( ratio > max_ratio )
            {
              max_ratio = ratio;
              index = ax;
            }
          else if ( 1.0 /ratio > max_ratio )
            {
              max_ratio = 1.0 /ratio;
              index = next_ax;
            }
        }
    return std::make_pair(index, max_ratio);
  }


  template<int dim, int spacedim>
  void
  remove_hanging_nodes( Triangulation<dim,spacedim> &tria,
                        const bool isotropic,
                        const unsigned int max_iterations)
  {
    unsigned int iter = 0;
    bool continue_refinement = true;

    typename Triangulation<dim, spacedim>::active_cell_iterator
    cell = tria.begin_active(),
    endc = tria.end();

    while ( continue_refinement && (iter < max_iterations) )
      {
        if (max_iterations != numbers::invalid_unsigned_int) iter++;
        continue_refinement = false;

        for (cell=tria.begin_active(); cell!= endc; ++cell)
          for (unsigned int j = 0; j < GeometryInfo<dim>::faces_per_cell; j++)
            if (cell->at_boundary(j)==false && cell->neighbor(j)->has_children())
              {
                if (isotropic)
                  {
                    cell->set_refine_flag();
                    continue_refinement = true;
                  }
                else
                  continue_refinement |= cell->flag_for_face_refinement(j);
              }

        tria.execute_coarsening_and_refinement();
      }
  }

  template<int dim, int spacedim>
  void
  remove_anisotropy(  Triangulation<dim,spacedim> &tria,
                      const double max_ratio,
                      const unsigned int max_iterations)
  {
    unsigned int iter = 0;
    bool continue_refinement = true;

    typename Triangulation<dim, spacedim>::active_cell_iterator
    cell = tria.begin_active(),
    endc = tria.end();

    while ( continue_refinement && (iter<max_iterations) )
      {
        iter++;
        continue_refinement = false;
        for (cell=tria.begin_active(); cell!= endc; ++cell)
          {
            std::pair<unsigned int, double> info = GridTools::get_longest_direction<dim,spacedim>(cell);
            if (info.second > max_ratio)
              {
                cell->set_refine_flag(RefinementCase<dim>::cut_axis(info.first));
                continue_refinement = true;
              }
          }
        tria.execute_coarsening_and_refinement ();
      }
  }

} /* namespace GridTools */


// explicit instantiations
#include "grid_tools.inst"

DEAL_II_NAMESPACE_CLOSE
