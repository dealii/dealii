// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2018 by the deal.II authors
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

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi.templates.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/shared_tria.h>

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q_generic.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/manifold.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/filtered_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.h>

#include <deal.II/numerics/matrix_tools.h>

#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

#include <array>
#include <cmath>
#include <numeric>
#include <list>
#include <set>
#include <tuple>
#include <unordered_map>
#include <iostream>

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
    unsigned int mapping_degree = 1;
    if (const auto *p = dynamic_cast<const MappingQGeneric<dim,spacedim>*>(&mapping))
      mapping_degree = p->get_degree();
    else if (const auto *p = dynamic_cast<const MappingQ<dim,spacedim>*>(&mapping))
      mapping_degree = p->get_degree();

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
  cell_measure<1>
  (const std::vector<Point<1> > &all_vertices,
   const unsigned int (&vertex_indices)[GeometryInfo<1>::vertices_per_cell])
  {
    return all_vertices[vertex_indices[1]][0]
           - all_vertices[vertex_indices[0]][0];
  }



  template <>
  double
  cell_measure<3>
  (const std::vector<Point<3> > &all_vertices,
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
  cell_measure<2>
  (const std::vector<Point<2> > &all_vertices,
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



  template <int dim, int spacedim>
  BoundingBox<spacedim>
  compute_bounding_box(const Triangulation<dim, spacedim> &tria)
  {
    using iterator = typename Triangulation<dim, spacedim>::active_cell_iterator;
    const auto predicate = [](const iterator &)
    {
      return true;
    };

    return compute_bounding_box(tria, std::function<bool(const iterator &)>(predicate));
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

  template <int dim>
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
    void laplace_solve (const SparseMatrix<double>                     &S,
                        const std::map<types::global_dof_index,double> &fixed_dofs,
                        Vector<double>                                 &u)
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

      SF.add_constraints(fixed_dofs);
      SF.apply_constraints (f, true);
      solver.solve(SF, u, f, PF);
    }
  }



  // Implementation for 1D only
  template <>
  void laplace_transform (const std::map<unsigned int,Point<1> > &,
                          Triangulation<1> &,
                          const Function<1> *,
                          const bool )
  {
    Assert(false, ExcNotImplemented());
  }


  // Implementation for dimensions except 1
  template <int dim>
  void
  laplace_transform (const std::map<unsigned int,Point<dim> > &new_points,
                     Triangulation<dim> &triangulation,
                     const Function<dim> *coefficient,
                     const bool solve_for_absolute_positions)
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

    // set up the boundary values for the laplace problem
    std::map<types::global_dof_index,double> fixed_dofs[dim];
    typename std::map<unsigned int,Point<dim> >::const_iterator map_end=new_points.end();

    // fill these maps using the data given by new_points
    typename DoFHandler<dim>::cell_iterator cell=dof_handler.begin_active(),
                                            endc=dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        // loop over all vertices of the cell and see if it is listed in the map
        // given as first argument of the function
        for (unsigned int vertex_no=0;
             vertex_no<GeometryInfo<dim>::vertices_per_cell; ++vertex_no)
          {
            const unsigned int vertex_index=cell->vertex_index(vertex_no);
            const Point<dim> &vertex_point=cell->vertex(vertex_no);

            const typename std::map<unsigned int,Point<dim> >::const_iterator map_iter
              = new_points.find(vertex_index);

            if (map_iter!=map_end)
              for (unsigned int i=0; i<dim; ++i)
                fixed_dofs[i].insert(std::pair<types::global_dof_index,double>
                                     (cell->vertex_dof_index(vertex_no, 0),
                                      (solve_for_absolute_positions ?
                                       map_iter->second(i) :
                                       map_iter->second(i) - vertex_point[i])
                                     ));
          }
      }

    // solve the dim problems with different right hand sides.
    Vector<double> us[dim];
    for (unsigned int i=0; i<dim; ++i)
      us[i].reinit (dof_handler.n_dofs());

    // solve linear systems in parallel
    Threads::TaskGroup<> tasks;
    for (unsigned int i=0; i<dim; ++i)
      tasks += Threads::new_task (&laplace_solve,
                                  S, fixed_dofs[i], us[i]);
    tasks.join_all ();

    // change the coordinates of the points of the triangulation
    // according to the computed values
    std::vector<bool> vertex_touched (triangulation.n_vertices(), false);
    for (cell=dof_handler.begin_active(); cell!=endc; ++cell)
      for (unsigned int vertex_no=0;
           vertex_no<GeometryInfo<dim>::vertices_per_cell; ++vertex_no)
        if (vertex_touched[cell->vertex_index(vertex_no)] == false)
          {
            Point<dim> &v = cell->vertex(vertex_no);

            const types::global_dof_index dof_index = cell->vertex_dof_index(vertex_no, 0);
            for (unsigned int i=0; i<dim; ++i)
              if (solve_for_absolute_positions)
                v(i) = us[i](dof_index);
              else
                v(i) += us[i](dof_index);

            vertex_touched[cell->vertex_index(vertex_no)] = true;
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
    std::vector<bool>   at_boundary (keep_boundary ? triangulation.n_vertices() : 0, false);
    // for parallel::shared::Triangulation we need to work on all vertices,
    // not just the ones related to loacally owned cells;
    const bool is_parallel_shared
      = (dynamic_cast<parallel::shared::Triangulation<dim,spacedim>*> (&triangulation) != nullptr);
    for (typename Triangulation<dim,spacedim>::active_cell_iterator
         cell=triangulation.begin_active(); cell!=triangulation.end(); ++cell)
      if (is_parallel_shared || cell->is_locally_owned())
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
                       const Point<spacedim>        &p,
                       const std::vector<bool>      &marked_vertices)
  {
    // first get the underlying
    // triangulation from the
    // mesh and determine vertices
    // and used vertices
    const Triangulation<dim, spacedim> &tria = mesh.get_triangulation();

    const std::vector< Point<spacedim> > &vertices = tria.get_vertices();

    Assert ( tria.get_vertices().size() == marked_vertices.size() || marked_vertices.size() ==0,
             ExcDimensionMismatch(tria.get_vertices().size(), marked_vertices.size()));

    // If p is an element of marked_vertices,
    // and q is that of used_Vertices,
    // the vector marked_vertices does NOT
    // contain unused vertices if p implies q.
    // I.e., if p is true q must be true
    // (if p is false, q could be false or true).
    // p implies q logic is encapsulated in ~p|q.
    Assert( marked_vertices.size()==0
            ||
            std::equal( marked_vertices.begin(),
                        marked_vertices.end(),
                        tria.get_used_vertices().begin(),
                        [](bool p, bool q)
    {
      return !p || q;
    }),
    ExcMessage("marked_vertices should be a subset of used vertices in the triangulation "
               "but marked_vertices contains one or more vertices that are not used vertices!") );

    // In addition, if a vector bools
    // is specified (marked_vertices)
    // marking all the vertices which
    // could be the potentially closest
    // vertex to the point, use it instead
    // of used vertices
    const std::vector<bool> &used     =
      (marked_vertices.size()==0) ? tria.get_used_vertices() : marked_vertices;

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



  template <int dim, template <int, int> class MeshType, int spacedim>
  unsigned int
  find_closest_vertex (const Mapping<dim,spacedim>  &mapping,
                       const MeshType<dim,spacedim> &mesh,
                       const Point<spacedim>        &p,
                       const std::vector<bool>      &marked_vertices)
  {
    // Take a shortcut in the simple case.
    if (mapping.preserves_vertex_locations() == true)
      return find_closest_vertex(mesh, p, marked_vertices);

    // first get the underlying
    // triangulation from the
    // mesh and determine vertices
    // and used vertices
    const Triangulation<dim, spacedim> &tria = mesh.get_triangulation();

    auto vertices = extract_used_vertices(tria, mapping);

    Assert ( tria.get_vertices().size() == marked_vertices.size() || marked_vertices.size() ==0,
             ExcDimensionMismatch(tria.get_vertices().size(), marked_vertices.size()));

    // If p is an element of marked_vertices,
    // and q is that of used_Vertices,
    // the vector marked_vertices does NOT
    // contain unused vertices if p implies q.
    // I.e., if p is true q must be true
    // (if p is false, q could be false or true).
    // p implies q logic is encapsulated in ~p|q.
    Assert( marked_vertices.size()==0
            ||
            std::equal( marked_vertices.begin(),
                        marked_vertices.end(),
                        tria.get_used_vertices().begin(),
                        [](bool p, bool q)
    {
      return !p || q;
    }),
    ExcMessage("marked_vertices should be a subset of used vertices in the triangulation "
               "but marked_vertices contains one or more vertices that are not used vertices!") );

    // Remove from the map unwanted elements.
    if (marked_vertices.size())
      for (auto it = vertices.begin(); it != vertices.end(); )
        {
          if (marked_vertices[it->first] == false)
            {
              vertices.erase(it++);
            }
          else
            {
              ++it;
            }
        }

    return find_closest_vertex(vertices, p);
  }



  template <int dim, template <int, int> class MeshType, int spacedim>
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
    template <int dim, template <int, int> class MeshType, int spacedim>
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

  template <int dim, template <int, int> class MeshType, int spacedim>
#ifndef _MSC_VER
  typename MeshType<dim, spacedim>::active_cell_iterator
#else
  typename dealii::internal::ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim> >::type
#endif
  find_active_cell_around_point (const MeshType<dim,spacedim> &mesh,
                                 const Point<spacedim>        &p,
                                 const std::vector<bool>      &marked_vertices)
  {
    return
      find_active_cell_around_point<dim,MeshType,spacedim>
      (StaticMappingQ1<dim,spacedim>::mapping,
       mesh, p, marked_vertices).first;
  }


  template <int dim, template <int, int> class MeshType, int spacedim>
#ifndef _MSC_VER
  std::pair<typename MeshType<dim, spacedim>::active_cell_iterator, Point<dim> >
#else
  std::pair<typename dealii::internal::ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim> >::type, Point<dim> >
#endif
  find_active_cell_around_point (const Mapping<dim,spacedim>  &mapping,
                                 const MeshType<dim,spacedim> &mesh,
                                 const Point<spacedim>        &p,
                                 const std::vector<bool>      &marked_vertices)
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
                                      find_closest_vertex(mapping, mesh, p, marked_vertices));

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

        // if the user provided a custom mask for vertices,
        // terminate the search without trying to expand the search
        // to all cells of the triangulation, as done below.
        if (marked_vertices.size() > 0)
          cells_searched = n_active_cells;

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



  template <int dim,int spacedim>
  std::vector<std::vector<Tensor<1,spacedim> > >
  vertex_to_cell_centers_directions(const Triangulation<dim,spacedim> &mesh,
                                    const std::vector<std::set<typename Triangulation<dim,spacedim>::active_cell_iterator> > &vertex_to_cells)
  {
    const std::vector<Point<spacedim> > &vertices = mesh.get_vertices();
    const unsigned int n_vertices = vertex_to_cells.size();

    AssertDimension(vertices.size(), n_vertices);


    std::vector<std::vector<Tensor<1,spacedim> > > vertex_to_cell_centers(n_vertices);
    for (unsigned int vertex=0; vertex<n_vertices; ++vertex)
      if (mesh.vertex_used(vertex))
        {
          const unsigned int n_neighbor_cells = vertex_to_cells[vertex].size();
          vertex_to_cell_centers[vertex].resize(n_neighbor_cells);

          typename std::set<typename Triangulation<dim,spacedim>::active_cell_iterator>::iterator it = vertex_to_cells[vertex].begin();
          for (unsigned int cell=0; cell<n_neighbor_cells; ++cell,++it)
            {
              vertex_to_cell_centers[vertex][cell] = (*it)->center() - vertices[vertex];
              vertex_to_cell_centers[vertex][cell] /= vertex_to_cell_centers[vertex][cell].norm();
            }
        }
    return vertex_to_cell_centers;
  }


  namespace
  {
    template <int spacedim>
    bool
    compare_point_association(const unsigned int a,
                              const unsigned int b,
                              const Tensor<1,spacedim> &point_direction,
                              const std::vector<Tensor<1,spacedim> > &center_directions)
    {
      const double scalar_product_a = center_directions[a] * point_direction;
      const double scalar_product_b = center_directions[b] * point_direction;

      // The function is supposed to return if a is before b. We are looking
      // for the alignment of point direction and center direction, therefore
      // return if the scalar product of a is larger.
      return (scalar_product_a > scalar_product_b);
    }
  }

  template <int dim, template <int, int> class MeshType, int spacedim>
#ifndef _MSC_VER
  std::pair<typename MeshType<dim, spacedim>::active_cell_iterator, Point<dim> >
#else
  std::pair<typename dealii::internal::ActiveCellIterator<dim, spacedim, MeshType<dim, spacedim> >::type, Point<dim> >
#endif
  find_active_cell_around_point (const Mapping<dim,spacedim>                                                    &mapping,
                                 const MeshType<dim,spacedim>                                                   &mesh,
                                 const Point<spacedim>                                                          &p,
                                 const std::vector<std::set<typename MeshType<dim,spacedim>::active_cell_iterator > > &vertex_to_cells,
                                 const std::vector<std::vector<Tensor<1,spacedim> > >                            &vertex_to_cell_centers,
                                 const typename MeshType<dim, spacedim>::active_cell_iterator                    &cell_hint ,
                                 const std::vector<bool>                                                         &marked_vertices)
  {
    std::pair<typename MeshType<dim, spacedim>::active_cell_iterator, Point<dim> > cell_and_position;
    // To handle points at the border we keep track of points which are close to the unit cell:
    std::pair<typename MeshType<dim, spacedim>::active_cell_iterator, Point<dim> > cell_and_position_approx;

    bool found_cell = false;
    bool approx_cell = false;

    unsigned int closest_vertex_index = 0;
    Tensor<1,spacedim> vertex_to_point;
    auto current_cell = cell_hint;

    while (found_cell == false)
      {
        // First look at the vertices of the cell cell_hint. If it's an
        // invalid cell, then query for the closest global vertex
        if (current_cell.state() == IteratorState::valid)
          {
            const unsigned int closest_vertex = find_closest_vertex_of_cell<dim,spacedim>(current_cell , p);
            vertex_to_point = p - current_cell ->vertex(closest_vertex);
            closest_vertex_index = current_cell ->vertex_index(closest_vertex);
          }
        else
          {
            closest_vertex_index = GridTools::find_closest_vertex(mesh,p,marked_vertices);
            vertex_to_point = p - mesh.get_vertices()[closest_vertex_index];
          }

        const double vertex_point_norm = vertex_to_point.norm();
        if (vertex_point_norm > 0)
          vertex_to_point /= vertex_point_norm;

        const unsigned int n_neighbor_cells = vertex_to_cells[closest_vertex_index].size();

        // Create a corresponding map of vectors from vertex to cell center
        std::vector<unsigned int> neighbor_permutation(n_neighbor_cells);

        for (unsigned int i=0; i<n_neighbor_cells; ++i)
          neighbor_permutation[i] = i;

        auto comp = [&](const unsigned int a, const unsigned int b) -> bool
        {
          return compare_point_association<spacedim>(a,b,vertex_to_point,vertex_to_cell_centers[closest_vertex_index]);
        };

        std::sort(neighbor_permutation.begin(),
                  neighbor_permutation.end(),
                  comp);
        // It is possible the vertex is close
        // to an edge, thus we add a tolerance
        // setting it initially to 1e-10
        // to keep also the "best" cell
        double best_distance = 1e-10;

        // Search all of the cells adjacent to the closest vertex of the cell hint
        // Most likely we will find the point in them.
        for (unsigned int i=0; i<n_neighbor_cells; ++i)
          {
            try
              {
                auto cell = vertex_to_cells[closest_vertex_index].begin();
                std::advance(cell,neighbor_permutation[i]);
                const Point<dim> p_unit = mapping.transform_real_to_unit_cell(*cell, p);
                if (GeometryInfo<dim>::is_inside_unit_cell(p_unit))
                  {
                    cell_and_position.first = *cell;
                    cell_and_position.second = p_unit;
                    found_cell = true;
                    approx_cell = false;
                    break;
                  }
                // The point is not inside this cell: checking how far outside it is
                // and whether we want to use this cell as a backup if we can't find
                // a cell within which the point lies.
                const double dist = GeometryInfo<dim>::distance_to_unit_cell(p_unit);
                if (dist < best_distance)
                  {
                    best_distance = dist;
                    cell_and_position_approx.first = *cell;
                    cell_and_position_approx.second = p_unit;
                    approx_cell = true;
                  }
              }
            catch (typename Mapping<dim>::ExcTransformationFailed &)
              {}
          }

        if (found_cell == true)
          return cell_and_position;
        else if (approx_cell == true)
          return cell_and_position_approx;

        // The first time around, we check for vertices in the hint_cell. If that
        // does not work, we set the cell iterator to an invalid one, and look
        // for a global vertex close to the point. If that does not work, we are in
        // trouble, and just throw an exception.
        //
        // If we got here, then we did not find the point. If the
        // current_cell.state() here is not IteratorState::valid, it means that
        // the user did not provide a hint_cell, and at the beginning of the
        // while loop we performed an actual global search on the mesh
        // vertices. Not finding the point then means the point is outside the
        // domain.
        AssertThrow(current_cell.state() == IteratorState::valid,
                    ExcPointNotFound<spacedim>(p));

        current_cell = typename MeshType<dim,spacedim>::active_cell_iterator();
      }
    return cell_and_position;
  }



  template <int dim, int spacedim>
  unsigned int
  find_closest_vertex_of_cell(const typename Triangulation<dim,spacedim>::active_cell_iterator &cell,
                              const Point<spacedim> &position)
  {
    double minimum_distance = position.distance_square(cell->vertex(0));
    unsigned int closest_vertex = 0;

    for (unsigned int v=1; v<GeometryInfo<dim>::vertices_per_cell; ++v)
      {
        const double vertex_distance = position.distance_square(cell->vertex(v));
        if (vertex_distance < minimum_distance)
          {
            closest_vertex = v;
            minimum_distance = vertex_distance;
          }
      }
    return closest_vertex;
  }



  namespace internal
  {
    namespace BoundingBoxPredicate
    {
      template < class MeshType >
      std::tuple< BoundingBox < MeshType::space_dimension >, bool >
      compute_cell_predicate_bounding_box
      (const typename MeshType::cell_iterator &parent_cell,
       const std::function<bool (const typename MeshType::active_cell_iterator &)> &predicate)
      {
        bool has_predicate = false; // Start assuming there's no cells with predicate inside
        std::vector< typename MeshType::active_cell_iterator > active_cells;
        if (parent_cell->active())
          active_cells = {parent_cell};
        else
          //Finding all active cells descendants of the current one (or the current one if it is active)
          active_cells = get_active_child_cells < MeshType > (parent_cell);

        const unsigned int spacedim = MeshType::space_dimension;

        // Looking for the first active cell which has the property predicate
        unsigned int i = 0;
        while ( i < active_cells.size() && !predicate(active_cells[i]) )
          ++i;

        // No active cells or no active cells with property
        if ( active_cells.size() == 0 || i == active_cells.size() )
          {
            BoundingBox<spacedim> bbox;
            return std::make_tuple(bbox, has_predicate);
          }

        // The two boundary points defining the boundary box
        Point<spacedim> maxp = active_cells[i]->vertex(0);
        Point<spacedim> minp = active_cells[i]->vertex(0);

        for (; i < active_cells.size() ; ++i)
          if ( predicate(active_cells[i]) )
            for (unsigned int v=0; v<GeometryInfo<spacedim>::vertices_per_cell; ++v)
              for ( unsigned int d=0; d<spacedim; ++d)
                {
                  minp[d] = std::min( minp[d], active_cells[i]->vertex(v)[d]);
                  maxp[d] = std::max( maxp[d], active_cells[i]->vertex(v)[d]);
                }

        has_predicate = true;
        BoundingBox < spacedim > bbox(std::make_pair(minp,maxp));
        return std::make_tuple(bbox, has_predicate);
      }
    }
  }



  template < class MeshType >
  std::vector< BoundingBox<MeshType::space_dimension> >
  compute_mesh_predicate_bounding_box
  (const MeshType                                                              &mesh,
   const std::function<bool (const typename MeshType::active_cell_iterator &)> &predicate,
   const unsigned int                                                           refinement_level,
   const bool                                                                   allow_merge,
   const unsigned int                                                           max_boxes)
  {
    // Algorithm brief description: begin with creating bounding boxes of all cells at
    // refinement_level (and coarser levels if there are active cells) which have the predicate
    // property. These are then merged

    Assert( refinement_level <= mesh.n_levels(),
            ExcMessage ( "Error: refinement level is higher then total levels in the triangulation!") );

    const unsigned int spacedim = MeshType::space_dimension;
    std::vector< BoundingBox < spacedim > > bounding_boxes;

    // Creating a bounding box for all active cell on coarser level

    for (unsigned int i=0; i < refinement_level; ++i)
      for (typename MeshType::cell_iterator cell: mesh.active_cell_iterators_on_level(i))
        {
          bool has_predicate = false;
          BoundingBox < spacedim > bbox;
          std::tie(bbox, has_predicate) =
            internal::BoundingBoxPredicate::compute_cell_predicate_bounding_box <MeshType> (cell, predicate);
          if (has_predicate)
            bounding_boxes.push_back(bbox);
        }

    // Creating a Bounding Box for all cells on the chosen refinement_level
    for (const typename MeshType::cell_iterator &cell: mesh.cell_iterators_on_level(refinement_level))
      {
        bool has_predicate = false;
        BoundingBox < spacedim > bbox;
        std::tie(bbox, has_predicate) =
          internal::BoundingBoxPredicate::compute_cell_predicate_bounding_box <MeshType> (cell, predicate);
        if (has_predicate)
          bounding_boxes.push_back(bbox);
      }

    if ( !allow_merge)
      // If merging is not requested return the created bounding_boxes
      return bounding_boxes;
    else
      {
        // Merging part of the algorithm
        // Part 1: merging neighbors
        // This array stores the indices of arrays we have already merged
        std::vector<unsigned int> merged_boxes_idx;
        bool found_neighbors = true;

        // We merge only neighbors which can be expressed by a single bounding box
        // e.g. in 1d [0,1] and [1,2] can be described with [0,2] without losing anything
        while (found_neighbors)
          {
            found_neighbors = false;
            for (unsigned int i=0; i<bounding_boxes.size()-1; ++i)
              {
                if ( std::find(merged_boxes_idx.begin(),merged_boxes_idx.end(),i) == merged_boxes_idx.end())
                  for (unsigned int j=i+1; j<bounding_boxes.size(); ++j)
                    if ( std::find(merged_boxes_idx.begin(),merged_boxes_idx.end(),j) == merged_boxes_idx.end()
                         && bounding_boxes[i].get_neighbor_type (bounding_boxes[j]) == NeighborType::mergeable_neighbors  )
                      {
                        bounding_boxes[i].merge_with(bounding_boxes[j]);
                        merged_boxes_idx.push_back(j);
                        found_neighbors = true;
                      }
              }
          }

        // Copying the merged boxes into merged_b_boxes
        std::vector< BoundingBox < spacedim > > merged_b_boxes;
        for (unsigned int i=0; i<bounding_boxes.size(); ++i)
          if (std::find(merged_boxes_idx.begin(),merged_boxes_idx.end(),i) == merged_boxes_idx.end())
            merged_b_boxes.push_back(bounding_boxes[i]);

        // Part 2: if there are too many bounding boxes, merging smaller boxes
        // This has sense only in dimension 2 or greater, since  in dimension 1,
        // neighboring intervals can always be merged without problems
        if ( (merged_b_boxes.size() > max_boxes) && (spacedim > 1) )
          {
            std::vector<double> volumes;
            for (unsigned int i=0; i< merged_b_boxes.size(); ++i)
              volumes.push_back(merged_b_boxes[i].volume());

            while ( merged_b_boxes.size() > max_boxes)
              {
                unsigned int min_idx = std::min_element(volumes.begin(),volumes.end()) -
                                       volumes.begin();
                volumes.erase(volumes.begin() + min_idx);
                //Finding a neighbor
                bool not_removed = true;
                for (unsigned int i=0; i<merged_b_boxes.size() && not_removed; ++i)
                  // We merge boxes if we have "attached" or "mergeable" neighbors, even though mergeable should
                  // be dealt with in Part 1
                  if ( i != min_idx &&
                       (merged_b_boxes[i].get_neighbor_type (merged_b_boxes[min_idx]) == NeighborType::attached_neighbors ||
                        merged_b_boxes[i].get_neighbor_type (merged_b_boxes[min_idx]) == NeighborType::mergeable_neighbors) )
                    {
                      merged_b_boxes[i].merge_with(merged_b_boxes[min_idx]);
                      merged_b_boxes.erase(merged_b_boxes.begin() + min_idx);
                      not_removed = false;
                    }
                Assert( !not_removed,
                        ExcMessage ( "Error: couldn't merge bounding boxes!") );
              }
          }
        Assert( merged_b_boxes.size() <= max_boxes,
                ExcMessage ( "Error: couldn't reach target number of bounding boxes!") );
        return merged_b_boxes;
      }
  }




  template <int spacedim>
  std::tuple< std::vector< std::vector< unsigned int > >,
      std::map< unsigned int, unsigned int>,
      std::map< unsigned int, std::vector< unsigned int > > >
      guess_point_owner (const std::vector< std::vector< BoundingBox<spacedim> > >
                         &global_bboxes,
                         const std::vector< Point<spacedim> >    &points)
  {
    unsigned int n_procs = global_bboxes.size();
    std::vector< std::vector< unsigned int > > point_owners(n_procs);
    std::map< unsigned int, unsigned int> map_owners_found;
    std::map< unsigned int, std::vector< unsigned int > > map_owners_guessed;

    unsigned int n_points = points.size();
    for (unsigned int pt=0; pt<n_points; ++pt)
      {
        // Keep track of how many processes we guess to own the point
        std::vector< unsigned int > owners_found;
        // Check in which other processes the point might be
        for (unsigned int rk=0; rk<n_procs; ++rk)
          {
            for (const BoundingBox<spacedim> &bbox: global_bboxes[rk])
              if (bbox.point_inside(points[pt]))
                {
                  point_owners[rk].emplace_back(pt);
                  owners_found.emplace_back(rk);
                  break; // We can check now the next process
                }
          }
        Assert(owners_found.size() > 0,
               ExcMessage("No owners found for the point " + std::to_string(pt)));
        if (owners_found.size()==1)
          map_owners_found[pt] = owners_found[0];
        else
          // Multiple owners
          map_owners_guessed[pt] = owners_found;
      }

    std::tuple< std::vector< std::vector< unsigned int > >,
        std::map< unsigned int, unsigned int>,
        std::map< unsigned int, std::vector< unsigned int > > >
        output_tuple;

    std::get<0>(output_tuple) = point_owners;
    std::get<1>(output_tuple) = map_owners_found;
    std::get<2>(output_tuple) = map_owners_guessed;

    return output_tuple;
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
    // use parallel::distributed::Triangulation in any meaningful
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
    std::map<types::subdomain_id,std::vector<std::tuple<types::global_vertex_index,
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
                        for (; adjacent_cell!=end_adj_cell; ++adjacent_cell)
                          if ((*adjacent_cell)->subdomain_id()!=cell->subdomain_id())
                            {
                              std::pair<types::subdomain_id,types::global_vertex_index>
                              tmp((*adjacent_cell)->subdomain_id(), cell->vertex_index(i));
                              if (vertices_added.find(tmp)==vertices_added.end())
                                {
                                  vertices_to_send[(*adjacent_cell)->subdomain_id()].emplace_back
                                  (i, cell->vertex_index(i), cell->id().to_string());
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
    int ierr = MPI_Allgather(&next_index, 1, DEAL_II_DOF_INDEX_MPI_TYPE, indices.data(),
                             indices.size(), DEAL_II_DOF_INDEX_MPI_TYPE, triangulation.get_communicator());
    AssertThrowMPI(ierr);
    const types::global_vertex_index shift = std::accumulate(indices.begin(),
                                                             indices.begin()+triangulation.locally_owned_subdomain(),
                                                             types::global_vertex_index(0));

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
             std::vector<std::tuple<types::global_vertex_index,
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
            vertices_send_buffers[i][2*j] = std::get<0>(vert_to_send_it->second[j]);
            vertices_send_buffers[i][2*j+1] =
              local_to_global_vertex_index[std::get<1>(vert_to_send_it->second[j])];
          }

        // Send the message
        ierr = MPI_Isend(&vertices_send_buffers[i][0],buffer_size,DEAL_II_VERTEX_INDEX_MPI_TYPE,
                         destination, 0, triangulation.get_communicator(), &first_requests[i]);
        AssertThrowMPI(ierr);
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
        ierr = MPI_Recv(&vertices_recv_buffers[i][0],buffer_size,DEAL_II_VERTEX_INDEX_MPI_TYPE,
                        source, 0, triangulation.get_communicator(), MPI_STATUS_IGNORE);
        AssertThrowMPI(ierr);
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
            std::string cell_id = std::get<2>(vert_to_send_it->second[j]);
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
        ierr = MPI_Isend(&cellids_send_buffers[i][0], buffer_size, MPI_CHAR,
                         destination, 0, triangulation.get_communicator(), &second_requests[i]);
        AssertThrowMPI(ierr);
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
        ierr = MPI_Recv(&cellids_recv_buffers[i][0],buffer_size, MPI_CHAR,
                        source, 0, triangulation.get_communicator(), MPI_STATUS_IGNORE);
        AssertThrowMPI(ierr);
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
  get_vertex_connectivity_of_cells_on_level (const Triangulation<dim,spacedim> &triangulation,
                                             const unsigned int                level,
                                             DynamicSparsityPattern            &cell_connectivity)
  {
    std::vector<std::vector<unsigned int> > vertex_to_cell(triangulation.n_vertices());
    for (typename Triangulation<dim,spacedim>::cell_iterator cell=
           triangulation.begin(level); cell != triangulation.end(level); ++cell)
      {
        for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
          vertex_to_cell[cell->vertex_index(v)].push_back(cell->index());
      }

    cell_connectivity.reinit (triangulation.n_cells(level),
                              triangulation.n_cells(level));
    for (typename Triangulation<dim,spacedim>::cell_iterator cell=
           triangulation.begin(level); cell != triangulation.end(level); ++cell)
      {
        for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
          for (unsigned int n=0; n<vertex_to_cell[cell->vertex_index(v)].size(); ++n)
            cell_connectivity.add(cell->index(), vertex_to_cell[cell->vertex_index(v)][n]);
      }
  }



  template <int dim, int spacedim>
  void
  partition_triangulation (const unsigned int           n_partitions,
                           Triangulation<dim,spacedim> &triangulation,
                           const SparsityTools::Partitioner  partitioner
                          )
  {
    Assert ((dynamic_cast<parallel::distributed::Triangulation<dim,spacedim>*>
             (&triangulation)
             == nullptr),
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
    DynamicSparsityPattern cell_connectivity;
    get_face_connectivity_of_cells (triangulation, cell_connectivity);

    SparsityPattern sp_cell_connectivity;
    sp_cell_connectivity.copy_from(cell_connectivity);
    partition_triangulation (n_partitions,
                             sp_cell_connectivity,
                             triangulation,
                             partitioner
                            );
  }



  template <int dim, int spacedim>
  void
  partition_triangulation (const unsigned int           n_partitions,
                           const SparsityPattern        &cell_connection_graph,
                           Triangulation<dim,spacedim>  &triangulation,
                           const SparsityTools::Partitioner  partitioner
                          )
  {
    Assert ((dynamic_cast<parallel::distributed::Triangulation<dim,spacedim>*>
             (&triangulation)
             == nullptr),
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
    SparsityTools::partition (cell_connection_graph, n_partitions,  partition_indices, partitioner);

    // finally loop over all cells and set the
    // subdomain ids
    for (typename dealii::internal::ActiveCellIterator<dim, spacedim, Triangulation<dim, spacedim> >::type
         cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell)
      cell->set_subdomain_id (partition_indices[cell->active_cell_index()]);
  }


  namespace
  {
    /**
     * recursive helper function for partition_triangulation_zorder
     */
    template <class IT>
    void set_subdomain_id_in_zorder_recursively(IT                 cell,
                                                unsigned int       &current_proc_idx,
                                                unsigned int       &current_cell_idx,
                                                const unsigned int n_active_cells,
                                                const unsigned int n_partitions)
    {
      if (cell->active())
        {
          while (current_cell_idx >= floor((long)n_active_cells*(current_proc_idx+1)/n_partitions))
            ++current_proc_idx;
          cell->set_subdomain_id(current_proc_idx);
          ++current_cell_idx;
        }
      else
        {
          for (unsigned int n=0; n<cell->n_children(); ++n)
            set_subdomain_id_in_zorder_recursively(cell->child(n),
                                                   current_proc_idx,
                                                   current_cell_idx,
                                                   n_active_cells,
                                                   n_partitions);
        }
    }
  }

  template <int dim, int spacedim>
  void
  partition_triangulation_zorder (const unsigned int          n_partitions,
                                  Triangulation<dim,spacedim> &triangulation)
  {
    Assert ((dynamic_cast<parallel::distributed::Triangulation<dim,spacedim>*>
             (&triangulation)
             == nullptr),
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

    // Duplicate the coarse cell reordoring
    // as done in p4est
    std::vector<types::global_dof_index> coarse_cell_to_p4est_tree_permutation;
    std::vector<types::global_dof_index> p4est_tree_to_coarse_cell_permutation;

    DynamicSparsityPattern cell_connectivity;
    GridTools::get_vertex_connectivity_of_cells_on_level (triangulation, 0, cell_connectivity);
    coarse_cell_to_p4est_tree_permutation.resize (triangulation.n_cells(0));
    SparsityTools::reorder_hierarchical (cell_connectivity,
                                         coarse_cell_to_p4est_tree_permutation);

    p4est_tree_to_coarse_cell_permutation
      = Utilities::invert_permutation (coarse_cell_to_p4est_tree_permutation);

    unsigned int current_proc_idx=0;
    unsigned int current_cell_idx=0;
    const unsigned int n_active_cells = triangulation.n_active_cells();

    // set subdomain id for active cell descendants
    // of each coarse cell in permuted order
    for (unsigned int idx=0; idx<triangulation.n_cells(0); ++idx)
      {
        const unsigned int coarse_cell_idx = p4est_tree_to_coarse_cell_permutation[idx];
        typename Triangulation<dim,spacedim>::cell_iterator
        coarse_cell (&triangulation, 0, coarse_cell_idx);

        set_subdomain_id_in_zorder_recursively(coarse_cell,
                                               current_proc_idx,
                                               current_cell_idx,
                                               n_active_cells,
                                               n_partitions);
      }

    // if all children of a cell are active (e.g. we
    // have a cell that is refined once and no part
    // is refined further), p4est places all of them
    // on the same processor. The new owner will be
    // the processor with the largest number of children
    // (ties are broken by picking the lower rank).
    // Duplicate this logic here.
    {
      typename Triangulation<dim,spacedim>::cell_iterator
      cell = triangulation.begin(),
      endc = triangulation.end();
      for (; cell!=endc; ++cell)
        {
          if (cell->active())
            continue;
          bool all_children_active = true;
          std::map<unsigned int, unsigned int> map_cpu_n_cells;
          for (unsigned int n=0; n<cell->n_children(); ++n)
            if (!cell->child(n)->active())
              {
                all_children_active = false;
                break;
              }
            else
              ++map_cpu_n_cells[cell->child(n)->subdomain_id()];

          if (!all_children_active)
            continue;

          unsigned int new_owner = cell->child(0)->subdomain_id();
          for (std::map<unsigned int, unsigned int>::iterator it = map_cpu_n_cells.begin();
               it != map_cpu_n_cells.end();
               ++it)
            if (it->second > map_cpu_n_cells[new_owner])
              new_owner = it->first;

          for (unsigned int n=0; n<cell->n_children(); ++n)
            cell->child(n)->set_subdomain_id(new_owner);
        }
    }
  }


  template <int dim, int spacedim>
  void
  partition_multigrid_levels (Triangulation<dim,spacedim> &triangulation)
  {
    unsigned int n_levels = triangulation.n_levels();
    for (int lvl = n_levels-1; lvl>=0; --lvl)
      {
        typename Triangulation<dim,spacedim>::cell_iterator
        cell = triangulation.begin(lvl),
        endc = triangulation.end(lvl);
        for (; cell!=endc; ++cell)
          {
            if (!cell->has_children())
              cell->set_level_subdomain_id(cell->subdomain_id());
            else
              {
                Assert(cell->child(0)->level_subdomain_id()
                       != numbers::artificial_subdomain_id, ExcInternalError());
                cell->set_level_subdomain_id(cell->child(0)->level_subdomain_id());
              }
          }
      }
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



  template <int dim, int spacedim>
  double
  minimal_cell_diameter (const Triangulation<dim, spacedim> &triangulation)
  {
    double min_diameter = std::numeric_limits<double>::max();
    for (const auto &cell: triangulation.active_cell_iterators())
      if (!cell->is_artificial())
        min_diameter = std::min (min_diameter,
                                 cell->diameter());

    double global_min_diameter = 0;

#ifdef DEAL_II_WITH_MPI
    if (const parallel::Triangulation<dim,spacedim> *p_tria
        = dynamic_cast<const parallel::Triangulation<dim,spacedim>*>(&triangulation))
      global_min_diameter = Utilities::MPI::min (min_diameter, p_tria->get_communicator());
    else
#endif
      global_min_diameter = min_diameter;

    return global_min_diameter;
  }



  template <int dim, int spacedim>
  double
  maximal_cell_diameter (const Triangulation<dim, spacedim> &triangulation)
  {
    double max_diameter = 0.;
    for (const auto &cell: triangulation.active_cell_iterators())
      if (!cell->is_artificial())
        max_diameter = std::max (max_diameter,
                                 cell->diameter());

    double global_max_diameter = 0;

#ifdef DEAL_II_WITH_MPI
    if (const parallel::Triangulation<dim,spacedim> *p_tria
        = dynamic_cast<const parallel::Triangulation<dim,spacedim>*>(&triangulation))
      global_max_diameter = Utilities::MPI::max (max_diameter, p_tria->get_communicator());
    else
#endif
      global_max_diameter = max_diameter;

    return global_max_diameter;
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
          = std::accumulate (parent_alternating_forms,
                             parent_alternating_forms + GeometryInfo<structdim>::vertices_per_cell,
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
                         std::integral_constant<int, 1>)
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
                         std::integral_constant<int, 2>)
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
                         std::integral_constant<int, 3>)
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
                                  std::integral_constant<int, structdim>())
                                 .distance (get_face_midpoint
                                            (object,
                                             e,
                                             std::integral_constant<int, structdim>())));

        return diameter;
      }



      /**
       * Try to fix up a single cell by moving around its midpoint. Return whether we succeeded with this.
       */
      template <typename Iterator>
      bool
      fix_up_object (const Iterator &object)
      {
        const unsigned int structdim = Iterator::AccessorType::structure_dimension;
        const unsigned int spacedim  = Iterator::AccessorType::space_dimension;

        // right now we can only deal with cells that have been refined
        // isotropically because that is the only case where we have a cell
        // mid-point that can be moved around without having to consider
        // boundary information
        Assert (object->has_children(), ExcInternalError());
        Assert (object->refinement_case() == RefinementCase<structdim>::isotropic_refinement,
                ExcNotImplemented());

        // get the current location of the object mid-vertex:
        Point<spacedim> object_mid_point
          = object->child(0)->vertex (GeometryInfo<structdim>::max_children_per_cell-1);

        // now do a few steepest descent steps to reduce the objective
        // function. compute the diameter in the helper function above
        unsigned int iteration = 0;
        const double diameter = minimal_diameter (object);

        // current value of objective function and initial delta
        double current_value = objective_function (object, object_mid_point);
        double initial_delta = 0;

        do
          {
            // choose a step length that is initially 1/4 of the child
            // objects' diameter, and a sequence whose sum does not converge
            // (to avoid premature termination of the iteration)
            const double step_length = diameter / 4 / (iteration + 1);

            // compute the objective function's derivative using a two-sided
            // difference formula with eps=step_length/10
            Tensor<1,spacedim> gradient;
            for (unsigned int d=0; d<spacedim; ++d)
              {
                const double eps = step_length/10;

                Tensor<1,spacedim> h;
                h[d] = eps/2;

                gradient[d] = (objective_function (object,
                                                   project_to_object(object, object_mid_point + h))
                               -
                               objective_function (object,
                                                   project_to_object(object, object_mid_point - h)))
                              /eps;
              }

            // there is nowhere to go
            if (gradient.norm() == 0)
              break;

            // We need to go in direction -gradient. the optimal value of the
            // objective function is zero, so assuming that the model is
            // quadratic we would have to go -2*val/||gradient|| in this
            // direction, make sure we go at most step_length into this
            // direction
            object_mid_point -= std::min(2 * current_value / (gradient*gradient),
                                         step_length / gradient.norm()) * gradient;
            object_mid_point = project_to_object(object, object_mid_point);

            // compute current value of the objective function
            const double previous_value = current_value;
            current_value = objective_function (object, object_mid_point);

            if (iteration == 0)
              initial_delta = (previous_value - current_value);

            // stop if we aren't moving much any more
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
                         std::integral_constant<int, 1>,
                         std::integral_constant<int, 1>)
      {
        // nothing to do for the faces of cells in 1d
      }



      // possibly fix up the faces of a cell by moving around its mid-points
      template <int dim, int spacedim>
      void fix_up_faces (const typename dealii::Triangulation<dim,spacedim>::cell_iterator &cell,
                         std::integral_constant<int, dim>,
                         std::integral_constant<int, spacedim>)
      {
        // see if we first can fix up some of the faces of this object. We can
        // mess with faces if and only if the neighboring cell is not even
        // more refined than we are (since in that case the sub-faces have
        // themselves children that we can't move around any more). however,
        // the latter case shouldn't happen anyway: if the current face is
        // distorted but the neighbor is even more refined, then the face had
        // been deformed before already, and had been ignored at the time; we
        // should then also be able to ignore it this time as well
        for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
          {
            Assert (cell->face(f)->has_children(), ExcInternalError());
            Assert (cell->face(f)->refinement_case() ==
                    RefinementCase<dim - 1>::isotropic_refinement,
                    ExcInternalError());

            bool subface_is_more_refined = false;
            for (unsigned int g=0; g<GeometryInfo<dim>::max_children_per_face; ++g)
              if (cell->face(f)->child(g)->has_children())
                {
                  subface_is_more_refined = true;
                  break;
                }

            if (subface_is_more_refined == true)
              continue;

            // we finally know that we can do something about this face
            fix_up_object (cell->face(f));
          }
      }
    } /* namespace FixUpDistortedChildCells */
  } /* namespace internal */


  template <int dim, int spacedim>
  typename Triangulation<dim,spacedim>::DistortedCellList
  fix_up_distorted_child_cells
  (const typename Triangulation<dim,spacedim>::DistortedCellList &distorted_cells,
   Triangulation<dim,spacedim> &/*triangulation*/)
  {
    typename Triangulation<dim,spacedim>::DistortedCellList unfixable_subset;

    // loop over all cells that we have to fix up
    for (typename std::list<typename Triangulation<dim,spacedim>::cell_iterator>::const_iterator
         cell_ptr = distorted_cells.distorted_cells.begin();
         cell_ptr != distorted_cells.distorted_cells.end(); ++cell_ptr)
      {
        const typename Triangulation<dim,spacedim>::cell_iterator
        cell = *cell_ptr;

        Assert(!cell->active(),
               ExcMessage("This function is only valid for a list of cells that "
                          "have children (i.e., no cell in the list may be active)."));

        internal::FixUpDistortedChildCells
        ::fix_up_faces (cell,
                        std::integral_constant<int, dim>(),
                        std::integral_constant<int, spacedim>());

        // If possible, fix up the object.
        if (!internal::FixUpDistortedChildCells::fix_up_object (cell))
          unfixable_subset.distorted_cells.push_back (cell);
      }

    return unfixable_subset;
  }



  template <int dim, int spacedim>
  void copy_boundary_to_manifold_id(Triangulation<dim, spacedim> &tria,
                                    const bool reset_boundary_ids)
  {
    const auto src_boundary_ids = tria.get_boundary_ids();
    std::vector<types::manifold_id> dst_manifold_ids(src_boundary_ids.size());
    auto m_it = dst_manifold_ids.begin();
    for (auto b : src_boundary_ids)
      {
        *m_it = static_cast<types::manifold_id>(b);
        ++m_it;
      }
    const std::vector<types::boundary_id> reset_boundary_id =
      reset_boundary_ids ?
      std::vector<types::boundary_id>(src_boundary_ids.size(), 0) : src_boundary_ids;
    map_boundary_to_manifold_ids(src_boundary_ids, dst_manifold_ids, tria, reset_boundary_id);
  }



  template <int dim, int spacedim>
  void map_boundary_to_manifold_ids(const std::vector<types::boundary_id> &src_boundary_ids,
                                    const std::vector<types::manifold_id> &dst_manifold_ids,
                                    Triangulation<dim, spacedim> &tria,
                                    const std::vector<types::boundary_id> &reset_boundary_ids_)
  {
    AssertDimension(src_boundary_ids.size(), dst_manifold_ids.size());
    const auto reset_boundary_ids = reset_boundary_ids_.size() ?
                                    reset_boundary_ids_ : src_boundary_ids;
    AssertDimension(reset_boundary_ids.size(), src_boundary_ids.size());

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
              {
                auto bid = cell->face(f)->line(e)->boundary_id();
                auto ind = std::find(src_boundary_ids.begin(), src_boundary_ids.end(), bid)-
                           src_boundary_ids.begin();
                if ((unsigned int)ind < src_boundary_ids.size())
                  cell->face(f)->line(e)->set_manifold_id(dst_manifold_ids[ind]);
              }

    // now do cells
    for (typename Triangulation<dim,spacedim>::active_cell_iterator
         cell=tria.begin_active();
         cell != tria.end(); ++cell)
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        if (cell->face(f)->at_boundary())
          {
            auto bid = cell->face(f)->boundary_id();
            auto ind = std::find(src_boundary_ids.begin(), src_boundary_ids.end(), bid)-
                       src_boundary_ids.begin();

            if ((unsigned int)ind < src_boundary_ids.size())
              {
                // assign the manifold id
                cell->face(f)->set_manifold_id(dst_manifold_ids[ind]);
                // then reset boundary id
                cell->face(f)->set_boundary_id(reset_boundary_ids[ind]);
              }

            if (dim >= 3)
              for (signed int e=0; e<static_cast<signed int>(GeometryInfo<dim>::lines_per_face); ++e)
                {
                  auto bid = cell->face(f)->line(e)->boundary_id();
                  auto ind = std::find(src_boundary_ids.begin(), src_boundary_ids.end(), bid)-
                             src_boundary_ids.begin();
                  if ((unsigned int)ind < src_boundary_ids.size())
                    cell->face(f)->line(e)->set_boundary_id(reset_boundary_ids[ind]);
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

  template <int dim, int spacedim>
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


  template <int dim, int spacedim>
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

  template <int dim, int spacedim>
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


  template <int dim, int spacedim>
  void regularize_corner_cells (Triangulation<dim,spacedim> &tria,
                                const double limit_angle_fraction)
  {
    if (dim == 1)
      return; // Nothing to do

    // Check that we don't have hanging nodes
    AssertThrow(!tria.has_hanging_nodes(), ExcMessage("The input Triangulation cannot "
                                                      "have hanging nodes."));


    bool has_cells_with_more_than_dim_faces_on_boundary = true;
    bool has_cells_with_dim_faces_on_boundary = false;

    unsigned int refinement_cycles = 0;

    while (has_cells_with_more_than_dim_faces_on_boundary)
      {
        has_cells_with_more_than_dim_faces_on_boundary = false;

        for (auto cell: tria.active_cell_iterators())
          {
            unsigned int boundary_face_counter = 0;
            for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
              if (cell->face(f)->at_boundary())
                boundary_face_counter++;
            if (boundary_face_counter > dim)
              {
                has_cells_with_more_than_dim_faces_on_boundary = true;
                break;
              }
            else if (boundary_face_counter == dim)
              has_cells_with_dim_faces_on_boundary = true;
          }
        if (has_cells_with_more_than_dim_faces_on_boundary)
          {
            tria.refine_global(1);
            refinement_cycles++;
          }
      }

    if (has_cells_with_dim_faces_on_boundary)
      {
        tria.refine_global(1);
        refinement_cycles++;
      }
    else
      {
        while (refinement_cycles>0)
          {
            for (auto cell: tria.active_cell_iterators())
              cell->set_coarsen_flag();
            tria.execute_coarsening_and_refinement();
            refinement_cycles--;
          }
        return;
      }

    std::vector<bool> cells_to_remove(tria.n_active_cells(), false);
    std::vector<Point<spacedim> > vertices = tria.get_vertices();

    std::vector<bool> faces_to_remove(tria.n_raw_faces(),false);

    std::vector<CellData<dim> > cells_to_add;
    SubCellData                 subcelldata_to_add;

    // Trick compiler for dimension independent things
    const unsigned int
    v0 = 0, v1 = 1,
    v2 = (dim > 1 ? 2:0), v3 = (dim > 1 ? 3:0);

    for (auto cell : tria.active_cell_iterators())
      {
        double angle_fraction = 0;
        unsigned int vertex_at_corner = numbers::invalid_unsigned_int;

        if (dim == 2)
          {
            Tensor<1,spacedim> p0;
            p0[spacedim > 1 ? 1 : 0] = 1;
            Tensor<1,spacedim> p1;
            p1[0] = 1;

            if (cell->face(v0)->at_boundary() && cell->face(v3)->at_boundary())
              {
                p0 = cell->vertex(v0) -  cell->vertex(v2);
                p1 = cell->vertex(v3) -  cell->vertex(v2);
                vertex_at_corner = v2;
              }
            else  if (cell->face(v3)->at_boundary() && cell->face(v1)->at_boundary())
              {
                p0 = cell->vertex(v2) -  cell->vertex(v3);
                p1 = cell->vertex(v1) -  cell->vertex(v3);
                vertex_at_corner = v3;
              }
            else  if (cell->face(1)->at_boundary() && cell->face(2)->at_boundary())
              {
                p0 = cell->vertex(v0) -  cell->vertex(v1);
                p1 = cell->vertex(v3) -  cell->vertex(v1);
                vertex_at_corner = v1;
              }
            else  if (cell->face(2)->at_boundary() && cell->face(0)->at_boundary())
              {
                p0 = cell->vertex(v2) -  cell->vertex(v0);
                p1 = cell->vertex(v1) -  cell->vertex(v0);
                vertex_at_corner = v0;
              }
            p0 /= p0.norm();
            p1 /= p1.norm();
            angle_fraction = std::acos(p0*p1)/numbers::PI;

          }
        else
          {
            Assert(false, ExcNotImplemented());
          }

        if (angle_fraction > limit_angle_fraction)
          {

            auto flags_removal = [&](unsigned int f1, unsigned int f2,
                                     unsigned int n1, unsigned int n2) -> void
            {
              cells_to_remove[cell->active_cell_index()] = true;
              cells_to_remove[cell->neighbor(n1)->active_cell_index()] = true;
              cells_to_remove[cell->neighbor(n2)->active_cell_index()] = true;

              faces_to_remove[cell->face(f1)->index()] = true;
              faces_to_remove[cell->face(f2)->index()] = true;

              faces_to_remove[cell->neighbor(n1)->face(f1)->index()] = true;
              faces_to_remove[cell->neighbor(n2)->face(f2)->index()] = true;
            };

            auto cell_creation = [&](
                                   const unsigned int vv0,
                                   const unsigned int vv1,
                                   const unsigned int f0,
                                   const unsigned int f1,

                                   const unsigned int n0,
                                   const unsigned int v0n0,
                                   const unsigned int v1n0,

                                   const unsigned int n1,
                                   const unsigned int v0n1,
                                   const unsigned int v1n1)
            {
              CellData<dim> c1, c2;
              CellData<1> l1, l2;

              c1.vertices[v0] = cell->vertex_index(vv0);
              c1.vertices[v1] = cell->vertex_index(vv1);
              c1.vertices[v2] = cell->neighbor(n0)->vertex_index(v0n0);
              c1.vertices[v3] = cell->neighbor(n0)->vertex_index(v1n0);

              c1.manifold_id = cell->manifold_id();
              c1.material_id = cell->material_id();

              c2.vertices[v0] = cell->vertex_index(vv0);
              c2.vertices[v1] = cell->neighbor(n1)->vertex_index(v0n1);
              c2.vertices[v2] = cell->vertex_index(vv1);
              c2.vertices[v3] = cell->neighbor(n1)->vertex_index(v1n1);

              c2.manifold_id = cell->manifold_id();
              c2.material_id = cell->material_id();

              l1.vertices[0] = cell->vertex_index(vv0);
              l1.vertices[1] = cell->neighbor(n0)->vertex_index(v0n0);

              l1.boundary_id = cell->line(f0)->boundary_id();
              l1.manifold_id = cell->line(f0)->manifold_id();
              subcelldata_to_add.boundary_lines.push_back(l1);

              l2.vertices[0] = cell->vertex_index(vv0);
              l2.vertices[1] = cell->neighbor(n1)->vertex_index(v0n1);

              l2.boundary_id = cell->line(f1)->boundary_id();
              l2.manifold_id = cell->line(f1)->manifold_id();
              subcelldata_to_add.boundary_lines.push_back(l2);

              cells_to_add.push_back(c1);
              cells_to_add.push_back(c2);
            };

            if (dim == 2)
              {
                switch (vertex_at_corner)
                  {
                  case 0:
                    flags_removal(0,2,3,1);
                    cell_creation(0,3, 0,2, 3,2,3,  1,1,3);
                    break;
                  case 1:
                    flags_removal(1,2,3,0);
                    cell_creation(1,2, 2,1, 0,0,2, 3,3,2);
                    break;
                  case 2:
                    flags_removal(3,0,1,2);
                    cell_creation(2,1, 3,0, 1,3,1, 2,0,1);
                    break;
                  case 3:
                    flags_removal(3,1,0,2);
                    cell_creation(3,0, 1,3, 2,1,0, 0,2,0);
                    break;
                  }
              }
            else
              {
                Assert(false, ExcNotImplemented());
              }
          }
      }

    // if no cells need to be added, then no regularization is necessary. Restore things
    // as they were before this function was called.
    if (cells_to_add.size() == 0)
      {
        while (refinement_cycles>0)
          {
            for (auto cell: tria.active_cell_iterators())
              cell->set_coarsen_flag();
            tria.execute_coarsening_and_refinement();
            refinement_cycles--;
          }
        return;
      }

    // add the cells that were not marked as skipped
    for (auto cell : tria.active_cell_iterators())
      {
        if (cells_to_remove[cell->active_cell_index()] == false)
          {
            CellData<dim> c;
            for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
              c.vertices[v] = cell->vertex_index(v);
            c.manifold_id = cell->manifold_id();
            c.material_id = cell->material_id();
            cells_to_add.push_back(c);
          }
      }

    // Face counter for both dim == 2 and dim == 3
    typename Triangulation<dim,spacedim>::active_face_iterator
    face = tria.begin_active_face(),
    endf = tria.end_face();
    for (; face != endf; ++face)
      if ( (face->at_boundary() || face->manifold_id() != numbers::invalid_manifold_id)
           && faces_to_remove[face->index()] == false)
        {
          for (unsigned int l=0; l<GeometryInfo<dim>::lines_per_face; ++l)
            {
              CellData<1> line;
              if (dim == 2)
                {
                  for (unsigned int v=0; v<GeometryInfo<1>::vertices_per_cell; ++v)
                    line.vertices[v] = face->vertex_index(v);
                  line.boundary_id = face->boundary_id();
                  line.manifold_id = face->manifold_id();
                }
              else
                {
                  for (unsigned int v=0; v<GeometryInfo<1>::vertices_per_cell; ++v)
                    line.vertices[v] = face->line(l)->vertex_index(v);
                  line.boundary_id = face->line(l)->boundary_id();
                  line.manifold_id = face->line(l)->manifold_id();
                }
              subcelldata_to_add.boundary_lines.push_back(line);
            }
          if (dim == 3)
            {
              CellData<2> quad;
              for (unsigned int v=0; v<GeometryInfo<2>::vertices_per_cell; ++v)
                quad.vertices[v] = face->vertex_index(v);
              quad.boundary_id = face->boundary_id();
              quad.manifold_id = face->manifold_id();
              subcelldata_to_add.boundary_quads.push_back(quad);
            }
        }
    GridTools::delete_unused_vertices(vertices, cells_to_add, subcelldata_to_add);
    GridReordering<dim,spacedim>::reorder_cells(cells_to_add, true);

    // Save manifolds
    auto manifold_ids = tria.get_manifold_ids();
    std::map<types::manifold_id, std::unique_ptr<Manifold<dim,spacedim> > > manifolds;
    // Set manifolds in new Triangulation
    for (auto manifold_id: manifold_ids)
      if (manifold_id != numbers::invalid_manifold_id)
        manifolds[manifold_id] = tria.get_manifold(manifold_id).clone();

    tria.clear();

    tria.create_triangulation(vertices, cells_to_add, subcelldata_to_add);

    // Restore manifolds
    for (auto manifold_id: manifold_ids)
      if (manifold_id != numbers::invalid_manifold_id)
        tria.set_manifold(manifold_id, *manifolds[manifold_id]);
  }



  template <int dim, int spacedim>
  std::tuple<
  std::vector<typename Triangulation<dim, spacedim>::active_cell_iterator >,
      std::vector< std::vector< Point<dim> > >,
      std::vector< std::vector<unsigned int> > >
      compute_point_locations(const Cache<dim,spacedim>                 &cache,
                              const std::vector<Point<spacedim> >       &points,
                              const typename Triangulation<dim, spacedim>::active_cell_iterator &cell_hint)
  {
    // How many points are here?
    const unsigned int np = points.size();

    std::tuple<
    std::vector<typename Triangulation<dim, spacedim>::active_cell_iterator >,
        std::vector< std::vector< Point<dim> > >,
        std::vector< std::vector<unsigned int> > >
        cell_qpoint_map;

    // Now the easy case.
    if (np==0) return cell_qpoint_map;

    // We begin by finding the cell/transform of the first point
    std::pair<typename Triangulation<dim, spacedim>::active_cell_iterator, Point<dim> >
    my_pair;
    if (cell_hint.state() == IteratorState::valid)
      my_pair = GridTools::find_active_cell_around_point
                (cache, points[0],cell_hint);
    else
      my_pair = GridTools::find_active_cell_around_point
                (cache, points[0]);

    std::get<0>(cell_qpoint_map).emplace_back(my_pair.first);
    std::get<1>(cell_qpoint_map).emplace_back(1, my_pair.second);
    std::get<2>(cell_qpoint_map).emplace_back(1, 0);

    // Now the second easy case.
    if (np==1) return cell_qpoint_map;
    // Computing the cell center and diameter
    Point<spacedim> cell_center = std::get<0>(cell_qpoint_map)[0]->center();
    double cell_diameter = std::get<0>(cell_qpoint_map)[0]->diameter()*
                           (0.5 + std::numeric_limits<double>::epsilon() );

    // Cycle over all points left
    for (unsigned int p=1; p< np; ++p)
      {
        // Checking if the point is close to the cell center, in which
        // case calling find active cell with a cell hint
        if ( cell_center.distance(points[p]) < cell_diameter )
          my_pair  = GridTools::find_active_cell_around_point
                     (cache, points[p],std::get<0>(cell_qpoint_map).back());
        else
          my_pair  = GridTools::find_active_cell_around_point
                     (cache, points[p]);

        // Assuming the cell is probably the last cell added
        if ( my_pair.first == std::get<0>(cell_qpoint_map).back() )
          {
            // Found in the last cell: adding the data
            std::get<1>(cell_qpoint_map).back().emplace_back(my_pair.second);
            std::get<2>(cell_qpoint_map).back().emplace_back(p);
          }
        else
          {
            // Check if it is in another cell already found
            typename std::vector<typename Triangulation<dim, spacedim>::active_cell_iterator>::iterator
            cells_it = std::find(std::get<0>(cell_qpoint_map).begin(),std::get<0>(cell_qpoint_map).end()-1,my_pair.first);

            if ( cells_it == std::get<0>(cell_qpoint_map).end()-1 )
              {
                // Cell not found: adding a new cell
                std::get<0>(cell_qpoint_map).emplace_back(my_pair.first);
                std::get<1>(cell_qpoint_map).emplace_back(1, my_pair.second);
                std::get<2>(cell_qpoint_map).emplace_back(1, p);
                // Updating center and radius of the cell
                cell_center = std::get<0>(cell_qpoint_map).back()->center();
                cell_diameter = std::get<0>(cell_qpoint_map).back()->diameter()*
                                (0.5 + std::numeric_limits<double>::epsilon() );
              }
            else
              {
                unsigned int current_cell = cells_it - std::get<0>(cell_qpoint_map).begin();
                // Cell found: just adding the point index and qpoint to the list
                std::get<1>(cell_qpoint_map)[current_cell].emplace_back(my_pair.second);
                std::get<2>(cell_qpoint_map)[current_cell].emplace_back(p);
              }
          }
      }

    // Debug Checking
    Assert(std::get<0>(cell_qpoint_map).size() == std::get<2>(cell_qpoint_map).size(),
           ExcDimensionMismatch(std::get<0>(cell_qpoint_map).size(), std::get<2>(cell_qpoint_map).size()));

    Assert(std::get<0>(cell_qpoint_map).size() == std::get<1>(cell_qpoint_map).size(),
           ExcDimensionMismatch(std::get<0>(cell_qpoint_map).size(), std::get<1>(cell_qpoint_map).size()));

#ifdef DEBUG
    unsigned int c = std::get<0>(cell_qpoint_map).size();
    unsigned int qps = 0;
    // The number of points in all
    // the cells must be the same as
    // the number of points we
    // started off from.
    for (unsigned int n=0; n<c; ++n)
      {
        Assert(std::get<1>(cell_qpoint_map)[n].size() ==
               std::get<2>(cell_qpoint_map)[n].size(),
               ExcDimensionMismatch(std::get<1>(cell_qpoint_map)[n].size(),
                                    std::get<2>(cell_qpoint_map)[n].size()));
        qps += std::get<1>(cell_qpoint_map)[n].size();
      }
    Assert(qps == np,
           ExcDimensionMismatch(qps, np));
#endif

    return cell_qpoint_map;
  }



  namespace internal
  {
    // Functions are needed for distributed compute point locations
    namespace distributed_cptloc
    {
      // Hash function for cells; needed for unordered maps/multimaps
      template < int dim, int spacedim>
      struct cell_hash
      {
        std::size_t operator()(const typename Triangulation<dim, spacedim>::active_cell_iterator &k) const
        {
          // Return active cell index, which is faster than CellId to compute
          return k->active_cell_index();
        }
      };



      // Compute point locations; internal version which returns an unordered map
      // The algorithm is the same as GridTools::compute_point_locations
      template <int dim, int spacedim>
      std::unordered_map< typename Triangulation<dim, spacedim>::active_cell_iterator,
          std::pair<std::vector<Point<dim> >,std::vector<unsigned int> >, cell_hash<dim,spacedim> >
          compute_point_locations_unmap(const GridTools::Cache<dim,spacedim>     &cache,
                                        const std::vector<Point<spacedim> >      &points)
      {
        // How many points are here?
        const unsigned int np = points.size();
        // Creating the output tuple
        std::unordered_map< typename Triangulation<dim, spacedim>::active_cell_iterator,
            std::pair<std::vector<Point<dim> >,std::vector<unsigned int> >, cell_hash<dim,spacedim> >
            cell_qpoint_map;

        // Now the easy case.
        if (np==0) return cell_qpoint_map;
        // We begin by finding the cell/transform of the first point
        auto my_pair  = GridTools::find_active_cell_around_point
                        (cache, points[0]);

        auto last_cell = cell_qpoint_map.emplace(
                           std::make_pair(my_pair.first, std::make_pair(
                                            std::vector<Point<dim> > {my_pair.second},
                                            std::vector<unsigned int> {0})));
        // Now the second easy case.
        if (np==1) return cell_qpoint_map;
        // Computing the cell center and diameter
        Point<spacedim> cell_center = my_pair.first->center();
        double cell_diameter = my_pair.first->diameter()*
                               (0.5 + std::numeric_limits<double>::epsilon() );

        // Cycle over all points left
        for (unsigned int p=1; p< np; ++p)
          {
            // Checking if the point is close to the cell center, in which
            // case calling find active cell with a cell hint
            if ( cell_center.distance(points[p]) < cell_diameter )
              my_pair  = GridTools::find_active_cell_around_point
                         (cache, points[p],last_cell.first->first);
            else
              my_pair  = GridTools::find_active_cell_around_point
                         (cache, points[p]);

            if ( last_cell.first->first == my_pair.first)
              {
                last_cell.first->second.first.emplace_back(my_pair.second);
                last_cell.first->second.second.emplace_back(p);
              }
            else
              {
                // Check if it is in another cell already found
                last_cell = cell_qpoint_map.emplace(std::make_pair(my_pair.first, std::make_pair(
                                                                     std::vector<Point<dim> > {my_pair.second},
                                                                     std::vector<unsigned int> {p})));

                if ( last_cell.second == false )
                  {
                    // Cell already present: adding the new point
                    last_cell.first->second.first.emplace_back(my_pair.second);
                    last_cell.first->second.second.emplace_back(p);
                  }
                else
                  {
                    // New cell was added, updating center and diameter
                    cell_center = my_pair.first->center();
                    cell_diameter = my_pair.first->diameter()*
                                    (0.5 + std::numeric_limits<double>::epsilon() );
                  }
              }
          }

#ifdef DEBUG
        unsigned int qps = 0;
        // The number of points in all
        // the cells must be the same as
        // the number of points we
        // started off from.
        for (const auto &m: cell_qpoint_map)
          {
            Assert(m.second.second.size() ==
                   m.second.first.size(),
                   ExcDimensionMismatch(m.second.second.size(),
                                        m.second.first.size()));
            qps += m.second.second.size();
          }
        Assert(qps == np,
               ExcDimensionMismatch(qps, np));
#endif
        return cell_qpoint_map;
      }



      // Merging the output means to add data to a previous output, here contained
      // in output unmap:
      // if the cell is already present: add information about new points
      // if the cell is not present: add the cell with all information
      //
      // Notice we call "information" the data associated with a point of the sort:
      // cell containing it, transformed point on reference cell, index,
      // rank of the owner etc.
      template <int dim, int spacedim>
      void
      merge_cptloc_outputs(
        std::unordered_map< typename Triangulation<dim, spacedim>::active_cell_iterator,
        std::tuple<
        std::vector< Point<dim> >,
        std::vector< unsigned int >,
        std::vector< Point<spacedim> >,
        std::vector< unsigned int >
        >,
        cell_hash<dim,spacedim>>                                                        &output_unmap,
        const std::vector<typename Triangulation<dim, spacedim>::active_cell_iterator > &in_cells,
        const std::vector< std::vector< Point<dim> > >                                  &in_qpoints,
        const std::vector< std::vector<unsigned int> >                                  &in_maps,
        const std::vector< std::vector< Point<spacedim> > >                             &in_points,
        const unsigned int                                                               in_rank
      )
      {
        // Adding cells, one by one
        for (unsigned int c=0; c< in_cells.size(); ++c)
          {
            // Attempt to add a new cell with its relative data
            auto current_c = output_unmap.emplace(
                               std::make_pair(in_cells[c],
                                              std::make_tuple(in_qpoints[c],
                                                              in_maps[c],
                                                              in_points[c],
                                                              std::vector<unsigned int>
                                                              (in_points[c].size(),in_rank))));
            // If the flag is false no new cell was added:
            if ( current_c.second == false )
              {
                // Cell in output map at current_c.first:
                // Adding the information to it
                auto &cell_qpts = std::get<0>(current_c.first->second);
                auto &cell_maps = std::get<1>(current_c.first->second);
                auto &cell_pts = std::get<2>(current_c.first->second);
                auto &cell_ranks = std::get<3>(current_c.first->second);
                cell_qpts.insert(cell_qpts.end(),
                                 in_qpoints[c].begin(),
                                 in_qpoints[c].end());
                cell_maps.insert(cell_maps.end(),
                                 in_maps[c].begin(),
                                 in_maps[c].end());
                cell_pts.insert(cell_pts.end(),
                                in_points[c].begin(),
                                in_points[c].end());
                std::vector< unsigned int > ranks_tmp(in_points[c].size(),in_rank);
                cell_ranks.insert(cell_ranks.end(),
                                  ranks_tmp.begin(),
                                  ranks_tmp.end());
              }
          }
      }



      // This function initializes the output by calling compute point locations
      // on local points; vector containing points which are probably local.
      // Its output is then sorted in the following manner:
      // - output unmap: points, with relative information, inside locally onwed cells,
      // - ghost loc pts: points, with relative information, inside ghost cells,
      // - classified pts: vector of all points returned in output map and ghost loc pts
      //   (these are stored as indices)
      template <int dim, int spacedim>
      void
      compute_and_classify_points(
        const GridTools::Cache<dim,spacedim>                              &cache,
        const std::vector<Point<spacedim> >                               &local_points,
        const std::vector< unsigned int >                                 &local_points_idx,
        std::unordered_map<
        typename Triangulation<dim, spacedim>::active_cell_iterator,
        std::tuple<
        std::vector< Point<dim> >,
        std::vector< unsigned int >,
        std::vector< Point<spacedim> >,
        std::vector< unsigned int >
        >,
        cell_hash<dim,spacedim>>                                          &output_unmap,
        std::map< unsigned int,
        std::tuple<
        std::vector< CellId >,
        std::vector< std::vector< Point<dim> > >,
        std::vector< std::vector< unsigned int > >,
        std::vector< std::vector< Point<spacedim> > >
        > >                                                               &ghost_loc_pts,
        std::vector< unsigned int >                                       &classified_pts
      )
      {
        auto cpt_loc_pts = compute_point_locations_unmap(cache,local_points);

        // Alayzing the output discarding artificial cell
        // and storing in the proper container locally owned and ghost cells
        for (auto const &cell_tuples : cpt_loc_pts)
          {
            auto &cell_loc = cell_tuples.first;
            auto &q_loc = std::get<0>(cell_tuples.second);
            auto &indices_loc = std::get<1>(cell_tuples.second);
            if (cell_loc->is_locally_owned() )
              {
                // Point inside locally owned cell: storing all its data
                std::vector < Point<spacedim> > cell_points(indices_loc.size());
                for (unsigned int i=0; i< indices_loc.size(); ++i)
                  {
                    // Adding the point to the cell points
                    cell_points[i] = local_points[indices_loc[i]];
                    // Storing the index: notice indices loc refer to the local points
                    // vector, but we need to return the index with respect of
                    // the points owned by the current process
                    classified_pts.emplace_back(local_points_idx[indices_loc[i]]);
                  }
                output_unmap.emplace(std::make_pair(cell_loc,
                                                    std::make_tuple(q_loc,
                                                                    indices_loc,
                                                                    cell_points,
                                                                    std::vector<unsigned int>
                                                                    (indices_loc.size(),cell_loc->subdomain_id()))));
              }
            else if ( cell_loc->is_ghost() )
              {
                // Point inside ghost cell: storing all its information and preparing
                // it to be sent
                std::vector < Point<spacedim> > cell_points(indices_loc.size());
                for (unsigned int i=0; i< indices_loc.size(); ++i)
                  {
                    cell_points[i] = local_points[indices_loc[i]];
                    classified_pts.emplace_back(local_points_idx[indices_loc[i]]);
                  }
                // Each key of the following map represent a process,
                // each mapped value is a tuple containing the information to be sent:
                // preparing the output for the owner, which has rank subdomain id
                auto &map_tuple_owner = ghost_loc_pts[cell_loc->subdomain_id()];
                // To identify the cell on the other process we use the cell id
                std::get<0>(map_tuple_owner).emplace_back(cell_loc->id());
                std::get<1>(map_tuple_owner).emplace_back(q_loc);
                std::get<2>(map_tuple_owner).emplace_back(indices_loc);
                std::get<3>(map_tuple_owner).emplace_back(cell_points);
              }
            // else: the cell is artificial, nothing to do
          }
      }



      // Given the map obtained from a communication, where the key is rank and the mapped
      // value is a pair of (points,indices), calls compute point locations; its output
      // is then merged with output tuple
      // if check_owned is set to true only points
      // lying inside locally onwed cells shall be merged, otherwise all points shall be merged.
      template <int dim, int spacedim>
      void
      compute_and_merge_from_map(
        const GridTools::Cache<dim,spacedim>                                         &cache,
        const   std::map< unsigned int,
        std::pair<
        std::vector < Point<spacedim> >,
        std::vector < unsigned int > >
        >                                                                            &map_pts,
        std::unordered_map< typename Triangulation<dim, spacedim>::active_cell_iterator,
        std::tuple<
        std::vector< Point<dim> >,
        std::vector< unsigned int >,
        std::vector< Point<spacedim> >,
        std::vector< unsigned int >
        >,
        cell_hash<dim,spacedim>>                                                     &output_unmap,
        const bool                                                                   &check_owned
      )
      {
        bool no_check = !check_owned;

        // rank and points is a pair: first rank, then a pair of vectors (points, indices)
        for (auto const &rank_and_points : map_pts)
          {
            // Rewriting the contents of the map in human readable format
            const auto &received_process = rank_and_points.first;
            const auto &received_points = rank_and_points.second.first;
            const auto &received_ranks = rank_and_points.second.second;

            // Initializing the vectors needed to store the result of compute point location
            std::vector< typename Triangulation<dim, spacedim>::active_cell_iterator > in_cell;
            std::vector< std::vector< Point<dim> > > in_qpoints;
            std::vector< std::vector< unsigned int > > in_maps;
            std::vector< std::vector< Point<spacedim> > > in_points;

            auto cpt_loc_pts = compute_point_locations_unmap(cache,rank_and_points.second.first);
            for (const auto &map_c_pt_idx: cpt_loc_pts)
              {
                // Human-readable variables:
                const auto &proc_cell = map_c_pt_idx.first;
                const auto &proc_qpoints = map_c_pt_idx.second.first;
                const auto &proc_maps = map_c_pt_idx.second.second;

                // This is stored either if we're not checking if the cell is owned or
                // if the cell is locally owned
                if ( no_check || proc_cell->is_locally_owned() )
                  {
                    in_cell.emplace_back(proc_cell);
                    in_qpoints.emplace_back(proc_qpoints);
                    // The other two vectors need to be built
                    unsigned int loc_size = proc_qpoints.size();
                    std::vector< unsigned int > cell_maps(loc_size);
                    std::vector< Point<spacedim> > cell_points(loc_size);
                    for (unsigned int pt=0; pt<loc_size; ++pt)
                      {
                        cell_maps[pt] = received_ranks[proc_maps[pt]];
                        cell_points[pt] = received_points[proc_maps[pt]];
                      }
                    in_maps.emplace_back(cell_maps);
                    in_points.emplace_back(cell_points);
                  }
              }

            // Merge everything from the current process
            internal::distributed_cptloc::merge_cptloc_outputs(output_unmap,
                                                               in_cell,
                                                               in_qpoints,
                                                               in_maps,
                                                               in_points,
                                                               received_process);
          }
      }
    } // namespace distributed_cptloc
  } // namespace internal



  template <int dim, int spacedim>
  std::tuple<
  std::vector< typename Triangulation<dim, spacedim>::active_cell_iterator >,
      std::vector< std::vector< Point<dim> > >,
      std::vector< std::vector< unsigned int > >,
      std::vector< std::vector< Point<spacedim> > >,
      std::vector< std::vector< unsigned int > >
      >
      distributed_compute_point_locations
      (const GridTools::Cache<dim,spacedim>                &cache,
       const std::vector<Point<spacedim> >                 &local_points,
       const std::vector< BoundingBox<spacedim> >          &local_bbox)
  {
#ifndef DEAL_II_WITH_MPI
    (void)cache;
    (void)local_points;
    (void)local_bbox;
    Assert(false, ExcMessage("GridTools::distributed_compute_point_locations() requires MPI."));
    std::tuple<
    std::vector< typename Triangulation<dim, spacedim>::active_cell_iterator >,
        std::vector< std::vector< Point<dim> > >,
        std::vector< std::vector< unsigned int > >,
        std::vector< std::vector< Point<spacedim> > >,
        std::vector< std::vector< unsigned int > >
        > tup;
    return tup;
#else
    // Recovering the mpi communicator used to create the triangulation
    const auto &tria_mpi =
      dynamic_cast< const parallel::Triangulation< dim, spacedim >*>(&cache.get_triangulation());
    // If the dynamic cast failed we can't recover the mpi communicator: throwing an assertion error
    Assert(tria_mpi, ExcMessage("GridTools::distributed_compute_point_locations() requires a parallel triangulation."));
    auto mpi_communicator = tria_mpi->get_communicator();
    // Preparing the output tuple
    std::tuple<
    std::vector< typename Triangulation<dim, spacedim>::active_cell_iterator >,
        std::vector< std::vector< Point<dim> > >,
        std::vector< std::vector< unsigned int > >,
        std::vector< std::vector< Point<spacedim> > >,
        std::vector< std::vector< unsigned int > >
        >                                                           output_tuple;

    // Preparing the temporary unordered map
    std::unordered_map< typename Triangulation<dim, spacedim>::active_cell_iterator,
        std::tuple<
        std::vector< Point<dim> >,
        std::vector< unsigned int >,
        std::vector< Point<spacedim> >,
        std::vector< unsigned int >
        >,
        internal::distributed_cptloc::cell_hash<dim,spacedim> >
        temporary_unmap;

    // Obtaining the global mesh description through an all to all communication
    std::vector< std::vector< BoundingBox<spacedim> > > global_bounding_boxes;
    global_bounding_boxes = Utilities::MPI::all_gather(mpi_communicator,local_bbox);

    // Step 1 (part 1): Using the bounding boxes to guess the owner of each points
    // in local_points
    unsigned int my_rank = Utilities::MPI::this_mpi_process(mpi_communicator);

    // Using global bounding boxes to guess/find owner/s of each point
    std::tuple< std::vector< std::vector< unsigned int > >, std::map< unsigned int, unsigned int >,
        std::map< unsigned int, std::vector< unsigned int > > > guessed_points;
    guessed_points =
      GridTools::guess_point_owner(global_bounding_boxes, local_points);

    // Preparing to call compute point locations on points which are/might be
    // local
    const auto &guess_loc_idx = std::get<0>(guessed_points)[my_rank];
    const unsigned int n_local_guess = guess_loc_idx.size();
    // Vector containing points which are probably local
    std::vector< Point<spacedim> > guess_local_pts(n_local_guess);
    for (unsigned int i=0; i<n_local_guess; ++i)
      guess_local_pts[i] = local_points[ guess_loc_idx[i] ];

    // Preparing the map with data on points lying on locally owned cells
    std::map< unsigned int,
        std::tuple<
        std::vector< CellId >,
        std::vector< std::vector< Point<dim> > >,
        std::vector< std::vector< unsigned int > >,
        std::vector< std::vector< Point<spacedim> > > > >  ghost_loc_pts;
    // Vector containing indices of points lying either on locally owned
    // cells or ghost cells, to avoid computing them more than once
    std::vector< unsigned int >                          classified_pts;

    // Thread used to call compute point locations on guess local pts
    Threads::Task<void>
    cpt_loc_tsk
      = Threads::new_task (
          &internal::distributed_cptloc::compute_and_classify_points<dim,spacedim>,
          cache,
          guess_local_pts,
          guess_loc_idx,
          temporary_unmap,
          ghost_loc_pts,
          classified_pts);

    // Step 1 (part 2): communicate point which are owned by a certain process
    // Preparing the map with points whose owner is known with certainty:
    const auto &other_owned_idx = std::get<1>(guessed_points);
    std::map<
    unsigned int,
             std::pair< std::vector<Point<spacedim>> , std::vector<unsigned int > > >
             other_owned_pts;

    for (const auto &indices: other_owned_idx)
      if (indices.second != my_rank)
        {
          // Finding/adding in the map the current process
          auto &current_pts = other_owned_pts[indices.second];
          current_pts.first.emplace_back(local_points[indices.first]);
          current_pts.second.emplace_back(indices.first);
        }

    // Communicating the points whose owner is sure
    auto owned_rank_pts = Utilities::MPI::some_to_some(mpi_communicator,other_owned_pts);
    // Waiting for part 1 to finish to avoid concurrency problems
    cpt_loc_tsk.join();

    // Step 2 (part 1): compute received points which are owned
    Threads::Task<void>
    owned_pts_tsk
      = Threads::new_task (&internal::distributed_cptloc::compute_and_merge_from_map<dim,spacedim>,
                           cache,
                           owned_rank_pts,
                           temporary_unmap,
                           false);

    // Step 2 (part 2): communicate info on points lying on ghost cells
    auto cpt_ghost = Utilities::MPI::some_to_some(mpi_communicator,ghost_loc_pts);

    // Step 3: construct vectors containing uncertain points i.e. those whose owner
    // is known among few guesses
    std::map<
    unsigned int,
             std::pair< std::vector < Point<spacedim> >,
             std::vector<unsigned int > > >
             other_check_pts;

    const auto &other_check_idx = std::get<2>(guessed_points);

    // Points in classified pts need not to be communicated;
    // sorting the array classified pts in order to use
    // binary search when checking if the points needs to be
    // communicated
    // Notice classified pts is a vector of integer indexes
    std::sort (classified_pts.begin(), classified_pts.end());

    for (const auto &pt_to_guesses: other_check_idx)
      {
        if ( !std::binary_search(
               classified_pts.begin(), classified_pts.end(),pt_to_guesses.first) )
          // The point wasn't found in ghost or locally owned cells: adding it to the map
          for (unsigned int rank=0; rank<pt_to_guesses.second.size(); ++rank)
            if (pt_to_guesses.second[rank] != my_rank)
              {
                auto &current_pts = other_check_pts[pt_to_guesses.second[rank]];
                current_pts.first.emplace_back(local_points[pt_to_guesses.first]);
                current_pts.second.emplace_back(pt_to_guesses.second[rank]);
              }
      }

    // Step 4: send around uncertain points
    auto check_pts = Utilities::MPI::some_to_some(mpi_communicator,other_check_pts);
    // Before proceeding, merging threads to avoid concurrency problems
    owned_pts_tsk.join();

    // Step 5: add the received ghost cell data to output
    for ( const auto &rank_vals: cpt_ghost)
      {
        // Transforming CellsIds into Tria iterators
        const auto &cell_ids = std::get<0>(rank_vals.second);
        unsigned int n_cells = cell_ids.size();
        std::vector< typename Triangulation<dim, spacedim>::active_cell_iterator >
        cell_iter(n_cells);
        for (unsigned int c=0; c<n_cells; ++c)
          cell_iter[c] = cell_ids[c].to_cell(cache.get_triangulation());

        internal::distributed_cptloc::merge_cptloc_outputs(temporary_unmap,
                                                           cell_iter,
                                                           std::get<1>(rank_vals.second),
                                                           std::get<2>(rank_vals.second),
                                                           std::get<3>(rank_vals.second),
                                                           rank_vals.first);
      }

    // Step 6: use compute point locations on the uncertain points and
    // merge output
    internal::distributed_cptloc::compute_and_merge_from_map(
      cache,
      check_pts,
      temporary_unmap,
      true);

    // Copying data from the unordered map to the tuple
    // and returning output
    unsigned int size_output = temporary_unmap.size();
    auto &out_cells   = std::get<0>(output_tuple);
    auto &out_qpoints = std::get<1>(output_tuple);
    auto &out_maps    = std::get<2>(output_tuple);
    auto &out_points  = std::get<3>(output_tuple);
    auto &out_ranks   = std::get<4>(output_tuple);

    out_cells.resize(size_output);
    out_qpoints.resize(size_output);
    out_maps.resize(size_output);
    out_points.resize(size_output);
    out_ranks.resize(size_output);

    unsigned int c = 0;
    for (const auto &rank_and_tuple: temporary_unmap)
      {
        out_cells[c]   = rank_and_tuple.first;
        out_qpoints[c] = std::get<0>(rank_and_tuple.second);
        out_maps[c]    = std::get<1>(rank_and_tuple.second);
        out_points[c]  = std::get<2>(rank_and_tuple.second);
        out_ranks[c]   = std::get<3>(rank_and_tuple.second);
        ++c;
      }

    return output_tuple;
#endif
  }


  template<int dim, int spacedim>
  std::map<unsigned int, Point<spacedim> >
  extract_used_vertices(const Triangulation<dim, spacedim> &container,
                        const Mapping<dim,spacedim> &mapping)
  {
    std::map<unsigned int, Point<spacedim> > result;
    for (const auto &cell : container.active_cell_iterators())
      {
        const auto vs = mapping.get_vertices(cell);
        for (unsigned int i=0; i<vs.size(); ++i)
          result[cell->vertex_index(i)]=vs[i];
      }
    Assert(result.size() == container.n_used_vertices(),
           ExcInternalError());
    return result;
  }


  template<int spacedim>
  unsigned int
  find_closest_vertex(const std::map<unsigned int,Point<spacedim> > &vertices,
                      const Point<spacedim> &p)
  {
    auto id_and_v = std::min_element(vertices.begin(), vertices.end(),
                                     [&](const std::pair<const unsigned int, Point<spacedim>> &p1,
                                         const std::pair<const unsigned int, Point<spacedim>> &p2) -> bool
    {
      return p1.second.distance(p) < p2.second.distance(p);
    }
                                    );
    return id_and_v->first;
  }


  template<int dim, int spacedim>
  std::pair<typename Triangulation<dim,spacedim>::active_cell_iterator, Point<dim> >
  find_active_cell_around_point(const Cache<dim,spacedim> &cache,
                                const Point<spacedim> &p,
                                const typename Triangulation<dim,spacedim>::active_cell_iterator &cell_hint,
                                const std::vector<bool> &marked_vertices)
  {
    const auto &mesh = cache.get_triangulation();
    const auto &mapping = cache.get_mapping();
    const auto &vertex_to_cells = cache.get_vertex_to_cell_map();
    const auto &vertex_to_cell_centers = cache.get_vertex_to_cell_centers_directions();

    return find_active_cell_around_point(mapping, mesh, p,
                                         vertex_to_cells,
                                         vertex_to_cell_centers,
                                         cell_hint,
                                         marked_vertices);
  }

  template<int spacedim>
  std::vector< std::vector< BoundingBox<spacedim> > >
  exchange_local_bounding_boxes(const std::vector< BoundingBox<spacedim> > &local_bboxes,
                                MPI_Comm                                    mpi_communicator)
  {
#ifndef DEAL_II_WITH_MPI
    (void)local_bboxes;
    (void)mpi_communicator;
    Assert(false, ExcMessage("GridTools::exchange_local_bounding_boxes() requires MPI."));
    return {};
#else
    // Step 1: preparing data to be sent
    unsigned int n_bboxes = local_bboxes.size();
    // Dimension of the array to be exchanged (number of double)
    int n_local_data = 2*spacedim*n_bboxes;
    // data array stores each entry of each point describing the bounding boxes
    std::vector<double> loc_data_array(n_local_data);
    for (unsigned int i=0; i<n_bboxes; ++i)
      for (unsigned int d=0; d < spacedim; ++d)
        {
          // Extracting the coordinates of each boundary point
          loc_data_array[2*i*spacedim + d] = local_bboxes[i].get_boundary_points().first[d];
          loc_data_array[2*i*spacedim + spacedim + d] = local_bboxes[i].get_boundary_points().second[d];
        }

    // Step 2: exchanging the size of local data
    unsigned int n_procs = Utilities::MPI::n_mpi_processes(mpi_communicator);

    // Vector to store the size of loc_data_array for every process
    std::vector<int> size_all_data(n_procs);

    // Exchanging the number of bboxes
    MPI_Allgather(&n_local_data, 1, MPI_INT,
                  &(size_all_data[0]), 1, MPI_INT,
                  mpi_communicator);

    // Now computing the the displacement, relative to recvbuf,
    // at which to store the incoming data
    std::vector<int> rdispls(n_procs);
    rdispls[0] = 0;
    for (unsigned int i=1; i < n_procs; ++i)
      rdispls[i] = rdispls[i-1] + size_all_data[i-1];

    // Step 3: exchange the data and bounding boxes:
    // Allocating a vector to contain all the received data
    std::vector<double> data_array(rdispls.back() + size_all_data.back());

    MPI_Allgatherv(&(loc_data_array[0]), n_local_data, MPI_DOUBLE,
                   &(data_array[0]), &(size_all_data[0]),
                   &(rdispls[0]), MPI_DOUBLE, mpi_communicator);

    // Step 4: create the array of bboxes for output
    std::vector< std::vector< BoundingBox<spacedim> > > global_bboxes(n_procs);
    unsigned int begin_idx = 0;
    for (unsigned int i=0; i < n_procs; ++i)
      {
        // Number of local bounding boxes
        unsigned int n_bbox_i = size_all_data[i]/(spacedim*2);
        global_bboxes[i].resize(n_bbox_i);
        for (unsigned int bbox=0; bbox<n_bbox_i; ++bbox)
          {
            Point<spacedim> p1,p2; // boundary points for bbox
            for (unsigned int d=0; d<spacedim; ++d)
              {
                p1[d] = data_array[ begin_idx + 2*bbox*spacedim + d];
                p2[d] = data_array[ begin_idx + 2*bbox*spacedim + spacedim + d];
              }
            BoundingBox<spacedim> loc_bbox(std::make_pair(p1,p2));
            global_bboxes[i][bbox] = loc_bbox;
          }
        // Shifting the first index to the start of the next vector
        begin_idx += size_all_data[i];
      }
    return global_bboxes;
#endif // DEAL_II_WITH_MPI
  }


} /* namespace GridTools */


// explicit instantiations
#include "grid_tools.inst"

DEAL_II_NAMESPACE_CLOSE
