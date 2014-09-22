// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2014 by the deal.II authors
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
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/multigrid/mg_dof_handler.h>

#include <cmath>
#include <numeric>


DEAL_II_NAMESPACE_OPEN


namespace GridTools
{

// This anonymous namespace contains utility functions to extract the
// triangulation from any container such as DoFHandler
// and the like
  namespace
  {
    template<int dim, int spacedim>
    const Triangulation<dim, spacedim> &
    get_tria(const Triangulation<dim, spacedim> &tria)
    {
      return tria;
    }

    template<int dim, int spacedim>
    const Triangulation<dim, spacedim> &
    get_tria(const parallel::distributed::Triangulation<dim, spacedim> &tria)
    {
      return tria;
    }

    template<int dim, template<int, int> class Container, int spacedim>
    const Triangulation<dim,spacedim> &
    get_tria(const Container<dim,spacedim> &container)
    {
      return container.get_tria();
    }


    template<int dim, int spacedim>
    Triangulation<dim, spacedim> &
    get_tria(Triangulation<dim, spacedim> &tria)
    {
      return tria;
    }

    template<int dim, int spacedim>
    Triangulation<dim, spacedim> &
    get_tria(parallel::distributed::Triangulation<dim, spacedim> &tria)
    {
      return tria;
    }

    template<int dim, template<int, int> class Container, int spacedim>
    const Triangulation<dim,spacedim> &
    get_tria(Container<dim,spacedim> &container)
    {
      return container.get_tria();
    }
  }


  template <int dim, int spacedim>
  double
  diameter (const Triangulation<dim, spacedim> &tria)
  {
    // we can't deal with distributed meshes
    // since we don't have all vertices
    // locally. there is one exception,
    // however: if the mesh has never been
    // refined. the way to test this is not to
    // ask tria.n_levels()==1, since this is
    // something that can happen on one
    // processor without being true on
    // all. however, we can ask for the global
    // number of active cells and use that
#ifdef DEAL_II_WITH_P4EST
    if (const parallel::distributed::Triangulation<dim,spacedim> *p_tria
        = dynamic_cast<const parallel::distributed::Triangulation<dim,spacedim>*>(&tria))
      Assert (p_tria->n_global_active_cells() == tria.n_cells(0),
              ExcNotImplemented());
#endif

    // the algorithm used simply
    // traverses all cells and picks
    // out the boundary vertices. it
    // may or may not be faster to
    // simply get all vectors, don't
    // mark boundary vertices, and
    // compute the distances thereof,
    // but at least as the mesh is
    // refined, it seems better to
    // first mark boundary nodes, as
    // marking is O(N) in the number of
    // cells/vertices, while computing
    // the maximal distance is O(N*N)
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

    // now traverse the list of
    // boundary vertices and check
    // distances. since distances are
    // symmetric, we only have to check
    // one half
    double max_distance_sqr = 0;
    std::vector<bool>::const_iterator pi = boundary_vertices.begin();
    const unsigned int N = boundary_vertices.size();
    for (unsigned int i=0; i<N; ++i, ++pi)
      {
        std::vector<bool>::const_iterator pj = pi+1;
        for (unsigned int j=i+1; j<N; ++j, ++pj)
          if ((*pi==true) && (*pj==true) &&
              ((vertices[i]-vertices[j]).square() > max_distance_sqr))
            max_distance_sqr = (vertices[i]-vertices[j]).square();
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

    FE_DGQ<dim,spacedim> dummy_fe(0);
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
    if (const parallel::distributed::Triangulation<dim,spacedim> *p_tria
        = dynamic_cast<const parallel::distributed::Triangulation<dim,spacedim>*>(&triangulation))
      global_volume = Utilities::MPI::sum (local_volume, p_tria->get_communicator());
    else
      global_volume = local_volume;
#else
    global_volume = local_volume;
#endif

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
    // first check which vertices are
    // actually used
    std::vector<bool> vertex_used (vertices.size(), false);
    for (unsigned int c=0; c<cells.size(); ++c)
      for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
        vertex_used[cells[c].vertices[v]] = true;

    // then renumber the vertices that
    // are actually used in the same
    // order as they were beforehand
    const unsigned int invalid_vertex = numbers::invalid_unsigned_int;
    std::vector<unsigned int> new_vertex_numbers (vertices.size(), invalid_vertex);
    unsigned int next_free_number = 0;
    for (unsigned int i=0; i<vertices.size(); ++i)
      if (vertex_used[i] == true)
        {
          new_vertex_numbers[i] = next_free_number;
          ++next_free_number;
        };

    // next replace old vertex numbers
    // by the new ones
    for (unsigned int c=0; c<cells.size(); ++c)
      for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
        cells[c].vertices[v] = new_vertex_numbers[cells[c].vertices[v]];

    // same for boundary data
    for (unsigned int c=0; c<subcelldata.boundary_lines.size(); ++c)
      for (unsigned int v=0; v<GeometryInfo<1>::vertices_per_cell; ++v)
        subcelldata.boundary_lines[c].vertices[v]
          = new_vertex_numbers[subcelldata.boundary_lines[c].vertices[v]];
    for (unsigned int c=0; c<subcelldata.boundary_quads.size(); ++c)
      for (unsigned int v=0; v<GeometryInfo<2>::vertices_per_cell; ++v)
        subcelldata.boundary_quads[c].vertices[v]
          = new_vertex_numbers[subcelldata.boundary_quads[c].vertices[v]];

    // finally copy over the vertices
    // which we really need to a new
    // array and replace the old one by
    // the new one
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
      new_vertex_numbers[i]=i;

    // if the considered_vertices vector is
    // empty, consider all vertices
    if (considered_vertices.size()==0)
      considered_vertices=new_vertex_numbers;

    // now loop over all vertices to be
    // considered and try to find an identical
    // one
    for (unsigned int i=0; i<considered_vertices.size(); ++i)
      {
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

    delete_unused_vertices(vertices,cells,subcelldata);
  }



// define some transformations in an anonymous namespace
  namespace
  {
    template <int spacedim>
    class ShiftPoint
    {
    public:
      ShiftPoint (const Point<spacedim> &shift)
        :
        shift(shift)
      {}
      Point<spacedim> operator() (const Point<spacedim> p) const
      {
        return p+shift;
      }
    private:
      const Point<spacedim> shift;
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


    template <int spacedim>
    class ScalePoint
    {
    public:
      ScalePoint (const double factor)
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
  shift (const Point<spacedim>   &shift_vector,
         Triangulation<dim, spacedim> &triangulation)
  {
    transform (ShiftPoint<spacedim>(shift_vector), triangulation);
  }



  void
  rotate (const double      angle,
          Triangulation<2> &triangulation)
  {
    transform (Rotate2d(angle), triangulation);
  }



  template <int dim, int spacedim>
  void
  scale (const double        scaling_factor,
         Triangulation<dim, spacedim> &triangulation)
  {
    Assert (scaling_factor>0, ExcScalingFactorNotPositive (scaling_factor));
    transform (ScalePoint<spacedim>(scaling_factor), triangulation);
  }


  template <int dim>
  void
  laplace_transform (const std::map<unsigned int,Point<dim> > &new_points,
                     Triangulation<dim> &triangulation,
                     const Function<dim> *coefficient)
  {
    //TODO: Move implementation of this function into the current
    // namespace
    GridGenerator::laplace_transformation(triangulation, new_points, coefficient);
  }



  template <int dim, int spacedim>
  void
  distort_random (const double        factor,
                  Triangulation<dim, spacedim> &triangulation,
                  const bool          keep_boundary)
  {
    //TODO: Move implementation of this function into the current
    // namespace
    triangulation.distort_random (factor, keep_boundary);
  }


  template <int dim, template <int, int> class Container, int spacedim>
  unsigned int
  find_closest_vertex (const Container<dim,spacedim> &container,
                       const Point<spacedim> &p)
  {
    // first get the underlying
    // triangulation from the
    // container and determine vertices
    // and used vertices
    const Triangulation<dim, spacedim> &tria = get_tria(container);

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
    double       best_dist   = (p - vertices[best_vertex]).square();

    // For all remaining vertices, test
    // whether they are any closer
    for (unsigned int j = best_vertex+1; j < vertices.size(); j++)
      if (used[j])
        {
          double dist = (p - vertices[j]).square();
          if (dist < best_dist)
            {
              best_vertex = j;
              best_dist   = dist;
            }
        }

    return best_vertex;
  }


  template<int dim, template<int, int> class Container, int spacedim>
  std::vector<typename Container<dim,spacedim>::active_cell_iterator>
  find_cells_adjacent_to_vertex(const Container<dim,spacedim> &container,
                                const unsigned int    vertex)
  {
    // make sure that the given vertex is
    // an active vertex of the underlying
    // triangulation
    Assert(vertex < get_tria(container).n_vertices(),
           ExcIndexRange(0,get_tria(container).n_vertices(),vertex));
    Assert(get_tria(container).get_used_vertices()[vertex],
           ExcVertexNotUsed(vertex));

    // use a set instead of a vector
    // to ensure that cells are inserted only
    // once
    std::set<typename Container<dim,spacedim>::active_cell_iterator> adjacent_cells;

    typename Container<dim,spacedim>::active_cell_iterator
    cell = container.begin_active(),
    endc = container.end();

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
      std::vector<typename Container<dim,spacedim>::active_cell_iterator>
      (adjacent_cells.begin(), adjacent_cells.end());
  }



  namespace
  {
    template <int dim, template<int, int> class Container, int spacedim>
    void find_active_cell_around_point_internal(const Container<dim,spacedim> &container,
                                                std::set<typename Container<dim,spacedim>::active_cell_iterator> &searched_cells,
                                                std::set<typename Container<dim,spacedim>::active_cell_iterator> &adjacent_cells)
    {
      typedef typename Container<dim,spacedim>::active_cell_iterator cell_iterator;

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
          get_active_neighbors<Container<dim, spacedim> >(*cell, active_neighbors);
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
          cell_iterator it = container.begin_active();
          for ( ; it!=container.end(); ++it)
            if (searched_cells.find(it) == searched_cells.end())
              {
                adjacent_cells.insert(it);
                break;
              }
        }
    }
  }

  template <int dim, template<int, int> class Container, int spacedim>
  typename Container<dim,spacedim>::active_cell_iterator
  find_active_cell_around_point (const Container<dim,spacedim>  &container,
                                 const Point<spacedim> &p)
  {
    return
      find_active_cell_around_point<dim,Container,spacedim>
      (StaticMappingQ1<dim,spacedim>::mapping,
       container, p).first;
  }


  template <int dim, template <int, int> class Container, int spacedim>
  std::pair<typename Container<dim,spacedim>::active_cell_iterator, Point<dim> >
  find_active_cell_around_point (const Mapping<dim,spacedim>   &mapping,
                                 const Container<dim,spacedim> &container,
                                 const Point<spacedim>     &p)
  {
    typedef typename Container<dim,spacedim>::active_cell_iterator active_cell_iterator;

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
      = find_cells_adjacent_to_vertex(container,
                                      find_closest_vertex(container, p));

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
    const unsigned int n_active_cells = get_tria(container).n_active_cells();
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
            catch (typename MappingQ1<dim,spacedim>::ExcTransformationFailed &)
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
            find_active_cell_around_point_internal<dim,Container,spacedim>
            (container, searched_cells, adjacent_cells);
          }
      }

    AssertThrow (best_cell.first.state() == IteratorState::valid,
                 ExcPointNotFound<spacedim>(p));

    return best_cell;
  }



  template <int dim, int spacedim>
  std::pair<typename hp::DoFHandler<dim,spacedim>::active_cell_iterator, Point<dim> >
  find_active_cell_around_point (const hp::MappingCollection<dim,spacedim>   &mapping,
                                 const hp::DoFHandler<dim,spacedim> &container,
                                 const Point<spacedim>     &p)
  {
    Assert ((mapping.size() == 1) ||
            (mapping.size() == container.get_fe().size()),
            ExcMessage ("Mapping collection needs to have either size 1 "
                        "or size equal to the number of elements in "
                        "the FECollection."));

    typedef typename hp::DoFHandler<dim,spacedim>::active_cell_iterator cell_iterator;

    std::pair<cell_iterator, Point<dim> > best_cell;
    //If we have only one element in the MappingCollection,
    //we use find_active_cell_around_point using only one
    //mapping.
    if (mapping.size() == 1)
      best_cell = find_active_cell_around_point(mapping[0], container, p);
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
        unsigned int vertex = find_closest_vertex(container, p);

        std::vector<cell_iterator> adjacent_cells_tmp =
          find_cells_adjacent_to_vertex(container, vertex);

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
        const unsigned int n_cells =get_tria(container).n_cells();
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
                catch (typename MappingQ1<dim,spacedim>::ExcTransformationFailed &)
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
                (container, searched_cells, adjacent_cells);
              }

          }
      }

    AssertThrow (best_cell.first.state() == IteratorState::valid,
                 ExcPointNotFound<spacedim>(p));

    return best_cell;
  }



  template <int dim, int spacedim>
  void
  get_face_connectivity_of_cells (const Triangulation<dim,spacedim> &triangulation,
                                  SparsityPattern          &cell_connectivity)
  {
    // as built in this function, we
    // only consider face neighbors,
    // which leads to a fixed number of
    // entries per row (don't forget
    // that each cell couples with
    // itself, and that neighbors can
    // be refined)
    cell_connectivity.reinit (triangulation.n_active_cells(),
                              triangulation.n_active_cells(),
                              GeometryInfo<dim>::faces_per_cell
                              * GeometryInfo<dim>::max_children_per_face
                              +
                              1);

    // create a map pair<lvl,idx> -> SparsityPattern index
    // TODO: we are no longer using user_indices for this because we can get
    // pointer/index clashes when saving/restoring them. The following approach
    // works, but this map can get quite big. Not sure about more efficient solutions.
    std::map< std::pair<unsigned int,unsigned int>, unsigned int >
    indexmap;
    unsigned int index = 0;
    for (typename Triangulation<dim,spacedim>::active_cell_iterator
         cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell, ++index)
      indexmap[std::pair<unsigned int,unsigned int>(cell->level(),cell->index())] = index;

    // next loop over all cells and
    // their neighbors to build the
    // sparsity pattern. note that it's
    // a bit hard to enter all the
    // connections when a neighbor has
    // children since we would need to
    // find out which of its children
    // is adjacent to the current
    // cell. this problem can be
    // omitted if we only do something
    // if the neighbor has no children
    // -- in that case it is either on
    // the same or a coarser level than
    // we are. in return, we have to
    // add entries in both directions
    // for both cells
    index = 0;
    for (typename Triangulation<dim,spacedim>::active_cell_iterator
         cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell, ++index)
      {
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

    // now compress the so-built connectivity pattern
    cell_connectivity.compress ();
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
        for (typename Triangulation<dim,spacedim>::active_cell_iterator
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
        for (typename Triangulation<dim,spacedim>::active_cell_iterator
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
    unsigned int index = 0;
    for (typename Triangulation<dim,spacedim>::active_cell_iterator
         cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell, ++index)
      cell->set_subdomain_id (partition_indices[index]);
  }



  template <int dim, int spacedim>
  void
  get_subdomain_association (const Triangulation<dim, spacedim>  &triangulation,
                             std::vector<types::subdomain_id> &subdomain)
  {
    Assert (subdomain.size() == triangulation.n_active_cells(),
            ExcDimensionMismatch (subdomain.size(),
                                  triangulation.n_active_cells()));
    unsigned int index = 0;
    for (typename Triangulation<dim, spacedim>::active_cell_iterator
         cell = triangulation.begin_active();
         cell!=triangulation.end(); ++cell, ++index)
      subdomain[index] = cell->subdomain_id();

    Assert (index == subdomain.size(), ExcInternalError());
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



  template <typename Container>
  std::list<std::pair<typename Container::cell_iterator,
      typename Container::cell_iterator> >
      get_finest_common_cells (const Container &mesh_1,
                               const Container &mesh_2)
  {
    Assert (have_same_coarse_mesh (mesh_1, mesh_2),
            ExcMessage ("The two containers must be represent triangulations that "
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
    std::list<std::pair<typename Container::cell_iterator,
        typename Container::cell_iterator> >
        CellList;

    CellList cell_list;

    // first push the coarse level cells
    typename Container::cell_iterator
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



  template <typename Container>
  bool
  have_same_coarse_mesh (const Container &mesh_1,
                         const Container &mesh_2)
  {
    return have_same_coarse_mesh (get_tria(mesh_1),
                                  get_tria(mesh_2));
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



  template <int dim, int spacedim>
  void
  create_union_triangulation (const Triangulation<dim, spacedim> &triangulation_1,
                              const Triangulation<dim, spacedim> &triangulation_2,
                              Triangulation<dim, spacedim>       &result)
  {
    // this function is deprecated. call the function that replaced it
    GridGenerator::create_union_triangulation (triangulation_1, triangulation_2, result);
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

                Point<spacedim> h;
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
              old_min_product = std::min (old_min_product,
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
              new_min_product = std::min (new_min_product,
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



  template <class Container>
  std::vector<typename Container::active_cell_iterator>
  get_patch_around_cell(const typename Container::active_cell_iterator &cell)
  {
    Assert (cell->is_locally_owned(),
            ExcMessage ("This function only makes sense if the cell for "
                        "which you are asking for a patch, is locally "
                        "owned."));

    std::vector<typename Container::active_cell_iterator> patch;
    patch.push_back (cell);
    for (unsigned int face_number=0; face_number<GeometryInfo<Container::dimension>::faces_per_cell; ++face_number)
      if (cell->face(face_number)->at_boundary()==false)
        {
          if (cell->neighbor(face_number)->has_children() == false)
            patch.push_back (cell->neighbor(face_number));
          else
            // the neighbor is refined. in 2d/3d, we can simply ask for the children
            // of the neighbor because they can not be further refined and,
            // consequently, the children is active
            if (Container::dimension > 1)
              {
                for (unsigned int subface=0; subface<cell->face(face_number)->n_children(); ++subface)
                  patch.push_back (cell->neighbor_child_on_subface (face_number, subface));
              }
            else
              {
                // in 1d, we need to work a bit harder: iterate until we find
                // the child by going from cell to child to child etc
                typename Container::cell_iterator neighbor
                  = cell->neighbor (face_number);
                while (neighbor->has_children())
                  neighbor = neighbor->child(1-face_number);

                Assert (neighbor->neighbor(1-face_number) == cell, ExcInternalError());
                patch.push_back (neighbor);
              }
        }
    return patch;
  }




  template <template <int,int> class Container, int dim, int spacedim>
#ifndef _MSC_VER
  std::map<typename Container<dim-1,spacedim>::cell_iterator,
      typename Container<dim,spacedim>::face_iterator>
#else
  typename ExtractBoundaryMesh<Container,dim,spacedim>::return_type
#endif
      extract_boundary_mesh (const Container<dim,spacedim> &volume_mesh,
                             Container<dim-1,spacedim>     &surface_mesh,
                             const std::set<types::boundary_id> &boundary_ids)
  {
    // this function is deprecated. call the one that replaced it
    return GridGenerator::extract_boundary_mesh (volume_mesh, surface_mesh, boundary_ids);
  }



  /*
   * Internally used in orthogonal_equality
   *
   * An orthogonal equality test for points:
   *
   * point1 and point2 are considered equal, if
   *    (point1 + offset) - point2
   * is parallel to the unit vector in <direction>
   */
  template<int spacedim>
  inline bool orthogonal_equality (const dealii::Point<spacedim> &point1,
                                   const dealii::Point<spacedim> &point2,
                                   const int                     direction,
                                   const dealii::Tensor<1,spacedim> &offset)
  {
    Assert (0<=direction && direction<spacedim,
            ExcIndexRange (direction, 0, spacedim));
    for (int i = 0; i < spacedim; ++i)
      {
        // Only compare coordinate-components != direction:
        if (i == direction)
          continue;

        if (fabs(point1(i) + offset[i] - point2(i)) > 1.e-10)
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
      return 4; // [true ,false,false]
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
                       const dealii::Tensor<1,FaceIterator::AccessorType::space_dimension> &offset)
  {
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
                                    direction, offset))
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
                       const dealii::Tensor<1,FaceIterator::AccessorType::space_dimension> &offset)
  {
    // Call the function above with a dummy orientation array
    std::bitset<3> dummy;
    return orthogonal_equality (dummy, face1, face2, direction, offset);
  }



  /*
   * Internally used in collect_periodic_faces
   */
  template<typename CellIterator>
  void
  match_periodic_face_pairs
  (std::set<std::pair<CellIterator, unsigned int> > &pairs1,
   std::set<std::pair<typename identity<CellIterator>::type, unsigned int> > &pairs2,
   const int direction,
   std::vector<PeriodicFacePair<CellIterator> > &matched_pairs,
   const dealii::Tensor<1,CellIterator::AccessorType::space_dimension> &offset)
  {
    static const int space_dim = CellIterator::AccessorType::space_dimension;
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
                                               direction, offset))
              {
                // We have a match, so insert the matching pairs and
                // remove the matched cell in pairs2 to speed up the
                // matching:
                const PeriodicFacePair<CellIterator> matched_face
                = {{cell1, cell2},{face_idx1, face_idx2}, orientation};
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



  template<typename CONTAINER>
  void
  collect_periodic_faces
  (const CONTAINER          &container,
   const types::boundary_id b_id1,
   const types::boundary_id b_id2,
   const int                direction,
   std::vector<PeriodicFacePair<typename CONTAINER::cell_iterator> > &matched_pairs,
   const dealii::Tensor<1,CONTAINER::space_dimension> &offset)
  {
    static const int dim = CONTAINER::dimension;
    static const int space_dim = CONTAINER::space_dimension;
    Assert (0<=direction && direction<space_dim,
            ExcIndexRange (direction, 0, space_dim));

    // Loop over all cells on the highest level and collect all boundary
    // faces belonging to b_id1 and b_id2:

    std::set<std::pair<typename CONTAINER::cell_iterator, unsigned int> > pairs1;
    std::set<std::pair<typename CONTAINER::cell_iterator, unsigned int> > pairs2;

    for (typename CONTAINER::cell_iterator cell = container.begin(0);
         cell != container.end(0); ++cell)
      {
        for (unsigned int i = 0; i < GeometryInfo<dim>::faces_per_cell; ++i)
          {
            const typename CONTAINER::face_iterator face = cell->face(i);
            if (face->at_boundary() && face->boundary_indicator() == b_id1)
              {
                const std::pair<typename CONTAINER::cell_iterator, unsigned int> pair1
                  = std::make_pair(cell, i);
                pairs1.insert(pair1);
              }

            if (face->at_boundary() && face->boundary_indicator() == b_id2)
              {
                const std::pair<typename CONTAINER::cell_iterator, unsigned int> pair2
                  = std::make_pair(cell, i);
                pairs2.insert(pair2);
              }
          }
      }

    Assert (pairs1.size() == pairs2.size(),
            ExcMessage ("Unmatched faces on periodic boundaries"));

    // and call match_periodic_face_pairs that does the actual matching:
    match_periodic_face_pairs(pairs1, pairs2, direction, matched_pairs, offset);
  }



  template<typename CONTAINER>
  void
  collect_periodic_faces
  (const CONTAINER          &container,
   const types::boundary_id b_id,
   const int                direction,
   std::vector<PeriodicFacePair<typename CONTAINER::cell_iterator> > &matched_pairs,
   const dealii::Tensor<1,CONTAINER::space_dimension> &offset)
  {
    static const int dim = CONTAINER::dimension;
    static const int space_dim = CONTAINER::space_dimension;
    Assert (0<=direction && direction<space_dim,
            ExcIndexRange (direction, 0, space_dim));

    Assert(dim == space_dim,
           ExcNotImplemented());

    // Loop over all cells on the highest level and collect all boundary
    // faces 2*direction and 2*direction*1:

    std::set<std::pair<typename CONTAINER::cell_iterator, unsigned int> > pairs1;
    std::set<std::pair<typename CONTAINER::cell_iterator, unsigned int> > pairs2;

    for (typename CONTAINER::cell_iterator cell = container.begin(0);
         cell != container.end(0); ++cell)
      {
        const typename CONTAINER::face_iterator face_1 = cell->face(2*direction);
        const typename CONTAINER::face_iterator face_2 = cell->face(2*direction+1);

        if (face_1->at_boundary() && face_1->boundary_indicator() == b_id)
          {
            const std::pair<typename CONTAINER::cell_iterator, unsigned int> pair1
              = std::make_pair(cell, 2*direction);
            pairs1.insert(pair1);
          }

        if (face_2->at_boundary() && face_2->boundary_indicator() == b_id)
          {
            const std::pair<typename CONTAINER::cell_iterator, unsigned int> pair2
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
    match_periodic_face_pairs(pairs1, pairs2, direction, matched_pairs, offset);

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

} /* namespace GridTools */


// explicit instantiations
#include "grid_tools.inst"

DEAL_II_NAMESPACE_CLOSE
