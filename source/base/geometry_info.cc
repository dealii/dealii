// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
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

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/tensor.h>

DEAL_II_NAMESPACE_OPEN


unsigned int
GeometryInfo<0>::n_children(const RefinementCase<0> &)
{
  return 0;
}



template <int dim> const unsigned int GeometryInfo<dim>::max_children_per_cell;
template <int dim> const unsigned int GeometryInfo<dim>::faces_per_cell;
template <int dim> const unsigned int GeometryInfo<dim>::max_children_per_face;
template <int dim> const unsigned int GeometryInfo<dim>::vertices_per_cell;
template <int dim> const unsigned int GeometryInfo<dim>::vertices_per_face;
template <int dim> const unsigned int GeometryInfo<dim>::lines_per_face;
template <int dim> const unsigned int GeometryInfo<dim>::quads_per_face;
template <int dim> const unsigned int GeometryInfo<dim>::lines_per_cell;
template <int dim> const unsigned int GeometryInfo<dim>::quads_per_cell;
template <int dim> const unsigned int GeometryInfo<dim>::hexes_per_cell;


using namespace numbers;

// make sure that also the icc compiler defines (and not only declares)
// these variables
namespace internal
{
  void foo (const unsigned int *) {}

  template <int dim>
  void define_variables ()
  {
    foo(&::dealii::GeometryInfo<dim>::vertices_per_cell);
  }

  template void define_variables<2> ();
  template void define_variables<3> ();
  template void define_variables<4> ();
}



template <>
const unsigned int
GeometryInfo<1>::unit_normal_direction[faces_per_cell]
  = { 0, 0 };

template <>
const unsigned int
GeometryInfo<2>::unit_normal_direction[faces_per_cell]
  = { 0, 0, 1, 1 };

template <>
const unsigned int
GeometryInfo<3>::unit_normal_direction[faces_per_cell]
  = { 0, 0, 1, 1, 2, 2 };

template <>
const unsigned int
GeometryInfo<4>::unit_normal_direction[faces_per_cell]
  = { 0, 0, 1, 1, 2, 2, 3, 3 };



template <>
const int
GeometryInfo<1>::unit_normal_orientation[faces_per_cell]
  = { -1, 1 };

template <>
const int
GeometryInfo<2>::unit_normal_orientation[faces_per_cell]
  = { -1, 1, -1, 1 };

template <>
const int
GeometryInfo<3>::unit_normal_orientation[faces_per_cell]
  = { -1, 1, -1, 1, -1, 1 };

template <>
const int
GeometryInfo<4>::unit_normal_orientation[faces_per_cell]
  = { -1, 1, -1, 1, -1, 1, -1, 1 };



template <>
const unsigned int
GeometryInfo<1>::opposite_face[faces_per_cell]
  = { 1, 0 };

template <>
const unsigned int
GeometryInfo<2>::opposite_face[faces_per_cell]
  = { 1, 0, 3, 2 };

template <>
const unsigned int
GeometryInfo<3>::opposite_face[faces_per_cell]
  = { 1, 0, 3, 2, 5, 4 };

template <>
const unsigned int
GeometryInfo<4>::opposite_face[faces_per_cell]
  = { 1, 0, 3, 2, 5, 4, 7, 6 };



template <>
const unsigned int GeometryInfo<1>::ucd_to_deal[GeometryInfo<1>::vertices_per_cell]
  = { 0, 1};

template <>
const unsigned int GeometryInfo<2>::ucd_to_deal[GeometryInfo<2>::vertices_per_cell]
  = { 0, 1, 3, 2};

template <>
const unsigned int GeometryInfo<3>::ucd_to_deal[GeometryInfo<3>::vertices_per_cell]
  = { 0, 1, 5, 4, 2, 3, 7, 6};

template <>
const unsigned int GeometryInfo<4>::ucd_to_deal[GeometryInfo<4>::vertices_per_cell]
  = {  invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int
    };


template <>
const unsigned int GeometryInfo<1>::dx_to_deal[GeometryInfo<1>::vertices_per_cell]
  = { 0, 1};

template <>
const unsigned int GeometryInfo<2>::dx_to_deal[GeometryInfo<2>::vertices_per_cell]
  = { 0, 2, 1, 3};

template <>
const unsigned int GeometryInfo<3>::dx_to_deal[GeometryInfo<3>::vertices_per_cell]
  = { 0, 4, 2, 6, 1, 5, 3, 7};

template <>
const unsigned int GeometryInfo<4>::dx_to_deal[GeometryInfo<4>::vertices_per_cell]
  = {  invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int,
       invalid_unsigned_int
    };

template <>
const unsigned int GeometryInfo<1>::vertex_to_face
[GeometryInfo<1>::vertices_per_cell][1]
= { { 0 },
  { 1 }
};

template <>
const unsigned int GeometryInfo<2>::vertex_to_face
[GeometryInfo<2>::vertices_per_cell][2]
= { { 0, 2 },
  { 1, 2 },
  { 0, 3 },
  { 1, 3 }
};

template <>
const unsigned int GeometryInfo<3>::vertex_to_face
[GeometryInfo<3>::vertices_per_cell][3]
= { { 0, 2, 4 },
  { 1, 2, 4 },
  { 0, 3, 4 },
  { 1, 3, 4 },
  { 0, 2, 5 },
  { 1, 2, 5 },
  { 0, 3, 5 },
  { 1, 3, 5 }
};

template <>
const unsigned int GeometryInfo<4>::vertex_to_face
[GeometryInfo<4>::vertices_per_cell][4]
= { { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
  { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
  { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
  { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
  { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
  { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
  { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
  { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
  { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
  { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
  { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
  { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
  { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
  { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
  { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int },
  { invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int, invalid_unsigned_int }
};


template<int dim>
unsigned int
GeometryInfo<dim>::n_children(const RefinementCase<dim> &ref_case)
{
  static const unsigned int n_children[RefinementCase<3>::cut_xyz+1]=
  {0, 2, 2, 4, 2, 4, 4, 8};

  return n_children[ref_case];
}


template<>
unsigned int
GeometryInfo<1>::n_subfaces(const internal::SubfaceCase<1> &)
{
  Assert(false, ExcImpossibleInDim(1));
  return 0;
}



template<>
unsigned int
GeometryInfo<2>::n_subfaces(const internal::SubfaceCase<2> &subface_case)
{
  return (subface_case == internal::SubfaceCase<2>::case_x) ? 2 : 0;
}



template<>
unsigned int
GeometryInfo<3>::n_subfaces(const internal::SubfaceCase<3> &subface_case)
{
  static const unsigned int nsubs[internal::SubfaceCase<3>::case_isotropic+1]=
  {0, 2, 3, 3, 4, 2, 3, 3, 4, 4};
  return nsubs[subface_case];
}


template<>
double
GeometryInfo<1>::subface_ratio(const internal::SubfaceCase<1> &,
                               const unsigned int)
{
  return 1;
}


template<>
double
GeometryInfo<2>::subface_ratio(const internal::SubfaceCase<2> &subface_case,
                               const unsigned int)
{
  const unsigned int dim=2;

  double ratio=1;
  switch (subface_case)
    {
    case internal::SubfaceCase<dim>::case_none:
      // Here, an
      // Assert(false,ExcInternalError())
      // would be the right
      // choice, but
      // unfortunately the
      // current function is
      // also called for faces
      // without children (see
      // tests/fe/mapping.cc).
//          Assert(false, ExcMessage("Face has no subfaces."));
      // Furthermore, assign
      // following value as
      // otherwise the
      // bits/volume_x tests
      // break
      ratio=1./GeometryInfo<dim>::max_children_per_face;
      break;
    case internal::SubfaceCase<dim>::case_x:
      ratio=0.5;
      break;
    default:
      // there should be no
      // cases left
      Assert(false, ExcInternalError());
      break;
    }

  return ratio;
}


template<>
double
GeometryInfo<3>::subface_ratio(const internal::SubfaceCase<3> &subface_case,
                               const unsigned int subface_no)
{
  const unsigned int dim=3;

  double ratio=1;
  switch (subface_case)
    {
    case internal::SubfaceCase<dim>::case_none:
      // Here, an
      // Assert(false,ExcInternalError())
      // would be the right
      // choice, but
      // unfortunately the
      // current function is
      // also called for faces
      // without children (see
      // tests/bits/mesh_3d_16.cc). Add
      // following switch to
      // avoid diffs in
      // tests/bits/mesh_3d_16
      ratio=1./GeometryInfo<dim>::max_children_per_face;
      break;
    case internal::SubfaceCase<dim>::case_x:
    case internal::SubfaceCase<dim>::case_y:
      ratio=0.5;
      break;
    case internal::SubfaceCase<dim>::case_xy:
    case internal::SubfaceCase<dim>::case_x1y2y:
    case internal::SubfaceCase<dim>::case_y1x2x:
      ratio=0.25;
      break;
    case internal::SubfaceCase<dim>::case_x1y:
    case internal::SubfaceCase<dim>::case_y1x:
      if (subface_no<2)
        ratio=0.25;
      else
        ratio=0.5;
      break;
    case internal::SubfaceCase<dim>::case_x2y:
    case internal::SubfaceCase<dim>::case_y2x:
      if (subface_no==0)
        ratio=0.5;
      else
        ratio=0.25;
      break;
    default:
      // there should be no
      // cases left
      Assert(false, ExcInternalError());
      break;
    }

  return ratio;
}



template<>
RefinementCase<0>
GeometryInfo<1>::face_refinement_case(const RefinementCase<1> &,
                                      const unsigned int,
                                      const bool,
                                      const bool,
                                      const bool)
{
  Assert(false, ExcImpossibleInDim(1));

  return RefinementCase<0>::no_refinement;
}


template<>
RefinementCase<1>
GeometryInfo<2>::face_refinement_case(const RefinementCase<2> &cell_refinement_case,
                                      const unsigned int face_no,
                                      const bool,
                                      const bool,
                                      const bool)
{
  const unsigned int dim=2;
  Assert(cell_refinement_case<RefinementCase<dim>::isotropic_refinement+1,
         ExcIndexRange(cell_refinement_case, 0, RefinementCase<dim>::isotropic_refinement+1));
  Assert(face_no<GeometryInfo<dim>::faces_per_cell,
         ExcIndexRange(face_no, 0, GeometryInfo<dim>::faces_per_cell));

  static const RefinementCase<dim-1>
  ref_cases[RefinementCase<dim>::isotropic_refinement+1][GeometryInfo<dim>::faces_per_cell/2]=
  {
    {
      RefinementCase<dim-1>::no_refinement,  // no_refinement
      RefinementCase<dim-1>::no_refinement
    },

    {
      RefinementCase<dim-1>::no_refinement,
      RefinementCase<dim-1>::cut_x
    },

    {
      RefinementCase<dim-1>::cut_x,
      RefinementCase<dim-1>::no_refinement
    },

    {
      RefinementCase<dim-1>::cut_x,          // cut_xy
      RefinementCase<dim-1>::cut_x
    }
  };

  return ref_cases[cell_refinement_case][face_no/2];
}


template<>
RefinementCase<2>
GeometryInfo<3>::face_refinement_case(const RefinementCase<3> &cell_refinement_case,
                                      const unsigned int face_no,
                                      const bool face_orientation,
                                      const bool /*face_flip*/,
                                      const bool face_rotation)
{
  const unsigned int dim=3;
  Assert(cell_refinement_case<RefinementCase<dim>::isotropic_refinement+1,
         ExcIndexRange(cell_refinement_case, 0, RefinementCase<dim>::isotropic_refinement+1));
  Assert(face_no<GeometryInfo<dim>::faces_per_cell,
         ExcIndexRange(face_no, 0, GeometryInfo<dim>::faces_per_cell));

  static const RefinementCase<dim-1>
  ref_cases[RefinementCase<dim>::isotropic_refinement+1][GeometryInfo<dim>::faces_per_cell/2]=
  {
    {
      RefinementCase<dim-1>::no_refinement,  // no_refinement
      RefinementCase<dim-1>::no_refinement,
      RefinementCase<dim-1>::no_refinement
    },

    {
      RefinementCase<dim-1>::no_refinement,  // cut_x
      RefinementCase<dim-1>::cut_y,
      RefinementCase<dim-1>::cut_x
    },

    {
      RefinementCase<dim-1>::cut_x,          // cut_y
      RefinementCase<dim-1>::no_refinement,
      RefinementCase<dim-1>::cut_y
    },

    {
      RefinementCase<dim-1>::cut_x,          // cut_xy
      RefinementCase<dim-1>::cut_y,
      RefinementCase<dim-1>::cut_xy
    },

    {
      RefinementCase<dim-1>::cut_y,          // cut_z
      RefinementCase<dim-1>::cut_x,
      RefinementCase<dim-1>::no_refinement
    },

    {
      RefinementCase<dim-1>::cut_y,          // cut_xz
      RefinementCase<dim-1>::cut_xy,
      RefinementCase<dim-1>::cut_x
    },

    {
      RefinementCase<dim-1>::cut_xy,         // cut_yz
      RefinementCase<dim-1>::cut_x,
      RefinementCase<dim-1>::cut_y
    },

    {
      RefinementCase<dim-1>::cut_xy,         // cut_xyz
      RefinementCase<dim-1>::cut_xy,
      RefinementCase<dim-1>::cut_xy
    },
  };

  const RefinementCase<dim-1> ref_case=ref_cases[cell_refinement_case][face_no/2];

  static const RefinementCase<dim-1> flip[4]=
  {
    RefinementCase<dim-1>::no_refinement,
    RefinementCase<dim-1>::cut_y,
    RefinementCase<dim-1>::cut_x,
    RefinementCase<dim-1>::cut_xy
  };

  // correct the ref_case for face_orientation
  // and face_rotation. for face_orientation,
  // 'true' is the default value whereas for
  // face_rotation, 'false' is standard. If
  // <tt>face_rotation==face_orientation</tt>,
  // then one of them is non-standard and we
  // have to swap cut_x and cut_y, otherwise no
  // change is necessary.  face_flip has no
  // influence. however, in order to keep the
  // interface consistent with other functions,
  // we still include it as an argument to this
  // function
  return (face_orientation==face_rotation) ? flip[ref_case] : ref_case;
}



template<>
RefinementCase<1>
GeometryInfo<1>::line_refinement_case(const RefinementCase<1> &cell_refinement_case,
                                      const unsigned int line_no)
{
  const unsigned int dim = 1;
  Assert(cell_refinement_case<RefinementCase<dim>::isotropic_refinement+1,
         ExcIndexRange(cell_refinement_case, 0, RefinementCase<dim>::isotropic_refinement+1));
  Assert(line_no<GeometryInfo<dim>::lines_per_cell,
         ExcIndexRange(line_no, 0, GeometryInfo<dim>::lines_per_cell));

  return cell_refinement_case;
}


template<>
RefinementCase<1>
GeometryInfo<2>::line_refinement_case(const RefinementCase<2> &cell_refinement_case,
                                      const unsigned int line_no)
{
  // Assertions are in face_refinement_case()
  return face_refinement_case(cell_refinement_case, line_no);
}


template<>
RefinementCase<1>
GeometryInfo<3>::line_refinement_case(const RefinementCase<3> &cell_refinement_case,
                                      const unsigned int line_no)
{
  const unsigned int dim=3;
  Assert(cell_refinement_case<RefinementCase<dim>::isotropic_refinement+1,
         ExcIndexRange(cell_refinement_case, 0, RefinementCase<dim>::isotropic_refinement+1));
  Assert(line_no<GeometryInfo<dim>::lines_per_cell,
         ExcIndexRange(line_no, 0, GeometryInfo<dim>::lines_per_cell));

  // array indicating, which simple refine
  // case cuts a line in dirextion x, y or
  // z. For example, cut_y and everything
  // containing cut_y (cut_xy, cut_yz,
  // cut_xyz) cuts lines, which are in y
  // direction.
  static const RefinementCase<dim>
  cut_one[dim] =
  {
    RefinementCase<dim>::cut_x,
    RefinementCase<dim>::cut_y,
    RefinementCase<dim>::cut_z
  };

  // order the direction of lines
  // 0->x, 1->y, 2->z
  static const unsigned int direction[lines_per_cell]=
  {1,1,0,0,1,1,0,0,2,2,2,2};

  return ((cell_refinement_case & cut_one[direction[line_no]]) ?
          RefinementCase<1>::cut_x : RefinementCase<1>::no_refinement);
}



template<>
RefinementCase<1>
GeometryInfo<1>::min_cell_refinement_case_for_face_refinement(const RefinementCase<0> &,
    const unsigned int,
    const bool,
    const bool,
    const bool)
{
  const unsigned int dim = 1;
  Assert(false, ExcImpossibleInDim(dim));

  return RefinementCase<dim>::no_refinement;
}


template<>
RefinementCase<2>
GeometryInfo<2>::min_cell_refinement_case_for_face_refinement(const RefinementCase<1> &face_refinement_case,
    const unsigned int face_no,
    const bool,
    const bool,
    const bool)
{
  const unsigned int dim = 2;
  Assert(face_refinement_case<RefinementCase<dim-1>::isotropic_refinement+1,
         ExcIndexRange(face_refinement_case, 0, RefinementCase<dim-1>::isotropic_refinement+1));
  Assert(face_no<GeometryInfo<dim>::faces_per_cell,
         ExcIndexRange(face_no, 0, GeometryInfo<dim>::faces_per_cell));

  if (face_refinement_case==RefinementCase<dim>::cut_x)
    return (face_no/2) ? RefinementCase<dim>::cut_x : RefinementCase<dim>::cut_y;
  else
    return RefinementCase<dim>::no_refinement;
}


template<>
RefinementCase<3>
GeometryInfo<3>::min_cell_refinement_case_for_face_refinement(const RefinementCase<2> &face_refinement_case,
    const unsigned int face_no,
    const bool face_orientation,
    const bool /*face_flip*/,
    const bool face_rotation)
{
  const unsigned int dim=3;
  Assert(face_refinement_case<RefinementCase<dim-1>::isotropic_refinement+1,
         ExcIndexRange(face_refinement_case, 0, RefinementCase<dim-1>::isotropic_refinement+1));
  Assert(face_no<GeometryInfo<dim>::faces_per_cell,
         ExcIndexRange(face_no, 0, GeometryInfo<dim>::faces_per_cell));

  static const RefinementCase<2> flip[4]=
  {
    RefinementCase<2>::no_refinement,
    RefinementCase<2>::cut_y,
    RefinementCase<2>::cut_x,
    RefinementCase<2>::cut_xy
  };

  // correct the face_refinement_case for
  // face_orientation and face_rotation. for
  // face_orientation, 'true' is the default
  // value whereas for face_rotation, 'false'
  // is standard. If
  // <tt>face_rotation==face_orientation</tt>,
  // then one of them is non-standard and we
  // have to swap cut_x and cut_y, otherwise no
  // change is necessary.  face_flip has no
  // influence. however, in order to keep the
  // interface consistent with other functions,
  // we still include it as an argument to this
  // function
  const RefinementCase<dim-1> std_face_ref = (face_orientation==face_rotation) ? flip[face_refinement_case] : face_refinement_case;

  static const RefinementCase<dim> face_to_cell[3][4]=
  {
    {
      RefinementCase<dim>::no_refinement,  // faces 0 and 1
      RefinementCase<dim>::cut_y,          // cut_x in face 0 means cut_y for the cell
      RefinementCase<dim>::cut_z,
      RefinementCase<dim>::cut_yz
    },

    {
      RefinementCase<dim>::no_refinement,  // faces 2 and 3 (note that x and y are "exchanged on faces 2 and 3")
      RefinementCase<dim>::cut_z,
      RefinementCase<dim>::cut_x,
      RefinementCase<dim>::cut_xz
    },

    {
      RefinementCase<dim>::no_refinement,  // faces 4 and 5
      RefinementCase<dim>::cut_x,
      RefinementCase<dim>::cut_y,
      RefinementCase<dim>::cut_xy
    }
  };

  return face_to_cell[face_no/2][std_face_ref];
}



template<>
RefinementCase<1>
GeometryInfo<1>::min_cell_refinement_case_for_line_refinement(const unsigned int line_no)
{
  Assert(line_no==0, ExcIndexRange(line_no,0,1));

  return RefinementCase<1>::cut_x;
}


template<>
RefinementCase<2>
GeometryInfo<2>::min_cell_refinement_case_for_line_refinement(const unsigned int line_no)
{
  const unsigned int dim = 2;
  Assert(line_no<GeometryInfo<dim>::lines_per_cell,
         ExcIndexRange(line_no, 0, GeometryInfo<dim>::lines_per_cell));

  return (line_no/2) ? RefinementCase<2>::cut_x : RefinementCase<2>::cut_y;
}


template<>
RefinementCase<3>
GeometryInfo<3>::min_cell_refinement_case_for_line_refinement(const unsigned int line_no)
{
  const unsigned int dim=3;
  Assert(line_no<GeometryInfo<dim>::lines_per_cell,
         ExcIndexRange(line_no, 0, GeometryInfo<dim>::lines_per_cell));

  static const RefinementCase<dim> ref_cases[6]=
  {
    RefinementCase<dim>::cut_y,  // lines  0 and  1
    RefinementCase<dim>::cut_x,  // lines  2 and  3
    RefinementCase<dim>::cut_y,  // lines  4 and  5
    RefinementCase<dim>::cut_x,  // lines  6 and  7
    RefinementCase<dim>::cut_z,  // lines  8 and  9
    RefinementCase<dim>::cut_z
  }; // lines 10 and 11

  return ref_cases[line_no/2];
}



template <>
unsigned int
GeometryInfo<3>::standard_to_real_face_vertex(const unsigned int vertex,
                                              const bool face_orientation,
                                              const bool face_flip,
                                              const bool face_rotation)
{
  Assert(vertex<GeometryInfo<3>::vertices_per_face,
         ExcIndexRange(vertex,0,GeometryInfo<3>::vertices_per_face));

  // set up a table to make sure that
  // we handle non-standard faces correctly
  //
  // so set up a table that for each vertex (of
  // a quad in standard position) describes
  // which vertex to take
  //
  // first index: four vertices 0...3
  //
  // second index: face_orientation; 0:
  // opposite normal, 1: standard
  //
  // third index: face_flip; 0: standard, 1:
  // face rotated by 180 degrees
  //
  // forth index: face_rotation: 0: standard,
  // 1: face rotated by 90 degrees

  static const unsigned int vertex_translation[4][2][2][2] =
  {
    { { { 0, 2 },  // vertex 0, face_orientation=false, face_flip=false, face_rotation=false and true
        { 3, 1 }
      }, // vertex 0, face_orientation=false, face_flip=true, face_rotation=false and true
      { { 0, 2 },  // vertex 0, face_orientation=true, face_flip=false, face_rotation=false and true
        { 3, 1 }
      }
    },// vertex 0, face_orientation=true, face_flip=true, face_rotation=false and true

    { { { 2, 3 },  // vertex 1 ...
        { 1, 0 }
      },
      { { 1, 0 },
        { 2, 3 }
      }
    },

    { { { 1, 0 },  // vertex 2 ...
        { 2, 3 }
      },
      { { 2, 3 },
        { 1, 0 }
      }
    },

    { { { 3, 1 },  // vertex 3 ...
        { 0, 2 }
      },
      { { 3, 1 },
        { 0, 2 }
      }
    }
  };

  return vertex_translation[vertex][face_orientation][face_flip][face_rotation];
}



template <int dim>
unsigned int
GeometryInfo<dim>::standard_to_real_face_vertex(const unsigned int vertex,
                                                const bool,
                                                const bool,
                                                const bool)
{
  Assert(dim>1, ExcImpossibleInDim(dim));
  Assert(vertex<GeometryInfo<dim>::vertices_per_face,
         ExcIndexRange(vertex,0,GeometryInfo<dim>::vertices_per_face));
  return vertex;
}



template <>
unsigned int
GeometryInfo<3>::real_to_standard_face_vertex(const unsigned int vertex,
                                              const bool face_orientation,
                                              const bool face_flip,
                                              const bool face_rotation)
{
  Assert(vertex<GeometryInfo<3>::vertices_per_face,
         ExcIndexRange(vertex,0,GeometryInfo<3>::vertices_per_face));

  // set up a table to make sure that
  // we handle non-standard faces correctly
  //
  // so set up a table that for each vertex (of
  // a quad in standard position) describes
  // which vertex to take
  //
  // first index: four vertices 0...3
  //
  // second index: face_orientation; 0:
  // opposite normal, 1: standard
  //
  // third index: face_flip; 0: standard, 1:
  // face rotated by 180 degrees
  //
  // forth index: face_rotation: 0: standard,
  // 1: face rotated by 90 degrees

  static const unsigned int vertex_translation[4][2][2][2] =
  {
    { { { 0, 2 },  // vertex 0, face_orientation=false, face_flip=false, face_rotation=false and true
        { 3, 1 }
      }, // vertex 0, face_orientation=false, face_flip=true, face_rotation=false and true
      { { 0, 1 },  // vertex 0, face_orientation=true, face_flip=false, face_rotation=false and true
        { 3, 2 }
      }
    },// vertex 0, face_orientation=true, face_flip=true, face_rotation=false and true

    { { { 2, 3 },  // vertex 1 ...
        { 1, 0 }
      },
      { { 1, 3 },
        { 2, 0 }
      }
    },

    { { { 1, 0 },  // vertex 2 ...
        { 2, 3 }
      },
      { { 2, 0 },
        { 1, 3 }
      }
    },

    { { { 3, 1 },  // vertex 3 ...
        { 0, 2 }
      },
      { { 3, 2 },
        { 0, 1 }
      }
    }
  };

  return vertex_translation[vertex][face_orientation][face_flip][face_rotation];
}



template <int dim>
unsigned int
GeometryInfo<dim>::real_to_standard_face_vertex(const unsigned int vertex,
                                                const bool,
                                                const bool,
                                                const bool)
{
  Assert(dim>1, ExcImpossibleInDim(dim));
  Assert(vertex<GeometryInfo<dim>::vertices_per_face,
         ExcIndexRange(vertex,0,GeometryInfo<dim>::vertices_per_face));
  return vertex;
}



template <>
unsigned int
GeometryInfo<3>::standard_to_real_face_line(const unsigned int line,
                                            const bool face_orientation,
                                            const bool face_flip,
                                            const bool face_rotation)
{
  Assert(line<GeometryInfo<3>::lines_per_face,
         ExcIndexRange(line,0,GeometryInfo<3>::lines_per_face));


  // make sure we handle
  // non-standard faces correctly
  //
  // so set up a table that for each line (of a
  // quad) describes which line to take
  //
  // first index: four lines 0...3
  //
  // second index: face_orientation; 0:
  // opposite normal, 1: standard
  //
  // third index: face_flip; 0: standard, 1:
  // face rotated by 180 degrees
  //
  // forth index: face_rotation: 0: standard,
  // 1: face rotated by 90 degrees

  static const unsigned int line_translation[4][2][2][2] =
  {
    { { { 2, 0 },  // line 0, face_orientation=false, face_flip=false, face_rotation=false and true
        { 3, 1 }
      }, // line 0, face_orientation=false, face_flip=true, face_rotation=false and true
      { { 0, 3 },  // line 0, face_orientation=true, face_flip=false, face_rotation=false and true
        { 1, 2 }
      }
    },// line 0, face_orientation=true, face_flip=true, face_rotation=false and true

    { { { 3, 1 },  // line 1 ...
        { 2, 0 }
      },
      { { 1, 2 },
        { 0, 3 }
      }
    },

    { { { 0, 3 },  // line 2 ...
        { 1, 2 }
      },
      { { 2, 0 },
        { 3, 1 }
      }
    },

    { { { 1, 2 },  // line 3 ...
        { 0, 3 }
      },
      { { 3, 1 },
        { 2, 0 }
      }
    }
  };

  return line_translation[line][face_orientation][face_flip][face_rotation];
}



template <int dim>
unsigned int
GeometryInfo<dim>::standard_to_real_face_line(const unsigned int line,
                                              const bool,
                                              const bool,
                                              const bool)
{
  Assert(false, ExcNotImplemented());
  return line;
}



template <>
unsigned int
GeometryInfo<3>::real_to_standard_face_line(const unsigned int line,
                                            const bool face_orientation,
                                            const bool face_flip,
                                            const bool face_rotation)
{
  Assert(line<GeometryInfo<3>::lines_per_face,
         ExcIndexRange(line,0,GeometryInfo<3>::lines_per_face));


  // make sure we handle
  // non-standard faces correctly
  //
  // so set up a table that for each line (of a
  // quad) describes which line to take
  //
  // first index: four lines 0...3
  //
  // second index: face_orientation; 0:
  // opposite normal, 1: standard
  //
  // third index: face_flip; 0: standard, 1:
  // face rotated by 180 degrees
  //
  // forth index: face_rotation: 0: standard,
  // 1: face rotated by 90 degrees

  static const unsigned int line_translation[4][2][2][2] =
  {
    { { { 2, 0 },  // line 0, face_orientation=false, face_flip=false, face_rotation=false and true
        { 3, 1 }
      }, // line 0, face_orientation=false, face_flip=true, face_rotation=false and true
      { { 0, 2 },  // line 0, face_orientation=true, face_flip=false, face_rotation=false and true
        { 1, 3 }
      }
    },// line 0, face_orientation=true, face_flip=true, face_rotation=false and true

    { { { 3, 1 },  // line 1 ...
        { 2, 0 }
      },
      { { 1, 3 },
        { 0, 2 }
      }
    },

    { { { 0, 3 },  // line 2 ...
        { 1, 2 }
      },
      { { 2, 1 },
        { 3, 0 }
      }
    },

    { { { 1, 2 },  // line 3 ...
        { 0, 3 }
      },
      { { 3, 0 },
        { 2, 1 }
      }
    }
  };

  return line_translation[line][face_orientation][face_flip][face_rotation];
}



template <int dim>
unsigned int
GeometryInfo<dim>::real_to_standard_face_line(const unsigned int line,
                                              const bool,
                                              const bool,
                                              const bool)
{
  Assert(false, ExcNotImplemented());
  return line;
}



template <>
unsigned int
GeometryInfo<1>::child_cell_on_face (const RefinementCase<1> &,
                                     const unsigned int face,
                                     const unsigned int subface,
                                     const bool, const bool, const bool,
                                     const RefinementCase<0> &)
{
  Assert (face<faces_per_cell, ExcIndexRange(face, 0, faces_per_cell));
  Assert (subface<max_children_per_face,
          ExcIndexRange(subface, 0, max_children_per_face));

  return face;
}



template <>
unsigned int
GeometryInfo<2>::child_cell_on_face (const RefinementCase<2> &ref_case,
                                     const unsigned int face,
                                     const unsigned int subface,
                                     const bool /*face_orientation*/,
                                     const bool face_flip,
                                     const bool /*face_rotation*/,
                                     const RefinementCase<1> &)
{
  Assert (face<faces_per_cell, ExcIndexRange(face, 0, faces_per_cell));
  Assert (subface<max_children_per_face,
          ExcIndexRange(subface, 0, max_children_per_face));

  // always return the child adjacent to the specified
  // subface. if the face of a cell is not refined, don't
  // throw an assertion but deliver the child adjacent to
  // the face nevertheless, i.e. deliver the child of
  // this cell adjacent to the subface of a possibly
  // refined neighbor. this simplifies setting neighbor
  // information in execute_refinement.
  static const unsigned int
  subcells[2][RefinementCase<2>::isotropic_refinement][faces_per_cell][max_children_per_face] =
  {
    {
      // Normal orientation (face_filp = false)
      {{0,0},{1,1},{0,1},{0,1}},          // cut_x
      {{0,1},{0,1},{0,0},{1,1}},          // cut_y
      {{0,2},{1,3},{0,1},{2,3}}           // cut_z
    },
    {
      // Flipped orientation (face_flip = true)
      {{0,0},{1,1},{1,0},{1,0}},          // cut_x
      {{1,0},{1,0},{0,0},{1,1}},          // cut_y
      {{2,0},{3,1},{1,0},{3,2}}           // cut_z
    }
  };

  return subcells[face_flip][ref_case-1][face][subface];
}



template <>
unsigned int
GeometryInfo<3>::child_cell_on_face (const RefinementCase<3> &ref_case,
                                     const unsigned int face,
                                     const unsigned int subface,
                                     const bool face_orientation,
                                     const bool face_flip,
                                     const bool face_rotation,
                                     const RefinementCase<2> &face_ref_case)
{
  const unsigned int dim = 3;

  Assert (ref_case>RefinementCase<dim-1>::no_refinement, ExcMessage("Cell has no children."));
  Assert (face<faces_per_cell, ExcIndexRange(face, 0, faces_per_cell));
  Assert (subface<GeometryInfo<dim-1>::n_children(face_ref_case) ||
          (subface==0 && face_ref_case==RefinementCase<dim-1>::no_refinement),
          ExcIndexRange(subface, 0, GeometryInfo<2>::n_children(face_ref_case)));

  // invalid number used for invalid cases,
  // e.g. when the children are more refined at
  // a given face than the face itself
  static const unsigned int e=invalid_unsigned_int;

  // the whole process of finding a child cell
  // at a given subface considering the
  // possibly anisotropic refinement cases of
  // the cell and the face as well as
  // orientation, flip and rotation of the face
  // is quite complicated. thus, we break it
  // down into several steps.

  // first step: convert the given face refine
  // case to a face refine case concerning the
  // face in standard orientation (, flip and
  // rotation). This only affects cut_x and
  // cut_y
  static const RefinementCase<dim-1> flip[4]=
  {
    RefinementCase<dim-1>::no_refinement,
    RefinementCase<dim-1>::cut_y,
    RefinementCase<dim-1>::cut_x,
    RefinementCase<dim-1>::cut_xy
  };
  // for face_orientation, 'true' is the
  // default value whereas for face_rotation,
  // 'false' is standard. If
  // <tt>face_rotation==face_orientation</tt>,
  // then one of them is non-standard and we
  // have to swap cut_x and cut_y, otherwise no
  // change is necessary.
  const RefinementCase<dim-1> std_face_ref = (face_orientation==face_rotation) ? flip[face_ref_case] : face_ref_case;

  // second step: convert the given subface
  // index to the one for a standard face
  // respecting face_orientation, face_flip and
  // face_rotation

  // first index:  face_ref_case
  // second index: face_orientation
  // third index:  face_flip
  // forth index:  face_rotation
  // fifth index:  subface index
  static const unsigned int subface_exchange[4][2][2][2][4]=
  {
    // no_refinement (subface 0 stays 0,
    // all others are invalid)
    { { { {0,e,e,e},
          {0,e,e,e}
        },
        { {0,e,e,e},
          {0,e,e,e}
        }
      },
      { { {0,e,e,e},
          {0,e,e,e}
        },
        { {0,e,e,e},
          {0,e,e,e}
        }
      }
    },
    // cut_x (here, if the face is only
    // rotated OR only falsely oriented,
    // then subface 0 of the non-standard
    // face does NOT correspond to one of
    // the subfaces of a standard
    // face. Thus we indicate the subface
    // which is located at the lower left
    // corner (the origin of the face's
    // local coordinate system) with
    // '0'. The rest of this issue is
    // taken care of using the above
    // conversion to a 'standard face
    // refine case')
    { { { {0,1,e,e},
          {0,1,e,e}
        },
        { {1,0,e,e},
          {1,0,e,e}
        }
      },
      { { {0,1,e,e},
          {0,1,e,e}
        },
        { {1,0,e,e},
          {1,0,e,e}
        }
      }
    },
    // cut_y (the same applies as for
    // cut_x)
    { { { {0,1,e,e},
          {1,0,e,e}
        },
        { {1,0,e,e},
          {0,1,e,e}
        }
      },
      { { {0,1,e,e},
          {1,0,e,e}
        },
        { {1,0,e,e},
          {0,1,e,e}
        }
      }
    },
    // cut_xyz: this information is
    // identical to the information
    // returned by
    // GeometryInfo<3>::real_to_standard_face_vertex()
    { { { {0,2,1,3},    // face_orientation=false, face_flip=false, face_rotation=false, subfaces 0,1,2,3
          {2,3,0,1}
        },   // face_orientation=false, face_flip=false, face_rotation=true,  subfaces 0,1,2,3
        { {3,1,2,0},    // face_orientation=false, face_flip=true,  face_rotation=false, subfaces 0,1,2,3
          {1,0,3,2}
        }
      },  // face_orientation=false, face_flip=true,  face_rotation=true,  subfaces 0,1,2,3
      { { {0,1,2,3},    // face_orientation=true,  face_flip=false, face_rotation=false, subfaces 0,1,2,3
          {1,3,0,2}
        },   // face_orientation=true,  face_flip=false, face_rotation=true,  subfaces 0,1,2,3
        { {3,2,1,0},    // face_orientation=true,  face_flip=true,  face_rotation=false, subfaces 0,1,2,3
          {2,0,3,1}
        }
      }
    }
  };// face_orientation=true,  face_flip=true,  face_rotation=true,  subfaces 0,1,2,3

  const unsigned int std_subface=subface_exchange
                                 [face_ref_case]
                                 [face_orientation]
                                 [face_flip]
                                 [face_rotation]
                                 [subface];
  Assert (std_subface!=e, ExcInternalError());

  // third step: these are the children, which
  // can be found at the given subfaces of an
  // isotropically refined (standard) face
  //
  // first index:  (refinement_case-1)
  // second index: face_index
  // third index:  subface_index (isotropic refinement)
  static const unsigned int
  iso_children[RefinementCase<dim>::cut_xyz][faces_per_cell][max_children_per_face] =
  {
    // cut_x
    { {0, 0, 0, 0},  // face 0, subfaces 0,1,2,3
      {1, 1, 1, 1},  // face 1, subfaces 0,1,2,3
      {0, 0, 1, 1},  // face 2, subfaces 0,1,2,3
      {0, 0, 1, 1},  // face 3, subfaces 0,1,2,3
      {0, 1, 0, 1},  // face 4, subfaces 0,1,2,3
      {0, 1, 0, 1}
    }, // face 5, subfaces 0,1,2,3
    // cut_y
    { {0, 1, 0, 1},
      {0, 1, 0, 1},
      {0, 0, 0, 0},
      {1, 1, 1, 1},
      {0, 0, 1, 1},
      {0, 0, 1, 1}
    },
    // cut_xy
    { {0, 2, 0, 2},
      {1, 3, 1, 3},
      {0, 0, 1, 1},
      {2, 2, 3, 3},
      {0, 1, 2, 3},
      {0, 1, 2, 3}
    },
    // cut_z
    { {0, 0, 1, 1},
      {0, 0, 1, 1},
      {0, 1, 0, 1},
      {0, 1, 0, 1},
      {0, 0, 0, 0},
      {1, 1, 1, 1}
    },
    // cut_xz
    { {0, 0, 1, 1},
      {2, 2, 3, 3},
      {0, 1, 2, 3},
      {0, 1, 2, 3},
      {0, 2, 0, 2},
      {1, 3, 1, 3}
    },
    // cut_yz
    { {0, 1, 2, 3},
      {0, 1, 2, 3},
      {0, 2, 0, 2},
      {1, 3, 1, 3},
      {0, 0, 1, 1},
      {2, 2, 3, 3}
    },
    // cut_xyz
    { {0, 2, 4, 6},
      {1, 3, 5, 7},
      {0, 4, 1, 5},
      {2, 6, 3, 7},
      {0, 1, 2, 3},
      {4, 5, 6, 7}
    }
  };

  // forth step: check, whether the given face
  // refine case is valid for the given cell
  // refine case. this is the case, if the
  // given face refine case is at least as
  // refined as the face is for the given cell
  // refine case

  // note, that we are considering standard
  // face refinement cases here and thus must
  // not pass the given orientation, flip and
  // rotation flags
  if ((std_face_ref & face_refinement_case(ref_case, face))
      == face_refinement_case(ref_case, face))
    {
      // all is fine. for anisotropic face
      // refine cases, select one of the
      // isotropic subfaces which neighbors the
      // same child

      // first index: (standard) face refine case
      // second index: subface index
      static const unsigned int equivalent_iso_subface[4][4]=
      {
        {0,e,e,e},                    // no_refinement
        {0,3,e,e},                    // cut_x
        {0,3,e,e},                    // cut_y
        {0,1,2,3}
      };                   // cut_xy

      const unsigned int equ_std_subface
        =equivalent_iso_subface[std_face_ref][std_subface];
      Assert (equ_std_subface!=e, ExcInternalError());

      return iso_children[ref_case-1][face][equ_std_subface];
    }
  else
    {
      // the face_ref_case was too coarse,
      // throw an error
      Assert(false,
             ExcMessage("The face RefineCase is too coarse "
                        "for the given cell RefineCase."));
    }
  // we only get here in case of an error
  return e;
}



template <>
unsigned int
GeometryInfo<4>::child_cell_on_face (const RefinementCase<4> &,
                                     const unsigned int,
                                     const unsigned int,
                                     const bool, const bool, const bool,
                                     const RefinementCase<3> &)
{
  Assert(false, ExcNotImplemented());
  return invalid_unsigned_int;
}



template <>
unsigned int
GeometryInfo<1>::line_to_cell_vertices (const unsigned int line,
                                        const unsigned int vertex)
{
  Assert (line<lines_per_cell, ExcIndexRange(line, 0, lines_per_cell));
  Assert (vertex<2, ExcIndexRange(vertex, 0, 2));

  return vertex;
}


template <>
unsigned int
GeometryInfo<2>::line_to_cell_vertices (const unsigned int line,
                                        const unsigned int vertex)
{
  return child_cell_on_face(RefinementCase<2>::isotropic_refinement, line, vertex);
}



template <>
unsigned int
GeometryInfo<3>::line_to_cell_vertices (const unsigned int line,
                                        const unsigned int vertex)
{
  Assert (line<lines_per_cell, ExcIndexRange(line, 0, lines_per_cell));
  Assert (vertex<2, ExcIndexRange(vertex, 0, 2));

  static const unsigned
  vertices[lines_per_cell][2] = {{0, 2},  // bottom face
    {1, 3},
    {0, 1},
    {2, 3},
    {4, 6},  // top face
    {5, 7},
    {4, 5},
    {6, 7},
    {0, 4},  // connects of bottom
    {1, 5},  //   top face
    {2, 6},
    {3, 7}
  };
  return vertices[line][vertex];
}


template <>
unsigned int
GeometryInfo<4>::line_to_cell_vertices (const unsigned int,
                                        const unsigned int)
{
  Assert(false, ExcNotImplemented());
  return invalid_unsigned_int;
}


template <>
unsigned int
GeometryInfo<1>::face_to_cell_lines (const unsigned int face,
                                     const unsigned int line,
                                     const bool, const bool, const bool)
{
  Assert (face+1<faces_per_cell+1, ExcIndexRange(face, 0, faces_per_cell));
  Assert (line+1<lines_per_face+1, ExcIndexRange(line, 0, lines_per_face));

  // There is only a single line, so
  // it must be this.
  return 0;
}



template <>
unsigned int
GeometryInfo<2>::face_to_cell_lines (const unsigned int face,
                                     const unsigned int line,
                                     const bool, const bool, const bool)
{
  Assert (face<faces_per_cell, ExcIndexRange(face, 0, faces_per_cell));
  Assert (line<lines_per_face, ExcIndexRange(line, 0, lines_per_face));

  // The face is a line itself.
  return face;
}



template <>
unsigned int
GeometryInfo<3>::face_to_cell_lines (const unsigned int face,
                                     const unsigned int line,
                                     const bool face_orientation,
                                     const bool face_flip,
                                     const bool face_rotation)
{
  Assert (face<faces_per_cell, ExcIndexRange(face, 0, faces_per_cell));
  Assert (line<lines_per_face, ExcIndexRange(line, 0, lines_per_face));

  static const unsigned
  lines[faces_per_cell][lines_per_face] = {{8,10, 0, 4}, // left face
    {9,11, 1, 5}, // right face
    {2, 6, 8, 9}, // front face
    {3, 7,10,11}, // back face
    {0, 1, 2, 3}, // bottom face
    {4, 5, 6, 7}
  };// top face
  return lines[face][real_to_standard_face_line(line,
                                                face_orientation,
                                                face_flip,
                                                face_rotation)];
}



template<int dim>
unsigned int
GeometryInfo<dim>::face_to_cell_lines (const unsigned int,
                                       const unsigned int,
                                       const bool, const bool, const bool)
{
  Assert(false, ExcNotImplemented());
  return invalid_unsigned_int;
}



template <int dim>
unsigned int
GeometryInfo<dim>::face_to_cell_vertices (const unsigned int face,
                                          const unsigned int vertex,
                                          const bool face_orientation,
                                          const bool face_flip,
                                          const bool face_rotation)
{
  return child_cell_on_face(RefinementCase<dim>::isotropic_refinement, face, vertex,
                            face_orientation, face_flip, face_rotation);
}



template <int dim>
Point<dim>
GeometryInfo<dim>::project_to_unit_cell (const Point<dim> &q)
{
  Point<dim> p = q;
  for (unsigned int i=0; i<dim; i++)
    if      (p[i] < 0.)  p[i] = 0.;
    else if (p[i] > 1.)  p[i] = 1.;

  return p;
}



template <int dim>
double
GeometryInfo<dim>::distance_to_unit_cell (const Point<dim> &p)
{
  double result = 0.0;

  for (unsigned int i=0; i<dim; i++)
    if ((-p[i]) > result)
      result = -p[i];
    else if ((p[i]-1.) > result)
      result = (p[i] - 1.);

  return result;
}



template <int dim>
double
GeometryInfo<dim>::
d_linear_shape_function (const Point<dim> &xi,
                         const unsigned int i)
{
  Assert (i < GeometryInfo<dim>::vertices_per_cell,
          ExcIndexRange (i, 0, GeometryInfo<dim>::vertices_per_cell));

  switch (dim)
    {
    case 1:
    {
      const double x = xi[0];
      switch (i)
        {
        case 0:
          return 1-x;
        case 1:
          return x;
        }
    }

    case 2:
    {
      const double x = xi[0];
      const double y = xi[1];
      switch (i)
        {
        case 0:
          return (1-x)*(1-y);
        case 1:
          return x*(1-y);
        case 2:
          return (1-x)*y;
        case 3:
          return x*y;
        }
    }

    case 3:
    {
      const double x = xi[0];
      const double y = xi[1];
      const double z = xi[2];
      switch (i)
        {
        case 0:
          return (1-x)*(1-y)*(1-z);
        case 1:
          return x*(1-y)*(1-z);
        case 2:
          return (1-x)*y*(1-z);
        case 3:
          return x*y*(1-z);
        case 4:
          return (1-x)*(1-y)*z;
        case 5:
          return x*(1-y)*z;
        case 6:
          return (1-x)*y*z;
        case 7:
          return x*y*z;
        }
    }

    default:
      Assert (false, ExcNotImplemented());
    }
  return -1e9;
}



template <>
Tensor<1,1>
GeometryInfo<1>::
d_linear_shape_function_gradient (const Point<1> &,
                                  const unsigned int i)
{
  Assert (i < GeometryInfo<1>::vertices_per_cell,
          ExcIndexRange (i, 0, GeometryInfo<1>::vertices_per_cell));

  switch (i)
    {
    case 0:
      return Point<1>(-1.);
    case 1:
      return Point<1>(1.);
    }

  return Point<1>(-1e9);
}



template <>
Tensor<1,2>
GeometryInfo<2>::
d_linear_shape_function_gradient (const Point<2> &xi,
                                  const unsigned int i)
{
  Assert (i < GeometryInfo<2>::vertices_per_cell,
          ExcIndexRange (i, 0, GeometryInfo<2>::vertices_per_cell));

  const double x = xi[0];
  const double y = xi[1];
  switch (i)
    {
    case 0:
      return Point<2>(-(1-y),-(1-x));
    case 1:
      return Point<2>(1-y,-x);
    case 2:
      return Point<2>(-y, 1-x);
    case 3:
      return Point<2>(y,x);
    }
  return Point<2> (-1e9, -1e9);
}



template <>
Tensor<1,3>
GeometryInfo<3>::
d_linear_shape_function_gradient (const Point<3> &xi,
                                  const unsigned int i)
{
  Assert (i < GeometryInfo<3>::vertices_per_cell,
          ExcIndexRange (i, 0, GeometryInfo<3>::vertices_per_cell));

  const double x = xi[0];
  const double y = xi[1];
  const double z = xi[2];
  switch (i)
    {
    case 0:
      return Point<3>(-(1-y)*(1-z),
                      -(1-x)*(1-z),
                      -(1-x)*(1-y));
    case 1:
      return Point<3>((1-y)*(1-z),
                      -x*(1-z),
                      -x*(1-y));
    case 2:
      return Point<3>(-y*(1-z),
                      (1-x)*(1-z),
                      -(1-x)*y);
    case 3:
      return Point<3>(y*(1-z),
                      x*(1-z),
                      -x*y);
    case 4:
      return Point<3>(-(1-y)*z,
                      -(1-x)*z,
                      (1-x)*(1-y));
    case 5:
      return Point<3>((1-y)*z,
                      -x*z,
                      x*(1-y));
    case 6:
      return Point<3>(-y*z,
                      (1-x)*z,
                      (1-x)*y);
    case 7:
      return Point<3>(y*z, x*z, x*y);
    }

  return Point<3> (-1e9, -1e9, -1e9);
}



template <int dim>
Tensor<1,dim>
GeometryInfo<dim>::
d_linear_shape_function_gradient (const Point<dim> &,
                                  const unsigned int)
{
  Assert (false, ExcNotImplemented());
  return Tensor<1,dim>();
}





namespace internal
{
  namespace GeometryInfo
  {
    // wedge product of a single
    // vector in 2d: we just have to
    // rotate it by 90 degrees to the
    // right
    inline
    Tensor<1,2>
    wedge_product (const Tensor<1,2> (&derivative)[1])
    {
      Tensor<1,2> result;
      result[0] = derivative[0][1];
      result[1] = -derivative[0][0];

      return result;
    }


    // wedge product of 2 vectors in
    // 3d is the cross product
    inline
    Tensor<1,3>
    wedge_product (const Tensor<1,3> (&derivative)[2])
    {
      Tensor<1,3> result;
      cross_product (result, derivative[0], derivative[1]);

      return result;
    }


    // wedge product of dim vectors
    // in dim-d: that's the
    // determinant of the matrix
    template <int dim>
    inline
    Tensor<0,dim>
    wedge_product (const Tensor<1,dim> (&derivative)[dim])
    {
      Tensor<2,dim> jacobian;
      for (unsigned int i=0; i<dim; ++i)
        jacobian[i] = derivative[i];

      return determinant (jacobian);
    }
  }
}


template <int dim>
template <int spacedim>
void
GeometryInfo<dim>::
alternating_form_at_vertices
#ifndef DEAL_II_CONSTEXPR_BUG
(const Point<spacedim> (&vertices)[vertices_per_cell],
 Tensor<spacedim-dim,spacedim> (&forms)[vertices_per_cell])
#else
(const Point<spacedim> *vertices,
 Tensor<spacedim-dim,spacedim> *forms)
#endif
{
  // for each of the vertices,
  // compute the alternating form
  // of the mapped unit
  // vectors. consider for
  // example the case of a quad
  // in spacedim==3: to do so, we
  // need to see how the
  // infinitesimal vectors
  // (d\xi_1,0) and (0,d\xi_2)
  // are transformed into
  // spacedim-dimensional space
  // and then form their cross
  // product (i.e. the wedge product
  // of two vectors). to this end, note
  // that
  //    \vec x = sum_i \vec v_i phi_i(\vec xi)
  // so the transformed vectors are
  //   [x(\xi+(d\xi_1,0))-x(\xi)]/d\xi_1
  // and
  //   [x(\xi+(0,d\xi_2))-x(\xi)]/d\xi_2
  // which boils down to the columns
  // of the 3x2 matrix \grad_\xi x(\xi)
  //
  // a similar reasoning would
  // hold for all dim,spacedim
  // pairs -- we only have to
  // compute the wedge product of
  // the columns of the
  // derivatives
  for (unsigned int i=0; i<vertices_per_cell; ++i)
    {
      Tensor<1,spacedim> derivatives[dim];

      for (unsigned int j=0; j<vertices_per_cell; ++j)
        {
          const Tensor<1,dim> grad_phi_j
            = d_linear_shape_function_gradient (unit_cell_vertex(i),
                                                j);
          for (unsigned int l=0; l<dim; ++l)
            derivatives[l] += vertices[j] * grad_phi_j[l];
        }

      forms[i] = internal::GeometryInfo::wedge_product (derivatives);
    }
}


template struct GeometryInfo<1>;
template struct GeometryInfo<2>;
template struct GeometryInfo<3>;
template struct GeometryInfo<4>;

template
void
GeometryInfo<1>::
alternating_form_at_vertices
#ifndef DEAL_II_CONSTEXPR_BUG
(const Point<1> (&)[vertices_per_cell],
 Tensor<1-1,1> (&)[vertices_per_cell])
#else
(const Point<1> *, Tensor<1-1,1> *)
#endif
;

template
void
GeometryInfo<1>::
alternating_form_at_vertices
#ifndef DEAL_II_CONSTEXPR_BUG
(const Point<2> (&)[vertices_per_cell],
 Tensor<2-1,2> (&)[vertices_per_cell])
#else
(const Point<2> *, Tensor<2-1,2> *)
#endif
;

template
void
GeometryInfo<2>::
alternating_form_at_vertices
#ifndef DEAL_II_CONSTEXPR_BUG
(const Point<2> (&vertices)[vertices_per_cell],
 Tensor<2-2,2> (&forms)[vertices_per_cell])
#else
(const Point<2> *, Tensor<2-2,2> *)
#endif
;

template
void
GeometryInfo<2>::
alternating_form_at_vertices
#ifndef DEAL_II_CONSTEXPR_BUG
(const Point<3> (&vertices)[vertices_per_cell],
 Tensor<3-2,3> (&forms)[vertices_per_cell])
#else
(const Point<3> *, Tensor<3-2,3> *)
#endif
;


template
void
GeometryInfo<3>::
alternating_form_at_vertices
#ifndef DEAL_II_CONSTEXPR_BUG
(const Point<3> (&vertices)[vertices_per_cell],
 Tensor<3-3,3> (&forms)[vertices_per_cell])
#else
(const Point<3> *, Tensor<3-3,3> *)
#endif
;


DEAL_II_NAMESPACE_CLOSE
