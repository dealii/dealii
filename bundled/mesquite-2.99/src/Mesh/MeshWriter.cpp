/* ***************************************************************** 
    MESQUITE -- The Mesh Quality Improvement Toolkit

    Copyright 2005 Lawrence Livermore National Laboratory.  Under 
    the terms of Contract B545069 with the University of Wisconsin -- 
    Madison, Lawrence Livermore National Laboratory retains certain
    rights in this software.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License 
    (lgpl.txt) along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

    kraftche@cae.wisc.edu    

  ***************************************************************** */

#ifndef MSQ_MESH_WRITER_CPP
#define MSQ_MESH_WRITER_CPP

#include "MeshWriter.hpp"
#include "Mesquite.hpp"
#include "MeshImpl.hpp"
#include "MsqError.hpp"
#include "PatchData.hpp"
#include "PlanarDomain.hpp"
#include "VtkTypeInfo.hpp"
#include "EdgeIterator.hpp"

#include <memory>
#include <limits>
#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include <iomanip>
using namespace std;

#include <stdio.h>

namespace MESQUITE_NS {

namespace MeshWriter {

/**\brief Transform from coordinates in the XY-plane 
 *        to graphics coordinates.
 */
class Transform2D
{
  public:
  
    Transform2D( PatchData* pd,
                 Projection& proj,
                 unsigned width, 
                 unsigned height,
                 bool flip_about_horizontal );

    Transform2D( const Vector3D* verts,
                 size_t num_vert,
                 Projection& projection,
                 unsigned width, 
                 unsigned height );
    
    void transform( const Vector3D& coords,
                    int& horizontal,
                    int& vertical ) const;
    
    int max_horizontal() const { return horizMax; }
    int max_vertical() const   { return vertMax; } 
    
  private:
  
    Projection& myProj;
    float myScale;
    int horizOffset, vertOffset;
    int horizMax, vertMax;
    int vertSign;
};
                             


/* Write VTK file
 *
 * Copied from src/Mesh/MeshSet.cpp and adapted for removal of 
 * MeshSet class by J.Kraftcheck on 2005-7-28.
 *
 * This code is provided mainly for debugging.  A more efficient
 * and complete writer implementation is provided in the MeshImpl
 * class for saving meshes that were read from a file initially.
 */
void write_vtk( Mesh* mesh, const char* out_filename, MsqError &err )
{
  if (MeshImpl* msq_mesh = dynamic_cast<MeshImpl*>(mesh)) {
    msq_mesh->write_vtk( out_filename, err );
    MSQ_CHKERR(err);
    return;
  }

    // loads a global patch
  PatchData pd;
  pd.set_mesh( mesh );
  pd.fill_global_patch( err ); MSQ_ERRRTN(err);
  
    // write mesh
  write_vtk( pd, out_filename, err ); MSQ_CHKERR(err);
}


void write_vtk( PatchData& pd, const char* out_filename, MsqError &err,
                const Vector3D* OF_gradient)
{
    // Open the file
  std::ofstream file(out_filename);
  if (!file)
  {
    MSQ_SETERR(err)(MsqError::FILE_ACCESS);
    return;
  }
    
    // Write a header
  file << "# vtk DataFile Version 2.0\n";
  file << "Mesquite Mesh " << out_filename << " .\n";
  file << "ASCII\n";
  file << "DATASET UNSTRUCTURED_GRID\n";
  
    // Write vertex coordinates
  file << "POINTS " << pd.num_nodes() << " float\n";
  size_t i;
  for (i = 0; i < pd.num_nodes(); i++)
  {
    file << pd.vertex_by_index(i)[0] << ' '
         << pd.vertex_by_index(i)[1] << ' '
         << pd.vertex_by_index(i)[2] << '\n';
  }
  
    // Write out the connectivity table
  size_t connectivity_size = 0;
  for (i = 0; i < pd.num_elements(); ++i)
    connectivity_size += pd.element_by_index(i).node_count()+1;
    
  file << "CELLS " << pd.num_elements() << ' ' << connectivity_size << '\n';
  for (i = 0; i < pd.num_elements(); i++)
  {
    std::vector<size_t> vtx_indices;
    pd.element_by_index(i).get_node_indices(vtx_indices);
    
      // Convert native to VTK node order, if not the same
    const VtkTypeInfo* info = VtkTypeInfo::find_type( pd.element_by_index(i).get_element_type(),
                                                      vtx_indices.size(),
                                                      err ); MSQ_ERRRTN(err);
    info->mesquiteToVtkOrder( vtx_indices );
     
    file << vtx_indices.size();
    for (std::size_t j = 0; j < vtx_indices.size(); ++j)
    {
      file << ' ' << vtx_indices[j];
    }
    file << '\n';
  }
  
    // Write out the element types
  file << "CELL_TYPES " << pd.num_elements() << '\n';
  for (i = 0; i < pd.num_elements(); i++)
  {
    const VtkTypeInfo* info = VtkTypeInfo::find_type( 
                               pd.element_by_index(i).get_element_type(),
                               pd.element_by_index(i).node_count(),
                               err ); MSQ_ERRRTN(err);
    file << info->vtkType << '\n';
  }
  
    // Write out which points are fixed.
  file << "POINT_DATA " << pd.num_nodes()
       << "\nSCALARS fixed int\nLOOKUP_TABLE default\n";
  for (i = 0; i < pd.num_nodes(); ++i)
  {
    if (pd.vertex_by_index(i).get_flags() & MsqVertex::MSQ_CULLED)
      file << "1\n";
    else
      file << "0\n";
  }
  file << "SCALARS culled short\nLOOKUP_TABLE default\n";
  for (i = 0; i < pd.num_nodes(); ++i)
  {
    if (pd.vertex_by_index(i).is_free_vertex())
      file << "0\n";
    else
      file << "1\n";
  }
  
  
  if (OF_gradient) {
    file << "VECTORS gradient double\n";
    for (i = 0; i < pd.num_free_vertices(); ++i) 
      file << OF_gradient[i].x() << " " << OF_gradient[i].y() << " " << OF_gradient[i].z() << "\n";
    for (i = pd.num_free_vertices(); i < pd.num_nodes(); ++i)
      file << "0.0 0.0 0.0\n";
  }
  
  
  
    // Close the file
  file.close();
}


void write_gnuplot( Mesh* mesh, const char* out_filebase, MsqError &err)
{
    // loads a global patch
  PatchData pd;
  pd.set_mesh( mesh );
  pd.fill_global_patch( err ); MSQ_ERRRTN(err);
  write_gnuplot( pd, out_filebase, err );
}

void write_gnuplot( Mesh* mesh, std::vector<Mesh::ElementHandle>& elems,
                    const char* out_filebase, MsqError &err)
{
    // loads a global patch
  PatchData pd;
  pd.set_mesh( mesh );
  std::vector<Mesh::VertexHandle> verts;
  pd.set_mesh_entities( elems, verts, err ); MSQ_ERRRTN(err);
  write_gnuplot( pd, out_filebase, err );
}


/*  Writes a gnuplot file directly from the MeshSet.
 *  This means that any mesh imported successfully into Mesquite
 *  can be outputed in gnuplot format.
 *
 *  Within gnuplot, use \b plot 'file1.gpt' w l, 'file2.gpt' w l  
 *   
 *  This is not geared for performance, since it has to load a global Patch from
 *  the mesh to write a mesh file. 
 *
 * Copied from src/Mesh/MeshSet.cpp and adapted for removal of 
 * MeshSet class by J.Kraftcheck on 2005-7-28.
 * 
 * Re-written to use EdgeIterator by J.Kraftcheck on 2009-6-11
*/
void write_gnuplot( PatchData& pd, const char* out_filebase, MsqError& err )
{
    // Open the file
  std::string out_filename = out_filebase;
  out_filename += ".gpt";
  std::ofstream file(out_filename.c_str());
  if (!file)
  {
    MSQ_SETERR(err)(MsqError::FILE_ACCESS);
    return;
  }

  EdgeIterator edges( &pd, err ); MSQ_ERRRTN(err);
   
    // Write a header
  file << "\n";
  
  while (!edges.is_at_end())
  {
    const Vector3D& s = edges.start();
    const Vector3D& e = edges.end();
    const Vector3D* m = edges.mid();
    
    file << s[0] << ' ' << s[1] << ' ' << s[2] << std::endl;
    if (m) 
      file << (*m)[0] << ' ' << (*m)[1] << ' ' << (*m)[2] << std::endl;
    file << e[0] << ' ' << e[1] << ' ' << e[2] << std::endl;
    file << std::endl << std::endl;
    
    edges.step(err); MSQ_ERRRTN(err);
  }
  
    // Close the file
  file.close();
}

/** Helper function for write_gnuplot_animator and write_gnuplot_overlay
 *
 * Read a set of input files to determine the bounding box
 * of the combined data.
 */
static void find_gnuplot_agregate_range( int count,
                                         const char* basename,
                                         Vector3D& min,
                                         Vector3D& max,
                                         MsqError& err )
{
    // read all input files to determine extents
  min = Vector3D( HUGE_VAL, HUGE_VAL, HUGE_VAL);
  max = Vector3D(-HUGE_VAL,-HUGE_VAL,-HUGE_VAL);
  for (int i = 0; i <= count; ++i) {
    stringstream s;
    s << basename << '.' << i << ".gpt";
    ifstream infile( s.str().c_str() );
    if (!infile) {
      MSQ_SETERR(err)(s.str(), MsqError::FILE_ACCESS);
      return;
    }
    double c[3];
    while (infile >> c[0] >> c[1] >> c[2]) {
      for (int j = 0; j < 3; ++j) {
        if (c[j] < min[j])
          min[j] = c[j];
        if (c[j] > max[j])
          max[j] = c[j];
      }
    }
  }
}
  
/** Write a GNU Plot script to produce an animation from a 
 *  sequence of data files 
 */
void write_gnuplot_animator( int count, const char* basename, MsqError& err )
{
  if (count <= 0)
    return;

  const int DELAY = 10;
  const int WIDTH = 640;
  const int HEIGHT = 480;
  
    // read all input files to determine extents
  Vector3D min, max;
  find_gnuplot_agregate_range( count, basename, min, max, err ); MSQ_ERRRTN(err);

    // chose coordinate plane to plot in
  Vector3D range = max - min;
  int haxis = 0, vaxis = 1;
  if (range[0] < range[1] && range[1] < range[2]) {
    haxis = 1;
    vaxis = 2;
  }
  else if (range[1] < range[2]) {
    vaxis = 2;
  }
  
    // open output file
  string base(basename);
  ofstream file( (string(basename) + ".gnuplot").c_str() );
  if (!file)
  {
    MSQ_SETERR(err)(MsqError::FILE_ACCESS);
    return;
  }
  
    // write header
  file << "#!gnuplot" << endl;
  file << "#" << endl;
  file << "# Mesquite Animation of " << basename << ".0 to " << basename << '.' << count << endl;
  file << "#" << endl;
  
    // write plot settings
  file << "set term gif animate transparent opt delay " << DELAY << " size " << WIDTH << "," << HEIGHT << endl;
  file << "set xrange [" << min[haxis] - 0.05 * range[haxis] << ":" << max[haxis] + 0.05 * range[haxis] << "]" << endl;
  file << "set yrange [" << min[vaxis] - 0.05 * range[vaxis] << ":" << max[vaxis] + 0.05 * range[vaxis] << "]" << endl;
  file << "set title '" << basename << "'" << endl;
  file << "set output '" << basename << ".gif'" << endl;
  
    // plot data
  for (int i = 0; i <= count; ++i) 
    file << "plot '" << basename << '.' << i << ".gpt'" 
         << " using " << haxis+1 << ":" << vaxis+1 << " w l" << endl;
}



static unsigned red( int i, int c )
{
  if (i * 4 <= c)
    return 255;
  else if (i * 4 >= 3 * c)
    return 0;
  else
    return 384 - i * 511 / c;
}

static unsigned green( int i, int c )
{
  if (i * 4 < c)
    return i * 510 / c;
  else if (i * 4 > 3 * c)
    return 1023 - i * 1023 / c;
  else
    return 255;
}

static unsigned blue( int i, int c )
{
  if (i * 4 <= c)
    return 0;
  else if (i * 4 >= 3 * c)
    return 255;
  else
    return i * 511 / c - 127;
}

/** Write a GNU Plot script to produce a single plot from a 
 *  sequence of data files 
 */
void write_gnuplot_overlay( int count, const char* basename, MsqError& err )
{
  if (count <= 0)
    return;

  const int WIDTH = 640;
  const int HEIGHT = 480;
  
    // read all input files to determine extents
  Vector3D min, max;
  find_gnuplot_agregate_range( count, basename, min, max, err ); MSQ_ERRRTN(err);

    // chose coordinate plane to plot in
  Vector3D range = max - min;
  int haxis = 0, vaxis = 1;
  if (range[0] < range[1] && range[1] < range[2]) {
    haxis = 1;
    vaxis = 2;
  }
  else if (range[1] < range[2]) {
    vaxis = 2;
  }
  
    // open output file
  string base(basename);
  FILE* file = fopen( (string(basename) + ".gnuplot").c_str(), "w" );
  if (!file)
  {
    MSQ_SETERR(err)(MsqError::FILE_ACCESS);
    return;
  }
  
    // write header
  fprintf( file, "#!gnuplot\n" );
  fprintf( file, "#\n" );
  fprintf( file, "# Mesquite Overlay of %s.0 to %s.%d\n", basename, basename, count );
  fprintf( file, "#\n" );
  
    // write plot settings
  fprintf( file, "set term gif size %d,%d\n", WIDTH, HEIGHT );
  fprintf( file, "set xrange [%f:%f]\n", min[haxis] - 0.05 * range[haxis], max[haxis] + 0.05 * range[haxis] );
  fprintf( file, "set yrange [%f:%f]\n", min[vaxis] - 0.05 * range[vaxis], max[vaxis] + 0.05 * range[vaxis] );
  fprintf( file, "set title '%s'\n", basename );
  fprintf( file, "set output '%s.gif'\n", basename );
  
    // plot data
  fprintf( file, "plot '%s.0.gpt' using %d:%d w l lc rgb '#%02x%02x%02x' title 't0'",
           basename, haxis+1, vaxis+1, red(0,count), green(0,count), blue(0,count) );
  for (int i = 1; i <= count; ++i) {
    fprintf( file, ", \\\n     '%s.%d.gpt' using %d:%d w l lc rgb '#%02x%02x%02x' title 't%d'",
             basename, i, haxis+1, vaxis+1, red(i,count), green(i,count), blue(i,count), i);
  }
  fprintf( file, "\n" );
  
  fclose(file);
}
  

void write_stl( Mesh* mesh, const char* filename, MsqError& err )
{
    // Open the file
  ofstream file(filename);
  if (!file)
  {
    MSQ_SETERR(err)(MsqError::FILE_ACCESS);
    return;
  }
  
    // Write header
  char header[70];
  sprintf( header, "Mesquite%d", rand() );
  file << "solid " << header << endl;
  
  MsqVertex coords[3];
  std::vector<Mesh::VertexHandle> verts(3);
  std::vector<size_t> offsets(2);
  
    // Iterate over all elements
  size_t count = 0;
  std::vector<Mesh::ElementHandle> elems;
  std::vector<Mesh::ElementHandle>::iterator iter;
  mesh->get_all_elements( elems, err ); MSQ_ERRRTN(err);
  for (iter = elems.begin(); iter != elems.end(); ++iter)
  {
      // Skip non-triangles
    Mesh::ElementHandle elem = *iter;
    EntityTopology type;
    mesh->elements_get_topologies( &elem, &type, 1, err ); MSQ_ERRRTN(err);
    if (type != TRIANGLE) 
      continue;
    ++count;
    
      // Get vertex coordinates
    mesh->elements_get_attached_vertices( &elem, 1, verts, offsets, err ); MSQ_ERRRTN(err);
    mesh->vertices_get_coordinates( arrptr(verts), coords, 3, err ); MSQ_ERRRTN(err);
    
      // Get triagnle normal
    Vector3D n = (coords[0] - coords[1]) * (coords[0] - coords[2]);
    n.normalize();
    
      // Write triangle
    file << "facet normal " << n.x() << " " << n.y() << " " << n.z() << endl;
    file << "outer loop" << endl;
    for (unsigned i = 0; i < 3; ++i)
      file << "vertex " << coords[i].x() << " " 
                        << coords[i].y() << " " 
                        << coords[i].z() << endl;
    file << "endloop" << endl;
    file << "endfacet" << endl;
  }
  
  file << "endsolid " << header << endl;
  
  file.close();
  if (count == 0)
  {
    std::remove(filename);
    MSQ_SETERR(err)("Mesh contains no triangles", MsqError::INVALID_STATE);
  }
}
  
  

Projection::Projection( PlanarDomain* domain )
{
  Vector3D normal;
  domain->vertex_normal_at( 0, normal );
  init( normal );
}

Projection::Projection( const Vector3D& n )
{ init( n ); }

Projection::Projection( const Vector3D& n, const Vector3D& up )
{ init( n, up ); }

Projection::Projection( Axis h, Axis v )
{
  Vector3D horiz(0,0,0), vert(0,0,0);
  horiz[h] = 1.0;
  vert[v] = 1.0;
  init( horiz * vert, vert );
}

void Projection::init( const Vector3D& n )
{
    // Choose an "up" direction

  Axis max = X;
  for (Axis i = Y; i <= Z; i = (Axis)(i+1))
    if (fabs(n[i]) > fabs(n[max]))
      max = i;
  
  Axis up;
  if (max == Y)
    up = Z;
  else 
    up = Y;
  
    // Initialize rotation matrix
  
  Vector3D up_vect(0,0,0);
  up_vect[up] = 1.0;
  init( n, up_vect );
}  

void Projection::init(  const Vector3D& n1, const Vector3D& up1 )
{
  MsqError err;
  const Vector3D n = n1/n1.length();
  const Vector3D u = up1/up1.length();
  
  // Rotate for projection
  const Vector3D z( 0., 0., 1. );
  Vector3D r = n * z;
  double angle = r.interior_angle( n, z, err );
  Matrix3D rot1 = rotation( r, angle );
   
  // In-plane rotation for up vector
  Vector3D pu = u - n * (n % u);
  Vector3D y( 0., 1., 0. );
  angle = z.interior_angle( pu, y, err );
  Matrix3D rot2 = rotation( z, angle );
  
  this->myTransform = rot1 * rot2;
}

Matrix3D Projection::rotation( const Vector3D& axis, double angle )
{
  const double c = cos( angle );
  const double s = sin( angle );
  const double x = axis.x();
  const double y = axis.y();
  const double z = axis.z(); 
  
  const double xform[9] = 
    {    c + x*x*(1-c),  -z*s + x*y*(1-c),  y*s + x*z*(1-c),
       z*s + x*y*(1-c),     c + y*y*(1-c), -x*s + y*z*(1-c),
      -y*s + x*z*(1-c),   x*s + y*z*(1-c),    c + z*z*(1-c) };
  return Matrix3D( xform );
}

void Projection::project( const Vector3D& p, float& h, float& v )
{
  Vector3D pr = myTransform * p;
  h = (float)pr.x();
  v = (float)pr.y();
}

Transform2D::Transform2D( PatchData* pd,
                          Projection& projection,
                          unsigned width, 
                          unsigned height,
                          bool flip )
  : myProj(projection),
    vertSign(flip ? -1 : 1)
{
    // Get the bounding box of the projected points
  float w_max, w_min, h_max, h_min;
  w_max = h_max = -std::numeric_limits<float>::max();
  w_min = h_min =  std::numeric_limits<float>::max();
  MsqError err;
  const MsqVertex* verts = pd->get_vertex_array( err );
  const size_t num_vert = pd->num_nodes();
  for (unsigned i = 0; i < num_vert; ++i)
  {
    float w, h;
    myProj.project( verts[i], w, h );
    if (w > w_max) w_max = w;
    if (w < w_min) w_min = w;
    if (h > h_max) h_max = h;
    if (h < h_min) h_min = h;
  }
  
    // Determine the scale factor
  const float w_scale = (float)width  / (w_max - w_min);
  const float h_scale = (float)height / (h_max - h_min);
  myScale = w_scale > h_scale ? h_scale : w_scale;
  
    // Determine offset
  horizOffset = -(int)(myScale * w_min);
  vertOffset  = -(int)(myScale * (flip ? -h_max : h_min));
  
    // Determine bounding box
  horizMax = (int)(                 w_max  * myScale) + horizOffset;
  vertMax  = (int)((flip ? -h_min : h_max) * myScale) +  vertOffset; 
    
}
 
Transform2D::Transform2D( const Vector3D* verts,
                          size_t num_vert,
                          Projection& projection,
                          unsigned width, 
                          unsigned height )
  : myProj(projection),
    vertSign(1)
{
    // Get the bounding box of the projected points
  float w_max, w_min, h_max, h_min;
  w_max = h_max = -std::numeric_limits<float>::max();
  w_min = h_min =  std::numeric_limits<float>::max();
  for (unsigned i = 0; i < num_vert; ++i)
  {
    float w, h;
    myProj.project( verts[i], w, h );
    if (w > w_max) w_max = w;
    if (w < w_min) w_min = w;
    if (h > h_max) h_max = h;
    if (h < h_min) h_min = h;
  }
  
    // Determine the scale factor
  const float w_scale = (float)width  / (w_max - w_min);
  const float h_scale = (float)height / (h_max - h_min);
  myScale = w_scale > h_scale ? h_scale : w_scale;
  
    // Determine offset
  horizOffset = -(int)(myScale * w_min);
  vertOffset  = -(int)(myScale * h_min);
  
    // Determine bounding box
  horizMax = (int)(w_max * myScale) + horizOffset;
  vertMax  = (int)(h_max * myScale) +  vertOffset; 
    
}
   
void Transform2D::transform( const Vector3D& coords,
                             int& horizontal,
                             int& vertical ) const
{
    float horiz, vert;
    myProj.project( coords, horiz, vert );
    horizontal =            (int)(myScale * horiz) + horizOffset;
    vertical   = vertSign * (int)(myScale *  vert) +  vertOffset;
}

/** Write quadratic edge shape in PostScript format.
 *
 * Given the three points composing a quadratic mesh edge,
 * write the cubic Bezier curve of the same shape in 
 * PostScript format.  The formulas for P1 and P2 
 * at the start of this function will result in the cubic
 * terms of the Bezier curve dropping out, leaving the
 * quadratic curve matching the edge shape function as 
 * described in Section 3.6 of Hughes.  (If you're attempting
 * to verify this, don't forget to adjust for the different
 * parameter ranges: \f$ \xi = 2 t - 1 \f$).
 */
static void write_eps_quadratic_edge( ostream &s,
                                      Transform2D& xform,
                                      Vector3D start,
                                      Vector3D mid, 
                                      Vector3D end )
{
  Vector3D P1 = 1./3 * (4 * mid - end);
  Vector3D P2 = 1./3 * (4 * mid - start);

  int x, y;
  xform.transform( start, x, y );
  s << x << ' ' << y << " moveto" << endl;
  xform.transform( P1, x, y );
  s << x << ' ' << y << ' ';
  xform.transform( P2, x, y );
  s << x << ' ' << y << ' ';
  xform.transform( end, x, y );
  s << x << ' ' << y << " curveto" << endl;
}

void write_eps( Mesh* mesh, 
                const char* filename, 
                Projection proj, 
                MsqError& err,
                int width, int height )
{
    // Get a global patch
  PatchData pd;
  pd.set_mesh( mesh );
  pd.fill_global_patch( err ); MSQ_ERRRTN(err);
  
  Transform2D transf( &pd, proj, width, height, false );
    
    // Open the file
  ofstream s(filename);
  if (!s)
  {
    MSQ_SETERR(err)(MsqError::FILE_ACCESS);
    return;
  }

    // Write header
  s << "%!PS-Adobe-2.0 EPSF-2.0"                      << endl;
  s << "%%Creator: Mesquite"                          << endl;
  s << "%%Title: Mesquite "                           << endl;
  s << "%%DocumentData: Clean7Bit"                    << endl;
  s << "%%Origin: 0 0"                                << endl;
  s << "%%BoundingBox: 0 0 " 
    << transf.max_horizontal() <<  ' ' 
    << transf.max_vertical()                          << endl;
  s << "%%Pages: 1"                                   << endl;
  
  s << "%%BeginProlog"                                << endl;
  s << "save"                                         << endl;
  s << "countdictstack"                               << endl;
  s << "mark"                                         << endl;
  s << "newpath"                                      << endl;
  s << "/showpage {} def"                             << endl;
  s << "/setpagedevice {pop} def"                     << endl;
  s << "%%EndProlog"                                  << endl;
  
  s << "%%Page: 1 1"                                  << endl;
  s << "1 setlinewidth"                               << endl;
  s << "0.0 setgray"                                  << endl;
  
    // Write mesh edges
  EdgeIterator iter( &pd, err );  MSQ_ERRRTN(err);
  while( !iter.is_at_end() )
  {
    int s_w, s_h, e_w, e_h;
    transf.transform( iter.start(), s_w, s_h );
    transf.transform( iter.end  (), e_w, e_h );
    
    s << "newpath"                                    << endl;
    s << s_w << ' ' << s_h << " moveto"               << endl;
    
    if (!iter.mid()) {
      s << e_w << ' ' << e_h << " lineto"               << endl;
    }
    else {
      write_eps_quadratic_edge( s, transf, iter.start(), *iter.mid(), iter.end() );
        // draw rings at mid-edge node location
      //transf.transform( *(iter.mid()), w1, h1 );
      //s << w1+2 << ' ' << h1 <<  " moveto"            << endl;
      //s << w1 << ' ' << h1 <<  " 2 0 360 arc"         << endl;
    }
    s << "stroke"                                       << endl;
    
    iter.step(err); MSQ_ERRRTN(err);
  }
  
    // Write footer
  s << "%%Trailer"                                    << endl;
  s << "cleartomark"                                  << endl;
  s << "countdictstack"                               << endl;
  s << "exch sub { end } repeat"                      << endl;
  s << "restore"                                      << endl;  
  s << "%%EOF"                                        << endl;
}

/** Quadratic triangle shape function for use in write_eps_triangle */
static double tN0( double r, double s ) { double t = 1 - r - s; return t*(2*t - 1); }
static double tN1( double r, double   ) { return r*(2*r - 1); }
static double tN2( double  , double s ) { return s*(2*s - 1); }
static double tN3( double r, double s ) { double t = 1 - r - s; return 4*r*t; }
static double tN4( double r, double s ) { return 4*r*s; }
static double tN5( double r, double s ) { double t = 1 - r - s; return 4*s*t; }
/** Quadratic triangle shape function for use in write_eps_triangle */
static Vector3D quad_tri_pt( double r, double s, const Vector3D* coords )
{
  Vector3D result = tN0(r,s) * coords[0];
  result += tN1(r,s) * coords[1];
  result += tN2(r,s) * coords[2];
  result += tN3(r,s) * coords[3];
  result += tN4(r,s) * coords[4];
  result += tN5(r,s) * coords[5];
  return result;
}

void write_eps_triangle( Mesh* mesh, 
                         Mesh::ElementHandle elem,
                         const char* filename, 
                         bool draw_iso_lines, 
                         bool draw_nodes,
                         MsqError& err,
                         int width, int height )
{
    // Get triangle vertices
  MsqVertex coords[6];
  EntityTopology type;
  mesh->elements_get_topologies( &elem, &type, 1, err ); MSQ_ERRRTN(err);
  if (type != TRIANGLE) {
    MSQ_SETERR(err)("Invalid element type", MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  std::vector<Mesh::VertexHandle> verts;
  std::vector<size_t> junk;
  mesh->elements_get_attached_vertices( &elem, 1, verts, junk, err );  MSQ_ERRRTN(err);
  if (verts.size() != 3 && verts.size() != 6) {
    MSQ_SETERR(err)("Invalid element type", MsqError::UNSUPPORTED_ELEMENT );
    return;
  }
  mesh->vertices_get_coordinates( arrptr(verts), coords, verts.size(), err ); MSQ_ERRRTN(err);
  
  Vector3D coords2[6];
  std::copy( coords, coords+verts.size(), coords2 );
  
  std::vector<bool> fixed(verts.size(), false);
  if (draw_nodes) {
    mesh->vertices_get_fixed_flag( arrptr(verts), fixed, verts.size(), err ); MSQ_ERRRTN(err);
  }
  write_eps_triangle( coords2, verts.size(), filename, draw_iso_lines, draw_nodes, err, fixed, width, height );
}

void write_eps_triangle( const Vector3D* coords,
                         size_t num_vtx,
                         const char* filename, 
                         bool draw_iso_lines, 
                         bool draw_nodes,
                         MsqError& err,
                         const std::vector<bool>& fixed,
                         int width, int height )
{
  const int PT_RAD = 3; // radius of circles for drawing nodes, in points
  const int PAD = PT_RAD + 2; // margin in points
  const double EDGE_GRAY  = 0.0; // color for triangle edges, 0.0 => black
  const double ISO_GRAY   = 0.7; // color for parameter iso-lines
  const double NODE_GRAY  = 0.0; // color for node circle 
  const double FIXED_GRAY = 1.0; // color to fill fixed nodes with, 1.0 => white
  const double FREE_GRAY  = 0.0; // color to fill free nodes with

  Projection proj( X, Y );
  Transform2D transf( coords, num_vtx, proj, width, height );
    
    // Open the file
  ofstream str(filename);
  if (!str)
  {
    MSQ_SETERR(err)(MsqError::FILE_ACCESS);
    return;
  }

    // Write header
  str << "%!PS-Adobe-2.0 EPSF-2.0"                      << endl;
  str << "%%Creator: Mesquite"                          << endl;
  str << "%%Title: Mesquite "                           << endl;
  str << "%%DocumentData: Clean7Bit"                    << endl;
  str << "%%Origin: 0 0"                                << endl;
  str << "%%BoundingBox: " << -PAD << ' ' << -PAD << ' ' 
      << transf.max_horizontal() + PAD <<  ' ' 
      << transf.max_vertical()   + PAD                  << endl;
  str << "%%Pages: 1"                                   << endl;
  
  str << "%%BeginProlog"                                << endl;
  str << "save"                                         << endl;
  str << "countdictstack"                               << endl;
  str << "mark"                                         << endl;
  str << "newpath"                                      << endl;
  str << "/showpage {} def"                             << endl;
  str << "/setpagedevice {pop} def"                     << endl;
  str << "%%EndProlog"                                  << endl;
  
  str << "%%Page: 1 1"                                  << endl;
  str << "1 setlinewidth"                               << endl;
  str << EDGE_GRAY << " setgray"                        << endl;

  const double h = 0.5, t = 1.0/3.0, w = 2.0/3.0, s = 1./6, f = 5./6;
  const int NUM_ISO = 15;
  const double iso_params[NUM_ISO][2][2] = 
    { { { h, 0 }, { h, h } },  // r = 1/2
      { { t, 0 }, { t, w } },  // r = 1/3
      { { w, 0 }, { w, t } },  // r = 2/3
      { { s, 0 }, { s, f } },  // r = 1/6
      { { f, 0 }, { f, s } },  // r = 5/6
      { { 0, h }, { h, h } },  // s = 1/2
      { { 0, t }, { w, t } },  // s = 1/3
      { { 0, w }, { t, w } },  // s = 2/3
      { { 0, s }, { f, s } },  // s = 1/6
      { { 0, f }, { s, f } },  // s = 5/6
      { { 0, h }, { h ,0 } },  // t = 1 - r - s = 1/2
      { { 0, w }, { w, 0 } },  // t = 1 - r - s = 1/3
      { { 0, t }, { t, 0 } },  // t = 1 - r - s = 2/3
      { { 0, f }, { f, 0 } },  // t = 1 - r - s = 1/6
      { { 0, s }, { s, 0 } }   // t = 1 - r - s = 5/6
     };

  if (num_vtx == 3) {
    int x[3], y[3];
    for (size_t i = 0; i < 3; ++i) 
      transf.transform( coords[i], x[i], y[i] );
    
    str << "newpath"                                      << endl;
    str << x[0] << ' ' << y[0] << " moveto"               << endl;
    str << x[1] << ' ' << y[1] << " lineto"               << endl;
    str << x[2] << ' ' << y[2] << " lineto"               << endl;
    str << x[0] << ' ' << y[0] << " lineto"               << endl;
    str << "stroke"                                       << endl;
    
    if (draw_iso_lines) {
      str << ISO_GRAY << " setgray"                       << endl;
      str << "newpath"                                    << endl;
      for (int i = 0; i < NUM_ISO; ++i) {
        double R[2] = { iso_params[i][0][0], iso_params[i][1][0] };
        double S[2] = { iso_params[i][0][1], iso_params[i][1][1] };
        double T[2] = { 1 - R[0] - S[0], 1 - R[1] - S[1] };
        Vector3D p[2] = { T[0] * coords[0] + R[0] * coords[1] + S[0] * coords[2],
                          T[1] * coords[0] + R[1] * coords[1] + S[1] * coords[2] };
        transf.transform( p[0], x[0], y[0] );
        transf.transform( p[1], x[1], y[1] );
        str << x[0] << ' ' << y[0] << " moveto"           << endl;
        str << x[1] << ' ' << y[1] << " lineto"           << endl;
        
      }
      str << "    stroke"                                 << endl;
    }
  }
  else if (num_vtx == 6) {
    str << "newpath"                                      << endl;
    write_eps_quadratic_edge( str, transf, coords[0], coords[3], coords[1] );
    write_eps_quadratic_edge( str, transf, coords[1], coords[4], coords[2] );
    write_eps_quadratic_edge( str, transf, coords[2], coords[5], coords[0] );
    str << "stroke"                                       << endl;
    
    if (draw_iso_lines) {
      str << ISO_GRAY << " setgray"                       << endl;
      str << "newpath"                                    << endl;
      for (int i = 0; i < NUM_ISO; ++i) {
        double R[3] = { iso_params[i][0][0], 0, iso_params[i][1][0] };
        double S[3] = { iso_params[i][0][1], 0, iso_params[i][1][1] };
        R[1] = 0.5*(R[0]+R[2]);
        S[1] = 0.5*(S[0]+S[2]);
        Vector3D p[3] = { quad_tri_pt( R[0], S[0], coords ),
                          quad_tri_pt( R[1], S[1], coords ),
                          quad_tri_pt( R[2], S[2], coords ) };
        write_eps_quadratic_edge( str, transf, p[0], p[1], p[2] );
      }
      str << "    stroke"                                 << endl;
    }
  }
  
  if (draw_nodes) {
    for (size_t i = 0; i < num_vtx; ++i) {
      int w, h;
        // fill interior with either white or black depending
        // on whether or not the vertex is fixed.
      if (fixed[i]) 
        str << FIXED_GRAY << " setgray"                     << endl;
      else
        str << FREE_GRAY << " setgray"                      << endl;
      transf.transform( coords[i], w, h );
      str << w+PT_RAD << ' ' << h << " moveto"              << endl;
      str << w << ' ' << h << ' ' << PT_RAD << " 0 360 arc" << endl;
      str << "closepath fill"                               << endl;
      str << NODE_GRAY << " setgray"                        << endl;
      str << "newpath"                                      << endl;
      str << w << ' ' << h << ' ' << PT_RAD << " 0 360 arc" << endl;
      str << "stroke"                                       << endl;
    }
  }
  
    // Write footer
  str << "%%Trailer"                                    << endl;
  str << "cleartomark"                                  << endl;
  str << "countdictstack"                               << endl;
  str << "exch sub { end } repeat"                      << endl;
  str << "restore"                                      << endl;  
  str << "%%EOF"                                        << endl;
}

void write_svg( Mesh* mesh, 
                const char* filename, 
                Projection proj, 
                MsqError& err )
{
    // Get a global patch
  PatchData pd;
  pd.set_mesh( mesh );
  pd.fill_global_patch( err ); MSQ_ERRRTN(err);
  
  Transform2D transf( &pd, proj, 400, 400, true );
    
    // Open the file
  ofstream file(filename);
  if (!file)
  {
    MSQ_SETERR(err)(MsqError::FILE_ACCESS);
    return;
  }

    // Write header
  file << "<?xml version=\"1.0\" standalone=\"no\"?>"                << endl;
  file << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" " 
       << "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">"    << endl;
  file <<                                                               endl;
  file << "<svg width=\"" << transf.max_horizontal() 
       << "\" height=\"" << transf.max_vertical()
       << "\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">" << endl;
  
    // Write mesh edges
  EdgeIterator iter( &pd, err );  MSQ_ERRRTN(err);
  while( !iter.is_at_end() )
  {
    int s_w, s_h, e_w, e_h;
    transf.transform( iter.start(), s_w, s_h );
    transf.transform( iter.end  (), e_w, e_h );
    
    file << "<line "
         << "x1=\"" << s_w << "\" "
         << "y1=\"" << s_h << "\" "
         << "x2=\"" << e_w << "\" "
         << "y2=\"" << e_h << "\" "
         << " style=\"stroke:rgb(99,99,99);stroke-width:2\""
         << "/>" << endl;
    
    iter.step( err ); MSQ_ERRRTN(err);
  }
  
    // Write footer
  file << "</svg>" << endl;
}

} // namespace MeshWriter

} // namespace Mesquite

#endif
