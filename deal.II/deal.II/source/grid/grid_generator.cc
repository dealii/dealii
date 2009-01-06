//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/quadrature_lib.h>
#include <base/thread_management.h>
#include <lac/vector.h>
#include <lac/vector_memory.h>
#include <lac/filtered_matrix.h>
#include <lac/precondition.h>
#include <lac/solver_cg.h>
#include <lac/sparse_matrix.h>
#include <grid/grid_generator.h>
#include <grid/grid_reordering.h>
#include <grid/grid_tools.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_constraints.h>
#include <fe/mapping_q1.h>
#include <fe/fe_q.h>
#include <numerics/matrices.h>

#include <iostream>
#include <cmath>
#include <limits>

DEAL_II_NAMESPACE_OPEN



namespace
{
#if deal_II_dimension == 3

				   // Corner points of the cube [-1,1]^3
  const Point<3> hexahedron[8] =
  {
	Point<3>(-1,-1,-1),
	Point<3>(+1,-1,-1),
	Point<3>(-1,+1,-1),
	Point<3>(+1,+1,-1),
	Point<3>(-1,-1,+1),
	Point<3>(+1,-1,+1),
	Point<3>(-1,+1,+1),
	Point<3>(+1,+1,+1)
  };

				   // Octahedron inscribed in the cube
				   // [-1,1]^3
  const Point<3> octahedron[6] =
  {
	Point<3>(-1, 0, 0),
	Point<3>( 1, 0, 0),
	Point<3>( 0,-1, 0),
	Point<3>( 0, 1, 0),
	Point<3>( 0, 0,-1),
	Point<3>( 0, 0, 1)
  };
  
#endif
}


template <int dim, int spacedim>
void
GridGenerator::hyper_rectangle (Triangulation<dim,spacedim> &tria,
				const Point<spacedim>   &p_1,
				const Point<spacedim>   &p_2,
				const bool          colorize)
{
				   // First, normalize input such that
				   // p1 is lower in all coordinate directions.
  Point<spacedim> p1(p_1);
  Point<spacedim> p2(p_2);
  
  for (unsigned int i=0;i<spacedim;++i)
    if (p1(i) > p2(i))
      std::swap (p1(i), p2(i));
  
  std::vector<Point<spacedim> > vertices (GeometryInfo<dim>::vertices_per_cell);
  switch (dim)
    {
      case 1:
	    vertices[0] = p1;
	    vertices[1] = p2;
	    break;
      case 2:
	    vertices[0] = vertices[1] = p1;
	    vertices[2] = vertices[3] = p2;
	    
	    vertices[1](0) = p2(0);
	    vertices[2](0) = p1(0);
	    break;
      case 3:
	    vertices[0] = vertices[1] = vertices[2] = vertices[3] = p1;
	    vertices[4] = vertices[5] = vertices[6] = vertices[7] = p2;
	    
	    vertices[1](0) = p2(0);
	    vertices[2](1) = p2(1);
	    vertices[3](0) = p2(0);
	    vertices[3](1) = p2(1);
	    
	    vertices[4](0) = p1(0);
	    vertices[4](1) = p1(1);
	    vertices[5](1) = p1(1);
	    vertices[6](0) = p1(0);
	    
	    break;
      default:
	    Assert (false, ExcNotImplemented ());
    }

				   // Prepare cell data
  std::vector<CellData<dim> > cells (1);
  for (unsigned int i=0;i<GeometryInfo<dim>::vertices_per_cell;++i)
    cells[0].vertices[i] = i;
  cells[0].material_id = 0;

  tria.create_triangulation (vertices, cells, SubCellData());

				   // Assign boundary indicators
  if (colorize && (spacedim == dim))
    colorize_hyper_rectangle (tria);
}



#if deal_II_dimension == 1

// Implementation for 1D only
template <int dim, int spacedim>
void
GridGenerator::colorize_hyper_rectangle (Triangulation<dim,spacedim> &)
{
				   // Nothing to do in 1D
}


#else

// Implementation for dimensions except 1
template <int dim, int spacedim>
void
GridGenerator::colorize_hyper_rectangle (Triangulation<dim,spacedim> &tria)
{
				   // there is only one cell, so
				   // simple task
  const typename Triangulation<dim,spacedim>::cell_iterator cell = tria.begin();
  for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
    cell->face(f)->set_boundary_indicator (f);
}

#endif


template <int dim, int spacedim>
void GridGenerator::hyper_cube (Triangulation<dim,spacedim> &tria,
			        const double        left,
			        const double        right)
{
  Assert (left < right,
	  ExcMessage ("Invalid left and right bounds of hypercube"));
  
  Point<spacedim> p1;
  Point<spacedim> p2;

  p1(spacedim-1) = 0;
  p2(spacedim-1) = 0;

  for (unsigned int i=0;i<dim;++i)
    {
      p1(i) = left;
      p2(i) = right;
    }
  hyper_rectangle (tria, p1, p2);
}


#if deal_II_dimension == 3

void
GridGenerator::moebius (
  Triangulation<3>&  tria,
  const unsigned int   n_cells,
  const unsigned int   n_rotations,
  const double         R,
  const double         r)
{
  const unsigned int dim=3;
  Assert (n_cells>4, ExcMessage("More than 4 cells are needed to create a moebius grid."));
  Assert (r>0 && R>0, ExcMessage("Outer and inner radius must be positive."));
  Assert (R>r, ExcMessage("Outer radius must be greater than inner radius."));
    
  
  std::vector<Point<dim> > vertices (4*n_cells);
  double beta_step=n_rotations*numbers::PI/2.0/n_cells;
  double alpha_step=2.0*numbers::PI/n_cells;
  
  for (unsigned int i=0; i<n_cells; ++i)
    for (unsigned int j=0; j<4; ++j)
      {
	vertices[4*i+j][0]=R*std::cos(i*alpha_step)+r*std::cos(i*beta_step+j*numbers::PI/2.0)*std::cos(i*alpha_step);
	vertices[4*i+j][1]=R*std::sin(i*alpha_step)+r*std::cos(i*beta_step+j*numbers::PI/2.0)*std::sin(i*alpha_step);
	vertices[4*i+j][2]=r*std::sin(i*beta_step+j*numbers::PI/2.0);
      }

  unsigned int offset=0;
  
  std::vector<CellData<dim> > cells (n_cells);
  for (unsigned int i=0; i<n_cells; ++i)
    {
      for (unsigned int j=0; j<2; ++j)
	{
	  cells[i].vertices[0+4*j]=offset+0+4*j;
	  cells[i].vertices[1+4*j]=offset+3+4*j;
	  cells[i].vertices[2+4*j]=offset+2+4*j;
	  cells[i].vertices[3+4*j]=offset+1+4*j;
	}
      offset+=4;
      cells[i].material_id=0;
    }
  
				   // now correct the last four vertices
  cells[n_cells-1].vertices[4]=(0+n_rotations)%4;
  cells[n_cells-1].vertices[5]=(3+n_rotations)%4;
  cells[n_cells-1].vertices[6]=(2+n_rotations)%4;
  cells[n_cells-1].vertices[7]=(1+n_rotations)%4;

  GridReordering<dim>::invert_all_cells_of_negative_grid(vertices,cells);
  tria.create_triangulation_compatibility (vertices, cells, SubCellData());
}

#endif



// Implementation for 2D only
template<int dim>
void
GridGenerator::parallelogram (
  Triangulation<dim>&  tria,
  const Tensor<2,dim>& corners,
  const bool      colorize)
{
  Assert (dim==2, ExcNotImplemented());
  
  std::vector<Point<dim> > vertices (GeometryInfo<dim>::vertices_per_cell);

  vertices[1] = corners[0];
  vertices[2] = corners[1];
  vertices[3] = vertices[1];
  vertices[3] += vertices[2];
				   // Prepare cell data
  std::vector<CellData<dim> > cells (1);
  for (unsigned int i=0;i<GeometryInfo<dim>::vertices_per_cell;++i)
    cells[0].vertices[i] = i;
  cells[0].material_id = 0;

  tria.create_triangulation (vertices, cells, SubCellData());

				   // Assign boundary indicators
  if (colorize)
    colorize_hyper_rectangle (tria);
}


template <int dim>
void
GridGenerator::subdivided_hyper_cube (Triangulation<dim> &tria,
                                      const unsigned int  repetitions,
                                      const double        left,
                                      const double        right)
{
  Assert (repetitions >= 1, ExcInvalidRepetitions(repetitions));
  Assert (left < right,
	  ExcMessage ("Invalid left and right bounds of hypercube"));
  
                                   // first generate the necessary
                                   // points
  const double delta = (right-left)/repetitions;
  std::vector<Point<dim> > points;
  switch (dim)
    {
      case 1:
            for (unsigned int x=0; x<=repetitions; ++x)
              points.push_back (Point<dim> (left+x*delta));
            break;

      case 2:
            for (unsigned int y=0; y<=repetitions; ++y)
              for (unsigned int x=0; x<=repetitions; ++x)
                points.push_back (Point<dim> (left+x*delta,
                                              left+y*delta));
            break;

      case 3:
            for (unsigned int z=0; z<=repetitions; ++z)
              for (unsigned int y=0; y<=repetitions; ++y)
                for (unsigned int x=0; x<=repetitions; ++x)
                  points.push_back (Point<dim> (left+x*delta,
                                                left+y*delta,
                                                left+z*delta));
            break;

      default:
            Assert (false, ExcNotImplemented());
    }

                                   // next create the cells
				   // Prepare cell data
  std::vector<CellData<dim> > cells;
				   // Define these as abbreviations
				   // for the step sizes below. The
				   // number of points in a single
				   // direction is repetitions+1
  const unsigned int dy = repetitions+1;
  const unsigned int dz = dy*dy;
  switch (dim)
    {
      case 1:
	cells.resize (repetitions);
	for (unsigned int x=0; x<repetitions; ++x)
	  {
	    cells[x].vertices[0] = x;
	    cells[x].vertices[1] = x+1;
	    cells[x].material_id = 0;
	  }
	break;
	
      case 2:
	cells.resize (repetitions*repetitions);
	for (unsigned int y=0; y<repetitions; ++y)
	  for (unsigned int x=0; x<repetitions; ++x)
	    {
	      const unsigned int c = x  +y*repetitions;
	      cells[c].vertices[0] = x  +y*dy;
	      cells[c].vertices[1] = x+1+y*dy;
	      cells[c].vertices[2] = x  +(y+1)*dy;
	      cells[c].vertices[3] = x+1+(y+1)*dy;
	      cells[c].material_id = 0;
	    }
	break;

      case 3:
	cells.resize (repetitions*repetitions*repetitions);
	for (unsigned int z=0; z<repetitions; ++z)
	  for (unsigned int y=0; y<repetitions; ++y)
	    for (unsigned int x=0; x<repetitions; ++x)
	      {
		const unsigned int c = x+y*repetitions
				       +z*repetitions*repetitions;
		cells[c].vertices[0] = x  +y*dy    +z*dz;
		cells[c].vertices[1] = x+1+y*dy    +z*dz;
		cells[c].vertices[2] = x  +(y+1)*dy+z*dz;
		cells[c].vertices[3] = x+1+(y+1)*dy+z*dz;
		cells[c].vertices[4] = x  +y*dy    +(z+1)*dz;
		cells[c].vertices[5] = x+1+y*dy    +(z+1)*dz;
		cells[c].vertices[6] = x  +(y+1)*dy+(z+1)*dz;
		cells[c].vertices[7] = x+1+(y+1)*dy+(z+1)*dz;
		cells[c].material_id = 0;
	      }
	break;

      default:
                                             // should be trivial to
                                             // do for 3d as well, but
                                             // am too tired at this
                                             // point of the night to
                                             // do that...
                                             //
                                             // contributions are welcome!
            Assert (false, ExcNotImplemented());
    }

  tria.create_triangulation (points, cells, SubCellData());  
}



template <int dim>
void
GridGenerator::subdivided_hyper_rectangle (
  Triangulation<dim>              &tria,
  const std::vector<unsigned int> &repetitions,
  const Point<dim>                &p_1,
  const Point<dim>                &p_2,
  const bool                       colorize)
{
				   // contributed by Joerg R. Weimar
				   // (j.weimar@jweimar.de) 2003
  Assert(repetitions.size() == dim, 
	 ExcInvalidRepetitionsDimension(dim));
				   // First, normalize input such that
				   // p1 is lower in all coordinate
				   // directions.
  Point<dim> p1(p_1);
  Point<dim> p2(p_2);
  
  for (unsigned int i=0;i<dim;++i)
    if (p1(i) > p2(i))
      std::swap (p1(i), p2(i));

				   // then check that all repetitions
				   // are >= 1, and calculate deltas
				   // convert repetitions from double
				   // to int by taking the ceiling.
  Point<dim> delta;
  
  for (unsigned int i=0; i<dim; ++i)
    {
      Assert (repetitions[i] >= 1, ExcInvalidRepetitions(repetitions[i]));
    		
      delta[i] = (p2[i]-p1[i])/repetitions[i];
    }
 
                                   // then generate the necessary
                                   // points
  std::vector<Point<dim> > points;
  switch (dim)
    {
      case 1:
            for (unsigned int x=0; x<=repetitions[0]; ++x)
              points.push_back (Point<dim> (p1[0]+x*delta[0]));
            break;

      case 2:
            for (unsigned int y=0; y<=repetitions[1]; ++y)
              for (unsigned int x=0; x<=repetitions[0]; ++x)
                points.push_back (Point<dim> (p1[0]+x*delta[0],
                                              p1[1]+y*delta[1]));
            break;

      case 3:
            for (unsigned int z=0; z<=repetitions[2]; ++z)
              for (unsigned int y=0; y<=repetitions[1]; ++y)
                for (unsigned int x=0; x<=repetitions[0]; ++x)
                  points.push_back (Point<dim> (p1[0]+x*delta[0],
                                                p1[1]+y*delta[1],
                                                p1[2]+z*delta[2]));
            break;

      default:
            Assert (false, ExcNotImplemented());
    }

                                   // next create the cells
				   // Prepare cell data
  std::vector<CellData<dim> > cells;
  switch (dim)
    {
      case 1:
      {
        cells.resize (repetitions[0]);
        for (unsigned int x=0; x<repetitions[0]; ++x)
          {
            cells[x].vertices[0] = x;
            cells[x].vertices[1] = x+1;
            cells[x].material_id = 0;
          }
        break;
      }
      
      case 2:
      {
        cells.resize (repetitions[1]*repetitions[0]);
        for (unsigned int y=0; y<repetitions[1]; ++y)
          for (unsigned int x=0; x<repetitions[0]; ++x)
            {
              const unsigned int c = x+y*repetitions[0];
              cells[c].vertices[0] = y*(repetitions[0]+1)+x;
              cells[c].vertices[1] = y*(repetitions[0]+1)+x+1;
              cells[c].vertices[2] = (y+1)*(repetitions[0]+1)+x;
              cells[c].vertices[3] = (y+1)*(repetitions[0]+1)+x+1;
              cells[c].material_id = 0;
            }
        break;
      }
      
      case 3:
      {
        const unsigned int n_x  = (repetitions[0]+1);
        const unsigned int n_xy = (repetitions[0]+1)*(repetitions[1]+1);
        
        cells.resize (repetitions[2]*repetitions[1]*repetitions[0]);
        for (unsigned int z=0; z<repetitions[2]; ++z)
          for (unsigned int y=0; y<repetitions[1]; ++y)
            for (unsigned int x=0; x<repetitions[0]; ++x)
              {
                const unsigned int c = x+y*repetitions[0] +
                                       z*repetitions[0]*repetitions[1];
                cells[c].vertices[0] = z*n_xy + y*n_x + x;
                cells[c].vertices[1] = z*n_xy + y*n_x + x+1;
                cells[c].vertices[2] = z*n_xy + (y+1)*n_x + x;
                cells[c].vertices[3] = z*n_xy + (y+1)*n_x + x+1;
                cells[c].vertices[4] = (z+1)*n_xy + y*n_x + x;
                cells[c].vertices[5] = (z+1)*n_xy + y*n_x + x+1;
                cells[c].vertices[6] = (z+1)*n_xy + (y+1)*n_x + x;
                cells[c].vertices[7] = (z+1)*n_xy + (y+1)*n_x + x+1;
                cells[c].material_id = 0;
              }
        break;
        
      }

      default:
            Assert (false, ExcNotImplemented());
    }

  tria.create_triangulation (points, cells, SubCellData());  

  if (colorize)
    {
				       // to colorize, run through all
				       // faces of all cells and set
				       // boundary indicator to the
				       // correct value if it was 0.
     
				       // use a large epsilon to
				       // compare numbers to avoid
				       // roundoff problems.
      const double epsilon
        = 0.01 * *std::min_element (&delta[0], &delta[dim]);
    
                                       // actual code is external since
                                       // 1-D is different from 2/3D.
      colorize_subdivided_hyper_rectangle (tria, p1, p2, epsilon);
    }
}



template <int dim>
void
GridGenerator::
subdivided_hyper_rectangle(Triangulation<dim>              &tria,
                           const std::vector<std::vector<double> > &step_sz,
                           const Point<dim>                &p_1,
                           const Point<dim>                &p_2,
                           const bool                       colorize)
{
				   // contributed by Joerg R. Weimar
				   // (j.weimar@jweimar.de) 2003
				   // modified by Yaqi Wang 2006
  Assert(step_sz.size() == dim, 
	 ExcInvalidRepetitionsDimension(dim));


				   // First, normalize input such that
				   // p1 is lower in all coordinate
				   // directions.

                                   // and check the consistency of
                                   // step sizes, i.e. that they all
                                   // add up to the sizes specified by
                                   // p_1 and p_2
  Point<dim> p1(p_1);
  Point<dim> p2(p_2);
  std::vector< std::vector<double> > step_sizes(step_sz);

  for (unsigned int i=0;i<dim;++i)
    {
      if (p1(i) > p2(i))
	{
	  std::swap (p1(i), p2(i));
	  std::reverse (step_sizes[i].begin(), step_sizes[i].end());
	}

      double x = 0;
      for (unsigned int j=0; j<step_sizes.at(i).size(); j++) 
	x += step_sizes[i][j];
      Assert(std::fabs(x - (p2(i)-p1(i))) <= 1e-12*std::fabs(x),
	     ExcInvalidRepetitions (i) );
    }

 
                                   // then generate the necessary
                                   // points
  std::vector<Point<dim> > points;
  switch (dim)
    {
      case 1:
      {
	double x=0;
	for (unsigned int i=0; ; ++i)
	  { 
	    points.push_back (Point<dim> (p1[0]+x));

					     // form partial sums. in
					     // the last run through
					     // avoid accessing
					     // non-existent values
					     // and exit early instead
	    if (i == step_sizes[0].size())
	      break;
	    
	    x += step_sizes[0][i];
	  }
	break;
      }
      
      case 2:
      {
	double y=0;
	for (unsigned int j=0; ; ++j)
	  {
	    double x=0;
	    for (unsigned int i=0; ; ++i)
	      {
		points.push_back (Point<dim> (p1[0]+x,
					      p1[1]+y));
		if (i == step_sizes[0].size())
		  break;
		
		x += step_sizes[0][i];
	      }
	    
	    if (j == step_sizes[1].size())
	      break;
	    
	    y += step_sizes[1][j];
	  }
	break;

      }
      case 3:
      {
	double z=0;
	for (unsigned int k=0; ; ++k)
	  {
	    double y=0;
	    for (unsigned int j=0; ; ++j)
	      {
		double x=0;
		for (unsigned int i=0; ; ++i)
		  {
		    points.push_back (Point<dim> (p1[0]+x,
						  p1[1]+y,
						  p1[2]+z));
		    if (i == step_sizes[0].size())
		      break;
		    
		    x += step_sizes[0][i];
		  }

		if (j == step_sizes[1].size())
		  break;
		
		y += step_sizes[1][j];
	      }
	    
	    if (k == step_sizes[2].size())
	      break;
	    
	    z += step_sizes[2][k];
	  }
	break;
      }
      
      default:
            Assert (false, ExcNotImplemented());
    }

                                   // next create the cells
				   // Prepare cell data
  std::vector<CellData<dim> > cells;
  switch (dim)
    {
      case 1:
      {
        cells.resize (step_sizes[0].size());
        for (unsigned int x=0; x<step_sizes[0].size(); ++x)
          {
            cells[x].vertices[0] = x;
            cells[x].vertices[1] = x+1;
            cells[x].material_id = 0;
          }
        break;
      }
      
      case 2:
      {
        cells.resize (step_sizes[1].size()*step_sizes[0].size());
        for (unsigned int y=0; y<step_sizes[1].size(); ++y)
          for (unsigned int x=0; x<step_sizes[0].size(); ++x)
            {
              const unsigned int c = x+y*step_sizes[0].size();
              cells[c].vertices[0] = y*(step_sizes[0].size()+1)+x;
              cells[c].vertices[1] = y*(step_sizes[0].size()+1)+x+1;
              cells[c].vertices[2] = (y+1)*(step_sizes[0].size()+1)+x;
              cells[c].vertices[3] = (y+1)*(step_sizes[0].size()+1)+x+1;
              cells[c].material_id = 0;
            }
        break;
      }
      
      case 3:
      {
        const unsigned int n_x  = (step_sizes[0].size()+1);
        const unsigned int n_xy = (step_sizes[0].size()+1)*(step_sizes[1].size()+1);
        
        cells.resize (step_sizes[2].size()*step_sizes[1].size()*step_sizes[0].size());
        for (unsigned int z=0; z<step_sizes[2].size(); ++z)
          for (unsigned int y=0; y<step_sizes[1].size(); ++y)
            for (unsigned int x=0; x<step_sizes[0].size(); ++x)
              {
                const unsigned int c = x+y*step_sizes[0].size() +
                                       z*step_sizes[0].size()*step_sizes[1].size();
                cells[c].vertices[0] = z*n_xy + y*n_x + x;
                cells[c].vertices[1] = z*n_xy + y*n_x + x+1;
                cells[c].vertices[2] = z*n_xy + (y+1)*n_x + x;
                cells[c].vertices[3] = z*n_xy + (y+1)*n_x + x+1;
                cells[c].vertices[4] = (z+1)*n_xy + y*n_x + x;
                cells[c].vertices[5] = (z+1)*n_xy + y*n_x + x+1;
                cells[c].vertices[6] = (z+1)*n_xy + (y+1)*n_x + x;
                cells[c].vertices[7] = (z+1)*n_xy + (y+1)*n_x + x+1;
                cells[c].material_id = 0;
              }
        break;
        
      }

      default:
            Assert (false, ExcNotImplemented());
    }

  tria.create_triangulation (points, cells, SubCellData());  

  if (colorize)
    {
				       // to colorize, run through all
				       // faces of all cells and set
				       // boundary indicator to the
				       // correct value if it was 0.
     
				       // use a large epsilon to
				       // compare numbers to avoid
				       // roundoff problems.
      double min_size = *std::min_element (step_sizes[0].begin(),
					   step_sizes[0].end());
      for (unsigned int i=1; i<dim; ++i)
	min_size = std::min (min_size,
			     *std::min_element (step_sizes[i].begin(),
						step_sizes[i].end()));
      const double epsilon = 0.01 * min_size;
    
                                       // actual code is external since
                                       // 1-D is different from 2/3D.
      colorize_subdivided_hyper_rectangle (tria, p1, p2, epsilon);
    }
}



#if deal_II_dimension == 1

template <>
void
GridGenerator::
subdivided_hyper_rectangle (Triangulation<1>&                             tria,
                            const std::vector< std::vector<double> >&     spacing,
                            const Point<1>&                               p,
			    const Table<1,unsigned char>&                 material_id,
                            const bool                                    colorize)
{
				   // contributed by Yaqi Wang 2006
  Assert(spacing.size() == 1, 
	 ExcInvalidRepetitionsDimension(1));

  const unsigned int n_cells = material_id.size(0);

  Assert(spacing[0].size() == n_cells, 
	 ExcInvalidRepetitionsDimension(1));

  double delta = std::numeric_limits<double>::max();
  for (unsigned int i=0; i<n_cells; i++) {
    Assert (spacing[0][i] >= 0, ExcInvalidRepetitions(-1));
    delta = std::min (delta, spacing[0][i]);
  }
  
                                   // generate the necessary points
  std::vector<Point<1> > points;
  double ax = p[0];
  for (unsigned int x=0; x<=n_cells; ++x) {
    points.push_back (Point<1> (ax));
    if (x<n_cells)
      ax += spacing[0][x];
  }
                                   // create the cells
  unsigned int n_val_cells = 0;
  for (unsigned int i=0; i<n_cells; i++)
    if (material_id[i]!=(unsigned char)(-1)) n_val_cells++;

  std::vector<CellData<1> > cells(n_val_cells);
  unsigned int id = 0;
  for (unsigned int x=0; x<n_cells; ++x)
    if (material_id[x] != (unsigned char)(-1))
      {
	cells[id].vertices[0] = x;
	cells[id].vertices[1] = x+1;
	cells[id].material_id = material_id[x];
	id++;
      }
                                   // create triangulation
  SubCellData t;
  GridTools::delete_unused_vertices (points, cells, t);

  tria.create_triangulation (points, cells, t);  

                                   // set boundary indicator
  if (colorize)
    Assert (false, ExcNotImplemented());
}

#endif

#if deal_II_dimension == 2

template <>
void
GridGenerator::
subdivided_hyper_rectangle (Triangulation<2>&                         tria,
                            const std::vector< std::vector<double> >&     spacing,
                            const Point<2>&                               p,
			    const Table<2,unsigned char>&                 material_id,
                            const bool                                    colorize)
{
				   // contributed by Yaqi Wang 2006
  Assert(spacing.size() == 2, 
	 ExcInvalidRepetitionsDimension(2));

  std::vector<unsigned int> repetitions(2);
  unsigned int n_cells = 1;
  double delta = std::numeric_limits<double>::max();
  for (unsigned int i=0; i<2; i++)
    {
      repetitions[i] = spacing[i].size();
      n_cells *= repetitions[i];
      for (unsigned int j=0; j<repetitions[i]; j++)
	{
	  Assert (spacing[i][j] >= 0, ExcInvalidRepetitions(-1));
	  delta = std::min (delta, spacing[i][j]);
	}
      Assert(material_id.size(i) == repetitions[i], 
	     ExcInvalidRepetitionsDimension(i));
    }
 
                                   // generate the necessary points
  std::vector<Point<2> > points;
  double ay = p[1];
  for (unsigned int y=0; y<=repetitions[1]; ++y)
    {
      double ax = p[0];
      for (unsigned int x=0; x<=repetitions[0]; ++x)
	{
	  points.push_back (Point<2> (ax,ay));
	  if (x<repetitions[0])
	    ax += spacing[0][x];
	}
      if (y<repetitions[1])
	ay += spacing[1][y];
    }

                                   // create the cells
  unsigned int n_val_cells = 0;
  for (unsigned int i=0; i<material_id.size(0); i++)
    for (unsigned int j=0; j<material_id.size(1); j++)
      if (material_id[i][j] != (unsigned char)(-1))
	n_val_cells++;

  std::vector<CellData<2> > cells(n_val_cells);
  unsigned int id = 0;
  for (unsigned int y=0; y<repetitions[1]; ++y)
    for (unsigned int x=0; x<repetitions[0]; ++x)
      if (material_id[x][y]!=(unsigned char)(-1))
	{
	  cells[id].vertices[0] = y*(repetitions[0]+1)+x;
	  cells[id].vertices[1] = y*(repetitions[0]+1)+x+1;
	  cells[id].vertices[2] = (y+1)*(repetitions[0]+1)+x;
	  cells[id].vertices[3] = (y+1)*(repetitions[0]+1)+x+1;
	  cells[id].material_id = material_id[x][y];
	  id++;
	}

                                    // create triangulation
  SubCellData t;
  GridTools::delete_unused_vertices (points, cells, t);

  tria.create_triangulation (points, cells, t);  

                                    // set boundary indicator
  if (colorize)
    {
      double eps = 0.01 * delta;
      Triangulation<2>::cell_iterator cell = tria.begin_raw(),
				      endc = tria.end();
      for (; cell !=endc; ++cell)
	{
	  Point<2> cell_center = cell->center();
	  for(unsigned int f=0; f<GeometryInfo<2>::faces_per_cell; ++f)
	    if (cell->face(f)->boundary_indicator() == 0)
	      {
		Point<2> face_center = cell->face(f)->center();
		for (unsigned int i=0; i<2; ++i)
		  {
		    if (face_center[i]<cell_center[i]-eps)
		      cell->face(f)->set_boundary_indicator(i*2);
		    if (face_center[i]>cell_center[i]+eps)
		      cell->face(f)->set_boundary_indicator(i*2+1);
		  }
	      }
	}
    }
}

#endif

#if deal_II_dimension == 3

template <>
void
GridGenerator::
subdivided_hyper_rectangle (Triangulation<3>&                           tria,
                            const std::vector< std::vector<double> >&     spacing,
                            const Point<3>&                             p,
			    const Table<3,unsigned char>&               material_id,
                            const bool                                    colorize)
{
				   // contributed by Yaqi Wang 2006
  const unsigned int dim = 3;
  
  Assert(spacing.size() == dim, 
	 ExcInvalidRepetitionsDimension(dim));

  std::vector<unsigned int> repetitions(dim);
  unsigned int n_cells = 1;
  double delta = std::numeric_limits<double>::max();
  for (unsigned int i=0; i<dim; i++)
    {
      repetitions[i] = spacing[i].size();
      n_cells *= repetitions[i];
      for (unsigned int j=0; j<repetitions[i]; j++)
	{
	  Assert (spacing[i][j] >= 0, ExcInvalidRepetitions(-1));
	  delta = std::min (delta, spacing[i][j]);
	}
      Assert(material_id.size(i) == repetitions[i], 
	     ExcInvalidRepetitionsDimension(i));
    }

                                   // generate the necessary points
  std::vector<Point<dim> > points;
  double az = p[2];
  for (unsigned int z=0; z<=repetitions[2]; ++z)
    {
      double ay = p[1];
      for (unsigned int y=0; y<=repetitions[1]; ++y)
	{
	  double ax = p[0];
	  for (unsigned int x=0; x<=repetitions[0]; ++x)
	    {
	      points.push_back (Point<dim> (ax,ay,az));
	      if (x<repetitions[0])
		ax += spacing[0][x];
	    }
	  if (y<repetitions[1])
	    ay += spacing[1][y];
	}
      if (z<repetitions[2])
	az += spacing[2][z];
    }

                                   // create the cells
  unsigned int n_val_cells = 0;
  for (unsigned int i=0; i<material_id.size(0); i++)
    for (unsigned int j=0; j<material_id.size(1); j++)
      for (unsigned int k=0; k<material_id.size(2); k++)
	if (material_id[i][j][k]!=(unsigned char)(-1))
	  n_val_cells++;

  std::vector<CellData<dim> > cells(n_val_cells);
  unsigned int id = 0;
  const unsigned int n_x  = (repetitions[0]+1);
  const unsigned int n_xy = (repetitions[0]+1)*(repetitions[1]+1);
  for (unsigned int z=0; z<repetitions[2]; ++z)
    for (unsigned int y=0; y<repetitions[1]; ++y)
      for (unsigned int x=0; x<repetitions[0]; ++x)
	if (material_id[x][y][z]!=(unsigned char)(-1))
	  {
	    cells[id].vertices[0] = z*n_xy + y*n_x + x;
	    cells[id].vertices[1] = z*n_xy + y*n_x + x+1;
	    cells[id].vertices[2] = z*n_xy + (y+1)*n_x + x;
	    cells[id].vertices[3] = z*n_xy + (y+1)*n_x + x+1;
	    cells[id].vertices[4] = (z+1)*n_xy + y*n_x + x;
	    cells[id].vertices[5] = (z+1)*n_xy + y*n_x + x+1;
	    cells[id].vertices[6] = (z+1)*n_xy + (y+1)*n_x + x;
	    cells[id].vertices[7] = (z+1)*n_xy + (y+1)*n_x + x+1;
	    cells[id].material_id = material_id[x][y][z];
	    id++;
	  }

				 // create triangulation
  SubCellData t;
  GridTools::delete_unused_vertices (points, cells, t);

  tria.create_triangulation (points, cells, t);  

                                  // set boundary indicator
  if (colorize)
    {
      double eps = 0.01 * delta;
      Triangulation<dim>::cell_iterator cell = tria.begin_raw(),
					endc = tria.end();
      for (; cell !=endc; ++cell)
	{
	  Point<dim> cell_center = cell->center();
	  for(unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
	    if (cell->face(f)->boundary_indicator() == 0)
	      {
		Point<dim> face_center = cell->face(f)->center();
		for (unsigned int i=0; i<dim; ++i)
		  {
		    if (face_center[i]<cell_center[i]-eps)
		      cell->face(f)->set_boundary_indicator(i*2);
		    if (face_center[i]>cell_center[i]+eps)
		      cell->face(f)->set_boundary_indicator(i*2+1);
		  }
	      }
	}
    }
}

#endif


#if deal_II_dimension == 1

// Implementation for 1D only
template <int dim>
void
GridGenerator::colorize_subdivided_hyper_rectangle (
  Triangulation<dim>& tria,
  const Point<dim>&,
  const Point<dim>&,
  const double)
{
  for (typename Triangulation<dim>::cell_iterator cell = tria.begin();
       cell != tria.end(); ++cell)
    if (cell->center()(0) > 0)
      cell->set_material_id(1);
				   // boundary indicators are set to
				   // 0 (left) and 1 (right) by default.
}

#else

// Implementation for dimensions except 1
template <int dim>
void
GridGenerator::colorize_subdivided_hyper_rectangle (Triangulation<dim> &tria,
						    const Point<dim>   &p1,
						    const Point<dim>   &p2,
						    const double        epsilon)
{
	
				   // run through all faces and check
				   // if one of their center coordinates matches
				   // one of the corner points. Comparisons
				   // are made using an epsilon which
				   // should be smaller than the smallest cell
				   // diameter.
				    
  typename Triangulation<dim>::raw_face_iterator face = tria.begin_raw_face(),
					      endface = tria.end_face();
  for (; face!=endface; ++face)
    {
      if (face->boundary_indicator() == 0)
	{
	  const Point<dim> center (face->center());
	  if (std::abs(center(0)-p1[0]) < epsilon)
	    face->set_boundary_indicator(0);
	  else if (std::abs(center(0) - p2[0]) < epsilon)
	    face->set_boundary_indicator(1);
	  else if (dim > 1 && std::abs(center(1) - p1[1]) < epsilon)
	    face->set_boundary_indicator(2);
	  else if (dim > 1 && std::abs(center(1) - p2[1]) < epsilon)
	    face->set_boundary_indicator(3);
	  else if (dim > 2 && std::abs(center(2) - p1[2]) < epsilon)
	    face->set_boundary_indicator(4);
	  else if (dim > 2 && std::abs(center(2) - p2[2]) < epsilon)
	    face->set_boundary_indicator(5);
	  else
					     // triangulation says it
					     // is on the boundary,
					     // but we could not find
					     // on which boundary.
	    Assert (false, ExcInternalError());
	  
	}
    }
  for (typename Triangulation<dim>::cell_iterator cell = tria.begin();
       cell != tria.end(); ++cell)
    {
      char id = 0;
      for (unsigned int d=0;d<dim;++d)
	if (cell->center()(d) > 0) id += 1 << d;
      cell->set_material_id(id);
    }
}

#endif



#if deal_II_dimension == 1

// The following implementations are for 1D only. See below for other
// dimensions.
template <int dim>
void GridGenerator::hyper_cube_slit (Triangulation<dim> &,
				     const double,
				     const double,
				     const bool)
{
  Assert (false, ExcNotImplemented());
}



// Implementation for 1D only
template <int dim>
void GridGenerator::enclosed_hyper_cube (Triangulation<dim>&,
					 const double,
					 const double,
					 const double,
					 const bool)
{
  Assert (false, ExcNotImplemented());
}



// Implementation for 1D only
template <int dim>
void GridGenerator::hyper_L (Triangulation<dim> &,
			     const double,
			     const double)
{
  Assert (false, ExcNotImplemented());
}



// Implementation for 1D only
template <int dim>
void GridGenerator::hyper_ball (Triangulation<dim> &,
				const Point<dim> &,
				const double)
{
  Assert (false, ExcNotImplemented());
}



// Implementation for 1D only
template <int dim>
void GridGenerator::cylinder (Triangulation<dim> &,
			      const double,
			      const double)
{
  Assert (false, ExcNotImplemented());
}



// Implementation for 1D only
template <int dim>
void GridGenerator::hyper_shell (Triangulation<dim> &,
				 const Point<dim> &,
				 const double,
				 const double,
				 const unsigned int,
				 const bool)
{
  Assert (false, ExcNotImplemented());
}

// Implementation for 1D only
template <int dim>
void GridGenerator::colorize_hyper_shell (Triangulation<dim> &,
					  const Point<dim> &,
					  const double,
					  const double)
{
  Assert (false, ExcNotImplemented());
}

// Implementation for 1D only
template <int dim>
void GridGenerator::cylinder_shell (Triangulation<dim>&,
                                    const double,
                                    const double,
                                    const double,
                                    const unsigned int,
                                    const unsigned int)
{
  Assert (false, ExcNotImplemented());
}


// Implementation for 1D only
template <int dim>
void
GridGenerator::half_hyper_ball (Triangulation<dim>&,
				const Point<dim>&,
				const double)
{
  Assert (false, ExcNotImplemented());
}


// Implementation for 1D only
template <int dim>
void
GridGenerator::half_hyper_shell (Triangulation<dim>&,
				 const Point<dim>&,
				 const double,
				 const double,
				 const unsigned int)
{
  Assert (false, ExcNotImplemented());
}

#endif



#if deal_II_dimension == 2

// Implementation for 2D only
template <int dim>
void GridGenerator::enclosed_hyper_cube (Triangulation<dim> &tria,
					 const double      l,
					 const double      r,
					 const double      d,
					 const bool        colorize)
{
  std::vector<Point<dim> > vertices(16);
  double coords[4];
  coords[0] = l-d;
  coords[1] = l;
  coords[2] = r;
  coords[3] = r+d;

  unsigned int k=0;
  for (unsigned int i0=0;i0<4;++i0)
    for (unsigned int i1=0;i1<4;++i1)
      vertices[k++] = Point<dim>(coords[i1], coords[i0]);

  const unsigned char materials[9] = { 5, 4, 6,
				       1, 0, 2,
				       9, 8,10
  };
  
  std::vector<CellData<dim> > cells(9);
  k = 0;
  for (unsigned int i0=0;i0<3;++i0)
    for (unsigned int i1=0;i1<3;++i1)
      {
	cells[k].vertices[0] = i1+4*i0;
	cells[k].vertices[1] = i1+4*i0+1;
	cells[k].vertices[2] = i1+4*i0+4;
	cells[k].vertices[3] = i1+4*i0+5;
	if (colorize)
	  cells[k].material_id = materials[k];
	++k;
      }
  tria.create_triangulation (vertices,
			     cells,
			     SubCellData());       // no boundary information
}



// Implementation for 2D only
template <int dim>
void
GridGenerator::hyper_cube_slit (Triangulation<dim> &tria,
				const double left,
				const double right,
				const bool colorize)
{
  const double rl2=(right+left)/2;
  const Point<dim> vertices[10] = { Point<dim>(left, left ),
				  Point<dim>(rl2,  left ),
				  Point<dim>(rl2,  rl2  ),
				  Point<dim>(left, rl2  ),
				  Point<dim>(right,left ),
				  Point<dim>(right,rl2  ),
				  Point<dim>(rl2,  right),
				  Point<dim>(left, right),
				  Point<dim>(right,right),
				  Point<dim>(rl2,  left ) };
  const int cell_vertices[4][4] = { { 0,1,3,2 },
				    { 9,4,2,5 },
				    { 3,2,7,6 },
				    { 2,5,6,8 } };
  std::vector<CellData<dim> > cells (4, CellData<dim>());
  for (unsigned int i=0; i<4; ++i)
    {
      for (unsigned int j=0; j<4; ++j)
	cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    };
  tria.create_triangulation (
    std::vector<Point<dim> >(&vertices[0], &vertices[10]),
    cells,
    SubCellData());       // no boundary information

  if (colorize)
    {
      typename Triangulation<dim>::cell_iterator cell = tria.begin();
      cell->face(1)->set_boundary_indicator(1);
      ++cell;
      cell->face(3)->set_boundary_indicator(2);
    }
}


//TODO: Colorize edges as circumference, left and right radius
// Implementation for 2D only
template <int dim>
void
GridGenerator::hyper_L (Triangulation<dim> &tria,
			const double a,
			const double b)
{
  const Point<dim> vertices[8] = { Point<dim> (a,a),
				   Point<dim> ((a+b)/2,a),
				   Point<dim> (b,a),
				   Point<dim> (a,(a+b)/2),
				   Point<dim> ((a+b)/2,(a+b)/2),
				   Point<dim> (b,(a+b)/2),
				   Point<dim> (a,b),
				   Point<dim> ((a+b)/2,b)  };
  const int cell_vertices[3][4] = {{0, 1, 3, 4},
				   {1, 2, 4, 5},
				   {3, 4, 6, 7}};

  std::vector<CellData<dim> > cells (3, CellData<dim>());
  
  for (unsigned int i=0; i<3; ++i) 
    {
      for (unsigned int j=0; j<4; ++j)
	cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    };
  
  tria.create_triangulation (
    std::vector<Point<dim> >(&vertices[0], &vertices[8]),
    cells,
    SubCellData());       // no boundary information
}



// Implementation for 2D only
template <int dim>
void
GridGenerator::hyper_ball (Triangulation<dim> &tria,
			   const Point<dim>   &p,
			   const double      radius)
{
				   // equilibrate cell sizes at
				   // transition from the inner part
				   // to the radial cells
  const double a = 1./(1+std::sqrt(2.0));
  const Point<dim> vertices[8] = { p+Point<dim>(-1,-1)*(radius/std::sqrt(2.0)),
				 p+Point<dim>(+1,-1)*(radius/std::sqrt(2.0)),
				 p+Point<dim>(-1,-1)*(radius/std::sqrt(2.0)*a),
				 p+Point<dim>(+1,-1)*(radius/std::sqrt(2.0)*a),
				 p+Point<dim>(-1,+1)*(radius/std::sqrt(2.0)*a),
				 p+Point<dim>(+1,+1)*(radius/std::sqrt(2.0)*a),
				 p+Point<dim>(-1,+1)*(radius/std::sqrt(2.0)),
				 p+Point<dim>(+1,+1)*(radius/std::sqrt(2.0)) };
  
  const int cell_vertices[5][4] = {{0, 1, 2, 3},
				   {0, 2, 6, 4},
				   {2, 3, 4, 5},
				   {1, 7, 3, 5},
				   {6, 4, 7, 5}};

  std::vector<CellData<dim> > cells (5, CellData<dim>());
  
  for (unsigned int i=0; i<5; ++i) 
    {
      for (unsigned int j=0; j<4; ++j)
	cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    };
  
  tria.create_triangulation (
    std::vector<Point<dim> >(&vertices[0], &vertices[8]),
    cells,
    SubCellData());       // no boundary information
}



// Implementation for 2D only
template <int dim>
void GridGenerator::hyper_shell (Triangulation<dim>   &tria,
				 const Point<dim>     &center,
				 const double        inner_radius,
				 const double        outer_radius,
				 const unsigned int  n_cells,
				 const bool colorize)
{
  Assert ((inner_radius > 0) && (inner_radius < outer_radius),
	  ExcInvalidRadii ());
  
  const double pi = numbers::PI;
  
				   // determine the number of cells
				   // for the grid. if not provided by
				   // the user determine it such that
				   // the length of each cell on the
				   // median (in the middle between
				   // the two circles) is equal to its
				   // radial extent (which is the
				   // difference between the two
				   // radii)
  const unsigned int N = (n_cells == 0 ?
			  static_cast<unsigned int>
			  (std::ceil((2*pi* (outer_radius + inner_radius)/2) /
				     (outer_radius - inner_radius))) :
			  n_cells);

				   // set up N vertices on the
				   // outer and N vertices on
				   // the inner circle. the
				   // first N ones are on the
				   // outer one, and all are
				   // numbered counter-clockwise
  std::vector<Point<dim> > vertices(2*N);
  for (unsigned int i=0; i<N; ++i)
    {
      vertices[i] = Point<dim>( std::cos(2*pi * i/N),
			      std::sin(2*pi * i/N)) * outer_radius;
      vertices[i+N] = vertices[i] * (inner_radius/outer_radius);

      vertices[i]   += center;
      vertices[i+N] += center;
    };

  std::vector<CellData<dim> > cells (N, CellData<dim>());
	
  for (unsigned int i=0; i<N; ++i) 
    {
      cells[i].vertices[0] = i;
      cells[i].vertices[1] = (i+1)%N;
      cells[i].vertices[2] = N+i;
      cells[i].vertices[3] = N+((i+1)%N);
	    
      cells[i].material_id = 0;
    };
  
  tria.create_triangulation (
    vertices, cells, SubCellData());

  if (colorize)
    colorize_hyper_shell(tria, center, inner_radius, outer_radius);
}


template<int dim>
void
GridGenerator::colorize_hyper_shell (
  Triangulation<dim>& tria,
  const Point<dim>&, const double, const double)
{
				   // Inspite of receiving geometrical
				   // data, we do this only based on
				   // topology.

				   // For the mesh based on  cube,
				   // this is highly irregular
  for (typename Triangulation<dim>::cell_iterator cell = tria.begin();
       cell != tria.end(); ++cell)
    {
      cell->face(2)->set_boundary_indicator(1);
    }
}



// Implementation for 2D only
template <int dim>
void
GridGenerator::cylinder (Triangulation<dim> &tria,
			 const double radius,
			 const double half_length)
{
  Point<dim> p1 (-half_length, -radius);
  Point<dim> p2 (half_length, radius);

  hyper_rectangle(tria, p1, p2, true);

  typename Triangulation<dim>::face_iterator f = tria.begin_face();
  typename Triangulation<dim>::face_iterator end = tria.end_face();
  while (f != end)
    {
      switch (f->boundary_indicator())
	{
	  case 0:
		f->set_boundary_indicator(1);
		break;
	  case 1:
		f->set_boundary_indicator(2);
		break;
	  default:
		f->set_boundary_indicator(0);
		break;	    
	}
      ++f;
    }
}



// Implementation for 2D only
template <int dim>
void GridGenerator::cylinder_shell (Triangulation<dim>&,
                                    const double,
                                    const double,
                                    const double,
                                    const unsigned int,
                                    const unsigned int)
{
  Assert (false, ExcNotImplemented());
}


template <int dim>
void
GridGenerator::half_hyper_ball (Triangulation<dim> &tria,
				const Point<dim>   &p,
				const double      radius)
{
				   // equilibrate cell sizes at
				   // transition from the inner part
				   // to the radial cells
  const double a = 1./(1+std::sqrt(2.0));
  const Point<dim> vertices[8] = { p+Point<dim>(0,-1)*radius,
				 p+Point<dim>(+1,-1)*(radius/std::sqrt(2.0)),
				 p+Point<dim>(0,-1)*(radius/std::sqrt(2.0)*a),
				 p+Point<dim>(+1,-1)*(radius/std::sqrt(2.0)*a),
				 p+Point<dim>(0,+1)*(radius/std::sqrt(2.0)*a),
				 p+Point<dim>(+1,+1)*(radius/std::sqrt(2.0)*a),
				 p+Point<dim>(0,+1)*radius,
				 p+Point<dim>(+1,+1)*(radius/std::sqrt(2.0)) };
  
  const int cell_vertices[5][4] = {{0, 1, 2, 3},
				   {2, 3, 4, 5},
				   {1, 7, 3, 5},
				   {6, 4, 7, 5}};

  std::vector<CellData<dim> > cells (4, CellData<dim>());
  
  for (unsigned int i=0; i<4; ++i) 
    {
      for (unsigned int j=0; j<4; ++j)
	cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    };
  
  tria.create_triangulation (
    std::vector<Point<dim> >(&vertices[0], &vertices[8]),
    cells,
    SubCellData());       // no boundary information

    typename Triangulation<dim>::cell_iterator cell = tria.begin();
    typename Triangulation<dim>::cell_iterator end = tria.end();


    while (cell != end)
    {
	for (unsigned int i=0;i<GeometryInfo<dim>::faces_per_cell;++i)
	{
	    if (cell->face(i)->boundary_indicator() == 255)
		continue;

	    // If x is zero, then this is part of the plane
	    if (cell->face(i)->center()(0) < p(0)+1.e-5)
		cell->face(i)->set_boundary_indicator(1);
	}
	++cell;
    }
}



// Implementation for 2D only
template <int dim>
void
GridGenerator::half_hyper_shell (Triangulation<dim>   &tria,
				 const Point<dim>     &center,
				 const double        inner_radius,
				 const double        outer_radius,
				 const unsigned int  n_cells)
{
  Assert ((inner_radius > 0) && (inner_radius < outer_radius),
	  ExcInvalidRadii ());
  
  const double pi     = numbers::PI;
				   // determine the number of cells
				   // for the grid. if not provided by
				   // the user determine it such that
				   // the length of each cell on the
				   // median (in the middle between
				   // the two circles) is equal to its
				   // radial extent (which is the
				   // difference between the two
				   // radii)
  const unsigned int N = (n_cells == 0 ?
			  static_cast<unsigned int>
			  (std::ceil((pi* (outer_radius + inner_radius)/2) /
				     (outer_radius - inner_radius))) :
			  n_cells);

				   // set up N+1 vertices on the
				   // outer and N+1 vertices on
				   // the inner circle. the
				   // first N+1 ones are on the
				   // outer one, and all are
				   // numbered counter-clockwise
  std::vector<Point<dim> > vertices(2*(N+1));
  for (unsigned int i=0; i<=N; ++i)
    {
				       // enforce that the x-coordinates
				       // of the first and last point of
				       // each half-circle are exactly
				       // zero (contrary to what we may
				       // compute using the imprecise
				       // value of pi)
      vertices[i] =  Point<dim>( ( (i==0) || (i==N) ?
				 0 :
				 std::cos(pi * i/N - pi/2) ),
			       std::sin(pi * i/N - pi/2)) * outer_radius;
      vertices[i+N+1] = vertices[i] * (inner_radius/outer_radius);

      vertices[i]     += center;
      vertices[i+N+1] += center;
    };


  std::vector<CellData<dim> > cells (N, CellData<dim>());
	
  for (unsigned int i=0; i<N; ++i) 
    {
      cells[i].vertices[0] = i;
      cells[i].vertices[1] = (i+1)%(N+1);
      cells[i].vertices[2] = N+1+i;
      cells[i].vertices[3] = N+1+((i+1)%(N+1));
	    
      cells[i].material_id = 0;
    };
  
  tria.create_triangulation (vertices, cells, SubCellData());
}



#endif


#if deal_II_dimension == 3


// Implementation for 3D only
template <int dim>
void GridGenerator::hyper_cube_slit (Triangulation<dim>& tria,
				     const double left,
				     const double right,
				     const bool colorize)
{
  const double rl2=(right+left)/2;
  const double len = (right-left)/2.;
  
  const Point<dim> vertices[20] = {
	Point<dim>(left, left , -len/2.),
	Point<dim>(rl2,  left , -len/2.),
	Point<dim>(rl2,  rl2  , -len/2.),
	Point<dim>(left, rl2  , -len/2.),
	Point<dim>(right,left , -len/2.),
	Point<dim>(right,rl2  , -len/2.),
	Point<dim>(rl2,  right, -len/2.),
	Point<dim>(left, right, -len/2.),
	Point<dim>(right,right, -len/2.),
	Point<dim>(rl2,  left , -len/2.),
	Point<dim>(left, left , len/2.),
	Point<dim>(rl2,  left , len/2.),
	Point<dim>(rl2,  rl2  , len/2.),
	Point<dim>(left, rl2  , len/2.),
	Point<dim>(right,left , len/2.),
	Point<dim>(right,rl2  , len/2.),
	Point<dim>(rl2,  right, len/2.),
	Point<dim>(left, right, len/2.),
	Point<dim>(right,right, len/2.),
	Point<dim>(rl2,  left , len/2.)
  };
  const int cell_vertices[4][8] = { { 0,1,3,2, 10, 11, 13, 12 },
				    { 9,4,2,5, 19,14, 12, 15 },
				    { 3,2,7,6,13,12,17,16 },
				    { 2,5,6,8,12,15,16,18 } };
  std::vector<CellData<dim> > cells (4, CellData<dim>());
  for (unsigned int i=0; i<4; ++i)
    {
      for (unsigned int j=0; j<8; ++j)
	cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    };
  tria.create_triangulation (
    std::vector<Point<dim> >(&vertices[0], &vertices[20]),
    cells,
    SubCellData());       // no boundary information

  if (colorize)
    {
      Assert(false, ExcNotImplemented());
      typename Triangulation<dim>::cell_iterator cell = tria.begin();
      cell->face(1)->set_boundary_indicator(1);
      ++cell;
      cell->face(3)->set_boundary_indicator(2);
    }
}



// Implementation for 3D only
template <int dim>
void GridGenerator::enclosed_hyper_cube (Triangulation<dim> &tria,
					 const double      l,
					 const double      r,
					 const double      d,
					 const bool        colorize)
{
  std::vector<Point<dim> > vertices(64);
  double coords[4];
  coords[0] = l-d;
  coords[1] = l;
  coords[2] = r;
  coords[3] = r+d;

  unsigned int k=0;
  for (unsigned int z=0;z<4;++z)
    for (unsigned int y=0;y<4;++y)
      for (unsigned int x=0;x<4;++x)
	vertices[k++] = Point<dim>(coords[x], coords[y], coords[z]);

  const unsigned char materials[27] = {
	21,20,22,
	17,16,18,
	25,24,26,
	5 , 4, 6,
	1 , 0, 2,
	9 , 8,10,
	37,36,38,
	33,32,34,
	41,40,42
  };
  
  std::vector<CellData<dim> > cells(27);
  k = 0;
  for (unsigned int z=0;z<3;++z)
    for (unsigned int y=0;y<3;++y)
      for (unsigned int x=0;x<3;++x)
	{
	  cells[k].vertices[0] = x+4*y+16*z;
	  cells[k].vertices[1] = x+4*y+16*z+1;
	  cells[k].vertices[2] = x+4*y+16*z+4;
	  cells[k].vertices[3] = x+4*y+16*z+5;
	  cells[k].vertices[4] = x+4*y+16*z+16;
	  cells[k].vertices[5] = x+4*y+16*z+17;
	  cells[k].vertices[6] = x+4*y+16*z+20;
	  cells[k].vertices[7] = x+4*y+16*z+21;
	  if (colorize)
	    cells[k].material_id = materials[k];
	  ++k;
	}
  tria.create_triangulation (
    vertices,
    cells,
    SubCellData());       // no boundary information
}



// Implementation for 3D only
template <int dim>
void
GridGenerator::hyper_L (Triangulation<dim> &tria,
			const double      a,
			const double      b)
{
				   // we slice out the top back right
				   // part of the cube
  const Point<dim> vertices[26]
    = {
					   // front face of the big cube
	  Point<dim> (a,      a,a),
	  Point<dim> ((a+b)/2,a,a),
	  Point<dim> (b,      a,a),
	  Point<dim> (a,      a,(a+b)/2),
	  Point<dim> ((a+b)/2,a,(a+b)/2),
	  Point<dim> (b,      a,(a+b)/2),
	  Point<dim> (a,      a,b),
	  Point<dim> ((a+b)/2,a,b),
	  Point<dim> (b,      a,b),
					   // middle face of the big cube
	  Point<dim> (a,      (a+b)/2,a),
	  Point<dim> ((a+b)/2,(a+b)/2,a),
	  Point<dim> (b,      (a+b)/2,a),
	  Point<dim> (a,      (a+b)/2,(a+b)/2),
	  Point<dim> ((a+b)/2,(a+b)/2,(a+b)/2),
	  Point<dim> (b,      (a+b)/2,(a+b)/2),
	  Point<dim> (a,      (a+b)/2,b),
	  Point<dim> ((a+b)/2,(a+b)/2,b),
	  Point<dim> (b,      (a+b)/2,b),
					   // back face of the big cube
					   // last (top right) point is missing
	  Point<dim> (a,      b,a),
	  Point<dim> ((a+b)/2,b,a),
	  Point<dim> (b,      b,a),
	  Point<dim> (a,      b,(a+b)/2),
	  Point<dim> ((a+b)/2,b,(a+b)/2),
	  Point<dim> (b,      b,(a+b)/2),
	  Point<dim> (a,      b,b),
	  Point<dim> ((a+b)/2,b,b)
    };
  const int cell_vertices[7][8] = {{0, 1, 9, 10, 3, 4, 12, 13},
				   {1, 2, 10, 11, 4, 5, 13, 14},
				   {3, 4, 12, 13, 6, 7, 15, 16},
				   {4, 5, 13, 14, 7, 8, 16, 17},
				   {9, 10, 18, 19, 12, 13, 21, 22},
				   {10, 11, 19, 20, 13, 14, 22, 23},
				   {12, 13, 21, 22, 15, 16, 24, 25}};

  std::vector<CellData<dim> > cells (7, CellData<dim>());
  
  for (unsigned int i=0; i<7; ++i) 
    {
      for (unsigned int j=0; j<8; ++j)
	cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    };

  tria.create_triangulation (
    std::vector<Point<dim> >(&vertices[0], &vertices[26]),
    cells,
    SubCellData());       // no boundary information
}



// Implementation for 3D only
template <int dim>
void
GridGenerator::hyper_ball (Triangulation<dim> &tria,
			   const Point<dim>   &p,
			   const double radius)
{
  const double a = 1./(1+std::sqrt(3.0)); // equilibrate cell sizes at transition
				          // from the inner part to the radial
				          // cells
  const unsigned int n_vertices = 16;
  const Point<dim> vertices[n_vertices]
    = {
					   // first the vertices of the inner
					   // cell
	  p+Point<dim>(-1,-1,-1)*(radius/std::sqrt(3.0)*a),
	  p+Point<dim>(+1,-1,-1)*(radius/std::sqrt(3.0)*a),
	  p+Point<dim>(+1,-1,+1)*(radius/std::sqrt(3.0)*a),
	  p+Point<dim>(-1,-1,+1)*(radius/std::sqrt(3.0)*a),
	  p+Point<dim>(-1,+1,-1)*(radius/std::sqrt(3.0)*a),
	  p+Point<dim>(+1,+1,-1)*(radius/std::sqrt(3.0)*a),
	  p+Point<dim>(+1,+1,+1)*(radius/std::sqrt(3.0)*a),
	  p+Point<dim>(-1,+1,+1)*(radius/std::sqrt(3.0)*a),
					   // now the eight vertices at
					   // the outer sphere
	  p+Point<dim>(-1,-1,-1)*(radius/std::sqrt(3.0)),
	  p+Point<dim>(+1,-1,-1)*(radius/std::sqrt(3.0)),
	  p+Point<dim>(+1,-1,+1)*(radius/std::sqrt(3.0)),
	  p+Point<dim>(-1,-1,+1)*(radius/std::sqrt(3.0)),
	  p+Point<dim>(-1,+1,-1)*(radius/std::sqrt(3.0)),
	  p+Point<dim>(+1,+1,-1)*(radius/std::sqrt(3.0)),
	  p+Point<dim>(+1,+1,+1)*(radius/std::sqrt(3.0)),
	  p+Point<dim>(-1,+1,+1)*(radius/std::sqrt(3.0)),
    };

				   // one needs to draw the seven cubes to
				   // understand what's going on here
  const unsigned int n_cells = 7;
  const int cell_vertices[n_cells][8] = {{0, 1, 4, 5, 3, 2, 7, 6}, // center
					 {8, 9, 12, 13, 0, 1, 4, 5}, // bottom
					 {9, 13, 1, 5, 10, 14, 2, 6}, // right
					 {11, 10, 3, 2, 15, 14, 7, 6}, // top
					 {8, 0, 12, 4, 11, 3, 15, 7}, // left
					 {8, 9, 0, 1, 11, 10, 3, 2}, // front
					 {12, 4, 13, 5, 15, 7, 14, 6}}; // back
  
  std::vector<CellData<dim> > cells (n_cells, CellData<dim>());
  
  for (unsigned int i=0; i<n_cells; ++i) 
    {
      for (unsigned int j=0; j<GeometryInfo<dim>::vertices_per_cell; ++j)
	cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    };

  tria.create_triangulation (
    std::vector<Point<dim> >(&vertices[0], &vertices[n_vertices]),
    cells,
    SubCellData());       // no boundary information
}



// Implementation for 3D only
template <int dim>
void
GridGenerator::cylinder (Triangulation<dim> &tria,
			 const double radius,
			 const double half_length)
{
  Assert (dim <= 3, ExcNotImplemented());
  
				   // Copy the base from hyper_ball<dim>
				   // and transform it to yz
  const double d = radius/std::sqrt(2.0);
  const double a = d/(1+std::sqrt(2.0));
  Point<dim> vertices[24] = {
	Point<dim>(-d, -half_length,-d),
	Point<dim>( d, -half_length,-d),
	Point<dim>(-a, -half_length,-a),
	Point<dim>( a, -half_length,-a),
	Point<dim>(-a, -half_length, a),
	Point<dim>( a, -half_length, a),
	Point<dim>(-d, -half_length, d),
	Point<dim>( d, -half_length, d),
	Point<dim>(-d, 0,-d),
	Point<dim>( d, 0,-d),
	Point<dim>(-a, 0,-a),
	Point<dim>( a, 0,-a),
	Point<dim>(-a, 0, a),
	Point<dim>( a, 0, a),
	Point<dim>(-d, 0, d),
	Point<dim>( d, 0, d),
	Point<dim>(-d, half_length,-d),
	Point<dim>( d, half_length,-d),
	Point<dim>(-a, half_length,-a),
	Point<dim>( a, half_length,-a),
	Point<dim>(-a, half_length, a),
	Point<dim>( a, half_length, a),
	Point<dim>(-d, half_length, d),
	Point<dim>( d, half_length, d),
  };
				   // Turn cylinder such that y->x
  for (unsigned int i=0;i<24;++i)
    {
      const double h = vertices[i](1);
      vertices[i](1) = -vertices[i](0);
      vertices[i](0) = h;
    }
  
  int cell_vertices[10][8] = {
	{0, 1, 8, 9, 2, 3, 10, 11},
	{0, 2, 8, 10, 6, 4, 14, 12},
	{2, 3, 10, 11, 4, 5, 12, 13},
	{1, 7, 9, 15, 3, 5, 11, 13},
	{6, 4, 14, 12, 7, 5, 15, 13}
  };
  for (unsigned int i=0;i<5;++i)
    for (unsigned int j=0;j<8;++j)
      cell_vertices[i+5][j] = cell_vertices[i][j]+8;
  
  std::vector<CellData<dim> > cells (10, CellData<dim>());
  
  for (unsigned int i=0; i<10; ++i) 
    {
      for (unsigned int j=0; j<8; ++j)
	cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    };

  tria.create_triangulation (
    std::vector<Point<dim> >(&vertices[0], &vertices[24]),
    cells,
    SubCellData());       // no boundary information

				   // set boundary indicators for the
				   // faces at the ends to 1 and 2,
				   // respectively. note that we also
				   // have to deal with those lines
				   // that are purely in the interior
				   // of the ends. we determine wether
				   // an edge is purely in the
				   // interior if one of its vertices
				   // is at coordinates '+-a' as set
				   // above
  typename Triangulation<dim>::cell_iterator cell = tria.begin();
  typename Triangulation<dim>::cell_iterator end = tria.end();
  
  for (; cell != end; ++cell)
    for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
      if (cell->at_boundary(i))
	{
	  if (cell->face(i)->center()(0) > half_length-1.e-5)
	    {
	      cell->face(i)->set_boundary_indicator(2);

	      for (unsigned int e=0; e<GeometryInfo<dim>::lines_per_face; ++e)
		if ((std::fabs(cell->face(i)->line(e)->vertex(0)[1]) == a) ||
		    (std::fabs(cell->face(i)->line(e)->vertex(0)[2]) == a) ||
		    (std::fabs(cell->face(i)->line(e)->vertex(1)[1]) == a) ||
		    (std::fabs(cell->face(i)->line(e)->vertex(1)[2]) == a))
		  cell->face(i)->line(e)->set_boundary_indicator(2);
	    }
	  else if (cell->face(i)->center()(0) < -half_length+1.e-5)
	    {
	      cell->face(i)->set_boundary_indicator(1);

	      for (unsigned int e=0; e<GeometryInfo<dim>::lines_per_face; ++e)
		if ((std::fabs(cell->face(i)->line(e)->vertex(0)[1]) == a) ||
		    (std::fabs(cell->face(i)->line(e)->vertex(0)[2]) == a) ||
		    (std::fabs(cell->face(i)->line(e)->vertex(1)[1]) == a) ||
		    (std::fabs(cell->face(i)->line(e)->vertex(1)[2]) == a))
		  cell->face(i)->line(e)->set_boundary_indicator(1);
	    }
	}
}



// Implementation for 3D only
template <int dim>
void
GridGenerator::half_hyper_ball (Triangulation<dim>& tria,
				const Point<dim>& center,
				const double radius)
{
    // These are for the two lower squares
    const double d = radius/std::sqrt(2.0);
    const double a = d/(1+std::sqrt(2.0));
    // These are for the two upper square
    const double b = a/2.0;
    const double c = d/2.0;
    // And so are these
    const double hb = radius*std::sqrt(3.0)/4.0;
    const double hc = radius*std::sqrt(3.0)/2.0;

    Point<dim> vertices[16] = {
	center+Point<dim>( 0,  d, -d),
	center+Point<dim>( 0, -d, -d),
	center+Point<dim>( 0,  a, -a),
	center+Point<dim>( 0, -a, -a),
	center+Point<dim>( 0,  a,  a),
	center+Point<dim>( 0, -a,  a),
	center+Point<dim>( 0,  d,  d),
	center+Point<dim>( 0, -d,  d),

	center+Point<dim>(hc,  c, -c),
	center+Point<dim>(hc, -c, -c),
	center+Point<dim>(hb,  b, -b),
	center+Point<dim>(hb, -b, -b),
	center+Point<dim>(hb,  b,  b),
	center+Point<dim>(hb, -b,  b),
	center+Point<dim>(hc,  c,  c),
	center+Point<dim>(hc, -c,  c),
    };

    int cell_vertices[6][8] = {
	{0, 1, 8, 9, 2, 3, 10, 11},
	{0, 2, 8, 10, 6, 4, 14, 12},
	{2, 3, 10, 11, 4, 5, 12, 13},
	{1, 7, 9, 15, 3, 5, 11, 13},
	{6, 4, 14, 12, 7, 5, 15, 13},
	{8, 10, 9, 11, 14, 12, 15, 13}
    };

    std::vector<CellData<dim> > cells (6, CellData<dim>());

    for (unsigned int i=0; i<6; ++i) 
    {
	for (unsigned int j=0; j<8; ++j)
	    cells[i].vertices[j] = cell_vertices[i][j];
	cells[i].material_id = 0;
    };

    tria.create_triangulation (
	    std::vector<Point<dim> >(&vertices[0], &vertices[16]),
	    cells,
	    SubCellData());       // no boundary information

    typename Triangulation<dim>::cell_iterator cell = tria.begin();
    typename Triangulation<dim>::cell_iterator end = tria.end();

    while (cell != end)
    {
	for (unsigned int i=0;i<GeometryInfo<dim>::faces_per_cell;++i)
	{
	    if (!cell->at_boundary(i))
		continue;

	    // If the center is on the plane x=0, this is a planar
	    // element
	    if (cell->face(i)->center()(0) < center(0)+1.e-5) {
		cell->face(i)->set_boundary_indicator(1);
		for (unsigned int j=0;j<GeometryInfo<dim>::lines_per_face;++j) 
		    cell->face(i)->line(j)->set_boundary_indicator(1);
	    }
	}
	// With this loop we restore back the indicator of the outer lines
	for (unsigned int i=0;i<GeometryInfo<dim>::faces_per_cell;++i)
	{
	    if (!cell->at_boundary(i))
		continue;

	    // If the center is not on the plane x=0, this is a curvilinear
	    // element
	    if (cell->face(i)->center()(0) > center(0)+1.e-5) {
		for (unsigned int j=0;j<GeometryInfo<dim>::lines_per_face;++j) 
		    cell->face(i)->line(j)->set_boundary_indicator(0);
	    }
	}
	++cell;
    }
}

// Implementation for 3D only
template <int dim>
void GridGenerator::hyper_shell (Triangulation<dim>& tria,
				 const Point<dim>& p,
				 const double inner_radius,
				 const double outer_radius,
				 const unsigned int n,
				 const bool colorize)
{
  Assert ((inner_radius > 0) && (inner_radius < outer_radius),
	  ExcInvalidRadii ());

  const double irad = inner_radius/std::sqrt(3.0);
  const double orad = outer_radius/std::sqrt(3.0);
  std::vector<Point<dim> > vertices;
  std::vector<CellData<dim> > cells;
  
				   // Start with the shell bounded by
				   // two nested cubes
  if (n <= 6)
    {
      for (unsigned int i=0;i<8;++i)
	vertices.push_back(p+hexahedron[i]*irad);
      for (unsigned int i=0;i<8;++i)
	vertices.push_back(p+hexahedron[i]*orad);
      
      const unsigned int n_cells = 6;
      const int cell_vertices[n_cells][8] =
	{{8, 9, 10, 11, 0, 1, 2, 3}, // bottom
	 {9, 11, 1, 3, 13, 15, 5, 7}, // right
	 {12, 13, 4, 5, 14, 15, 6, 7}, // top
	 {8, 0, 10, 2, 12, 4, 14, 6}, // left
	 {8, 9, 0, 1, 12, 13, 4, 5}, // front
	 {10, 2, 11, 3, 14, 6, 15, 7}}; // back
      
      cells.resize(n_cells, CellData<dim>());
      
      for (unsigned int i=0; i<n_cells; ++i) 
	{
	  for (unsigned int j=0; j<GeometryInfo<dim>::vertices_per_cell; ++j)
	    cells[i].vertices[j] = cell_vertices[i][j];
	  cells[i].material_id = 0;
	}
    }
				   // A more regular subdivision can
				   // be obtained by two nested
				   // rhombic dodecahedra
  else if (n <= 12)
    {
      for (unsigned int i=0;i<8;++i)
	vertices.push_back(p+hexahedron[i]*irad);
      for (unsigned int i=0;i<6;++i)
	vertices.push_back(p+octahedron[i]*inner_radius);
      for (unsigned int i=0;i<8;++i)
	vertices.push_back(p+hexahedron[i]*orad);
      for (unsigned int i=0;i<6;++i)
	vertices.push_back(p+octahedron[i]*outer_radius);

      const unsigned int n_cells = 12;
      const unsigned int rhombi[n_cells][4] =
	{{ 10,  4,  0,  8},
	 {  4, 13,  8,  6},
	 { 10,  5,  4, 13},
	 {  1,  9, 10,  5},
	 {  9,  7,  5, 13},
	 {  7, 11, 13,  6},
	 {  9,  3,  7, 11},
	 {  1, 12,  9,  3},
	 { 12,  2,  3, 11},
	 {  2,  8, 11,  6},
	 { 12,  0,  2,  8},
	 {  1, 10, 12,  0}};
      
      cells.resize(n_cells, CellData<dim>());
      
      for (unsigned int i=0; i<n_cells; ++i) 
	{
	  for (unsigned int j=0; j<4; ++j)
	    {
	      cells[i].vertices[j  ] = rhombi[i][j];
	      cells[i].vertices[j+4] = rhombi[i][j] + 14;
	    }
	  cells[i].material_id = 0;
	}
    }
  else
    {
      Assert(false, ExcIndexRange(n, 1, 7));
    }
  
  tria.create_triangulation (vertices, cells,
			     SubCellData());       // no boundary
						   // information
  
  if (colorize)
    colorize_hyper_shell(tria, p, inner_radius, outer_radius);
}


template<int dim>
void
GridGenerator::colorize_hyper_shell (
  Triangulation<dim>& tria,
  const Point<dim>&, const double, const double)
{
				   // Inspite of receiving geometrical
				   // data, we do this only based on
				   // topology.

				   // For the mesh based on  cube,
				   // this is highly irregular
  if (tria.n_cells() == 6)
    {
      typename Triangulation<dim>::cell_iterator cell = tria.begin();
      cell->face(4)->set_boundary_indicator(1);
      (++cell)->face(2)->set_boundary_indicator(1);
      (++cell)->face(2)->set_boundary_indicator(1);
      (++cell)->face(0)->set_boundary_indicator(1);
      (++cell)->face(2)->set_boundary_indicator(1);
      (++cell)->face(0)->set_boundary_indicator(1);      
    }
  else
				     // For higher polyhedra, this is regular.
    {
      for (typename Triangulation<dim>::cell_iterator cell = tria.begin();
	   cell != tria.end(); ++cell)
	cell->face(5)->set_boundary_indicator(1);
    }
}


// Implementation for 3D only
template <int dim>
void
GridGenerator::half_hyper_shell (Triangulation<dim>&,
				 const Point<dim>&,
				 const double,
				 const double,
				 const unsigned int)
{
  Assert (false, ExcNotImplemented());
}



// Implementation for 3D only
template <int dim>
void GridGenerator::cylinder_shell (Triangulation<dim>   &tria,
                                    const double        length,
                                    const double        inner_radius,
                                    const double        outer_radius,
                                    const unsigned int  n_radial_cells,
                                    const unsigned int  n_axial_cells)
{
  Assert ((inner_radius > 0) && (inner_radius < outer_radius),
	  ExcInvalidRadii ());
  
  const double pi = numbers::PI;
  
				   // determine the number of cells
				   // for the grid. if not provided by
				   // the user determine it such that
				   // the length of each cell on the
				   // median (in the middle between
				   // the two circles) is equal to its
				   // radial extent (which is the
				   // difference between the two
				   // radii)
  const unsigned int N_r = (n_radial_cells == 0 ?
                            static_cast<unsigned int>
                            (std::ceil((2*pi* (outer_radius + inner_radius)/2) /
                                       (outer_radius - inner_radius))) :
                            n_radial_cells);
  const unsigned int N_z = (n_axial_cells == 0 ?
                            static_cast<unsigned int>
                            (std::ceil (length /
                                        (2*pi*(outer_radius + inner_radius)/2/N_r))) :
                            n_axial_cells);

				   // set up N vertices on the
				   // outer and N vertices on
				   // the inner circle. the
				   // first N ones are on the
				   // outer one, and all are
				   // numbered counter-clockwise
  std::vector<Point<2> > vertices_2d(2*N_r);
  for (unsigned int i=0; i<N_r; ++i)
    {
      vertices_2d[i] = Point<2>( std::cos(2*pi * i/N_r),
                                 std::sin(2*pi * i/N_r)) * outer_radius;
      vertices_2d[i+N_r] = vertices_2d[i] * (inner_radius/outer_radius);
    };

  std::vector<Point<3> > vertices_3d;
  vertices_3d.reserve (2*N_r*(N_z+1));
  for (unsigned int j=0; j<=N_z; ++j)
    for (unsigned int i=0; i<2*N_r; ++i)
      {
        const Point<3> v (vertices_2d[i][0],
                          vertices_2d[i][1],
                          j*length/N_z);
        vertices_3d.push_back (v);
      }
                            
  std::vector<CellData<dim> > cells (N_r*N_z, CellData<dim>());
  
  for (unsigned int j=0; j<N_z; ++j)
    for (unsigned int i=0; i<N_r; ++i) 
      {
        cells[i+j*N_r].vertices[0] = i + (j+1)*2*N_r;
        cells[i+j*N_r].vertices[1] = (i+1)%N_r + (j+1)*2*N_r;
        cells[i+j*N_r].vertices[2] = i + j*2*N_r;
        cells[i+j*N_r].vertices[3] = (i+1)%N_r + j*2*N_r;

        cells[i+j*N_r].vertices[4] = N_r+i + (j+1)*2*N_r;
        cells[i+j*N_r].vertices[5] = N_r+((i+1)%N_r) + (j+1)*2*N_r;
        cells[i+j*N_r].vertices[6] = N_r+i + j*2*N_r;
        cells[i+j*N_r].vertices[7] = N_r+((i+1)%N_r) + j*2*N_r;
        
        cells[i+j*N_r].material_id = 0;
      }
  
  tria.create_triangulation (
    vertices_3d, cells, SubCellData());
}



#endif


// make the following function inline. this is so that it becomes marked
// internal/weak for the linker and we don't get multiply defined errors
// when linking with more than one dimension at a time. Usually we used
// the trick of putting these functions in a .all_dimensions.cc file, but
// this is not necessary here as this is an internal only function.
inline
void GridGenerator::laplace_solve (const SparseMatrix<double> &S,
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


#if deal_II_dimension == 1

// Implementation for 1D only
template <int dim>
void GridGenerator::laplace_transformation (Triangulation<dim> &,
					    const std::map<unsigned int,Point<dim> > &)
{
  Assert(false, ExcNotImplemented());
}

#else

// Implementation for dimensions except 1
template <int dim>
void GridGenerator::laplace_transformation (Triangulation<dim> &tria,
					    const std::map<unsigned int,Point<dim> > &new_points)
{
				   // first provide everything that is
				   // needed for solving a Laplace
				   // equation.  
  MappingQ1<dim> mapping_q1;  
  FE_Q<dim> q1(1);
  
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(q1);
  SparsityPattern sparsity_pattern (dof_handler.n_dofs (), dof_handler.n_dofs (),
				    dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress ();
  
  SparseMatrix<double> S(sparsity_pattern);
  
  QGauss4<dim> quadrature;
  
  MatrixCreator::create_laplace_matrix(mapping_q1, dof_handler, quadrature, S);

				   // set up the boundary values for
				   // the laplace problem
  std::vector<std::map<unsigned int,double> > m(dim);
  typename std::map<unsigned int,Point<dim> >::const_iterator map_iter;
  typename std::map<unsigned int,Point<dim> >::const_iterator map_end=new_points.end();

				   // fill these maps using the data
				   // given by new_points
  typename DoFHandler<dim>::cell_iterator cell=dof_handler.begin_active(),
					  endc=dof_handler.end();
  typename DoFHandler<dim>::face_iterator face;
  for (; cell!=endc; ++cell)
    {
      if (cell->at_boundary())
	for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
	  {
	    face=cell->face(face_no);
	    if (face->at_boundary())
	      for (unsigned int vertex_no=0;
		   vertex_no<GeometryInfo<dim>::vertices_per_face; ++vertex_no)
		{
		  const unsigned int vertex_index=face->vertex_index(vertex_no);
		  map_iter=new_points.find(vertex_index);

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
  Threads::ThreadGroup<> threads;
  for (unsigned int i=0; i<dim; ++i)
    threads += Threads::spawn (&GridGenerator::laplace_solve)(S, m[i], us[i]);
  threads.join_all ();
  
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


#endif


#if deal_II_dimension == 1

template<int dim>
void GridGenerator::hyper_cube_with_cylindrical_hole (Triangulation<dim> &,
				   const double,
				   const double,
				   const double,
				   const unsigned int,
				   bool){
  Assert(false, ExcNotImplemented());
}

#endif


#if deal_II_dimension == 2

template<int dim>
void GridGenerator::hyper_cube_with_cylindrical_hole (Triangulation<dim> &triangulation, 
				   const double inner_radius,
				   const double outer_radius,
				   const double, // width,
				   const unsigned int, // width_repetition,
				   bool colorize)
{
  Assert(inner_radius < outer_radius,
	 ExcMessage("outer_radius has to be bigger than inner_radius."));

  Point<dim> center;
  // We create an hyper_shell in two dimensions, and then we modify it.
  GridGenerator::hyper_shell (triangulation,
			      center, inner_radius, outer_radius,
                              8);
  typename Triangulation<dim>::active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();
  std::vector<bool> treated_vertices(triangulation.n_vertices(), false);
  for(; cell != endc; ++cell) {
    for(unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      if(cell->face(f)->at_boundary()) {
	for(unsigned int v=0; v < GeometryInfo<dim>::vertices_per_face; ++v) {
	  unsigned int vv = cell->face(f)->vertex_index(v);
	  if(treated_vertices[vv] == false) {
	    treated_vertices[vv] = true;
	    switch(vv) {
	    case 1:
	      cell->face(f)->vertex(v) = center+Point<dim>(outer_radius,outer_radius);
	      break;
	    case 3:
	      cell->face(f)->vertex(v) = center+Point<dim>(-outer_radius,outer_radius);
	      break;
	    case 5:
	      cell->face(f)->vertex(v) = center+Point<dim>(-outer_radius,-outer_radius);
	      break;
	    case 7:
	      cell->face(f)->vertex(v) = center+Point<dim>(outer_radius,-outer_radius);
	    default:
	      break;
	    }
	  }
	}
      } 
  }
  double eps = 1e-3 * outer_radius;
  cell = triangulation.begin_active();
   for(; cell != endc; ++cell) {
    for(unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      if(cell->face(f)->at_boundary()) {
	double dx = cell->face(f)->center()(0) - center(0);
	double dy = cell->face(f)->center()(1) - center(1);
	if(colorize) {
	  if(std::abs(dx + outer_radius) < eps)
	    cell->face(f)->set_boundary_indicator(0);
	  else if(std::abs(dx - outer_radius) < eps)
	    cell->face(f)->set_boundary_indicator(1);
	  else if(std::abs(dy + outer_radius) < eps)
	    cell->face(f)->set_boundary_indicator(2);
	  else if(std::abs(dy - outer_radius) < eps)
	    cell->face(f)->set_boundary_indicator(3);
	  else 
	    cell->face(f)->set_boundary_indicator(4);
	} else {
	  double d = (cell->face(f)->center() - center).norm();
	  if(d-inner_radius < 0)
	    cell->face(f)->set_boundary_indicator(1);
	  else
	    cell->face(f)->set_boundary_indicator(0);
	}
      }
   }	
}

#endif

#if deal_II_dimension == 3

template<int dim>
void GridGenerator::hyper_cube_with_cylindrical_hole(Triangulation<dim> &triangulation, 
						 const double inner_radius,
						 const double outer_radius,
						 const double L,
						 const unsigned int Nz,
						 bool colorize)
{
  Assert(inner_radius < outer_radius,
	 ExcMessage("outer_radius has to be bigger than inner_radius."));
  Assert(L > 0, 
	 ExcMessage("Must give positive extension L"));
  Assert(Nz >= 1, ExcLowerRange(1, Nz));
  
  GridGenerator::cylinder_shell (triangulation,
				 L, inner_radius, outer_radius,
				 8,
				 Nz);
  
  typename Triangulation<dim>::active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();
  std::vector<bool> treated_vertices(triangulation.n_vertices(), false);
  for(; cell != endc; ++cell) {
    for(unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      if(cell->face(f)->at_boundary()) {
	for(unsigned int v=0; v < GeometryInfo<dim>::vertices_per_face; ++v) {
	  unsigned int vv = cell->face(f)->vertex_index(v);
	  if(treated_vertices[vv] == false) {
	    treated_vertices[vv] = true;
	    for(unsigned int i=0; i<=Nz; ++i) {
	      double d = ((double) i)*L/((double) Nz);
	      switch(vv-i*16) {
	      case 1:
		cell->face(f)->vertex(v) = Point<dim>(outer_radius,outer_radius,d);
		break;
	      case 3:
		cell->face(f)->vertex(v) = Point<dim>(-outer_radius,outer_radius,d);
		break;
	      case 5:
		cell->face(f)->vertex(v) = Point<dim>(-outer_radius,-outer_radius,d);
		break;
	      case 7:
		cell->face(f)->vertex(v) = Point<dim>(outer_radius,-outer_radius,d);
		break;
	      default:
		break;
	      }
	    }
	  }
	}
      } 
  }
  double eps = 1e-3 * outer_radius;
  cell = triangulation.begin_active();
   for(; cell != endc; ++cell) {
    for(unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
      if(cell->face(f)->at_boundary()) {
	double dx = cell->face(f)->center()(0);
	double dy = cell->face(f)->center()(1);
	double dz = cell->face(f)->center()(2);
	
	if(colorize) {
	  if(std::abs(dx + outer_radius) < eps)
	    cell->face(f)->set_boundary_indicator(0);
	  
	  else if(std::abs(dx - outer_radius) < eps)
	    cell->face(f)->set_boundary_indicator(1);
	  
	  else if(std::abs(dy + outer_radius) < eps)
	    cell->face(f)->set_boundary_indicator(2);
	  
	  else if(std::abs(dy - outer_radius) < eps)
	    cell->face(f)->set_boundary_indicator(3);
	  
	  else if(std::abs(dz) < eps)
	    cell->face(f)->set_boundary_indicator(4);

	  else if(std::abs(dz - L) < eps)
	    cell->face(f)->set_boundary_indicator(5);
	  
	  else {
	    cell->face(f)->set_boundary_indicator(6);
	    for(unsigned int l=0; l<GeometryInfo<dim>::lines_per_face; ++l)
	      cell->face(f)->line(l)->set_boundary_indicator(6);
	  }

	} else {
	  Point<dim> c = cell->face(f)->center();
	  c(2) = 0;
	  double d = c.norm();
	  if(d-inner_radius < 0) {
	    cell->face(f)->set_boundary_indicator(1);
	    for(unsigned int l=0; l<GeometryInfo<dim>::lines_per_face; ++l)
	      cell->face(f)->line(l)->set_boundary_indicator(1);
	  } else
	    cell->face(f)->set_boundary_indicator(0);
	}
      }
   }	
}

#endif


// explicit instantiations
template void
GridGenerator::hyper_cube<deal_II_dimension> (
  Triangulation<deal_II_dimension> &, const double, const double);

template void
GridGenerator::subdivided_hyper_cube<deal_II_dimension> (
  Triangulation<deal_II_dimension> &,
  const unsigned int, const double, const double);

template void
GridGenerator::hyper_rectangle<deal_II_dimension> (
  Triangulation<deal_II_dimension> &,
  const Point<deal_II_dimension>&, const Point<deal_II_dimension>&,
  const bool);

template void
GridGenerator::subdivided_hyper_rectangle<deal_II_dimension> 
(Triangulation<deal_II_dimension> &,
 const std::vector<unsigned int>&,
 const Point<deal_II_dimension>&,
 const Point<deal_II_dimension>&, bool);

template void
GridGenerator::subdivided_hyper_rectangle<deal_II_dimension> 
(Triangulation<deal_II_dimension> &,
 const std::vector<std::vector<double> >&,
 const Point<deal_II_dimension>&,
 const Point<deal_II_dimension>&, bool);

template void
GridGenerator::parallelogram<deal_II_dimension> (
  Triangulation<deal_II_dimension> &,
  const Tensor<2,deal_II_dimension>&,
  const bool);

template void
GridGenerator::enclosed_hyper_cube (
  Triangulation<deal_II_dimension>&, double, double, double, bool);

template void
GridGenerator::hyper_ball (
  Triangulation<deal_II_dimension>&,
  const Point<deal_II_dimension>&, double);

template void
GridGenerator::cylinder (
  Triangulation<deal_II_dimension>&, double, double);


template void
GridGenerator::hyper_L (
  Triangulation<deal_II_dimension>&, double, double);

template void
GridGenerator::hyper_cube_slit (
  Triangulation<deal_II_dimension>&, double, double, bool);

template void
GridGenerator::hyper_shell (
  Triangulation<deal_II_dimension>&,
  const Point<deal_II_dimension>&, double, double, unsigned int, bool);


template void
GridGenerator::cylinder_shell (
  Triangulation<deal_II_dimension>&,
  double, double, double, unsigned int, unsigned int);

template void
GridGenerator::half_hyper_ball (
  Triangulation<deal_II_dimension>&, const Point<deal_II_dimension>&, double);

template void
GridGenerator::half_hyper_shell (
  Triangulation<deal_II_dimension>&,
  const Point<deal_II_dimension>&, double, double, unsigned int);


template void 
GridGenerator::hyper_cube_with_cylindrical_hole (
  Triangulation<deal_II_dimension> &,
  const double, const double, const double, const unsigned int, bool);

template void
GridGenerator::colorize_hyper_shell(
  Triangulation<deal_II_dimension>& tria,
  const Point<deal_II_dimension>& center,
  const double inner_radius, const double outer_radius);

template void
GridGenerator::
laplace_transformation<deal_II_dimension> (Triangulation<deal_II_dimension> &,
					   const std::map<unsigned int,Point<deal_II_dimension> > &);

#if deal_II_dimension != 3

template void
GridGenerator::hyper_cube<deal_II_dimension, deal_II_dimension+1> (
  Triangulation<deal_II_dimension,deal_II_dimension+1> &, const double, const double);
template void
GridGenerator::hyper_rectangle<deal_II_dimension,deal_II_dimension+1> (
  Triangulation<deal_II_dimension,deal_II_dimension+1> &,
  const Point<deal_II_dimension+1>&, const Point<deal_II_dimension+1>&,
  const bool);

#endif


DEAL_II_NAMESPACE_CLOSE
