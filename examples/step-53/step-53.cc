/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2014 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Authors: Wolfgang Bangerth, Texas A&M University, 2014
 *          Luca Heltai, SISSA, 2014
 *          D. Sarah Stamps, MIT, 2014
 */

// Let us start with the include files we need here. Obviously, we need the
// ones that describe the triangulation (<code>tria.h</code>), and that allow
// us to create and output triangulations (<code>grid_generator.h</code> and
// <code>grid_out.h</code>). Furthermore, we need the header file that
// declares the Manifold and ChartManifold classes that we will need to
// describe the geometry (<code>manifold.h</code>). We will then also need
// the GridTools::transform() function from the last of the following header
// files; the purpose for this function will become discussed at the point
// where we use it.
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold.h>
#include <deal.II/grid/grid_tools.h>

// The remainder of the include files relate to reading the topography data.
// As explained in the introduction, we will read it from a file and then
// use the Functions::InterpolatedUniformGridData class that is declared in the
// first of the following header files. Because the data is large, the file
// we read from is stored as gzip compressed data and we make use of
// some BOOST-provided functionality to read directly from gzipped data.
#include <deal.II/base/function_lib.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file.hpp>

#include <fstream>


// The final part of the top matter is to open a namespace into which to put
// everything, and then to import the dealii namespace into it.
namespace Step53
{
  using namespace dealii;


  // @sect3{Describing topography: AfricaTopography}
  //
  // The first significant part of this program is the class that describes
  // the topography $h(\hat phi,\hat \theta)$ as a function of longitude
  // and latitude. As discussed in the introduction, we will make our life
  // a bit easier here by not writing the class in the most general way
  // possible but by only writing it for the particular purpose we are
  // interested in here: interpolating data obtained from one very specific
  // data file that contains information about a particular area of the
  // world for which we know the extents.
  //
  // The general layout of the class has been discussed already above.
  // Following is its declaration, including three static member functions
  // that we will need in initializing the <code>topography_data</code>
  // member variable.
  class AfricaTopography
  {
  public:
    AfricaTopography ();

    double value (const double lon,
                  const double lat) const;

  private:
    const Functions::InterpolatedUniformGridData<2> topography_data;

    static std_cxx11::array<std::pair<double,double>,2> get_endpoints ();
    static std_cxx11::array<unsigned int,2>             n_intervals ();
    static std::vector<double>                           get_data ();
  };


  // Let us move to the implementation of the class. The interesting parts
  // of the class are the constructor and the <code>value()</code> function.
  // The former initializes the Functions::InterpolatedUniformGridData member
  // variable and we will use the constructor that requires us to pass in
  // the end points of the 2-dimensional data set we want to interpolate
  // (which are here given by the intervals $[-6.983333, 11.98333]$,
  // using the trick of switching end points discussed in the introduction,
  // and $[25, 35.983333]$, both given in degrees), the number of intervals
  // into which the data is split (379 in latitude direction and 219 in
  // longitude direction, for a total of $380\times 220$ data points), and
  // a Table object that contains the data. The data then of course has
  // size $380\times 220$ and we initialize it by providing an iterator
  // to the first of the 83,600 elements of a std::vector object returned
  // by the <code>get_data()</code> function below. Note that all of the
  // member functions we call here are static because (i) they do not
  // access any member variables of the class, and (ii) because they are
  // called at a time when the object is not initialized fully anyway.
  AfricaTopography::AfricaTopography ()
    :
    topography_data (get_endpoints(),
                     n_intervals(),
                     Table<2,double> (380, 220,
                                      get_data().begin()))
  {}


  double
  AfricaTopography::value (const double lon,
                           const double lat) const
  {
    return topography_data.value (Point<2>(-lat * 180/numbers::PI,
                                           lon * 180/numbers::PI));
  }


  std_cxx11::array<std::pair<double,double>,2>
  AfricaTopography::get_endpoints ()
  {
    std_cxx11::array<std::pair<double,double>,2> endpoints;
    endpoints[0] = std::make_pair (-6.983333, 11.966667);
    endpoints[1] = std::make_pair (25, 35.95);
    return endpoints;
  }


  std_cxx11::array<unsigned int,2>
  AfricaTopography::n_intervals ()
  {
    std_cxx11::array<unsigned int,2> endpoints;
    endpoints[0] = 379;
    endpoints[1] = 219;
    return endpoints;
  }


  // The only other function of greater interest is the <code>get_data()</code>
  // function. It returns a temporary vector that contains all 83,600 data
  // points describing the altitude and is read from the file
  // <code>topography.txt.gz</code>. Because the file is compressed by gzip,
  // we cannot just read it through an object of type std::ifstream, but
  // there are convenient methods in the BOOST library (see
  // http://www.boost.org) that allows us to read from compressed files
  // without first having to uncompress it on disk. The result is, basically,
  // just another input stream that, for all practical purposes, looks just like
  // the ones we always use.
  //
  // When reading the data, we read the three columns but throw ignore the
  // first two. The datum in the last column is appended to an array that we
  // the return and that will be copied into the table from which
  // <code>topography_data</code> is initialized. Since the BOOST.iostreams
  // library does not provide a very useful exception when the input file
  // does not exist, is not readable, or does not contain the correct
  // number of data lines, we catch all exceptions it may produce and
  // create our own one. To this end, in the <code>catch</code>
  // clause, we let the program run into an <code>AssertThrow(false, ...)</code>
  // statement. Since the condition is always false, this always triggers an
  // exception. In other words, this is equivalent to writing
  // <code>throw ExcMessage("...")</code> but it also fills certain fields
  // in the exception object that will later be printed on the screen
  // identifying the function, file and line where the exception happened.
  std::vector<double>
  AfricaTopography::get_data ()
  {
    std::vector<double> data;

    // create a stream where we read from gzipped data
    boost::iostreams::filtering_istream in;
    in.push(boost::iostreams::basic_gzip_decompressor<>());
    in.push(boost::iostreams::file_source("topography.txt.gz"));

    for (unsigned int line=0; line<83600; ++line)
      {
        try
          {
            double lat, lon, elevation;
            in >> lat >> lon >> elevation;

            data.push_back (elevation);
          }
        catch (...)
          {
            AssertThrow (false,
                         ExcMessage ("Could not read all 83,600 data points "
                                     "from the file <topography.txt.gz>!"));
          }
      }

    return data;
  }


  // @sect3{Describing the geometry: AfricaGeometry}
  //
  // The following class is then the main one of this program. Its structure
  // has been described in much detail in the introduction and does not need
  // much introduction any more.
  class AfricaGeometry : public ChartManifold<3,3>
  {
  public:
    virtual
    Point<3>
    pull_back(const Point<3> &space_point) const;

    virtual
    Point<3>
    push_forward(const Point<3> &chart_point) const;

  private:
    static const double    R;
    static const double    ellipticity;

    const AfricaTopography topography;

    Point<3> push_forward_wgs84 (const Point<3> &phi_theta_d) const;
    Point<3> pull_back_wgs84 (const Point<3> &x) const;

    Point<3> push_forward_topo (const Point<3> &phi_theta_d_hat) const;
    Point<3> pull_back_topo (const Point<3> &phi_theta_d) const;
  };


  const double AfricaGeometry::R           = 6378137;
  const double AfricaGeometry::ellipticity = 8.1819190842622e-2;


  // The implementation, as well, is pretty straightforward if you have
  // read the introduction. In particular, both of the pull back and
  // push forward functions are just concatenations of the respective
  // functions of the WGS 84 and topography mappings:
  Point<3>
  AfricaGeometry::pull_back(const Point<3> &space_point) const
  {
    return pull_back_topo (pull_back_wgs84 (space_point));
  }

  Point<3>
  AfricaGeometry::push_forward(const Point<3> &chart_point) const
  {
    return push_forward_wgs84 (push_forward_topo (chart_point));
  }


  // The following two functions then define the forward and inverse
  // transformations that correspond to the WGS 84 reference shape of
  // Earth. The forward transform follows the formula shown in the
  // introduction. The inverse transform is significantly more complicated
  // and is, at the very least, not intuitive. It also suffers from the
  // fact that it returns an angle that at the end of the function we
  // need to clip back into the interval $[0,2\pi]$ if it should have
  // escaped from there.
  Point<3>
  AfricaGeometry::push_forward_wgs84(const Point<3> &phi_theta_d) const
  {
    const double phi   = phi_theta_d[0];
    const double theta = phi_theta_d[1];
    const double d     = phi_theta_d[2];

    const double R_bar = R / std::sqrt(1 - (ellipticity * ellipticity *
                                            std::sin(theta) * std::sin(theta)));

    return Point<3> ((R_bar + d) * std::cos(phi) * std::cos(theta),
                     (R_bar + d) * std::sin(phi) * std::cos(theta),
                     ((1 - ellipticity * ellipticity) * R_bar + d) * std::sin(theta));
  }

  Point<3>
  AfricaGeometry::pull_back_wgs84(const Point<3> &x) const
  {
    const double b     = std::sqrt(R * R * (1 - ellipticity * ellipticity));
    const double ep    = std::sqrt((R * R - b * b) / (b * b));
    const double p     = std::sqrt(x(0) * x(0) + x(1) * x(1));
    const double th    = std::atan2(R * x(2), b * p);
    const double phi   = std::atan2(x(1), x(0));
    const double theta = std::atan2(x(2) + ep * ep * b * std::pow(std::sin(th),3),
                                    (p - (ellipticity * ellipticity * R  * std::pow(std::cos(th),3))));
    const double R_bar = R / (std::sqrt(1 - ellipticity * ellipticity * std::sin(theta) * std::sin(theta)));
    const double R_plus_d = p / std::cos(theta);

    Point<3> phi_theta_d;
    if (phi < 0)
      phi_theta_d[0] = phi + 2*numbers::PI;
    else if (phi > 2*numbers::PI)
      phi_theta_d[0] = phi - 2*numbers::PI;
    else
      phi_theta_d[0] = phi;
    phi_theta_d[1] = theta;
    phi_theta_d[2] = R_plus_d - R_bar;
    return phi_theta_d;
  }


  // In contrast, the topography transformations follow exactly the
  // description in the introduction. There is not consequently not
  // much to add:
  Point<3>
  AfricaGeometry::push_forward_topo(const Point<3> &phi_theta_d_hat) const
  {
    const double d_hat = phi_theta_d_hat[2];
    const double h     = topography.value(phi_theta_d_hat[0],
                                          phi_theta_d_hat[1]);
    const double d = d_hat + (d_hat + 500000)/500000*h;
    const Point<3> phi_theta_d (phi_theta_d_hat[0],
                                phi_theta_d_hat[1],
                                d);
    return phi_theta_d;
  }

  Point<3>
  AfricaGeometry::pull_back_topo(const Point<3> &phi_theta_d) const
  {
    const double d = phi_theta_d[2];
    const double h = topography.value(phi_theta_d[0],
                                      phi_theta_d[1]);
    const double d_hat = 500000 * (d-h)/(500000+h);
    const Point<3> phi_theta_d_hat (phi_theta_d[0],
                                    phi_theta_d[1],
                                    d_hat);
    return phi_theta_d_hat;
  }


  // @sect3{Creating the mesh}
  //
  // Having so described the properties of the geometry, not it is
  // time to deal with the mesh used to discretize it. To this end,
  // we create objects for the geometry and triangulation, and then
  // proceed to create a $1\times 2\times 1$ rectangular mesh that
  // corresponds to the reference domain
  // $\hat U=[26,35]\times[-10,5]\times[-500000,0]$. We choose
  // this number of subdivisions because it leads to cells that
  // are roughly like cubes instead of stretched in one direction or
  // another.
  //
  // Of course, we are not actually interested in meshing the
  // reference domain. We are interested in meshing the real domain.
  // Consequently, we will use the GridTools::transform() function
  // that simply moves every point of a triangulation according to
  // a given transformation. The transformation function it wants is
  // a function that takes as its single argument a point in the reference
  // domain and returns the corresponding location in the domain that we
  // want to map to. This is, of course, exactly the push forward
  // function of the geometry we use. However,
  // <code>AfricaGeometry::push_forward()</code> requires two arguments:
  // the <code>AfricaGeometry</code> object to work with via its implicit
  // <code>this</code> pointer, and the point. We bind the first of these
  // to the geometry object we have created at the top of the function
  // and leave the second one open, obtaining the desired object to
  // do the transformation.
  void run ()
  {
    AfricaGeometry   geometry;
    Triangulation<3> triangulation;

    {
      const Point<3> corner_points[2] = { Point<3>(26*numbers::PI/180,
                                                   -10*numbers::PI/180,
                                                   -500000),
                                          Point<3>(35*numbers::PI/180,
                                                   5*numbers::PI/180,
                                                   0)
                                        };
      std::vector<unsigned int> subdivisions(3);
      subdivisions[0] = 1;
      subdivisions[1] = 2;
      subdivisions[2] = 1;
      GridGenerator::subdivided_hyper_rectangle (triangulation, subdivisions,
                                                 corner_points[0], corner_points[1],
                                                 true);

      GridTools::transform (std_cxx11::bind(&AfricaGeometry::push_forward,
                                            std_cxx11::cref(geometry),
                                            std_cxx11::_1),
                            triangulation);
    }

    // The next step is to explain to the triangulation to use our geometry
    // object whenever a new point is needed upon refining the mesh. We do
    // this by telling the triangulation to use our geometry for everythin
    // that has manifold indicator zero, and then proceed to mark all cells
    // and their bounding faces and edges with manifold indicator zero. This
    // ensures that the triangulation consults our geometry object everytime
    // a new vertex is needed. Since manifold indicators are inherited from
    // mother to children, this also happens after several recursive
    // refinement steps.
    triangulation.set_manifold(0, geometry);
    for (Triangulation<3>::active_cell_iterator cell=triangulation.begin_active();
         cell!=triangulation.end(); ++cell)
      cell->set_all_manifold_ids(0);

    // The last step is to refine the mesh beyond its initial $1\times 2\times 1$
    // coarse mesh. We could just refine globally a number of times, but since for
    // the purpose of this tutorial program we're really only interested in what
    // is happening close to the surface, we just refine 6 times all of the cells
    // that have a face at a boundary with indicator 5. Looking this up in the
    // documentation of the GridGenerator::subdivided_hyper_rectangle() function
    // we have used above reveals that boundary indicator 5 corresponds to the top
    // surface of the domain (and this is what the last <code>true</code> argument
    // in the call to GridGenerator::subdivided_hyper_rectangle() above meant: to
    // "color" the boundaries by assigning each boundary a unique boundary indicator).
    for (unsigned int i=0; i<6; ++i)
      {
        for (Triangulation<3>::active_cell_iterator cell=triangulation.begin_active();
             cell!=triangulation.end(); ++cell)
          for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f)
            if (cell->face(f)->boundary_indicator() == 5)
              {
                cell->set_refine_flag ();
                break;
              }
        triangulation.execute_coarsening_and_refinement();

        std::cout << "Refinement step " << i+1 << ": "
                  << triangulation.n_active_cells() << " cells, "
                  << GridTools::minimal_cell_diameter (triangulation)/1000
                  << "km minimal cell diameter"
                  << std::endl;
      }

    // Having done this all, we can now output the mesh into a file of its own:
    const std::string filename = "mesh.vtu";
    std::ofstream out (filename.c_str());
    GridOut grid_out;
    grid_out.write_vtu (triangulation, out);
  }
}



// @sect3{The main function}

// Finally, the main function, which follows the same scheme used in all
// tutorial programs starting with step-6. There isn't much to do here, only
// to call the single <code>run()</code> function.
int main ()
{
  try
    {
      Step53::run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
}

