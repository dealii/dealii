// $Id$
// (c) Guido Kanschat

// A little program reading a grid *.inp and writing it to *.eps.
// Some more functionality should be added som time.

#include <grid/tria.h>
#include <grid/grid_in.h>
#include <grid/grid_out.h>

#include <iostream>
#include <fstream>
#include <string>

#include <unistd.h>

using namespace std;

template <int dim>
void convert(const char* infile,
	     const char* outfile,
	     GridOut::OutputFormat oformat)
{
  Triangulation<dim> tr;
  GridIn<dim> gin;
  gin.attach_triangulation(tr);
  gin.read(infile);
  
  GridOut gout;
  GridOutFlags::DX dx_flags(true, true, true, false, true);
  gout.set_flags(dx_flags);
  
  ofstream out(outfile);
  gout.write(tr, out, oformat);
}


int main(int argc, char** argv)
{
  if (argc<4)
    {
      cerr << "Usage: " << argv[0]
	   << " dim infile outfile" << endl;
      exit(1);
    }

  const unsigned int d = atoi(argv[1]);
  
  std::string outstring(argv[3]);
  GridOut::OutputFormat oformat = GridOut::eps;
  
  const unsigned int dotpos = outstring.find_last_of('.');
  if (dotpos < outstring.length())
    {
      std::string ext = outstring.substr(dotpos+1);
      if (ext == "inp")
	ext = "ucd";
      
      oformat = GridOut::parse_output_format(ext);
    }

  if (d == 2)
    convert<2>(argv[2], argv[3], oformat);
  else if (d == 3)
    convert<3>(argv[2], argv[3], oformat);
}
