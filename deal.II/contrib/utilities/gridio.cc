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

int main(int argc, char** argv)
{
  if (argc<2)
    {
      cerr << "Usage: " << argv[0]
	   << " -i informat -o outformat <filename>" << endl;
      exit(1);
    }

  GridOut::OutputFormat oformat;
  int c;
  const char* optstring="i:o:";
  while((c=getopt(argc, argv, optstring)) != -1)
  {
      switch (c)
      {
	  case 'i':
	      break;
	  case 'o':
	      oformat = GridOut::parse_output_format(optarg);
	      break;
      }
  }
  string iname =argv[optind];
//  if (iname.find_last_of('.') <= 0)
//    iname.append(".inp");
  
  GridOutFlags::XFig xfig_flags;
  xfig_flags.level_depth = false;

  Triangulation<2> tr;

   ifstream in (iname.c_str());
   AssertThrow(in, ExcFileNotOpen(argv[1]));

   GridIn<2> gin;
   gin.attach_triangulation(tr);
   gin.read_ucd(in);

   GridOut gout;
   gout.set_flags(xfig_flags);
   
   string oname(iname, 0, iname.find_last_of('.'));  
   oname.append(GridOut::default_suffix(oformat));
   cout << iname << " -> " << oname << endl;
   ofstream out(oname.c_str());

   gout.write(tr, out, oformat);
   
//  GridOut
}
