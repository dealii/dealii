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

using namespace std;

int main(int argc, const char** argv)
{
  if (argc<2)
    {
      cerr << "Usage: " << argv[0] << " <filename>" << endl;
      exit(1);
    }

  string iname =argv[1];
//  if (iname.find_last_of('.') <= 0)
//    iname.append(".inp");
  
  Triangulation<2> tr;

   ifstream in (iname.c_str());
   AssertThrow(in, ExcFileNotOpen(argv[1]));

   GridIn<2> gin;
   gin.attach_triangulation(tr);
   gin.read_ucd(in);

   GridOut gout;
   
   string oname(iname, 0, iname.find_last_of('.'));  
   oname.append(".eps");
   cout << iname << " -> " << oname << endl;
   ofstream out(oname.c_str());

   gout.write_eps(tr, out);
   
//  GridOut
}
