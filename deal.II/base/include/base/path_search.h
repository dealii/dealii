//-----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------
#ifndef __deal2__path_search_h
#define __deal2__path_search_h


#include <string>
#include <fstream>
#include <map>
#include <vector>
#include <memory>

#include <base/config.h>
#include <base/exceptions.h>

/**
 * Support for searching files in a list of paths and with a list of
 * suffixes.
 *
 * A list of search paths is maintained for each file class supported.
 * A file class is defined by s unique string. The classes provided are
 * <dl>
 * <dt> MESH <dd> mesh input files in various formats (see GridIn)
 * <dt> PARAMETER <dd> Parameter files (<tt>.prm</tt>)
 * </dl>
 *
 * Additional file classes can be added easily by using add_class().
 *
 * Usage: First, you construct a PathSearch object for a certain file class,
 * e.g. meshes. Then, you use the open() method to obtain a reference
 * to an <tt>istream</tt> object.
 * @verbatim
 * #include <base/path_search.h>
 *
 * PathSearch search("MESH");
 * istream& in = search.open("grid");
 * ...
 * @endverbatim
 *
 * This piece of code will first traverse all paths in the list set up
 * for file class <TT>MESH</tt>. If it manages to open a file, it
 * returns the <tt>istream</tt> object. If not, it will try to append
 * the first suffix of the suffix list and do the same. And so on. If
 * no file is found in the end, an exception is thrown.
 *
 * If you want to restrict your search to a certain mesh format,
 * <tt>.inp</tt> for instance, then either use <tt>"grid.inp"</tt> in
 * the code above or use the alternative open function
 * @verbatim
 * istream& in = search.open("grid", ".inp");
 * @endverbatim
 *
 * Path lists are by default starting with the current directory
 * (<tt>""</tt>), followed optionally by a standard directory of
 * deal.II. Use show() to find out the path list
 * for a given class. Paths and suffixes can be added using the
 * functions add_path() and add_suffix(), repectively.
 *
 * @note Directories in the path list should always end with a
 * trailing <tt>"/"</tt>, while suffixes should always start with a
 * dot. These characters are not added automatically (allowing you to
 * do some real file name editing).
 *
 * @todo Add support for environment variables like in kpathsea.
 *
 * @author Guido Kanschat, 2005
 */

class PathSearch
{
  public:
				     /**
				      * Position for adding a new item to a list.
				      */
    enum Position
    {
					   /// Add new item at end of list
	  back,
					   /// Add new item at front of list
	  front,
					   /// Add in path list after empty element
	  after_none
    };
    
				     /**
				      * Constructor. The first argument is a
				      * string identifying the class of files
				      * to be searched for.
				      *
				      * The debug argument determines
				      * the verbosity of this class.
				      */
    PathSearch(const std::string& cls, const unsigned int debug=0);

				     /**
				      * Open a file in the class
				      * specified by the constructor.
				      */
    std::istream& open(const std::string& filename);

				     /**
				      * Open a file in the class
				      * specified by the
				      * constructor. Do not use the
				      * standard suffix list, but only
				      * try to apply the suffix given.
				      */
    std::istream& open(const std::string& filename,
		       const std::string& suffix);
    
				     /**
				      * Show the paths and suffixes
				      * used for this object.
				      */
    template <class STREAM>
    void show(STREAM& stream) const;
    
				     /**
				      * Add a new class.
				      */
    static void add_class (const std::string& cls);
    
				     /**
				      * Add a path to the current
				      * class. See
				      * PathSearch::Position for
				      * possible position arguments.
				      */
    void add_path (const std::string& path, Position pos = back);
    
				     /**
				      * Add a path to the current
				      * class. See
				      * PathSearch::Position for
				      * possible position arguments.
				      */
    void add_suffix (const std::string& suffix, Position pos = back);

				     /**
				      * This class was not
				      * registered in the path
				      * search mechanism.
				      */
    DeclException1(ExcNoClass,
		   std::string,
		   << "The class "
		   << arg1
		   << " must be registered before referring it in PathSearch");
    
    DeclException2(ExcFileNotFound,
		   std::string, std::string,
		   << "The file \"" << arg1
		   << "\"was not found in the path for files of class "
		   << arg2);
    
  private:
				     /**
				      * Type of values in the class maps.
				      */
    typedef std::map<std::string, std::vector<std::string> >::value_type map_type;
    
				     /**
				      * Initialize the static list objects for further use.
				      */
    static void initialize_classes();
    
				     /**
				      * Get path list for a certain
				      * class. Used to set up
				      * #my_path_list in constructor.
				      */
    static std::vector<std::string>& get_path_list(const std::string& cls);
    
				     /**
				      * Get suffix list for a certain
				      * class. Used to set up
				      * #my_suffix_list in constructor.
				      */
    static std::vector<std::string>& get_suffix_list(const std::string& cls);
    
				     /**
				      * The flie class handled by this object.
				      */
    const std::string cls;
    
				     /**
				      * All path lists for all
				      * classes, such that we can
				      * build them only once.
				      */
    static std::map<std::string, std::vector<std::string> > path_lists;
    
				     /**
				      * List of suffixes for each class.
				      */
    static std::map<std::string, std::vector<std::string> > suffix_lists;
    
				     /**
				      * Path list for the class this object belongs to.
				      */
    std::vector<std::string>& my_path_list;
    
				     /**
				      * Suffix list for the class this object belongs to.
				      */
    std::vector<std::string>& my_suffix_list;
    
				     /**
				      * The file stream after open was called.
				      */
    std::auto_ptr<std::ifstream> stream;

				     /**
				      * Debug flag. No output if zero.
				      */
    const unsigned int debug;
    
				     /**
				      * The empty string.
				      */
    static std::string empty;
};


template <class STREAM>
void
PathSearch::show(STREAM& out) const
{
  out << "DEAL_II_" << cls << "PATH=\"";
  bool first = true;
  for (std::vector<std::string>::iterator p = my_path_list.begin();
       p != my_path_list.end(); ++p)
    {
      if (!first)
	out << ':';
      out << *p;
      first = false;
    }
  out << '"' << std::endl << " Suffixes";
  for (std::vector<std::string>::iterator s = my_suffix_list.begin();
       s != my_suffix_list.end(); ++s)
    out << " \"" << *s << '"';
  out << std::endl;
}


#endif

