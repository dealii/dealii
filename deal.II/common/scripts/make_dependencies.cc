#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <set>
#include <vector>
#include <cassert>

                                 // base path for object files
std::string basepath;

                                 // list of include directories
std::vector<std::string> include_directories;

                                 // for each file that we ever visit,
                                 // store the set of other files it
                                 // includes directly
std::map<std::string,std::set<std::string> > direct_includes;


                                 // for the given file, fill a
                                 // respective entry in the "direct_includes"
                                 // map listing the names of those
                                 // files that are directly included
void determine_direct_includes (const std::string &file) 
{
                                   // if this file has already been
                                   // treated, then leave it at this
  if (direct_includes.find (file) != direct_includes.end())
    return;

                                   // otherwise, open the file and go
                                   // through it line by line to
                                   // search for other includes. we
                                   // will have to have the path to
                                   // the present file later, so get
                                   // it already here
  const std::string present_path (file.find ('/') != std::string::npos
                                  ?
                                  std::string (file.begin(),
                                               file.begin()+file.rfind ('/')+1)
                                  :
                                  "");
  std::ifstream in(file.c_str());
  assert (in);

  std::string line;
  while (in) 
    {
                                       // get one line, eat whitespace
                                       // at the beginning and see
                                       // whether the first
                                       // non-whitespace is a #
                                       // character
      getline (in, line);
      unsigned int pos=0;
      for (; pos<line.length(); ++pos)
        if ((line[pos] != ' ') && (line[pos] != '\t'))
          break;

                                       // if no non-whitespace, or
                                       // something other than #: next
                                       // line
      if ((pos == line.length()) || (line[pos] != '#'))
        continue;

                                       // ok, this is a preprocessor
                                       // line. eat pound sign and
                                       // again the next couple of
                                       // whitespaces
      ++pos;
      for (; pos<line.length(); ++pos)
        if ((line[pos] != ' ') && (line[pos] != '\t'))
          break;

                                       // and let's see whether the
                                       // following is the word
                                       // include
      if ((line.length() < pos+7)
          ||
          (! ((line[pos+0] == 'i') &&
              (line[pos+1] == 'n') &&
              (line[pos+2] == 'c') &&
              (line[pos+3] == 'l') &&
              (line[pos+4] == 'u') &&
              (line[pos+5] == 'd') &&
              (line[pos+6] == 'e'))))
        continue;

                                       // ok, word found. advance pos
                                       // and eat more whitespace
      pos += 7;
      for (; pos<line.length(); ++pos)
        if ((line[pos] != ' ') && (line[pos] != '\t'))
          break;

                                       // check that the next char is
                                       // either '<' or '"'
      if ((line[pos] != '"') && (line[pos] != '<'))
        continue;

                                       // copy out name
      std::string included_file;
      for (unsigned int endpos=pos+1; endpos<line.length(); ++endpos)
        if ((line[endpos]=='"') || (line[endpos] == '>'))
          {
            included_file = std::string (line.begin()+pos+1,
                                         line.begin()+endpos);
            break;
          }
      assert (included_file.length() > 0);

                                       // next try to locate the file
                                       // in absolute paths. this is
                                       // easy if it was included via
                                       // "...", but for <...> we have
                                       // to work a little harder
      if (included_file[0] != '/')
        {
          if (line[pos] == '"')
            included_file = present_path+included_file;
          else
            for (std::vector<std::string>::const_iterator
                   include_dir=include_directories.begin();
                 include_dir!=include_directories.end(); ++include_dir)
              if (std::ifstream((*include_dir+included_file).c_str()))
                {
                  included_file = *include_dir+included_file;
                  break;
                }
        }
      
                                       // make sure the file
                                       // exists, otherwise just
                                       // ignore the line
      if (!std::ifstream(included_file.c_str()))
        continue;

                                       // ok, so we did find an
                                       // appropriate file. add it to
                                       // the correct list
      direct_includes[file].insert (included_file);

                                       // work on the include file
                                       // recursively. note that the
                                       // first line of this file
                                       // saves us from infinite
                                       // recursions in case of
                                       // include loops
      determine_direct_includes (included_file);
    }
}



                                 // return the set of all included
                                 // files, directly or indirectly, for
                                 // the given file
std::set<std::string>
get_all_includes (const std::string &name) 
{
                                   // start with direct includes
  std::set<std::string> all_includes = direct_includes[name];

  std::set<std::string> next_level_includes = all_includes;
  while (true)
    {
                                       // traverse all next level
                                       // includes and get their
                                       // direct include files. the
                                       // set makes sure that
                                       // duplicates are removed
      std::set<std::string> second_next;
      for (std::set<std::string>::const_iterator
             next=next_level_includes.begin();
           next!=next_level_includes.end(); ++next)
        second_next.insert (direct_includes[*next].begin(),
                            direct_includes[*next].end());

                                       // for each of them, if it
                                       // hasn't been treated then add
                                       // it to the files of the next
                                       // level and later add it to
                                       // the all_includes
      next_level_includes.clear ();
      for (std::set<std::string>::const_iterator f=second_next.begin();
           f!=second_next.end(); ++f)
        if (all_includes.find(*f) == all_includes.end())
          next_level_includes.insert (*f);

                                       // if no new includes found no
                                       // more, then quit
      if (next_level_includes.size() == 0)
        return all_includes;
      
                                       // otherwise, copy over and
                                       // start over on the next level
                                       // of the tree
      all_includes.insert (next_level_includes.begin(),
                           next_level_includes.end());
    }
}



int main (int argc, char **argv) 
{
  std::vector<std::string> filenames;

                                   // parse all arguments (except the
                                   // name of the executable itself)
  for (unsigned int c=1; c<argc; ++c)
    {
      const std::string arg = argv[c];

                                       // if string starts with -I,
                                       // take this as an include path
      if ((arg.length()>2) && (arg[0]=='-') && (arg[1]=='I'))
        {
          std::string dir (arg.begin()+2, arg.end());

                                           // append a slash if not
                                           // already there
          if (dir[dir.length()-1] != '/')
            dir += '/';

                                           // drop initial ./ if this
                                           // is there
          if ((dir[0]=='.') && (dir[1]=='/'))
            dir = std::string(dir.begin()+2, dir.end());
          
          include_directories.push_back (dir);
        }
                                       // if string starts with -B,
                                       // then this is the base name
                                       // for object files
      else if ((arg.length()>2) && (arg[0]=='-') && (arg[1]=='B'))
        basepath = std::string(arg.begin()+2, arg.end());

                                       // otherwise assume that this
                                       // is one of the files for
                                       // input
      else
        {
          assert (arg.size()>=1);
          assert (arg[0] != '-');
          
          filenames.push_back (arg);
        }
    }

                                   // next iterate through all files
                                   // and figure out which other files
                                   // they include
  for (std::vector<std::string>::const_iterator file=filenames.begin();
       file != filenames.end(); ++file)
    determine_direct_includes (*file);


                                   // now we have all files that are
                                   // directly or indirectly included
                                   // into the files given on the
                                   // command lines. for each of them,
                                   // we have recorded which files
                                   // they include themselves. for
                                   // each file on the command line,
                                   // we can thus form a complete
                                   // include tree. do exactly this by
                                   // populating the all_includes
                                   // variable. we walk over the list
                                   // of all these files, but since
                                   // the tree is constructed
                                   // recursively, it may happen that
                                   // for some later filenames the
                                   // all_includes is already
                                   // built. the function then returns
                                   // immediately
//   for (std::map<std::string,std::set<std::string> >::const_iterator
//          file = direct_includes.begin(); file!=direct_includes.end(); ++file)
//     complete_tree (file->first);

  for (std::vector<std::string>::const_iterator file=filenames.begin();
       file != filenames.end(); ++file)
    {
                                       // get base of filename by
                                       // chipping away .cc extension
                                       // as well as path
      std::string basename;
      if (file->find (".cc") != std::string::npos)
        basename = std::string (file->begin(),
                                file->begin()+file->find (".cc"));
      else if (file->find (".cpp") != std::string::npos)
        basename = std::string (file->begin(),
                                file->begin()+file->find (".cpp"));
      else
        basename = *file;

      if (basename.rfind ("/") != std::string::npos)
        basename = std::string (basename.begin()+basename.rfind("/")+1,
                                basename.end());
      
      std::cout << basepath << "/" << basename << ".o "
                << basepath << "/" << basename << ".g.o: \\"
                << std::endl
                << "\t\t" << *file;
      
      const std::set<std::string> includes = get_all_includes (*file);
      for (std::set<std::string>::const_iterator i=includes.begin();
           i!=includes.end(); ++i)
        std::cout << "\\\n\t\t" << *i;
      std::cout << std::endl;
    }
}

