// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// This is the program that we use to generate explicit instantiations for a
// variety of template arguments. It takes two kinds of input files. The first
// is given as arguments on the command line and contains entries of the
// following form:
// --------------------
// REAL_SCALARS    := { double; float}
// COMPLEX_SCALARS := { std::complex<double>;
//                      std::complex<float>}
// VECTORS := { Vector<double>; Vector<float>}
// --------------------
//
// The input file is typically located in share/deal.II/template-arguments in
// the build directory and it is built from cmake/config/template-arguments.in
// to contain the list of vectors etc. that make sense for the current
// configuration. For example, the list of VECTORS is going to contain PETSc
// vectors if so configured.
//
// The second input is read from the command line and consists of a sequence
// of statements of the following form:
// --------------------
// for (u,v:VECTORS; z:SCALARS) { f(u, z, const v &); }
// --------------------
// Here, everything between {...} will be copied as many times as there are
// combinations of arguments u,v in the list of substitutions given by
// VECTORS. For each copy, the arguments u,v will be replaced by one of these
// combinations.

// Author: Wolfgang Bangerth, 2007


#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <string>

// Returns a map from the keys in the expansion lists to the list itself. For
// instance, the example above will lead to the entry
//      get_expansion_lists()[REAL_SCALARS] = (double, float)
// in this map, among others
std::map<std::string, std::list<std::string>> &
get_expansion_lists()
{
  static std::map<std::string, std::list<std::string>> expansion_lists;
  return expansion_lists;
}



/* ======================== auxiliary functions ================= */

// Return whether or not one string starts with a given prefix
bool
has_prefix(const std::string &base, const std::string &prefix)
{
  if (prefix.size() > base.size())
    return false;
  else
    return std::equal(prefix.begin(), prefix.end(), base.begin());
}


// replace all occurrences of 'pattern' by 'substitute' in 'in', and
// return the result
std::string
replace_all(const std::string &in,
            const std::string &pattern,
            const std::string &substitute)
{
  std::string x = in;
  while (x.find(pattern) != std::string::npos)
    x.replace(x.find(pattern), pattern.size(), substitute);
  return x;
}


// extract from the start of #in the part of the string that ends with one
// of the characters in #delim_list. The extracted part is deleted from #in
// and returned. We skip characters in delim_list if they are preceded by a
// backslash
std::string
get_substring_with_delim(std::string &in, const std::string &delim_list)
{
  std::string x;
  while (in.size() != 0)
    {
      // stop copying to the result if the current character is a
      // delimiter, but only if the previous character was not a backslash
      if ((delim_list.find(in[0]) != std::string::npos) &&
          !((x.size() > 0) && (x[x.size() - 1] == '\\')))
        break;

      x += in[0];
      in.erase(0, 1);
    }

  // We often end up with the '}' delimiter on a separate line, but
  // not in the first column of the line. Since whitespace isn't
  // harmful, that isn't a problem in itself, except that it makes
  // producing nicely formatted output a bit harder than necessary. It
  // would be nice if we could just end the text we read with the last
  // newline in such cases. To this end, just trim trailing
  // whitespace.
  while ((x.size() > 0) && (x.back() == ' '))
    x.erase(x.size() - 1, 1);

  return x;
}


// delete all whitespace at the beginning of the given argument
void
skip_space(std::string &in)
{
  while ((in.size() != 0) &&
         ((in[0] == ' ') || (in[0] == '\t') || (in[0] == '\n')))
    in.erase(0, 1);
}


std::string
remove_comments(std::string line)
{
  const std::string::size_type double_slash_comment = line.find("//");
  if (double_slash_comment != std::string::npos)
    line.erase(double_slash_comment, std::string::npos);

  const std::string::size_type slash_star_comment_begin = line.find("/*");
  if (slash_star_comment_begin != std::string::npos)
    {
      const std::string::size_type slash_star_comment_end = line.find("*/");
      if (slash_star_comment_end == std::string::npos)
        {
          std::cerr << "The program can currently only handle /* block */"
                    << "comments that start and end within the same line."
                    << std::endl;
          std::exit(1);
        }
      line.erase(slash_star_comment_begin,
                 slash_star_comment_end - slash_star_comment_begin + 2);
    }

  return line;
}


// read the whole file specified by the stream given as argument into a string
// for simpler parsing, and return it
std::string
read_whole_file(std::istream &in)
{
  std::string whole_file;
  while (in)
    {
      std::string line;
      getline(in, line);

      whole_file += remove_comments(line);
      whole_file += '\n';
    }
  // substitute tabs by spaces
  std::replace(whole_file.begin(), whole_file.end(), '\t', ' ');
  // substitute multiple spaces by single ones
  std::size_t position = 0;
  while ((position = whole_file.find("  ", position)) != std::string::npos)
    whole_file.replace(position, 2, 1, ' ');

  return whole_file;
}



// split a given string assumed to consist of a list of substrings
// delimited by a particular character into its components
std::list<std::string>
split_string_list(const std::string &s, const char delimiter)
{
  std::string            tmp = s;
  std::list<std::string> split_list;

  // split the input list
  while (tmp.length() != 0)
    {
      std::string name;
      name = tmp;

      if (name.find(delimiter) != std::string::npos)
        {
          name.erase(name.find(delimiter), std::string::npos);
          tmp.erase(0, tmp.find(delimiter) + 1);
        }
      else
        tmp = "";

      skip_space(name);

      while ((name.size() != 0) && (name[name.length() - 1] == ' ' ||
                                    name[name.length() - 1] == '\n'))
        name.erase(name.length() - 1, 1);

      split_list.push_back(name);
    }

  return split_list;
}



// return the given list but without empty entries
std::list<std::string>
delete_empty_entries(const std::list<std::string> &list)
{
  std::list<std::string> return_list;
  for (const auto &entry : list)
    if (!entry.empty())
      return_list.push_back(entry);

  return return_list;
}



// determine whether a given substring at position #pos and length #length
// in the string #text is a real token, i.e. not just part of another word
bool
is_real_token(const std::string           &text,
              const std::string::size_type pos,
              const std::string::size_type length)
{
  static const std::string token_chars("abcdefghijklmnopqrstuvwxyz"
                                       "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
                                       "0123456789"
                                       "_");
  if ((pos != 0) && (token_chars.find(text[pos - 1]) != std::string::npos))
    return false;

  if ((pos + length < text.size()) &&
      (token_chars.find(text[pos + length]) != std::string::npos))
    return false;

  return true;
}


// substitute all occurrences of #token in #text by #substitute. because a
// replacement token could be a templated class like std::complex<double>
// and because the token to the substituted may be a template argument
// itself, we surround the substitution by a space which shouldn't matter
// in C++
std::string
substitute_tokens(const std::string &text,
                  const std::string &token,
                  const std::string &substitute)
{
  std::string            x_text = text;
  std::string::size_type pos    = 0;
  while ((pos = x_text.find(token, pos)) != std::string::npos)
    {
      if (is_real_token(x_text, pos, token.size()))
        {
          x_text.replace(pos,
                         token.size(),
                         std::string(" ") + substitute + std::string(" "));
          pos += substitute.size() + 2;
        }
      else
        ++pos;
    }

  return x_text;
}



/* ======================== the main functions ================= */


// read and parse the expansion lists like
//   REAL_SCALARS    := { double; float}
// as specified at the top of the file and store them in the static
// expansion_lists variable
void
read_expansion_lists(const std::string &filename)
{
  std::ifstream in(filename.c_str());

  if (!in)
    {
      std::cerr << "Instantiation list file can not be read!" << std::endl;
      std::exit(1);
    }

  // read the entire file into a string for simpler processing. replace
  // end-of-line characters by spaces
  std::string whole_file = read_whole_file(in);

  skip_space(whole_file);

  // now process entries of the form
  //   NAME := { class1; class2; ...}.
  while (whole_file.size() != 0)
    {
      const std::string name = get_substring_with_delim(whole_file, " :");

      skip_space(whole_file);
      if (!has_prefix(whole_file, ":="))
        {
          std::cerr << "Invalid entry <" << name << '>' << std::endl;
          std::exit(1);
        }
      whole_file.erase(0, 2);
      skip_space(whole_file);
      if (whole_file[0] != '{')
        {
          std::cerr << "Invalid entry <" << name << '>' << std::endl;
          std::exit(1);
        }
      whole_file.erase(0, 1);
      skip_space(whole_file);

      std::string expansion = get_substring_with_delim(whole_file, "}");

      if (whole_file[0] != '}')
        {
          std::cerr << "Invalid entry <" << name << '>' << std::endl;
          std::exit(1);
        }
      whole_file.erase(0, 1);
      skip_space(whole_file);

      // assign but remove empty entries; this may happen if an expansion
      // list ends in a semicolon (then we get an empty entry at the end),
      // or if there are multiple semicolons after each other (this may
      // happen if, for example, we have "Vector<double>; TRILINOS_VECTOR;"
      // and if TRILINOS_VECTOR is an empty expansion after running
      // ./configure)
      get_expansion_lists()[name] =
        delete_empty_entries(split_string_list(expansion, ';'));
    }
}



// produce all combinations of substitutions of the tokens given in the
// #substitutions list in #text and output it to std::cout
void
substitute(const std::string                                    &text,
           const std::list<std::pair<std::string, std::string>> &substitutions)
{
  // do things recursively: if the list of substitutions has a single
  // entry, then process all of them. otherwise, process the first in the
  // list and call the function recursively with the rest of the
  // substitutions
  if (substitutions.size() > 1)
    {
      // do the first substitution, then call function recursively
      const std::string name    = substitutions.front().first,
                        pattern = substitutions.front().second;

      if (get_expansion_lists().find(pattern) == get_expansion_lists().end())
        {
          std::cerr << "could not find pattern '" << pattern << "'"
                    << std::endl;
          std::exit(1);
        }


      const std::list<std::pair<std::string, std::string>>
        rest_of_substitutions(++substitutions.begin(), substitutions.end());

      for (std::list<std::string>::const_iterator expansion =
             get_expansion_lists()[pattern].begin();
           expansion != get_expansion_lists()[pattern].end();
           ++expansion)
        {
          std::string new_text = substitute_tokens(text, name, *expansion);

          substitute(new_text, rest_of_substitutions);
        }
    }
  else if (substitutions.size() == 1)
    {
      // do the substitutions
      const std::string name    = substitutions.front().first,
                        pattern = substitutions.front().second;

      if (get_expansion_lists().find(pattern) == get_expansion_lists().end())
        {
          std::cerr << "could not find pattern '" << pattern << "'"
                    << std::endl;
          std::exit(1);
        }

      for (std::list<std::string>::const_iterator expansion =
             get_expansion_lists()[pattern].begin();
           expansion != get_expansion_lists()[pattern].end();
           ++expansion)
        {
          // surround each block in the for loop with an if-def hack
          // that allows us to split instantiation files into several
          // chunks to be used in different .cc files (to reduce
          // compiler memory usage).
          // Just define SPLIT_INSTANTIATIONS_COUNT to a positive number (number
          // of sections) to split the definitions into and
          // SPLIT_INSTANTIATIONS_INDEX as a number between 0 and
          // SPLIT_INSTANTIATIONS_COUNT-1 to get the instantiations of that
          // particular chunk.
          static unsigned int counter = 0;
          std::cout << "#if (SPLIT_INSTANTIATIONS_CHECK(" << counter++ << "))"
                    << std::endl;
          std::cout << substitute_tokens(text, name, *expansion);
          std::cout << "#endif" << std::endl << std::endl;
        }
    }
  else
    {
      std::cout << text << std::endl;
    }
}



// process the list of instantiations given in the form
//   for (u,v:VECTORS; z:SCALARS) { f(u, z, const v &); }
void
process_instantiations()
{
  std::string whole_file = read_whole_file(std::cin);

  // process entries of the form
  //   for (X:Y; A:B) { INST }
  while (whole_file.size() != 0)
    {
      // skip space, tabs, comments:
      skip_space(whole_file);

      // output preprocessor defines as is:
      if (has_prefix(whole_file, "#"))
        {
          std::cout << get_substring_with_delim(whole_file, "\n") << '\n';
          skip_space(whole_file);
          continue;
        }

      if (!has_prefix(whole_file, "for"))
        {
          std::cerr << "Invalid instantiation list: missing 'for'" << std::endl;
          std::exit(1);
        }
      whole_file.erase(0, 3);
      skip_space(whole_file);
      if (whole_file[0] != '(')
        {
          std::cerr << "Invalid instantiation list: missing '('" << std::endl;
          std::exit(1);
        }
      whole_file.erase(0, 1);
      skip_space(whole_file);

      const std::list<std::string> substitutions_list =
        split_string_list(get_substring_with_delim(whole_file, ")"), ';');
      if (whole_file[0] != ')')
        {
          std::cerr << "Invalid instantiation list: missing ')'" << std::endl;
          std::exit(1);
        }
      whole_file.erase(0, 1);
      skip_space(whole_file);

      // process the header
      std::list<std::pair<std::string, std::string>> substitutions;

      for (const auto &substitution : substitutions_list)
        {
          const std::list<std::string> names_and_type =
            split_string_list(substitution, ':');
          if (names_and_type.size() != 2)
            {
              std::cerr << "Invalid instantiation header: '" << substitution
                        << "'" << std::endl;
              std::exit(1);
            }

          const std::list<std::string> names =
            split_string_list(names_and_type.front(), ',');

          for (const auto &name : names)
            substitutions.emplace_back(name, names_and_type.back());
        }

      // now read the part in {...}
      skip_space(whole_file);
      if (whole_file[0] != '{')
        {
          std::cerr << "Invalid substitution text" << std::endl;
          std::exit(1);
        }
      whole_file.erase(0, 1);
      skip_space(whole_file);
      const std::string text_to_substitute =
        get_substring_with_delim(whole_file, "}");
      whole_file.erase(0, 1);
      skip_space(whole_file);

      // now produce the substitutions. first replace all occurrences of
      // "\{" by "{"
      substitute(replace_all(replace_all(text_to_substitute, "\\{", "{"),
                             "\\}",
                             "}"),
                 substitutions);
    }
}



int
main(int argc, char **argv)
{
  if (argc < 2)
    {
      std::cerr
        << "Usage: " << std::endl
        << "  expand_instantiations class_list_files < in_file > out_file"
        << std::endl;
      std::exit(1);
    }

  for (int i = 1; i < argc; ++i)
    read_expansion_lists(argv[i]);

  // write header:
  std::cout
    << "// This file is automatically generated from corresponding .inst.in, do not edit."
    << std::endl
    << std::endl;

  // Make sure SPLIT_INSTANTIATIONS_* is working correctly if the user doesn't
  // use it. The defaults will not split the instantiations. This logic is
  // somewhat tricky to get right for two reasons: 1. icc 14 will not allow an
  // expressition like "#if !defined(B) || (A % B == C)" 2. we have .cc files
  // where more than one .inst is included and splitting is only required in
  // one of them. So we need to handle the case where _COUNT is undefined but
  // _INDEX is defined, which might be needed later.
  std::cout
    << "#ifdef SPLIT_INSTANTIATIONS_COUNT" << std::endl
    << "  #define SPLIT_INSTANTIATIONS_CHECK(C) (((C) % SPLIT_INSTANTIATIONS_COUNT) == SPLIT_INSTANTIATIONS_INDEX)"
    << std::endl
    << "#else" << std::endl
    << "  #define SPLIT_INSTANTIATIONS_CHECK(C) (1)" << std::endl
    << "#endif" << std::endl
    << std::endl;

  process_instantiations();

  // undefine the macro to avoid issues when more than one .inst file is
  // included in a single .cc
  std::cout << std::endl << "#undef SPLIT_INSTANTIATIONS_CHECK" << std::endl;
}
