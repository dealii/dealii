## ---------------------------------------------------------------------
##
## Copyright (C) 2017 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
##
## ---------------------------------------------------------------------

#
# Usage:
#   setup_source_list("${_list_of_unity_source_files}"
#                     "${_list_of_non_unity_source_files}"
#                     _includes_per_unity_file
#                     _source_files)
#
# This macro sets up the list of source files for the current directory. If
# DEAL_II_UNITY_BUILD=ON then this script calls SETUP_UNITY_TARGET to do so;
# otherwise, if DEAL_II_UNITY_BUILD=OFF, then all source files are added to the
# variable _source_files.
#
macro(setup_source_list _unity_include_src _separate_src _n_includes_per_unity_file _output_src)
  if(DEAL_II_UNITY_BUILD)
    setup_unity_target("${_unity_include_src}" ${_n_includes_per_unity_file} ${_output_src})
    set(${_output_src}
      ${${_output_src}}
      ${_separate_src}
      )
  else()
    set(${_output_src}
      ${${_output_src}}
      ${_unity_include_src}
      ${_separate_src}
      )
  endif()
endmacro()
