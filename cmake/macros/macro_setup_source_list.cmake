## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2017 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

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
