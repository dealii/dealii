## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------


# This file defines some common infrastructure for the two scripts
# that convert header and source files into module partition
# units. Specifically, it contains lists of header files for each
# external dependency of deal.II, and compiles them into regular
# expressions that can be matched against #include statements.


import re

# Define lists of (potentially wildcarded) names of header files for
# deal.II and all other external dependencies we use.. The list of C++
# headers is from
# https://stackoverflow.com/questions/2027991/list-of-standard-header-files-in-c-and-c
dealii_includes = ["deal.II/(.*).h"]

adolc_includes = ["adolc/.*"]

boost_includes = ["boost/.*"]

cgal_includes = ["CGAL/.*"]

hdf5_includes = ["hdf5.h"]

kokkos_includes = ["Kokkos.*"]

metis_includes = ["metis.h"]

mpi_includes = ["mpi.h"]

mumps_includes = ["mumps_.*.h"]

muparser_includes = ["muParser.h"]

opencascade_includes = [
    "Adaptor3d_Curve.hxx",
    "Adaptor3d_HCurve.hxx",
    "BRep.*.hxx",
    "Geom.*.hxx",
    "Handle.*.hxx",
    "IGESControl.*.hxx",
    "IFSelect_ReturnStatus.hxx",
    "IntCurve.*.hxx",
    "Poly_Triangulation.hxx",
    "ShapeAnalysis.*.hxx",
    "Standard_Transient.hxx",
    "STEPControl.*.hxx",
    "StlAPI.*.hxx",
    "TCol.*.hxx",
    "Top.*.hxx",
    "TopoDS.*.hxx",
    "gp_.*.hxx",
]

p4est_includes = ["p4est.*", "p8est.*", "sc_containers.h"]

slepc_includes = ["slepc.*.h"]

petsc_includes = [
    "petscconf.h",
    "petscdm.h",
    "petscerror.h",
    "petscis.h",
    "petscistypes.h",
    "petscksp.h",
    "petscmat.h",
    "petscpc.h",
    "petscsf.h",
    "petscsnes.h",
    "petscsys.h",
    "petscts.h",
    "petscvec.h",
    "petsc/.*.h",
]

std_includes = [
    "assert.h",
    "limits.h",
    "signal.h",
    "stdlib.h",
    "ctype.h",
    "locale.h",
    "stdarg.h",
    "string.h",
    "errno.h",
    "math.h",
    "stddef.h",
    "time.h",
    "float.h",
    "setjmp.h",
    "stdio.h",
    "algorithm",
    "format",
    "new",
    "stdexcept",
    "any",
    "forward_list",
    "numbers",
    "stdfloat",
    "array",
    "fstream",
    "numeric",
    "stop_token",
    "atomic",
    "functional",
    "optional",
    "streambuf",
    "barrier",
    "future",
    "ostream",
    "string",
    "bit",
    "generator",
    "print",
    "string_view",
    "bitset",
    "initializer_list",
    #                    "queue",   # excluded for now due to https://github.com/llvm/llvm-project/issues/138558
    "strstream",
    "charconv",
    "iomanip",
    #                    "random",   # excluded for now due to https://github.com/llvm/llvm-project/issues/138558
    "syncstream",
    "chrono",
    "ios",
    "ranges",
    "system_error",
    "codecvt",
    "iosfwd",
    "ratio",
    "thread",
    "compare",
    "iostream",
    "regex",
    "tuple",
    "complex",
    "istream",
    "scoped_allocator",
    "type_traits",
    "concepts",
    "iterator",
    "semaphore",
    "typeindex",
    "condition_variable",
    "latch",
    "set",
    "typeinfo",
    "coroutine",
    "limits",
    "shared_mutex",
    "unordered_map",
    "deque",
    "list",
    "source_location",
    "unordered_set",
    "exception",
    "locale",
    "span",
    "utility",
    "execution",
    "map",
    "spanstream",
    "valarray",
    "expected",
    "mdspan",
    "sstream",
    "variant",
    "filesystem",
    "memory",
    #                    "stack",   # excluded for now due to https://github.com/llvm/llvm-project/issues/138558
    "vector",
    "flat_map",
    "memory_resource",
    "stacktrace",
    "version",
    "flat_set",
    "mutex",
    "cassert",
    "cfenv",
    "climits",
    "csetjmp",
    "cstddef",
    "cstdlib",
    "cuchar",
    "cctype",
    "cfloat",
    "clocale",
    "csignal",
    "cstdint",
    "cstring",
    "cwchar",
    "cerrno",
    "cinttypes",
    "cmath",
    "cstdarg",
    "cstdio",
    "ctime",
    "cwctype",
]

sundials_includes = [
    "arkode/arkode.*",
    "ida/ida.h",
    "idas/idas.h",
    "kinsol/kinsol.*",
    "nvector/nvector_parallel.h",
    "nvector/nvector_serial.h",
    "sundials/sundials_.*",
    "sunlinsol/sunlinsol_.*",
    "sunmatrix/sunmatrix_.*",
    "sunnonlinsol/sunnonlinsol_.*",
]

taskflow_includes = ["taskflow/.*"]

tbb_includes = ["tbb/.*"]

trilinos_includes = [
    "Amesos.*",
    "AztecOO.h",
    "Belos.*",
    "Epetra.*",
    "Ifpack2.*",
    "NOX.*",
    "Rol.*",
    "ROL.*",
    "Sacado.*",
    "Teuchos.*",
    "Tpetra.*",
    "zoltan_cpp.h",
    "exodusII.h",
]

umfpack_includes = ["umfpack.h"]

vtk_includes = ["vtk.*"]

zlib_includes = ["zlib.h"]

# Define regexs to match header #includes for each of the categories above
match_dealii_includes = re.compile("# *include *<" + "|".join(dealii_includes) + ">")

match_adolc_includes = re.compile("# *include *<(" + "|".join(adolc_includes) + ")>")
match_boost_includes = re.compile("# *include *<(" + "|".join(boost_includes) + ")>")
match_cgal_includes = re.compile("# *include *<(" + "|".join(cgal_includes) + ")>")
match_hdf5_includes = re.compile("# *include *<(" + "|".join(hdf5_includes) + ")>")
match_kokkos_includes = re.compile("# *include *<(" + "|".join(kokkos_includes) + ")>")
match_metis_includes = re.compile("# *include *<(" + "|".join(metis_includes) + ")>")
match_mpi_includes = re.compile("# *include *<(" + "|".join(mpi_includes) + ")>")
match_mumps_includes = re.compile("# *include *<(" + "|".join(mumps_includes) + ")>")
match_muparser_includes = re.compile(
    "# *include *<(" + "|".join(muparser_includes) + ")>"
)
match_opencascade_includes = re.compile(
    "# *include *<(" + "|".join(opencascade_includes) + ")>"
)
match_p4est_includes = re.compile("# *include *<(" + "|".join(p4est_includes) + ")>")
match_petsc_includes = re.compile("# *include *<(" + "|".join(petsc_includes) + ")>")
match_slepc_includes = re.compile("# *include *<(" + "|".join(slepc_includes) + ")>")
match_std_includes = re.compile("# *include *<(" + "|".join(std_includes) + ")>")
match_sundials_includes = re.compile(
    "# *include *<(" + "|".join(sundials_includes) + ")>"
)
match_taskflow_includes = re.compile(
    "# *include *<(" + "|".join(taskflow_includes) + ")>"
)
match_tbb_includes = re.compile("# *include *<(" + "|".join(tbb_includes) + ")>")
match_trilinos_includes = re.compile(
    "# *include *<(" + "|".join(trilinos_includes) + ")>"
)
match_umfpack_includes = re.compile(
    "# *include *<(" + "|".join(umfpack_includes) + ")>"
)
match_vtk_includes = re.compile("# *include *<(" + "|".join(vtk_includes) + ")>")
match_zlib_includes = re.compile("# *include *<(" + "|".join(zlib_includes) + ")>")

external_package_headers_regex_map = {
    "adolc": match_adolc_includes,
    "boost": match_boost_includes,
    "cgal": match_cgal_includes,
    "hdf5": match_hdf5_includes,
    "kokkos": match_kokkos_includes,
    "metis": match_metis_includes,
    "mpi": match_mpi_includes,
    "mumps": match_mumps_includes,
    "muparser": match_muparser_includes,
    "opencascade": match_opencascade_includes,
    "p4est": match_p4est_includes,
    "petsc": match_petsc_includes,
    "slepc": match_slepc_includes,
    "std": match_std_includes,
    "sundials": match_sundials_includes,
    "taskflow": match_taskflow_includes,
    "tbb": match_tbb_includes,
    "trilinos": match_trilinos_includes,
    "umfpack": match_umfpack_includes,
    "vtk": match_vtk_includes,
    "zlib": match_zlib_includes,
}


# A function that given a list of strings (specifically: the
# individual lines of a file) returns the list of deal.II header files
# that are being #included, excluding "config.h" and "exception_macros.h".
def get_dealii_includes(lines):
    dealii_include_list = []
    for line in lines:
        m = match_dealii_includes.match(line)
        if m:
            if (m.group(1) != "base/config") and (
                m.group(1) != "base/exception_macros"
            ):
                dealii_include_list.append(m.group(1))
    return dealii_include_list


# Return whether a line in a header or source file matches an #include
# statement for any of the recognized external projects that we wrap
# in module partitions:
def matches_external_header_include(line):
    for package, regex in external_package_headers_regex_map.items():
        if regex.match(line):
            return True
    return False


# A function that given a list of strings (specifically: the
# individual lines of a file) returns which of the external projects
# we wrap in module partitions are actually being used by this one
# file (as indicated by this file #including any of the corresponding
# header files of that external project).
def get_used_external_projects(lines):
    used_external_projects = set()
    my_external_package_headers_regex_map = external_package_headers_regex_map.copy()

    for line in lines:
        for package, regex in my_external_package_headers_regex_map.items():
            if regex.match(line):
                # The current line matches an external package's
                # headers. Once we have determined that a file uses an
                # external package, we no longer need to test for that and
                # can remove that entry from the list of regexes. Note that
                # for this to work, we need to use .copy() above when creating
                # the my_... dictionary, as we would otherwise simply be
                # deleting the entry in the original dictionary.
                used_external_projects.add(package)
                del my_external_package_headers_regex_map[package]
                break

    return used_external_projects
