#####
##
## Copyright (C) 2012 by the deal.II authors
##
## This file is part of the deal.II library.
##
## <TODO: Full License information>
## This file is dual licensed under QPL 1.0 and LGPL 2.1 or any later
## version of the LGPL license.
##
#####

#
# Check for various compiler bugs:
#


#
# Some compiler versions, notably ICC, have trouble with the
# following code in which we explicitly call a destructor.
# This has to be worked around with a typedef. The problem is
# that the workaround fails with some other compilers, so that
# we can not unconditionally use the workaround...
#
# - Matthias Maier, rewritten 2012
#
CHECK_CXX_COMPILER_BUG(
  "
  namespace dealii
  {
    namespace FEValuesViews
    {
      template <int dim, int spacedim> struct Scalar {};
    }

    template <int dim, int spacedim>
    struct X
    {
        FEValuesViews::Scalar<dim,spacedim> scalars[dim*spacedim];

        void f()
          {
            scalars[0].dealii::FEValuesViews::Scalar<dim,spacedim>::~Scalar ();
          }
    };

    template struct X<2,2>;
  }
  int main() { return 0; }
  "
  DEAL_II_EXPLICIT_DESTRUCTOR_BUG)


#
# On some gcc 4.3 snapshots, a 'const' qualifier on a return type triggers a
# warning. This is unfortunate, since we happen to stumble on this
# in some of our template trickery with iterator classes. If necessary,
# do not use the relevant warning flag
#
# - Matthias Maier, rewritten 2012
#
PUSH_TEST_FLAG("-Wreturn-type")
PUSH_TEST_FLAG("-Werror")
CHECK_CXX_COMPILER_BUG(
  "
  const double foo() { return 1.; }
  int main() { return 0; }
  "
  DEAL_II_WRETURN_TYPE_CONST_QUALIFIER_BUG)
POP_TEST_FLAG()
POP_TEST_FLAG()

IF(DEAL_II_WRETURN_TYPE_CONST_QUALIFIER_BUG)
  ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS -Wno-return-type)
  ENABLE_IF_SUPPORTED(CMAKE_C_FLAGS -Wno-return-type)
ENDIF()


#
# gcc 4.4 has an interesting problem in that it doesn't
# care for one of BOOST signals2's header files and produces
# dozens of pages of error messages of the form
#   warning: invoking macro BOOST_PP_CAT argument 1: \
#   empty macro arguments are undefined in ISO C90 and ISO C++98
# This can be avoided by not using -pedantic for this compiler.
# For all other versions, we use this flag, however.
#
# - Matthias Maier, rewritten 2012
#
IF(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND
   CMAKE_CXX_COMPILER_VERSION MATCHES "4.4.")
  STRIP_FLAG(CMAKE_CXX_FLAGS "-pedantic")
  STRIP_FLAG(CMAKE_C_FLAGS "-pedantic")
ENDIF()


#
# Some gcc compiler versions have a problem when using an unsigned count
# in the std::advance function. Unfortunately, this also happens
# occasionally from within the standard library, so we can't prevent the
# warning messages. Since this is annoying, switch of the flag -W which
# causes this.
#
# - Matthias Maier, rewritten 2012
#

# TODO: We use the mpi.h header file for this check. We should test this
# bug with another header file than mpi.h ...
CHECK_INCLUDE_FILE_CXX("mpi.h" HAVE_MPI_H)

IF(HAVE_MPI_H)
  PUSH_TEST_FLAG("-Wunused-parameter")
  PUSH_TEST_FLAG("-Werror")
  CHECK_CXX_COMPILER_BUG(
    "
    #include <mpi.h>
    int main() { return 0; }
    "
    DEAL_II_ADVANCE_WARNING_BUG)
  POP_TEST_FLAG()
  POP_TEST_FLAG()

  IF(DEAL_II_ADVANCE_WARNING_BUG)
    ENABLE_IF_SUPPORTED(CMAKE_CXX_FLAGS "-Wno-unused-parameter")
    ENABLE_IF_SUPPORTED(CMAKE_C_FLAGS "-Wno-unused-parameter")
  ENDIF()
ENDIF()


#
# Disable shared libraries on cygwin targets for the moment. Shared
# libraries on Windows are a bit buggy atm.
#
# - Matthias Maier, 2012
#
IF(CMAKE_SYSTEM_NAME MATCHES "CYGWIN")
  MESSAGE(WARNING "\n"
    "BUILD_SHARED_LIBS forced to OFF\n\n"
    )
  SET(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)
ENDIF()

