#
# I'm sure we will have to split this file in various sensible pieces...
#


#
# Check for various compiler bugs:
#



#
# TODO: Obsolete. Remove and clean source.
# Versions of GCC before 3.0 had a problem with the explicit
# instantiation of member templates when the member was in fact
# an operator. In that case, they needed the "template" keyword,
# which is actually not allowed at this place. Test case is
#
# struct X
# {
#     template <typename T2>
#     X operator = (T2 &) { return X(); };
# };
#
# template X X::operator=<float> (float &);
#
# The compiler only groks this if the "operator=" is prepended
# by "template". We detect this, and either set the
# DEAL_II_MEMBER_OP_TEMPLATE_INST to "template" or nothing, so
# that it gets expanded to the right string needed in this place.
#
CHECK_CXX_SOURCE_COMPILES(
  "
  struct X
  {
      template <typename T2>
      X operator = (T2 &) { return X(); }
  };

  template X X::operator=<float> (float &);
  int main(){return 0;}
  "
  DEAL_II_MEMBER_OP_TEMPLATE_INST_OK)

IF(DEAL_II_MEMBER_OP_TEMPLATE_INST_OK)
  SET(DEAL_II_MEMBER_OP_TEMPLATE_INST "")
ELSE()
  SET(DEAL_II_MEMBER_OP_TEMPLATE_INST "template")
ENDIF()



#
# Some compiler versions, notably ICC, have trouble with the
# following code in which we explicitly call a destructor.
# This has to be worked around with a typedef. The problem is
# that the workaround fails with some other compilers, so that
# we can not unconditionally use the workaround...
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
# gcc 4.4 has an interesting problem in that it doesn't
# care for one of BOOST signals2's header files and produces
# dozens of pages of error messages of the form
#   warning: invoking macro BOOST_PP_CAT argument 1: \
#   empty macro arguments are undefined in ISO C90 and ISO C++98
# This can be avoided by not using -pedantic for this compiler.
# For all other versions, we use this flag, however.
#
IF(CMAKE_CXX_COMPILER_ID MATCHES "GNU" AND
   CMAKE_CXX_COMPILER_VERSION MATCHES "4.4.")
  STRIP_FLAG(CMAKE_C_FLAGS "-pedantic")
  STRIP_FLAG(CMAKE_CXX_FLAGS "-pedantic")
ENDIF()



