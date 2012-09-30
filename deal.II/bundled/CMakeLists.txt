
#
# TODO:
# The corresponding <feature>_folder and <feature_HAVE_BUNDLED variables
# are set under cmake/configure/configure_*. This is a bit unfortunate and
# should be corrected.
#

IF(FEATURE_TBB_BUNDLED_CONFIGURED)
  ADD_SUBDIRECTORY(${tbb_folder}/src)
ENDIF()


IF(FEATURE_BOOST_BUNDLED_CONFIGURED)
  INSTALL(DIRECTORY ${boost_folder}/include/boost
    DESTINATION ${DEAL_II_INCLUDE_RELDIR}/deal.II/bundled
    COMPONENT library
    FILES_MATCHING PATTERN "*.hh"
    PATTERN ".svn" EXCLUDE
    )

  ADD_SUBDIRECTORY(${boost_folder}/libs/serialization/src)

  IF( DEAL_II_USE_MT AND NOT DEAL_II_CAN_USE_CXX11)
    #
    # If the C++ compiler doesn't completely support the C++11 standard
    # (and consequently we can't use std::thread, std::mutex, etc), then
    # include all the files that form BOOST's thread implementation so that
    # we don't have to build BOOST itself only to get at this small part of
    # it. it also ensures that we use the correct compiler and flags
    #
    ADD_SUBDIRECTORY(${boost_folder}/libs/thread/src)
  ENDIF()
ENDIF()


IF(FEATURE_FUNCTIONPARSER_BUNDLED_CONFIGURED)
  ADD_SUBDIRECTORY(${functionparser_folder})
ENDIF()


IF(FEATURE_UMFPACK_BUNDLED_CONFIGURED)
  ADD_SUBDIRECTORY(${umfpack_folder}/UMFPACK/Source)
  ADD_SUBDIRECTORY(${umfpack_folder}/AMD/Source)
ENDIF()

