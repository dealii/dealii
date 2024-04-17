## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2012 - 2023 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

set(src_boost_serialization
  archive_exception.cpp
  basic_archive.cpp
  basic_iarchive.cpp
  basic_iserializer.cpp
  basic_oarchive.cpp
  basic_oserializer.cpp
  basic_pointer_iserializer.cpp
  basic_pointer_oserializer.cpp
  basic_serializer_map.cpp
  basic_text_iprimitive.cpp
  basic_text_oprimitive.cpp
  basic_text_wiprimitive.cpp
  basic_text_woprimitive.cpp
  basic_xml_archive.cpp
  binary_iarchive.cpp
  binary_oarchive.cpp
  binary_wiarchive.cpp
  binary_woarchive.cpp
  codecvt_null.cpp
  extended_type_info.cpp
  extended_type_info_no_rtti.cpp
  extended_type_info_typeid.cpp
  polymorphic_iarchive.cpp
  polymorphic_oarchive.cpp
  stl_port.cpp
  text_iarchive.cpp
  text_oarchive.cpp
  text_wiarchive.cpp
  text_woarchive.cpp
  utf8_codecvt_facet.cpp
  void_cast.cpp
  xml_archive_exception.cpp
  xml_grammar.cpp
  xml_iarchive.cpp
  xml_oarchive.cpp
  xml_wgrammar.cpp
  xml_wiarchive.cpp
  xml_woarchive.cpp
  )
enable_if_supported(DEAL_II_WARNING_FLAGS -Wno-c11-extensions)
enable_if_supported(DEAL_II_WARNING_FLAGS -Wno-deprecated-copy)
enable_if_supported(DEAL_II_WARNING_FLAGS -Wno-uninitialized)
enable_if_supported(DEAL_II_WARNING_FLAGS -Wno-deprecated-builtins)
enable_if_supported(DEAL_II_WARNING_FLAGS -Wno-deprecated-declarations)
enable_if_supported(DEAL_II_WARNING_FLAGS -Wno-suggest-override)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../../include)
define_object_library(bundled_boost_serialization OBJECT ${src_boost_serialization})
