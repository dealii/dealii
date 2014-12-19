## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2012 - 2014 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

INCLUDE_DIRECTORIES(${BOOST_BUNDLED_INCLUDE_DIRS})

SET(src_boost_serialization
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
  basic_xml_archive.cpp
  binary_iarchive.cpp
  binary_oarchive.cpp
  extended_type_info.cpp
  extended_type_info_typeid.cpp
  extended_type_info_no_rtti.cpp
  polymorphic_iarchive.cpp
  polymorphic_oarchive.cpp
  stl_port.cpp
  text_iarchive.cpp
  text_oarchive.cpp
  void_cast.cpp
  archive_exception.cpp
  xml_grammar.cpp
  xml_iarchive.cpp
  xml_oarchive.cpp
  xml_archive_exception.cpp
  )

DEAL_II_ADD_LIBRARY(obj_boost_serialization OBJECT ${src_boost_serialization})

