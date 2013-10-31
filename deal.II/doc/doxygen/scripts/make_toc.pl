## ---------------------------------------------------------------------
## $Id$
##
## Copyright (C) 2006 - 2013 by the deal.II authors
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

$/ = undef;

# Read source for web page
open TOC, "<toc.html.in";
$toc = <TOC>;
close TOC;

# Read generated map file
open MAP, "<steps.cmapx";
$map = <MAP>;
close MAP;

# Insert contents of map file for @@MAP@@
$toc =~ s/\@\@MAP\@\@/$map/;

print $toc;
