## ---------------------------------------------------------------------
##
## Copyright (C) 2013, 2015 by the deal.II authors
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

if ($#ARGV != 1) {
  print "\nUsage: make_gallery.pl gallery cmake_source_dir\n";
  exit;
}

$gallery=$ARGV[0];
$gallery_underscore=$gallery;
$gallery_underscore=~ s/-/_/;

$cmake_source_dir=$ARGV[1];

print
"/**
  * \@page code_gallery_$gallery_underscore The $gallery code gallery program
\@htmlonly
<table class=\"tutorial\" width=\"50%\">
<tr><th colspan=\"2\"><b><small>Table of contents</small></b></th></tr>
<tr><td width=\"50%\" valign=\"top\">
\@endhtmlonly
";

print
"*/
";
