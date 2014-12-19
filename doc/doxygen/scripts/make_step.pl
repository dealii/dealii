## ---------------------------------------------------------------------
##
## Copyright (C) 2013 by the deal.II authors
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
  print "\nUsage: make_step.pl step cmake_source_dir\n";
  exit;
}

$step=$ARGV[0];
$step_underscore=$step;
$step_underscore=~ s/-/_/;

$cmake_source_dir=$ARGV[1];

print
"/**
  * \@page $step_underscore The $step tutorial program
\@htmlonly
<table class=\"tutorial\" width=\"50%\">
<tr><th colspan=\"2\"><b><small>Table of contents</small></b></th></tr>
<tr><td width=\"50%\" valign=\"top\">
<ol>
  <li> <a href=\"#Intro\" class=bold>Introduction</a>
";

system $^X, "$cmake_source_dir/doc/doxygen/scripts/intro2toc", "$cmake_source_dir/examples/$step/doc/intro.dox";

print "  <li> <a href=\"#CommProg\" class=bold>The commented program</a>\n";

system $^X, "$cmake_source_dir/doc/doxygen/scripts/program2toc", "$cmake_source_dir/examples/$step/$step.cc";

print
"</ol></td><td width=\"50%\" valign=\"top\"><ol>
  <li value=\"3\"> <a href=\"#Results\" class=bold>Results</a>
";

system $^X, "$cmake_source_dir/doc/doxygen/scripts/intro2toc", "$cmake_source_dir/examples/$step/doc/results.dox";

print
"  <li> <a href=\"#PlainProg\" class=bold>The plain program</a>
</ol> </td> </tr> </table>
\@endhtmlonly
";

system $^X, "$cmake_source_dir/doc/doxygen/scripts/create_anchors", "$cmake_source_dir/examples/$step/doc/intro.dox";

system $^X, "$cmake_source_dir/doc/doxygen/scripts/program2doxygen", "$cmake_source_dir/examples/$step/$step.cc";

system $^X, "$cmake_source_dir/doc/doxygen/scripts/create_anchors", "$cmake_source_dir/examples/$step/doc/results.dox";

print
"<a name=\"PlainProg\"></a>
<h1> The plain program</h1>
\@include \"$step.cc\"
 */
";
