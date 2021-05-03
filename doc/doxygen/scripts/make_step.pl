## ---------------------------------------------------------------------
##
## Copyright (C) 2013 - 2020 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
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
\@page $step_underscore The $step tutorial program
";

open BF, "$cmake_source_dir/examples/$step/doc/builds-on"
    or die "Can't open builds-on file $cmake_source_dir/examples/$step/doc/builds-on";
my $buildson = <BF>;
close BF;
chop $buildson;

# At the very top, print which other programs this one builds on. The
# filter script will replace occurrences of step-XX by the appropriate
# links.
if ($buildson ne "")
{
    $buildson =~ s/ /, /g;
    print "This tutorial depends on $buildson.\n\n";
}

# then show the table of contents
print
"\@htmlonly
<table class=\"tutorial\" width=\"50%\">
<tr><th colspan=\"2\"><b><small>Table of contents</small></b></th></tr>
<tr><td width=\"50%\" valign=\"top\">
<ol>
  <li> <a href=\"#Intro\" class=bold>Introduction</a>
";

system $^X, "$cmake_source_dir/doc/doxygen/scripts/intro2toc", "$cmake_source_dir/examples/$step/doc/intro.dox";

print "  <li> <a href=\"#CommProg\" class=bold>The commented program</a>\n";

my $file_extension;

if (-f "$cmake_source_dir/examples/$step/$step.cc")
{
  $file_extension = cc;
}

if (-f "$cmake_source_dir/examples/$step/$step.cu")
{
  $file_extension = cu;
}

system $^X, "$cmake_source_dir/doc/doxygen/scripts/program2toc", "$cmake_source_dir/examples/$step/$step.$file_extension";

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


# Start the commented program by writing two empty lines. We have had
# cases where the end of the intro.dox was missing a newline, and in
# that case doxygen might get confused about what is being added here
# to the end of an existing line. So add a newline.
#
# But then we also had a situation where doxygen was confused about a
# line starting with an anchor (see #9357). It's not clear what the
# cause is, but making sure that there is an empty line in between
# solves the problem -- so a second newline character.
print " *\n";
print " *\n";
print " * <a name=\"CommProg\"></a>\n";
print " * <h1> The commented program</h1>\n";

system $^X, "$cmake_source_dir/doc/doxygen/scripts/program2doxygen", "$cmake_source_dir/examples/$step/$step.$file_extension";

system $^X, "$cmake_source_dir/doc/doxygen/scripts/create_anchors", "$cmake_source_dir/examples/$step/doc/results.dox";


# Move to the stripped, plain program. The same principle as above
# applies for newlines.
print " *\n";
print " *\n";
print
"<a name=\"PlainProg\"></a>
<h1> The plain program</h1>
\@include \"$step.$file_extension\"
*/
";
