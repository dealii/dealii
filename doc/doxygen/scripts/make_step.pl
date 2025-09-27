## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2013 - 2024 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

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
<table class=\"tutorial\" width=\"75%\">
<tr><th colspan=\"2\"><b><small>Table of contents</small></b></th></tr>
<tr><td width=\"50%\" valign=\"top\">
<ol>
  <li> <a href=\"#$step_underscore-Intro\" class=bold>Introduction</a>
";

system $^X, "$cmake_source_dir/doc/doxygen/scripts/intro2toc.pl", "$cmake_source_dir/examples/$step/doc/intro.dox" , "--prefix=$step_underscore";

print "  <li> <a href=\"#$step_underscore-CommProg\" class=bold>The commented program</a>\n";

system $^X, "$cmake_source_dir/doc/doxygen/scripts/program2toc.pl", "$cmake_source_dir/examples/$step/$step.cc" , "--prefix=$step_underscore";

print
"</ol></td><td width=\"50%\" valign=\"top\"><ol>
  <li value=\"3\"> <a href=\"#$step_underscore-Results\" class=bold>Results</a>
";

system $^X, "$cmake_source_dir/doc/doxygen/scripts/intro2toc.pl", "$cmake_source_dir/examples/$step/doc/results.dox" , "--prefix=$step_underscore";

print
"  <li> <a href=\"#$step_underscore-PlainProg\" class=bold>The plain program</a>
</ol> </td> </tr> </table>
\@endhtmlonly

";

system $^X, "$cmake_source_dir/doc/doxygen/scripts/create_anchors.pl", "$cmake_source_dir/examples/$step/doc/intro.dox" , "--prefix=$step_underscore";


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
print " * <a name=\"$step_underscore-CommProg\"></a>\n";
print " * <h1> The commented program</h1>\n";

system $^X, "$cmake_source_dir/doc/doxygen/scripts/program2doxygen.pl", "$cmake_source_dir/examples/$step/$step.cc" , "--prefix=$step_underscore";

system $^X, "$cmake_source_dir/doc/doxygen/scripts/create_anchors.pl", "$cmake_source_dir/examples/$step/doc/results.dox" , "--prefix=$step_underscore";


# Move to the stripped, plain program. The same principle as above
# applies for newlines.
print " *\n";
print " *\n";
print
"<a name=\"$step_underscore-PlainProg\"></a>
<h1> The plain program</h1>
\@include \"$step.cc\"
*/
";
