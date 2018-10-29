## ---------------------------------------------------------------------
##
## Copyright (C) 2006 - 2015 by the deal.II authors
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

use strict;

my $tutorial_file = shift;
open TUTORIAL, "<$tutorial_file";

# Print the first part of tutorial.h.in up until the point where we
# find the line with '@@TUTORIAL_MAP@@'
while (my $line = <TUTORIAL>)
{
  last if($line =~ m/\@\@TUTORIAL_MAP\@\@/);
  print $line;
}

# List of additional node attributes to highlight purpose and state of the example
my %style = (
 "basic"          => ',height=.8,width=.8,shape="octagon",fillcolor="green"',
 "techniques"     => ',height=.35,width=.35,fillcolor="orange"',
 "fluids"         => ',height=.25,width=.25,fillcolor="yellow"',
 "solids"         => ',height=.25,width=.25,fillcolor="lightblue"',
 "time dependent" => ',height=.25,width=.25,fillcolor="blue"',
 "unfinished"     => ',height=.25,width=.25,style="dashed"',
 "code-gallery"   => ',height=.08,width=.125,shape="circle"',
    );

# Print a preamble setting common attributes
print << 'EOT'
digraph StepsMap
{
  overlap=false;
  edge [fontname="FreeSans",
        fontsize="10",
        labelfontname="FreeSans",
        labelfontsize="10",
        color="black",
        style="solid"];
  node [fontname="FreeSans",
        fontsize="10",
        shape="rectangle",
        height=0.2,
        width=0.4,
        color="black",
        fillcolor="white",
        style="filled"];
EOT
    ;

# Print all nodes of the graph by looping over the remaining
# command line arguments denoting the tutorial programs

my $step;
foreach $step (@ARGV)
{
    # read first line of tooltip file
    open TF, "$step/doc/tooltip"
        or die "Can't open tooltip file $step/doc/tooltip";
    my $tooltip = <TF>;
    close TF;
    chop $tooltip;

    # read first line of 'kind' file if it is a step;
    # otherwise assume it is a code gallery program. for
    # each of them, output something for 'dot' to generate
    # the dependencies graph from
    if ($step =~ /step-/
        &&
        !($step =~ /code-gallery/))
    {
      open KF, "$step/doc/kind"
          or die "Can't open kind file $step/doc/kind";
      my $kind = <KF>;
      chop $kind;
      close KF;

      die "Unknown kind '$kind' in file $step/doc/kind" if (! defined $style{$kind});

      my $number = $step;
      $number =~ s/^.*-//;

      printf "  Step$number [label=\"$number\", URL=\"\\ref step_$number\", tooltip=\"$tooltip\"";
      print "$style{$kind}";
    }
    else
    {
      # get at the name of the program; also create something
      # that can serve as a tag without using special characters
      my $name = $step;
      $name =~ s/^.*code-gallery\///;
      my $tag = $name;
      $tag =~ s/[^a-zA-Z]/_/g;

      printf "  code_gallery_$tag [label=\"\", URL=\"\\ref code_gallery_$tag\", tooltip=\"$tooltip\"";
      my $kind = "code-gallery";
      print "$style{$kind}";
    }

    print "];\n";
}

# Print all edges by going over the same list of tutorials again.
# Keep sorted by second node on edge!

my $step;
foreach $step (@ARGV)
{
    # read first line of dependency file
    open BF, "$step/doc/builds-on"
        or die "Can't open builds-on file $step/doc/builds-on";
    my $buildson = <BF>;
    close BF;
    chop $buildson;

    my $destination;
    if ($step =~ /step-/
        &&
        !($step =~ /code-gallery/))
    {
      my $number = $step;
      $number =~ s/^.*-//;
      $destination = "Step$number";
    }
    else
    {
      my $name = $step;
      $name =~ s/^.*code-gallery\///;
      my $tag = $name;
      $tag =~ s/[^a-zA-Z]/_/g;
      $destination = "code_gallery_$tag";
    }

    my $source;
    foreach $source (split ' ', $buildson) {
        $source =~ s/step-/Step/g;
        print "  $source -> $destination";
        if ($destination =~ /code_gallery/)
        {
            print " [style=\"dashed\", arrowhead=\"empty\"]";
        }
        print "\n";
    }
}

print "}\n";

# Copy that part of tutorial.h.in up until the point where we
# find the line with '@@TUTORIAL_LEGEND@@'
while (my $line = <TUTORIAL>)
{
  last if($line =~ m/\@\@TUTORIAL_LEGEND\@\@/);
  print $line;
}

# Print a preamble setting common attributes
print << 'EOT'
graph StepsDescription
{
  overlap=false;
  edge [fontname="FreeSans",
        fontsize="10",
        labelfontname="FreeSans",
        labelfontsize="10",
        color="black",
        style="solid"];
  node [fontname="FreeSans",
        fontsize="10",
        shape="rectangle",
        height=0.2,
        width=0.4,
        color="black",
        fillcolor="white",
        style="filled"];
EOT
    ;

my %kind_descriptions = (
 "basic"          => 'Basic techniques',
 "techniques"     => 'Advanced techniques',
 "fluids"         => 'Fluid dynamics',
 "solids"         => 'Solid mechanics',
 "time dependent" => 'Time dependent problems',
 "unfinished"     => 'Unfinished codes',
 "code-gallery"   => 'Code gallery',
    );

# for each kind, print a box in the same style as used in
# the connections graph; also print a fake box with a
# description of what each kind is. then connect these
my $kind;
foreach $kind (keys %style)
{
    my $escaped_kind = $kind;
    $escaped_kind =~ s/[^a-zA-Z]/_/g;
    printf "  $escaped_kind [label=\"\" $style{$kind}];\n";
    printf "  fake_$escaped_kind [label=\"$kind_descriptions{$kind}\", shape=plaintext];\n";
    printf "  $escaped_kind -- fake_$escaped_kind [style=dotted, arrowhead=odot, arrowsize=1];\n";
}
# now add connections to make sure they appear nicely next to each other
# in the legend
print "  basic -- techniques -- fluids -- solids -- time_dependent -- unfinished -- code_gallery;\n";

# we need to tell 'dot' that all of these are at the same
# rank to ensure they appear next to (as opposed to atop)
# each other
print "  {rank=same; basic, techniques, fluids, solids, time_dependent, unfinished, code_gallery}";

# end the graph
print "}\n";



# Then print the rest of tutorial.h.in
while (my $line = <TUTORIAL>)
{
  print $line;
}
close TUTORIAL;
