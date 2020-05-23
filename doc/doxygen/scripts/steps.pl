## ---------------------------------------------------------------------
##
## Copyright (C) 2006 - 2020 by the deal.II authors
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

# List of additional node and edge attributes to highlight purpose and state of
# a tutorial or code gallery program. For a list of colors, take a look here:
#   https://www.graphviz.org/doc/info/colors.html
my %colors = (
 "basic"          => 'green',
 "techniques"     => 'orange',
 "fluids"         => 'yellow2',
 "solids"         => 'lightblue',
 "time dependent" => 'dodgerblue1',
 "unfinished"     => 'white',
 "code-gallery"   => 'white',
    );

my %style = (
 "basic"          => ',height=.8,width=.8,shape="octagon"',
 "techniques"     => ',height=.35,width=.35',
 "fluids"         => ',height=.25,width=.25',
 "solids"         => ',height=.25,width=.25',
 "time dependent" => ',height=.25,width=.25',
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
my %kind_map;
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

      $kind_map{"Step$number"} = $kind;

      printf "  Step$number [label=\"$number\", URL=\"\\ref step_$number\", tooltip=\"$tooltip\"";
      print "$style{$kind},fillcolor=\"$colors{$kind}\"";
    }
    else
    {
      # get at the name of the program; also create something
      # that can serve as a tag without using special characters
      my $name = $step;
      $name =~ s/^.*code-gallery\///;
      my $tag = $name;
      $tag =~ s/[^a-zA-Z]/_/g;

      $kind_map{"code_gallery_$tag"} = "code-gallery";

      printf "  code_gallery_$tag [label=\"\", URL=\"\\ref code_gallery_$tag\", tooltip=\"$tooltip\"";
      my $kind = "code-gallery";
      print "$style{$kind},fillcolor=\"$colors{$kind}\"";
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

        # We want to treat the step-6 -> step-40 edge differently. If
        # it is printed like any other edge, i.e., with step-40
        # directly below step-6, then we end up with a tangle of lines
        # because there are so many sub-graphs that originate from
        # step-6 (with the next node at the same level as step-40) but
        # where step-40 then feeds with a long line into one of the
        # programs further down. It looks better if we place step-40 a
        # couple of levels further down in the graph so that the
        # higher up parts of the graph consists of the sequential
        # programs and the lower-down parts to the parallel ones.
        #
        # The way to do this is to insert a few invisible nodes (with
        # invisible edges) between step-6 and step-40.
        if ($source eq "Step6" && $destination eq "Step40")
        {
            print "  Step40a [style=\"invis\"];";
            print "  Step40b [style=\"invis\"];";
            print "  Step40c [style=\"invis\"];";

            print "  Step6 -> Step40a [style=\"invis\"];";
            print "  Step40a -> Step40b [style=\"invis\"];";
            print "  Step40b -> Step40c [style=\"invis\"];";
            print "  Step40c -> Step40 [style=\"invis\"];";

            print "  Step6 -> Step40 [weight=100,color=\"$colors{$kind_map{$source}}\"];";
        }

        # All other edges in the graph
        else
        {
            print "  $source -> $destination";

            my $edge_attributes = "";

            # Determine the style of the arrow that connects
            # the two nodes. If the two nodes are of the same
            # kind, use the same color as the nodes as this makes
            # reading the flow of the graph a bit easier. Furthermore,
            # set the edge weight to 5 (instead of the default of 1)
            # to try and keep programs of the same kind together. The
            # exception is the "basic" tutorial programs: these are
            # going to be connected by edges of weight 100, ensuring
            # that they are all essentially aligned vertically.
            if ($kind_map{$source} eq $kind_map{$destination})
            {
                $edge_attributes = "color=\"$colors{$kind_map{$source}}\",";
                if ($kind_map{$source} eq "basic")
                {
                    $edge_attributes .= "weight=100,";
                }
                else
                {
                    $edge_attributes .= "weight=5,";
                }
            }

            # If the destination is a code gallery program, used a dashed line
            if ($destination =~ /code_gallery/)
            {
                $edge_attributes .= "style=\"dashed\", arrowhead=\"empty\",";
            }
            print " [$edge_attributes];\n";
        }
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
    printf "  $escaped_kind [label=\"\" $style{$kind}, fillcolor=\"$colors{$kind}\"];\n";
    printf "  fake_$escaped_kind [label=\"$kind_descriptions{$kind}\", shape=plaintext];\n";
    printf "  $escaped_kind -- fake_$escaped_kind [style=dotted, arrowhead=odot, arrowsize=1];\n";
}
# now add connections to make sure they appear nicely next to each other
# in the legend
print "  basic -- techniques [style=invis];\n";
print "  techniques -- fluids [style=invis];\n";
print "  fluids -- solids [style=invis];\n";
print "  solids -- time_dependent [style=invis];\n";
print "  time_dependent -- unfinished [style=invis];\n";
print "  unfinished -- code_gallery [style=invis];\n";

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
