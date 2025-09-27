## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2009 - 2025 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

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
 "solids"         => 'coral',
 "time dependent" => 'darkolivegreen1',
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
  bgcolor=transparent;
  overlap=false;
  edge [fontname="FreeSans",
        fontsize="10",
        labelfontname="FreeSans",
        labelfontsize="10",
        color="cornflowerblue",
        style="solid"];
  node [fontname="FreeSans",
        fontsize="10",
        shape="rectangle",
        height=0.2,
        width=0.4,
        color="cornflowerblue",
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
    # we need to escape any double-quotes here before putting it in the dot file
    $tooltip =~ s/"/\\"/g;

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
      $tag =~ s/[^a-zA-Z_0-9]/_/g;

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
      $tag =~ s/[^a-zA-Z_0-9]/_/g;
      $destination = "code_gallery_$tag";
    }

    my $source;
    foreach $source (split ' ', $buildson) {
        $source =~ s/step-/Step/g;

        # Some tutorial programs have a large number of ancestors: this gets
        # particularly noticeable when we switch from sequential to distributed
        # programs. Moving step-40 down by adding some invisible nodes helps a
        # lot.
        if ($source eq "Step6" && $destination eq "Step40")
        {
            print "  Step40a [style=\"invis\"];";
            print "  Step40b [style=\"invis\"];";
            print "  Step40c [style=\"invis\"];";
            print "  Step40d [style=\"invis\"];";
            print "  Step40e [style=\"invis\"];";
            print "  Step40f [style=\"invis\"];";

            print "  Step6 -> Step40a [style=\"invis\"];";
            print "  Step40a -> Step40b [style=\"invis\"];";
            print "  Step40b -> Step40c [style=\"invis\"];";
            print "  Step40c -> Step40d [style=\"invis\"];";
            print "  Step40d -> Step40e [style=\"invis\"];";
            print "  Step40e -> Step40f [style=\"invis\"];";
            print "  Step40f -> Step40 [style=\"invis\"];";

            print "  Step6 -> Step40 [weight=100,color=\"$colors{$kind_map{$source}}\"];";
        }

        elsif ($source eq "Step6" && $destination eq "Step8")
        {
            print "  Step8a [style=\"invis\"];";

            print "  Step6 -> Step8a [style=\"invis\"];";
            print "  Step8a -> Step8 [style=\"invis\"];";

            print "  Step6 -> Step8";
        }
        elsif ($source eq "Step6" && $destination eq "Step15")
        {
            print "  Step15a [style=\"invis\"];";

            print "  Step6 -> Step15a [style=\"invis\"];";
            print "  Step15a -> Step15 [style=\"invis\"];";

            print "  Step6 -> Step15";
        }
        elsif ($source eq "Step16" && $destination eq "Step37")
        {
            print "  Step37a [style=\"invis\"];";
            print "  Step37b [style=\"invis\"];";
            print "  Step37c [style=\"invis\"];";
            print "  Step37d [style=\"invis\"];";

            print "  Step16 -> Step37a [style=\"invis\"];";
            print "  Step37a -> Step37b [style=\"invis\"];";
            print "  Step37b -> Step37c [style=\"invis\"];";
            print "  Step37c -> Step37d [style=\"invis\"];";
            print "  Step37d -> Step37 [style=\"invis\"];";

            print "  Step16 -> Step37 [weight=100,color=\"$colors{$kind_map{$source}}\"];";
        }
        # This one looks better if it is down a level due to the other
        # artificial nodes.
        elsif ($source eq "Step21" && $destination eq "Step22")
        {
            print "  Step22a [style=\"invis\"];";
            print "  Step22b [style=\"invis\"];";

            print "  Step21 -> Step22a [style=\"invis\"];";
            print "  Step22a -> Step22b [style=\"invis\"];";
            print "  Step22b -> Step22 [style=\"invis\"];";

            print "  Step21 -> Step22 [weight=100,color=\"$colors{$kind_map{$source}}\"];";
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
            # to try and keep programs of the same kind together.
            #
            # There are two exceptions:
            # * The "basic" tutorial programs: these are
            #   going to be connected by edges of weight 100, ensuring
            #   that they are all essentially aligned vertically.
            # * Code gallery programs: Here, we don't care much where
            #   they are placed, and so don't treat edges between these kinds
            #   of programs as special
            if ($kind_map{$source} eq $kind_map{$destination}
                &&
                ! ($kind_map{$source} eq "code-gallery"))
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

            # If the destination is a code gallery program, use a dashed line
            if ($kind_map{$destination} eq "code-gallery")
            {
                $edge_attributes .= "style=\"dashed\", arrowhead=\"empty\", color=\"gray\",";
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
  bgcolor=transparent;
  overlap=false;
  edge [fontname="FreeSans",
        fontsize="12",
        labelfontname="FreeSans",
        labelfontsize="10",
        color="cornflowerblue",
        style="solid"];
  node [fontname="FreeSans",
        fontsize="11",
        shape="rectangle",
        height=0.2,
        width=0.4,
        color="cornflowerblue"];
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
    printf "  $escaped_kind [label=\"\" $style{$kind}, style=\"filled\" fillcolor=\"$colors{$kind}\"];\n";
    printf "  fake_$escaped_kind [label=\"$kind_descriptions{$kind}\", fontcolor=\"cornflowerblue\" fontsize=12 shape=plaintext];\n";
    printf "  $escaped_kind -- fake_$escaped_kind [style=\"bold,dotted\", arrowhead=odot, arrowsize=1];\n";
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
