######################################################################
# $Id$
######################################################################

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

# Find tooltip files
$list = `ls ../../../examples/step-*/doc/tooltip`;
@list = split "\n", $list;

foreach (@list)
{
    # Only the first line of the tooltip file is used
    open TF, "<$_";
    $tooltip = <TF>;
    close TF;
    chop $tooltip;
    
    m/step-(\d+)/;
    $n = $1;
    $toc =~ s/\@step$n\@/$tooltip/;
}

print $toc;
