#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
package kdocHTML;

## KDOC HTML output
#
# From main:
# children classRef classSource dontDoPrivate 
# generatedByText headerList quiet

require Ast;

BEGIN 
{
	$lib = "";
	$outputdir = ".";
	$linkRef = 0;
	$modname = "kdocHTML";
	%classList = ();
}


sub dumpDoc
{
	$genText = $main::generatedByText;

	($lib, $root, $outputdir) = @_;

	print "Generating HTML documentation.\n" unless $main::quiet;
	writeList();
	writeHier();
	writeHeaderList();

	print "Generating HTMLized headers.\n" unless $main::quiet;
	foreach $header ( @main::headerList ) {
		markupCxxHeader( $header );
	}
}

sub writeHeaderList
{

	open(HDRIDX, ">$outputdir/header-list.html") 
	|| die "Couldn't create $outputdir/header-list.html\n";

	print HDRIDX<<EOF;
<HTML><HEAD><TITLE>$lib Header Index</TITLE></HEAD><BODY bgcolor="#ffffff">
	<H1>$lib Header Index</H1>
<p>
[<A HREF="index.html">Index</A>] [<A HREF="hier.html">Hierarchy</A>] [Headers]
</p>
<HR>
<UL>
EOF

	foreach $header ( sort @main::headerList ) {
		$_ = $header;
# convert dashes to double dash, convert path to dash
		s/-/--g/g;
		s/\/|\./-/g;
		
		# never let a file name start with a '-' sign, since this makes
		# some programs think the file name really was a parameter to the
		# program. Filenames starting with '-' usually happen if you let
		# kdoc work on files like '../../include/foo/bar.h', given the
		# above substitution rules. Therefore: if a filename starts with
		# '-', prepend a '_'
		s/^-/_-/;

		print HDRIDX  "\t<LI><A HREF=\"$_.html\">$header</A></LI>\n";
	}
print HDRIDX "</UL>\n<HR>\n<address>$genText",
			"</address>\n</BODY>\n</HTML>\n";
}

sub escape
{
	my( $str ) = @_;

	$str =~ s/&/\&amp;/g;
	$str =~ s/</\&lt;/g;
	$str =~ s/>/\&gt;/g;
	$str =~ s/"/\&quot;/g;

	return $str;
}

# In file names, we don't want to have the special characters be
# replaced by the rules '>' -> '&gt;' etc., so we revert the
# replacement. The reason is that when giving Netscape the tag
# <a href="abc&lt;2&gt;.html"> xyz </a>, it looks for the file
# "abc<2>.html". Strange, eh? We'd like to replace filenames
# back to the original names, which is not a good idea also, since
# filenames containing <>&" are no good. Therefore, we stay with
# the HTML tags, but replace the '&' and ';' by ',' which should
# not be a character which appears in class and file names, but
# which has no special meaning for the shell either.
sub filenameUnescape
{
	my( $str ) = @_;

	$str =~ s/\&amp;/,amp,/g;
	$str =~ s/\&lt;/,lt,/g;
	$str =~ s/\&gt;/,gt,/g;
	$str =~ s/\&quot;/,quot,/g;

	return $str;
}

sub writeList()
{

	$root->Visit( $modname );

	foreach $class ( @{$classes} )
	{
		$classList{ $class->{'astNodeName'} } = $class;
	}

	open ( LIBIDX, ">".$outputdir."/index.html" )
		|| die "Couldn't write to ".$outputdir."/index.html" ;

	print LIBIDX<<EOF;
<HTML><HEAD><TITLE>$lib Class Index</TITLE></HEAD><BODY bgcolor="#ffffff">
	<H1>$lib Class Index</H1>
<p>
[Index] [<A HREF="hier.html">Hierarchy</A>] [<A HREF="header-list.html">Headers</A>]
</p>
<HR>
<TABLE BORDER="0" WIDTH="100%">
EOF

	foreach $className ( sort keys %classList ) {
		$class = $classList { $className };

		$class->Visit( $modname );

		$className = $astNodeName;

		# escape special characters in class names (e.g. <>&)
		# however, we need to unescape them when using the
		# class name as a file name
		$astNodeName = escape( $astNodeName );

		print LIBIDX "<TR><TD ALIGN=\"LEFT\" VALIGN=\"TOP\">",
		"<a href=\"", filenameUnescape($astNodeName), ".html\">",
		$astNodeName, "</a></TD>\n", "\t<TD>", 
		makeReferencedText( $ClassShort ), "</TD></TR>\n";	

		writeClass();

		Ast::UnVisit();
		
	}

	Ast::UnVisit();

print LIBIDX<<EOF;
</TABLE>
<HR>
<TABLE WIDTH="100%" BORDER="0">
<TR>
<TD ALIGN="LEFT" VALIGN="TOP"><address>$genText</address></TD>
<TD ALIGN="RIGHT" VALIGN="TOP">
<b>K</b><i>doc</i>
</TD>
</TR>
</TABLE>
</BODY></HTML>
EOF

}

###########################################################
#
# Function
#  makeReferencedText( text )
#
###########################################################

sub makeReferencedText
{
    local( $t ) = @_;

    while ( $t =~ /(\@ref\s+)([\w_\#~\d\:]+\W)/ )
    {##
	$pat = $1.$2;
	$ref = $2;
	chop $ref;
	$repl = findClassReference( $ref ) . " ";
	$pat =~ s/\(/\\(/g;
	$pat =~ s/\)/\\)/g;
	$t =~ s/$pat/$repl/;
    }

    $t =~ s/</&lt;/g;
    $t =~ s/>/&gt;/g;
    $t =~ s/&lt;\/p&gt;&lt;p&gt;/<\/p><p>/g;
    $t =~ s/\#([^\#]*)\#/<CODE>$1<\/CODE>/g;
    $t =~ s/\\begin\{verbatim\}/<pre>/g;
    $t =~ s/\\end\{verbatim\}/<\/pre>/g;
    $t =~ s/\\begin\{itemize\}/<UL>/g;
    $t =~ s/\\end\{itemize\}/<\/UL>/g;
    $t =~ s/\\item/<LI>/g;
    $t =~ s/\{\\bf([^}]*)}/<b>$1<\/b>/g;
    $t =~ s/\\section\{([^}]*)\}/<H3>$1<\/H3>/g;
    $t =~ s/\\subsection\{([^}]*)\}/<H4>$1<\/H4>/g;
    $t =~ s/\\subsubsection\{([^}]*)\}/<H5>$1<\/H5>/g;

#while($t =~ m/(\#([^\#]*)\#)/)
#{
#    $full = $1;
#    $text = $2;
#    $text =~ s/</&lt;/g;
#    $text =~ s/>/&gt;/g;
#    $text =~ s/\&/&amp;/g;
#    $t =~ s/$full/<CODE>$text<\/CODE>/;
#}
    return $t;
}

###########################################################
#
# Function
#  findClassReference( class )
#
###########################################################

sub findClassReference
{
    my( $c ) = @_;
    my( $old );


    if ( $c eq "" ) { return ""; }
    if ( $c eq "#" ) { return "#"; }

	$c =~ s/\#/::/g;

    $old = $c;

    if ( $c =~ /^\s*(\#|::)(.*)$/ )
    {
	$c = "$className\:\:$2";
	$pref = $1;
	$old =~ s/$pref//;
    }

    if ( !defined $main::classRef{ $c } )
    {
	if ( !defined $main::noErr{ $c } )
	{
#	    print main::ERROUT "$filename: '$c' undefined Reference\n\n"; 
	}
	return $c;
    }
    else
    {
	if ( $c =~ /\:\:/ )
	{
	    return "<a href=\"".escape($main::classRef{ $c })."\">".
	    	escape($old)."</a>";
	}
	else
	{
	    return "<a href=\"".escape($main::classRef{ $c })."\">".
	    	escape($c)."</a>";
	}
    }
}

#############
# 
# writeHier -- Write HTML Class Hierarchy
#
#############

sub writeHier
{

open( HIER, ">".$outputdir."/hier.html" ) 
	||die "Cannot write to $outputdir/hier.html";

print HIER <<EOF;
<HTML><HEAD><TITLE>$lib Class Hierarchy</TITLE></HEAD>
<BODY BGCOLOR="#ffffff">
<H1>$lib Class Hierarchy</H1>
<p>
[<A HREF="index.html">Index</A>] [Hierarchy] [<A HREF="header-list.html">Headers</A>] 
</p>
<HR>
EOF

printHierarchy( "<>" );

print HIER<<EOF;
<HR>
<TABLE WIDTH="100%" BORDER="0">
<TR>
<TD ALIGN="LEFT" VALIGN="TOP"><address>$genText</address></TD>
<TD ALIGN="RIGHT" VALIGN="TOP">
<b>K</b><i>doc</i>
</TD>
</TR>
</TABLE>
</BODY></HTML>
EOF

### here we are!!
    # build a list of parents and children for each class
    # use it to generate a tree of classes suitable for
    # a tree view of the inheritance lattice

    # we may not use the %classList for the children information
    # since some classes are parents but not in %classList (those
    # classes that are found as base classes but are not documented
    # within this module
#    my (%kids);
#    my (%parents);
#    my (%classLevel);
#    my ($classNames);
#
#    foreach $className (keys %main::children)
#    {
#	$kids{$className} = $main::children{$className};
#    }
#
#    foreach $className (keys %kids)
#    {
#	foreach $c (split /;/, $main::children{$className})
#	{
#	    $parents{$c} = $parents{$c} . $className . ";"
#		if !($className eq "<>");
#	}
#    }
#
#    # build a list of all class names, irrespective of whether
#    # they are documented or not.
#    foreach $class (keys %kids, keys %parents, keys %classList)
#    {
#	$classNames = $classNames . $class . ";"
#	    if (! ($classNames =~ /$class;/) &&
#		! ($class eq "<>"));
#    }
#
#
#    # now evaluate the list to get the proper hierarchy level
#    # for each class. Loop until all classes are processed,
#    # which may take some iterations since a class can only
#    # be processed if all of its parents have been processed
#    my ($unprocessedClasses);
#    my ($classLevel);
#    my ($maxParentLevel);
#    $unprocessedClasses = 1;
#    while ($unprocessedClasses == 1)
#    {
#	$unprocessedClasses = 0;
#	foreach $class (split /;/, $classNames)
#	{
#	    if (! defined $classLevel{$class})
#	    {
#		# if class has no children or no parents: assign
#		# level zero. for classes with no parents, this
#		# decision may later be revised
#		if ((! $parents{$class}) ||
#		    (! defined $parents{$class}))
#		{
#		    $classLevel{$class} = 0;
#		}
#		else
#		{
#		    foreach $parent (split /;/, $parents{$class})
#		    {
#			# find the maximum level any of the
#			# children has
#			$maxParentLevel = 0;
#
#			# if the parent has not been processed,
#			# then don't do so for this class also
#			if (! defined $classLevel{$parent})
#			{
#			    $unprocessedClasses = 1;
#			    last;
#			}
#			else
#			{
#			    if ($classLevel{$parent} > $maxParentLevel)
#			    {
#				$maxParentLevel = $classLevel{$parent};
#			    }
#			}
#			$classLevel{$class} = $maxParentLevel+1;
#		    }
#		}
#	    }
#	}
#    }
#
#
#    foreach $class (sort split /;/, $classNames )
#    {
#	print "--> ", $class, ": $classLevel{$class} <c=", $kids{$class}, "><p=",
#	$parents{$class}, ">\n";
#
#    }
}


######
#
# Dump the hierarchy for a particular class
#
######

sub printHierarchy
{
	my( $className ) = @_;

	print HIER "\t<LI>", findClassReference($className),"\n" 
		if $className ne "<>";

	if( defined $main::classSource{ $className } ) {
		print HIER "<small> (",$main::classSource{$className},")</small>\n";
	}

	if ( defined $main::children{$className} ) {
	
		print HIER "\t<UL>\n";

	        foreach $kid ( sort split /;/, $main::children{$className} )
	        {
	                printHierarchy( $kid );
		}
		print HIER "\t</UL>\n";
       	 }

	 print HIER "\t</LI>" if $className ne "<>";
}

# #############
#
# Write out the documentation for a single class.
#
# #############

sub writeClass()
{
	open( CLASS, ">" .
              $outputdir."/".filenameUnescape($astNodeName) .
              ".html"
            ) || die "Couldn't write to file " .
                     $outputdir."/".filenameUnescape($astNodeName).".html";

	print CLASS<<EOF;
<HTML><HEAD><TITLE>$astNodeName Class</TITLE></HEAD>
<BODY bgcolor="#ffffff\">\n
<H1>$astNodeName Class Reference</H1>
<p>
[<A HREF="index.html">$lib Index</A>] [<A HREF="hier.html">$lib Hierarchy</A>]
[<A HREF="header-list.html">Headers</A>]
</p>
<HR>
EOF
	
	print CLASS "<h3>This class is Deprecated.</h3>\n" if $Deprecated;
	print CLASS "<h3>Internal class - general use not recommended.</h3>\n" 
		if $Internal;

	print CLASS "<P>", makeReferencedText( $ClassShort ),
    			"  <a href=\"#short\">More...</a></P>\n" 
			if $ClassShort ne "";

	$_ = $Header;
	s/-/--/g;
	s/\/|\./-/g;

        # never let a file name start with a '-' sign, since this makes
        # some programs think the file name really was a parameter to the
        # program. Filenames starting with '-' usually happen if you let
        # kdoc work on files like '../../include/foo/bar.h', given the
        # above substitution rules. Therefore: if a filename starts with
        # '-', prepend a '_'
        s/^-/_-/;
    
	$escapedpath = $_;

# include file

	print CLASS<<EOF;
<P>
<code>
	#include &lt;<a href=\"$escapedpath.html\">$Header</a>&gt;
</code>\n
</P>
EOF


# Template

	if( $class->{ "TmplArgs" } ne "" ) {

		print CLASS "\n<p>Template Form: <code>",
		"template &lt; ", escape( $class->{ "TmplArgs" }), 
		" &gt; $astNodeName</code> </p>\n";
	}

# Inheritance
	my ( $pcount ) = 0;

	if( defined $class->{ "Ancestors" } ) {
		print CLASS "\n<P>\nInherits: ";
		
		foreach $foreparent ( @{$Ancestors} ) {
			print CLASS ", " if $pcount != 0;
			$pcount = 1;

			print CLASS findClassReference( 
				$foreparent->{"astNodeName"} );

			if( defined $main::classSource{$foreparent->{
					"astNodeName"}} ) {
				print CLASS " (",$main::classSource{
					$foreparent->{"astNodeName"}},")";
			}
		}

		print CLASS "\n<\P>";
	}

# Now list members.

	$linkRef = 0;

	listMethods( "Public Members", $public );
	listMethods( "Public Slots", $public_slots );
	listMethods( "Protected Members", $protected );
	listMethods( "Protected Slots", $protected_slots );
	listMethods( "Signals", $signals );

	if( $main::dontDoPrivate == 0 ) {
		listMethods( "Private Members", $private );
		listMethods( "Private Slots", $private_slots );
	}
print CLASS<<EOF;
<HR>
EOF

# Dump Description.
	Ast::UnVisit();
	$class->Visit($modname);
	writeClassDescription();

# Document members.
	
	$linkRef = 0;

	writeMemberDoc( "public",	$public );
	writeMemberDoc( "public slot",	$public_slots );
	writeMemberDoc( "protected",	$protected );
	writeMemberDoc( "protected slot", $protected_slots );
	writeMemberDoc( "signal", $signals );

	if( $main::dontDoPrivate == 0 ) {
		writeMemberDoc( "private", $private );
		writeMemberDoc( "private slot", $private_slots );
	}

# File Generation info
	print CLASS "<HR>\n<TABLE WIDTH=\"100%\">",
			"<TR><TD ALIGN=\"left\" VALIGN=\"top\">\n";

	if( $Author ne "" || $Version ne "" ) {
		print CLASS "\n<UL>";

		print CLASS "<LI><I>Author</I>: ", escape($Author), "</LI>\n"
			if $Author ne "";

		print CLASS "<LI><I>Version</I>: ", escape($Version), "</LI>\n"
			if $Version ne "";
		print CLASS "<LI>$genText</LI>\n";

		print CLASS "</UL>";
	} else {
		print CLASS "<address>$genText<address>\n";
	}

            print CLASS<<EOF;
</TD><TD ALIGN="RIGHT" VALIGN="TOP">
<b>K</b><i>doc</i>
</TD>
</TR></TABLE></BODY></HTML>\n
EOF

}

######
# Lists all methods of a particular visibility.
#####

sub listMethods()
{
	my( $desc, $members ) = @_;
	
	return if !defined $members || $#{$members} == -1;

print CLASS<<EOF;

<H2>$desc</H2>
<UL>
EOF

	foreach $member ( @{$members} ) {
		$member->Visit($modname);

		print CLASS "<LI>";

		if( $Description eq "" && $See eq "" ) {
			$link = "name=\"";
		} else {
			$link = "href=\"\#";
		}

		if( $Keyword eq "property" ) {
			print CLASS escape($Type), 
				"<b><a ",$link,"ref$linkRef\">",
				escape($astNodeName), "</a></b>",
			        escape($Array),
			        "\n";
		}
		elsif( $Keyword eq "method" ) {
			print CLASS escape($ReturnType),
				" <b><a ", $link, "ref$linkRef\">",
				escape($astNodeName), "</a></b> (",
				escape($Parameters),
				") ",$Const,"\n";
		}
		elsif( $Keyword eq "enum" ) {
			print CLASS "enum <b><a ", $link, "ref$linkRef\">",
				escape($astNodeName),"</a></b> {",
				escape($Constants),"}\n";
		}
		elsif( $Keyword eq "typedef" ) {
			print CLASS "typedef ", escape($Type), " <b><a ", $link, 
				"ref$linkRef\">", escape($astNodeName),
			        "</a>", escape($Array), "</b>";
		}

		print CLASS "</LI>\n";

		$linkRef += 1;

		Ast::UnVisit();
	}

print CLASS<<EOF;
</UL>
EOF

}

sub writeClassDescription
{
	my( $term ) = "";
	my( $once ) = 0;

	if( $Description ne "" ) {
		print CLASS "<H2><a name=\"short\">Detailed Description</a>",
			"</H2>\n<P>\n";

		$term = "<HR>\n";

		$Description =~ s/\n(\s*\n)+/\n<\/p><p>\n/g;

		print CLASS "\n", makeReferencedText($Description) ,"\n</P>";
	}
	
	if( $ClassSee ne "" ) {
		print CLASS "<p><b>See Also</b>: ";

		$once = 0;
		$term = "<HR>" if $term eq "";

		foreach $item ( split /[\s,]+/, $ClassSee ) {

			if( $once ) { print CLASS ", "; } 
			else { $once = 1; }

			print CLASS findClassReference( $item ); 
		}

		print CLASS "</p>\n";
	}

	print CLASS $term if $term ne "";
}

sub writeMemberDoc
{
	my( $visibility, $node ) = @_;
	my( $once ) = 0;
	my( $myvis ) = "";

	return if !defined $node || $#{$node} == -1 ;

	foreach $member ( @{$node} ) {
		$member->Visit($modname);

		if( $Description eq "" && $See eq "" ) {
			$linkRef++;
			Ast::UnVisit();
			next;
		}
			
		if( $Keyword eq "property" ) {
			print CLASS "<H4><b>",refString(escape($Type)),
			"<a name=\"ref",$linkRef,"\"></a>",
			"<a name=\"",$astNodeName,"\">",
			escape($astNodeName), "</a>",
			refString(escape($Array)),
			"   </b><font color=gray><code>[", $visibility, "]</code>",
			"</font></H4>\n";	
		}
		elsif( $Keyword eq "method" ) {
			$myvis = $visibility;

			if ( $ReturnType =~ /virtual/ ) {
				$myvis .= " virtual";
				$ReturnType =~ s/\s*virtual\s*//g;
			}

			if ( $ReturnType =~ /static/ ) {
				$myvis .= " static";
				$ReturnType =~ s/\s*static\s*//g;
			}

			print CLASS "<H4><b>",
			refString(escape($ReturnType)), 
			" <a name=\"ref",$linkRef,"\"></a>", 
			"<a name=\"",$astNodeName,"\">", 
			escape($astNodeName), "</a>(",
			refString(escape($Parameters)),
			") $Const </b><font color=gray>",
			"<code>[", $myvis, "]</code></font></H4>\n";
		}
		elsif( $Keyword eq "enum" ) {
			print CLASS "<H4><b>enum <a name=\"ref", 
			$linkRef,"\"></a>",
			"<a name=\"",$astNodeName, "\"></a>",
			escape($astNodeName), " (",
			$Constants,") </b><font color=gray>",
			"<code>[", $visibility, "]</code></font></H4>\n";
		}
		elsif( $Keyword eq "typedef" ) {
			print CLASS "<H4><b>typedef ",refString(escape($Type)),
			" <a name=\"ref", $linkRef,"\"></a>",
			escape($astNodeName),
			refString(escape($Array)),
			"  </b><font color=gray>",
			"<code>[", $visibility, "]</code></font></H4>";
		}

		print CLASS "<p><strong>Deprecated member.</strong></p>\n" 
			if $MethDeprecated;
		print CLASS "<p><strong>For internal use only.</strong></p>\n" 
			if $MethInternal;

		$Description =~ s/\n(\s*\n)+/\n<\/p><p>\n/g;

		print CLASS "<p>",makeReferencedText($Description),"</p>\n";
		

		if( $Keyword eq "method" ) {

			if( $#{$ParamDoc} != -1 ) {
				print CLASS "<dl><dt><b>Parameters</b>:<dd>\n",
					"<table width=\"100%\" border=\"0\">\n";

				foreach $parameter ( @{$ParamDoc} ) {
					print CLASS 
					"<tr><td align=\"left\" valign=\"top\">\n",
					$parameter->{"astNodeName"};

					print CLASS
					"</td><td align=\"left\" valign=\"top\">\n",
					makeReferencedText(
					$parameter->{"Description"}),
					"</td></tr>\n";
				}

				print CLASS "</table>\n</dl>\n";

			}

			print CLASS "<dl><dt><b>Returns</b>:<dd>\n",
				makeReferencedText( $Returns ),
				"</dl>\n" if $Returns ne "";

			print CLASS "<dl><dt><b>Throws</b>:<dd>\n",
				makeReferencedText( $Exceptions ),
				"</dl>\n" if $Exceptions ne "";
		}

	if( $See ne "" ) {
		print CLASS "<dl><dt><b>See Also</b>:<dd>";

		$once = 0;

		foreach $item ( split /[\s,]+/, $See ) {

			if( $once ) { print CLASS ", "; } 
			else { $once = 1; }

			print CLASS findClassReference( $item ); 
		}

		print CLASS "</dl>\n";
	}

	print CLASS $term if $term ne "";

		Ast::UnVisit();

		$linkRef++;
	}
}

#
# Check and markup words that are in the symbol table.
# Use sparingly, since EVERY WORD is looked up in the table.
#

sub refString
{
	my( $source ) = @_;

	foreach $werd ( split /[^\w><]+/, $source ) {

		$ref = findClassReference( $werd );
		if( $ref ne $werd ) {
			$source =~ s/$werd/$ref/g;
		}
	}

	return $source;
}

######################################################################
# markupCxxHeader -- Converts the header into a semi-fancy HTML doc
######################################################################

$outputFilename = "";

# tokens
#$nonIdentifier="[^\w]+";
#$keywords="if|union|class|struct|operator|for|while|do|goto|new|delete|friend|typedef|return|switch|const|inline|asm|sizeof|virtual|static";
#$types="void|int|long|float|double|unsigned|char|bool|short";

# styles
#$keywordStyle="strong";
#$typeStyle="u";
#$stringStyle="i";

# Takes the path to a C++ source file,
# spits out a marked-up and cross-referenced HTML file.

sub markupCxxHeader
{
	my( $filename ) = @_;
	$className = "";
	my( $reference );
	my( @inheritance );
	my( $word );

	open( HFILE, $filename ) || die "Couldn't open $filename to read.\n";

	$_ = $filename;
	# convert dashes to double dash, convert path to dash
	s/-/--g/g;
	s/\/|\./-/g;

        # never let a file name start with a '-' sign, since this makes
        # some programs think the file name really was a parameter to the
        # program. Filenames starting with '-' usually happen if you let
        # kdoc work on files like '../../include/foo/bar.h', given the
        # above substitution rules. Therefore: if a filename starts with
        # '-', prepend a '_'
        s/^-/_-/;

	$outputFilename = $_;
	$outputFilename = $outputdir."/".$outputFilename;

	open( HTMLFILE, ">$outputFilename.html" ) 
		|| die "Couldn't open $outputFilename to read.\n";

	print HTMLFILE "<HTML>\n<HEAD><TITLE>$outputFilename</TITLE></HEAD>\n".
		"<BODY BGCOLOR=\"#ffffff\">\n<PRE>";

	while( <HFILE> )
	{

		s/</\&lt;/g;
		s/>/\&gt;/g;
		s/"/\&quot;/g;

		if( /^\s*(template.*\s+)?(class|struct)/ ) {
			$_ = refString($_);
		}

		print HTMLFILE;
	}

	print HTMLFILE "</PRE>\n<HR>\n<address>$genText",
			"</address>\n</BODY>\n</HTML>\n";
}

#sub Linkize
#{
#	my( $str ) = @_;
#
#	$str =~ s/ /%20/g;
#	$str =~ s/</%3C/g;
#	$str =~ s/>/%3E/g;
#
#	return $str;
#}

1;

