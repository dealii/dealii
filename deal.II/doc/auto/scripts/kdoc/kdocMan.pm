#	kdocMan	-- man page output for kdoc.
#	Copyright(c) 1998, Sirtaj Singh Kang (taj@kde.org).

#	$Id$

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

package kdocMan;

require Ast;

BEGIN {
	$modname   = "kdocMan";
	$outputdir = ".";
	$genText   = "";
	$className = "";
}

sub dumpDoc
{
	( $lib, $root, $outputdir ) = @_;

	print "Generating man pages.\n" unless $main::quiet;
	$genText = $main::generatedByText;
	$datestamp = `date +%D`;
	chop $datestamp;

	writeIndex();
	writeClasses();
}

sub writeIndex()
{
	@classNames	= ();
	%classNodes	= ();

	foreach $class ( @{$root->{'classes'}} ) {
		push( @classNames, $class->{'astNodeName'} );
		$classNodes{ $class->{'astNodeName'} } = $class;
	}

	@classNames = sort @classNames;

	open( MANPAGE, ">$outputdir/$lib.1$lib" ) 
		|| die "Couldn't write to $outputdir/$lib.1$lib\n";

print MANPAGE<<EOF;
.TH $lib 1 "$datestamp" "KDOC"
.SH NAME
$lib - C++ Class Library
.SH SYNOPSIS
The $lib library consists of the following classes:
.PP
EOF
	writeHierarchy("<>");

print MANPAGE<<EOF;	
.SH DESCRIPTION
EOF

	foreach $class ( @classNames ) {
		$short = $classNodes{$class}->{ "ClassShort" };

		print MANPAGE "\n.TP 4\n.BI \"$class\"\n$short\n.PP\n";
	}


print MANPAGE<<EOF;
.SH AUTHOR
$genText
.PP
EOF
	close MANPAGE;

}

sub writeClasses
{
	foreach $classNode ( @{$root->{'classes'}} ) {
		$classNode->Visit( $modname );
		$className=$astNodeName;
		writeClass();
		Ast::UnVisit();
	}
}

sub writeClass
{
	($fname = $astNodeName) =~ tr/[A-Z]/[a-z]/;
	$fname .= ".3$lib";

	$className = $astNodeName;
	$Description = makeReferencedText( $Description );
	$Description =~ s/\n\s+/\n/g;

	open ( CLASSMAN, ">$outputdir/$fname" )
		|| die "couldn't open $outputdir/$fname for writing.\n";
	
	open( REDIR, ">$outputdir/$astNodeName.3$lib" )
		||die "couldn't open $outputdir/$astNodeName.3$lib ".
		" for writing.\n";

	print REDIR ".so $fname\n";
	close REDIR;


print CLASSMAN<<EOF;
.TH $className 3 "$datestamp\" "KDOC" "$lib Class Library"
.SH NAME
$className - $ClassShort

.SH SYNOPSIS 
#include <$Header>

EOF

# Internal or Deprecated

print ".ti +1c\n.b \"Deprecated Class\"\n" if $Deprecated; 
print ".ti +1c\n.b \"Internal use only\"\n" if $Internal; 


print CLASSMAN<<EOF;
Inherits: 
EOF

# Inheritance list
	my($pcount) = 0;
	if (defined @{$Ancestors} ) {
		foreach $foreparent ( @{$Ancestors} ) {
			print CLASSMAN ", " if $pcount != 0;
			$pcount = 1;

			print CLASSMAN $foreparent->{"astNodeName"};

			if( defined $main::classSource{
				$foreparent->{"astNodeName"}} ) {
				print CLASSMAN " (",
				$main::classSource{$foreparent->{"astNodeName"}},")";
			}
		}
	}
	else {
		print CLASSMAN "nothing";
	}

print CLASSMAN<<EOF;

.PP
EOF
	listMembers( "Public Members", $public );
        listMembers( "Public Slots", $public_slots );
        listMembers( "Protected Members", $protected );
        listMembers( "Protected Slots", $protected_slots );
        listMembers( "Signals", $signals );

	if( $main::dontDoPrivate == 0 ) {
                listMembers( "Private Members", $private );
                listMembers( "Private Slots", $private_slots );
        } 
                                   

print CLASSMAN<<EOF;

.SH DESCRIPTION
$Description
.SH MEMBER DOCUMENTATION
EOF

	documentMembers( $public );
        documentMembers( $public_slots );
        documentMembers( $protected );
        documentMembers( $protected_slots );
        documentMembers( $signals );

	if( $main::dontDoPrivate == 0 ) {
                documentMembers( $private );
                documentMembers( $private_slots );
        } 

	if( $ClassSee ne "" ) {
		$ClassSee = ", ".$ClassSee
	}

print CLASSMAN<<EOF;
.SH "SEE ALSO"
$lib(1)$ClassSee
EOF

print CLASSMAN<<EOF;
.SH AUTHOR
EOF

	if( ! ($Author =~ /^\s*$/) ) {
		print CLASSMAN "$Author\n.PP\n";
	}

print CLASSMAN<<EOF;
$genText
.SH VERSION
$Version
EOF

	close CLASSMAN;
}

sub writeHierarchy
{
        local( $node ) = @_;
        local( @kids );

# Display class
	if( $node ne "<>" ) {
		print MANPAGE "$node";
		if( defined $main::classSource{ $node } ) {
			print MANPAGE " (from ",$main::classSource{$node},")";
		}
		print MANPAGE "\n";
	}

# Recurse for each child class
        if ( defined $main::children{$node} ) {
		print MANPAGE ".in +1c\n";
                foreach $kid ( sort split /;/, $main::children{$node} )
                {
			print MANPAGE ".br\n";
                        writeHierarchy( $kid );
                }

		print MANPAGE ".in -1c\n";
         }
}
 
sub makeReferencedText
{
	my( $str ) = @_;

	$str =~ s/\@ref(\s+[-:#\w\._]+)/\n.BI "$1"\n/g;

	return $str;
}

sub listMembers
{
	my( $access, $nodelist ) = @_;
	my( $node );

	return if !defined $nodelist || $#{$nodelist} == -1;

print CLASSMAN<<EOF;
.PP
.BI "$access"
.in +2c
EOF
	
	foreach $node ( @{$nodelist} ) {
		$node->Visit( $modname );

		print CLASSMAN ".ti -1c\n";

		if( $Keyword eq "property" ) {
			$Type =~ s/^\s+//g;
			$Type =~ s/\s+$//g;
			print CLASSMAN ".BI \"$Type $astNodeName\"\n.br\n";
		}
		elsif( $Keyword eq "method" ) {
			$ReturnType =~ s/^\s+//g;
			$ReturnType =~ s/\s+$//g;
			$ReturnType .= " " if $ReturnType ne "";

			print CLASSMAN ".BI \"$ReturnType$astNodeName ",
			"($Parameters) $Const\"\n.br\n";
		}
		elsif( $Keyword eq "enum" ) {
			print CLASSMAN ".BI \"",
			"enum $astNodeName {$Constants}\"\n.br\n";
		}
		elsif( $Keyword eq "typedef" ) {
			print CLASSMAN ".BI \"",
			"typedef $astNodeName\"\n.br\n";
		}
		Ast::UnVisit();
	}

print CLASSMAN<<EOF;
.in -2c
EOF
}

sub documentMembers
{
	my( $nodelist ) = @_;

	return if !defined $nodelist || $#{$nodelist} == -1;

	foreach $member ( @{$nodelist} ) {
		$member->Visit( $modname );
 
		if( $Description eq "" && $See eq "" ) {
			Ast::UnVisit();
			next;
		}

		if( $Keyword eq "property" ) {
			$Type =~ s/^\s+//g;
			print CLASSMAN ".SH \"$Type $className","::",
				"$astNodeName\"\n";
		}
		elsif( $Keyword eq "method" ) {
			$ReturnType =~ s/^\s+//g;
			$ReturnType =~ s/\s+$//g;
			$ReturnType .= " " if $ReturnType ne "";

			print CLASSMAN ".SH \"$ReturnType$className",
				"::","$astNodeName ","($Parameters) $Const\"\n";
		}
		elsif( $Keyword eq "enum" ) {
			print CLASSMAN ".SH \"",
				"enum $className","::",
				"$astNodeName {$Constants}\"\n";
		}
		elsif( $Keyword eq "typedef" ) {
			print CLASSMAN ".SH \"", "typedef $className",
				"::","$astNodeName\"\n";
		}

# Internal or Deprecated

		print ".ti +1c\n.b \"Deprecated member\"\n" if $MethDeprecated; 
		print ".ti +1c\n.b \"Internal use only\"\n" if $MethInternal; 

# Description
		print CLASSMAN makeReferencedText($Description),
			"\n.br\n" if $Description ne "";

# Parameters
		if( $Keyword eq "method" ) {
			if( $#{$ParamDoc} != -1 ) {
				print CLASSMAN "\n.br\n.in +1c\n",
					"Parameters\n.in +1c\n.";
				
				foreach $parameter ( @{$ParamDoc} ) {
print CLASSMAN ".br\n.BI \"", $parameter->{"astNodeName"}, 
	"\"\n.in +1c\n",
		$parameter->{"Description"}, "\n.in -1c\n";
				}
print CLASSMAN<<EOF;
.in -1c
.in -1c
EOF
			}

			print CLASSMAN ".in +1c\nReturns\n.in +1c\n",
				"$Returns\n.in -1c\n.in -1c\n" if $Returns ne "";
			print CLASSMAN ".in +1c\nThrows\n.in +1c\n",
				"$Exceptions\n.in -1c\n.in -1c\n" 
				if $Exceptions ne "";
		}

		print CLASSMAN ".in +1c\nSee Also\n.in +1c\n$See\n",
			".in -1c\n.in -1c\n" if $See ne "";

		Ast::UnVisit();
	}
}

1;
