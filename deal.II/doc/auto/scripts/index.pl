@entries = ();

while(<>)
{
    next if (m/DeclException/);
    chop;
    s/\</\&lt\;/g;
    s/\>/\&gt\;/g;
    if (m/(.+)=([^=]+)=?$/)
    {
	@n = split "::", $1;

	$entry = pop @n;
	while ($#n>=0)
	{
	    $entry .= "+" . pop @n;
	}
	$entry .= ',' . $2 . ',' . $library;
	push @entries, $entry;
    }
    else
    {
	$library = $_;
    }
}

print << 'EOT'
<HTML><HEAD><TITLE>DEAL Name Index</TITLE></HEAD>
<BODY><H1>Name Index of DEAL</H1><hr>
Remark: Lowercase is sorted behind uppercase and destructors are isolated at
the end.
<hr>
<TABLE CELLPADDING=5 CELLSPACING=0 RULES=ROWS BORDER=2>
<TR><TH WIDTH=300><B>function/class</B><TH ALIGN=CENTER WIDTH=300><B>in class</B></TR>
</TABLE>
EOT
    ;
foreach $entry (sort @entries)
{
    @l = split ",", $entry;
    @c = split '\+', @l[0];
    @f = split "#", @l[1];

    print '<TABLE CELLPADDING=5 CELLSPACING=0 RULES=ROWS BORDER=2>';
    print "\n<TR><TH ALIGN=left WIDTH=300><A HREF = \"$l[2]/@f[0]\">@c[0]</A><TD WIDTH=300>@c[1]</TR></TABLE>\n";
}
print << 'EOT'
</BODY>
</HTML>
EOT
    ;
