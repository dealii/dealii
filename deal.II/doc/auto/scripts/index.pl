@entries = ();

while(<>)
{
    next if (m/DeclException/);
    chop;
    if (m/(.+)=(.+)/)
    {
	@n = split "::", $1;

	$entry = pop @n;
	while ($#n>=0)
	{
	    $entry .= "+" . pop @n;
	}
	$entry .= '=' . $2 . '=' . $library;
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
<TABLE CELLPADDING=9 CELLSPACING=0 RULES=ROWS BORDER=6>
<TR><TH><B>function/class</B><TD ALIGN=CENTER><B>in class</B></TR>
EOT
    ;
foreach $entry (sort @entries)
{
    @l = split "=", $entry;
    @c = split '\+', @l[0];
    @f = split "#", @l[1];

    print "<TR><TH ALIGN=left><A HREF = \"$l[2]/@f[0]\">@c[0]</A><TD>@c[1]</TR>\n";
}
print << 'EOT'
</TABLE>
</BODY>
</HTML>
EOT
    ;
