<Title>
deal.II Mailinglist
</Title>

<!-- General markup first -->


<!-- Start of main index page -->

<IdxPgBegin>
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML//EN">
<HTML>
  <HEAD>
    <TITLE>$IDXTITLE$</TITLE>
    <link href="../screen.css" rel="StyleSheet" title="deal.II Homepage" media="screen">
    <link href="../print.css" rel="StyleSheet" title="deal.II Homepage" media="print">
    <meta name="keywords" content="deal.II">
  </HEAD>
  <BODY style="background-image:url(../pictures/title-background.jpg);" lang="en">
    <H1>$IDXTITLE$</H1>
</IdxPgBegin>

<!-- style="background-image:url(../pictures/title-background.jpg);" lang="en" -->

<!-- End of main index page -->

<IdxPgEnd>
<ADDRESS><A HREF="mailto:deal@iwr.uni-heidelberg.de">The deal.II mailing list</ADDRESS>
</BODY>
</HTML>
</IdxPgEnd>


<!-- Start of thread index page -->

<TIdxPgBegin>
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML//EN">
<HTML>
  <HEAD>
    <TITLE>$TIDXTITLE$</TITLE>
    <link href="../screen.css" rel="StyleSheet" title="deal.II Homepage" media="screen">
    <link href="../print.css" rel="StyleSheet" title="deal.II Homepage" media="print">
    <meta name="keywords" content="deal.II">
  </HEAD>
  <BODY style="background-image:url(../pictures/title-background.jpg);" lang="en">
  <H1>$TIDXTITLE$</H1>
</TIdxPgBegin>

<!-- End of thread index page -->

<TIdxPgEnd>
<ADDRESS><A HREF="mailto:deal@iwr.uni-heidelberg.de">The deal.II mailing list</ADDRESS>
</BODY>
</HTML>
</TIdxPgEnd>

<!-- Start of message page -->

<MsgPgBegin>
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML//EN">
<HTML>
  <HEAD>
    <TITLE>$SUBJECTNA:72$</TITLE>
    <link href="../screen.css" rel="StyleSheet" title="deal.II Homepage" media="screen">
    <link href="../print.css" rel="StyleSheet" title="deal.II Homepage" media="print">
  </HEAD>
  <BODY>
</MsgPgBegin>

<!-- End of message page -->

<MsgPgEnd>
<HR>
<ADDRESS><A HREF="mailto:deal@iwr.uni-heidelberg.de">The deal.II mailing list</ADDRESS>
</BODY>
</HTML>
</MsgPgEnd>




<!-- Sorting and filtering -->


<Sort>
<Reverse>
<TSort>
<TReverse>
<SpamMode>
<FromFields>
reply-to:from:return-path:apparently-from:sender:resent-sender
</FromFields>

<!-- List styles -->


<ListBegin>
<UL>
<LI><A HREF="$TIDXFNAME$">Thread Index</A></LI>
</UL>
<HR>
</ListBegin>

<DayBegin>
<p><strong>$MSGGMTDATE(CUR;%Y - %b - %d)$</strong><br>
<table border=0 width=100%>
</DayBegin>

<DayEnd>
</table>
</DayEnd>

<!--    Define LITEMPLATE to display the time of day the message was
        sent, message subject, author, and any annotation for the
        message.
  -->
<LiTemplate>
<tr valign="top">
<td align="left">$SUBJECT$</td>
<td align="right">$FROMNAME$</td>
</tr>
</LiTemplate>

<!--    Define LISTEND to close table  -->

<ListEnd>
</ListEnd>


<!-- Message formatting -->


SubjectHeader>
<H1>$SUBJECTNA$</H1>
</SubjectHeader>

<FieldsBeg>
</FieldsBeg>

<FieldsEnd>
</FieldsEnd>

<LabelBeg>
 <!--
</LabelBeg>

<FldEnd>
 -->
</FldEnd>

<HeadBodySep>
<p>
<strong>$FROMNAME$</strong> ($FROMADDRNAME$ <I>at</I> $FROMADDRDOMAIN$)
<BR>
<strong>$GMTDATE$</strong></p>
<HR>
</HeadBodySep>