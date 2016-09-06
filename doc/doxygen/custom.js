$( document ).ready( function() {
    // replace the text inside every div with class 'doxygen-generated-exception-message' to be more readable:
    $('.doxygen-generated-exception-message').each(function(i, obj) {
	var s = $(this).html();
	s=s.replace(/\\n|std::endl/g,"<br/>");
	s=s.replace(/\\"/g,"'");
	s=s.replace(/"|&lt;&lt;/g,"");
	s=s.replace(/arg1/g,"<i>arg1</i>");
	s=s.replace(/arg2/g,"<i>arg2</i>");
	s=s.replace(/arg3/g,"<i>arg3</i>");
	s=s.replace(/arg4/g,"<i>arg4</i>");
	s=s.replace(/arg5/g,"<i>arg5</i>");
	$(this).html(s);
    });
});
