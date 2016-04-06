#!/usr/bin/python

# run this script on all headers of deal.II to fix comment line wrapping for
# doxygen comments.
# Example:
# cd include
# find . -name "*h" -print | while read file;do ../contrib/utilities/wrapcomments.py $file >temp;mv temp $file;done

from __future__ import print_function
import textwrap
import sys, re
wrapper = textwrap.TextWrapper()

# take an array of lines and wrap them to 78 columns and let each line start
# with @p startwith
def wrap_block(lines, startwith):
    longline = " ".join(lines)
    wrapper.initial_indent = startwith
    wrapper.subsequent_indent = startwith
    wrapper.break_long_words = False
    wrapper.width = 78
    return wrapper.wrap(longline)

# strips whitespace and leading "*" from each lines
def remove_junk(lines):
    out = []
    for line in lines:
        line = line.strip()
        if line.startswith("*"):
            line = line[1:].strip()
        out.append(line)
    return out

# returns True if at least one entry in @p list is contained in @p str
def one_in(list, str):
    for li in list:
        if li in str:
            return True
    return False

# returns True if @p str starts with one of the entries in @p list
def starts_with_one(list, str):
    for li in list:
        if str.startswith(li):
            return True
    return False

# take a doxygen comment block "/** bla ... */" given as a list of lines and
# format it in a pretty way rewrapping lines while taking care to keep certain
# structure.
def format_block(lines, infostr=""):
    if len(lines)==1 and "/*" in lines[0] and "*/" in lines[0]:
        return lines

    if (not "/**" in lines[0]
          or not "*/" in lines[-1]):
        print("%s not a code block"%infostr, file=sys.stderr)
        return lines

    for line in lines[1:-1]:
        assert(not "/*" in line)
        assert(not "*/" in line)

    lines[-1]=lines[-1].replace("**/","*/")

    if not lines[0].strip().startswith("/**"):
        print ("%s error, ignoring code block with junk in same line before"%infostr, file=sys.stderr)
        return lines
    if not lines[-1].strip().endswith("*/"):
        print ("%s error, ignoring code block not ending at end of line"%infostr, file=sys.stderr)
        return lines

    if lines[0].strip()!="/**":
        #print ("%s warning code block not starting in separate line"%infostr, file=sys.stderr)
        idx = lines[0].find("/**")
        temp = [lines[0][0:idx+3], lines[0][idx+3:]]
        temp.extend(lines[1:])
        lines = temp
    if lines[-1].strip()!="*/":
        #print ("%s warning code block not ending in separate line"%infostr, file=sys.stderr)
        idx = lines[-1].find("*/")
        temp = lines[0:-1]
        temp.append(lines[-1][0:idx])
        temp.append(lines[-1][idx:])
        lines = temp

    idx = lines[0].find("/**")
    start = lines[0][:idx]+" * "
    
    out = [lines[0].rstrip()]
    idx = 1
    endidx = len(lines)-1
    curlines = []

    ops_startline = ["<li>", "@param", "@returns", "@warning", "@ingroup", "@author", "@date", "@related", "@relates", "@relatesalso", "@deprecated", "@image", "@return", "@brief", "@attention", "@copydoc", "@addtogroup", "@todo", "@tparam", "@see", "@note", "@skip", "@skipline", "@until", "@line", "@dontinclude", "@include"]

    # subset of ops_startline that does not want stuff from the next line appended
    # to this.
    ops_also_end_paragraph = ["@image", "@skip", "@skipline", "@until", "@line", "@dontinclude", "@include"]

    # stuff handled in the while loop down: @code, @verbatim, @f @ref

    ops_separate_line = ["<ol>", "</ol>", "<ul>", "</ul>", "@{", "@}", "<br>"]

    # separate and do not break:
    ops_title_line = ["@page", "@name"]

    #todo:
    #  @arg @c @cond  @em @endcond @f{ @internal @name @post @pre  @sa 

    while idx<endidx:
        if one_in(ops_separate_line, lines[idx]):
            if curlines!=[]:
                out.extend(wrap_block(remove_junk(curlines), start))
                curlines=[]
            thisline = remove_junk([lines[idx]])[0]
            for it in ops_separate_line:
                if it in thisline and thisline!=it:
                    print ("%s warning %s not in separate line"%(infostr, it), file=sys.stderr)
            out.append(start + thisline)
        elif re.match(r'\*\s+- ',lines[idx].strip()):
            # bullet ('-') list
            if curlines!=[]:
                out.extend(wrap_block(remove_junk(curlines), start))
                curlines=[]
            thisline = lines[idx].strip()[2:]
            out.append(start + thisline)
        elif lines[idx].strip().startswith("* ") and re.match(r'\s*\d+.',lines[idx][3:]):
            # numbered list
            if curlines!=[]:
                out.extend(wrap_block(remove_junk(curlines), start))
                curlines=[]
            thisline = lines[idx].strip()[2:]
            out.append(start + thisline)
            
        elif one_in(ops_title_line, lines[idx]):
            # do not break @page, etc.
            if curlines!=[]:
                out.extend(wrap_block(remove_junk(curlines), start))
                curlines=[]
            thisline = remove_junk([lines[idx]])[0]
            
            if not thisline.split(" ")[0] in ops_title_line:
                print ("%s warning title not at start of line"%(infostr), file=sys.stderr)
            out.append(start + thisline.strip())

        elif "@ref" in lines[idx]:
            # @ref link "some long description"
            # is special, and we mustn't break it
            if curlines!=[]:
                out.extend(wrap_block(remove_junk(curlines), start))
                curlines=[]
            thisline = remove_junk([lines[idx]])[0]
            if not thisline.startswith("@ref") and not thisline.startswith("(@ref"):
                print ("%s warning %s not at start of line"%(infostr, "@ref"), file=sys.stderr)

            # format:
            # @ref name "some text" blub
            # or @ref name blurb
            withquotes = thisline.split('"')
            if len(withquotes)==3:
                thisline = withquotes[0]+'"'+withquotes[1]+'"'
                remain = withquotes[2]
                for st in [')', '.', ',', ':']:
                    if remain.startswith(st):
                        thisline = thisline + st
                        remain = remain[1:]

                # do not wrap the @ref line:
                out.append(start + thisline.strip())
                if len(withquotes[0].strip().split(' '))!=2:
                    print ("%s warning @ref line looks broken"%(infostr), file=sys.stderr)     
            elif len(withquotes)==1:
                words = thisline.strip().split(" ")
                if len(words)<2 or len(words[0])==0 or len(words[1])==0:
                    print ("%s warning @ref line looks broken"%(infostr), file=sys.stderr)     
                thisline = words[0] + ' ' + words[1]
                out.append(start + thisline.strip())
                remain = " ".join(words[2:])
            else:
                print ("%s warning @ref quotes are not in single line"%(infostr), file=sys.stderr)
                remain = ''

            if len(remain)>0:
                curlines.append(remain)
            
        elif one_in(ops_startline, lines[idx]):
            if curlines!=[]:
                out.extend(wrap_block(remove_junk(curlines), start))
                curlines=[]
            thisline = remove_junk([lines[idx]])[0]
            if not starts_with_one(ops_startline, thisline):
                for it in ops_startline:
                    if it in thisline:
                        print ("%s warning %s not at start of line"%(infostr, it), file=sys.stderr)
            if one_in(ops_also_end_paragraph, lines[idx]):
                out.append(lines[idx].rstrip())
            else:
                curlines.append(lines[idx])
        elif one_in(["@code", "@verbatim", "@f[", "@f{"], lines[idx]):
            if "@f{" in lines[idx]:
                if not lines[idx].endswith("}{"):
                    print ("%s warning malformed @f{*}{"%(infostr), file=sys.stderr)
            if curlines!=[]:
                out.extend(wrap_block(remove_junk(curlines), start))
                curlines=[]
            while True:
                thisline = lines[idx].rstrip()
                if thisline.strip()=="*":
                    thisline = start
                elif thisline.strip()=="@code" or thisline.strip()=="@endcode":
                    thisline = start + thisline.strip()
                elif thisline.strip()[0:2]!="* ":
                    if thisline[0:len(start)].strip()=="":
                        # just a missing *, so keep old indentation
                        thisline = start + thisline[len(start):]
                    else:
                        # no way to recover indentation:
                        print ("%s Error: wrong formatting inside @code block"%infostr, file=sys.stderr)
                        thisline = start + thisline.strip()
                else:
                    thisline = start + thisline.strip()[2:]
                out.append(thisline.rstrip())
                if one_in(["@endcode", "@endverbatim", "@f]", "@f}"], lines[idx]):
                    break
                idx += 1
        elif lines[idx].strip()=="*":
            if curlines!=[]:
                out.extend(wrap_block(remove_junk(curlines), start))
                curlines=[]
            out.append(start[:-1]) #skip whitespace at the end
        else:
            curlines.append(lines[idx])
        idx += 1

    if curlines!=[]:
        out.extend(wrap_block(remove_junk(curlines), start))        

    out.append(start[0:-2] + lines[-1].strip())

    return out


# test the routines:
lineI = ["   * blub", \
         "   * two three ", \
         "   * four"]
lineO = ["blub", \
         "two three", \
         "four"]
assert(remove_junk(lineI)==lineO)

lineI = ["blub", \
         "two three", \
         "four"]
lineO = ["   * blub two three four"]
assert(wrap_block(lineI,"   * ")==lineO)

lineI = [" * 1 2 3 4 5 6 7 9 0 1 2 3 4 5 6 7 9 0 1 2 3 4 5 6 7 9 0 1 2 3 4 5 6 7 9 0 1 2",\
         " * A 4 5 6 7 9 0 1 2 3 4 5 6 7 9 0"]
assert(wrap_block(remove_junk(lineI)," * ")==lineI)

lineI = [" * Structure which is passed to Triangulation::create_triangulation. It",\
         " * contains all data needed to construct a cell, namely the indices of the",\
         " * vertices and the material indicator."]
assert(wrap_block(remove_junk(lineI)," * ")==lineI)

lineI = ["  /**", \
         "   * blub", \
         "   * two three", \
         "   * four", \
         "   */"]
lineO = ["  /**", \
         "   * blub two three four", \
         "   */"]
assert(format_block(lineI)==lineO)

lineI = ["  /**", \
         "   * blub", \
         "   * two three", \
         "   * ", \
         "   * four", \
         "   */"]
lineO = ["  /**", \
         "   * blub two three", \
         "   *", \
         "   * four", \
         "   */"]
assert(format_block(lineI)==lineO)

lineI = ["  /**", \
         "   * blub", \
         "   * @code", \
         "   *   two three", \
         "   * @endcode", \
         "   * four", \
         "   */"]
assert(format_block(lineI)==lineI)

lineI = ["  /**", \
         "   * blub", \
         "   * @code ", \
         "   *   two three", \
         "   *   two three  ", \
         "   * ", \
         "   *   two three", \
         "   * @endcode ", \
         "   * four", \
         "   */"]
lineO = ["  /**", \
         "   * blub", \
         "   * @code", \
         "   *   two three", \
         "   *   two three", \
         "   *", \
         "   *   two three", \
         "   * @endcode", \
         "   * four", \
         "   */"]
assert(format_block(lineI)==lineO)
 
lineI = ["  /**", \
         "   * blub", \
         "       @code ", \
         "       two three", \
         "      @endcode ", \
         "   * four", \
         "   */"]
lineO = ["  /**", \
         "   * blub", \
         "   * @code", \
         "   *   two three", \
         "   * @endcode", \
         "   * four", \
         "   */"]
assert(format_block(lineI)==lineO)


lineI = ["  /**", \
         "   * blub", \
         "   * @code ", \
         "    *   two three", \
         "   *   two three  ", \
         " * ", \
         "       two three", \
         "    ", \
         "   * @endcode ", \
         "   * four", \
         "   */"]
lineO = ["  /**", \
         "   * blub", \
         "   * @code", \
         "   *   two three", \
         "   *   two three", \
         "   *", \
         "   *   two three", \
         "   *", \
         "   * @endcode", \
         "   * four", \
         "   */"]
assert(format_block(lineI)==lineO)



lineI = ["  /**", \
         "   * blub", \
         "   * <ul>", \
         "   * <li> bla", \
         "   * <li> blub", \
         "   * </ul>", \
         "   */"]
assert(format_block(lineI)==lineI)
 
lineI = ["    /** @addtogroup Exceptions", \
         "     * @{ */"]
lineO = ["    /**", \
         "     * @addtogroup Exceptions", \
         "     * @{", \
         "     */"]
assert(format_block(lineI)==lineO)

lineI = [" /** ", \
         "  *   bla", \
         "  * @image testing.png", \
         "  *  blub", \
         "  */"]
lineO = [" /**", \
         "  * bla", \
         "  * @image testing.png", \
         "  * blub", \
         "  */"]
assert(format_block(lineI)==lineO)

lineI = [" /** ", \
         "  * @ref a b c d e", \
         "  * @ref a \"b c\" d dd", \
         "  */"]
lineO = [" /**", \
         "  * @ref a", \
         "  * b c d e", \
         "  * @ref a \"b c\"", \
         "  * d dd", \
         "  */"]
assert(format_block(lineI)==lineO)

long = "long "*20
lineI = [" /** ", \
         "  * @ref a \""+long+"\" d", \
         "  */"]
lineO = [" /**", \
         "  * @ref a \""+long+"\"", \
         "  * d", \
         "  */"]
assert(format_block(lineI)==lineO)

lineI = [" /** ", \
         "  * @ref a. c", \
         "  * @ref a \"b c\". c2", \
         "  */"]
lineO = [" /**", \
         "  * @ref a.", \
         "  * c", \
         "  * @ref a \"b c\".", \
         "  * c2", \
         "  */"]
assert(format_block(lineI)==lineO)

# do not break @page:
longtext = "bla bla"*20
lineI = [" /**", \
         "  * @page " + longtext, \
         "  * hello", \
         "  */"]
assert(format_block(lineI)==lineI)

# do not break $very_long_formula_without_spacing$:
longtext = "blabla"*20
lineI = [" /**", \
         "  * a $" + longtext + "$", \
         "  */"]
lineO = [" /**", \
         "  * a", \
         "  * $" + longtext + "$", \
         "  */"]
assert(format_block(lineI)==lineO)

# nested lists
lineI = [" /**", \
         "  * Hello:", \
         "  * - A", \
         "  * - B", \
         "  *   - C", \
         "  *   - D", \
         "  * - E", \
         "  *         - very indented", \
         "  */"]
lineO = lineI
assert(format_block(lineI)==lineO)

# @f{}
lineI = [" /**", \
         "  * Hello:", \
         "  * @f{aligned*}{", \
         "  *   A\\\\", \
         "  *   B", \
         "  * @f}", \
         "  * bla", \
         "  */"]
lineO = lineI
assert(format_block(lineI)==lineO)

# @until 
lineI = [" /**", \
         "  * Hello:", \
         "  * @include a", \
         "  * bla", \
         "  * @dontinclude a", \
         "  * bla", \
         "  * @line a", \
         "  * bla", \
         "  * @skip a", \
         "  * bla", \
         "  * @until a", \
         "  * bla", \
         "  */"]
lineO = lineI
assert(format_block(lineI)==lineO)

# lists
lineI = [" /**", \
         "  * Hello:", \
         "  *  - a", \
         "  *  - b", \
         "  * the end.", \
         "  */"]
lineO = lineI
assert(format_block(lineI)==lineO)

# numbered lists
lineI = [" /**", \
         "  * Hello:", \
         "  * 1. a", \
         "  * 2. b", \
         "  * the end.", \
         "  */"]
lineO = lineI
assert(format_block(lineI)==lineO)

# do not break @name:
longtext = "bla bla"*20
lineI = [" /**", \
         "  * @name " + longtext, \
         "  * hello", \
         "  */"]
assert(format_block(lineI)==lineI)


#print (lineI)
#print (format_block(lineI))




# now open the file and do the work


args=sys.argv
args.pop(0)

if len(args)!=1:
    print("Usage: wrapcomments.py infile >outfile")
    exit(0)

f = open(args[0])
lines = f.readlines()
f.close()

out = []
cur = []
inblock = False
lineidx = 0
for line in lines:
    lineidx += 1
    line = line.replace("\n","")
    if not inblock and "/**" in line:
        inblock = True
        cur = []
    if inblock:
        cur.append(line)
    else:
        out.append(line)
    if inblock and "*/" in line:
        out.extend(format_block(cur, args[0]+":%d"%lineidx))
        cur = []
        inblock = False

assert(cur==[])

for line in out:
    print (line)





