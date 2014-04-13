#!/usr/bin/python

# run this script on all headers of deal.II to fix comment line wrapping for
# doxygen comments.
# Example:
# cd include
# find . -name "*h" -print | while read file;do ../contrib/utilities/wrapcomments.py $file >temp;mv temp $file;done

from __future__ import print_function
import textwrap
import sys
import string
wrapper = textwrap.TextWrapper()

# take an array of lines and wrap them to 78 columns and let each lines start
# with @p startwith
def wrap_block(lines, startwith):
    longline = " ".join(lines)
    wrapper.initial_indent = startwith
    wrapper.subsequent_indent = startwith
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
        idx = string.find(lines[0],"/**")
        temp = [lines[0][0:idx+3], lines[0][idx+3:]]
        temp.extend(lines[1:])
        lines = temp
    if lines[-1].strip()!="*/":
        #print ("%s warning code block not ending in separate line"%infostr, file=sys.stderr)
        idx = string.find(lines[-1],"*/")
        temp = lines[0:-1]
        temp.append(lines[-1][0:idx])
        temp.append(lines[-1][idx:])
        lines = temp

    idx = string.find(lines[0],"/**")
    start = lines[0][:idx]+" * "
    
    out = [lines[0]]
    idx = 1
    endidx = len(lines)-1
    curlines = []

    ops_startline = ["<li>", "@param", "@returns", "@warning", "@ingroup", "@author", "@date", "@related", "@deprecated"]

    ops_separate_line = ["<ol>", "</ol>", "<ul>", "</ul>", "@{", "@}", "<br>"]

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
        elif one_in(ops_startline, lines[idx]):
            if curlines!=[]:
                out.extend(wrap_block(remove_junk(curlines), start))
                curlines=[]
            thisline = remove_junk([lines[idx]])[0]
            if not starts_with_one(ops_startline, thisline):
                for it in ops_startline:
                    if it in thisline:
                        print ("%s warning %s not at start of line"%(infostr, it), file=sys.stderr)
            curlines.append(lines[idx])
        elif one_in(["@code", "@verbatim", "@f["], lines[idx]):
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
                if one_in(["@endcode", "@endverbatim", "@f]"], lines[idx]):
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





