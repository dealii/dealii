#!/usr/bin/python

## ---------------------------------------------------------------------
##
## Copyright (C) 2016 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

#
# Written by timo.heister@gmail.com
#
# This is an (incomplete) script that parses test output from ctest's xml
# files and generates a static html page from it.
#
# Execute in the directory build/Testing

import xml.etree.ElementTree as ET
import glob

class Group:
    def __init__(self, name):
        self.name = name
        self.n_tests = 0
        self.n_fail = 0
        self.fail = []
        self.fail_text = {}
        self.fail_status = {}
        self.n_status = [0,0,0,0,0]

class Revision:
    def __init__(self):
        self.groups = {}
        self.number = -1
        self.name = ''
        self.n_tests = 0
        self.n_fail = 0

def parse_revision(dirname):
    rev = Revision()

    if len(glob.glob(dirname+'/Update.xml'))>0:
        #new format
        tree = ET.parse(dirname+'/Update.xml')
        rev.name = tree.getroot().find('BuildName').text
        rev.number = tree.getroot().find('Revision').text
    elif len(glob.glob(dirname+'/Notes.xml'))>0:
        #old format
        tree = ET.parse(dirname+'/Notes.xml')
        rev.name = tree.getroot().attrib['BuildName']
        number = rev.name.split('-')[-1]
        rev.number = number[1:]
    else:
        return None
        
    print dirname, "BUILD: ", rev.name

    #now Test.xml:
    tree = ET.parse(dirname+'/Test.xml')
    root = tree.getroot()
    testing = root.find('Testing')

    for test in testing.findall("Test"):
        fail=False
        if test.attrib['Status']=="failed": fail=True
        name = test.find('Name').text
        group = name.split('/')[0]
        status = 4
        if fail:
            text = test.find('Results').find('Measurement').find('Value').text
            if text == None: text=""
            failtext = text.encode('utf-8')
            failtextlines = failtext.replace('"','').split('\n')
            failstatustxt = failtextlines[0].split(' ')[-1]
            for i in range(0,len(failtextlines)):
                failtextlines[i] = failtextlines[i][0:80]
                if failtextlines[i].startswith('FAILED: '): failtextlines[i]='FAILED: ...';
            failtext = '\n'.join(failtextlines[4:min(25,len(failtext))])
            statuslist=['CONFIGURE','BUILD','RUN','DIFF']
            if failstatustxt in statuslist:
                status = statuslist.index(failstatustxt)
            else:
                print "unknown status '%s' in test %s "% (failstatustxt,name)
                status=0           

        if not group in rev.groups:
            rev.groups[group]= Group(group)
            
        rev.groups[group].n_tests += 1
        rev.n_tests += 1
        rev.groups[group].n_status[status] += 1
        if fail: 
            rev.groups[group].n_fail += 1
            rev.n_fail += 1
            rev.groups[group].fail.append(name)
            rev.groups[group].fail_text[name]=failtext
            rev.groups[group].fail_status[name]=status
        
    for g in sorted(rev.groups):
        g = rev.groups[g]
        #print g.name, g.n_tests, g.n_fail, g.fail
        
    return rev




#from xml.dom import minidom


n=glob.glob("*/Build.xml")
n.sort(reverse=True)
numberofrevisions=10
n = n[0:min(10,len(n))]

revs = []

allgroups = set()

for f in n:
    dirname = f.replace('/Build.xml','')
    rev = parse_revision(dirname)
    if rev!=None:
        revs.append(rev)
        for gr in rev.groups:
            allgroups.add(gr)

revs.sort(key=lambda x: x.number, reverse=True)
    
allgroups = sorted(allgroups)

f = open('tests.html', 'w')
f.write('<html><head></head>')
f.write("""<style type="text/css">
table {
border-collapse:collapse;
}
.fail {background-color:red;}
.togglebody tr {background-color:#CCC}
.test4 {
    background-color: #90ff80;
text-align: right;
}

.test3 {
    background-color: #FFFF00;
text-align: right;
}

.test2 {
    background-color: #FFA000;
text-align: right;
}

.test1 {
    background-color: #FF2020;
text-align: right;
}

.test0 {
    background-color: #C030D0;
text-align: right;
}

thead{
border: 2px solid;
}

.groupALL {
text-align: center;
}

.colgroup {
border: 2px solid;
}

.onerow {background: #FFE}
.otherrow {background: #EEE}

</style>""")
f.write("""<script>
function toggle_id(obj)
{
var e = document.getElementById(obj);
if(e.style.display != 'none')
 e.style.display = 'none';
else
 e.style.display = '';
}
</script>
<body>""")
f.write('<table>')

f.write('<colgroup span="1" class="colgroup""/>')
for rev in revs:    
    f.write('<colgroup span="5" class="colgroup"/>')
f.write('\n')


f.write('<thead><tr>')
f.write('<th style="width:250px">&nbsp;</th>')

for rev in revs:
    f.write('<th colspan="5"><a href="http://www.dealii.org/websvn/revision.php?repname=deal.II+Repository&rev=%s">r%s</th>'%(rev.number,rev.number))
f.write('</tr></thead>\n')

f.write('<tbody><tr>')
f.write('<td>ALL</td>')
for rev in revs:
    if (rev.n_fail>0):
        f.write('<td colspan="5" class="groupALL"><span class="fail">' + str(rev.n_fail) + '</span> / ' + str(rev.n_tests) + '</td>')
    else:
        f.write('<td colspan="5" class="groupALL">' + str(rev.n_fail) + ' / ' + str(rev.n_tests) + '</td>')
f.write('</tr></tbody>\n')

#second header
f.write('<tbody><tr style="border-bottom: 2px solid">')
f.write('<td></td>')
for rev in revs:
    for c in range(0,5):
        
        titles=['Configure','Build','Run','Diff','Pass']
        caption=['C','B','R','D','P']
        f.write('<td title="%s" class="test%d">%s</td>'%(titles[c],c,caption[c]))
f.write('</tr></tbody>\n')

counter=0
for group in allgroups:
    counter+=1
    if counter % 2==0:
        f.write('<tbody class="onerow"><tr>')
    else:
        f.write('<tbody class="otherrow"><tr>')

    failing = set()
    for rev in revs:
        if group in rev.groups:
            failing |= set(rev.groups[group].fail)
    failing = sorted(failing)

    if (len(failing)>0):
        f.write('<td><a href="#" onclick="toggle_id(\'group:%s\');return false">%s</a></td>'%(group,group))
    else:
        f.write('<td>' + group + '</td>')

    for rev in revs:
        if group not in rev.groups:
            for c in range(0,5):
                f.write('<td></td>')
        else:
            gr = rev.groups[group]
            for c in range(0,5):
                if gr.n_status[c]==0:
                    f.write('<td></td>')
                else:
                    f.write('<td class="test%d">%d</td>'%(c,gr.n_status[c]))

    f.write('</tr></tbody>\n')

                
    #failing tests in group:
    if len(failing)>0:
        f.write('<tbody class="togglebody" style="display:none" id="group:%s">'%group)
        for fail in failing:
            f.write('<tr>')
            name = fail[len(group):] # cut off group name
            f.write('<td>&nbsp;' + name + '</td>')
            for rev in revs:
                if group in rev.groups and fail in rev.groups[group].fail:
                    status=rev.groups[group].fail_status[fail]
                    text=rev.groups[group].fail_text[fail]
                    for c in range(0,5):
                        if c==status:
                            f.write('<td class="test%d" title="%s">X</td>'%(c,text))
                        else:
                            f.write('<td></td>')

                else:
                    f.write('<td colspan="5"></td>')

            f.write('</tr>\n')
        f.write('</tbody>\n')
                
    f.write('\n\n')


f.write('</table>')


f.write('</body></html>')


