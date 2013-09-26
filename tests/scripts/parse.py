import xml.etree.ElementTree as ET
import glob

class Group:
    def __init__(self, name):
        self.name = name
        self.n_tests = 0
        self.n_fail = 0
        self.fail = []
        self.fail_text = {}

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
        status = test.attrib['Status']
        fail=False
        if status=="failed": fail=True
        name = test.find('Name').text
        group = name.split('/')[0]

        if not group in rev.groups:
            rev.groups[group]= Group(group)
            
        rev.groups[group].n_tests += 1
        rev.n_tests += 1
        if fail: 
            rev.groups[group].n_fail += 1
            rev.n_fail += 1
            rev.groups[group].fail.append(name)
            failtext = test.find('Results').find('Measurement').find('Value').text.encode('utf-8')
            failtext = failtext.replace('"','').split('\n')
            failtext = '\n'.join(failtext[4:min(20,len(failtext))])
            rev.groups[group].fail_text[name]=failtext
        
    for g in sorted(rev.groups):
        g = rev.groups[g]
        #print g.name, g.n_tests, g.n_fail, g.fail
        
    return rev




#from xml.dom import minidom


n=glob.glob("*/Build.xml")
n.sort(reverse=True)
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
table,td,th {
border:1px solid black; 
}
.fail {background-color:red;}
.togglebody tr {background-color:#eeeeee}

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

f.write('<thead><tr>')
f.write('<th style="width:250px">&nbsp;</th>')

for rev in revs:
    f.write('<th>r'+rev.number+'</th>')
f.write('</tr></thead>\n')

f.write('<tbody><tr>')
f.write('<td>ALL</td>')
for rev in revs:
    if (rev.n_fail>0):
        f.write('<td><span class="fail">' + str(rev.n_fail) + '</span> / ' + str(rev.n_tests) + '</td>')
    else:
        f.write('<td>' + str(rev.n_fail) + ' / ' + str(rev.n_tests) + '</td>')
f.write('</tr></tbody>\n')


for group in allgroups:
    f.write('<tbody><tr>')

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
            f.write('<td>-</td>')
        else:
            gr = rev.groups[group]
            if (gr.n_fail>0):
                f.write('<td><span class="fail">' + str(gr.n_fail) + '</span> / ' + str(gr.n_tests) + '</td>')
            else:
                f.write('<td>' + str(gr.n_fail) + ' / ' + str(gr.n_tests) + '</td>')

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
                    f.write('<td><span class="fail" title="%s">1</span></td>'%rev.groups[group].fail_text[fail])
                else:
                    f.write('<td>-</td>')

            f.write('</tr>\n')
        f.write('</tbody>\n')
                
    f.write('\n\n')


f.write('</table>')


f.write('</body></html>')


