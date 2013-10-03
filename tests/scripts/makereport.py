# accepts 0,1, or 2 arguments. If a string starting with a number is handed in, it is assumed to be a subdirectory of the current directory to run on. If not specified, the newest build is used. Any other string is taken as the branch name for this test (or treated as mainline). Order of the arguments does not matter.
# for questions: Timo Heister

import xml.etree.ElementTree as ET
import glob
import sys
from datetime import datetime
import subprocess

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

branch=''
args=sys.argv
args.pop(0)

dirname=""
while len(args)>0:
    if args[0].startswith("20"): #I hope this script is not used in the year 2100
        dirname=args[0].replace('/','')
    else:
        branch=args[0].replace('/','')+'/'
    args.pop(0)

if len(glob.glob(dirname+'/Update.xml'))>0:
    #new format
    tree = ET.parse(dirname+'/Update.xml')
    name = tree.getroot().find('BuildName').text
    number = tree.getroot().find('Revision').text
    date = datetime.strptime(dirname,"%Y%m%d-%H%M")
else:
    #old format
    tree = ET.parse(dirname+'/Notes.xml')
    name = tree.getroot().attrib['BuildName']
    number = name.split('-')[-1]
    number = number[1:]
    date = datetime.strptime(dirname,"%Y%m%d-%H%M")

header = "Revision: %s"%number + "\n"
header += "Date: %s"%(date.strftime("%Y %j  %F  %U-%w")) + '\n'
id = subprocess.check_output(["id","-un"])+'@'+subprocess.check_output(["hostname"])
id=id.replace('\n','')
header += "Id:  %s"%id

#now Test.xml:
tree = ET.parse(dirname+'/Test.xml')
root = tree.getroot()
testing = root.find('Testing')

tests={}

for test in testing.findall("Test"):
    status = test.attrib['Status']
    fail=False
    if status=="failed": fail=True
    name = test.find('Name').text
    group = name.split('/')[0]

    if fail:
        line = "%s  3   %s%s"%(date,branch,name)
    else:
        line = "%s   +  %s%s"%(date,branch,name)

    if group not in tests: tests[group]=[]
    tests[group].append( line )

for g in sorted(tests):
    group = tests[g]
    print header
    for l in group:
        print l



