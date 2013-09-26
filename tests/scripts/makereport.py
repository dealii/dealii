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



if len(sys.argv)==2:
    dirname=sys.argv[1].replace('/','')
else:
    n=glob.glob("*/Build.xml")
    n.sort(reverse=True)
    dirname = n[0].replace('/Build.xml','')

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

print "Revision: %s"%number
print "Date: %s"%(date.strftime("%Y %j  %F  %U-%w"))
#print "Id:  %s"%name
id = subprocess.check_output(["id","-un"])+'@'+subprocess.check_output(["hostname"])
id=id.replace('\n','')
print "Id:  %s"%id

#now Test.xml:
tree = ET.parse(dirname+'/Test.xml')
root = tree.getroot()
testing = root.find('Testing')

for test in testing.findall("Test"):
    status = test.attrib['Status']
    fail=False
    if status=="failed": fail=True
    name = test.find('Name').text

    if fail:
        print "%s  3   %s"%(date,name)
    else:
        print "%s   +  %s"%(date,name)


