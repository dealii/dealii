#!/bin/bash
cp *png ~/public_html/bench/
chmod a+r ~/public_html/bench/*png
#rsync -a *png root@dealii.org:/home/archiver/public_html/
#ssh root@dealii.org chmod a+r /home/archiver/public_html/*png
