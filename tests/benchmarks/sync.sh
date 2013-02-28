#!/bin/bash
echo "copying files to public_html/"
cp *png ~/public_html/bench/
cp index.html ~/public_html/bench/
chmod a+r ~/public_html/bench/*png ~/public_html/bench/*html

