#!/bin/bash

if [ ! -f index.html ]
then
  echo "Please run this script in the doc output directory (install/doc)"
  exit 1
fi

mkdir images >/dev/null

echo "Downloading images (press ctrl-c to cancel) ..."
cd images
{
trap "echo \"(skipping)\"" SIGINT
wget -q -nd -A png,gif -m -l 1 -np  http://www.dealii.org/images/steps/developer/
}
cd ..

echo "Patching html files ..."
sed -i 's#"http://www.dealii.org/images/steps/developer/\(step-.*\)"#"images/\1"#g' step_*.html

echo "all done!"