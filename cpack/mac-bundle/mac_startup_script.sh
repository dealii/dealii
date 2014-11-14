#!/bin/sh

function relocate_if_necessary {
    # If necessary, we relocate the libraries and binaries to the new place, 
    # by calling install_name_tool on each library in the installation tree
    OLD=`cat $1/Contents/Resources/etc/dealii.location`
    if [ ! "$OLD"=="$1" ];then
	echo $1 > $1/Contents/Resources/etc/dealii.location
	echo Needs relocation from $OLD to $1!
    fi  
    return
}

if [ "$BASH_SOURCE" == "$0" ]
then
  export DEAL_II_BUNDLE=`echo "$0" | sed -e 's|/Contents/MacOS/.*||'`
  export DEAL_II_RESOURCES=$DEAL_II_BUNDLE/Contents/Resources
  open -a /Applications/Utilities/Terminal.app $DEAL_II_RESOURCES/bin/dealii-terminal
else
  export DEAL_II_BUNDLE=`echo "$BASH_SOURCE" | sed -e 's|/Contents/MacOS/.*||'`
  export DEAL_II_RESOURCES=$DEAL_II_BUNDLE/Contents/Resources
  source $DEAL_II_RESOURCES/etc/dealii.conf
fi
