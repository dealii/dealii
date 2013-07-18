#!/bin/bash

FILE=$1

FIRST_YEAR=$(grep Copyright $FILE | perl -pi -e 's/.*?(\d{4}).*/\1/')
if [ "$FIRST_YEAR" = "2013" ] || (! [[ "$FIRST_YEAR" =~ "^[0-9]+$" ]])
then
  YEAR="2013"
else if
  YEAR="$FIRST_YEAR - 2013"
fi

sed -e "s/@YEAR@/$YEAR/" doc/license/header-template > /tmp/foobar22.temp

header_lines=1
on_whitespace=true
whitespace_detected=false
while read line
do
  if $on_whitespace && [[ $line =~ "^$" ]]
  then
    whitespace_detected=true
    header_lines=$[$header_lines+1]
    continue
  else
    on_whitespace=false
  fi

  if ! [[ $line =~ "^//.*$" ]]
  then
    break
  else
    header_lines=$[$header_lines+1]
  fi
done < $FILE

if $whitespace_detected
then
  echo "$FILE --> whitespace"
  echo "$FILE --> $header_lines"
else
  if [ $header_lines -lt 10 ]
  then
    echo "$FILE --> $header_lines"
  fi
fi

tail -n+$header_lines $FILE >> /tmp/foobar22.temp

cp /tmp/foobar22.temp $FILE

