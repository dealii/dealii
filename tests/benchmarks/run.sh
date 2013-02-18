#!/bin/bash


while :
do
   ./bench.sh

   ./doplots.sh
   
   ./sync.sh

   sleep 3600
done