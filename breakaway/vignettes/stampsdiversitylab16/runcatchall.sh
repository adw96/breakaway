#!/bin/bash
# A script to run CatchAll on every csv file in (only) the current directory

mylist=$(find . -maxdepth 1 -name "*.csv" -type f) 
for i in $mylist; do
  j=${i%.csv};
  mono ./CatchAllCmdL.exe $i $j 1
done
