#!/bin/bash

while read -r line
do
  without_version=${line%%.*}
  url="http://www.uniprot.org/uniprot/?sort=score&desc=&query=database:(type:embl%20$without_version)&fil=&force=no&format=list"
  curl $url 2> /dev/null
done
