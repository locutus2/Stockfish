#!/bin/bash

x=-56
while [ $x -le 56 ];
do
    ./sf bench $x 16 1 13 pos1000.fen 2>&1 | grep 'Nodes searched' | cut -d':' -f2|sed 's/ //g'
    x=$(($x+4))
done
