#!/bin/bash

CMD=./sf_ev3a
H=64
T=1
D=20
F=pos1000.fen

S=0
while true;
do
	MAX=$S
	MIN=$((-$S))
	X1=$MIN
	while [ $X1 -le $MAX ];
	do
		X2=$MIN
		while [ $X2 -le $MAX ];
		do
			X3=$MIN
			while [ $X3 -le $MAX ];
			do
				if [[ $X1 == $MIN || $X1 == $MAX || $X2 == $MIN || $X2 == $MAX || $X3 == $MIN || $X3 == $MAX ]]
				then
					echo -n $X1 $X2 $X3
					if [ -f $CMD ];
					then
						$CMD bench $X1 $X2 $X3 $H $T $D $F 2>&1 | grep '^Nodes searched'| cut -d':' -f2
					else
						exit 1;
					fi
				fi
				X3=$(($X3+1))
			done
			X2=$(($X2+1))
		done
		X1=$(($X1+1))
	done
	S=$((S+1))
done
