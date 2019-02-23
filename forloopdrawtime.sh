#!/bin/bash
if [ $# -eq 0 ]; then
let st=1
let N=100
else
let st=$1
let N=$2
fi

for (( i=$st; i <= $N; ++i ))
do
	if [ -d "Job$i" ]; then
		./resdraw.sh $i
	else
		echo "Job$i non-exist!"
	fi
done
