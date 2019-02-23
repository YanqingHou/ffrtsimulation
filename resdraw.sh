#!/bin/bash
module load python/3.5.1 gnuplot/4.6.4
if [ ! -n "$1" ]
then
   cp ./DiaryLogs/Taskdata1.txt  ./DiaryLogs/Task1.diary.txt
   echo "Using previous file!"
else
   cp /scratch/tb2n14/hyq/Job$1/Task1.diary.txt ./DiaryLogs/Taskdata${1}.txt 
   cp /scratch/tb2n14/hyq/Job$1/Task1.diary.txt ./DiaryLogs/Task1.diary.txt 
fi
python3 drawprog.py

if test -e "Taskdata.txt"
then
cat Taskdata.txt | gnuplot -e "set terminal dumb; plot '<cat'  title 'Time for each block'"
rm  Taskdata.txt 
fi
#cat Taskdata.txt | gnuplot -e "set terminal dumb; set style data histogram; set style histogram cluster gap 1; plot '<cat'  title 'Time hist'"
#cat Taskdata.txt | gnuplot -e "set terminal dumb; set style data histograms; plot '<cat'  title 'Time hist'"
#cat Taskdata.txt | gnuplot -e "set terminal dumb; plot '<cat' smooth freq with boxes title 'Time hist'"
