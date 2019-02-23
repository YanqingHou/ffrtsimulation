#!/bin/bash
if qstat -r | grep -q R
then
printf "running \a"
fi
