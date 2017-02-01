#!/bin/bash

for i in {12..40}
do
    root -l -b -q "analysis.C+(${i},30000,100000,0,30)" | tee plots/log_${i}.log
done
