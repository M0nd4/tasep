#!/bin/bash

begin=$1
end=$2
for i in $(seq $begin $end); do
  epoch=$((2**$i))
  echo "RUNNING WITH epoch = $epoch"
  ./tasep_playground --rates ./data/rates6.txt  --epoch $epoch -m
done

