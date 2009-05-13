#!/bin/bash
for ((i=1;i<=10;i+=1)); do 
let nodes=$i+1; 

for ((run=1;run<=10;run+=1)); do
#echo "hops: $i run: $run";
ns mesh.tcl \
-random-nodeid off \
-run $run \
-out out-tmp-$run \
-topology chain \
-n $nodes \
-trfsrc 0 \
-trfdst $i \
-trftype cbr \
-trfnsrc 1 \
-trfprio 1 \
-cbr-rate 12000000 \
-cbr-pkt 500 \
-prfall 2 \
-cbr-rnd 1 \
-duration 15 \
-warm 5 \
> /dev/null 2>&1; 
done

cat out-tmp-* > out-tptperhop-$i;
rm out-tmp-*;
echo "";
echo "hops: $i"; 
recover out-tptperhop-$i | grep -A 2 e2e_tpt;
done
