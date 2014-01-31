#!/bin/bash
for x in {1..80}
do
     nohup wq sub -r "job_name:mat_fisher; priority:low" -c "source ~/.bashrc;source ~/.bash_profile;python2.7 cluster_fisher.py -v 3 -n 500000 -t 200 -i $x" >/dev/null &
     sleep 0.2
done
#EOF
