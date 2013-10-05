#!/bin/bash
for x in {1..100}
do
     nohup wq sub -c "source ~/.bashrc;source ~/.bash_profile;python2.7 toymodel_mod.py -n 200 -1 -0.01 -2 0.02 -e 0.05 -s 0.3 -t 1.49e-03 -i $x" >/dev/null &
done
#EOF
