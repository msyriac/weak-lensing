#!/bin/bash
for x in {1..350}
do
     nohup wq sub -c "source ~/.bashrc;source ~/.bash_profile;python2.7 toy.py -n 1000 -p -e 0.05 -s 0.3 -t 1.49e-06 -i $x" >/dev/null &
     sleep 0.2
done
#EOF
