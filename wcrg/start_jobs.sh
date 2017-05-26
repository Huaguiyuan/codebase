#!/bin/bash
ssh -i chrisstartscloud.pem ubuntu@$1 -t 'nohup mpirun -np 16 --host host0,host1 ./frg > out 2 > out.err < /dev/null &'
#ssh -i chrisstartscloud.pem ubuntu@$1 -t 'nohup ./copy_res.sh 1000 < /dev/null > outshell 2>&1 &'

