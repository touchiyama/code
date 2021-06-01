#!/home/touchiyama/.pyenv/versions/anaconda3-2019.03/bin/python

import sys
import subprocess
from concurrent import futures

args=sys.argv
f_name=args[1];
N_cpu=int(args[2]);

f=open(f_name,"r")
array=[]
for line in f:
  array.append(line.rstrip('\n'))
f.close()

def run(cmd):
  subprocess.call(cmd,shell=True)
  print('%s' % cmd)

future_list=[]
with futures.ThreadPoolExecutor(max_workers=N_cpu) as executer:
  future=[executer.submit(run,cmd) for cmd in array]
  future_list.append(future)
  _ = futures.as_completed(fs=future_list)

print("finished")

