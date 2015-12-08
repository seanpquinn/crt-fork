import numpy as np
import scipy.stats as ss
import multiprocessing as mp
import subprocess as sp
import os
import sys

rndseed = np.random.randint(0,10**7)

n=100

l=306.59
lat=-7.45
E=50.73
eid=152691828300

FNULL = open(os.devnull,'w')

a1,b1=(13.8-18.)/2.5,(20.-18.)/2.5
a2,b2=(0.-0.2)/0.12,(np.inf-0.2)/0.12

def backtrack(counter):
  #First determine how many params to vary
  #num_params=np.random.randint(1,20,size=1)[0]
  #Draw num_params from param list without replacement
  #params_used=np.random.choice(np.arange(1,21),num_params,replace=False)
  np.random.seed(counter+rndseed)
  params=np.array([np.random.normal(0.1,1.8),
    np.random.normal(3.0,0.6),
    np.random.normal(-0.9,0.8),
    np.random.normal(-0.8,0.3),
    np.random.normal(-2.0,0.1),
    np.random.normal(-4.2,0.5),
    np.random.normal(0.0,1.8),
    np.random.normal(0.1,0.1),
    np.random.normal(0.4,0.03),
    np.random.normal(0.27,0.08),
    np.random.normal(1.4,0.1),
    np.random.normal(-1.1,0.1),
    np.random.normal(9.22,0.08),
    ss.truncnorm.rvs(a1,b1,18.,2.5,size=1),
    ss.truncnorm.rvs(a2,b2,0.2,0.12,size=1),
    np.random.normal(5.3,1.6),
    np.random.normal(4.6,0.3),
    np.random.normal(49.,1.),
    np.random.normal(4.8,0.2),
    np.random.normal(2.9,0.1)
    ])
  line1=bytearray('C %.2f %.2f %.2f 1 -1\n' %(E,l,lat),'ascii')
  line2=bytearray('F jf2012 0 0'+''.join([' %.6f']*len(params)) % tuple(params),'ascii')
  #proc=sp.Popen(['./bin/CRT','-backtrack','-rawinp'],stdin=sp.PIPE,stdout=sp.PIPE,stderr=sp.PIPE,close_fds=True)
  proc=sp.Popen(['./bin/CRT','-backtrack','-rawinp'],stdin=sp.PIPE,stdout=sp.PIPE,stderr=sp.PIPE,close_fds=True)
  proc.stdin.write(line1)
  proc.stdin.write(line2)
  out,err=proc.communicate()
  out=out.decode('ascii').split('\t')
  #print(out)
  return float(out[21]),float(out[22])

pool = mp.Pool(processes=4)
srcl,srcb=zip(*pool.map(backtrack,range(n)))
src=np.zeros((len(srcl),2))
src[:,0]=srcl
src[:,1]=srcb
np.save('evt_152691828300',src)