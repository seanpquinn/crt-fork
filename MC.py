import numpy as np
import scipy.stats as ss
import multiprocessing as mp
import subprocess as sp
import os

np.random.seed(2368687)

n=5*10**3

l=306.59
lat=-7.45
E=50.73
eid=152691828300

FNULL = open(os.devnull,'w')

def backtrack(counter):
  #First determine how many params to vary
  num_params=np.random.randint(1,20,size=1)[0]
  #Draw num_params from param list without replacement
  params_used=np.random.choice(np.arange(1,21),num_params,replace=False)
  p1=np.random.normal(0.1,1.8)
  p2=np.random.normal(3.0,0.6)
  p3=np.random.normal(-0.9,0.8)
  p4=np.random.normal(-0.8,0.3)
  p5=np.random.normal(-2.0,0.1)
  p6=np.random.normal(-4.2,0.5)
  p7=np.random.normal(0.0,1.8)
  p9=np.random.normal(0.1,0.1)
  p10=np.random.normal(0.4,0.03)
  p11=np.random.normal(0.27,0.08)
  p12=np.random.normal(1.4,0.1)
  p13=np.random.normal(-1.1,0.1)
  p14=np.random.normal(9.22,0.08)
  a,b=(13.8-18.)/2.5,(20.-18.)/2.5
  p15=ss.truncnorm.rvs(a,b,18.,2.5,size=1)
  a,b=(0.-0.2)/0.12,(np.inf-0.2)/0.12
  p16=ss.truncnorm.rvs(a,b,0.2,0.12,size=1)
  p17=np.random.normal(5.3,1.6)
  p18=np.random.normal(4.6,0.3)
  p19=np.random.normal(49.,1.)
  p20=np.random.normal(4.8,0.2)
  p21=np.random.normal(2.9,0.1)
  line1='C %.2f %.2f %.2f 1 -1\n' %(E,l,lat)
  line2='F jf2012 0 0 %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f' %(p1,p2,p3,p4,p5,p6,p7,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18,p19,p20,p21)
  f2=open('mc_infile_%i'%counter,'w')
  f2.write(line1+line2)
  f2.close()
  res=sp.Popen(['./bin/CRT','-backtrack','infile=./mc_infile_%i'%counter],stdout=sp.PIPE,stderr=sp.PIPE,close_fds=True)
  out,err=res.communicate()
  out=out.decode('ascii').split('\t')
  sp.call(['rm','mc_infile_%i'%counter]) #Clean up
  counter+=1
  return float(out[21]),float(out[22])

pool = mp.Pool(processes=4)
srcl,srcb=zip(*pool.map(backtrack,range(n)))
src=np.zeros((len(srcl),2))
src[:,0]=srcl
src[:,1]=srcb
np.save('evt_152691828300',src)