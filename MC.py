import numpy as np
import scipy.stats as ss
import multiprocessing as mp
import subprocess as sp
import datetime as dt
import os
import sys

n=70

#l=348.63
#lat=-21.11
#E=113.06
#eid=152364323000

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
  line1=bytearray('C %.2f %.2f %.2f 1 -1\n' %(E,lon,lat),'ascii')
  line2=bytearray('F jf2012 0 0'+''.join([' %.6f']*len(params)) % tuple(params),'ascii')
  #proc=sp.Popen(['./bin/CRT','-backtrack','-rawinp'],stdin=sp.PIPE,stdout=sp.PIPE,stderr=sp.PIPE,close_fds=True)
  proc=sp.Popen(['./bin/CRT','-backtrack','-rawinp'],stdin=sp.PIPE,stdout=sp.PIPE,stderr=sp.PIPE,close_fds=True)
  proc.stdin.write(line1)
  proc.stdin.write(line2)
  out,err=proc.communicate()
  out=out.decode('ascii').split('\t')
  #print(out)
  return float(out[21]),float(out[22])

idarr,Earr,lonarr,latarr=np.loadtxt('data/events50_parsed.txt',
  dtype={'names':('ID','E','lon','lat'),
  'formats':('int','float','float','float')},skiprows=7,unpack=True)

nevents=len(idarr)
begin=dt.datetime.now()
beginstr=dt.datetime.now().strftime("%c")
print("Starting program: %s" %beginstr)
for j in range(nevents):
  #Globals rndseed, E, lon, and lat must be specified for each loop iteration
  rndseed=np.random.randint(0,9*10**7)
  E=Earr[j]
  lon=lonarr[j]
  lat=latarr[j]
  eid=idarr[j]
  start=dt.datetime.now()
  startstr=dt.datetime.now().strftime("%c")
  print("\nStart sampling for Event %i: %s" %(eid,startstr))
  pool = mp.Pool(processes=4)
  srcl,srcb=zip(*pool.map(backtrack,range(n)))
  src=np.zeros((len(srcl),2))
  src[:,0]=srcl
  src[:,1]=srcb
  np.save('evt_%i' %eid,src)
  pool.close()
  stop=dt.datetime.now()
  stopstr=dt.datetime.now().strftime("%c")
  deltat=stop-start
  hours,mins=divmod(deltat.seconds,3600)
  mins=int(mins/60)
  print("Finished: %s  Elapsed time: %02i hours %02i minutes" %(stopstr,hours,mins))
  
stop=dt.datetime.now()
stopstr=dt.datetime.now().strftime("%c")
deltat=stop-begin
days=deltat.days
hours,mins=divmod(deltat.seconds,3600)
mins=int(mins/60)
print("\n\nEnding program: %s Elapsed time: %02i:%02i:%02i" %(stopstr,days,hours,mins))