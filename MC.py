import numpy as np
import scipy.stats as ss
import multiprocessing as mp
import subprocess as sp
import datetime as dt
import os
import sys

n=2*10**5 #num samples

#Truncated Gaussian interval params
a1,b1=(12.6-18.25)/2.75,(19.5-18.25)/2.75
a2,b2=(0.-0.2)/0.12,(1-0.2)/0.12

def backtrack(counter):
  """Performs pRNG, calls CRT, then saves results in memory"""
  np.random.seed(counter+rndseed) #Start seed is randomly assigned
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
    ss.truncnorm.rvs(a1,b1,18.25,2.75,size=1),
    ss.truncnorm.rvs(a2,b2,0.2,0.12,size=1),
    np.random.normal(5.3,1.6),
    np.random.normal(4.6,0.3),
    np.random.normal(49.,1.),
    np.random.normal(4.8,0.2),
    np.random.normal(2.9,0.1)
    ])
  #Text data that is usually read from a file is instead
  #converted to a byte array and piped directly to stdin
  line1=bytearray('C %.2f %.2f %.2f 1 -1\n' %(E,lon,lat),'ascii')
  line2=bytearray('F jf2012 0 0'+''.join([' %.6f']*len(params)) % tuple(params),'ascii')
  #Call CRT using timeout to avoid spurious hangs
  proc=sp.Popen(['timeout','3','./bin/CRT','-backtrack','-rawinp'],stdin=sp.PIPE,stdout=sp.PIPE,stderr=sp.PIPE,close_fds=True)
  proc.stdin.write(line1)
  proc.stdin.write(line2)
  #Save text output to list in memory
  out,err=proc.communicate()
  out=out.decode('ascii').split('\t')
  #Try except to record output long, lat, timeout flag
  try:
    return float(out[21]),float(out[22]),int(out[-1][0])
  except:
    print("Timeout condition")
    return 0,0,1

#Pre-processed Herald data file
idarr,Earr,lonarr,latarr=np.loadtxt('data/events50_parsed.txt',
  dtype={'names':('ID','E','lon','lat'),
  'formats':('int','float','float','float')},skiprows=166,unpack=True)

nevents=len(idarr)
#Record date info for run duration
begin=dt.datetime.now()
beginstr=dt.datetime.now().strftime("%c")
print("Starting program: %s" %beginstr)
#Main loop
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
  #Take advantage of multiple cores. My machine has 4
  pool = mp.Pool(processes=4)
  srcl,srcb,status=zip(*pool.map(backtrack,range(n)))
  src=np.zeros((len(srcl),3))
  src[:,0]=srcl
  src[:,1]=srcb
  src[:,2]=status
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