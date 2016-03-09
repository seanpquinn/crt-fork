import numpy as np
import matplotlib.pyplot as plt
import os
import itertools
import scipy.stats as ss
import scipy.optimize as so
import numpy.linalg as linalg
import scipy.spatial.distance as ssd

#Load pre-processed data
darray=np.genfromtxt('events50_parsed.txt',dtype=None,names=('id','E','lon','lat'))

#Script excepts .npy arrival direction data in
#it's own directory. List of file names is generated using os.listdir
filelist=[]
datadir='./raw_arrival/'
dirlist=os.listdir(datadir)
for s in dirlist:
  if ('evt' in s) and ('.npy' in s):#Grab only event data files
    filelist.append(datadir+s)

#Compute area of 2D polygon
#Reference: http://stackoverflow.com/questions/19873596/convex-hull-area-in-python
def PolyArea2D(pts):
    lines = np.hstack([pts,np.roll(pts,-1,axis=0)])
    area = 0.5*abs(sum(x1*y2-x2*y1 for x1,y1,x2,y2 in lines))
    return area

def onesig(x):
  """Numerically solve for one sigma contour"""
  return f_kde_flat[f_kde_flat>x].sum()-0.68*N_const

def twosig(x):
  """Numerically solve for two sigma contour"""
  return f_kde_flat[f_kde_flat>x].sum()-0.95*N_const

def threesig(x):
  """Numerically solve for three sigma contour"""
  return f_kde_flat[f_kde_flat>x].sum()-0.997*N_const

def bbox(points):
  a=np.zeros((2,2))
  a[:,0]=np.min(points,axis=0)
  a[:,1]=np.max(points,axis=0)
  return a

def fitEllipse(x,y):
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    E, V =  linalg.eig(np.dot(linalg.inv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:,n]
    return a

def ellipse_center(a):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return np.array([x0,y0])

def ellipse_angle_of_rotation( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    return 0.5*np.arctan(2*b/(a-c))

def ellipse_axis_length( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=np.sqrt(up/down1)
    res2=np.sqrt(up/down2)
    return np.array([res1, res2])

def find_ellipse(x, y):
    xmean = x.mean()
    ymean = y.mean()
    x -= xmean
    y -= ymean
    a = fitEllipse(x,y)
    center = ellipse_center(a)
    center[0] += xmean
    center[1] += ymean
    phi = ellipse_angle_of_rotation(a)
    axes = ellipse_axis_length(a)
    x += xmean
    y += ymean
    return center, phi, axes

pfile_exist=False
if 'analysis_parameters.txt' in os.listdir('./'):
  pfile_exist=True

paramfile=open('analysis_parameters.txt','ab')
#Write header if opening new file only
if pfile_exist==False:
  paramfile.write(bytes('#ID E meanl meanb ecc num_modes semimajor semiminor lmincut lmaxcut bmincut bmaxcut\n'))
  
for evt in filelist:
  rawdata=np.load(evt)
  rawdata=rawdata[:,:2]
  evtnum=int(''.join([k for k in evt if k.isdigit()]))
  #pos_ang=np.load('./pos_ang/pos_ang_%i.npy' %evtnum)
  obs_ind=np.where(darray['id']==evtnum)[0][0]
  obs_E=darray['E'][obs_ind]
  #max_pos_ang=0.5*0.06*1e4/(10.*obs_E)
  #outlier_index=np.where(pos_ang>max_pos_ang)[0]
  dvarinv=linalg.inv(np.cov(rawdata,rowvar=0))
  mahal=np.zeros(len(rawdata))
  sampmean=np.array([rawdata[:,0].mean(),rawdata[:,1].mean()])
  plt.subplot(211)
  plt.scatter(rawdata[:,0],rawdata[:,1],s=3)
  for i in range(len(rawdata)):
    mahal[i]=ssd.mahalanobis(np.array([rawdata[i,0],rawdata[i,1]]),sampmean,dvarinv)
  outlier_index=np.where(mahal<mahal.mean()+3*mahal.std())[0]
  plt.subplot(212)
  plt.scatter(rawdata[outlier_index,0],rawdata[outlier_index,1],marker='.',s=3)
  plt.savefig('mahal_cut_%i.png' %evtnum)
  plt.clf()
  data=rawdata[outlier_index]
  kernel=ss.gaussian_kde(data.T)
  xmin,xmax,ymin,ymax=data[:,0].min(),data[:,0].max(),data[:,1].min(),data[:,1].max()
  # Contour plotting routine from 
  # http://stackoverflow.com/questions/30145957/plotting-2d-kernel-density-estimation-with-python
  xx,yy=np.mgrid[xmin:xmax:150j,ymin:ymax:150j]
  pos = np.vstack([xx.ravel(), yy.ravel()])
  f_kde=np.reshape(kernel(pos).T, xx.shape) # Expensive
  #Remove numerical artifacts
  f_kde=np.ma.array(f_kde,mask=f_kde<1e-5)
  f_kde_flat=f_kde.flatten()
  f_kde_flat=f_kde_flat[~f_kde_flat.mask].data
  N_const=f_kde_flat.sum()
  cset=plt.contourf(xx,yy,f_kde,cmap='Blues')
  x1=cset.levels[2]
  x2=x1/10.
  x3=x1/100.
  l1=so.fsolve(onesig,x1)[0]
  l2=so.fsolve(twosig,x2)[0]
  l3=so.fsolve(threesig,x3)[0]
  plt.figure(1,figsize=(8,6))
  cmap=plt.cm.OrRd
  cmap.set_under("w")
  plt.contourf(xx,yy,f_kde,cmap=cmap,fignum=1)
  plt.colorbar()
  cset=plt.contour(xx,yy,f_kde,[l1,l2,l3],colors='k',fignum=1)
  ell_size=[]
  ps=cset.collections[2].get_paths()
  for arr in ps:
    ell_size.append(len(arr))
  ps=cset.collections[2].get_paths()[np.argmax(ell_size)]
  vert=ps.vertices
  np.save('./ellipse_coords/ellipse_data_%i' %evtnum,vert)
  xsig3=vert[:,0]
  ysig3=vert[:,1]
  center,phi,axes=find_ellipse(xsig3,ysig3)
  ecc=np.sqrt(1-(np.min(axes)/np.max(axes))**2)
  ellarea=PolyArea2D(vert)
  fmt={}
  strs=['68%','95%','99.7%']
  for l, s in zip(cset.levels,strs):
    fmt[l]=s
  plt.clabel(cset,fmt=fmt)
  plt.title(r'ID=%i $\bar{\ell}=$%.2f $\bar{b}=$%.2f $e=$%.3f' %(evtnum,center[0],
      center[1],ecc))
  fig = plt.gcf()
  fig.savefig("lat_lon_contour_%i.png" %evtnum)
  plt.clf()
  paramsout=np.array([evtnum,obs_E,center[0],center[1],ecc,
    np.max(axes),np.min(axes),xmin,xmax,ymin,ymax])
  paramsout=np.column_stack(paramsout)
  np.savetxt(paramfile,paramsout,fmt="%i %.2f %.4f %.4f %.5f %.4f %.4f %i %i %i %i")

paramfile.close()
