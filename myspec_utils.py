import numpy as np
from scipy.interpolate import splev,splrep,splint


def get_dx(x):
  dx1 = x[1:]-x[:-1]
  dx = np.hstack([dx1[0],dx1[1:]+dx1[:-1],dx1[-1]])
  return dx

def average_nbins(nbin,x,y):
  if nbin==1:
    return x,y
  else:
    nfin = len(x)//nbin*nbin
    dx = get_dx(x)
    xx = np.sum([(x*dx)[0+ii:nfin:nbin] for ii in range(nbin)],axis=0)/\
      np.sum([dx[0+ii:nfin:nbin] for ii in range(nbin)],axis=0)
    yy = np.sum([(y*dx)[0+ii:nfin:nbin] for ii in range(nbin)],axis=0)/\
      np.sum([dx[0+ii:nfin:nbin] for ii in range(nbin)],axis=0)
    return xx,yy

def rebin(x,y,xnew,conserve_count=True): 
  dx = get_dx(x) 
  dxnew = get_dx(xnew) 
  if conserve_count: # count conserved (input is per pix)
      spl = splrep(x,y/dx,k=1,task=0,s=0) 
      return np.array([splint(xn-0.5*dxn,xn+0.5*dxn,spl) \
        for xn,dxn in zip(xnew,dxnew)]) 
  else: #flux conserved (input is in physical unit)
      spl = splrep(x,y,k=1,task=0,s=0) 
      return np.array([splint(xn-0.5*dxn,xn+0.5*dxn,spl)/dxn \
        for xn,dxn in zip(xnew,dxnew)]) 

def x_sorted(xx,yy):
  argxx = np.argsort(xx)
  try:
    len(yy[0])
    return (xx[argxx],)+tuple([y[argxx] for y in yy])
  except:
    return xx[argxx],yy[argxx]