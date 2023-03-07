# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 16:25:22 2023

@author: CYR
"""

import numpy as np
import xarray as xr
import os
from tools import eof_cal
from defmap2 import def_map, def_fig

root = r"E:\atmdata\workdata"+'\\'
path = r'E:\python_code\atm\work'
os.chdir(path)

ds = xr.open_dataset(root+'ERA5-1981-2022.nc')
ds = ds.sel(expver=1)
#10hPa EOF
z_arr = ds.z.isel(time=slice(0,-31))
z = z_arr.sel(latitude=slice(90,20),level=10)
z1 = z.isel(time=285)
lons = z1.longitude
lats = z1.latitude

e_lon=1.0
e_lat=1.0
dy = np.deg2rad(e_lat)


def find_nearest(array,value):
    #array 1D array
    array = np.asarray(array)
    ind = np.abs(array-value).argmin()
    return ind

def Vortex_boundary(z,lats,zb_lat):
    #z: 3D numpy array(t,lat,lon)
    #given latitude, calculate Multi-year average Z10 along the given latitude
    #return 1 float
    i = find_nearest(lats, zb_lat)
    zb = z[:,i,:].mean()
    return zb.values
    

def Vortex_Geo(z,lons,lats,zb,k,l):
    #Z: 2D numpy array (lat,lon)
    ss = 0
    N, M = z.shape
    for i in range(N):
        #method1
        #R = np.deg2rad(90-lats[i])
        #method2
        #R = np.sqrt(2-2*abs(np.sin(np.deg2rad(lats[i]))))
        #method3
        R = np.cos(np.deg2rad(lats[i]))/(1+np.sin(np.deg2rad(lats[i])))
        
        dx = np.deg2rad(e_lon)*R
        for j in range(M):
            x = R*np.cos(np.deg2rad(lons[j]))
            y = R*np.sin(np.deg2rad(lons[j]))
            z_2 = abs(z[i][j]-zb) if z[i][j]-zb<0 else 0
            ss += (z_2)*pow(x,k)*pow(y,l)#*dx*dy
    return ss

def Vortex_Geo_centralized(z,lons,lats,zb,k,l,x_,y_):
    #Z: 2D numpy array (lat,lon)
    jj = 0
    N, M = z.shape
    for i in range(N):
        #method1
        #R = np.deg2rad(90-lats[i])
        #method2
        #R = np.sqrt(2-2*abs(np.sin(np.deg2rad(lats[i]))))
        #method3
        R = np.cos(np.deg2rad(lats[i]))/(1+np.sin(np.deg2rad(lats[i])))
        
        dx = np.deg2rad(e_lon)*R
        for j in range(M):
            xa = R*np.cos(np.deg2rad(lons[j]))-x_
            ya = R*np.sin(np.deg2rad(lons[j]))-y_
            z_2 = abs(z[i][j]-zb) if z[i][j]-zb<0 else 0
            jj += (z_2)*pow(xa,k)*pow(ya,l)#*dx*dy
    return jj

def Centroid(z,lons,lats,zb):
    ss00 = Vortex_Geo(z.values, lons.values, lats.values, 
                      zb, 0, 0)
    ss10 = Vortex_Geo(z.values, lons.values, lats.values, 
                      zb, 1, 0)
    ss01 = Vortex_Geo(z.values, lons.values, lats.values, 
                      zb, 0, 1)
    
    x_, y_ = ss10/ss00, ss01/ss00
    lon_ = np.rad2deg(np.arctan(ss01/ss10))
    if x_<0:
        lon_ = lon_+180
    #method1
    #R = x_/np.cos(np.deg2rad(lon_))
    #lat_ = 90-np.rad2deg(R)
    #method2
    #R2 = pow(x_,2) + pow(y_,2)
    #lat_ = np.rad2deg(np.arcsin(1-R2/2))
    #method3
    R2 = pow(x_,2) + pow(y_,2)
    sinlat_ = (1-R2)/(1+R2)
    lat_ = np.rad2deg(np.arcsin(sinlat_))
    
    return (lon_,lat_), (x_,y_)

zb = Vortex_boundary(z,lats,60)
(lon_,lat_),(x_,y_) = Centroid(z1,lons,lats,zb)
print(lon_,lat_)

jj11 = Vortex_Geo_centralized(z1.values, lons.values, lats.values,
                              zb, 1, 1, x_, y_)
jj20 = Vortex_Geo_centralized(z1.values, lons.values, lats.values,
                              zb, 2, 0, x_, y_)
jj02 = Vortex_Geo_centralized(z1.values, lons.values, lats.values,
                              zb, 0, 2, x_, y_)  
psi = np.rad2deg(0.5*np.arctan(2*jj11/(jj20-jj02)))

Ja = jj20+jj02
Jb2 = 4*pow(jj11,2)
Jc2 = pow(jj20-jj02,2)
r = np.sqrt(abs((Ja+np.sqrt(Jb2+Jc2))/(Ja-np.sqrt(Jb2+Jc2))))
print(psi,r)












