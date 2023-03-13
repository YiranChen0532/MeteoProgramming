# -*- coding: utf-8 -*-
"""
Created on Tue Mar  7 19:44:39 2023

@author: CYR
"""
# Copyright (c) 2023 Yiran Chen

#VortexPara.py:
#the code is used to calculate vortex moment diagnostics on Cartesian field  
#the Equations of calculating vortex moments diagnostics can be in Matthewman et al. (2009, J. Climate) and Mithell et al. (2011, J. Atmos. Sci.) 
# Copyright (c) 2023 Yiran Chen

#input: 
#1. z = Geopotential on pressure level(usually 10hPa) or PV on isentropic surface(usually 850K) of the north-hemisphere(2D-ndarray) 
#2. lons = longitude (1D-ndarray)
#3. lats = latitude (1D-ndarray)
#4. zb = boundary contour (float) (can be calculated by the function Vortex_boundary)
#return:
#Centroid:
#1. (lon_,lat_) = the centriod of vortex in Polar coordinates (lon_ in {0,360) lat_ in (0,90))
#2. (x_,y_) = the centriod of vortex in Cartsian coordinates
#Angel_AspectRatio:
#1. psi = the angle between x-axis and major axis of ellipse (in {0,180))
#2. r = aspect ratio, the ratio of major axis to minor axis (>1, usually 1.0-2.0, see Seviour(2013, GRL))
#VortexArea:
#1. area = the area of the vortex

import numpy as np
import xarray as xr

#e_lon=1.0
#e_lat=1.0
#dy = np.deg2rad(e_lat)

def find_nearest(array,value):
    #array 1D array
    array = np.asarray(array)
    ind = np.abs(array-value).argmin()
    return ind

def Vortex_boundary(z,lats,zb_lat):
    #z: 3D numpy array(t,lat,lon)
    #given latitude, calculate Multi-year average Z10 along the given latitude
    i = find_nearest(lats, zb_lat)
    zb = z[:,i,:].mean()
    
    #return 1 float: the boundary of vortex
    return zb

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
        
        #dx = np.deg2rad(e_lon)*R
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
        
        #dx = np.deg2rad(e_lon)*R
        for j in range(M):
            xa = R*np.cos(np.deg2rad(lons[j]))-x_
            ya = R*np.sin(np.deg2rad(lons[j]))-y_
            z_2 = abs(z[i][j]-zb) if z[i][j]-zb<0 else 0
            jj += (z_2)*pow(xa,k)*pow(ya,l)#*dx*dy
    return jj

def VortexArea(z,lons,lats,zb):
    #Z: 2D numpy array (lat,lon)
    ss00 = Vortex_Geo(z, lons, lats, zb, 0, 0)
    area = ss00/zb
    return area

def Centroid(z,lons,lats,zb):
    #Z: 2D numpy array (lat,lon)
    ss00 = Vortex_Geo(z, lons, lats, zb, 0, 0)
    ss10 = Vortex_Geo(z, lons, lats, zb, 1, 0)
    ss01 = Vortex_Geo(z, lons, lats, zb, 0, 1)
    
    x_, y_ = ss10/ss00, ss01/ss00
    lon_ = np.rad2deg(np.arctan(ss01/ss10))
    if x_<0:
        lon_ = lon_+180
    elif x_>0 and y_<0:
        lon_ = lon_+360
    elif x_>0 and y_>0:
        pass
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
    
    #The centroid of Votrex at the given moment
    #(lon_,lat_,) = in polar coord; (x_,y_) = in Cartsian coord
    return (lon_,lat_), (x_,y_) 

def Angel_AspectRatio(z,lons,lats,zb, x_, y_):
    #Z: 2D numpy array (lat,lon)
    jj11 = Vortex_Geo_centralized(z, lons, lats, zb, 1, 1, x_, y_)
    jj20 = Vortex_Geo_centralized(z, lons, lats, zb, 2, 0, x_, y_)
    jj02 = Vortex_Geo_centralized(z, lons, lats, zb, 0, 2, x_, y_)
    jjdif = jj20-jj02
    psi = np.rad2deg(0.5*np.arctan(2*jj11/jjdif))
    if jjdif>0 and jj11>0:
        pass
    elif jjdif>0 and jj11<0:
        psi += 180
    elif jjdif<0:
        psi += 90

    Ja = jj20+jj02
    Jb2 = 4*pow(jj11,2)
    Jc2 = pow(jj20-jj02,2)
    r = np.sqrt(abs((Ja+np.sqrt(Jb2+Jc2))/(Ja-np.sqrt(Jb2+Jc2))))
    
    #psi = the angle between x-axis and major axis of ellipse (in {0,180))
    #r = aspect ratio, the ratio of major axis to minor axis (>1)
    return psi,r
