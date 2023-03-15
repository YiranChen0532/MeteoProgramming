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

#usage:
#1.get zb using Vortex_boundary()
#2.get xy_map, dR_ls, Rdtheta_ls using Cart_coord()
#3.get lon_, lat_, x_, y_ using Centroid(); get area using VortexArea()
#4.get psi, r using Angel_AspectRatio()

#return:
#Centroid:
#1. (lon_,lat_) = the centriod of vortex in Polar coordinates (lon_ in {0,360) lat_ in (0,90))
#2. (x_,y_) = the centriod of vortex in Cartsian coordinates
#Angel_AspectRatio:
#1. psi = the angle between x-axis and major axis of ellipse (in {0,180))
#2. r = aspect ratio, the ratio of major axis to minor axis (>1, usually 1.0-2.0, see Seviour(2013, GRL))
#VortexArea:
#1. area = the area of the vortex
#Major_Minor:(optional)
#1. major = major axis of ellipse(2a)
#2. minor = minor axis of eliipse(2b)

import numpy as np
import xarray as xr

e_lon=1.0 #resolution of longitude/latitude (deg)
e_lat=1.0

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

def Cart_coord(lons,lats,e_lon,e_lat):
    #lons, lats: 1D numpy array
    nlon = len(lons)
    nlat = len(lats)
    xy_map = np.zeros((nlat,nlon,2))
    dR_ls = np.zeros((nlat))
    Rdtheta_ls = np.zeros((nlat))
    for i in range(nlat):
        #method1
        #R = np.deg2rad(90-lats[i])
        #method2 lambert equal-area, Matthewman et al., 2009 J.Climate
        R = np.sqrt(2-2*abs(np.sin(np.deg2rad(lats[i]))))
        dR_ls[i] = (np.cos(np.deg2rad(lats[i]))*np.deg2rad(e_lat))/R if R>0.0000005 else 0
        #method3 huang,2019 lanzhou Univ.
        #R = np.cos(np.deg2rad(lats[i]))/(1+np.sin(np.deg2rad(lats[i])))
        #dR_ls[i] = 1/(1+np.sin(np.deg2rad(lats[i])))*np.deg2rad(e_lat)
        ########################################################################
        Rdtheta_ls[i] = np.deg2rad(e_lon)*R
        for j in range(nlon):
            xy_map[i][j][0] = R*np.cos(np.deg2rad(lons[j])) #x
            xy_map[i][j][1] = R*np.sin(np.deg2rad(lons[j])) #y
    
    #return: 
    #xy_map: 3D numpy array(nlon,nlat,2), Cartesian coord for every points in Polor coord
    #As dxdy in Cartesian coord is converted to dR*Rdθ,
    #dR:1D numpy array(nlat), derivative in the radial direction
    #Rdtheta: 1D numpy array(nlat), derivative in the normal direction
    return xy_map, dR_ls, Rdtheta_ls

def Vortex_Geo(z,zb,xy_map,dR_ls,Rdtheta_ls,k,l):
    #Z: 2D numpy array (lat,lon)
    ss = 0
    nlat, nlon = z.shape
    for i in range(nlat):
        dR = dR_ls[i]
        Rdtheta = Rdtheta_ls[i]
        for j in range(nlon):
            x = xy_map[i][j][0]
            y = xy_map[i][j][1]
            z_2 = abs(z[i][j]-zb) if z[i][j]-zb<0 else 0
            ss += (z_2)*pow(x,k)*pow(y,l)*dR*Rdtheta
    return ss

def Vortex_Geo_centralized(z,zb,xy_map,dR_ls,Rdtheta_ls,k,l,x_,y_):
    #Z: 2D numpy array (lat,lon)
    jj = 0
    nlat, nlon = z.shape
    for i in range(nlat):
        dR = dR_ls[i]
        Rdtheta = Rdtheta_ls[i]
        for j in range(nlon):
            x = xy_map[i][j][0]
            y = xy_map[i][j][1]
            xa = x-x_
            ya = y-y_
            z_2 = abs(z[i][j]-zb) if z[i][j]-zb<0 else 0
            jj += (z_2)*pow(xa,k)*pow(ya,l)*dR*Rdtheta
    return jj

def VortexArea(z,zb,xy_map,dR_ls,Rdtheta_ls):
    #Z: 2D numpy array (lat,lon)
    ss00 = Vortex_Geo(z,zb,xy_map,dR_ls,Rdtheta_ls, 0, 0)
    area = ss00/zb
    return area

def Centroid(z,zb,xy_map,dR_ls,Rdtheta_ls):
    #Z: 2D numpy array (lat,lon)
    ss00 = Vortex_Geo(z,zb,xy_map,dR_ls,Rdtheta_ls, 0, 0)
    ss10 = Vortex_Geo(z,zb,xy_map,dR_ls,Rdtheta_ls, 1, 0)
    ss01 = Vortex_Geo(z,zb,xy_map,dR_ls,Rdtheta_ls, 0, 1)
    
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
    R2 = pow(x_,2) + pow(y_,2)
    lat_ = np.rad2deg(np.arcsin(1-R2/2))
    #method3
    #R2 = pow(x_,2) + pow(y_,2)
    #sinlat_ = (1-R2)/(1+R2)
    #lat_ = np.rad2deg(np.arcsin(sinlat_))
    
    #return:
    #The centroid of Votrex at the given moment
    #(lon_,lat_,) = in polar coord; (x_,y_) = in Cartsian coord
    return (lon_,lat_), (x_,y_) 

def Angel_AspectRatio(z,zb,xy_map,dR_ls,Rdtheta_ls, x_, y_):
    #Z: 2D numpy array (lat,lon)
    jj11 = Vortex_Geo_centralized(z,zb,xy_map,dR_ls,Rdtheta_ls, 1, 1, x_, y_)
    jj20 = Vortex_Geo_centralized(z,zb,xy_map,dR_ls,Rdtheta_ls, 2, 0, x_, y_)
    jj02 = Vortex_Geo_centralized(z,zb,xy_map,dR_ls,Rdtheta_ls, 0, 2, x_, y_)
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
    
    #return
    #psi = the angle between x-axis and major axis of ellipse (in {0,180))
    #r = aspect ratio, the ratio of major axis to minor axis (>1)
    return psi,r

def Major_Minor_axis(area,r):
    minor = np.sqrt(area/(np.pi*r))*2
    major = minor*r
    return major, minor
    #major = 2a, minor = 2b
    
    
    
    