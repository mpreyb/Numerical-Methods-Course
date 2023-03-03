# -*- coding: utf-8 -*-
"""
    File name: punto3_tarea2.py
    Author: Felipe Otalvaro Agudelo
    Date created: 25/09/2021
    Date last modified: 10/11/2021
    Python Version: 3.8
"""
#se llaman los paquetes
import meshio 
import numpy as np
import math as mt
from scipy.special import jv


#se lee el archivo obtenido en gmsh
membrana_malla = meshio.read("punto3_tarea2.msh")
pts =membrana_malla.points
tets = membrana_malla.cells[0].data

#Matriz nxm donde n es la raíz y m es la función de Bessel
J=np.array([[2.04048,3.8317,5.1356], \
            [5.5201,7.0156,8.4172], \
            [8.6537,10.1735,11.6198]]) 

#Función solución
n=2 
m=2
a=0.1

#se cambia de coordernadas cartesianas a cilindrícas
def coordenadas_car_cilin(nodo):
    x = nodo[0]
    y = nodo[1]
    rho = mt.sqrt((x**2 + y**2))
    phi = mt.atan(y/x)
    return rho, phi


#se define la función de Bessel
def bessel (v,x):
    v=jv(v,x)
    #print(v)
    return v


#se define la función de los modos de vibración
def modos_vibracion(nodo,n,m,a,t):
    rho, phi = coordenadas_car_cilin(nodo)
    lbd=J[n,m]/a
    u=(mt.cos(m*phi))*bessel(m,(lbd*rho))*mt.sin(lbd*t)
    return u

T=np.linspace(0.1,1,500) #función para generar 500 números entre 0.1 y 1
p=1
for t in T:
    i=0
    for pt in pts:
        z=modos_vibracion(pt,n,m,a,t)*0.1
        pts[i,2]=z
        i=i+1
    malla = meshio.Mesh(
    pts,[("triangle",tets)],
    point_data = {"zp":pts[:,2]}
    # cell_data = {} 
    )
    filename='membrana_t'+'_modo_'+str(n)+'x'+str(m)+' '+str(p)+'.vtk'
    meshio.write(filename, malla,"vtk") #se crean los archivos VTK para visualizar en Paraview
    p=int(p)+1


