"""
    File name: punto2_tarea2_2.py
    Author: Maria Paula Rey
    Date created: 12/09/2021
    Date last modified: 28/09/2021
    Python Version: 3.8
"""

import meshio as mes
import numpy as np
import sympy as sym
from sympy import Matrix
from sympy import *
from sympy import Array
from sympy import diff, exp
    
mesh=mes.read('concentric_cylinders.vtk')          # Loading the mesh
points=mesh.points                                 # Storing points
cell=mesh.cells                                    # Storing elements (squares)
cell2=cell[0]
elements=cell2[1]                                  # Contains the number of the elements and their connections
Potential=mesh.point_data
Potential=Potential['electrostatic_potential']
 

# Constants
epsilon=8.8541878128E-12     # Vacuum permittivity
deltaV=1                     # Potential difference
Utot=0                       # Energy

# ---------------------------------------Base functions------------------------------------------------
# symbolic variables are defined
r = sym.symbols("r")
s = sym.symbols("s")
N0 = (1/4)*(1-r)*(1-s)  
N1 = (1/4)*(1+r)*(1-s)
N2 = (1/4)*(1+r)*(1+s)
N3 = (1/4)*(1-r)*(1+s)
Nt = Array([N0, N1, N2, N3]) 
print(Nt)

dNR_0 = N0.diff(r);
dNR_1 = N1.diff(r);
dNR_2 = N2.diff(r);
dNR_3 = N3.diff(r);

dNS_0 = N0.diff(s);
dNS_1 = N1.diff(s);
dNS_2 = N2.diff(s);
dNS_3 = N3.diff(s);

dNR = Array([dNR_0, dNR_1, dNR_2, dNR_3])
#print(dNR)
dNS = Array([dNS_0, dNS_1, dNS_2, dNS_3])

# Mesh definition (with respect to points and nodes)
def mapping(l):
    nodes=elements[l]
    pts0=points[nodes[0]]
    pts1=points[nodes[1]]
    pts2=points[nodes[2]]
    pts3=points[nodes[3]]
    ptsx=Array([pts0[0], pts1[0], pts2[0], pts3[0]])
    ptsy=Array([pts0[1], pts1[1], pts2[1], pts3[1]])
    X=np.dot(Nt,ptsx)
    Y=np.dot(Nt,ptsy)
    return ptsx, ptsy, X, Y, nodes

# Jacobian Calculation 
def Jacobian(ptsx,ptsy):
    J0=np.dot(dNR,ptsx)
    J1=np.dot(dNR,ptsy)
    J2=np.dot(dNS,ptsx)
    J3=np.dot(dNS,ptsy)
    Jp=J=Array([J0, J1, J2, J3])
    J=Array([[J0, J1],[J2, J3]])
    detJ=J0*J3-J1*J2 ## Determinante del jacobiano 
    return Jp, J, detJ

# Electric potential
def potential_energy(nodes):
    P0=Potential[nodes[0]]
    P1=Potential[nodes[1]]
    P2=Potential[nodes[2]]
    P3=Potential[nodes[3]]
    P = Array([P0, P1, P2, P3])
    return P

# Gradient calculation
def gradient (Jp, dNR, dNS, pot, detJ):
    Grax=pot[0]*(dNR[0]*Jp[3] - dNS[0]*Jp[1]) + pot[1]*(dNR[1]*Jp[3] - dNS[1]*Jp[1]) + pot[2]*(dNR[2]*Jp[3] - dNS[2]*Jp[1]) + pot[3]*(dNR[3]*Jp[3] - dNS[3]*Jp[1])
    Gray= - pot[0]*(dNR[0]*Jp[2] - dNS[0]*Jp[0]) -pot[1]*(dNR[1]*Jp[2] - dNS[1]*Jp[0]) - pot[2]*(dNR[2]*Jp[2] - dNS[2]*Jp[0]) - pot[3]*(dNR[3]*Jp[2] - dNS[3]*Jp[0])
    Grax=(1/detJ)*Grax
    Gray=(1/detJ)*Gray
    Gradiente=Array([Grax,Gray])
    return Gradiente

# Electric field, using the gradient  
def Electric_Field (Gradiente):
    E=-Matrix(Gradiente)
    magnitudE = (E[0]**2 + E[1]**2)**(0.5)
    magnitudE2=magnitudE**2
    return magnitudE2

# Sum of the Gaussian square, to find the electric potential energy
def Gaussian_S(magnitudE2, detJ):
    R = 0.5773502691896258
    S=R
    gs=magnitudE2 * detJ
    G0=gs.subs([(r,-R),(s,-S)])
    G1=gs.subs([(r,R),(s,-S)])
    G2=gs.subs([(r,R),(s,S)])
    G3=gs.subs([(r,-R),(s,S)])
    G=G0+G1+G2+G3
    U=0.5*epsilon*G
    return U

# Total energy calculation
for elementos in range(0, len(elements)):
    mapp= mapping(elementos)
    Jac= Jacobian(mapp[0], mapp[1])
    POT = potential_energy(mapp[4])
    Grad = gradient(Jac[0], dNR, dNS, POT, Jac[2])
    magnitudE= Electric_Field(Grad)
    U= Gaussian_S( magnitudE, Jac[2])
    Utot= U +Utot
    
# ---------------------------------------Results------------------------------------------------
capacitance=2*Utot/deltaV                              # Capacitance calculation 
Cap_T = 2.266180070913597 * epsilon                    # Theoretical calculation
error = abs((Cap_T - capacitance) / Cap_T )* 100       # Error Calculation 
print('La capacitancia obtenida es: ', capacitance,'F.', 'La capacitancia teórica es: ', Cap_T, 'F.',' Porcentaje de error: ', error )

