# -*- coding: utf-8 -*-
"""
    File name: caballo_tarea1.py
    Author: Maria Paula Rey
    Date created: 12/08/2021
    Date last modified: 19/08/2021
    Python Version: 3.8
"""

# Importing necessary libraries
import numpy as np 
import meshio

mesh = meshio.read("knight.msh")                 # Loading the mesh
pts = mesh.points                                # Storing points
tets = mesh.cells[0].data                        # Storing elements (tetrahedral)

# ---------------------------------------Volume calculation------------------------------------------------
V_total = 0.0                                    # Setting the volume to 0.
tets_rows = tets.shape[0]                        # Number of rows
det_tets = np.ones((4,4), dtype = "float64")     # 4X4 array filled with ones (Double precision floats). para sacar el determinante de tets


# the matrix is filled with the information of each tetrahedron and its respective nodes, thus being able to find the volume of the n-th tetrahedron
for i in range (0,tets_rows):                    # Cycle through the number of elements
    for r in range (0,4):                        # Cycle through the number of rows (0 to 3)
        for c in range (0,4):                    # Cycle through the number of columns (0 to 3)
        
            if c != 3:                           # pts matrix has only 3 columns. (if column /= 3 (last column))                    
                det_tets[r,c] = pts[tets[i,r],c]
    
    det = np.linalg.det(det_tets)                # Compute determinant
    
    # V=volume of a square prism. V/6 is the volume of a tetrahedron
    V = (1/6)*abs(det)                           # Absolute value of determinant divided by 6.
    V_total = V_total + V                        # Sum of all volumes (of each tetrahedron)

# We print the total volume of the knight
print('The total volume is',V_total,'mm^3')