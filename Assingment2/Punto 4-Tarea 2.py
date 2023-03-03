# -*- coding: utf-8 -*-
"""
Punto 4-Tarea Métodos Numéricos
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from copy import deepcopy
import math

fig, ax = plt.subplots()


p = 2   #Densidad volumétrica de la masa
E = 3   #Módulo de Young 
v = math.sqrt(E/p)
x_left = 0.
x_right = 1.
t = 0.
#Se van a discretizar 101 puntos
num_points = 101
dx = (x_right-x_left)/(num_points-1)
#Reserve one extra for Neumann boundary
x = np.linspace(x_left, x_right+dx, num_points+1)
# u = x/Initialize at 1000
u = 1e-3*x
#Condiciones de frontera de Neumann
u[-1] = u[-2]

u_before = deepcopy(u)
u_next = deepcopy(u)

#Acumula el valor de u a intervalos regulares
u_at_t = [deepcopy(u)]


dt = 0.1 * dx/c

count = 0


while t < 0.5:
    print(np.diff(u,2))
    u_next[1:num_points] = 2.*u[1:num_points] - u_before[1:num_points] + c*c*dt*dt/(dx*dx) * np.diff(u, 2)

    #Condiciones de frontera de Neumann
    u_next[-1] = u_next[-2]
    u_next[0] = 0.
    u_before = deepcopy(u)
    u = deepcopy(u_next)
    u_next = np.zeros(u.shape)
    t += dt
    count += 1
    #Once every 16 steps
    if count == 16:
        u_at_t.append(deepcopy(u))
        count = 0


def update(u, x):
    plt.clf()
    plt.xlim(x_left, x_right)
    plt.ylim(-1e-3, 1e-3)
    plt.plot(x, u)


#Animación
anim = FuncAnimation(fig, update, fargs=(
    x,), frames=u_at_t, interval=10)

plt.show()

