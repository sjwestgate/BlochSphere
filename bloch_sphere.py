# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 13:03:28 2026

@author: swest
"""

from qutip import *
from scipy import *
from numpy import *

def simulate_rabi_bloch(gamma, omega_drive, omega21, psi0, tlist,
                        gamma1=0.0, gamma2=0.0):
    sx_op = sigmax()
    sy_op = sigmay()
    sz_op = sigmaz()

    # Detuning
    delta = omega_drive - omega21


    # Drive term
    delta = omega_drive - omega21

    H = (delta/2)*sz_op + (gamma/2)*sx_op


    # Collapse operators
    c_ops = []

    if gamma1 > 0:
        c_ops.append(sqrt(gamma1) * sigmam())

    if gamma2 > 0:
        c_ops.append(sqrt(gamma2) * sz_op)

    # Solve master equation
    result = mesolve(
        H,
        psi0,
        tlist,
        c_ops,
        [sx_op, sy_op, sz_op]
    )

    return result.expect
    
omega_d = 2.5
omega_res = 2.5
gam=1

sx, sy, sz = simulate_rabi_bloch(
    gamma=gam,
    omega_drive=omega_d,
    omega21=omega_res,
    psi0=basis(2,0),
    tlist=linspace(0,20,200),
    gamma1 = 0.2,
    gamma2 = 0.05
)

from pylab import *
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
sphere = Bloch(axes=ax)

def animate(i):
    ax.cla()
    sphere = Bloch(axes=ax)

    #sphere.add_vectors([np.sin(theta), 0, np.cos(theta)])
    sphere.add_points([sx[:i+1], sy[:i+1], sz[:i+1]])
    Omega_R = gam
    Delta = omega_d - omega_res

    axis = np.array([Omega_R, 0, Delta])
    axis = axis / np.linalg.norm(axis)

    sphere.add_vectors(axis)
    sphere.make_sphere()

    return []

ani = animation.FuncAnimation(
    fig,
    animate,
    frames=len(sx),
    interval=500,
    blit=False
)


ani.save("bloch_sphere.gif", writer="pillow", fps=20)

plt.close(fig)

tlist = linspace(0,20,200)
plt.plot(tlist,sz)
plt.xlabel('Time /s')
plt.ylabel('S_z')
plt.show()
