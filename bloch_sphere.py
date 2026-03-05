# -*- coding: utf-8 -*-
"""
Created on Thu Mar  5 13:03:28 2026

@author: swest
"""

from qutip import *
from scipy import *
from numpy import *

def qubit_integrate(w, theta, gamma1, gamma2, psi0, tlist):
    # operators and the hamiltonian
    sx = sigmax(); sy = sigmay(); sz = sigmaz(); sm = sigmam()
    H = w * (cos(theta) * sz + sin(theta) * sx)
    # collapse operators
    c_op_list = []
    n_th = 0.5 # temperature
    rate = gamma1 * (n_th + 1)
    if rate > 0.0: c_op_list.append(sqrt(rate) * sm)
    rate = gamma1 * n_th
    if rate > 0.0: c_op_list.append(sqrt(rate) * sm.dag())
    rate = gamma2
    if rate > 0.0: c_op_list.append(sqrt(rate) * sz)


    # evolve and calculate expectation values
    output = mesolve(H, psi0, tlist, c_op_list, [sx, sy, sz])  
    return output.expect[0], output.expect[1], output.expect[2]


def rabi_osc(gamma, omega, omega_21, t):
    rabi_frequency = sqrt(gamma ** 2 + (omega - omega_21) ** 2 / 4)
    f = 1j * (omega - omega_21) / 2
    c1 = (exp(f * t) * cos(rabi_frequency * t) - (1j * (omega-omega_21) * t 
         / (2*rabi_frequency)) * exp(f * t) * sin(rabi_frequency * t))
    c2 = (-1j * gamma / rabi_frequency) * exp(-f * t) * sin(rabi_frequency * t)
    psi = c1 * basis(2, 0) + c2 * basis(2, 1)
    return psi


def simulate_rabi_bloch(gamma, omega_drive, omega21, psi0, tlist,
                        gamma1=0.0, gamma2=0.0):
    sx_op = sigmax()
    sy_op = sigmay()
    sz_op = sigmaz()

    # Detuning
    delta = omega_drive - omega21

    # Static Hamiltonian part (rotating-frame form)
    H0 = delta/2 * sz_op

    # Drive term
    H1 = gamma * sx_op
    H = [H0, [H1, lambda t, args: cos(omega_drive * t)]]

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
    
#gamma = 0.1
#omega = 10
#omega21 = 9

#tlist = linspace(0, 50, 400)

#sx = []
#sy = []
#sz = []

#for t in tlist:
 #   psi = rabi_osc(gamma, omega, omega21, t)

  #  sx.append(expect(sigmax(), psi))
   # sy.append(expect(sigmay(), psi))
    #sz.append(expect(sigmaz(), psi))

#sx = array(sx)
#sy = array(sy)
#sz = array(sz)

sx, sy, sz = simulate_rabi_bloch(
    gamma=10,
    omega_drive=2.5,
    omega21=2.5,
    psi0=basis(2,0),
    tlist=linspace(0,2,400)
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