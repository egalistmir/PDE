"""CS tutorial demo program: Heat equation with Dirichlet conditions.
Test problem is chosen to give an exact solution at all nodes of the mesh.
  u'= Laplace(u) + f  in the unit square
  u = u_D             on the boundary
  u = u_0             at t = 0
  u = 1 + x^2 + alpha*y^2 + beta*t
  f = beta - 2 - 2*alpha
"""
from fenics import *
import numpy as np
import matplotlib.pyplot as plt
#import time

T = 2.0            # final time
num_steps = 50     # number of time steps
dt = T / num_steps # time step size
alpha = 3          # parameter alpha
#beta = 1.2         # parameter beta

# Create mesh and define function space
nx = ny = 30
mesh = RectangleMesh(Point(-2, -2), Point(2, 2), nx, ny)
V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition
def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, Constant(0), boundary)

# Define initial value
u_0 = Expression('exp(-a*pow(x[0], 2) - a * pow(x[1], 2))',
                 degree=2, a = 5)

u_n = interpolate(u_0, V)
#u_n = project(u_D, V)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(0)

F = u*v*dx + dt*dot(grad(u), grad(v))*dx - (u_n + dt*f)*v*dx
a, L = lhs(F), rhs(F)

# Time-stepping
u = Function(V)
t = 0
for n in range(num_steps):

    # Update current time
    t += dt
    #u_D.t = t

    # Compute solution
    solve(a == L, u, bc)

    # Plot solution
    plot(u)
    plt.savefig('ft03_heat_diffusion%.d.png' %n)

    # Compute error at vertices
    #u_e = interpolate(u_D, V)
    #error = np.abs(u_e.vector().array() - u.vector().array()).max()
    #print('t = %.2f: error = %.3g' % (t, error))

    # Update previous solution
    u_n.assign(u)


