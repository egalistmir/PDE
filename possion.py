"""CS tutorial demo program: Poisson equation with Dirichlet conditions.
Test problem is chosen to give an exact solution at all nodes of the mesh.
  -Laplace(u) = f    in the unit square
            u = u_D  on the boundary
  u_D = 1 + x^2 + 2y^2
    f = -6
"""

from __future__ import print_function
import fenics as fex
import matplotlib.pyplot as plt

# Create mesh and define function space
mesh = fex.UnitSquareMesh(8, 8)
V = fex.FunctionSpace(mesh, 'P', 1)

# Define boundary condition
u_D = fex.Expression('1 + x[0]*x[0] + 2*x[1]*x[1]', degree=2)

def boundary(x, on_boundary):
    return on_boundary

bc = fex.DirichletBC(V, u_D, boundary)

# Define variational problem
u = fex.TrialFunction(V)
v = fex.TestFunction(V)
f = fex.Constant(-6.0)
a = fex.dot(fex.grad(u), fex.grad(v))*fex.dx
L = f*v*fex.dx

# Compute solution
u = fex.Function(V)
fex.solve(a == L, u, bc)
print(u)
# Plot solution and mesh
plt.subplot(2, 1, 1)
fex.plot(u)
#define the boundary conditions
plt.subplot(2, 1, 2)
fex.plot(mesh)
plt.savefig('possion_fex2.png')
#define the boundary conditions
# Save solution to file in VTK format
vtkfile = fex.File('poisson/solution.pvd')
vtkfile << u

# Compute error in L2 norm
error_L2 = fex.errornorm(u_D, u, 'L2')

# Compute maximum error at vertices
vertex_values_u_D = u_D.compute_vertex_values(mesh)
vertex_values_u = u.compute_vertex_values(mesh)
import numpy as np
error_max = np.max(np.abs(vertex_values_u_D - vertex_values_u))

# Print errors
print('error_L2  =', error_L2)
print('error_max =', error_max)

# Hold plot
#plt.show()
