#import modules

import fenics as fex
import matplotlib.pyplot as plt
import numpy as np


#initializing variables
#declaring  parameters
alpha = 3
beta = 1.2
u_D = fex.Expression('1 + x[0] * x[0] + alpha * x[1] * x[1] + beta * t', degree = 2, alpha = alpha, beta = beta, t =0)

u_D.t = t

def boundary(x, on_boundary):
	return on_boundary

bc = fex.DirichletBC(V, u_D, on_boundary)

u_n = project(u_D, V)
u_n = interpolate(u_D< V)

u = fex.TrialFunction(V)
v = fex.TestFunctin(V)
f = fex.Constant(beta - 2 - 2 * alpha)

F = u * v * dx + dt * fex.dot(fex.grad(u) ,  fex.grad(v)) * dx - (u_n + dt * f) * dx
a, L = fex.lhs(F), fex.rhs(F)

u = fex.Function(V)
t = 0
for n in range(num_steps):
	
	#Update current time
	t += dt
	u_D.t = t
	
	#Solve variational problem
	fex.solve(a == L, u, bc)
	
	#Update previous solution
	u_n.assign(u)

u_e = interpolate(u_D, V)
error = np.abs(u_e.vector().array() - u.vector().array())
print("t = %.2f: error = %.3g" %(t, erro))












