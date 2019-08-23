import fenics as fex
import numpy as np
import matplotlib.pyplot as plt
import mshr

mesh = fex.UnitSquareMesh()
V = fex.

u 
u_D = fex.Expression()

bc = fex.Dirichlect()

def boundary(x, on_boundary)
	return on_boundary

fex.solve(u, x, bc)
fex.plot(u)
fex.plot(mesh)

