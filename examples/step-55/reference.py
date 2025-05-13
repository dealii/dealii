from sympy import *
from sympy.printing import print_ccode
from sympy.physics.vector import ReferenceFrame, gradient, divergence
from sympy.vector import CoordSysCartesian

R = ReferenceFrame("R")
x = R[0]
y = R[1]

a = -0.5
b = 1.5
visc = 1e-1
lambda_ = 1 / (2 * visc) - sqrt(1 / (4 * visc ** 2) + 4 * pi ** 2)
print(" visc=%f" % visc)

u = [0, 0]
u[0] = 1 - exp(lambda_ * x) * cos(2 * pi * y)
u[1] = lambda_ / (2 * pi) * exp(lambda_ * x) * sin(2 * pi * y)
p = (exp(3 * lambda_) - exp(-lambda_)) / (8 * lambda_) - exp(2 * lambda_ * x) / 2
p = p - integrate(p, (x, a, b))

grad_p = gradient(p, R).to_matrix(R)
f0 = -divergence(visc * gradient(u[0], R), R) + grad_p[0]
f1 = -divergence(visc * gradient(u[1], R), R) + grad_p[1]
f2 = divergence(u[0] * R.x + u[1] * R.y, R)

print("\n * RHS:")
print(ccode(f0, assign_to="values[0]"))
print(ccode(f1, assign_to="values[1]"))
print(ccode(f2, assign_to="values[2]"))


print("\n * ExactSolution:")
print(ccode(u[0], assign_to="values[0]"))
print(ccode(u[1], assign_to="values[1]"))
print(ccode(p, assign_to="values[2]"))

print("")
print("pressure mean:", N(integrate(p, (x, a, b))))
