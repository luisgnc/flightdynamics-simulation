import sympy as sp
from sympy import pretty_print


# Constantes
a, b, c, d, e = 5, 1, -5, 8, 7
x1, x2, x3, x4 = sp.symbols('x1 x2 x3 x4')
u1, u2 = sp.symbols('u1 u2')

#Funções Vel.
dx1 = x1**2 * sp.cos(x2**3) + x3 + a*x1*x4**2 + b*x2*u1
dx2 = x1 + x2**2 + x3 + c*x4
dx3 = d*x1*x2*x3 + u1 + u2
dx4 = x1**2 + x2*x3 + e*x4 + u1

x_state = sp.Matrix([x1, x2, x3, x4])
u_state = sp.Matrix([u1, u2])
dx_state = sp.Matrix([dx1, dx2, dx3, dx4])

#Vetor de Estado de Equilíbrio
eq_state_vector = {x1:1,
                   x2:0,
                   x3:-1,
                   x4: 0}

control_state_vector = {u1:-1,
                        u2:1}

x_jacobian = sp.Matrix([dx1, dx2, dx3, dx4]).jacobian(x_state)
u_jacobian = sp.Matrix([dx1, dx2, dx3, dx4]).jacobian(u_state)

#Avaliação das Expressões or matrix A, B
x_jacob_eval = x_jacobian.subs(eq_state_vector).subs(control_state_vector)
u_jacob_eval = u_jacobian.subs(eq_state_vector).subs(control_state_vector)

pretty_print(x_jacob_eval)
pretty_print(u_jacob_eval)
