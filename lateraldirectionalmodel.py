import sympy as sp
import math

#Flags
i = 0
#Simbolos
beta, phi, p, r, delta_a, delta_r = sp.symbols('beta, phi, p, r, delta_a, delta_r')
x1_dot, x2_dot, x3_dot, x4_dot = sp.symbols('x1_dot, x2_dot, x3_dot, x4_dot')

x_symbol = sp.Matrix([beta, phi, p, r])
u_symbol = sp.Matrix([delta_a, delta_r])
y_symbol = sp.Matrix([beta, phi])

#Matrizes Constantes
A = sp.Matrix([[-0.19, 0.14, 0.065, -0.16],
              [0, 0, 1, 0.014],
              [-4.58, 0, -3.5, 0.87],
              [2.041, 0, -0.3, -0.72]])

B = sp.Matrix([[0.0673, 0.004],
              [0, 0],
              [0.21, 0.0058],
              [0.001, -0.134]])

C = sp.Matrix([[1, 0, 0, 0],
              [0, 1, 0, 0]])

vp_A = sp.Matrix(A).eigenvals()
# Lista do tipo [Real1, Imaginario1][Real2, Imaginario2]
vp_A_lista = [[-3.39334582954197, 0],[-0.465628973363563, -0.678870059572729], [-0.0853962237309086,0],
        [-0.465628973363563, 0.678870059572729]]
print(vp_A_lista)

# Dicionario tipo Tpi: Value, w_ni: Value, zeta_i: Value
controlval = dict()
for ind in vp_A_lista:
    w = (ind[0] ** 2 + ind[1]** 2) ** 0.5
    zeta = abs(ind[0])/w
    if zeta < 1:
        tp = 2*math.pi/(w*math.sqrt(1 - zeta**2))
    else:
        tp = 'undefined'
    print(f"{tp} , {w} , {zeta}")
    controlval[f'Tp_{i}'] = tp
    controlval[f'w_n{i}'] = w
    controlval[f'zeta_{i}'] = zeta
    i += 1
x_dot = sp.Matrix([x1_dot, x2_dot, x3_dot, x4_dot])
x_dot = A*x_symbol + B*u_symbol
print(controlval)

# Matriz de Contrabilidade
controllability_blocks = [B]
for i in range(1, A.shape[0]):
    controllability_blocks.append(A @ controllability_blocks[-1])
Delta = sp.Matrix.hstack(*controllability_blocks)
# Matriz de Observabilidade
observability_blocks = [C]
for i in range(1, A.shape[0]):
    observability_blocks.append(observability_blocks[-1] @ A)
Gamma = sp.Matrix.vstack(*observability_blocks)

print("Matriz de Controlabilidade")
sp.pretty_print(Delta)
print("\n\n")
print("Matriz de Observabilidade")
sp.pretty_print(Gamma)

if A.shape[0] == Delta.rank():
    print(f"\nO sistema e controlavel. n = {A.shape[0]} , , caract(Delta) = {Delta.rank()}")
else:
    print(f"\nO sistema nao e controlavel.")

if A.shape[0] == Gamma.rank():
    print(f"\nO sistema e observavel. n = {A.shape[0]} , , caract(Gamma) = {Gamma.rank()}")
else:
    print(f"\nO sistema nao e observavel.")
