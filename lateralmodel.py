import sympy as sp
import math

#Flags
i = 0
#Simbolos
u, alpha, theta, q, delta_t, delta_e = sp.symbols('u, alpha, theta, q, delta_t, delta_e')
x1_dot, x2_dot, x3_dot, x4_dot = sp.symbols('x1_dot, x2_dot, x3_dot, x4_dot')

x_symbol = sp.Matrix([u, alpha, theta, q])
u_symbol = sp.Matrix([delta_t, delta_e])
y_symbol = sp.Matrix([u, theta, -alpha, q])

#Matrizes Constantes
A = sp.Matrix([[0.017, 19.03, -32.193, 0],
              [-0.0005, -1.56, 0.0056, 0.612],
              [0, 0, 0, 1],
              [0.0002, -3.80, -0.0034, -2.51]])
B = sp.Matrix([[0.0673, 0.0218],
              [0, -0.0012],
              [0, 0],
              [0, -0.0374]])
C = sp.Matrix([[1, 0, 0, 0],
              [0, 1, 0, 0]])
vp_A = sp.Matrix(A).eigenvals()

# Lista do tipo [Real1, Imaginario1][Real2, Imaginario 2]
vp_A_lista = list()
for _ in vp_A.keys():
    new = list()
    new.append(sp.re(_))
    new.append(sp.im(_))
    vp_A_lista.append(new)
# Dicionario tipo Tpi: Value, w_ni: Value, zeta_i: Value
controlval = dict()
for chave, valor in vp_A.items():
    if i % 2 == 0:
        w = (sp.re(chave)** 2 + sp.im(chave)** 2) ** 0.5
        zeta = abs(-sp.re(chave))/w
        tp = 2*math.pi/(w*math.sqrt(1 - zeta**2))

        controlval[f'Tp_{i}'] = tp
        controlval[f'w_n{i}'] = w
        controlval[f'zeta_{i}'] = zeta
    i += 1
for k, v in controlval.items():
    if 'Tp' in k:
        if v > 10:
            print('\nMovimento Fugóide')
        else:
            print('\nMovimento de Período Curto')
    print(f"{k} = {v}")

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
print("Matriz de Observabilidade")
sp.pretty_print(Gamma)

if A.shape[0] == Delta.rank():
    print(f"\nO sistema é controlável. n = {A.shape[0]} , , caract(Delta) = {Delta.rank()}")
else:
    print(f"\nO sistema não é controlável.")
