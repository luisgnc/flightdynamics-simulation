import sympy as sp
import math

#Flags
i = 0


#Simbolos
u, alpha, theta, q, delta_t, delta_e = sp.symbols('u, alpha, theta, q, delta_t, delta_e')
du, dalpha, dtheta, dq, d_delta_t,d_delta_e = sp.symbols('du, dalpha, dtheta, dq, d_delta_t, d_delta_e')
rho, S, L, C_L, C_D, T_max, epsilon_t, m = sp.symbols('rho, S, L, C_L, C_D, T_max, epsilon_t, m')
c_media, C_M, I_yy = sp.symbols('c_media, C_M, I_yy')


#cte SI UNITS
u_0 = 0
alpha_0 = 0
theta_0 = 0
q_0 = 0
x_0 = sp.Matrix([u_0, alpha_0, theta_0, q_0])
u_0 = sp.Matrix([0, 0])
u_1 = sp.Matrix([0.5, 0])
u_2 = sp.Matrix([0.5, 0.1])
u_3 = sp.Matrix([0.5, 0])
u_4 = sp.Matrix([0.25, 0])


g = 9.80665
rho_ref = 1.225
h = 400

#ISA MODEL
if h < 11000:
    rho_h = rho_ref * (1 - 2.2558 * h * 10 ** -5) ** 4.256060537
else:
    rho_h = 0.29707 * rho_ref * math.e **((1.576939464 * 10 ** -4) * (11000 - h))


#Matrizes Constantes

A = sp.Matrix([[0.017, 19.03, -32.193, 0],
              [-0.0005, -1.56, 0.0056, 0.612],
              [0, 0, 0, 1],
              [0.0002, -3.80, -0.0034, -2.51]])

B = sp.Matrix([[0.0673, 0.0218],
              [0, -0.0012],
              [0, 0],
              [0, -0.0374]])

vp_A = sp.Matrix(A).eigenvals()

# Lista do tipo [Real1, Imaginario1][Real2, Imaginario 2]
vp_A_lista = list()

# Dicionario tipo Tpi: Value, w_ni: Value, zeta_i: Value
controlval = dict()

for chave, valor in vp_A.items():
    new = list()
    new.append(sp.re(chave))
    new.append(sp.im(chave))
    vp_A_lista.append(new)

    if i % 2 == 0:
        w = (sp.re(chave)** 2 + sp.im(chave)** 2) ** 0.5
        zeta = -sp.re(chave)/w
        tp = 2*math.pi/(w*math.sqrt(1 - zeta**2))

        controlval[f'Tp_{i}'] = tp
        controlval[f'w_n{i}'] = w
        controlval[f'zeta_{i}'] = zeta

    i += 1

print(controlval)



'''
def longitudinal_model(x, u):
    du = (1/m)*(0.5*rho*S * (u**2/(sp.cos(alpha))**2) * (C_L*sp.sin(alpha) - C_D*sp.cos(alpha)) +
                delta_t*T_max*sp.cos(epsilon_t)) - g*sp.sin(theta) - q*u*sp.tan(alpha)

    dalpha = q + g*sp.cos(alpha)*sp.cos(theta-alpha)/u - sp.cos(alpha)/(m*u)*(L -
            delta_t*T_max*sp.cos(epsilon_t)*sp.sin(epsilon_t)*sp.sin(alpha))

    dtheta = q

    dq = rho*u**2*S*c_media*C_M/(2*I_yy*(sp.cos(alpha))**2)
    

x_state = sp.Matrix([u, alpha, theta, q])
y_state = sp.Matrix([u, theta, -alpha, q])
u_state = sp.Matrix([delta_t, delta_e])
'''
