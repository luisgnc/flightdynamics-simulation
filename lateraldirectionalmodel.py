import sympy as sp

#Matrizes Constantes

A = sp.Matrix([[-0.19, 0.14, 0.065, -0.16],
              [0, 0, 1, 0.14],
              [-4.58, 0, -3.5, 0.87],
              [2.041, 0, -0.3, -0.72]])

B = sp.Matrix([[0, 0.004],
              [0, 0],
              [0.021, 0.0058],
              [0.0010, -0.134]])

vp_A = sp.Matrix(A).eigenvals()

vp_A_mod = dict()

# Há valores demasiados baixos(O(-64)) onde... Im(H) ~= 0, logo nega-se esses valores
for chave, valor in vp_A.items():
    if abs(sp.im(chave)) < 1e-10:
        chave = sp.re(chave)
        vp_A_mod[chave] = valor
    else:
        vp_A_mod[chave] = valor

print(vp_A_mod)
