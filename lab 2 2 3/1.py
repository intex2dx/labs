import numpy as np
import matplotlib.pyplot as plt

def xi2(x,y, y_err):
    #y = a + bx
    N = len(x)
    avx = sum(x/(y_err**2))/sum(1/(y_err**2))
    avy = sum(y/(y_err**2))/sum(1/(y_err**2))
    av2x = sum(x*x/(y_err**2))/sum(1/(y_err**2))
    av2y = sum(y*y/(y_err**2))/sum(1/(y_err**2))
    R_xy = (sum(x*y/(y_err**2))/sum(1/(y_err**2))-avx*avy)
    S_x = av2x - avx**2
    S_y = av2y - avy**2
    b = R_xy/S_x
    a = avy - b * avx
    r = R_xy/(np.sqrt(S_x*S_y))
    delta_b = np.sqrt(1/(N-2)) * np.sqrt(S_y/S_x - b**2)
    delta_a = delta_b * np.sqrt(av2x)
    f = a + b * x
    xi = sum(((y-f)/y_err)**2)
    return a, b, delta_a, delta_b, xi,  r

r0 = np.array([0.02012, 0.02048, 0.02161, 0.02302, 0.02409])
r0 *= 1000
T = np.array([24, 29.9, 45.1, 65.1, 80.1])
dr0 = np.array([4.46, 4.72, 5.31, 3.21, 2.80])/1000
dT = np.array([0.1] * 5)
plt.errorbar(T, r0, dr0, dT, ls="none", marker="o", markersize=5, color="black")
a, b, delta_a, delta_b, xi,  r = xi2(T, r0, dr0)
xdof = xi/3
print(a, b, delta_a, delta_b, xi,  r)
plt.plot([23.2, 80.8], [a + b*23.2, a + b*80.8], color="red", label=r'''$R_0(t) = (18.39 \pm 0.3) + (0.0712 \pm 0.0004)t$
$r \approx 0.99995$
''')
plt.legend()
plt.grid(True)
plt.xlabel(r"$t, \ ^{\circ}C$")
plt.ylabel(r"$R_0, Ом$")
plt.show()
