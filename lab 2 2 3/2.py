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


k = np.array([0.030779491, 0.031353089, 0.032977581, 0.034662892, 0.036246012]) * 1000
T = np.array([24, 29.9, 45.1, 65.1, 80.1])

plt.plot(T, k, ".", color="black", markersize="8", label=r"$\kappa(t)$")
plt.xlabel(r"$t, \ ^{\circ}C$")
plt.ylabel(r"$\kappa, мВт/(Ом \cdot К)$")
plt.legend()



plt.show()