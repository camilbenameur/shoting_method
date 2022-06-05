from typing import Callable
import numpy as np
from scipy.interpolate import lagrange


# Permet d'intégrer numériquement l'équation différentielle à l'aide de la méthode d'Euler
def Euler(f: Callable, xa: float, xb: float, ya: float, y1a: float, n: int):
    y_list = []
    h = (xb - xa) / float(n) # V
    
    x_list = np.linspace(xa,xb,n,endpoint=True)  

    x = x_list[0]
    y = ya
    zn = y1a
    z_temp = zn


    for i in x_list:
        z_temp = zn
        zn += h * f(x, y, zn)
        y += h * z_temp
        x=i
        y_list.append(y)

       
    return ({"yb": y, "x_array":  x_list, "y_array": np.array(y_list)})


def Lagrange(x_array, y_array):
    def poly(x):
        # Effectue une interpolation de Lagrange
        return lagrange(x_array, y_array)(x)
    return poly
