from math import pi
import sys
from types import NoneType
from typing import Callable
from matplotlib import markers, pyplot as plt
from numpy import around
from functions import *

# Définition de l'équation et des conditions initiales

def f(x, y, y1): # Équation différentielle, de la forme y''=f(x,y,y')
    return y+1

xa=1 # Position de la première condition initiale xa
xb=2 # Position de la seconde condition initiale xb

ya=5 # Valeur de y(xa)
yb=25 # Valeur de y(xb)

dy1a=abs(ya-yb)/5 # Pas de tir, si la cible est manqué ou est inatteignable 
#modifier le pas de tir peut permettre l'obtention de meilleurs résultats

n= 1000 # Permet de définir le pas utilisé dans l'algorithme d'Euler h = (xb-xa)/n

#####################################################################################


def shot(f: Callable, xa: float, xb: float, ya: float, yb: float, dy1a: int, n: int):  # dy1a = m1-m0 /  n =
    tirs = []
    y1a = (yb-ya)/(xb-xa)   # initialisation de m0
    y1 = y1a
    yb_shot = Euler(f, xa, xb, ya, y1, n)["yb"]
    list_m = []
    list_yb = []

    preshot1 = Euler(f, xa, xb, ya, y1, n)["yb"]
    preshot2 = Euler(f, xa, xb, ya, y1+1/2, n)["yb"]
    is_growing = preshot2-preshot1>=0


    if yb_shot < yb and is_growing:  # on vérifie si le tir initial est trop court ou trop long
        while yb_shot < yb:  # on effectue des tirs par incrément tant que le tir
            #   est trop court en s'assurant de bien dépassé la cible
            new_shot = Euler(f, xa, xb, ya, y1, n)
            yb_shot = new_shot["yb"]
            tirs.append(new_shot)
            list_yb.append(yb_shot)
            list_m.append(y1)
            y1 += dy1a  # On incrémente la valeur de mk
    elif yb_shot > yb or not is_growing:
        while yb_shot > yb:
            new_shot = Euler(f, xa, xb, ya, y1, n)
            yb_shot = new_shot["yb"]
            tirs.append(new_shot)
            list_yb.append(yb_shot)
            list_m.append(y1)
            y1 -= dy1a  # On incrémente la valeur de mk

    return ({"tirs": tirs, "list_m": list_m, 'list_yb': list_yb, 'y1a': y1a})


new_shot = shot(f, xa,xb, ya, yb, dy1a, n)

tirs, list_m, list_yb, y1a = new_shot.values()


# Création avec matplotlib.pyplot de deux graphiques sur une même fenêtre


figure, (graph1,graph2) = plt.subplots(2)


# Trace la courbe de chaque tir

for i, tir in enumerate(tirs):
    graph1.plot(tir["x_array"],tir["y_array"], label = 'Tir '+str(i+1))



# Nous traçons la droite passant par le point de départ du tir et la cible 

x = np.linspace(list_m[0], list_m[-1], n*5, endpoint=True) # Intervalle de définition [m0,...,mk] du polynôme interpolateur de Lagrange
x_droite = np.linspace(xa, xb, n*5, endpoint=True) # on définit l'intervalle de définition [xa,xb] de la droite reliant (xa,ya) et (xb,yb)


graph1.plot(x_droite, y1a*(x_droite-xa) + ya, label='Droite reliant les deux point',c="red") 
graph2.plot(x_droite, y1a*(x_droite-xa) + ya, label='Droite reliant les deux point',c="red")


interpolation_m = lagrange(list_m, list_yb) # Polynôme de Lagrange avec noeuds: liste des m et image: liste des yb

best_y1a=None
try:
    for k in x:
    
        if np.around(interpolation_m(k),2) == np.around(yb,2): # on cherche k tel que le polynôme est égal à yb au 100eme près
            best_y1a = k
    
    if type(best_y1a)==NoneType: raise Exception("Une erreur est survenue. Essayez d'ajuster vos paramètres.")

except Exception as e:
    print(e)
    sys.exit()
    
            
            

best_shot = Euler(f, xa, xb, ya, best_y1a, n)

arounded_y1a = around(best_y1a,3)



# Définition des titres des graphs 

graph1.set_title('Graph des tirs effectués')
graph2.set_title('Graph du meilleur tir pour m='+str(arounded_y1a))

graph2.plot(best_shot["x_array"],best_shot["y_array"],label="Meilleur tir pour m=" + str(arounded_y1a))



# On place les points de départ et d'arrivée des tirs 

graph1.scatter(xa,ya, c="red",label="Point de départ du tir en ("+str(np.around(xa,2))+","+str(np.around(ya,2))+")")
graph2.scatter(xa,ya, c="red",label="Point de départ du tir en ("+str(np.around(xa,2))+","+str(np.around(ya,2))+")")

graph1.scatter(xb,yb, c="red",marker="x",label="Cible en ("+str(np.around(xb,2))+","+str(np.around(yb,2))+")")
graph2.scatter(xb,yb, c="red",marker="x",label="Cible en ("+str(np.around(xb,2))+","+str(np.around(yb,2))+")")

# Affichage des légendes et des graphiques

graph2.legend()
graph1.legend()


print('θ(m_k):' + str(np.around(list_yb, 2)))
print("La meilleure approximation de y'(a) est " + str(best_y1a))

plt.show()



