# Sorbonne Université LU3IN024 2021-2022
# TME 5 : Cryptographie à base de courbes elliptiques
#
# Etudiant.e 1 : Aime Cesaire Mugishawayo , 21340522
# Etudiant.e 2 : Alex Maloigne, 21105949
# Etudiante 3: Saarah Khodadin
import math
import time
from math import sqrt
import matplotlib.pyplot as plt
from random import randint


# Fonctions utiles

def exp(a, N, p):
    """Renvoie a**N % p par exponentiation rapide."""

    def binaire(N):
        L = list()
        while (N > 0):
            L.append(N % 2)
            N = N // 2
        L.reverse()
        return L

    res = 1
    for Ni in binaire(N):
        res = (res * res) % p
        if (Ni == 1):
            res = (res * a) % p
    return res


def factor(n):
    """ Return the list of couples (p, a_p) where p is a prime divisor of n and
    a_p is the p-adic valuation of n. """

    def factor_gen(n):
        j = 2
        while n > 1:
            for i in range(j, int(sqrt(n)) + 1):
                if n % i == 0:
                    n //= i
                    j = i
                    yield i
                    break
            else:
                if n > 1:
                    yield n
                    break

    factors_with_multiplicity = list(factor_gen(n))
    factors_set = set(factors_with_multiplicity)

    return [(p, factors_with_multiplicity.count(p)) for p in factors_set]


def inv_mod(x, p):
    """Renvoie l'inverse de x modulo p."""
    return exp(x, p - 2, p)


def racine_carree(a, p):
    """Renvoie une racine carrée de a mod p si p = 3 mod 4."""
    assert p % 4 == 3, "erreur: p != 3 mod 4"

    return exp(a, (p + 1) // 4, p)


""" OUR HELPEER METHODS"""


def delta(E):
    """Calculer le delta"""
    p, a, b = E

    disc = (4 * exp(a, 3, p) + 27 * exp(b, 2, p)) % p

    delt = (-16 * disc) % p

    if delt < 0:
        delt += p

    return delt


def mod(num, p):
    """Calculer la mod and accounts for negatives"""

    ret = num % p

    if ret < 0:
        ret += p

    return ret


def left(y, p):
    """ Helper function pour calculer le cote gauche de l'elliptique'"""
    return exp(y, 2, p)


def right(x, E):
    """ Helper function pour calculer le cote gauche de l'elliptique"""
    p, a, b = E
    return (exp(x, 3, p) + (a * x) + b) % p


def hasseBounds(p):
    """ Helper function pour calculer les bounds du theoreme de Hasse"""
    lower_bound = math.ceil(p + 1 - 2 * sqrt(p))
    upper_bound = math.floor(p + 1 + 2 * sqrt(p))

    return lower_bound, upper_bound


def est_inverse_additif(P1, P2, p):
    """ Helper function pour savoir si P2 est l'inverse additif de P1"""
    if P1 == () or P2 == ():
        return False

    return P1 == moins(P2, p)


def toBinaire(N):
    """ Helper function pour obtenir une representation binaire de N comme liste de bits."""
    L = list()
    while (N > 0):
        L.append(N & 1)
        N = N >> 1
    L.reverse()
    return L


# Fonctions demandées dans le TME


def est_elliptique(E):
    """
    Renvoie True si la courbe E est elliptique et False sinon.

    E : un triplet (p, a, b) représentant la courbe d'équation
    y^2 = x^3 + ax + b sur F_p, p > 3
    """

    """ Calcule delta et verifie qu'il est different de zero"""

    return delta(E) != 0


def point_sur_courbe(P, E):
    """Renvoie True si le point P appartient à la courbe E et False sinon."""
    if P == ():
        return True
    p, a, b = E
    x, y = P

    return left(y, p) == right(x, E)


def symbole_legendre(a, p):
    """Renvoie le symbole de Legendre de a mod p."""

    """
    Cette relation est guarantie par le critere de Euler;.
    Je l'ai de mon cours de Number Theory.
    
    PROOF:
    cas 1. p|a.  => p| (a ^ (p-1)/2 )
                donc, a ^ (p-1)/2 = 0 mod p = leg(a,p)
    caa 2. leg(a,p) = 1. donc, a = g ^ 2k mod p
                                a ^ (p-1)/2 = g^(2k(p-1))/2 mod p
                                            = g ^ (p-1)k mod p
                                            = 1 mod p, par le theorem de Fermat (ou par Z/pZ* est un groupe d' ordre p-1)
                                            = leg(a,p)
    cas 3. leg(a,p) = -1. donc a = g ^ (2k + 1) mod p
                            a ^ (p-1)/2 = (g^(p-1)k)  *  (g ^ (p-1)/2)
                                        = g ^ (p-1)/2 mod p
                            
                            Square both sides, donc a^(p-1) = 1 mod p
                            mais a ^ (p-1)/2 = g ^ (p-1)/2 mod p != 1 car ord(g) = p-1. (Z/pZ* est un groupe d' ordre p-1)
                            
                            donc a ^ (p-1)/2 = -1
                                            = leg(a,p)
    """

    leg = exp(a, (p - 1) // 2, p)

    return leg


def cardinal(E):
    """Renvoie le cardinal du groupe de points de la courbe E."""
    """
        La courbe est definie sur le corps Fp = Z/pZ. Donc on teste tout les p-1 elements de ce corps. 
        Si le cote droit de l'equation les deux racines carre, les deux solutions sont +- y.
        Si le cote droit egal a zero, alors y = 0
    """
    p, a, b = E
    card = 1  # l'element neutre

    for x in range(p):
        rhs = right(x, E)

        if symbole_legendre(rhs, p) == 1:
            card += 2
        elif rhs == 0:
            card += 1
    return card


def liste_points(E):
    """Renvoie la liste des points de la courbe elliptique E."""
    """
    On teste chaque x dans le corps sur lequel la courbe est definie.
    Si on trouve un x pour lequel le cote droit a les deux racine carres,
    on trouve les racines et les ajoutent a notre liste.
     Si le cote droit egal a zero, alors on ajoute (x,0)    
    """

    p, a, b = E

    assert p % 4 == 3, "erreur: p n'est pas congru à 3 mod 4."

    liste = [()]

    for x in range(p):
        rhs = right(x, E)
        if symbole_legendre(rhs, p) == 1:

            y = racine_carree(rhs, p)
            """ le fait que la racine carre de a mod p est (a ^ (p + 1)//4 ) mod p
            est utilisee dasns la fonction racine_carree"""

            liste.append((x, y))
            liste.append((x, -y))
        elif rhs == 0:
            liste.append((x, 0))

    return liste


# Verification de liste_points()
# Les tests dans test-4-liste-points.py font deja le travail, mais on peut aussi utiliser cette fonction
def verif_liste_points(E):
    """ Renvoie True si tous les points calculees par liste_points sont confirmes dans la liste"""
    liste = liste_points(E)

    for point in liste:
        if not point_sur_courbe(point, E):
            # print("[-]" + point + " n'est pas sur la courbe")
            return False

    return True


def cardinaux_courbes(p):
    """
    Renvoie la distribution des cardinaux des courbes elliptiques définies sur F_p.

    Renvoie un dictionnaire D où D[i] contient le nombre de courbes elliptiques
    de cardinal i sur F_p.
    """
    """
    On ajoute les cles possibles (= les cardinaux possibles) par le theoreme de Hasse.
    Pour chaque paire possible (a,b) sur Fp x Fp, on teste si E = (p,a,b) est une courbe elliptique.
    Si c'est une courbe elliptique, on enregistre son cardinal sur notre histogramme
    """
    D = {}
    lb, ub = hasseBounds(p)

    for possible in range(lb, ub + 1):
        D[possible] = 0

    for a in range(p):
        for b in range(p):
            EC = (p, a, b)
            if est_elliptique(EC):
                D[cardinal(EC)] += 1

    return D


def dessine_graphe(p):
    """Dessine le graphe de répartition des cardinaux des courbes elliptiques définies sur F_p."""
    bound = int(2 * sqrt(p))
    C = [c for c in range(p + 1 - bound, p + 1 + bound + 1)]
    D = cardinaux_courbes(p)

    plt.bar(C, [D[c] for c in C], color='b')
    plt.show()


def moins(P, p):
    """Retourne l'opposé du point P mod p."""
    x, y = P
    return x, -y % p


def est_egal(P1, P2, p):
    """Teste l'égalité de deux points mod p."""
    if est_zero(P1) and est_zero(P2): return True
    if est_zero(P1) or est_zero(P2): return False

    x1, y1 = P1
    x2, y2 = P2

    return mod(x1, p) == mod(x2, p) and mod(y1, p) == mod(y2, p)


def est_zero(P):
    """Teste si un point est égal au point à l'infini."""

    return P == ()


def addition(P1, P2, E):
    """Renvoie P1 + P2 sur la courbe E."""
    """
    Il y a 3 cas (non-trivials) possibles:
        - Si P1 = P2, la pente(=gradient) est la pente de la tangente au point P1. On utilise la differentiation pour
        l' obtenir.
        - Si P1 = -P2, P1 + P2 = point a l'infini.
        - Sinon, la pente de la ligne qui connecte P et Q est calculee normallement, avec modulo bien sur.
    Enfin, plug la pente dans la formule pour l'intersection pou trouver les coordones de P1 + P2.
    """

    if est_zero(P1) and not est_zero(P2):
        return P2
    elif est_zero(P2) and not est_zero(P1):
        return P1
    elif est_zero(P1) and est_zero(P2):
        return ()

    x1, y1 = P1
    x2, y2 = P2
    p, a, b = E

    if est_inverse_additif(P1, P2, p):
        return ()

    if est_egal(P1, P2, p):
        gradient = ((3 * exp(x1, 2, p) + a) * inv_mod(2 * y1, p)) % p

    else:
        gradient = ((y2 - y1) * inv_mod(x2 - x1, p)) % p

    x3 = (exp(gradient, 2, p) - x1 - x2) % p
    y3 = (gradient * (x1 - x3) - y1) % p

    return x3, y3


def multiplication_scalaire(k, P, E):
    """Renvoie la multiplication scalaire k*P sur la courbe E."""

    """
    On fait k additions, en se servant du code deja utilise pour l'exponentiation modulaire,
    mais on additione au lieu de multiplier.
    """
    product = ()

    if k == 0:
        return ()
    if k == 1:
        return P

    p = E[0]

    if k == -1:
        return moins(P, p)

    k_binary = toBinaire(abs(k))

    for bit in k_binary:
        product = addition(product, product, E)
        if bit == 1:
            product = addition(product, P, E)

    if k < 0:
        product = moins(product, p)

    return product


def ordre(N, factors_N, P, E):
    """Renvoie l'ordre du point P dans les points de la courbe E mod p. 
    N est le nombre de points de E sur Fp.
    factors_N est la factorisation de N en produit de facteurs premiers."""

    """ N a deja ete factorisee en factors_N donc pas besoin de le refactoriser.
    factors_N est une liste de tuples de forme (prime, exponent) dans la prime factorisation de N
    
    Explication de l'algo:
    On se sert  du theoreme de structure E = Z/cZ x Z/dZ ou c divise d et d divise p-1.
    Grace a ce theoreme, on peut se servir du Chinese remainder theoreme pour calculer l'orre.
    
         - L'ordre est au maximum N, et pour tout P, ord(P) divise N. N est donc un bon point de depart
         - Pour chaque factur (prime, exponent), on l'enleve de n. On multiplie tous les facteurs restants
          pour P pour trouver S.
         - Apres on reintroduit le prime du facteur qu'on avait enleve, mais from small to big.
          On guarde le minimum exposant qui nous ramene au point a l'infini.
          
        -On enleve donc le plus de termes(soit facteurs tout entiers, ou exposants de p non necessaires)
         possibles de N pour faire k mais en guarantissant que kP = O.
    """

    if P == ():
        return 1

    order = N

    for (prime, exponent) in factors_N:
        order = order // pow(prime, exponent)
        # Pour chaque facteur on enleve le current prime factor et son expsant

        new_point = multiplication_scalaire(order, P, E)
        # Multiplie toute la factorization restante par P

        for power in range(exponent + 1):
            # Pour le facteur qu;on a enleve, on ne garde que p ^ e, e etant le plus petit necessaire pour que
            # nous revenons au point aa l'infini.
            if new_point == ():
                break

            order = order * prime
            new_point = multiplication_scalaire(prime, new_point, E)
            # Cet exposant n'a pas marche, on ajoute p (ou ajoute + 1 sur l'exposant) sur notre soulution et on reessaie
    return order


def point_aleatoire_naif(E):
    """Renvoie un point aléatoire (différent du point à l'infini) sur la courbe E."""

    """
    Complexite worst case O(p^2) et avergae case O((p^2)/|E|). Explications donnees ci dessous.
    """
    p = E[0]
    x, y = 0, 0

    while not point_sur_courbe((x, y), E):
        x, y = randint(0, p - 1), randint(0, p - 1)
    return x, y


# E = (360040014289779780338359, 117235701958358085919867, 18575864837248358617992)
# t1_naif = time.time_ns()
# p_naif_test = point_aleatoire_naif(E)
# t2_naif = time.time_ns()
# print("Calculee dans: " + str(t2_naif-t1_naif) + " nanosecndes")

# print(p_naif_test)

""" Le test sur la courbe donnee prend trop de temps runtime (I had to terminate it). 
Ca c'est parceque on tire x et y aleatoirement sans guarantir quelconque relation
entre eux, et donc la probabilite de n'est pas etre sur la courbe est eleve.
Vu que x et y sont tire de Fp x Fp, ont peut tirer au pire des cas p^2 points.
Mais en moyenne (avergae case), on tire jusqu' a (p^2)/|E|. Le plus
|E| grandit, plus de chance on a pour trouver un (x,y) sur la courbe. """


def point_aleatoire(E):
    """Renvoie un point aléatoire (différent du point à l'infini) sur la courbe E."""

    """ Meme logique pour la complexite ici aussi.
    Dans le pire des cas, on choisit tous les valeurs de x en Fp, donc O(p) (= On assume bien sure,
    une distribution egale de la probabilite d'etre choisi).
    En average case, O(p/|E|). 
    """

    p, a, b = E

    if mod(p, 4) != 3:
        return point_aleatoire_naif(E)

    x = randint(0, p - 1)
    rhs = right(x, E)

    if symbole_legendre(rhs, p) == -1:
        return point_aleatoire(E)
    else:
        return x, racine_carree(rhs, p)


# E = (360040014289779780338359, 117235701958358085919867, 18575864837248358617992)
# t1 = time.time_ns()
# p_test = point_aleatoire(E)
# t2 = time.time_ns()
# print(p_test)
# print("Calculee dans: " + str(t2-t1) + " nanosecndes")

""" Calculee dans: 0 nanosecndes """


def point_ordre(E, N, factors_N, n):
    """Renvoie un point aléatoire d'ordre n sur la courbe E.
    Ne vérifie pas que n divise N."""

    randomPoint = point_aleatoire(E)

    if ordre(N, factors_N, randomPoint, E) == n:
        return randomPoint
    else:
        return point_ordre(E, N, factors_N, n)


def keygen_DH(P, E, n):
    """Génère une clé publique et une clé privée pour un échange Diffie-Hellman.
    P est un point d'ordre n sur la courbe E.
    """

    sec = randint(0, n - 1)
    pub = multiplication_scalaire(sec, P, E)

    return sec, pub


def echange_DH(sec_A, pub_B, E):
    """Renvoie la clé commune à l'issue d'un échange Diffie-Hellman.
    sec_A est l'entier secret d'Alice et pub_b est l'entier public de Bob."""

    return multiplication_scalaire(sec_A, pub_B, E)


""" Question 15. """


def diffie_hellman(Ec, Num, factors_Num):
    """Renvoie un point d'ordre N qui va etre utilise pour un echange Diffie-Hellman"""

    """ 
    On choisi le premier qui est le plus premier comme notre ordre de <g>. Iterate sur tous les facteurs pour le trouver.
    Le point est bien car l'ordre etant maximal rend sa tache d'attaquer notre systeme la plus dure possible pour ce N.
    Notez que c'est ne pas le meme point qui est choisi.
    """
    order = 1
    for fac in factors_Num:
        order = max(order, fac[0])

    return point_ordre(Ec, Num, factors_Num, order)


p = 248301763022729027652019747568375012323
E = (p, 1, 0)
N = 248301763022729027652019747568375012324
facteurs = [(2, 2), (62075440755682256913004936892093753081, 1)]

pointDHM = diffie_hellman(E, N, facteurs)

# print(" [+] Le bon point choisi est: ", pointDHM)

assert (ordre(N, facteurs, pointDHM, E) == facteurs[1][0])  # Sanity check here
