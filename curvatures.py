from sympy import *
import mpmath

def tangent_vectors(x, u, v):
    '''Computes tangent vectors of surface, returns x_u, x_v as tuple'''
    xu = simplify(diff(x, u))
    xv = simplify(diff(x, v))
    return xu, xv

def unit_normal(x, u, v, simple=True):
    '''Computes unit normal vector (Gauss map) of surface

    Optional argument simple (True/False) determines whether to simple expression for normal vector not.'''
    xu, xv = tangent_vectors(x, u, v)
    normal = xu.cross(xv)
    norm = normal.norm()
    if simple:
        return simplify(normal / norm)
    else:
        return normal/norm

def first_ff(x, u, v):
    '''Computes first fundamental form of surface, returns E,F,G as tuple'''
    xu, xv = tangent_vectors(x, u, v)
    E = simplify(xu.dot(xu))
    F = simplify(xu.dot(xv))
    G = simplify(xv.dot(xv))
    return E, F, G

def second_ff(x, u, v):
    '''Computes second fundamental form of surface, returns e,f,g as tuple'''
    xu, xv = tangent_vectors(x, u, v)
    N = unit_normal(x, u, v)
    xuu = simplify(diff(xu, u))
    xuv = simplify(diff(xu, v))
    xvv = simplify(diff(xv, v))

    e = simplify(N.dot(xuu))
    f = simplify(N.dot(xuv))
    g = simplify(N.dot(xvv))

    return e, f, g

def dN_p(x, u, v):
    '''Computes differential of the Gauss map in 2x2 matrix form'''
    E, F, G = first_ff(x, u, v)
    e, f, g = second_ff(x, u, v)
    fffmatrix = Matrix([[E,F],[F,G]])
    sffmatrix = Matrix([[e,f],[f,g]])
    dNtranspose = simplify(-sffmatrix * (fffmatrix**-1))
    return dNtranspose.T

def Gaussian_curvature(x, u, v):
    '''Computes Gaussian curvature of surface'''
    E, F, G = first_ff(x, u, v)
    e, f, g = second_ff(x, u, v)
    return simplify((e*g-f**2) / (E*G-F**2))

def mean_curvature(x, u, v):
    '''Computes mean curvature of surface'''
    E, F, G = first_ff(x, u, v)
    e, f, g = second_ff(x, u, v)
    return simplify((1/2)*(e*G-2*f*F+g*E) / (E*G-F**2))

def principal_curvatures(x, u, v):
    '''Computes principal curvatures of surface;
    returns two values k_1, k_2 as tuple if k_1 and k_2 are different, otherwise if k_1=k_2 then returns a single value'''
    dN = dN_p(x, u, v)
    try:
        k1, k2 = dN.eigenvals()
        return simplify(k1), simplify(k2)
    except:
        k1 = dN.eigenvals()
        return simplify(k1)

def compute_curvatures(x, u, v):
    '''All-in-one function that computes first and second fundamental forms, Gaussian and mean curvatures,
    and principal curvatures.

    This function reuses code from the individual functions, but is more efficient than running each
    of the above functions separately.

    May require more refined methods of simplification in order to get things to turn out nicely.'''

    xu, xv = tangent_vectors(x, u, v)

    E = simplify(xu.dot(xu))
    F = simplify(xu.dot(xv))
    G = simplify(xv.dot(xv))
    symE = symbols('E')
    symF = symbols('F')
    symG = symbols('G')
    print('First fundamental form:')
    display(Eq(symE, E))
    display(Eq(symF, F))
    display(Eq(symG, G))

    N = unit_normal(x, u, v)
    xuu = simplify(diff(xu, u))
    xuv = simplify(diff(xu, v))
    xvv = simplify(diff(xv, v))

    e = simplify(N.dot(xuu))
    f = simplify(N.dot(xuv))
    g = simplify(N.dot(xvv))
    syme = symbols('e')
    symf = symbols('f')
    symg = symbols('g')
    print('Second fundamental form:')
    display(Eq(syme, e))
    display(Eq(symf, f))
    display(Eq(symg, g))

    K = simplify((e*g-f**2) / (E*G-F**2))
    print('Gaussian curvature:')
    symK = symbols('K')
    display(Eq(symK, K))

    H = simplify(Rational(1,2)*(e*G-2*f*F+g*E) / (E*G-F**2))
    print('Mean curvature:')
    symH = symbols('H')
    display(Eq(symH, H))

    k1 = H + sqrt(H**2-K)
    k2 = H - sqrt(H**2-K)
    symk1 = symbols('k_1')
    symk2 = symbols('k_2')
    print('Principal curvatures:')
    if simplify(k1-k2) == 0:
        display(Eq(symk1, Eq(symk2,k1)))
    else:
        display(Eq(symk1, k1))
        display(Eq(symk2, k2))

    return E, F, G, e, f, g, K, H, k1, k2
