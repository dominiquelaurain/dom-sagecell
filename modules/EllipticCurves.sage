print "ElliptiCurves (EC) module loading ..."

# Elliptic Curves module is made to set helping functions about elliptic curves

# EC_Squared_y_Parameters returns equation coefficient after a change of var leaving a square equality
# Equation : y^2 + a1*x*y + a3*y = x^3 + a2*x^2 + a4*x + a6
# New equation after change of var :  x = u and y = 1/2*v - a1/2*u - a3/2
# v^2 = 4*u^3 + b2*u^2 + 2*b4*u + b6
def EC_Squared_y_Parameters(a1,a2,a3,a4,a6):
    b2 = a1^2 + 4*a2
    b4 = 2*a4 + a1*a3
    b6 = a3^2 + 4*a6
    b8 = a1^2*a6 + 4*a2*a6 - a1*a3*a4 + a2*a3^2 - a4^2
    return [b2,b4,b6,b8]

# EC_Parameters returns A,B for the equation  w^2 = t^3 + p*t + q
# after the change of var :  t = 4*u + b2/3   and w = 4*v
def EC_Parameters_from_bi(b2,b4,b6):
    p = -1/3*b2^2 + 8*b4
    q = 2/27*b2^3 - 8/3*b2*b4 + 16*b6
    return [p,q]

def EC_Parameters(a1,a2,a3,a4,a6):
    [b2,b4,b6,b8] = EC_Squared_y_Parameters(a1,a2,a3,a4,a6)
    [p,q] = EC_Parameters_from_bi(b2,b4,b6)
    return [p,q]

# EC_Discriminant returns discriminant of elliptic curve
# disc = - b2^2*b8 - 8*b4^3 - 27*b6^2 + 9*b2*b4*b6 # disc = -1/256*(4*p^3 + 27*q^2)
def EC_Discriminant_from_pq(p,q):
    return -1/256*(4*p^3 + 27*q^2)

def EC_Discriminant(a1,a2,a3,a4,a6):
    [p,q] = EC_Parameters(a1,a2,a3,a4,a6)
    disc = EC_Discriminant_from_pq(p,q)
    return  disc

# EC_CurveEq returns equation x^3 + p*x + q
def EC_CurveEq(p,q):
    x = var('x')
    g(x) = x^3 + p*x + q
    return g

# EC_PointsInflectionEq returns equation 3*x^4 + 6*p*x^2 + 12*q*x - p^2
def EC_PointsInflectionEq(p,q):
    x = var('x')
    i(x) = 3*x^4 + 6*p*x^2 + 12*q*x - p^2
    return i

# EC_Negate returns negated point for (x,y)
def EC_Negate(a1,a2,a3,a4,a6,x,y):
    x1 = x
    y1 = -a1*x - a3 - y
    return [x1,y1]

# EC_double returns double point (x,y) + (x,y)
def EC_Double(a1,a2,a3,a4,a6,x,y):
    k = (3*x^2 + 2*a2*x - a1*y + a4)/(2*y + a1*x + a3)
    x1 = k^2 + k*a1 - a2 -2*x
    y1 = - a1*x1 - a3 - k*x1 + k*x - y
    return [x1,y1]

# EC_add returns added point (x1,y1) + (x2,y2)
def EC_Add(a1,a2,a3,a4,a6,x1,y1,x2,y2):
    k = (y2 - y1)/(x2 - x1)
    x3 = k^2 + a1*k - a2 - x1 - x2
    y3 = -a1*x3 - a3 - k*x3 + k*x1 - y1
    return [x3,y3]



