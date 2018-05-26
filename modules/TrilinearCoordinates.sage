print "TrilinearCoordinates (TC) module loading ..."

# Trilinear coordinates are coordinates with respect to a reference triangle ABC
# Edge lengths are a = BC, b = AC, c = AB

# For ETC search (http://faculty.evansville.edu/ck6/encyclopedia/Search_13_6_9.html) we define ETC_search function
def ETC_Search(a,b,c,t):
    [ta,tb,tc] = t
    S = TriArea(a,b,c)
    k = 2*S/(a*ta + b*tb + c*tc)
    kx = k*ta
    return kx

def ETC_Search_13_6_9(t):
    a = 13 ; b = 6 ; c = 9
    kx = ETC_Search(a,b,c,t)
    return kx

def ETC_Search_6_9_13(t):
    a = 6 ; b = 9 ; c = 13
    kx = ETC_Search(a,b,c,t)
    return kx

# Evaluate list
def TriEvalRealList(le):
    return [e.n() for e in le]

# General function about trilinear, in order to evaluate real
def TriEvalReal(t):
    [ta,tb,tc] = t
    return [ta.n(),tb.n(),tc.n()]

def TriEvalRealArray(lt):
    return [TriEvalReal(t) for t in lt]

# TriEvalRealComplex returns only the real part of the complex values ta,tb,tc
def TriEvalRealComplex(t):
    [ta,tb,tc] = t
    return [ta.n().real(),tb.n().real(),tc.n().real()]

# TriSafeFactor handles exception for zero factoring
def TriSafeFactor(e):
    try:
        ev = e.factor()
    except ArithmeticError :
        ev = 0
    return ev

# General function about trilinear, in order to factorize symbolic coordinates
def TriFactor(t):
    [ta,tb,tc] = t
    ta = TriSafeFactor(ta)
    tb = TriSafeFactor(tb)
    tc = TriSafeFactor(tc)
    t = [ta,tb,tc]
    return t

# General function about trilinear, in order to gcd integer coordinates
def TriGcd(t):
    [ta,tb,tc] = t
    m = gcd(ta,gcd(tb,tc))
    return [ta/m,tb/m,tc/m]

# TriHomogeneousSimplify simplify result of equation system by setting free variables with 1
def TriHomogeneousSimplify(res,unknown_vars):
    # Compute set of free vars
    free_vars = []
    for sol in res:
        for sdef in sol:
            for v in sdef.rhs().variables():
                if v not in free_vars and v not in unknown_vars:
                    # select only free variables starting with 'r'
                    v_str = latex(v)
                    if v_str[0] == 'r': #
                        free_vars.append(v)
    # Substitue value of free vars with 1
    #print "free_vars = ",free_vars
    res1 = []
    for sol in res:
        sol1 = []
        for sdef in sol:
            for v in free_vars:
                sdef  = sdef.subs({v : 1})
            sol1.append(sdef)
        res1.append(sol1)
    #
    return res1

# General functions in order to compute a tangential quadrilateral
#  TriHalfCoTangent(x) returns cotg(alpha/2) is x = cotg(alpha)
#  TriDoubleCoTangent(x) returns cotg(2*alpha) is x = cotg(alpha)
#  TriTangentialQuadrilateralLastCotangent returns second "swing" angle cotangent for triangle constructed on BC
#  edge of reference triangle ABC
#    x1 = cotg(B) ; x2 = cotg(C) ; x3 = ct = cotangent first "swing angle"
def TriHalfCoTangent(x):
    return 1/( -x + sqrt(x^2 + 1))

def TriDoubleCoTangent(y):
    return 1/( 2*(1/y)/(1 - (1/y)^2))

def TriTangentialQuadrilateralLastCotangent(x1,x2,x3):
    y1 = TriHalfCoTangent(ctB)  # x1 = cotg(B/2)
    y2 = TriHalfCoTangent(ctC)  # x2 = cotg(C/2)
    y3 = TriHalfCoTangent(ct)   # x3 = ct
    # Iosifescu's characteristics for tangential quadrilateral
    # y1*y4 = y2*y3 => y4 = (y2*y3)/y1
    y4 = (y2*y3)/y1
    x4 = TriDoubleCoTangent(y4)
    return x4

# TriCuttingRatios returns points D,E,F on edges cutting in ratio d : AD/AB = d, BE/BC = d, CF/CA = d
def TriCuttingRatios(a,b,c,d):
    t1 = [(1-d)*b,d*a,0]
    t2 = [0,(1-d)*c,d*b]
    t3 = [d*c,0,(1-d)*a]
    return [t1,t2,t3]

# TriBarycentric returns barycentric coordinates of point P given by trilinear coordinates x:y:z with respect to
# triangle ABC.
# With arbitrary origin and vector notation P = ax/(ax+by+cz) A + by/(ax+by+cz) B + cz/(ax+by+cz) C
def TriBarycentric(a,b,c,t):
    [ta,tb,tc] = t
    n = a*ta + b*tb + c*tc
    return [a*ta/n,b*tb/n,c*tc/n]

# TriBarycentric2ArealCoordinates convert barycentric coordinates into areal coordinates (areas of sub-triangles) by multiplying them with
# the area S of the reference triangle
def TriBarycentric2ArealCoordinates(S,m):
    [ma,mb,mc] = m
    return [S*ma,S*mb,S*mc]

# TriAreal2BarycentricCoordinates convert areal coordinates into barycentric coordinates by dividing them with
# the area S of the reference triangle
def TriAreal2BarycentricCoordinates(S,u):
    [ua,ub,uc] = u
    return [ua/S,ub/S,uc/S]

# TriLinear2BarycentricCoordinates returns homogeneous barycentric coordinates (ma,mb) for point M = ma A + mb B + (1-ma-mb) C given by trilinear coordinates
def TriLinear2BarycentricCoordinates(a,b,c,t):
    [ma,mb,mc] = TriBarycentric(a,b,c,t)
    return [ma,mb]

# TriBarycentric2LinearCoordinates returns trilinear coordinates reference triangle ABC from homogeneous barycentric coordinates A (ma) and B (mb) and C (1-ma-mb)
# M = ma A + mb B + (1-ma-mb) C
# Example : m = [0,7/12]) for M = 0 A + (7/12) B + (5/12) C or BI:IC = 5:7
def TriBarycentric2LinearCoordinates(a,b,c,m):
    [ma,mb] = m
    ta = ma/a ; tb = mb/b ; tc = (1-ma-mb)/c
    return [ta,tb,tc]

# Convert cartesian coordinates to barycentric coordinates
def TriCartesian2Barycentric(p1,p2,p3,p4):
    [x1,y1] = p1 ; [x2,y2] = p2 ; [x3,y3] = p3 ; [x4,y4] = p4
    l1 = ( (y2 - y3)*(x4 - x3) + (x3 - x2)*(y4 - y3) ) / ( (y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3) )
    l2 = ( (y3 - y1)*(x4 - x3) + (x1 - x3)*(y4 - y3) ) / ( (y2 - y3)*(x1 - x3) + (x3 - x2)*(y1 - y3) )
    l3 = 1 - l1 - l2
    return [l1,l2,l3]

# Convert barycentric coordinates to trilinear coordinates
def TriBarycentric2Trilinear(a,b,c,l):
    [la,lb,lc] = l
    return [la/a,lb/b,lc/c]

# TriBarycentric2Points returns trilinear coordinates of barycentric of two points
def TriBarycentric2Points(a,b,c,w1,w2,t1,t2):
    w = w1 + w2
    k1 = TriBarycentric(a,b,c,t1)
    k2 = TriBarycentric(a,b,c,t2)
    k3 = [w1/w*k1[0]+w2/w*k2[0],w1/w*k1[1]+w2/w*k2[1],w1/w*k1[2]+w2/w*k2[2]]
    t3 = TriBarycentric2Trilinear(a,b,c,k3)
    return t3

# TriBarycentric3Points returns trilinear coordinates of barycentric of three points
def TriBarycentric3Points(a,b,c,w1,w2,w3,t1,t2,t3):
    w = w1 + w2 + w3
    k1 = TriBarycentric(a,b,c,t1)
    k2 = TriBarycentric(a,b,c,t2)
    k3 = TriBarycentric(a,b,c,t3)
    k4 = [w1/w*k1[0]+w2/w*k2[0]+w3/w*k3[0],w1/w*k1[1]+w2/w*k2[1]+w3/w*k3[1],w1/w*k1[2]+w2/w*k2[2]+w3/w*k3[2]]
    t4 = TriBarycentric2Trilinear(a,b,c,k4)
    return t4

# OLD : check it and if not used, delet function
# TriReflection returns reflection of point
#def TriReflection(t1,t2):
#    t3 = TriBarycentric2Points(a,b,c,1,-2,t2,t1)
#    return t3

# TriHarmonicConjugates returns the harmonic conjugates t3,t4 of two points t1,t2
# so they are on line t1t2 and such as t1t3/t2t3 = t2t4/t1t4
def TriHarmonicConjugates(t1,t2):
    [ta1,tb1,tc1] = t1; [ta2,tb2,tc2] = t2
    ta3 = ta1 + ta2; tb3 = tb1 + tb2; tc3 = tc1 + tc2; t3 = [ta3,tb3,tc3]
    ta4 = ta1 - ta2; tb4 = tb1 - tb2; tc4 = tc1 - tc2; t4 = [ta4,tb4,tc4]
    return [t3,t4]

# Convert cartesian coordinates to trilinear coordinates
def TriCartesian2Trilinear(a,b,c,p1,p2,p3,p4):
    l = TriCartesian2Barycentric(p1,p2,p3,p4)
    t = TriBarycentric2Trilinear(a,b,c,l)
    return t

# Remove one point from list
def TriRemoveCopy(pA,lp):
    lp1 = []
    for v in lp:
        p = (v[0].right().n(),v[1].right().n())
        if RT_Quadrance(p,pA) > 0.0001:
            lp1.append(p)
    return lp1

# TriSquaredArea returns squared of area of reference triangle, using Heron formula
def TriSquaredArea(a,b,c):
    return (a+b+c)*(-a+b+c)*(a-b+c)*(a+b-c)/16

# TriSquaredAreaFromQuadrances returns squared of area of reference triangle, using Heron formula
def TriSquaredAreaFromQuadrances(qa,qb,qc):
    return 1/16*((qa + qb + qc)^2 - 2*(qa^2 + qb^2 + qc^2))

# TriArea returns area of reference triangle, using Heron formula
def TriArea(a,b,c):
    S2 = TriSquaredArea(a,b,c)
    return sqrt(S2)

# TriSquaredDiameter return square of diameter D = abc/2S where S is area
def TriSquaredDiameter(a,b,c):
    S2 = TriSquaredArea(a,b,c)
    d2 = (a*b*c/2)^2/S2
    return d2

# TriDiameter returns diameter or circumcircle
def TriDiameter(a,b,c):
    d2 = TriSquaredDiameter(a,b,c)
    return sqrt(d2)

# TriSquaredCircumradius return square of circumradius
def TriSquaredCircumradius(a,b,c):
    return TriSquaredDiameter(a,b,c) / 4

# TriDiameter returns circumradius
def TriCircumradius(a,b,c):
    d = TriDiameter(a,b,c)
    R = d/2
    return R

# TriCircumradiusFromS returns circumradius R from formula 2*R = abc/2S
def TriCircumradiusFromS(a,b,c):
    return a*b*c/(4*S)


# TriInradius returns inradius
def TriInradius(a,b,c):
    S = RT_TriangleArea(a,b,c)
    s = RT_TriangleSemiperimeter(a,b,c)
    r = RT_TriangleInradius(S,s)
    return r

# TriEulerFormula returns quadrance OI^2 = d^2 between circumcenter I and incenter O because from Euler's formula
#  d^2 = R(R-2r) or  1/(R+d) + 1/(R-d) = 1/r   (r as an harmonic mean)
def TriEulerFormula(a,b,c):
    R = TriCircumradius(a,b,c)
    r = TriInradius(a,b,c)
    q = R*(R-2*r)
    return q

# TriConwayParameters return Conway parameters for reference triangle
# S_A = 2S cotg A = bc cos A ; S_B = 2S cotg B = ca cos B ; S_C = 2S cotg C = ab cos C
# where S is area triangle ABC
def TriConwayParameters(a,b,c):
    S_A = (-a^2 + b^2 + c^2)/2 ; S_B = (a^2 - b^2 + c^2)/2 ; S_C = (a^2 + b^2 - c^2)/2
    return [S_A,S_B,S_C]

def TriConwayParametersFromQuadrances(qa,qb,qc):
    S_A = (-qa + qb + qc)/2 ; S_B = (qa - qb + qc)/2 ; S_C = (qa + qb - qc)/2
    return [S_A,S_B,S_C]


# TriConwayQuadrances return quadrances for edge (squared edge lenghts)
def TriConwayQuadrances(S_A,S_B,S_C):
    qa = S_B + S_C ; qb = S_C + S_A ; qc = S_A + S_B
    return [qa,qb,qc]

# TriConwayTriangleNotation return S(theta) = 2S cotg(theta)  where S is area of
# triangle ABC with side lengths a,b,c
# Reference https://en.wikipedia.org/wiki/Conway_triangle_notation
#  ctx = cotg(theta) where theta is angle at vertex X
#  S_X is Conway parameter at vertex X
def TriConwayTriangleNotation(S,ctx):
    S_X = 2*S*ctx
    return S_X

# TriConwayCotangents2Quadrances returns quadrances edges lengths for triangle area S = 1 given
#  u = cotg(A)  ; v = cotg(B) and w = (1 - uv)/(u + v) = cotg(C)
def TriConwayCotangents2Quadrances(u,v,w):
    S = 1
    S_A = TriConwayTriangleNotation(S,u)
    S_B = TriConwayTriangleNotation(S,v)
    S_C = TriConwayTriangleNotation(S,w)
    [qa,qb,qc] = TriConwayQuadrances(S_A,S_B,S_C)
    return [qa,qb,qc]


# TriConwayTriangleCotangent returns ct = cotg(theta) = S(theta)/(2*S) where S is area of
# triangle ABC with side lengths a,b,c
# It's inverse function for TriConwayTriangleNotation
def TriConwayTriangleCotangent(S,S_X):
    return S_X/(2*S)

# TriConwayTriangleCotangents returns cotangents internal angles
def TriConwayTriangleCotangents(a,b,c,S):
    [S_A,S_B,S_C] = TriConwayParameters(a,b,c)
    ctA = TriConwayTriangleCotangent(S,S_A)
    ctB = TriConwayTriangleCotangent(S,S_B)
    ctC = TriConwayTriangleCotangent(S,S_C)
    return [ctA,ctB,ctC]

# TriConwayQuadrances2Cotangents returns cotangents internal angles for any triangle
def TriConwayQuadrances2Cotangents(qa,qb,qc):
    S2 = 1/16*( (qa + qb + qc)^2 - 2*(qa^2 + qb^2 + qc^2) ) ; S = sqrt(S2) # QuadArea formula
    [S_A,S_B,S_C] = TriConwayParametersFromQuadrances(qa,qb,qc)
    ctA = TriConwayTriangleCotangent(S,S_A)
    ctB = TriConwayTriangleCotangent(S,S_B)
    ctC = TriConwayTriangleCotangent(S,S_C)
    return [ctA,ctB,ctC]

# TriConwayQuadrancesArea2Cotangents returns cotangents internal angles for any triangle
def TriConwayQuadrancesArea2Cotangents(S,qa,qb,qc):
    [S_A,S_B,S_C] = TriConwayParametersFromQuadrances(qa,qb,qc)
    ctA = TriConwayTriangleCotangent(S,S_A)
    ctB = TriConwayTriangleCotangent(S,S_B)
    ctC = TriConwayTriangleCotangent(S,S_C)
    return [ctA,ctB,ctC]

# TriConwayTriangleTangents returns tangents internal angles
def TriConwayTriangleTangents(a,b,c,S):
    [ctA,ctB,ctC] = TriConwayTriangleCotangents(a,b,c,S)
    return [1/ctA,1/ctB,1/ctC]

# TriConwayFormula_A returns trilinear coordinates of a point P apex of triangle BPC
# erected on edge BC
# CBP is angle theta and PCB is angle phi  : "swings" angles
def TriConwayFormula_A(a,b,c,S,ct,cu):
    [S_A,S_B,S_C] = TriConwayParameters(a,b,c)
    S_t = TriConwayTriangleNotation(S,ct) ; S_u = TriConwayTriangleNotation(S,cu)
    return [(-a^2)/a,(S_C+S_u)/b,(S_B+S_t)/c]

# TriConwayFormula_B returns trilinear coordinates of a point P apex of triangle CPA
# erected on edge CA
# ACP is angle theta and PAC is angle phi  : "swings" angles
def TriConwayFormula_B(a,b,c,S,ct,cu):
    [S_A,S_B,S_C] = TriConwayParameters(a,b,c)
    S_t = TriConwayTriangleNotation(S,ct) ; S_u = TriConwayTriangleNotation(S,cu)
    return [(S_C+S_t)/a,(-b^2)/b,(S_A+S_u)/c]

# TriConwayFormula_C returns trilinear coordinates of a point P apex of triangle APB
# erected on edge AB
# BAP is angle theta and PBA is angle phi  : "swings" angles
def TriConwayFormula_C(a,b,c,S,ct,cu):
    [S_A,S_B,S_C] = TriConwayParameters(a,b,c)
    S_t = TriConwayTriangleNotation(S,ct) ; S_u = TriConwayTriangleNotation(S,cu)
    return [(S_B+S_u)/a,(S_A+S_t)/b,(-c^2)/c]

# TriConwayFoot_A return foot of A on BC
def TriConwayFoot_A(a,b,c):
    return [0, (a^2 + b^2 - c^2)*c, (a^2 - b^2 + c^2)*b]

# TriConwayFoot_B return foot of B on CA
def TriConwayFoot_B(a,b,c):
    return [(a^2 + b^2 - c^2)*c, 0,(- a^2 + b^2 + c^2)*a]

# TriConwayFoot_C return foot of C on AB
def TriConwayFoot_C(a,b,c):
    return [(a^2 - b^2 + c^2)*b, (- a^2 + b^2 + c^2)*a, 0]

# TriConwayBrocard return S(w) where w is Brocard angle
# cotg(w) = cotg(alpha) + cotg(beta) + cotg(gamma) where alpha,beta,gamma are internal
def TriConwayBrocard(a,b,c):
    S_w = (a^2 + b^2 + c^2)/2
    return S_w

# TriPolar2Barycentric returns barycentric coordinates of two points De with e = +1 or -1 such as ADe:BDe:CDe = u:v:w
def TriPolar2Barycentric(a,b,c,u,v,w):
    # Computing for triangle A,B,C
    S = RT_TriangleArea(a,b,c)
    [ctA,ctB,ctC] = TriConwayTriangleCotangents(a,b,c,S)
    # Computing for triangle edge lengths au, bv, cw
    S1 = RT_TriangleArea(a*u,b*v,c*w)
    [ctA1,ctB1,ctC1] = TriConwayTriangleCotangents(a*u,b*v,c*w,S1)
    # Computing barycentric coordinates of two points whose middle is circumcenter o
    [ma1,mb1,mc1] = [a^2*(ctA - ctA1),b^2*(ctB - ctB1),c^2*(ctC - ctC1)] ; m1 = ma1 + mb1 + mc1
    [ma2,mb2,mc2] = [a^2*(ctA + ctA1),b^2*(ctB + ctB1),c^2*(ctC + ctC1)] ; m2 = ma2 + mb2 + mc2
    # Returning them
    return [ [ma1/m1,mb1/m1], [ma2/m2,mb2/m2] ]

# TriPolar2Midpoint returns trilinear coordinates of midpoint D-1 and D1
def TriPolar2Midpoint(a,b,c,u,v,w):
    ta = ( 2*a^2*v^2*w^2 + b^2*v^4 + c^2*w^4  - ( a^2*u^2*(v^2 + w^2) + (b^2 + c^2)*v^2*w^2 ) )*b*c
    tb = ( 2*b^2*u^2*w^2 + c^2*w^4 + a^2*u^4  - ( b^2*v^2*(w^2 + u^2) + (c^2 + a^2)*w^2*u^2 ) )*c*a
    tc = ( 2*c^2*u^2*v^2 + a^2*u^4 + b^2*v^4  - ( c^2*w^2*(u^2 + v^2) + (a^2 + b^2)*u^2*v^2 ) )*a*b
    return [ta,tb,tc]

# TriEuler_Tripolar_ functions returns quadrance from A,B or C for a point P given the other quadrances
#
# the functions come from the Euler equation for a reference triangle ABC and for tripolar coordinates of a point P
# (qa + qb - qc)*(qx*qy + qc*qz) + (- qa + qb + qc)*(qy*qz + qa*qx) + (qa - qb + qc)*(qz*qx + qb*qy) - (qa*qx^2 + qb*qy^2 + qc*qz^2) - qa*qb*qc = 0
#
# there are two quadrances solutions as there are two points as intersection of two circles centered at two vertices of ABC

def TriEuler_Tripolar_A(qa,qb,qc,qy,qz):
    disc = (qa^2 + qb^2 + qc^2 - 2*qa*qb - 2*qa*qc - 2*qb*qc)*(qa^2 - 2*qa*qy + qy^2 - 2*qa*qz - 2*qy*qz + qz^2)
    k = qa^2 - qa*qb - qa*qc - qa*qy - qb*qy + qc*qy - qa*qz + qb*qz - qc*qz
    qx1 = -1/2*(k - sqrt(disc))/qa ; qx2 = -1/2*(k + sqrt(disc))/qa
    return [qx1,qx2]

def TriEuler_Tripolar_B(qa,qb,qc,qx,qz):
    disc = (qa^2 + qb^2 + qc^2 - 2*qa*qb - 2*qa*qc - 2*qb*qc)*(qb^2 - 2*qb*qx + qx^2 - 2*qb*qz - 2*qx*qz + qz^2)
    k =  qb^2 - qb*qc - qb*qa - qb*qz - qc*qz + qa*qz - qb*qx + qc*qx - qa*qx
    qy1 = -1/2*(k - sqrt(disc))/qb ; qy2 = -1/2*(k + sqrt(disc))/qb
    return [qy1,qy2]

def TriEuler_Tripolar_C(qa,qb,qc,qx,qy):
    disc = (qa^2 + qb^2 + qc^2 - 2*qa*qb - 2*qa*qc - 2*qb*qc)*(qc^2 - 2*qc*qx + qx^2 - 2*qc*qy - 2*qx*qy + qy^2)
    k =  qc^2 - qc*qa - qc*qb  - qc*qx - qa*qx + qb*qx  - qc*qy + qa*qy - qb*qy
    qz1 = -1/2*(k - sqrt(disc))/qc ; qz2 = -1/2*(k + sqrt(disc))/qc
    return [qz1,qz2]

# TriBrocardPoint return the (first) Brocard point of triange ABC
# Brocard point is point P inside triangle such as cevian angles are same
# NOTE : Brocard point is NOT a triangle center, because coordinates are not cyclic in term
#        of edge lengths a,b,c
# NOTE : second Brocard point is Brocard point of triangle ACB
def TriBrocardPoint(a,b,c):
    return [c/b,a/c,b/a]

def TriSecondBrocardPoint(a,b,c):
    return [b/c,c/a,a/b]

# TriAcuteness returns sign of conway parameters , if positive or zero angle is acute
def TriAcuteness(a,b,c):
    [S_A,S_B,S_C] = TriConwayParameters(a,b,c)
    return [sign(S_A),sign(S_B),sign(S_C)]

# TriIsAcuteCheck_X check whether angle at X is acute
def TriAcuteCheck_A(a,b,c):
    [S_A,S_B,S_C] = TriConwayParameters(a,b,c)
    return (S_A >= 0)

def TriAcuteCheck_B(a,b,c):
    [S_A,S_B,S_C] = TriConwayParameters(a,b,c)
    return (S_B >= 0)

def TriAcuteCheck_C(a,b,c):
    [S_A,S_B,S_C] = TriConwayParameters(a,b,c)
    return (S_C >= 0)

# TriIsRectangularCheck_X check whether angle at X is rectangular (pi/2)
def TriRectangularCheck_A(a,b,c):
    [S_A,S_B,S_C] = TriConwayParameters(a,b,c)
    return (S_A == 0)

def TriRectangularCheck_B(a,b,c):
    [S_A,S_B,S_C] = TriConwayParameters(a,b,c)
    return (S_B == 0)

def TriRectangularCheck_C(a,b,c):
    [S_A,S_B,S_C] = TriConwayParameters(a,b,c)
    return (S_C == 0)

# TriCosines returns the three cosines of internal angles at A,B,C for triangle ABC with a = BC, b = AC, c = AB
def TriCosines(a,b,c):
    [S_A,S_B,S_C] = TriConwayParameters(a,b,c)
    cs_A = S_A/(b*c) ; cs_B = S_B/(a*c) ; cs_C = S_C/(a*b)
    return [cs_A,cs_B,cs_C]

# TriSecants returns the three secants of internal angles at A,B,C for triangle ABC with a = BC, b = AC, c = AB
# secant is inverse of cosine
def TriSecants(a,b,c):
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    sec_A = 1/cs_A ; sec_B = 1/cs_B ; sec_C = 1/cs_C
    return [sec_A,sec_B,sec_C]

# TriCosecants returns the three cosecants of internal angles at A,B,C for triangle ABC with a = BC, b = AC, c = AB
# cosecant is inverse of sine
def TriCosecants(a,b,c):
    [sn_A,sn_B,sn_C] = TriTriangleSines(a,b,c)
    csc_A = 1/sn_A ; csc_B = 1/sn_B ; csc_C = 1/sn_C
    return [csc_A,csc_B,csc_C]


# TriSquaredCosinesBisector returns cos^2((B-C)/2), cos^2((C-A)/2), cos^2((A-B)/2)
def TriSquaredCosinesBisector(a,b,c):
    d2 = TriSquaredDiameter(a,b,c)
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    sq_A = (1 + cs_B*cs_C + b*c/d2)/2
    sq_B = (1 + cs_C*cs_A + c*a/d2)/2
    sq_C = (1 + cs_A*cs_B + a*b/d2)/2
    return [sq_A,sq_B,sq_C]

# TriExactCoordinates returns the exact (=actual distance) coordinates x':y',z' for trilinear coordinates x:y:z
def TriExactCoordinates(a,b,c,t):
    S = RT_TriangleArea(a,b,c)
    [ta,tb,tc] = t
    k = (2*S)/(a*ta+b*tb+c*tc)
    ea = k*ta ; eb = k*tb ; ec = k*tc
    return [ea,eb,ec]

# TriReflectionPoint returns reflection point P' of point P in point M with respect to triangle ABC
# P' is such as vector PM = MP  ; M is midpoint of PP'
def TriReflectionPoint(a,b,c,t1,t2):
    [ta1,tb1,tc1] = t1; [ta2,tb2,tc2] = t2
    rta = 2*ta1*(b*tb2+c*tc2)+ta2*(a*ta1-b*tb1-c*tc1)
    rtb = 2*tb1*(a*ta2+c*tc2)+tb2*(-a*ta1+b*tb1-c*tc1)
    rtc = 2*tc1*(a*ta2+b*tb2)+tc2*(-a*ta1-b*tb1+c*tc1)
    return [rta,rtb,rtc]

# TriComplementPoint returns complement point P' of point P with respect to triangle ABC
# P' is such as vector PG = 2 GP'  ; G centroid is barycentric of (P,1) and (P',2)
def TriComplementPoint(a,b,c,t):
    [ta,tb,tc] = t
    cta = (b*tb + c*tc)/a ; ctb = (a*ta + c*tc)/b ; ctc = (a*ta + b*tb)/c
    return [cta,ctb,ctc]

# TriIsogonalConjugate returns point at intersection of reflected cevians of a point P with respect
# to the respective internal angle bissectors
def TriIsogonalConjugatePoint(a,b,c,t):
    [ta,tb,tc] = t
    return [1/ta,1/tb,1/tc]

# TriHaimovPoint returns Haimov point for a point P
# If QaQbQc is Haimov triangle of point P, then AQa,BQb, CQc intersect at Haimov point because
# it is perspector of ABC and QaQbQc triangles
def TriHaimovPoint(a,b,c,t):
    t1 = TriComplementPoint(a,b,c,t)
    t2 = TriIsogonalConjugatePoint(a,b,c,t1)
    return t2

# TriAntiComplementPoint returns anticomplement point P' of point P with respect to triangle ABC
# P' is such as vector P'G = 2 GP  ; G centroid is barycentric of (P',1) and (P,2)
def TriAntiComplementPoint(a,b,c,t):
    [ta,tb,tc] = t
    acta = (-a*ta + b*tb + c*tc)/a ; actb = (a*ta -b*tb + c*tc)/b ; actc = (a*ta + b*tb - c*tc)/c
    return [acta,actb,actc]

# TriSines returns the three sines of internal angles at A,B,C for triangle ABC with a = BC, b = AC, c = AB
# and r is inradius, s = (a+b+c)/2 is semiperimeter and R = abc / 4rs is circumradius
# law of sines : a/ sinA = b / sinB = c / sinD = 2R
def TriSines(a,b,c,r):
    s = (a+b+c)/2
    R = a*b*c/(4*r*s)
    sn_A = a / (2*R) ; sn_B = b / (2*R) ; sn_C = c / (2*R)
    return [sn_A,sn_B,sn_C]

# TriTriangleSinesFromR returns angles from law of sines given circumradius R
def TriTriangleSinesFromR(a,b,c,R):
    sn_A = a / (2*R) ; sn_B = b / (2*R) ; sn_C = c / (2*R)
    return [sn_A,sn_B,sn_C]

# TriTriangleSines returns angles using Heron formula to get circumradius
def TriTriangleSines(a,b,c):
    s = (a+b+c)/2
    R = a*b*c/(4*sqrt(s*(s-a)*(s-b)*(s-c)))
    return TriTriangleSinesFromR(a,b,c,R)


# TriTriangleQuadranceTwoPoints returns quadrance for two points
def TriTriangleQuadranceTwoPoints(a,b,c,S,cs_A,cs_B,cs_C,t1,t2):
    [ta1,tb1,tc1] = t1 ; [ta2,tb2,tc2] = t2
    k1 = (2*S)/(a*ta1+b*tb1+c*tc1) ; k2 = (2*S)/(a*ta2+b*tb2+c*tc2)
    q = ((k1*ta1 - k2*ta2)^2 + (k1*tb1 - k2*tb2)^2 + 2*(k1*ta1 - k2*ta2)*(k1*tb1 - k2*tb2)*cs_C)/(1 - cs_C^2)
    return q

# TriQuadranceTwoPoints returns quadrance between two points
def TriQuadranceTwoPoints(a,b,c,t1,t2):
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    S = RT_TriangleArea(a,b,c)
    q = TriTriangleQuadranceTwoPoints(a,b,c,S,cs_A,cs_B,cs_C,t1,t2)
    return q


# TriSquaredCrossRatioThreePoints returns squared cross-ratio for three points given by trilinear coordinates
def TriSquaredCrossRatioThreePoints(a,b,c,t1,t2,t3):
    q13 = TriQuadranceTwoPoints(a,b,c,t1,t3)
    q23 = TriQuadranceTwoPoints(a,b,c,t2,t3)
    return (q13/q23)

# TriSquaredCrossRatioFourPoints returns squared cross-ratio for four points given by trilinear coordinates
def TriSquaredCrossRatioFourPoints(a,b,c,t1,t2,t3,t4):
    q13 = TriQuadranceTwoPoints(a,b,c,t1,t3)
    q23 = TriQuadranceTwoPoints(a,b,c,t2,t3)
    q14 = TriQuadranceTwoPoints(a,b,c,t1,t4)
    q24 = TriQuadranceTwoPoints(a,b,c,t2,t4)
    return (q13/q23)/(q14/q24)

# SquaredCrossRatioThreePoints returns squared cross-ratio for three points given by homogenous barycentric coordinates xi:yi:1-xi-yi
def SquaredCrossRatioThreePoints(a,b,c,x1,y1,x2,y2,x3,y3):
    t1 = [x1/a,y1/b,(1-x1-y1)/c]
    t2 = [x2/a,y2/b,(1-x2-y2)/c]
    t3 = [x3/a,y3/b,(1-x3-y3)/c]
    return TriSquaredCrossRatioThreePoints(a,b,c,t1,t2,t3)

# SquaredCrossRatioFourPoints returns cross-ratio for four points given by homogenous barycentric coordinates xi:yi:1-xi-yi
def SquaredCrossRatioFourPoints(a,b,c,x1,y1,x2,y2,x3,y3,x4,y4):
    t1 = [x1/a,y1/b,(1-x1-y1)/c]
    t2 = [x2/a,y2/b,(1-x2-y2)/c]
    t3 = [x3/a,y3/b,(1-x3-y3)/c]
    t4 = [x4/a,y4/b,(1-x4-y4)/c]
    return TriSquaredCrossRatioFourPoints(a,b,c,t1,t2,t3,t4)

# Functions Tri2Cart_x and Tri2Cart_y map trilinear coordinates (ta,tb,tc) to cartesian coordinates with respect
# to a triangle ABC with edge lengths a = BC b = AC c = AB, inradius r and semi-perimeter s
# Triangle has incenter I at origin (x,y) = (0,0) and BC edge segment y = r
# BE CAREFUL : some trilinear coordinates leads to a zero divide because a*ta + b*tb + c*tc near 0
def Tri2Cart_x(a,b,c,r,s,ta,tb,tc):
    return 1/2*(a*b*ta - b^2*ta - a*c*ta + c^2*ta + 2*b^2*tb - 2*b*s*tb - 2*c^2*tc + 2*c*s*tc)/(a*ta + b*tb + c*tc)

def Tri2Cart_y(a,b,c,r,s,ta,tb,tc):
    return (b*ta + c*ta - b*tb - c*tc)*r/(a*ta + b*tb + c*tc)

# point T with trilinear coordinates [ta,tb,tc] has barycentric coordinates [a*ta,b*tb,c*tc]
# I = [0,0]
# A = [(a-b-c)*(b-c)/(2*a) , r*(b+c)/a ]
# B = [-(s-b),-r]
# C = [s-c,-r]; pC = vector(C)
# IT = (a*ta IA + b*tb IB + c*tc IC)/(a*ta + b*tb + c*tc)

# Tri2Cart functions return a point as cartesian coordinates given t trilinear coordinates and parameters of reference triangle
def Tri2Cart(a,b,c,r,s,t):
    [ta,tb,tc] = t
    x = Tri2Cart_x(a,b,c,r,s,ta,tb,tc) ; y = Tri2Cart_y(a,b,c,r,s,ta,tb,tc)
    return [x,y]

# Function TriToCart returns (x,y) cartesian coordinates from trilinear coordinates t
def TriToCart(a,b,c,r,t):
    [ta,tb,tc] = t
    if ta != 0:
        U = tb/ta ; V = tc/ta
        x = 1/2*( (b - c)*(a - b - c) + b*(- a  + b - c)*U + c*(a + b - c)*V) /(a + b*U + c*V)
        y = ((b + c) - b*U  - c*V)*r/(a + b*U + c*V)
    else:
        if tb != 0:
            U = ta/tb ; V = tc/tb
            x = 1/2*( (b - c)*(a - b - c)*U+ b*(- a  + b - c) + c*(a + b - c)*V) /(a*U + b + c*V)
            y = ((b + c)*U - b  - c*V)*r/(a*U + b + c*V)
        else:
            U = ta/tc ; V = tb/tc
            x = 1/2*( (b - c)*(a - b - c)*U + b*(- a  + b - c)*V + c*(a + b - c)) /(a*U + b*V + c)
            y = ((b + c)*U - b*V  - c)*r/(a*U + b*V + c)
    return [x,y]

# Function TriFromCart returns trilinear coordinates from cartesian coordinates (x,y)
# inverse function of TriToCart
def TriFromCart(a,b,c,r,x,y):
    ta = 2*a*b*c*(y + r)
    tb = (- 2*(a + b + c)*r*x - (a^2 + b^2 - c^2)*y + 2*a*b*r)*c
    tc = (  2*(a + b + c)*r*x - (a^2 - b^2 + c^2)*y + 2*a*c*r)*b
    t = [ta,tb,tc]
    return t

# Functions Multisection_ta, _tb and _tc retrieves trilinear coordinate of the firts point among n on segment given
# by trilinear coordinates ta1:tb1:tc1 and ta2:tb2:tc2 of end points
def Multisection_ta(n,a,b,c,ta1,tb1,tc1,ta2,tb2,tc2):
    return n*a*ta1*ta2 + b*( ta2*tb1 + (n-1)*ta1*tb2 ) + c*( ta2*tc1 + (n-1)*ta1*tc2 )

def Multisection_tb(n,a,b,c,ta1,tb1,tc1,ta2,tb2,tc2):
    return a*( ta1*tb2 + (n-1)*ta2*tb1 ) + n*b*tb1*tb2 + c*( tb2*tc1 + (n-1)*tb1*tc2 )

def Multisection_tc(n,a,b,c,ta1,tb1,tc1,ta2,tb2,tc2):
    return a*( ta1*tc2 + (n-1)*ta2*tc1 ) + b*( tb1*tc2 + (n-1)*tb2*tc1 ) + n*c*tc1*tc2

# Multisection returns trilinear coordinates as vector for first point among n on segment
def RT_Multisection(n,a,b,c,t1,t2):
    [ta1,tb1,tc1] = t1 ; [ta2,tb2,tc2] = t2
    ta = Multisection_ta(n,a,b,c,ta1,tb1,tc1,ta2,tb2,tc2)
    tb = Multisection_tb(n,a,b,c,ta1,tb1,tc1,ta2,tb2,tc2)
    tc = Multisection_tc(n,a,b,c,ta1,tb1,tc1,ta2,tb2,tc2)
    t = [ta,tb,tc]
    return t

# ItemMultisection returns the k-th point among n on segment
def RT_ItemMultisection(k,n,a,b,c,t1,t2):
    for i in range(0,k):
        t = RT_Multisection(n-i,a,b,c,t1,t2)
        t1 = t
    return t

# Multisection returns trilinear coordinates as vector for first point among n on segment
def Multisection(n,a,b,c,t1,t2):
    [ta1,tb1,tc1] = t1 ; [ta2,tb2,tc2] = t2
    ta = Multisection_ta(n,a,b,c,ta1,tb1,tc1,ta2,tb2,tc2)
    tb = Multisection_tb(n,a,b,c,ta1,tb1,tc1,ta2,tb2,tc2)
    tc = Multisection_tc(n,a,b,c,ta1,tb1,tc1,ta2,tb2,tc2)
    t = TriEvalReal([ta,tb,tc])
    return t

# ItemMultisection returns the k-th point among n on segment
def ItemMultisection(k,n,a,b,c,t1,t2):
    for i in range(0,k):
        t = Multisection(n-i,a,b,c,t1,t2)
        t1 = t
    return t


# TriLine defines a line as the "orthogonal" vector (see TriOnLine)
# Line given in trilinear coordinates : determinant is 0 given two points A1 and A2
#  |  ta1 tb1 tc1 |     point A1 = ta1:tb1:tc1
#  |  ta2 tb2 tc2 |     point A2 = ta2:tb2:tc2
#  |  ta  tb  tc  |     any point M = ta:tb:tc on the line
# (tb1 tc2 - tc1 tb2)ta + (tc1 ta2 - ta1 tc2)tb + (ta1 tb2 - tb1 ta2)tc = 0
def TriLine(a,b,c,t1,t2):
    [ta1,tb1,tc1] = t1 ; [ta2,tb2,tc2] = t2
    la = tb1*tc2 - tc1*tb2 ; lb = tc1*ta2 - ta1*tc2 ; lc = ta1*tb2 - tb1*ta2
    u = a*la + b*lb + c*lc
    if u == 0:
        l = [la,lb,lc]
    else:
        l = [la/u,lb/u,lc/u]
    return l

def TriColinear(a,b,c,t1,t2,t3):
    [ta,tb,tc] = t3
    [la,lb,lc] = TriLine(a,b,c,t1,t2)
    return la*ta + lb*tb + lc*tc

# TriParallel returns parallel line thru point ta:tb:tc
def TriParallel(a,b,c,l,t):
    [la,lb,lc] = l
    [ta,tb,tc] = t
    z = - (la*ta + lb*tb + lc*tc)/(a*ta + b*tb + c*tc)  # denominator not zero because def of trilinear coordinates
    m = [la + a*z,lb + b*z,lc + c*z]
    return m

# TriMidpoint returns middle point of two trilinear coordinates points
def TriMidpoint(a,b,c,t1,t2):
    [ta1,tb1,tc1] = t1 ; [ta2,tb2,tc2] = t2
    ta = 2*a*ta1*ta2+b*(ta2*tb1+ta1*tb2)+c*(ta2*tc1+ta1*tc2)
    tb = a*(ta2*tb1+ta1*tb2)+2*b*tb1*tb2+c*(tb2*tc1+tb1*tc2)
    tc = a*(ta2*tc1+ta1*tc2)+b*(tb2*tc1+tb1*tc2)+2*c*tc1*tc2
    u = a*ta + b*tb + c*tc # u = 2 area(ABC) never 0
    t = [ta/u,tb/u,tc/u]
    return t

# TriOrthogonalLine returns orthogonal line of line at a crossing point
# line l = [la,lb,lc] has equation la ta + lb tb + lc tc = 0 where [ta,tb,tc] are trilinear coordinates of point on that line
# and orthogonal line has equation la ma + lb mb + lc mc - (lc mb + lb mc)cs_A - (la mc + lc ma)cs_B - (lb ma + la mb)cs_C
# OLD :
#def TriOrthogonalLine(a,b,c,l,t):
#    [la,lb,lc] = l
#    [ta,tb,tc] = t
#    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
#    ma,mb,mc = var('ma,mb,mc')
#    f1(ma,mb,mc) = la*ma + lb*mb + lc*mc - (lc*mb + lb*mc)*cs_A - (la*mc + lc*ma)*cs_B - (lb*ma + la*mb)*cs_C
#    f2(ma,mb,mc) = ta*ma + tb*mb + tc*mc
#    res = solve([f1(ma,mb,mc) == 0,f2(ma,mb,mc) == 0],ma,mb,mc)
    # Simplifying result : removing free vars
#    res = TriHomogeneousSimplify(res,[a,b,c,ma,mb,mc])
#    m = [res[0][0].right(),res[0][1].right(),res[0][2].right()]
#    return m

def TriOrthogonalLine(a,b,c,l,t):
    [la,lb,lc] = l
    [ta,tb,tc] = t
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    ma = (cs_B*la + cs_A*lb - lc)*tb - (cs_C*la + cs_A*lc - lb)*tc
    mb = (cs_C*lb + cs_B*lc - la)*tc - (cs_B*la + cs_A*lb - lc)*ta
    mc = (cs_C*la + cs_A*lc - lb)*ta - (cs_C*lb + cs_B*lc - la)*tb
    return [ma,mb,mc]


# TriFootOnLine returns orthogonal projection of point t3 on line defined by t1-t2
def TriFootOnLine(a,b,c,t1,t2,t3):
    l = TriLine(a,b,c,t1,t2)
    m = TriOrthogonalLine(a,b,c,l,t3)
    t = TriLineIntersect(a,b,c,l,m)
    return t

# TriFeetOnEdges returns feets on edges reference triangle ABC
def TriFeetOnEdges(a,b,c,t):
    [tri_A,tri_B,tri_C] = TriMainTriangle(a,b,c)
    t1 =  TriFootOnLine(a,b,c,tri_B,tri_C,t)
    t2 =  TriFootOnLine(a,b,c,tri_C,tri_A,t)
    t3 =  TriFootOnLine(a,b,c,tri_A,tri_B,t)
    return [t1,t2,t3]

# TriReflectionLine returns reflection of point t3 with respect to line defined by t1-t2
def TriReflectionLine(a,b,c,t1,t2,t3):
    t4 = TriFootOnLine(a,b,c,t1,t2,t3)
    t5 = TriReflectionPoint(a,b,c,t4,t3)
    return t5

# TriReflectionLinePoint returns reflection of point t with respect to line l
def TriReflectionLinePoint(a,b,c,l,t):
    m = TriOrthogonalLine(a,b,c,l,t)
    t1 = TriLineIntersect(a,b,c,l,m)  # foot of t on line l
    t2 = TriReflectionPoint(a,b,c,t1,t)
    return t2

# TriReflectedLine returns reflection of line t3-t4 with respect to line defined by t1-t2, as two points
def TriReflectedLine(a,b,c,t1,t2,t3,t4):
    t5 = TriReflectionLine(a,b,c,t1,t2,t3)
    t6 = TriReflectionLine(a,b,c,t1,t2,t4)
    return [t5,t6]

# TriReflectionQuadTriangle returns triangle from reflection
#   quad vertices t1, t2, t3, t4
#   other point (eg incenter for a the tangential quad) t5
def TriReflectionQuadTriangle(a,b,c,t1,t2,t3,t4,t5):
    t6 = TriFootOnLine(a,b,c,t1,t2,t5)
    t7 = TriFootOnLine(a,b,c,t2,t3,t5)
    t8 = TriFootOnLine(a,b,c,t3,t4,t5)
    t9 = TriFootOnLine(a,b,c,t4,t1,t5)
    [t10,t11] = TriReflectedLine(a,b,c,t7,t8,t1,t2)
    [t12,t13] = TriReflectedLine(a,b,c,t8,t9,t1,t2)
    [t14,t15] = TriReflectedLine(a,b,c,t9,t6,t1,t2)
    l10_11 = TriLine(a,b,c,t10,t11)
    l12_13 = TriLine(a,b,c,t12,t13)
    l14_15 = TriLine(a,b,c,t14,t15)
    t16 = TriLineIntersect(a,b,c,l10_11,l12_13)
    t17 = TriLineIntersect(a,b,c,l12_13,l14_15)
    t18 = TriLineIntersect(a,b,c,l14_15,l10_11)
    return [t16,t17,t18]


# TriRatioSegment is a function returning trilinear coordinates t5 of P such as UP:UV = v:w  with w odd
# U,V are given by trilinear coordinates t1,t2
# Starting point t3 is chosen such as discriminant of t1,t2,t3 is not 0
#
# When w is even : divide UV segment in two segments UW and WV.
#  - If v lesser than w/2 then dividing UV in v:w ratio is same than dividing UW in v:w1 ratio with w1 = w/2.
#  - If v greater than w/2 then dividing UV in v:w ratio is same than dividing WV in 2v-w1:w1 ratio with w1 = s/2.

def TriRatioSegment(a,b,c,t1,t2,v,w):
    n = Mod(2,w).multiplicative_order() # compute n such as 2^b = 1 mod s
    H = Integer((2^n - 1)*v/w)   # set the "distance" between U to P on UV mapped to [0,1] segment
    Hn = H.bits() ; m = len(Hn)  # binary decomposition (2-adic) of r/s = 0,(bn...b1) ; Hn = [b1,...,bn]
    ta = t1[1]*t2[2] - t1[2]*t2[1] ; tb = t1[2]*t2[0] - t1[0]*t2[2] ; tc = t1[0]*t2[1] - t1[1]*t2[0]
    t4 = t3 = [ta,tb,tc] # Any point not on UV
    for k in xrange(n/2+1):
        for l in xrange(2):
            i = 2*k + l
            if i < n:
                if i < m:
                    if Hn[i] == 0:
                        t5 = TriMidpoint(a,b,c,t1,t4)
                    else:
                        t5 = TriMidpoint(a,b,c,t2,t4)
                else:
                    t5 = TriMidpoint(a,b,c,t1,t4)
                t4 = TriEvalReal(t5)
    # Intersect lines to get P
    l12 = TriEvalReal(TriLine(a,b,c,t1,t2)); l34 = TriEvalReal(TriLine(a,b,c,t3,t4))
    t5 = TriEvalReal(TriLineIntersect(a,b,c,l12,l34))
    return t5

def TriRT_RatioSegment(a,b,c,t1,t2,v,w):
    n = Mod(2,w).multiplicative_order() # compute n such as 2^b = 1 mod s
    H = Integer((2^n - 1)*v/w)   # set the "distance" between U to P on UV mapped to [0,1] segment
    Hn = H.bits() ; m = len(Hn)  # binary decomposition (2-adic) of r/s = 0,(bn...b1) ; Hn = [b1,...,bn]
    ta = t1[1]*t2[2] - t1[2]*t2[1] ; tb = t1[2]*t2[0] - t1[0]*t2[2] ; tc = t1[0]*t2[1] - t1[1]*t2[0]
    t4 = t3 = [ta,tb,tc] # Any point not on UV
    for k in xrange(n/2+1):
        for l in xrange(2):
            i = 2*k + l
            if i < n:
                if i < m:
                    if Hn[i] == 0:
                        t5 = TriMidpoint(a,b,c,t1,t4)
                    else:
                        t5 = TriMidpoint(a,b,c,t2,t4)
                else:
                    t5 = TriMidpoint(a,b,c,t1,t4)
                t4 = TriFactor(t5)
    # Intersect lines to get P
    l12 = TriFactor(TriLine(a,b,c,t1,t2)); l34 = TriFactor(TriLine(a,b,c,t3,t4))
    t5 = TriFactor(TriLineIntersect(a,b,c,l12,l34))
    return t5


# TriCircleInversion returns image of point t2 with respect to circle inversion centered t1 radius k
#
# following code was computed as inverse of the function TriCirclesExternalSimilitudeCenter retrieving the
# external similitude center of two circles radius r and R = r - 1 centered at (ta,tb,tc) and t2
def RT_TriCircleInversion(a,b,c,t1,t2,k2):
    [ta1,tb1,tc1] = t1 ; [ta2,tb2,tc2] = t2
    # Convert to barycentric coordinates
    m1 = a*ta1 + b*tb1 + c*tc1 ; x1 = a*ta1/m1 ; y1 = b*tb1/m1  # barycentric coordinates of inversion center
    m2 = a*ta2 + b*tb2 + c*tc2 ; x2 = a*ta2/m2 ; y2 = b*tb2/m2  # barycentric coordinates of point M
    # Compute barycentric coordinates of image M'
    r1 = b^2*(x1 - x2)^2  + (a^2 + b^2 - c^2)*(x1 - x2)*(y1 - y2) + a^2*(y1 - y2)^2
    r2 = 1 - k2/r1
    r = 1/r2
    x = ((r - 1)*x2 + x1)/r
    y = ((r - 1)*y2 + y1)/r
    z = 1 - x - y
    # And compute trilinear coordinates from barycentric coordinates
    ta = x/a ; tb = y/b ; tc = z/c
    return [ta,tb,tc]

def TriCircleInversion(a,b,c,t1,t2,k):
    return RT_TriCircleInversion(a,b,c,t1,t2,k^2)

# TriAtQuadrance returns point at a given quadrance q from t1 on a line t1t2
# we use RT_TriCircleInversion for this and so we don't need a circle-line intersection
def TriAtQuadrance(a,b,c,t1,t2,q):
    q12 = TriQuadranceTwoPoints(a,b,c,t1,t2)
    k2 = sqrt(q12*q) # t1t2.t1t3 = k^2 => q12.q = k^4
    t3 = RT_TriCircleInversion(a,b,c,t1,t2,k2)
    return t3

# TriScaledPoint returns point at a scaled distance from t1 on t1t2
def TriScaledPoint(a,b,c,t1,t2,scale):
    q12 = TriQuadranceTwoPoints(a,b,c,t1,t2)
    t = TriAtQuadrance(a,b,c,t1,t2,q12*scale^2)
    return t

# TriTriangleOrthocenter returns orthocenter of triangle t1-t2-t3 as intersection of 2 altitude lines
def TriTriangleOrthocenter(a,b,c,t1,t2,t3):
    t4 = TriFootOnLine(a,b,c,t1,t2,t3)
    t5 = TriFootOnLine(a,b,c,t2,t3,t1)
    l34 = TriLine(a,b,c,t3,t4)
    l15 = TriLine(a,b,c,t1,t5)
    t = TriLineIntersect(a,b,c,l34,l15)
    return t

# TriTriangleNinePointCenter returns nine point center, midpoint of circumcenter and orthocenter
def TriTriangleNinePointCenter(a,b,c,t1,t2,t3):
    t4 =  TriTriangleCircumcenter(a,b,c,t1,t2,t3)
    t5 =  TriTriangleOrthocenter(a,b,c,t1,t2,t3)
    t = TriMidpoint(a,b,c,t4,t5)
    return t


# TriMedialTriangle returns midpoints of edges : medial triangle
def TriMedialTriangle(a,b,c):
    t_BC = [0,c,b]
    t_CA = [c,0,a]
    t_AB = [b,a,0]
    return [t_BC,t_CA,t_AB]

# TriTriangleCentroid returns centroid of triangle t1-t2-t3 as intersection of 2 median lines
def TriTriangleCentroid(a,b,c,t1,t2,t3):
    [t4,t5,t6] = TriMedialTriangle(a,b,c)
    l34 = TriLine(a,b,c,t3,t4)
    l15 = TriLine(a,b,c,t1,t5)
    t = TriLineIntersect(a,b,c,l34,l15)
    return t

# TriAngleBisectorPoints returns intersection angles bisector with opposite edges
def TriAngleBisectorPoints(a,b,c):
    [tri_A,tri_B,tri_C] = TriMainTriangle(a,b,c)
    tri_I = TriIncenter(a,b,c)
    lAB = TriLine(a,b,c,tri_A,tri_B) ; lBC = TriLine(a,b,c,tri_B,tri_C) ; lAC = TriLine(a,b,c,tri_A,tri_C)
    lAI = TriLine(a,b,c,tri_A,tri_I) ; lBI = TriLine(a,b,c,tri_B,tri_I) ; lCI = TriLine(a,b,c,tri_C,tri_I)
    t1 = TriLineIntersect(a,b,c,lAI,lBC) ; t2 = TriLineIntersect(a,b,c,lBI,lAC) ; t3 = TriLineIntersect(a,b,c,lCI,lAB)
    return [t1,t2,t3]

# TriTriangle120CirclesCenters return center of circles thru bisector points (see function TriAngleBisectorPoints) and two vertices
def TriTriangle120CirclesCenters(a,b,c):
    t1 = [(2*a + b + c), (a + c), (a + b)] # when angle at A is 120 degrees
    t2 = [(b + c), (a + 2*b + c), (a + b)] # when angle at B is 120 degrees
    t3 = [(b + c), (a + c), (a + b + 2*c)] # when angle at C is 120 degrees
    return [t1,t2,t3]

# TriAltitudePoints returns altitude points on line l2 at quadrance q to line l1
def TriAltitudePoints(a,b,c,S,q,l1,l2):
    [la1,lb1,lc1] = l1
    [la2,lb2,lc2] = l2
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    x,y,z = var('x,y,z')
    f1(x,y,z) = la2*x + lb2*y + lc2*z
    f2(x,y,z) = 4*S^2*(la1*x + lb1*y + lc1*z)^2 - q*(la1^2 + lb1^2 + lc1^2 - 2*lb1*lc1*cs_A - 2*lc1*la1*cs_B - 2*la1*lb1*cs_C)*(a*x + b*y + c*z)^2
    res = solve([f1(x,y,z) == 0,f2(x,y,z) == 0],x,y,z)
    # Simplifying result : removing free vars
    res = TriHomogeneousSimplify(res,[a,b,c,S])
    t1 = [res[0][0].right(),res[0][1].right(),res[0][2].right()]
    t2 = [res[1][0].right(),res[1][1].right(),res[1][2].right()]
    return [t1,t2]


# TriTriangleCircumcenter returns circumcenter of triangle t1-t2-t3 as intersection of 2 segment bissectors
def TriTriangleCircumcenter(a,b,c,t1,t2,t3):
    l12 = TriSegmentBissectorLine(a,b,c,t1,t2)
    l23 = TriSegmentBissectorLine(a,b,c,t2,t3)
    t = TriLineIntersect(a,b,c,l12,l23)
    return t

# TriTriangleApexVertices returns apex vertices t4,t5 such as t1-t2-t4 (or t5) is isosceles and t4 (or t5) at altitude
# d = sqrt(q) from t1-t2 line
def TriTriangleApexVertices(a,b,c,S,q,t1,t2):
    t3 = TriMidpoint(a,b,c,t1,t2)
    l12 = TriLine(a,b,c,t1,t2)
    l3 = TriOrthogonalLine(a,b,c,l12,t3)
    [t4,t5] = TriAltitudePoints(a,b,c,S,q,l12,l3)
    return [t4,t5]

# TriTriangleApexVertex returns apex vertex t such as t1-t2-t  is isosceles and t at altitude
# d = sqrt(q) from t1-t2 line
# Use TriReflectionPoint to reflect on midpoint of t1t2 to get symmetric point
def TriTriangleApexVertex(a,b,c,q,t1,t2):
    t3 = TriMidpoint(a,b,c,t1,t2)
    l12 = TriLine(a,b,c,t1,t2)
    l3 = TriOrthogonalLine(a,b,c,l12,t3)
    # Choose BC if not parallel with l3
    if c*l3[1] - b*l3[2] == 0:
        l = [0,1,0]
    else:
        l = [1,0,0]
    t4 = TriLineIntersect(a,b,c,l3,l)
    k = sqrt(sqrt(TriQuadranceTwoPoints(a,b,c,t3,t4)*q))
    t = TriCircleInversion(a,b,c,t3,t4,k)
    return t

# TriTriangleFermatPoint returns Fermat point of triangle t1-t2-t3 as intersection of Simson lines
def TriTriangleFermatPoint(a,b,c,t1,t2,t3):
    t4 = TriTriangleCircumcenter(a,b,c,t1,t2,t3)
    q12 = TriQuadranceTwoPoints(a,b,c,t1,t2)
    q31 = TriQuadranceTwoPoints(a,b,c,t3,t1)
    # Compute Simson line t3-t6
    t12 = TriMidpoint(a,b,c,t1,t2)
    k = sqrt(sqrt(TriQuadranceTwoPoints(a,b,c,t12,t4)*q12*3/4))
    t5 = TriCircleInversion(a,b,c,t12,t4,k)
    t6 = TriReflectionPoint(a,b,c,t12,t5)
    # Compute Simson line t2-t8
    t31 = TriMidpoint(a,b,c,t3,t1)
    k = sqrt(sqrt(TriQuadranceTwoPoints(a,b,c,t31,t4)*q31*3/4))
    t7 = TriCircleInversion(a,b,c,t31,t4,k)
    t8 = TriReflectionPoint(a,b,c,t31,t7)
    # Intersection Simson lines to get Fermat point
    l1 = TriLine(a,b,c,t3,t6)
    l2 = TriLine(a,b,c,t2,t8)
    t = TriLineIntersect(a,b,c,l1,l2)
    return t

# TriTrilliumTriangle returns trilinear coordinates of apex points of isosceles triangles on edges (see Trillium theorem)
# They are centers Oa,Ob,Oc of circles through  two vertices and the incenter
def TriTrilliumTriangle(a,b,c):
    t_BC =  [-a, (b + c), (b + c)]
    t_CA =  [(a + c), -b, (a + c)]
    t_AB =  [(a + b), (a + b), -c]
    return [t_BC,t_CA,t_AB]

# TriTrilliumQuadrances returns quadrances between apex points of isoscles triangles on edges to vertices A,B,C
def TriTrilliumQuadrances(a,b,c):
    qA = a^2*b*c/((a + b + c)*(-a + b + c))
    qB = a*b^2*c/((a + b + c)*(a - b + c))
    qC = a*b*c^2/((a + b + c)*(a + b - c))
    return [qA,qB,qC]

# TriTrilliumBisectors returns perpendicular bisectors lines of edges ABC triangle : lines from midpoints of edges to apex isosceles on edges
def TriTrilliumBisectors(a,b,c):
    lA = [(b + c)*(b - c), a*b, -a*c]
    lB = [-b*a,(c + a)*(c - a) , b*c]
    lC = [c*a, -c*b, (a + b)*(a - b)]
    return [lA,lB,lC]

# TriTrilliumAltitudes returns distances between apex points and midpoints of edges
# if triangle is not degenerated to flat segment, denominators are positive different from zero
def TriTrilliumAltitudes(a,b,c,S):
    hA = (2*a)/((b + c)^2 - a^2)*S
    hB = (2*b)/((c + a)^2 - b^2)*S
    hC = (2*c)/((a + b)^2 - c^2)*S
    return [hA,hB,hC]

# TriTrilliumCotangents returns cotangent of swing angles on edges
def TriTrilliumCotangents(a,b,c,S):
    [hA,hB,hC] = TriTrilliumAltitudes(a,b,c,S)
    return [a/(2*hA),b/(2*hB),c/(2*hC)]

# TriTrilliumSines returns sines of angles BOA2, COB2, AOC2 where O is circumcenter and A2,B2,C2 apex points
def TriTrilliumSines(a,b,c,S):
    return[2*S/(b*c),2*S/(c*a),2*S/(a*b)]

# TriTrilliumCosines returns cosines of angles BOA2, COB2, AOC2 where O is circumcenter and A2,B2,C2 apex points
def TriTrilliumCosines(a,b,c):
    return[(-a^2 + b^2 + c^2)/(2*b*c),(a^2 - b^2 + c^2)/(2*c*a),(a^2 + b^2 - c^2)/(2*a*b)]

# TriTrilliumCirclesFunctions returns functions of circles centered at Oa,Ob,Oc through two vertices and the incenter
def TriTrilliumCirclesFunctions(a,b,c):
    fOA =  [-1, 0, 0]
    fOB =  [0, -1, 0]
    fOC =  [0, 0, -1]
    return [fOA,fOB,fOC]

# TriEdgesPerpendicularPoints returns trilinear coordinates of colinear points P,Q,R (on BC,CA,AB extended edges) given a point O trilinear t, such as
# OA perpendicular with OP, OB perpendicular with OQ,   OC perpendicular with OR
# It's a Felix Laroche's colineairty theore;
def TriEdgesPerpendicularPoints(a,b,c,t):
    [ta,tb,tc] = t ; x = tb/ta ; y = tc/ta
    [S_A,S_B,S_C] = TriConwayParameters(a,b,c) # S_A = (-a^2 + b^2 + c^2)/2 ; S_B = (a^2 - b^2 + c^2)/2 ; S_C = (a^2 + b^2 - c^2)/2
    u = b*S_B*x^2 - c*S_C*x*y - a*S_A*x - a*b*c*y
    v = b*S_B*x*y - c*S_C*y^2 + a*b*c*x + a*S_A*y
    w = a*b*c*x*y + b*S_B*x + c*S_C*y - a*S_A
    t_P = [0,u,v] ; t_Q = [w,0,-v] ; t_R = [w,u,0]
    return [t_P,t_Q,t_R]

# TriCircumcevianTriangle returns circumcevian triangle of a point P with trilinear t : a triangle VA,VB,VP
# such as intersection of VA-P line with circumcircle of ABC are (A,VA) points (and the same for VB,VC)
def TriCircumcevianTriangle(a,b,c,t):
    [ta,tb,tc] = t
    t_VA = [-a*tb*tc,(b*tc+c*tb)*tb,(b*tc+c*tb)*tc]
    t_VB = [(c*ta+a*tc)*ta,-b*tc*ta,(c*ta+a*tc)*tc]
    t_VC = [(a*tb+b*ta)*ta,(a*tb+b*ta)*tb,-c*ta*tb]
    return [t_VA,t_VB,t_VC]

# TriMainTriangle returns trilinear coordinates of vertices
def TriMainTriangle(a,b,c):
    t_A = [1,0,0] ; t_B = [0,1,0] ; t_C = [0,0,1]
    return [t_A,t_B,t_C]

# TriOrthicTriangle returns trilinear coordinates of vertices of the orthic triangle (HA,HB,HC)
# Vertices of orthic triangle are orthogonal projection of vertices to opposite edges
def TriOrthicTriangle(a,b,c):
    [sec_A,sec_B,sec_C] = TriSecants(a,b,c)
    t_HA = [0,sec_B,sec_C]
    t_HB = [sec_A,0,sec_C]
    t_HC = [sec_A,sec_B,0]
    return [t_HA,t_HB,t_HC]

# TriOrthocenter returns orthocenter H = X(4) intersection of orthic lines A_HA, B_HB, C_HC
def TriOrthoCenter(a,b,c):
    [S_A,S_B,S_C] = TriConwayParameters(a,b,c)
    [f_A,f_B,f_C] = [a*S_A,b*S_B,c*S_C]
    t_H = [f_B*f_C,f_A*f_C,f_A*f_B]  # t_H = TriSecants(a,b,c)
    return t_H

# TriOrthocenterQuadrances returns squared distances from orthocenter H to vertices A,B,C
def TriOrthocenterQuadrances(a,b,c,S):
    qHA = (-a^2 + b^2 + c^2)^2*a^2/(16*S^2)
    qHB = ( a^2 - b^2 + c^2)^2*b^2/(16*S^2)
    qHC = ( a^2 + b^2 - c^2)^2*c^2/(16*S^2)
    return [qHA,qHB,qHC]

# TriCircumcenter returns circumcenter O
def TriCircumCenter(a,b,c):
    [S_A,S_B,S_C] = TriConwayParameters(a,b,c)
    [f_A,f_B,f_C] = [a*S_A,b*S_B,c*S_C]  # = [a*(-a^2+b^2+c^2),b*(a^2-b^2+c^2),c*(a^2+b^2-c^2)] = [a ctA, b ctB, c ctC]
    t_O = [f_A,f_B,f_C]  # t_O = TriCosines(a,b,c)
    return t_O

# TriIncenter returns incenter I
def TriIncenter(a,b,c):
    t_I = [1,1,1]
    return t_I

# TriExcenters returns excenters xA,xB,xC
# excenters are centers of circles tangent to two edges lines internally and one edge line externally
def TriExcenters(a,b,c):
    t_xA = [-1,1,1]
    t_xB = [1,-1,1]
    t_xC = [1,1,-1]
    return [t_xA,t_xB,t_xC]

# TriTangentialTriangle returns vertices of the tangential triangle (intersection of tangents to circumcircle at A,B,C vertices)
def TriTangentialTriangle(a,b,c):
    t_TA = [-a,b,c]
    t_TB = [a,-b,c]
    t_TC = [a,b,-c]
    return [t_TA,t_TB,t_TC]

# TriTangentialTriangleQuadrances returns quadrances (square of edges lengths) for tangential triangle, given quadrances for ABC
def TriTangentialTriangleQuadrances(qa,qb,qc):
    qta = (4*qa^3*qb*qc)/(qa^2 - (qb - qc)^2)^2
    qtb = (4*qb^3*qc*qa)/(qb^2 - (qc - qa)^2)^2
    qtc = (4*qc^3*qa*qb)/(qc^2 - (qa - qb)^2)^2
    return [qta,qtb,qtc]

# TriTangentialTriangleIncenter returns incenter of tangential triangle of an A-obtuse triangle ABC
def TriTangentialTriangleIncenter(a,b,c):
    t_J = [-a*(a^2 + b^2 + c^2), b*(a^2 + b^2 - c^2), c*(a^2 - b^2 + c^2)]
    return t_J

# TriTangentialTriangleFootEulerLineIncenter returns foot on Euler line for incenter of tangential triangle of an A-obtuse triangle ABC
def TriTangentialTriangleFootEulerLineIncenter(a,b,c):
    t_H =  [-a*(a^8 - 2*a^4*b^4 + b^8 + a^4*b^2*c^2 + a^2*b^4*c^2 - 2*b^6*c^2 - 2*a^4*c^4 + a^2*b^2*c^4 + 2*b^4*c^4 - 2*b^2*c^6 + c^8),
          b*(a^8 - 2*a^4*b^4 + b^8 - 2*a^6*c^2 + a^4*b^2*c^2 + 3*a^2*b^4*c^2 - 2*b^6*c^2 + 2*a^4*c^4 - 3*a^2*b^2*c^4 + 2*b^2*c^6 - c^8),
          c*(a^8 - 2*a^6*b^2 + 2*a^4*b^4 - b^8 + a^4*b^2*c^2 - 3*a^2*b^4*c^2 + 2*b^6*c^2 - 2*a^4*c^4 + 3*a^2*b^2*c^4 - 2*b^2*c^6 + c^8)]
    return t_H

# TriTangentialTriangleQuadranceEulerLineIncenter returns JH^2 quadrance to Euler line from incenter of tangential triangle of an A-obtuse triangle ABC
def TriTangentialTriangleQuadranceEulerLineIncenter(qa,qb,qc):
    return (4*qa^2*qb^2*qc^2*(qb - qc)^2)/((qa - qb + qc)^2*(qa + qb - qc)^2*(qa^3 + qb^3 + qc^3 + 3*qa*qb*qc - qa^2*(qb + qc) - qb^2*(qc + qa) - qc^2*(qa + qb)))

# TriTangentialTriangleQuadranceEulerLineExcenter returns OH^2 quadrance to Euler line from incenter of tangential triangle (O its Q-excenter) of an A-obtuse triangle ABC
def TriTangentialTriangleQuadranceEulerLineExcenter(qa,qb,qc):
    return -4*(2*qa^2 - qa*qb - qb^2 - qa*qc + 2*qb*qc - qc^2)^2*qa^2*qb^2*qc^2/((qa^3 - qa^2*qb - qa*qb^2 + qb^3 - qa^2*qc + 3*qa*qb*qc - qb^2*qc - qa*qc^2 - qb*qc^2 + qc^3)*(qa + qb - qc)^2*(qa - qb + qc)^2*(qa^2 - 2*qa*qb + qb^2 - 2*qa*qc - 2*qb*qc + qc^2))

# TriTangentialTriangleQuadranceExcenterCirumradius returns circumradius of circle thru cirumcenter O and TB,TC vertices of the tangential triangle
def TriTangentialTriangleQuadranceExcenterCircumradius(qa,qb,qc,S):
    return 1/2*sqrt(qa^3)*qb*qc/(S*(qa + qb - qc)*(qa - qb + qc))

# TriTangentialTriangleTACircleCenter returns center OA of circle thru vertices TB,TC of tangential triangles and O A-excenter of TATBTC
def TriTangentialTriangleTACircleCenter(a,b,c):
    t_OA = [-a*(a^6 - a^4*b^2 - a^2*b^4 + b^6 - a^4*c^2 - 2*a^2*b^2*c^2 - b^4*c^2 - a^2*c^4 - b^2*c^4 + c^6),
            b*(a^4 - 2*a^2*b^2 + b^4 - 2*b^2*c^2 + c^4)*(a^2 + b^2 - c^2),
            c*(a^4 + b^4 - 2*a^2*c^2 - 2*b^2*c^2 + c^4)*(a^2 - b^2 + c^2)]
    return t-OA

# TriContactTriangle, returns contact triangle, intersectin of incircle with edges
def TriContactTriangle(a,b,c):
    t_TA = [0,a*c/(a-b+c),a*b/(a+b-c)]    # intersect with BC line
    t_TB = [b*c/(-a+b+c),0,a*b/(a+b-c)]   # intersect with AC ine
    t_TC = [b*c/(-a+b+c),a*c/(a-b+c),0]   # intersect with AB line
    return [t_TA,t_TB,t_TC]

# TriIntouchTriangle, returns intouch triangle = contact triangle, but here trilinear
# coordinares are given using semiperimeter
def TriIntouchTriangle(a,b,c):
    s = RT_TriangleSemiperimeter(a,b,c)
    t_IA = [0,(s-c)/b,(s-b)/c]   # intersect with BC line
    t_IB = [(s-c)/a,0,(s-a)/c]   # intersect with AC ine
    t_IC = [(s-b)/a,(s-a)/b,0]   # intersect with AB line
    return [t_IA,t_IB,t_IC]

# TriAExtouchtriangle, TriBExtouchtriangle, TriBExtouchtriangle returns extouch triangle for A,B or C vertice
def TriAExtouchTriangle(a,b,c):
    s = RT_TriangleSemiperimeter(a,b,c)
    t_XA = [0,(s-b)/b,(s-c)/c]   # intersect with BC line
    t_XB = [-(s-b)/a, 0,s/c]     # intersect with AC line
    t_XC = [-(s-c)/a,s/b,0]      # intersect with AB line
    return [t_XA,t_XB,t_XC]

def TriBExtouchTriangle(a,b,c):
    s = RT_TriangleSemiperimeter(a,b,c)
    t_XA = [0,-(s-a)/b,s/c]      # intersect with BC line
    t_XB = [(s-a)/a,0,(s-c)/c]   # intersect with AC line
    t_XC = [s/a,-(s-c)/b,0]      # intersect with AB line
    return [t_XA,t_XB,t_XC]

def TriCExtouchTriangle(a,b,c):
    s = RT_TriangleSemiperimeter(a,b,c)
    t_XA = [0,s/b,-(s-a)/c]      # intersect with BC line
    t_XB = [s/a,0,-(s-b)/c]      # intersect with AC line
    t_XC = [(s-a)/a,(s-b)/b,0]   # intersect with AB line
    return [t_XA,t_XB,t_XC]

# TriABCtouchTriangle is triangle with same area than intouch triangle, which one vertice for each extouch triangle
def TriABCtouchTriangle(a,b,c):
    s = RT_TriangleSemiperimeter(a,b,c)
    t_XA = [0,(s-b)/b,(s-c)/c]   # intersect with BC line
    t_XB = [(s-a)/a,0,(s-c)/c]   # intersect with AC line
    t_XC = [(s-a)/a,(s-b)/b,0]   # intersect with AB line
    return [t_XA,t_XB,t_XC]

# TriHaimovTriangle return Haimov triangle QaQbQc of a point P
# If PaPbPc is Ceva triangle of P
#  Qa is intersection (different of A) of circle throuh A,B,Pb and circle throuh A,C,Pc
#  Qb is intersection (different of B) of circle throuh B,C,Pc and circle throuh B,A,Pa
#  Qc is intersection (different of C) of circle throuh C,A,Pc and circle throuh C,B,Pb
def TriHaimovTriangle(a,b,c,t):
    [ta,tb,tc] = t ; k = a*ta + b*tb + c*tc
    t_QA = [-a^3*ta^2 + a*b^2*ta^2 + a*c^2*ta^2 - a^2*b*ta*tb + b^3*ta*tb - a^2*c*ta*tc + c^3*ta*tc - a*b*c*tb*tc, k*(a*ta + b*tb)*b, k*(a*ta + c*tc)*c]
    t_QB = [k*(a*ta + b*tb)*a, a^3*ta*tb - a*b^2*ta*tb + a^2*b*tb^2 - b^3*tb^2 + b*c^2*tb^2 - a*b*c*ta*tc - b^2*c*tb*tc + c^3*tb*tc, k*(b*tb + c*tc)*c]
    t_QC = [k*(a*ta + c*tc)*a, k*(b*tb + c*tc)*b, -a*b*c*ta*tb + a^3*ta*tc - a*c^2*ta*tc + b^3*tb*tc - b*c^2*tb*tc + a^2*c*tc^2 + b^2*c*tc^2 - c^3*tc^2]
    return [t_QA,t_QB,t_QC]

# TriCentroidHaimovTriangle returns Haimov triangle for centroid
def TriCentroidHaimovTriangle(a,b,c):
    t_QA = [-2*a^2 + b^2 + c^2, 3*a*b, 3*a*c]
    t_QB = [ 3*a*b, a^2 - 2*b^2 + c^2, 3*b*c]
    t_QC = [ 3*a*c, 3*b*c, a^2 + b^2 - 2*c^2]
    return [t_QA,t_QB,t_QC]

# TriIncenterHaimovTriangle returns Haimov triangle for incenter
def TriIncenterHaimovTriangle(a,b,c):
    t_QA =  [-a^2 + b^2 - b*c + c^2, (a + b)*b, (a + c)*c]
    t_QB =  [(a + b)*a, a^2 - b^2 - a*c + c^2, (b + c)*c]
    t_QC =  [(a + c)*a, (b + c)*b, a^2 - a*b + b^2 - c^2]
    return [t_QA,t_QB,t_QC]

# TriCircumcenterHaimovTriangle returns Haimov triangle for incenter
def TriCircumcenterHaimovTriangle(a,b,c):
    t_QA =  [(a^2 - b^2 + b*c - c^2)*(- a^2 + b^2 + b*c + c^2)*a, (a^4 - 2*a^2*b^2 + b^4 - a^2*c^2 - b^2*c^2)*b, (a^4 - a^2*b^2 - 2*a^2*c^2 - b^2*c^2 + c^4)*c]
    t_QB =  [(-a^4 + 2*a^2*b^2 - b^4 + a^2*c^2 + b^2*c^2)*a, (a^2 - b^2 + a*c + c^2)*(a^2 - b^2 - a*c + c^2)*b, (a^2*b^2 - b^4 + a^2*c^2 + 2*b^2*c^2 - c^4)*c]
    t_QC =  [(-a^4 + a^2*b^2 + 2*a^2*c^2 + b^2*c^2 - c^4)*a, (a^2*b^2 - b^4 + a^2*c^2 + 2*b^2*c^2 - c^4)*b, (a^2 + a*b + b^2 - c^2)*(a^2 - a*b + b^2 - c^2)*c]
    return [t_QA,t_QB,t_QC]

# TriNagelPointHaimovTriangle returns Haimov triangle for Nagel point
def TriNagelPointHaimovTriangle(a,b,c):
    t_QA =  [-2*a + b + c, a, a]
    t_QB =  [ b, a - 2*b + c, b]
    t_QC =  [ c, c, a + b - 2*c]
    return [t_QA,t_QB,t_QC]

# TriOrthoCenterHaimovTriangle returns Haimov triangle for Orthocenter
def TriOrthocenterHaimovTriangle(a,b,c):
    t_QA =  [0, (a^2 + b^2 - c^2)*c, (a^2 - b^2 + c^2)*b]
    t_QB =  [(a^2 + b^2 - c^2)*c, 0, (-a^2 + b^2 + c^2)*a]
    t_QC =  [(a^2 - b^2 + c^2)*b, (-a^2 + b^2 + c^2)*a, 0]
    return [t_QA,t_QB,t_QC]

# TriNinePointCenter returns Nine Point Center M, center of Nine Points Circle
def TriNinePointCenter(a,b,c):
    t_M = [  b*c*(a^2*(b^2+c^2)-(b^2-c^2)^2) , c*a*(b^2*(c^2+a^2)-(c^2-a^2)^2) , a*b*(c^2*(a^2+b^2)-(a^2-b^2)^2) ]
    return t_M

# TriCentroid returns centroid G
def TriCentroid(a,b,c):
    t_G = [1/a,1/b,1/c]
    return t_G

# TriSpiekerCenter returns Spieker Center Sp
def TriSpiekerCenter(a,b,c):
    t_Sp = [b*c*(b+c),a*c*(a+c),a*b*(a+b)]
    return t_Sp

# TriGergonnePoint returns Gergonne point Ge
def TriGergonnePoint(a,b,c):
    t_Ge = [1/(a*(-a+b+c)),1/(b*(a-b+c)),1/(c*(a+b-c))]
    return t_Ge

# TriMittelpunkt returns Mittelpunkt Mp
def TriMittelpunkt(a,b,c):
    t_Mp = [(-a+b+c),(a-b+c),(a+b-c)]
    return t_Mp

# TriSymmedianPoint returns Symmedian point (= Lemoine point X6) K
def TriSymmedianPoint(a,b,c):
    t_K = [a,b,c]
    return t_K

# TriNagelPoint X8
def TriNagelPoint(a,b,c):
    t_N = [(b+c-a)/a,(c+a-b)/b,(a+b-c)/c]
    return t_N

# TriNagelLine returns the Nagel line
def TriNagelLine(a,b,c):
    l = [a*(b-c),b*(c-a),c*(a-b)]
    return l

# TriNagelTriangle returns the intersection of Nagel line with edges lines
def TriNagelTriangle(a,b,c):
    t_NA = [0,c*(a-b), b*(a-c)]  # intersect with BC line
    t_NB = [c*(a-b),0,a*(c-b)]   # intersect with AC ine
    t_NC = [b*(a-c),a*(b-c),0]   # intersect with AB line
    return [t_NA,t_NB,t_NC]

# TriToricelliPoint returns 1st Fermat point, isogonic center or Toricelli point X13
# It's point minimizing the sum of absolute distance to edges of triangle
# I don't use Kimberling definition, but mine as the intersection point of Simpson lines
def TriToricelliPoint(a,b,c):
    sn30 = 1/2 ; cs30 = sqrt(3)/2  # sin(pi/6) and cos(pi/6)
    [sn_A,sn_B,sn_C] = TriTriangleSines(a,b,c)
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    return[(a*cs30*cs_B + b*sn30*sn_A)*(b*cs30*cs_C + c*sn30*sn_B), (b*cs30*cs_A + a*sn30*sn_B)*(a*cs30*cs_C + c*sn30*sn_A), (b*cs30*cs_A + a*sn30*sn_B)*(a*cs30*cs_B + b*sn30*sn_A)]

def TriFirstFermatPoint(a,b,c):
    return TriToricelliPoint(a,b,c)

def RT_TriToricelliPoint(a,b,c,S):
    ks = 4/sqrt(3)*S
    kabc = a^2 + b^2 - c^2 + ks
    kbca = b^2 + c^2 - a^2 + ks
    kcab = c^2 + a^2 - b^2 + ks
    return [kabc*kcab/a,kabc*kbca/b,kcab*kbca/c]

# TriOuterNapoleonTriangle returns the Outer Napoleon triangle where vertices are
# centers of equilateral triangles constucted on edges
def TriOuterNapoleonTriangle(a,b,c):
    sn30 = 1/2 ; cs30 = sqrt(3)/2  # sin(pi/6) and cos(pi/6)
    [sn_A,sn_B,sn_C] = TriTriangleSines(a,b,c)
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    t_ONA = [           -1/2,         sn_C*cs30 + cs_C*sn30, sn_B*cs30 + cs_B*sn30 ]
    t_ONB = [sn_C*cs30 + cs_C*sn30,          -1/2,           sn_A*cs30 + cs_A*sn30 ]
    t_ONC = [sn_B*cs30 + cs_B*sn30,   sn_A*cs30 + cs_A*sn30,         -1/2          ]
    return [t_ONA,t_ONB,t_ONC]

def TriOuterNapoleonTriangleArea(a,b,c):
    return (sqrt(3)/24)*(a^2+b^2+c^2) + (1/2)*TriArea(a,b,c)

# TriFirstNapoleonPoint returns X17 center
# center function is csc(A+pi/6) = 1/sin(A+pi/6)
def TriFirstNapoleonPoint(a,b,c):
    sn30 = 1/2 ; cs30 = sqrt(3)/2  # sin(pi/6) and cos(pi/6)
    [sn_A,sn_B,sn_C] = TriTriangleSines(a,b,c)
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    return [1/(sn_A*cs30 + cs_A*sn30),1/(sn_B*cs30 + cs_B*sn30),1/(sn_C*cs30 + cs_C*sn30)]

# TriFirstNapoleonTriangle returns vertices making equilateral triangles with edges
# they are point on Simpson lines too, and intersecion of (A,EQA) (B,EQB) (C,EQC) is Toricelli (1st Fermat) point
def RT_TriFirstNapoleonTriangle(a,b,c,sn_A,cs_A,sn_B,cs_B,sn_C,cs_C,sn30,cs30):
    t_EQA = [       -a*cs30,             c*sn_A*sn30 + a*cs_C*cs30,   b*sn_A*sn30 + a*cs_B*cs30 ]
    t_EQB = [c*sn_B*sn30 + b*cs_C*cs30,          -b*cs30,             a*sn_B*sn30 + b*cs_A*cs30 ]
    t_EQC = [b*sn_C*sn30 + c*cs_B*cs30,  a*sn_C*sn30 + c*cs_A*cs30,          -c*cs30            ]
    return [t_EQA,t_EQB,t_EQC]

def TriFirstNapoleonTriangle(a,b,c):
    sn30 = 1/2 ; cs30 = sqrt(3)/2  # sin(pi/6) and cos(pi/6)
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    [sn_A,sn_B,sn_C] = TriTriangleSines(a,b,c)
    return  RT_TriFirstNapoleonTriangle(a,b,c,sn_A,cs_A,sn_B,cs_B,sn_C,cs_C,sn30,cs30)

def TriLemoineProblemSolution(a,b,c):
    [tri_A,tri_B,tri_C] = TriMainTriangle(a,b,c)
    [t_EQA,t_EQB,t_EQC] = TriFirstNapoleonTriangle(a,b,c)
    t_LA =  TriMidpoint(a,b,c,tri_A,t_EQA)
    t_LB =  TriMidpoint(a,b,c,tri_B,t_EQB)
    t_LC =  TriMidpoint(a,b,c,tri_C,t_EQC)
    return [t_LA,t_LB,t_LC]

# TriQuadranceSimpsonSegments returns quadrance of segments on Simpson lines from triangle vertices ABC to equilateral
# triangles vertices
def TriQuadranceSimpsonSegments(a,b,c):
    return (1/2)*(a^2+b^2+c^2) + 2*sqrt(3)*TriArea(a,b,c)

# TriInnerNapoleonTriangle returns the Outer Napoleon triangle where vertices are
# centers of equilateral triangles constucted on edges
def RT_TriInnerNapoleonTriangle(a,b,c,sn_A,cs_A,sn_B,cs_B,sn_C,cs_C,sn30,cs30):
    t_INA =  [b*c, (2*c*cs_B*sn30 + 2*c*cs30*sn_B - a)*c, (2*b*cs_C*sn30 + 2*b*cs30*sn_C - a)*b]
    t_INB =  [(2*c*cs_A*sn30 + 2*c*cs30*sn_A - b)*c, a*c, (2*a*cs_C*sn30 + 2*a*cs30*sn_C - b)*a]
    t_INC =  [(2*b*cs_A*sn30 + 2*b*cs30*sn_A - c)*b, (2*a*cs_B*sn30 + 2*a*cs30*sn_B - c)*a, a*b]
    return [t_INA,t_INB,t_INC]

def TriInnerNapoleonTriangle(a,b,c):
    sn30 = 1/2 ; cs30 = sqrt(3)/2  # sin(pi/6) and cos(pi/6)
    [sn_A,sn_B,sn_C] = TriTriangleSines(a,b,c)
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    return RT_TriInnerNapoleonTriangle(a,b,c,sn_A,cs_A,sn_B,cs_B,sn_C,cs_C,sn30,cs30)

def InnerNapoleonTriangleArea(a,b,c):
    return (sqrt(3)/24)*(a^2+b^2+c^2) - (1/2)*TriArea(a,b,c)

# TriSecondNapoleonPoint returns X18  point
# center function is csc(A-pi/6) = 1/sin(A-pi/6)
def TriSecondNapoleonPoint(a,b,c):
    sn30 = 1/2 ; cs30 = sqrt(3)/2  # sin(pi/6) and cos(pi/6)
    [sn_A,sn_B,sn_C] = TriTriangleSines(a,b,c)
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    return [1/(sn_A*cs30 - cs_A*sn30),1/(sn_B*cs30 - cs_B*sn30),1/(sn_C*cs30 - cs_C*sn30)]

# TriSecondNapoleonTriangle returns vertices making equilateral triangles with edges
# they are point on Simpson lines too, and intersecion of (A,FQA) (B,FQB) (C,FQC) is 2nd Fermat point
def RT_TriSecondNapoleonTriangle(a,b,c,sn_A,cs_A,sn_B,cs_B,sn_C,cs_C,sn30,cs30):
    t_FQA = [a*b*c*cs30, (a*c*cs30*cs_B + b*c*sn30*sn_A - a^2*cs30)*c, (a*b*cs30*cs_C + b*c*sn30*sn_A - a^2*cs30)*b]
    t_FQB = [(b*c*cs30*cs_A + a*c*sn30*sn_B - b^2*cs30)*c, a*b*c*cs30, (a*b*cs30*cs_C + a*c*sn30*sn_B - b^2*cs30)*a]
    t_FQC = [(b*c*cs30*cs_A + a*b*sn30*sn_C - c^2*cs30)*b, (a*c*cs30*cs_B + a*b*sn30*sn_C - c^2*cs30)*a, a*b*c*cs30]
    return [t_FQA,t_FQB,t_FQC]

def TriSecondNapoleonTriangle(a,b,c):
    sn30 = 1/2 ; cs30 = sqrt(3)/2  # sin(pi/6) and cos(pi/6)
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    [sn_A,sn_B,sn_C] = TriTriangleSines(a,b,c)
    return RT_TriSecondNapoleonTriangle(a,b,c,sn_A,cs_A,sn_B,cs_B,sn_C,cs_C,sn30,cs30)

# TriSecondFermatPoint returns 2nd Fermat point X14
# I don't use Kimberling definition, but mine as the intersection point of Simpson lines
def TriSecondFermatPoint(a,b,c):
    sn30 = 1/2 ; cs30 = sqrt(3)/2  # sin(pi/6) and cos(pi/6)
    [sn_A,sn_B,sn_C] = TriTriangleSines(a,b,c)
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    return [(b*c*cs30*cs_A + a*c*sn30*sn_B - b^2*cs30)*(a*b*cs30*cs_C + b*c*sn30*sn_A - a^2*cs30)*b*c, (a*c*cs30*cs_B + b*c*sn30*sn_A - a^2*cs30)*(a*b*cs30*cs_C + a*c*sn30*sn_B - b^2*cs30)*a*c, (a*b*cs30*cs_C + b*c*sn30*sn_A - a^2*cs30)*(a*b*cs30*cs_C + a*c*sn30*sn_B - b^2*cs30)*a*b]

# RT_TriSecondToricelliPoint(a,b,c,S) is same than TriSecondFermatPoint but using ETC definition with barycentric
def RT_SecondTriToricelliPoint(a,b,c,S):
    ks = 4/sqrt(3)*S
    kabc = a^2 + b^2 - c^2 - ks
    kbca = b^2 + c^2 - a^2 - ks
    kcab = c^2 + a^2 - b^2 - ks
    return [kabc*kcab/a,kabc*kbca/b,kcab*kbca/c]

# TriQuadranceSimpsonSegments2 returns quadrance of segments on Simpson lines from triangle vertices ABC to equilateral
# triangles vertices
def TriQuadranceSimpsonSegments2(a,b,c):
    return (1/2)*(a^2+b^2+c^2) - 2*sqrt(3)*TriArea(a,b,c)

# TriMorleyTriangle returns trilinear coordinates of first Morley triangle, given csArd = cos(A/3) ; csBrd = cos(B/3) ; csCrd = cos(C/3)
def TriMorleyTriangle(csArd,csBrd,csCrd):
    t_MTA = [1,2*csCrd,2*csBrd]
    t_MTB = [2*csCrd,1,2*csArd]
    t_MTC = [2*csBrd,2*csArd,1]
    return [t_MTA,t_MTB,t_MTC]

# TriTriangleCotangents returns cotangents of triangle LMN with ratio edge lengths l:m:n
def TriTriangleCotangents(l,m,n):
    S = TriArea(l,m,n)
    [S_L,S_M,S_N] = TriConwayParameters(l,m,n)
    ct_L = TriConwayTriangleCotangent(S,S_L)
    ct_M = TriConwayTriangleCotangent(S,S_M)
    ct_N = TriConwayTriangleCotangent(S,S_N)
    return [ct_L,ct_M,ct_N]

# TriConwayTriangle returns S_L,S_M,S_N for triangle LMN with ratio edge lengths l:m:n according to
# reference triangle ABC area S
def TriConwayTriangle(S,l,m,n):
    [ct_L,ct_M,ct_N] = TriTriangleCotangents(l,m,n)
    S_L = TriConwayTriangleNotation(S,ct_L)
    S_M = TriConwayTriangleNotation(S,ct_M)
    S_N = TriConwayTriangleNotation(S,ct_N)
    return [S_L,S_M,S_N]

# TriFirstGeneralizedFermatPoint returns 1st Generalized Fermat Point, point P minimizing l PA + m PB + n PC
# P is intersection point of lines AA1, BB1, CC1  where A1,B1,C1 are apex vertex of triangles constructed
# outward ABC on edges
def TriFirstGeneralizedFermatPoint(a,b,c,S,l,m,n):
    [S_A,S_B,S_C] = TriConwayParameters(a,b,c)
    [S_L,S_M,S_N] = TriConwayTriangle(S,l,m,n)
    tri_A1 = [    (-a^2)/a,        (S_C + S_N)/b,     (S_B + S_M)/c    ]
    tri_B1 = [  (S_C + S_N)/a,       (-b^2)/b,        (S_A + S_L)/c    ]
    tri_C1 = [  (S_B + S_M)/a,     (S_A + S_L)/b,        (-c^2)/c      ]
    tri_P  = [ (1/a)/(S_A + S_L), (1/b)/(S_B + S_M), (1/c)/(S_C + S_N) ]
    return [tri_P,tri_A1,tri_B1,tri_C1]

# TriSecondGeneralizedFermatPoint returns P intersection point of lines AA1, BB1, CC1  where A1,B1,C1 are
# apex vertex of triangles constructed outward ABC on edges
def TriSecondGeneralizedFermatPoint(a,b,c,S,l,m,n):
    [S_A,S_B,S_C] = TriConwayParameters(a,b,c)
    [S_L,S_M,S_N] = TriConwayTriangle(S,l,m,n)
    tri_A1 = [    (-a^2)/a,        (S_C - S_N)/b,     (S_B - S_M)/c    ]
    tri_B1 = [  (S_C - S_N)/a,       (-b^2)/b,        (S_A - S_L)/c    ]
    tri_C1 = [  (S_B - S_M)/a,     (S_A - S_L)/b,        (-c^2)/c      ]
    tri_P  = [ (1/a)/(S_A - S_L), (1/b)/(S_B - S_M), (1/c)/(S_C - S_N) ]
    return [tri_P,tri_A1,tri_B1,tri_C1]

# TriFeuerbachTriangle returns the intersection of Nine points circle and excircles
def TriFeuerbachTriangle(a,b,c):
    [sq_A,sq_B,sq_C] = TriSquaredCosinesBisector(a,b,c)
    t_FA = [sq_A - 1, sq_B, sq_C]
    t_FB = [sq_A, sq_B -1, sq_C]
    t_FC = [sq_A, sq_B, sq_C - 1]
    return [t_FA,t_FB,t_FC]

# TriFeuerbachPoint returns Feuerbach point X11, intersection of nine points circle and incircle
def TriFeuerbachPoint(a,b,c):
    return [b*c*(-a+b+c)*(b-c)^2, c*a*(a-b+c)*(c-a)^2, a*b*(a+b-c)*(a-b)^2]

# TriFirstIsodynamicPoint returns first isodynamic point X15, isogonal conjugate of the first Fermat Point
# r is inradius
def TriFirstIsodynamicPoint(a,b,c,r):
    [sn_A,sn_B,sn_C] = TriSines(a,b,c,r)
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    sn60 = sqrt(3)/2 ; cs60 = 1/2  # sin(pi/3) and cos(pi/3)
    return [sn_A*cs60 + sn60*cs_A,sn_B*cs60 + sn60*cs_B,sn_C*cs60 + sn60*cs_C]

# TriSecondIsodynamicPoint returns second isodynamic point X16, isogonal conjugate of the second Fermat Point
# r is inradius
def TriSecondIsodynamicPoint(a,b,c,r):
    [sn_A,sn_B,sn_C] = TriSines(a,b,c,r)
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    sn60 = sqrt(3)/2 ; cs60 = 1/2  # sin(pi/3) and cos(pi/3)
    return [sn_A*cs60 - sn60*cs_A,sn_B*cs60 - sn60*cs_B,sn_C*cs60 - sn60*cs_C]


# TriQuadranceCircumcenterIsodynamic returns  quadrance of segment from circumcenter to first isodynamic point
def TriQuadranceCircumcenterIsodynamic(a,b,c):
    [S_A,S_B,S_C] = TriConwayParameters(a,b,c)
    S = RT_TriangleArea(a,b,c)
    return (S_A^2 - S_A*S_B + S_B^2 - S_A*S_C - S_B*S_C + S_C^2)*(S_A + S_B)*(S_A + S_C)*(S_B + S_C)*S^2/((2*(S_A + S_B + S_C)*S + sqrt(3)*(S_A*S_B + S_A*S_C + S_B*S_C))^2*(S_A*S_B + S_A*S_C + S_B*S_C))

# TriQuadranceCircumcenter2ndIsodynamic returns  quadrance of segment from circumcenter to second isodynamic point
def TriQuadranceCircumcenter2ndIsodynamic(a,b,c):
    [S_A,S_B,S_C] = TriConwayParameters(a,b,c)
    S = RT_TriangleArea(a,b,c)
    return (S_A^2 - S_A*S_B + S_B^2 - S_A*S_C - S_B*S_C + S_C^2)*(S_A + S_B)*(S_A + S_C)*(S_B + S_C)*S^2/((-2*(S_A + S_B + S_C)*S + sqrt(3)*(S_A*S_B + S_A*S_C + S_B*S_C))^2*(S_A*S_B + S_A*S_C + S_B*S_C))

# TriIsodynamicEquation returns equation for inverse distances from O to isodynamic points
def TriIsodynamicEquation(a,b,c,x):
    S = RT_TriangleArea(a,b,c)
    [ctA,ctB,ctC] = TriConwayTriangleCotangents(a,b,c,S)
    [S_A,S_B,S_C] = TriConwayParameters(a,b,c)
    alpha = sqrt((S_A^2 - S_A*S_B + S_B^2 - S_A*S_C - S_B*S_C + S_C^2)*(S_A + S_B)*(S_A + S_C)*(S_B + S_C)/(S_A*S_B + S_A*S_C + S_B*S_C))
    beta = 4*(S_A + S_B + S_C)
    gamma = 4*(S_A*S_B + S_A*S_C + S_B*S_C)/((S_A + S_B)*(S_A + S_C)*(S_B + S_C))
    return x^2 - (beta/alpha)*x + (gamma)

# IsodynamicEquilateralTriangleArea returns area of equilateral triangle after inversion
def IsodynamicEquilateralTriangleArea(a,b,c,k):
    q = TriQuadranceCircumcenterIsodynamic(a,b,c)
    S = RT_TriangleArea(a,b,c)
    R = a*b*c/(4*S)
    area = (3*sqrt(3)/4)*(R^2/(q - R^2)^2)*k^4
    return area

# IsodynamicEquilateral2ndTriangleArea returns area of equilateral triangle after inversion
def IsodynamicEquilateral2ndTriangleArea(a,b,c,k):
    q = TriQuadranceCircumcenter2ndIsodynamic(a,b,c)
    S = RT_TriangleArea(a,b,c)
    R = a*b*c/(4*S)
    area = (3*sqrt(3)/4)*(R^2/(q - R^2)^2)*k^4
    return area

# IsodynamicSquaredRatio returns k^2 such as inversion leaves area unchanged
def IsodynamicSquaredRatio(a,b,c):
    q = TriQuadranceCircumcenterIsodynamic(a,b,c)
    S = RT_TriangleArea(a,b,c)
    k2 = 1/6*abs(a^2*b^2*c^2 - 16*S^2*q)*sqrt(sqrt(3)/S)/(a*b*c)
    return k2

# Isodynamic2ndSquaredRatio returns k^2 such as inversion leaves area unchanged
def Isodynamic2ndSquaredRatio(a,b,c):
    q = TriQuadranceCircumcenter2ndIsodynamic(a,b,c)
    S = RT_TriangleArea(a,b,c)
    k2 = 1/6*abs(a^2*b^2*c^2 - 16*S^2*q)*sqrt(sqrt(3)/S)/(a*b*c)
    return k2

# TriBevanPoint returns Bevan point X40, center of Bevan circle
def TriBevanPoint(a,b,c):
    return [b/(c+a-b)+c/(a+b-c)-a/(b+c-a),c/(a+b-c)+a/(b+c-a)-b/(c+a-b),a/(b+c-a)+b/(c+a-b)-c/(a+b-c)]

# TriFifthStevanovicPoint returns 5Th Stevanovic point X3021, on the incircle
def TriFifthStevanovicPoint(a,b,c):
    return [b*c*(-a+b+c)*(2*a^2-a*(b+c)+(b-c)^2)^2,c*a*(a-b+c)*(2*b^2-b*(c+a)+(c-a)^2)^2,a*b*(a+b-c)*(2*c^2-c*(a+b)*(a-b)^2)^2]

# TriSteinerPoint returns Steiner point X99
def TriSteinerPoint(a,b,c):
    t_S = [b*c*(a^2-b^2)*(a^2-c^2),c*a*(b^2-c^2)*(b^2-a^2),a*b*(c^2-a^2)*(c^2-b^2)]
    return t_S

# Tri_X110 returns triangle center X110,focus of Kieper parabola, named too Euler reflection point because it is
# intersection of reflection of Euler line on edges
def Tri_X110(a,b,c):
    t_X110 = [a/(b^2-c^2),b/(c^2-a^2),c/(a^2-b^2)]
    return t_X110

# Tri_X182 returns the midpoint of Brocard diaeter X3-X6 where X3 is circumcenter and X6 is Symmedian (Lemoine,Grebe) point
def Tri_X182(a,b,c):
    t_O = TriCircumCenter(a,b,c)
    t_K = TriSymmedianPoint(a,b,c)
    t_X182 =  TriMidpoint(a,b,c,t_O,t_K)
    return t_X182

# Tri_X351 returns triangle center X351, center of Parry circle
def Tri_X351(a,b,c):
    t_X351 = [a*(b^2-c^2)*(-2*a^2+b^2+c^2),b*(c^2-a^2)*(-2*b^2+c^2+a^2),c*(a^2-b^2)*(-2*c^2+a^2+b^2)]
    return t_X351

# TriParryPoint returns triangle center X111, ParryPoint
def TriParryPoint(a,b,c):
    t_P = [a/(2*a^2-b^2-c^2),b/(2*b^2-c^2-a^2),c/(2*c^2-a^2-b^2)]
    return t_P

# TriFarOutPoint return Far-Out point inverse of centroid G under inversion of circumcircle
def TriFarOutPoint(a,b,c):
    return [a*(b^4 + c^4 - a^4 - b^2*c^2),b*(c^4 + a^4 - b^4 - c^2*a^2),c*(a^4 + b^4 - c^4 - a^2*b^2)]

# TriParryCircleTriangle returns triangle FarOutPoint, Centroid
def TriParryCircleTriangle(a,b,c):
    t_G = TriCentroid(a,b,c)
    t_P = TriParryPoint(a,b,c)
    t_F = TriFarOutPoint(a,b,c)
    return [t_G,t_P,t_F]

# TriSchouteCenter return Schoute center X187 inverse of symmedian  K under inversion of circumcircle
def TriSchouteCenter(a,b,c):
    return [a*(2*a^2-b^2-c^2),b*(2*b^2-c^2-a^2),c*(2*c^2-a^2-b^2)]

# Tri_X190 returns center X190 Brianchon point of the Yff parabola
def Tri_X190(a,b,c):
    return [(a-b)*(c-a)/a,(a-b)*(b-c)/b,(b-c)*(c-a)/c]

# Tri_X648 returns center X648 trilinear pole of Euler line
def Tri_X648(a,b,c):
    return [(c^2-a^2)*(c^2+a^2-b^2)*(a^2-b^2)*(a^2+b^2-c^2)/a,(b^2-c^2)*(b^2+c^2-a^2)*(a^2-b^2)*(a^2+b^2-c^2)/b,(b^2-c^2)*(b^2+c^2-a^2)*(c^2-a^2)*(c^2+a^2-b^2)/c]


# TriLemoineAxis return Lemoine Axis, line through X351 center or Parry circle and X187 Schoute Center
def TriLemoineAxis(a,b,c):
    t_X357 = Tri_X351(a,b,c)
    t_Sc = TriSchouteCenter(a,b,c)
    l = TriLine(a,b,c,t_X357,t_Sc)
    return l

# TriLemoineTriangle returns triangle with circumcenter and two points intersection of circumcircle with LemoineAxis
def TriLemoineTriangle(a,b,c):
    k1 =  a^4 + b^4 + c^4 - a^2*b^2  - a^2*c^2 - b^2*c^2
    k2 = sqrt(k1)
    k3 = a^2 - b^2 + k2
    k4 = b^2 - c^2 - k2
    k5 = a^2 - b^2 - k2
    k6 = b^2 - c^2 + k2
    k7 = a^2 - c^2
    t_O = TriCircumCenter(a,b,c)
    t_L1 = [ a*k7*k4, -b*k3*k4, c*k3*k7 ]
    t_L2 = [ a*k7*k6, -b*k5*k6, c*k5*k7 ]
    return [t_O,t_L1,t_L2]

# TriKiepertHyperbolaAxis returns axis Kieper Hyperbola as two sets of three points on axis
# points are intersection of Simson lines with edge lines, for the two intersection points of
# Lemoine axis with circumcircle ABC
def TriKiepertHyperbolaAxis(a,b,c):
    k1 = sqrt(a^4 - a^2*b^2 + b^4 + c^4 - (a^2 + b^2)*c^2)
    k2 = a^2*b^2 - b^4 + a^2*c^2 - c^4 + k1*(-a^2 + b^2 + c^2)
    k3 = a^4 - a^2*b^2 + 2*b^4 - 2*a^2*c^2 - b^2*c^2 + c^4 - 2*k1*b^2
    k4 = a^4 + b^4 - a^2*c^2 - b^2*c^2 -k1*(a^2 + b^2 -c^2)
    k5 = b^2 - c^2 - k1
    k6 = a^2 - b^2 + k1
    k7 = a^2 - c^2
    k8 = b^2 - c^2 + k1
    k9 = a^2 - b^2 - k1
    k10 = a^4 - a^2*b^2 + 2*b^4 - 2*a^2*c^2 - b^2*c^2 + c^4 + 2*k1*b^2
    k11 = a^2*b^2 - b^4 + a^2*c^2 - c^4 - k1*(-a^2 + b^2 + c^2)
    k12 = a^4 + b^4 - a^2*c^2 - b^2*c^2 + k1*(a^2 + b^2 -c^2)
    #
    L1_A =  [0, c*k3*k5, b*k2*k7]
    L1_B =  [c*k4*k5, 0, a*k2*k6]
    L1_C =  [b*k4*k7, -a*k3*k6, 0]
    #
    L2_A =  [0, c*k10*k8, b*k11*k7]
    L2_B =  [c*k12*k8, 0, a*k11*k9]
    L2_C =  [b*k12*k7, -a*k10*k9, 0]
    #
    return [ [L1_A,L1_B,L1_C] , [L2_A,L2_B,L2_C] ]

# Tri_X3035 returns triangle center X3035, which is complement point of Feuerbach point
def Tri_X3035(a,b,c):
    #t_F = TriFeuerbachPoint(a,b,c)
    #t_X3035 = TriComplementPoint(a,b,c,t_F)
    t_X3035 = [b*c*((a-b+c)*(a-c)^2+(a+b-c)*(a-b)^2),c*a*((b-c+a)*(b-a)^2+(b+c-a)*(b-c)^2),a*b*((c-a+b)*(c-b)^2+(c+a-b)*(c-a)^2)]
    return t_X3035

# Tri_X3042 returns triangle center X3042 on Spieker circle
def Tri_X3042(a,b,c):
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    t_X3042 = [b*c*((a-b+c)*(a-c)^2*cs_B^2+(a+b-c)*(a-b)^2*cs_C^2),c*a*((b-c+a)*(b-a)^2*cs_C^2+(b+c-a)*(b-c)^2*cs_A^2),a*b*((c-a+b)*(c-b)^2*cs_A^2+(c+a-b)*(c-a)^2*cs_B^2)]
    return t_X3042

# TriSpiekerCircleTriangle returns triangle X3035, reflection of X3035 in Spieker center, X3042
# we choose these center triangles because defined for any values a,b,c among centers X3035 to X3042 on Spieker circle
def TriSpiekerCircleTriangle(a,b,c):
    t_S1 = Tri_X3035(a,b,c)
    t_Sp = TriSpiekerCenter(a,b,c)
    t_S2 = TriReflectionPoint(a,b,c,t_Sp,t_S1)
    t_S3 = Tri_X3042(a,b,c)
    return [t_S1,t_S2,t_S3]


# TriLineEuler returns the Euler line
# S_A, S_B, S_C is Conway Trianglenotation
def TriLineEuler(a,b,c):
    [S_A,S_B,S_C] = TriConwayParameters(a,b,c)
    l = [a*(b^2-c^2)*S_A,b*(c^2-a^2)*S_B,c*(a^2-b^2)*S_C]
    return l

# TriEulerLine returns Euler line OH
def TriEulerLine(a,b,c):
    l = [a*(b^2 - c^2)*(-a^2 + b^2 + c^2),b*(c^2 - a^2)*(a^2 - b^2 + c^2),c*(a^2 - b^2)*(a^2 + b^2 - c^2)]
    return l

# TriEulerLineCircumOrthoCenters returns two points O (circumcenter) and H (orthocenter) on the Euler Line for triangle t1-t2-t3
def TriEulerLineCircumOrthoCenters(a,b,c,t1,t2,t3):
    t4 =  TriTriangleCircumcenter(a,b,c,t1,t2,t3)
    t5 =  TriTriangleOrthocenter(a,b,c,t1,t2,t3)
    return [t4,t5]


# TriEulerTriangle returns the intersection of Euler line with edges lines
def TriEulerTriangle(a,b,c):
    [S_A,S_B,S_C] = TriConwayParameters(a,b,c)
    [f_A,f_B,f_C] = [a*S_A,b*S_B,c*S_C]
    t_EA = [0,(a^2-b^2)*f_C,(a^2-c^2)*f_B]   # intersect with BC line
    t_EB = [(a^2-b^2)*f_C,0,(c^2-b^2)*f_A]   # intersect with AC ine
    t_EC = [(c^2-a^2)*f_B,(c^2-b^2)*f_A,0]   # intersect with AB line
    return [t_EA,t_EB,t_EC]

# TriExcirclesCenters, returns trilinear coordinates of excircles for reference triangle ABC
# radii of excircles are 2S/(-a+b+c), 2S/(a-b+c) and 2S/(a+b-c) where S is area triangle ABC
# Excentral triangle
def TriExcirclesCenters(a,b,c):
    t_JA = [-1,1,1]   # center excircle opposite A, center JA
    t_JB = [1,-1,1]   # center excircle opposite B, center JB
    t_JC = [1,1,-1]   # center excircle opposite C, center JC
    return [t_JA,t_JB,t_JC]

# TriExtouchTriangle, returns contact triangle, intersection of excircles with edges
# It's contact of Mandard inellipse too
def TriExtouchTriangle(a,b,c):
    t_XA = [0,(-a+b+c)/a,(-a+b+c)/a] # intersect with BC line
    t_XB = [(a-b+c)/b,0,(a-b+c)/b]   # intersect with AC ine
    t_XC = [(a+b-c)/c,(a+b-c)/c,0]   # intersect with AB line
    return [t_XA,t_XB,t_XC]

# TriBevanTriangle, returns Bevan point of triangle, center of circumcircle excenters
def TriBevanTriangle(a,b,c):
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    return [cs_B+cs_C-cs_A-1,cs_C+cs_A-cs_B-1,cs_A+cs_B-cs_C-1]

# TriAntiAltimedialTriangle_A returns trilinear coordinates of AntiAltimedial triangle for A point
#  AA is reflection of A in BC line
#  BA is reflection of AA in B
#  CA is reflection of AA in C
def TriAntiAltimedialTriangle_A(a,b,c):
    t_AA =  [-a*b*c, (a^2 + b^2 - c^2)*c, (a^2 - b^2 + c^2)*b]
    t_BA =  [a*b*c, (a^2 - b^2 + c^2)*c, -(a^2 - b^2 + c^2)*b]
    t_CA =  [a*b*c, -(a^2 + b^2 - c^2)*c, (a^2 + b^2 - c^2)*b]
    return [t_AA,t_BA,t_CA]

# TriAntiAltimedialTriangle_B returns trilinear coordinates of AntiAltimedial triangle for B point
#  BB is reflection of B in CA line
#  CB is reflection of BB in C
#  AB is reflection of BB in A
def TriAntiAltimedialTriangle_B(a,b,c):
    t_BB =  [(a^2 + b^2 - c^2)*c, -a*b*c, -(a^2 - b^2 - c^2)*a]
    t_CB =  [-(a^2 + b^2 - c^2)*c, a*b*c, (a^2 + b^2 - c^2)*a]
    t_AB =  [-(a^2 - b^2 - c^2)*c, a*b*c, (a^2 - b^2 - c^2)*a]
    return [t_BB,t_CB,t_AB]

# TriAntiAltimedialTriangle_C returns trilinear coordinates of AntiAltimedial triangle for C point
#  CC is reflection of C in AB line
#  AC is reflection of CC in A
#  BC is reflection of CC in B
def TriAntiAltimedialTriangle_C(a,b,c):
    t_CC =  [(a^2 - b^2 + c^2)*b, -(a^2 - b^2 - c^2)*a, -a*b*c]
    t_AC =  [-(a^2 - b^2 - c^2)*b, (a^2 - b^2 - c^2)*a, a*b*c]
    t_BC =  [-(a^2 - b^2 + c^2)*b, (a^2 - b^2 + c^2)*a, a*b*c]
    return [t_CC,t_AC,t_BC]

# TriAntiAltimedialCubic returns coefficients of equation thru the three vertices of the three antialtimedial triangles
#   k3(a,b,c)*x^3 + k3(b,c,a)*y^3 + k3(c,a,b)*z^3 + k21(a,b,c)*x^2*y + k21(a,c,b)*x^2*z + k21(b,a,c)*y^2*x + k21(b,c,a)*y^2*z + k21(c,a,b)*z^2*x + k21(c,b,a)*z^2*y + k0(a,b,c)*x*y*z = 0
# we define helping functions TriAntiAltimedialCubic_k0,TriAntiAltimedialCubic_k3,TriAntiAltimedialCubic_k21

def TriAntiAltimedialCubic_k0(a,b,c):
    k0 = (3*a^36 - 20*a^34*b^2 + 59*a^32*b^4 - 108*a^30*b^6 + 148*a^28*b^8 - 184*a^26*b^10 + 324*a^24*b^12 - 788*a^22*b^14 + 1514*a^20*b^16 - 1896*a^18*b^18 + 1514*a^16*b^20 - 788*a^14*b^22 + 324*a^12*b^24 - 184*a^10*b^26 + 148*a^8*b^28 - 108*a^6*b^30 + 59*a^4*b^32 - 20*a^2*b^34 + 3*b^36 - 20*a^34*c^2 + 96*a^32*b^2*c^2 - 176*a^30*b^4*c^2 + 256*a^28*b^6*c^2 - 628*a^26*b^8*c^2 + 1092*a^24*b^10*c^2 - 576*a^22*b^12*c^2 - 656*a^20*b^14*c^2 + 612*a^18*b^16*c^2 + 612*a^16*b^18*c^2 - 656*a^14*b^20*c^2 - 576*a^12*b^22*c^2 + 1092*a^10*b^24*c^2 - 628*a^8*b^26*c^2 + 256*a^6*b^28*c^2 - 176*a^4*b^30*c^2 + 96*a^2*b^32*c^2 - 20*b^34*c^2 + 59*a^32*c^4 - 176*a^30*b^2*c^4 - 131*a^28*b^4*c^4 + 868*a^26*b^6*c^4 - 10*a^24*b^8*c^4 - 3000*a^22*b^10*c^4 + 3227*a^20*b^12*c^4 + 2116*a^18*b^14*c^4 - 5906*a^16*b^16*c^4 + 2116*a^14*b^18*c^4 + 3227*a^12*b^20*c^4 - 3000*a^10*b^22*c^4 - 10*a^8*b^24*c^4 + 868*a^6*b^26*c^4 - 131*a^4*b^28*c^4 - 176*a^2*b^30*c^4 + 59*b^32*c^4 - 108*a^30*c^6 + 256*a^28*b^2*c^6 + 868*a^26*b^4*c^6 - 4136*a^24*b^6*c^6 + 5076*a^22*b^8*c^6 + 1872*a^20*b^10*c^6 - 10084*a^18*b^12*c^6 + 6256*a^16*b^14*c^6 + 6256*a^14*b^16*c^6 - 10084*a^12*b^18*c^6 + 1872*a^10*b^20*c^6 + 5076*a^8*b^22*c^6 - 4136*a^6*b^24*c^6 + 868*a^4*b^26*c^6 + 256*a^2*b^28*c^6 - 108*b^30*c^6 + 148*a^28*c^8 - 628*a^26*b^2*c^8 - 10*a^24*b^4*c^8 + 5076*a^22*b^6*c^8 - 12141*a^20*b^8*c^8 + 8896*a^18*b^10*c^8 + 8147*a^16*b^12*c^8 - 18976*a^14*b^14*c^8 + 8147*a^12*b^16*c^8 + 8896*a^10*b^18*c^8 - 12141*a^8*b^20*c^8 + 5076*a^6*b^22*c^8 - 10*a^4*b^24*c^8 - 628*a^2*b^26*c^8 + 148*b^28*c^8 - 184*a^26*c^10 + 1092*a^24*b^2*c^10 - 3000*a^22*b^4*c^10 + 1872*a^20*b^6*c^10 + 8896*a^18*b^8*c^10 - 20616*a^16*b^10*c^10 + 11940*a^14*b^12*c^10 + 11940*a^12*b^14*c^10 - 20616*a^10*b^16*c^10 + 8896*a^8*b^18*c^10 + 1872*a^6*b^20*c^10 - 3000*a^4*b^22*c^10 + 1092*a^2*b^24*c^10 - 184*b^26*c^10 + 324*a^24*c^12 - 576*a^22*b^2*c^12 + 3227*a^20*b^4*c^12 - 10084*a^18*b^6*c^12 + 8147*a^16*b^8*c^12 + 11940*a^14*b^10*c^12 - 25929*a^12*b^12*c^12 + 11940*a^10*b^14*c^12 + 8147*a^8*b^16*c^12 - 10084*a^6*b^18*c^12 + 3227*a^4*b^20*c^12 - 576*a^2*b^22*c^12 + 324*b^24*c^12 - 788*a^22*c^14 - 656*a^20*b^2*c^14 + 2116*a^18*b^4*c^14 + 6256*a^16*b^6*c^14 - 18976*a^14*b^8*c^14 + 11940*a^12*b^10*c^14 + 11940*a^10*b^12*c^14 - 18976*a^8*b^14*c^14 + 6256*a^6*b^16*c^14 + 2116*a^4*b^18*c^14 - 656*a^2*b^20*c^14 - 788*b^22*c^14 + 1514*a^20*c^16 + 612*a^18*b^2*c^16 - 5906*a^16*b^4*c^16 + 6256*a^14*b^6*c^16 + 8147*a^12*b^8*c^16 - 20616*a^10*b^10*c^16 + 8147*a^8*b^12*c^16 + 6256*a^6*b^14*c^16 - 5906*a^4*b^16*c^16 + 612*a^2*b^18*c^16 + 1514*b^20*c^16 - 1896*a^18*c^18 + 612*a^16*b^2*c^18 + 2116*a^14*b^4*c^18 - 10084*a^12*b^6*c^18 + 8896*a^10*b^8*c^18 + 8896*a^8*b^10*c^18 - 10084*a^6*b^12*c^18 + 2116*a^4*b^14*c^18 + 612*a^2*b^16*c^18 - 1896*b^18*c^18 + 1514*a^16*c^20 - 656*a^14*b^2*c^20 + 3227*a^12*b^4*c^20 + 1872*a^10*b^6*c^20 - 12141*a^8*b^8*c^20 + 1872*a^6*b^10*c^20 + 3227*a^4*b^12*c^20 - 656*a^2*b^14*c^20 + 1514*b^16*c^20 - 788*a^14*c^22 - 576*a^12*b^2*c^22 - 3000*a^10*b^4*c^22 + 5076*a^8*b^6*c^22 + 5076*a^6*b^8*c^22 - 3000*a^4*b^10*c^22 - 576*a^2*b^12*c^22 - 788*b^14*c^22 + 324*a^12*c^24 + 1092*a^10*b^2*c^24 - 10*a^8*b^4*c^24 - 4136*a^6*b^6*c^24 - 10*a^4*b^8*c^24 + 1092*a^2*b^10*c^24 + 324*b^12*c^24 - 184*a^10*c^26 - 628*a^8*b^2*c^26 + 868*a^6*b^4*c^26 + 868*a^4*b^6*c^26 - 628*a^2*b^8*c^26 - 184*b^10*c^26 + 148*a^8*c^28 + 256*a^6*b^2*c^28 - 131*a^4*b^4*c^28 + 256*a^2*b^6*c^28 + 148*b^8*c^28 - 108*a^6*c^30 - 176*a^4*b^2*c^30 - 176*a^2*b^4*c^30 - 108*b^6*c^30 + 59*a^4*c^32 + 96*a^2*b^2*c^32 + 59*b^4*c^32 - 20*a^2*c^34 - 20*b^2*c^34 + 3*c^36)*a*b*c
    return k0

def TriAntiAltimedialCubic_k3(a,b,c):
    k3 = (a^24 - 4*a^22*b^2 + 4*a^20*b^4 + 8*a^18*b^6 - 31*a^16*b^8 + 48*a^14*b^10 - 44*a^12*b^12 + 32*a^10*b^14 - 33*a^8*b^16 + 36*a^6*b^18 - 24*a^4*b^20 + 8*a^2*b^22 - b^24 - 4*a^22*c^2 + 12*a^20*b^2*c^2 - 12*a^18*b^4*c^2 + 6*a^16*b^6*c^2 + 6*a^14*b^8*c^2 - 34*a^12*b^10*c^2 + 22*a^10*b^12*c^2 + 70*a^8*b^14*c^2 - 138*a^6*b^16*c^2 + 102*a^4*b^18*c^2 - 34*a^2*b^20*c^2 + 4*b^22*c^2 + 4*a^20*c^4 - 12*a^18*b^2*c^4 + 16*a^16*b^4*c^4 - 34*a^14*b^6*c^4 + 92*a^12*b^8*c^4 - 82*a^10*b^10*c^4 - 98*a^8*b^12*c^4 + 230*a^6*b^14*c^4 - 144*a^4*b^16*c^4 + 26*a^2*b^18*c^4 + 2*b^20*c^4 + 8*a^18*c^6 + 6*a^16*b^2*c^6 - 34*a^14*b^4*c^6 - 24*a^12*b^6*c^6 + 16*a^10*b^8*c^6 + 172*a^8*b^10*c^6 - 220*a^6*b^12*c^6 + 20*a^4*b^14*c^6 + 96*a^2*b^16*c^6 - 40*b^18*c^6 - 31*a^16*c^8 + 6*a^14*b^2*c^8 + 92*a^12*b^4*c^8 + 16*a^10*b^6*c^8 - 213*a^8*b^8*c^8 + 92*a^6*b^10*c^8 + 184*a^4*b^12*c^8 - 224*a^2*b^14*c^8 + 105*b^16*c^8 + 48*a^14*c^10 - 34*a^12*b^2*c^10 - 82*a^10*b^4*c^10 + 172*a^8*b^6*c^10 + 92*a^6*b^8*c^10 - 276*a^4*b^10*c^10 + 128*a^2*b^12*c^10 - 156*b^14*c^10 - 44*a^12*c^12 + 22*a^10*b^2*c^12 - 98*a^8*b^4*c^12 - 220*a^6*b^6*c^12 + 184*a^4*b^8*c^12 + 128*a^2*b^10*c^12 + 172*b^12*c^12 + 32*a^10*c^14 + 70*a^8*b^2*c^14 + 230*a^6*b^4*c^14 + 20*a^4*b^6*c^14 - 224*a^2*b^8*c^14 - 156*b^10*c^14 - 33*a^8*c^16 - 138*a^6*b^2*c^16 - 144*a^4*b^4*c^16 + 96*a^2*b^6*c^16 + 105*b^8*c^16 + 36*a^6*c^18 + 102*a^4*b^2*c^18 + 26*a^2*b^4*c^18 - 40*b^6*c^18 - 24*a^4*c^20 - 34*a^2*b^2*c^20 + 2*b^4*c^20 + 8*a^2*c^22 + 4*b^2*c^22 - c^24)*(a^2 + b^2 - c^2)*(a^2 - b^2 + b*c - c^2)*(a^2 - b^2 - b*c - c^2)*(a^2 - b^2 + c^2)*(a + b)*(a - b)*(a + c)*(a - c)*a^3
    return k3

def TriAntiAltimedialCubic_k21(a,b,c):
    k21 = (2*a^36 - 15*a^34*b^2 + 47*a^32*b^4 - 70*a^30*b^6 + 6*a^28*b^8 + 186*a^26*b^10 - 370*a^24*b^12 + 302*a^22*b^14 + 6*a^20*b^16 - 192*a^18*b^18 + 72*a^16*b^20 + 70*a^14*b^22 + 10*a^12*b^24 - 162*a^10*b^26 + 186*a^8*b^28 - 110*a^6*b^30 + 40*a^4*b^32 - 9*a^2*b^34 + b^36 - 12*a^34*c^2 + 75*a^32*b^2*c^2 - 205*a^30*b^4*c^2 + 337*a^28*b^6*c^2 - 329*a^26*b^8*c^2 - 109*a^24*b^10*c^2 + 1083*a^22*b^12*c^2 - 1675*a^20*b^14*c^2 + 751*a^18*b^16*c^2 + 789*a^16*b^18*c^2 - 879*a^14*b^20*c^2 - 413*a^12*b^22*c^2 + 1229*a^10*b^24*c^2 - 943*a^8*b^26*c^2 + 385*a^6*b^28*c^2 - 105*a^4*b^30*c^2 + 25*a^2*b^32*c^2 - 4*b^34*c^2 + 25*a^32*c^4 - 118*a^30*b^2*c^4 + 220*a^28*b^4*c^4 - 387*a^26*b^6*c^4 + 1107*a^24*b^8*c^4 - 1978*a^22*b^10*c^4 + 742*a^20*b^12*c^4 + 2891*a^18*b^14*c^4 - 4415*a^16*b^16*c^4 + 918*a^14*b^18*c^4 + 3368*a^12*b^20*c^4 - 3621*a^10*b^22*c^4 + 1249*a^8*b^24*c^4 + 202*a^6*b^26*c^4 - 250*a^4*b^28*c^4 + 45*a^2*b^30*c^4 + 2*b^32*c^4 - 8*a^30*c^6 + 250*a^26*b^4*c^6 - 721*a^24*b^6*c^6 - 190*a^22*b^8*c^6 + 3958*a^20*b^10*c^6 - 6830*a^18*b^12*c^6 + 2683*a^16*b^14*c^6 + 5840*a^14*b^16*c^6 - 8698*a^12*b^18*c^6 + 3462*a^10*b^20*c^6 + 2115*a^8*b^22*c^6 - 2898*a^6*b^24*c^6 + 1236*a^4*b^26*c^6 - 202*a^2*b^28*c^6 + 3*b^30*c^6 - 56*a^28*c^8 + 202*a^26*b^2*c^8 - 529*a^24*b^4*c^8 + 1850*a^22*b^6*c^8 - 3680*a^20*b^8*c^8 + 1197*a^18*b^10*c^8 + 7422*a^16*b^12*c^8 - 12576*a^14*b^14*c^8 + 6070*a^12*b^16*c^8 + 4563*a^10*b^18*c^8 - 8154*a^8*b^20*c^8 + 4956*a^6*b^22*c^8 - 1332*a^4*b^24*c^8 + 16*a^2*b^26*c^8 + 51*b^28*c^8 + 96*a^26*c^10 - 320*a^24*b^2*c^10 - 107*a^22*b^4*c^10 + 255*a^20*b^6*c^10 + 4113*a^18*b^8*c^10 - 9888*a^16*b^10*c^10 + 5383*a^14*b^12*c^10 + 7535*a^12*b^14*c^10 - 12887*a^10*b^16*c^10 + 7498*a^8*b^18*c^10 - 1593*a^6*b^20*c^10 - 533*a^4*b^22*c^10 + 703*a^2*b^24*c^10 - 255*b^26*c^10 + 28*a^24*c^12 + 514*a^22*b^2*c^12 - 18*a^20*b^4*c^12 - 3957*a^18*b^6*c^12 + 3836*a^16*b^8*c^12 + 6966*a^14*b^10*c^12 - 14736*a^12*b^12*c^12 + 7293*a^10*b^14*c^12 + 3072*a^8*b^16*c^12 - 3948*a^6*b^18*c^12 + 1440*a^4*b^20*c^12 - 1119*a^2*b^22*c^12 + 629*b^24*c^12 - 388*a^22*c^14 - 664*a^20*b^2*c^14 + 2344*a^18*b^4*c^14 + 2373*a^16*b^6*c^14 - 10358*a^14*b^8*c^14 + 6284*a^12*b^10*c^14 + 7140*a^10*b^12*c^14 - 10141*a^8*b^14*c^14 + 2742*a^6*b^16*c^14 + 990*a^4*b^18*c^14 + 768*a^2*b^20*c^14 - 1117*b^22*c^14 + 812*a^20*c^16 + 144*a^18*b^2*c^16 - 3713*a^16*b^4*c^16 + 4220*a^14*b^6*c^16 + 5110*a^12*b^8*c^16 - 12145*a^10*b^10*c^16 + 3590*a^8*b^12*c^16 + 3870*a^6*b^14*c^16 - 2804*a^4*b^16*c^16 - 502*a^2*b^18*c^16 + 1553*b^20*c^16 - 948*a^18*c^18 + 634*a^16*b^2*c^18 + 965*a^14*b^4*c^18 - 6679*a^12*b^6*c^18 + 5003*a^10*b^8*c^18 + 6792*a^8*b^10*c^18 - 5195*a^6*b^12*c^18 + 383*a^4*b^14*c^18 + 463*a^2*b^16*c^18 - 1643*b^18*c^18 + 702*a^16*c^20 - 338*a^14*b^2*c^20 + 2712*a^12*b^4*c^20 + 1925*a^10*b^6*c^20 - 8272*a^8*b^8*c^20 - 592*a^6*b^10*c^20 + 2380*a^4*b^12*c^20 + 223*a^2*b^14*c^20 + 1297*b^16*c^20 - 400*a^14*c^22 - 832*a^12*b^2*c^22 - 2862*a^10*b^4*c^22 + 2987*a^8*b^6*c^22 + 4688*a^6*b^8*c^22 - 1354*a^4*b^10*c^22 - 1074*a^2*b^12*c^22 - 815*b^14*c^22 + 296*a^12*c^24 + 1342*a^10*b^2*c^24 + 681*a^8*b^4*c^24 - 3184*a^6*b^6*c^24 - 930*a^4*b^8*c^24 + 920*a^2*b^10*c^24 + 465*b^12*c^24 - 280*a^10*c^26 - 864*a^8*b^2*c^26 + 483*a^6*b^4*c^26 + 1183*a^4*b^6*c^26 - 191*a^2*b^8*c^26 - 237*b^10*c^26 + 204*a^8*c^28 + 294*a^6*b^2*c^28 - 306*a^4*b^4*c^28 - 109*a^2*b^6*c^28 + 87*b^8*c^28 - 100*a^6*c^30 - 72*a^4*b^2*c^30 + 28*a^2*b^4*c^30 - 23*b^6*c^30 + 34*a^4*c^32 + 23*a^2*b^2*c^32 + 10*b^4*c^32 - 8*a^2*c^34 - 5*b^2*c^34 + c^36)*a^2*b
    return k21

def TriAntiAltimedialCubic(a,b,c):
    k3_abc = TriAntiAltimedialCubic_k3(a,b,c) ; k3_bca = TriAntiAltimedialCubic_k3(b,c,a) ; k3_cab = TriAntiAltimedialCubic_k3(c,a,b)
    k21_abc = TriAntiAltimedialCubic_k21(a,b,c) ; k21_acb = TriAntiAltimedialCubic_k21(a,c,b)
    k21_bac = TriAntiAltimedialCubic_k21(b,a,c) ; k21_bca = TriAntiAltimedialCubic_k21(b,c,a)
    k21_cab = TriAntiAltimedialCubic_k21(c,a,b) ; k21_cba = TriAntiAltimedialCubic_k21(c,b,a)
    k0_abc = TriAntiAltimedialCubic_k0(a,b,c)
    return [k3_abc,k3_bca,k3_cab,k21_abc,k21_acb,k21_bac,k21_bca,k21_cab,k21_cba,k0_abc]

def TriAntiAltimedialRationalCubic(a,b,c):
    [k3_abc,k3_bca,k3_cab,k21_abc,k21_acb,k21_bac,k21_bca,k21_cab,k21_cba,k0_abc] = TriAntiAltimedialCubic(a,b,c)
    m = gcd(gcd(gcd(gcd(gcd(gcd(gcd(gcd(gcd(ZZ(k3_abc),ZZ(k3_bca)),ZZ(k3_cab)),ZZ(k21_abc)),ZZ(k21_acb)),ZZ(k21_bac)),ZZ(k21_bca)),
                    ZZ(k21_cab)),ZZ(k21_cba)),ZZ(k0_abc))
    k3_abc /= m ; k3_bca /= m ; k3_cab /= m
    k21_abc /= m ; k21_acb /= m ; k21_bac /= m ; k21_bca /= m ; k21_cab /= m ; k21_cba /= m
    k0_abc /= m
    return [k3_abc,k3_bca,k3_cab,k21_abc,k21_acb,k21_bac,k21_bca,k21_cab,k21_cba,k0_abc]

# TriAntiAltimedialCubicEq returns cubic equation function h(x,y,z) = 0 with
# h(x,y,z) = k3(a,b,c)*x^3 + k3(b,c,a)*y^3 + k3(c,a,b)*z^3 + k21(a,b,c)*x^2*y + k21(a,c,b)*x^2*z + k21(b,a,c)*y^2*x + k21(b,c,a)*y^2*z + k21(c,a,b)*z^2*x + k21(c,b,a)*z^2*y + k0(a,b,c)*x*y*z
def TriAntiAltimedialCubicEq(a,b,c,t):
    [x,y,z] = t
    [k3_abc,k3_bca,k3_cab,k21_abc,k21_acb,k21_bac,k21_bca,k21_cab,k21_cba,k0_abc] = TriAntiAltimedialRationalCubic(a,b,c)
    eq = k3_abc*x^3 + k3_bca*y^3 + k3_cab*z^3 + k21_abc*x^2*y + k21_acb*x^2*z + k21_bac*y^2*x + k21_bca*y^2*z + k21_cab*z^2*x + k21_cba*z^2*y + k0_abc*x*y*z
    return eq

# TriCongruentIsoscelizerPoint returns the congruent isoscelizer point
def TriCongruentIsoscelizerPoint(a,b,c):
    [S_A,S_B,S_C] = TriConwayParameters(a,b,c)
    [f_A,f_B,f_C] = [a*S_A,b*S_B,c*S_C]  # TriCosines(a,b,c)
    [hf_A,hf_B,hf_C] = [sqrt(1+ f_A)/2,sqrt(1+ f_B)/2,sqrt(1+ f_C)/2]
    t_Ci = [hf_B+hf_C-hf_A,hf_A+hf_C-hf_B,hf_A+hf_B-hf_C]
    return t_Ci

# TriInnerSoddyCenter returns the inner Soddy center
def TriInnerSoddyCenter(a,b,c,r,s):
    R = a*b*c/(4*r*s) # circumradius of reference triangle
    t_Si = [1+b*c/(2*R*(-a+b+c)),1+a*c/(2*R*(a-b+c)),1+a*b/(2*R*(a+b-c))]
    return t_Si

# TriKyteCentroid returns centroid K of a kyte lamina ABDC where BC is symmetry segment a = BC, b = CA, c = AB
# computed simply are foot on BC triangle ABC centroid G,
def TriKyteCentroid(a,b,c):
    t_K = [0,c*(2*a^2 + b^2 - c^2), b*(2*a^2 - b^2 + c^2)]
    return t_K

# TriOuterSoddyCenter returns the outer Soddy center
def TriOuterSoddyCenter(a,b,c,r,s):
    R = a*b*c/(4*r*s) # circumradius of reference triangle
    t_So = [-1+b*c/(2*R*(-a+b+c)),-1+a*c/(2*R*(a-b+c)),-1+a*b/(2*R*(a+b-c))]
    return t_So

# TriOnLine check (if 0 is return) if one point given by trilinear coordinates on a line
# It's same as scalar product on trilinear coordinates
def TriOnLine(t,l):
    [la,lb,lc] = l
    [ta,tb,tc] = t
    on = (la*ta + lb*tb + lc*tc)
    return on

# TriLineIntersect return trilinear coordinates of point at intersection of two lines
def TriLineIntersect(a,b,c,l1,l2):
    [la1,lb1,lc1] = l1 ; [la2,lb2,lc2] = l2
    ta = lb1*lc2 - lb2*lc1 ; tb = lc1*la2 - lc2*la1 ; tc = la1*lb2 - la2*lb1
    u = a*ta + b*tb + c*tc
    t = [ta/u,tb/u,tc/u]
    return t

# TriAreLinesParallel returns 0 when lines are parallel
def TriAreLinesParallel(l1,l2):
    [la1,lb1,lc1] = l1 ; [la2,lb2,lc2] = l2
    u = a*(lb1*lc2 - lb2*lc1) + b*(lc1*la2 - lc2*la1) + c*(la1*lb2 - la2*lb1)
    return u

# TriExact converts trilinear coordinates to exact distances of the point to edge, where a,b,c,S are edge lengths
# and area of reference triangle
def TriExact(a,b,c,S,t):
    [ta,tb,tc] = t
    k = (2*S)/(a*ta+b*tb+c*tc)
    ea = k*ta ; eb = k*tb ; ec = k*tc
    e = [ea,eb,ec]
    return e

# TriQuadranceToLine returns quadrance of point t to line l
def TriQuadranceToLine(a,b,c,t,l):
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    [ta,tb,tc] = t
    [la,lb,lc] = l
    [ea,eb,ec] = TriExactCoordinates(a,b,c,t)
    q = (la*ea+lb*eb+lc*ec)^2/(la^2+lb^2+lc^2-2*lb*lc*cs_A-2*la*lc*cs_B-2*la*lb*cs_C)
    return q

# TriLineThreePoints returns three points on line l
def TriLineThreePoints(a,b,c,l):
    [tri_A,tri_B,tri_C] = TriMainTriangle(a,b,c)
    lAB = TriLine(a,b,c,tri_A,tri_B) ; lAC = TriLine(a,b,c,tri_A,tri_C) ; lBC = TriLine(a,b,c,tri_B,tri_C)
    tri_1 = TriLineIntersect(a,b,c,l,lAB)
    tri_2 = TriLineIntersect(a,b,c,l,lBC)
    tri_3 = TriLineIntersect(a,b,c,l,lAC)
    return [tri_1,tri_2,tri_3]

# TriLineParallelThruPoint returns point making with point t3, a parallel with line t1-t2
def TriLineParallelThruPoint(a,b,c,t1,t2,t3):
    t4 = TriReflectionPoint(a,b,c,t1,t2)
    t5 = TriReflectionPoint(a,b,c,t1,t3)
    t6 = TriReflectionPoint(a,b,c,t4,t5)
    return t6

#TriTriangleArea returns area of a sub-triangle of reference triangle, given its vertices as trilinear coordinates t1,t2,t3
def TriTriangleArea(a,b,c,S,t1,t2,t3):
    [ea1,eb1,ec1] = TriExact(a,b,c,S,t1) ; [ea2,eb2,ec2] = TriExact(a,b,c,S,t2) ; [ea3,eb3,ec3] = TriExact(a,b,c,S,t3)
    S123 = (a*b*c)*(ea1*(eb2*ec3-eb3*ec2) - eb1*(ea2*ec3-ea3*ec2) + ec1*(ea2*eb3-ea3*eb2))/(8*S^2)
    return S123

# TriArePerpendicularLines check (if 0 is return) if two lines are perpendicular
def TriArePerpendicularLines(a,b,c,l1,l2):
    [la1,lb1,lc1] = l1 ; [la2,lb2,lc2] = l2
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    perpendicular = (la1*la2 + lb1*lb2 + lc1*lc2) - (lb1*lc2 + lb2*lc1)*cs_A - (lc1*la2 + lc2*la1)*cs_B - (la1*lb2 + la2*lb1)*cs_C
    return perpendicular

# TriSegmentBissectorLine returns l perpendicular line to middle of segment (segment bissector)
def TriSegmentBissectorLine(a,b,c,t1,t2):
    t_12 = TriMidpoint(a,b,c,t1,t2)
    l_12 = TriLine(a,b,c,t1,t2)
    l = TriOrthogonalLine(a,b,c,l_12,t_12)
    return l

# TriLineParallel_BC returns one point parameter k on the line parallel with BC going thru point t
def TriLineParallel_BC(a,b,c,t,k):
    [ta,tb,tc] = t
    [sn_A,sn_B,sn_C] = TriTriangleSines(a,b,c)
    tka = ta
    tkb = tb - k*sn_C
    tkc = tc + k*sn_B
    return [tka,tkb,tkc]

# TriLineParallel_CA returns one point parameter k on the line parallel with CA going thru point t
def TriLineParallel_CA(a,b,c,t,k):
    [ta,tb,tc] = t
    [sn_A,sn_B,sn_C] = TriTriangleSines(a,b,c)
    tkb = tb
    tkc = tc - k*sn_A
    tka = ta + k*sn_C
    return [tka,tkb,tkc]

# TriLineParallel_AB returns one point parameter k on the line parallel with AB going thru point t
def TriLineParallel_AB(a,b,c,t,k):
    [ta,tb,tc] = t
    [sn_A,sn_B,sn_C] = TriTriangleSines(a,b,c)
    tkc = tc
    tka = ta - k*sn_B
    tkb = tb + k*sn_A
    return [tka,tkb,tkc]

# TriLineParallels returns three points parameter k on lines paralles with edges going thru point t
def TriLineParallels(a,b,c,t,k):
    tBC = TriLineParallel_BC(a,b,c,t,k)
    tCA = TriLineParallel_CA(a,b,c,t,k)
    tAB = TriLineParallel_AB(a,b,c,t,k)
    return [tBC,tCA,tAB]

# TriAtDistance_B_BC returns tc such as [0,-1,tc] at distance d from B
# Obsolete code (before november 2015), divide by zero
#def TriAtDistance_B_BC(a,b,c,r,d):
#    tc1 = b*d/(a*c + c*d)
#    tc2 = -b*d/(a*c - c*d)
#    return [tc1,tc2]
def TriAtDistance_B_BC(a,b,c,r,d):
    tb1 = a*c + c*d
    tc1 = -b*d
    tb2 = a*c - c*d
    tc2 = b*d
    return [[0,tb1,tc1],[0,tb2,tc2]]

# TriAtDistance_BC returns tc such as [0,-1,tc] at distance d from [0,-1,t]
# Obsolete code (before november 2015), divide by zero
#def TriAtDistance_BC(a,b,c,r,t,d):
#    tc1 = -(b^2*d + (a*b*c - b*c*d)*t)/(c^2*d*t - a*b*c - b*c*d)
#    tc2 = -(b^2*d - (a*b*c + b*c*d)*t)/(c^2*d*t + a*b*c - b*c*d)
#    return [tc1,tc2]
def TriAtDistance_BC(a,b,c,r,t,d):
    tb1 = c^2*d*t - a*b*c - b*c*d
    tc1 = b^2*d + (a*b*c - b*c*d)*t
    tb2 = c^2*d*t + a*b*c - b*c*d
    tc2 = b^2*d - (a*b*c + b*c*d)*t
    return [[0,tb1,tc1],[0,tb2,tc2]]

# TriAtDistance_C_CA returns tc such as [ta,0,-1] at distance d from C
def TriAtDistance_C_CA(a,b,c,r,d):
    A = a^4 + 2*a^2*b^2 + b^4 - 2*a^2*c^2 - 2*b^2*c^2 + c^4 - 4*a^2*d^2 + 4*a^2*r^2 + 8*a*b*r^2 + 4*b^2*r^2 + 8*a*c*r^2 + 8*b*c*r^2 + 4*c^2*r^2
    B = -4*a*c*d^2
    C = -2*c*d
    D = a^4 + 2*a^2*b^2 + b^4 - 2*a^2*c^2 - 2*b^2*c^2 + c^4 + 4*a^2*r^2 + 8*a*b*r^2 + 4*b^2*r^2 + 8*a*c*r^2 + 8*b*c*r^2 + 4*c^2*r^2
    # Obsolete code (before november 2015), divide by zero
    #ta1 = (B - C*sqrt(D))/A ; ta2 = (B + C*sqrt(D))/A
    ta1 = B - C*sqrt(D)
    ta2 = B + C*sqrt(D)
    tc = -A
    return [[ta1,0,tc],[ta2,0,tc]]

# TriAtDistance_CA returns ta such as [ta,0,-1] at distance d from [t,0,-1]
def TriAtDistance_CA(a,b,c,r,t,d):
    A = 4*a^4*d^2*t^2 - 8*a^3*c*d^2*t - a^4*c^2 - 2*a^2*b^2*c^2 - b^4*c^2 + 2*a^2*c^4 + 2*b^2*c^4 - c^6 + 4*a^2*c^2*d^2 - 4*a^2*c^2*r^2 - 8*a*b*c^2*r^2 - 4*b^2*c^2*r^2 - 8*a*c^3*r^2 - 8*b*c^3*r^2 - 4*c^4*r^2
    B = (4*a^3*d^2*t^2 - a^4*c*t - 2*a^2*b^2*c*t - b^4*c*t + 2*a^2*c^3*t + 2*b^2*c^3*t - c^5*t - 8*a^2*c*d^2*t - 4*a^2*c*r^2*t - 8*a*b*c*r^2*t - 4*b^2*c*r^2*t - 8*a*c^2*r^2*t - 8*b*c^2*r^2*t - 4*c^3*r^2*t + 4*a*c^2*d^2)*c
    C = 2*(a*t - c)^2*c*d
    D = a^4 + 2*a^2*b^2 + b^4 - 2*a^2*c^2 - 2*b^2*c^2 + c^4 + 4*a^2*r^2 + 8*a*b*r^2 + 4*b^2*r^2 + 8*a*c*r^2 + 8*b*c*r^2 + 4*c^2*r^2
    # Obsolete code (before november 2015), divide by zero
    #ta1 = (B + C*sqrt(D))/A ; ta2 = (B - C*sqrt(D))/A
    ta1 = B + C*sqrt(D)
    ta2 = B - C*sqrt(D)
    tc = -A
    return [[ta1,0,tc],[ta2,0,tc]]

# TriAtDistance_A_AB returns tb such as [-1,tb,0] at distance d from A
def TriAtDistance_A_AB(a,b,c,r,d):
    A = (a^4 - 2*a^2*b^2 + b^4 + 2*a^2*c^2 - 2*b^2*c^2 + c^4 - 4*a^2*d^2 + 4*a^2*r^2 + 8*a*b*r^2 + 4*b^2*r^2 + 8*a*c*r^2 + 8*b*c*r^2 + 4*c^2*r^2)*b
    B = -4*a^3*d^2
    C = -2*a^2*d
    D = a^4 - 2*a^2*b^2 + b^4 + 2*a^2*c^2 - 2*b^2*c^2 + c^4 + 4*a^2*r^2 + 8*a*b*r^2 + 4*b^2*r^2 + 8*a*c*r^2 + 8*b*c*r^2 + 4*c^2*r^2
    # Obsolete code (before november 2015), divide by zero
    #tb1 = (B - C*sqrt(D))/A ; tb2 = (B + C*sqrt(D))/A
    tb1 = B - C*sqrt(D)
    tb2 = B + C*sqrt(D)
    ta = -A
    return [[ta,tb1,0],[ta,tb2,0]]

# TriAtDistance_AB returns tb such as [-1,tb,0] at distance d from [-1,t,0]
def TriAtDistance_AB(a,b,c,r,t,d):
    A =  (4*b^2*d^2*t^2 - 8*a*b*d^2*t - a^4 + 2*a^2*b^2 - b^4 - 2*a^2*c^2 + 2*b^2*c^2 - c^4 + 4*a^2*d^2 - 4*a^2*r^2 - 8*a*b*r^2 - 4*b^2*r^2 - 8*a*c*r^2 - 8*b*c*r^2 - 4*c^2*r^2)*b
    B =  4*a*b^2*d^2*t^2 - a^4*b*t + 2*a^2*b^3*t - b^5*t - 2*a^2*b*c^2*t + 2*b^3*c^2*t - b*c^4*t - 8*a^2*b*d^2*t - 4*a^2*b*r^2*t - 8*a*b^2*r^2*t - 4*b^3*r^2*t - 8*a*b*c*r^2*t - 8*b^2*c*r^2*t - 4*b*c^2*r^2*t + 4*a^3*d^2
    C =  2*(b^2*d*t^2 - 2*a*b*d*t + a^2*d)
    D = a^4 - 2*a^2*b^2 + b^4 + 2*a^2*c^2 - 2*b^2*c^2 + c^4 + 4*a^2*r^2 + 8*a*b*r^2 + 4*b^2*r^2 + 8*a*c*r^2 + 8*b*c*r^2 + 4*c^2*r^2
    # Obsolete code (before november 2015), divide by zero
    #tb1 = (B - C*sqrt(D))/A ; tb2 = (B + C*sqrt(D))/A
    tb1 = B - C*sqrt(D)
    tb2 = B + C*sqrt(D)
    ta = -A
    return [[ta,tb1,0],[ta,tb2,0]]

# Limiting points for distant d = a circles , center B radius R = rb and center C radius r = rc
# See Wolfram - limiting points for algorithm
def TriLimitingPoints_BC(a,b,c,r,rb,rc):
    # Compute point M at distance (a^2-rc^2+rb^2) / (2*d) from B , with a = BC
    du = (a^2-rc^2+rb^2) / (2*a)
    [tu2,tu1] = TriAtDistance_B_BC(a,b,c,r,du)
    #print TriQuadrance(a,b,c,r,(a+b+c)/2,[0,1,0],[0,-1,tu1]) ,"vs", du^2
    # Computing points at quadrance q = ((a^2-rc^2+rb^2)^2 - 4a^2rb^2)/4a^2 from M
    #q = ((a^2-rc^2+rb2)^2 - 4*a^2*rb^2)/(4*a^2)  = du^2 - rb^2
    dv = sqrt(du^2 - rb^2)
    [tc1,tc2] = TriAtDistance_BC(a,b,c,r,tu1,dv)
    return [[0,-1,tc1],[0,-1,tu1],[0,-1,tc2]]

def TriLimitingPoints_CA(a,b,c,r,rc,ra):
    # Compute point M at distance (b^2-ra^2+rc^2) / (2*b) from C with b = CA
    du = (b^2-ra^2+rc^2) / (2*b)
    [tu2,tu1] = TriAtDistance_C_CA(a,b,c,r,du)
    #print TriQuadrance(a,b,c,r,(a+b+c)/2,[0,0,1],[tu1,0,-1]) ,"vs ", du^2
    # Computing points at quadrance q = ((b^2-ra^2+rc^2)^2 - 4b^2rc^2)/4b^2 from M
    #q = ((b^2-ra^2+rc^2)^2 - 4*b^2*rc^2)/(4*b^2)  = du^2 - rc^2
    dv = sqrt(du^2 - rc^2)
    [ta1,ta2] = TriAtDistance_CA(a,b,c,r,tu1,dv)
    return [[ta1,0,-1],[tu1,0,-1],[ta2,0,-1]]

def TriLimitingPoints_AB(a,b,c,r,ra,rb):
    # Compute point M at distance (c^2-ra^2+rc^2) / (2*c) from A with c = AB
    du = (c^2-rb^2+ra^2) / (2*c)
    [tu2,tu1] = TriAtDistance_A_AB(a,b,c,r,du)
    # Computing points at quadrance q = ((c^2-rb^2+ra^2)^2 - 4c^2ra^2)/4c^2 from M
    #q = ((c^2-rb^2+ra^2)^2 - 4*c^2*ra^2)/(4*c^2)  = du^2 - ra^2
    dv = sqrt(du^2 - ra^2)
    [tb1,tb2] = TriAtDistance_AB(a,b,c,r,tu1,dv)
    return [[-1,tb1,0],[-1,tu1,0],[-1,tb2,0]]

# Equation of circumcircle of reference triangle ABC
def TriOnCircumcircle(a,b,c,t):
    [ta,tb,tc] = t
    return c*ta*tb + a*tb*tc + b*tc*ta

# Points on circumcircle at same distance from point t than the vertices A,B,C
def TriCircumcircleMap(a,b,c,t):
    [ta,tb,tc] = t
    t_PA = [(a*c*ta - b*c*tb + a^2*tc - b^2*tc)*(a*b*ta + a^2*tb - c^2*tb - b*c*tc)*a,
            (a^2*c*tb + b^2*c*tb - c^3*tb - a^2*b*tc + b^3*tc - b*c^2*tc)*(a*b*ta + a^2*tb - c^2*tb - b*c*tc),
            -(a^2*c*tb + b^2*c*tb - c^3*tb - a^2*b*tc + b^3*tc - b*c^2*tc)*(a*c*ta - b*c*tb + a^2*tc - b^2*tc)]
    t_PB = [(a^2*c*ta + b^2*c*ta - c^3*ta + a^3*tc - a*b^2*tc - a*c^2*tc)*(b^2*ta - c^2*ta + a*b*tb - a*c*tc),
            -(a*c*ta - b*c*tb + a^2*tc - b^2*tc)*(b^2*ta - c^2*ta + a*b*tb - a*c*tc)*b,
            (a*c*ta - b*c*tb + a^2*tc - b^2*tc)*(a^2*c*ta + b^2*c*ta - c^3*ta + a^3*tc - a*b^2*tc - a*c^2*tc)]
    t_PC = [-(b^2*ta - c^2*ta + a*b*tb - a*c*tc)*(a^2*b*ta - b^3*ta + b*c^2*ta + a^3*tb - a*b^2*tb - a*c^2*tb),
            (a*b*ta + a^2*tb - c^2*tb - b*c*tc)*(a^2*b*ta - b^3*ta + b*c^2*ta + a^3*tb - a*b^2*tb - a*c^2*tb),
            (a*b*ta + a^2*tb - c^2*tb - b*c*tc)*(b^2*ta - c^2*ta + a*b*tb - a*c*tc)*c]
    return [t_PA,t_PB,t_PC]

# TriSimsonLine returns Simson line for a point on circumcircle of reference triangle
def TriSimsonLine(a,b,c,t):
    if TriOnCircumcircle(a,b,c,t) != 0:
        la = lb = lc = 0 # no line
    else:
        [ta,tb,tc] = t
        [sn_A,sn_B,sn_C] = TriTriangleSines(a,b,c)
        [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
        la = ta*(tb + tc*cs_A)/sn_A
        lb = tb*(tc + ta*cs_B)/sn_B
        lc = tc*(ta + tb*cs_C)/sn_C
    return [la,lb,lc]

# TriCevianTriangle returns cevian triangle for a point P with respect to reference triangle
# cevian triangle has vertices intersection of cevian lines (AP,BP,CP) with edges (BC,AC,AB)
def TriCevianTriangle(a,b,c,t):
    [ta,tb,tc] = t
    t1 = [0,            tb,           tc          ]  # on BC line
    t2 = [ta,           0,            tc          ]  # on CA line
    t3 = [ta,           tb,           0           ]  # on AB line
    return [t1,t2,t3]

# TriAntiCevianTriangle returns canti-cevian triangle for a point P with respect to reference triangle
# anticevian triangle is triangle with ABC as cevian triangle for a given point p
def TriAntiCevianTriangle(a,b,c,t):
    [ta,tb,tc] = t
    t1 = [-ta,          tb,           tc          ]
    t2 = [ta,          -tb,           tc          ]
    t3 = [ta,           tb,          -tc          ]
    return [t1,t2,t3]



# TriPedalTriangle returns pedal triangle for a point P with respect to reference triangle
# when point P is on the circumcircle of reference triangle, then the three points returned are colinear
# points are the feet of the point P on edges lines
def TriPedalTriangle(a,b,c,t):
    [ta,tb,tc] = t
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    t1 = [0,            tb + ta*cs_C, tc + ta*cs_B]  # on BC line
    t2 = [ta + tb*cs_C, 0,            tc + tb*cs_A]  # on CA line
    t3 = [ta + tc*cs_B, tb + tc*cs_A, 0           ]  # on AB line
    return [t1,t2,t3]

# TriAntiPedalTriangle returns antipedal triangle for a point P with respect to reference triangle
# points are the intersection of two perpendicular lines with PA,PB,PC thru A,B,C
def TriAntiPedalTriangle(a,b,c,t):
    [ta,tb,tc] = t
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    t1 = [-(tb + ta*cs_C)*(tc + ta*cs_B), (tc + ta*cs_B)*(ta + tb*cs_C), (tb + ta*cs_C)*(ta + tc*cs_B)]
    t2 = [ (tc + tb*cs_A)*(tb + ta*cs_C),-(tc + tb*cs_A)*(ta + tb*cs_C), (ta + tb*cs_C)*(tb + tc*cs_A)]
    t3 = [ (tb + tc*cs_A)*(tc + ta*cs_B), (ta + tc*cs_B)*(tc + tb*cs_A),-(ta + tc*cs_B)*(tb + tc*cs_A)]
    return [t1,t2,t3]

# TriPedalHomotheticCenter returns homothetic center when P doesn't lie on extended edges lines
def TriPedalHomotheticCenter(a,b,c,t):
    [ta,tb,tc] = t
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    return [a*ta*(ta + tb*cs_C)*(ta + tc*cs_B),b*tb*(tb + tc*cs_A)*(tb + ta*cs_C),c*tc*(tc + ta*cs_B)*(tb + tc*cs_A)]

# TriExcirclesRadicalCircle returns radical circle with circle function  (a^3 -ab^2 + 2abc -ac^2)/(4abc)
def TriExcirclesRadicalCircle(a,b,c):
    ga = (a^3 -a*b^2 + 2*a*b*c -a*c^2)/(4*a*b*c)
    gb = (b^3 -b*c^2 + 2*b*c*a -b*a^2)/(4*a*b*c)
    gc = (c^3 -c*a^2 + 2*c*a*b -c*b^2)/(4*a*b*c)
    g = [ga,gb,gc]
    return g

# TriExactCircle returns trilinear equation of circle center u and squared radius q
# equation is ka*(eta - eua)^2 + kb*(etb - eub)^2 + kc*(etc - euc)^2 = 1
# where [eta,etb,etc] and [eua,eub,euc] are "exact" coordinates (= exact distances)
# respectively of any point [ta,tb,tc] for center [ua,ub,uc]
def TriExactCircle(a,b,c,r,q):
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    [sn_A,sn_B,sn_C] = TriSines(a,b,c,r)
    ka = cs_A / (q*sn_B*sn_C) ; kb = cs_B / (q*sn_A*sn_C) ; kc = cs_C / (q*sn_A*sn_B)
    k = [ka,kb,kc] ; # equation coefficients
    return k

def TriExactOnCircle(a,b,c,S,u,k,t):
    [eua,eub,euc] = TriExact(a,b,c,S,u) ; # center point
    [ka,kb,kc] = k ; # circle
    [eta,etb,etc] = TriExact(a,b,c,S,t) ; # point to check
    return ka*(eta - eua)^2 + kb*(etb - eub)^2 + kc*(etc - euc)^2 - 1


# Equation of circumcircle three points in trilinear coordinates x:y:z
def TriCircumcircle(a,b,c,t1,t2,t3):
    [x1,y1,z1] = t1 ; [x2,y2,z2] = t2 ; [x3,y3,z3] = t3
    x,y,z = var('x,y,z')
    m = matrix(SR,4,4)
    m[0,0] = (a*y*z + b*z*x +c*x*y)/(a*x + b*y + c*z) ; m[0,1] = x ; m[0,2] = y ; m[0,3] = z
    m[1,0] = (a*y1*z1 + b*z1*x1 +c*x1*y1)/(a*x1 + b*y1 + c*z1) ; m[1,1] = x1 ; m[1,2] = y1 ; m[1,3] = z1
    m[2,0] = (a*y2*z2 + b*z2*x2 +c*x2*y2)/(a*x2 + b*y2 + c*z2) ; m[2,1] = x2 ; m[2,2] = y2 ; m[2,3] = z2
    m[3,0] = (a*y3*z3 + b*z3*x3 +c*x3*y3)/(a*x3 + b*y3 + c*z3) ; m[3,1] = x3 ; m[3,2] = y3 ; m[3,3] = z3
    return (m.determinant() * (a*x + b*y + c*z)).factor()

# Circumcircle functions l,m,n such as equation of circle : (lx + my + nz)(ax + by + cz) + (ayz + bxz + cxy) = 0
# function TriCircumcircleFactor returns p when equation is not normalized (l'x + m'y + n'z)(ax + by + cz) + p(ayz + bxz + cxy) = 0
# with l' = lp , m' = mp , n' = np
def TriCircumcircleFactor(a,b,c,t1,t2,t3):
    [x1,y1,z1] = t1 ; [x2,y2,z2] = t2 ; [x3,y3,z3] = t3
    mp = matrix(SR,3,3)
    mp[0,0] = x1 ; mp[0,1] = y1 ; mp[0,2] = z1
    mp[1,0] = x2 ; mp[1,1] = y2 ; mp[1,2] = z2
    mp[2,0] = x3 ; mp[2,1] = y3 ; mp[2,2] = z3
    p = mp.determinant()
    return p

# TriCircumcircleRefMidArcs return points at middle or ABC circumcircle arcs : BAC, CBA, ACB
def TriCircumcircleRefMidArcs(a,b,c):
    t_BAC = [a, c - b, b - c]
    t_CBA = [c - a, b, a - c]
    t_ACB = [b - a, a - b, c]
    return [t_BAC,t_CBA,t_ACB]

# TriCircumcircleOppositeRefMidArcs return points at middle or ABC circumcircle arcs : CAB, ABC, BCA
# they are reflected poins on circumcenter of points returned by TriCircumcircleRefMidArcs
def TriCircumcircleOppositeRefMidArcs(a,b,c):
    t_CAB = [-a/(b + c),1,1]
    t_ABC = [1,-b/(a + c),1]
    t_BCA = [1,1, -c/(a + b)]
    return [t_CAB,t_ABC,t_BCA]

# TriBrokenChordsHalver returns points splitting broken chords of ABC circumcircle in half length
def TriBrokenChordsHalver(a,b,c):
    t_BAC =  [(b + c)*c, 0, a*(b - c)] ; t_CAB = [(b + c)*b, a*(c - b), 0]  # Use t_BAC when c < b and t_CAB when b < c
    t_CBA =  [b*(c - a), (c + a)*a, 0] ; t_ABC = [0, (c + a)*c, b*(a - c)]  # Use t_CBA when a < c and t_ABC when c < a
    t_ACB =  [0, c*(a - b), (a + b)*b] ; t_BCA = [c*(b - a), 0, (a + b)*a]  # Use t_ACB when b < a and t_BCA when a < b
    return [t_BAC,t_CAB,t_CBA,t_ABC,t_ACB,t_BCA]

# TriCircumcircleRefFunctions return functions for special case cirumcircle of ABC triangle
def TriCircumcircleRefFunctions(a,b,c):
    return [0,0,0]

# TriCircumcircleFunctions returns functions for circumcircle of three points
def TriCircumcircleFunctions(a,b,c,t1,t2,t3):
    [x1,y1,z1] = t1 ; [x2,y2,z2] = t2 ; [x3,y3,z3] = t3
    p = TriCircumcircleFactor(a,b,c,t1,t2,t3)
    ml = matrix(SR,3,3)
    ml[0,0] = (a*y1*z1 + b*z1*x1 +c*x1*y1)/(a*x1 + b*y1 + c*z1) ; ml[0,1] = y1 ; ml[0,2] = z1
    ml[1,0] = (a*y2*z2 + b*z2*x2 +c*x2*y2)/(a*x2 + b*y2 + c*z2) ; ml[1,1] = y2 ; ml[1,2] = z2
    ml[2,0] = (a*y3*z3 + b*z3*x3 +c*x3*y3)/(a*x3 + b*y3 + c*z3) ; ml[2,1] = y3 ; ml[2,2] = z3
    l = ((-1/p)*ml.determinant())
    mm = matrix(SR,3,3)
    mm[0,0] = (a*y1*z1 + b*z1*x1 +c*x1*y1)/(a*x1 + b*y1 + c*z1) ; mm[0,1] = x1 ; mm[0,2] = z1
    mm[1,0] = (a*y2*z2 + b*z2*x2 +c*x2*y2)/(a*x2 + b*y2 + c*z2) ; mm[1,1] = x2 ; mm[1,2] = z2
    mm[2,0] = (a*y3*z3 + b*z3*x3 +c*x3*y3)/(a*x3 + b*y3 + c*z3) ; mm[2,1] = x3 ; mm[2,2] = z3
    m = ((1/p)*mm.determinant())
    mn = matrix(SR,3,3)
    mn[0,0] = (a*y1*z1 + b*z1*x1 +c*x1*y1)/(a*x1 + b*y1 + c*z1) ; mn[0,1] = x1 ; mn[0,2] = y1
    mn[1,0] = (a*y2*z2 + b*z2*x2 +c*x2*y2)/(a*x2 + b*y2 + c*z2) ; mn[1,1] = x2 ; mn[1,2] = y2
    mn[2,0] = (a*y3*z3 + b*z3*x3 +c*x3*y3)/(a*x3 + b*y3 + c*z3) ; mn[2,1] = x3 ; mn[2,2] = y3
    n = ((-1/p)*mn.determinant())
    return [l,m,n]

# TriConicSectionFunctions returns functions (u,v,w,f,g,h) of conic section defined by u x^2 + v y^2 + w z^2 + 2f yz + 2g zx + 2h xy = 0
# it is defined from five points on the conic
def TriConicSectionFunctions(a,b,c,t1,t2,t3,t4,t5):
    [x1,y1,z1] = t1 ; [x2,y2,z2] = t2 ; [x3,y3,z3] = t3 ; [x4,y4,z4] = t4 ; [x5,y5,z5] = t5
    mu = matrix(SR,5,5)
    mu[0,0] = y1^2 ; mu[0,1] = z1^2 ; mu[0,2] = y1*z1 ; mu[0,3] = z1*x1 ; mu[0,4] = x1*y1
    mu[1,0] = y2^2 ; mu[1,1] = z2^2 ; mu[1,2] = y2*z2 ; mu[1,3] = z2*x2 ; mu[1,4] = x2*y2
    mu[2,0] = y3^2 ; mu[2,1] = z3^2 ; mu[2,2] = y3*z3 ; mu[2,3] = z3*x3 ; mu[2,4] = x3*y3
    mu[3,0] = y4^2 ; mu[3,1] = z4^2 ; mu[3,2] = y4*z4 ; mu[3,3] = z4*x4 ; mu[3,4] = x4*y4
    mu[4,0] = y5^2 ; mu[4,1] = z5^2 ; mu[4,2] = y5*z5 ; mu[4,3] = z5*x5 ; mu[4,4] = x5*y5
    u = mu.determinant()
    mv = matrix(SR,5,5)
    mv[0,0] = x1^2 ; mv[0,1] = z1^2 ; mv[0,2] = y1*z1 ; mv[0,3] = z1*x1 ; mv[0,4] = x1*y1
    mv[1,0] = x2^2 ; mv[1,1] = z2^2 ; mv[1,2] = y2*z2 ; mv[1,3] = z2*x2 ; mv[1,4] = x2*y2
    mv[2,0] = x3^2 ; mv[2,1] = z3^2 ; mv[2,2] = y3*z3 ; mv[2,3] = z3*x3 ; mv[2,4] = x3*y3
    mv[3,0] = x4^2 ; mv[3,1] = z4^2 ; mv[3,2] = y4*z4 ; mv[3,3] = z4*x4 ; mv[3,4] = x4*y4
    mv[4,0] = x5^2 ; mv[4,1] = z5^2 ; mv[4,2] = y5*z5 ; mv[4,3] = z5*x5 ; mv[4,4] = x5*y5
    v = (-1)*mv.determinant()
    mw = matrix(SR,5,5)
    mw[0,0] = x1^2 ; mw[0,1] = y1^2 ; mw[0,2] = y1*z1 ; mw[0,3] = z1*x1 ; mw[0,4] = x1*y1
    mw[1,0] = x2^2 ; mw[1,1] = y2^2 ; mw[1,2] = y2*z2 ; mw[1,3] = z2*x2 ; mw[1,4] = x2*y2
    mw[2,0] = x3^2 ; mw[2,1] = y3^2 ; mw[2,2] = y3*z3 ; mw[2,3] = z3*x3 ; mw[2,4] = x3*y3
    mw[3,0] = x4^2 ; mw[3,1] = y4^2 ; mw[3,2] = y4*z4 ; mw[3,3] = z4*x4 ; mw[3,4] = x4*y4
    mw[4,0] = x5^2 ; mw[4,1] = y5^2 ; mw[4,2] = y5*z5 ; mw[4,3] = z5*x5 ; mw[4,4] = x5*y5
    w = mw.determinant()
    mf = matrix(SR,5,5)
    mf[0,0] = x1^2 ; mf[0,1] = y1^2 ; mf[0,2] = z1^2 ; mf[0,3] = z1*x1 ; mf[0,4] = x1*y1
    mf[1,0] = x2^2 ; mf[1,1] = y2^2 ; mf[1,2] = z2^2 ; mf[1,3] = z2*x2 ; mf[1,4] = x2*y2
    mf[2,0] = x3^2 ; mf[2,1] = y3^2 ; mf[2,2] = z3^2 ; mf[2,3] = z3*x3 ; mf[2,4] = x3*y3
    mf[3,0] = x4^2 ; mf[3,1] = y4^2 ; mf[3,2] = z4^2 ; mf[3,3] = z4*x4 ; mf[3,4] = x4*y4
    mf[4,0] = x5^2 ; mf[4,1] = y5^2 ; mf[4,2] = z5^2 ; mf[4,3] = z5*x5 ; mf[4,4] = x5*y5
    f = (-1)*mf.determinant()
    mg = matrix(SR,5,5)
    mg[0,0] = x1^2 ; mg[0,1] = y1^2 ; mg[0,2] = z1^2 ; mg[0,3] = y1*z1 ; mg[0,4] = x1*y1
    mg[1,0] = x2^2 ; mg[1,1] = y2^2 ; mg[1,2] = z2^2 ; mg[1,3] = y2*z2 ; mg[1,4] = x2*y2
    mg[2,0] = x3^2 ; mg[2,1] = y3^2 ; mg[2,2] = z3^2 ; mg[2,3] = y3*z3 ; mg[2,4] = x3*y3
    mg[3,0] = x4^2 ; mg[3,1] = y4^2 ; mg[3,2] = z4^2 ; mg[3,3] = y4*z4 ; mg[3,4] = x4*y4
    mg[4,0] = x5^2 ; mg[4,1] = y5^2 ; mg[4,2] = z5^2 ; mg[4,3] = y5*z5 ; mg[4,4] = x5*y5
    g = mg.determinant()
    mh = matrix(SR,5,5)
    mh[0,0] = x1^2 ; mh[0,1] = y1^2 ; mh[0,2] = z1^2 ; mh[0,3] = y1*z1 ; mh[0,4] = z1*x1
    mh[1,0] = x2^2 ; mh[1,1] = y2^2 ; mh[1,2] = z2^2 ; mh[1,3] = y2*z2 ; mh[1,4] = z2*x2
    mh[2,0] = x3^2 ; mh[2,1] = y3^2 ; mh[2,2] = z3^2 ; mh[2,3] = y3*z3 ; mh[2,4] = z3*x3
    mh[3,0] = x4^2 ; mh[3,1] = y4^2 ; mh[3,2] = z4^2 ; mh[3,3] = y4*z4 ; mh[3,4] = z4*x4
    mh[4,0] = x5^2 ; mh[4,1] = y5^2 ; mh[4,2] = z5^2 ; mh[4,3] = y5*z5 ; mh[4,4] = z5*x5
    h = (-1)*mh.determinant()
    return [u,v,w,f,g,h]

# TriConicCenter returns center of conic (intersection of axes)
# coordinates are trilinear coordinates up to a constant factor, because equation of conic
# is up to a constant
def TriConicCenter(a,b,c,u,v,w,f,g,h):
    ta = a*f^2 - b*f*g + c*v*g - a*v*w - c*f*h + b*w*h
    tb = b*g^2 - c*g*h + a*w*h - b*w*u - a*g*f + c*u*f
    tc = c*h^2 - a*h*f + b*u*f - c*u*v - b*h*g + a*v*g
    t = [ta,tb,tc]
    return t

# TriOnCircle check if point t is on circle given by functions f
def TriOnCircle(a,b,c,t,f):
    [ta,tb,tc] = t ; [fl,fm,fn] = f
    return (fl*ta+fm*tb+fn*tc)*(a*ta+b*tb+c*tc)+(a*tb*tc+b*tc*ta+c*ta*tb)

# TriCircleCenter returns center of circle given by functions f
def TriCircleCenter(a,b,c,f):
    [fl,fm,fn] = f
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    ta = fl + cs_A - fn*cs_B - fm*cs_C
    tb = fm - fn*cs_A + cs_B - fl*cs_C
    tc = fn - fm*cs_A - fl*cs_B + cs_C
    return [ta,tb,tc]

# TriCirclesEdgesFunctions returns functions for the circles thru every couple of vertices reference triangle and a point t
def TriCirclesEdgesFunctions(a,b,c,t):
    [ta,tb,tc] = t
    fAB = [0, 0, -(c*ta*tb + b*ta*tc + a*tb*tc)/((a*ta + b*tb + c*tc)*tc)]
    fBC = [-(c*ta*tb + b*ta*tc + a*tb*tc)/((a*ta + b*tb + c*tc)*ta), 0, 0]
    fCA = [0, -(c*ta*tb + b*ta*tc + a*tb*tc)/((a*ta + b*tb + c*tc)*tb), 0]
    return [fAB,fBC,fCA]

# TriCircleLineIntersect returns two chord points, intersection of circle with line (if they exist)
#def TriCircleLineIntersect_Old(a,b,c,f,l):
#    ta,tb,tc = var('ta,tb,tc')
#    t = [ta,tb,tc]
#    ta1 = 1
#    res = solve([TriOnCircle(a,b,c,t,f) == 0,TriOnLine(t,l) == 0],tb,tc)
#    if len(res) == 2 :
#        ub = ((res[0][0].right() + res[1][0].right())/(2*ta)).factor() ; print "ub = ",ub  # removing ta
#        uc = ((res[0][1].right() + res[1][1].right())/(2*ta)).factor() ; print "uc = ",uc  # removing ta
#        vb = (res[0][0].right()/ta).factor() - ub ; sq_vb = vb^2
#        vc = (res[0][1].right()/ta).factor() - uc ; sq_vc = vc^2
#    else:
#        if len(res) == 1 :
#            ub = (res[0][0].right()/ta).factor()
#            uc = (res[0][1].right()/ta).factor()
#            sq_vb = 0
#            sq_vc = 0
#        else :
#            ub = uc = sq_vb = sq_vc = ta1 = 0 # returns [0,0,0] when no intersection
#    tb1 = ub - sqrt(sq_vb) ; tc1 = uc - sqrt(sq_vc)
#    tb2 = ub + sqrt(sq_vb) ; tc2 = uc + sqrt(sq_vc)
#    return[ [ta1,tb1,tc1], [ta1,tb2,tc2] ]

def TriCircleLineIntersect(a,b,c,f,l):
    ta,tb,tc = var('ta,tb,tc')
    t = [ta,tb,tc]
    res = solve([TriOnCircle(a,b,c,t,f) == 0,TriOnLine(t,l) == 0],tb,tc)
    if len(res) == 2 :
        ta1 = 1 ; tb1 = res[0][0].right()/ta ; tc1 = res[0][1].right()/ta
        ta2 = 1 ; tb2 = res[1][0].right()/ta ; tc2 = res[1][1].right()/ta
    else:
        if len(res) == 1 :
            ta1 = 1 ; tb1 = res[0][0].right()/ta ; tc1 = res[0][1].right()/ta
            ta2 = ta1 ; tb2 = tb1 ; tc2 = tc1
        else:
            ta1 = tb1 = tc1 = ta2 = tb2 = tc2 = 0
    return [ [ta1,tb1,tc1], [ta1,tb2,tc2] ]

# TriCirclesRadicaline returns radical line for two circles
# if circles intersect, then radical line goes thru a chord
def TriCirclesRadicalLine(a,b,c,f1,f2):
    [f1l,f1m,f1n] = f1 ; [f2l,f2m,f2n] = f2
    return [f1l-f2l,f1m-f2m,f1n-f2n]

# TriCirclesIntersect returns intersection points of two circles
def TriCirclesIntersect(a,b,c,f1,f2):
    l = TriCirclesRadicalLine(a,b,c,f1,f2)
    [tri1,tri2] = TriCircleLineIntersect(a,b,c,f1,l)
    return [tri1,tri2]

# TriCirclesExternalSimilitudeCenter returns point of external similitude for two circles centered at t1,t2
# and with radii r1,r2
def TriCirclesExternalSimilitudeCenter(a,b,c,t1,t2,r1,r2):
    [ta1,tb1,tc1] = t1 ; [ta2,tb2,tc2] = t2
    u1 = a*ta1 + b*tb1 + c*tc1 ; u2 = a*ta2 + b*tb2 + c*tc2
    ta = ( a*r1*u1*ta2 - a*r2*u2*ta1 )/a
    tb = ( b*r1*u1*tb2 - b*r2*u2*tb1 )/b
    tc = ( c*r1*u1*tc2 - c*r2*u2*tc1 )/c
    return [ta,tb,tc]

# TriCirclesInternalSimilitudeCenter returns point of internal similitude for two circles centered at t1,t2
# and with radii r1,r2
def TriCirclesInternalSimilitudeCenter(a,b,c,t1,t2,r1,r2):
    [ta1,tb1,tc1] = t1 ; [ta2,tb2,tc2] = t2
    u1 = a*ta1 + b*tb1 + c*tc1 ; u2 = a*ta2 + b*tb2 + c*tc2
    ta = ( a*r1*u1*ta2 + a*r2*u2*ta1 )/a
    tb = ( b*r1*u1*tb2 + b*r2*u2*tb1 )/b
    tc = ( c*r1*u1*tc2 + c*r2*u2*tc1 )/c
    return [ta,tb,tc]

# TriCircumConicPoints returns two points on circumconic equation fyz + gzx + hxy = 0 image of isogonal conjugate of
# a line equation fx + gy + hz = 0  where t1 and t2 are two points on that line
# three other points are vertices A,B,C of reference triangle
def TriCircumConicPoints(a,b,c,t1,t2):
    tri_1 =  TriIsogonalConjugatePoint(a,b,c,t1)
    tri_2 =  TriIsogonalConjugatePoint(a,b,c,t2)
    return [tri_1,tri_2]


# TriInconicFunctions returns functions f of inconic given Brianchon point coordinates
def TriInconicFunctions(a,b,c,t):
    [ta,tb,tc] = t
    fl = 1/ta ; fm = 1/tb ; fn = 1/tc
    return [fl,fm,fn]

# Common inconics
def TriBrocardInellipse(a,b,c):
    t_K = TriSymmedianPoint(a,b,c)  # Brocard point is symmedian point
    f = TriInconicFunctions(a,b,c,t_K)
    return f

def TriKiepertParabola(a,b,c):
    t_S = TriSteinerPoint(a,b,c)  # Brocard point is Steiner point
    f = TriInconicFunctions(a,b,c,t_S)
    return f

def TriMandartInellipse(a,b,c):
    t_N = TriNagelPoint(a,b,c)  # Brocard point is Nagel point
    f = TriInconicFunctions(a,b,c,t_N)
    return f

def TriSteinerInellipse(a,b,c):
    t_G = TriCentroid(a,b,c)  # Brocard point is Centroid
    f = TriInconicFunctions(a,b,c,t_G)
    return f

# TriIsInconicParabola check if inconic given by function f is parabola  (returns 0 in that case)
def TriIsInconicParabola(a,b,c,f):
    [fl,fm,fn] = f
    return fl/a + fm/b + fn/c

# TriInconicCenter returns center of inconic given by functions f
def TriInconicCenter(a,b,c,f):
    [fl,fm,fn] = f
    ta = c*fm + b*fn ; tb = a*fn + c*fl ; tc = b*fl + a*fm
    return [ta,tb,tc]

# TriInconicFunctions returns functions f from center t
# inverse of TriInconicCenter
# ta =  c fm + b fn => a ta = ac fm + ab fn =>  (a ta - b tb)/c = a fm - b fl
# tb =  a fn + c fl => b tb = ab fn + bc fl
# tc =  b fl + a fm
# => 2 a fm = (a ta - b tb)/c + tc => fm = (a ta - b tb + c tc) / (2ac)
# => 2 b fl = (-a ta + b tb)/c + tc => fl = (-a ta + b tb + c tc) / (2bc)
# fn = (a ta + b tb - c tc) / (2ab)

def TriInconicFunctions(a,b,c,t):
    [ta,tb,tc] = t
    fl = (-a*ta + b*tb+ c*tc)/(2*b*c) ; fm = (-b*tb + c*tc + a*ta)/(2*c*a) ; fn = (-c*tc + a*ta + b*tb)/(2*a*b)
    return [fl,fm,fn]

# TriOnInconic check if point t is on inconic given by functions f
def TriOnInconic(a,b,c,t,f):
    [ta,tb,tc] = t ; [fl,fm,fn] = f
    return fl^2*ta^2 + fm^2*tb^2 + fn^2*tc^2 - 2*fm*fn*tb*tc - 2*fn*fl*tc*ta - 2*fl*fm*ta*tb

# TriInconicLineIntersect returns two chord points, intersection of inconic with line (if they exist)
def TriInconicLineIntersect(a,b,c,f,l):
    ta,tb,tc = var('ta,tb,tc')
    t = [ta,tb,tc]
    ta1 = 1
    res = solve([TriOnInconic(a,b,c,t,f) == 0,TriOnLine(t,l) == 0],tb,tc)
    if len(res) == 2 :
        ub = (res[0][0].right() + res[1][0].right())/(2*ta)
        uc = (res[0][1].right() + res[1][1].right())/(2*ta)
        sq_vb = res[0][0].right()/ta - ub
        sq_vc = res[0][1].right()/ta - uc
    else:
        if len(res) == 1 :
            ub = res[0][0].right()/ta
            uc = res[0][1].right()/ta
            sq_vb = 0
            sq_vc = 0
        else :
            ub = uc = sq_vb = sq_vc = ta1 = 0 # returns [0,0,0] when no intersection
    tb1 = ub - sqrt(sq_vb)
    tc1 = uc - sqrt(sq_vc)
    tb2 = ub + sqrt(sq_vb)
    tc2 = uc + sqrt(sq_vc)
    return[ [ta1,tb1,tc1], [ta1,tb2,tc2] ]

# TriConicPoint returns point on conic according to Colin MacLaurin construction of conic
#  tri_1,tri_2,tri_3 : three fixed given points
#  tri_4,tri_5 : two points defining line (tri_4,tri_5)
#  tri_6 : a varying vertice on a second line (tri_4,tri_6) of the varying triangle
def TriConicPoint(a,b,c,tri_1,tri_2,tri_3,tri_4,tri_5,tri_6):
    l45 = TriLine(a,b,c,tri_4,tri_5)                                       # given line
    l16 = TriLine(a,b,c,tri_1,tri_6) ; tri_7 = TriLineIntersect(a,b,c,l45,l16)   # set a line thru varying vertice to get new vertice
    l26 = TriLine(a,b,c,tri_2,tri_6) ; l37 = TriLine(a,b,c,tri_3,tri_7)          # intersect two lines to get last vertice
    tri_8 = TriLineIntersect(a,b,c,l26,l37)
    return [tri_8,tri_7,tri_6]

# Steiner circumellipse equation is 1/(ax) + 1/(by) + (1/cz) = 0 where x:y:z are trilinear coordinates
def TriIsOnSteinerCircumellipse(a,b,c,t):
    [ta,tb,tc] = t
    return 1/(a*ta) + 1/(b*tb) + 1/(c*tc)

# Steiner circumellipse center is centroid G
def TriSteinerCircumellipseCenter(a,b,c,t):
    return TriCentroid(a,b,c)

# TriSteinerCircumellipseLineIntersect returns intersection of circumellipse with line
def TriSteinerCircumellipseLineIntersectLine(a,b,c,l):
    [la,lb,lc] = l
    K = sqrt(b^2*c^2*la^2 - 2*a*b*c^2*la*lb + a^2*c^2*lb^2 - 2*a*b^2*c*la*lc - 2*a^2*b*c*lb*lc + a^2*b^2*lc^2)
    ta1 = -2*b*c*lb*lc
    tb1 = lc*(b*c*la + a*c*lb - a*b*lc + K)
    tc1 = lb*(b*c*la - a*c*lb + a*b*lc - K)
    ta2 = -2*b*c*lb*lc
    tb2 = lc*(b*c*la + a*c*lb - a*b*lc - K)
    tc2 = lb*(b*c*la - a*c*lb + a*b*lc + K)
    return [[ta1,tb1,tc1],[ta2,tb2,tc2]]

# TriTriangleSteinerCircumellipseBickartPoints returns Bickart points (= foci) of ABC circumellipse
# OBSOLETE
#def TriTriangleSteinerCircumellipseBickartPoints(a,b,c):
#    qa = a^2 ; qb = b^2 ; qc = c^2
#    qd = sqrt( 1/2*( (qb - qc)^2 + (qc - qa)^2 + (qa - qb)^2) )
#    p2 = (1/9)*(qa + qb + qc + 2*qd) #  squared semimajor axis length
#    q2 = (1/9)*(qa + qb + qc - 2*qd) #  squared semiminor axis length
#    if qb + qc - 2*qa == 0:
#        t = pi/2
#    else:
#        t = (1/2)*atan(sqrt(3)*(qc - qb)/(qb + qc - 2*qa))
#        if qb + qc - 2*qa < 0:
#            t = t + pi/2
#    if q2 < p2:
#        e = sqrt(1 - q2/p2)
#    else:
#        e = sqrt(1 - p2/q2)
#    mf1a = 1/2 + e*cos(t)
#    mf1b = 1/2 - e*cos(t - pi/3)
#    mf1c = 1/2 - e*cos(t + pi/3)
#    mf1 = mf1a + mf1b + mf1c
#    tri_F1 = TriBarycentric2LinearCoordinates(a,b,c,[mf1a/mf1,mf1b/mf1])
#    tri_G = TriCentroid(a,b,c)
#    tri_F2 = TriReflectionPoint(a,b,c,tri_G,tri_F1)
#    return [tri_F1,tri_F2]

# TriTriangleSteinerCircumellipseFoci returns Bickart points of ABC circumellipse when k = 2
# I removed obsolete TriTriangleSteinerCircumellipseBickartPoints(a,b,c)
def TriTriangleSteinerCircumellipseFocus(a,b,c,k):
    qa = a^2 ; qb = b^2 ; qc = c^2
    Sw = 1/2*(qa + qb + qc) ; S2 = 4*TriSquaredAreaFromQuadrances(qa,qb,qc)
    Z = sqrt(qa^2 + qb^2 + qc^2 - qa*qb - qb*qc - qc*qa)
    Q = sqrt(2*(qa*qb*qc*Z^3) - 2*(qb^2*qc^2*(2*S2 - qa*Sw) + qc^2*qa^2*(2*S2 - qb*Sw) + qa^2*qb^2*(2*S2 - qc*Sw)))
    ta = (Q - k*(qb - qc)*(qa^2 - qb*qc - qa*Z))/a
    tb = (Q - k*(qc - qa)*(qb^2 - qc*qa - qb*Z))/b
    tc = (Q - k*(qa - qb)*(qc^2 - qa*qb - qc*Z))/c
    return [ta,tb,tc]

#TriTriangleSteinerCircumellipseBickartPoints returns Bickart points (= foci) of ABC circumellipse
def TriTriangleSteinerCircumellipseBickartPoints(a,b,c):
    tri_G = TriCentroid(a,b,c)
    tri_F1 = TriTriangleSteinerCircumellipseFocus(a,b,c,2)
    tri_F2 = TriReflectionPoint(a,b,c,tri_G,tri_F1)
    return [tri_F1,tri_F2]

# TriSteinerCircumellipsePoints returns 5 points on the conic : A,B,C and centroid
# or  point [c*a + b*c,c*a + b*c, -a*b] or [a*b + b*c, -c*a,  a*b + b*c]
# or [ -b*c, b*a + a*c, b*a + a*c] and reflection of last point with respect to G

def TriIsDoublePoint(tri_1,tri_2,tri_3,tri_4):
    q14 = TriQuadranceTwoPoints(a,b,c,tri_1,tri_4)
    q24 = TriQuadranceTwoPoints(a,b,c,tri_2,tri_4)
    q34 = TriQuadranceTwoPoints(a,b,c,tri_3,tri_4)
    return  q14 == 0 or q24 == 0 or q34 == 0

def TriSteinerCircumellipsePoints(a,b,c):
    [tri_1,tri_2,tri_3] = TriMainTriangle(a,b,c) # A,B,C
    tri_G = TriCentroid(a,b,c)
    tri_4 = TriReflectionPoint(a,b,c,tri_G,tri_2)
    tri_5 = TriReflectionPoint(a,b,c,tri_G,tri_3)
    return [tri_1,tri_2,tri_3,tri_4,tri_5]

# TriSteinerCircumellipseTangential returns points A1,B1,C1 of tangential triangle at A,B,C
def TriSteinerCircumellipseTangential(a,b,c):
    return [ [-1/a,1/b,1/c], [1/a,-1/b,1/c], [1/a,1/b,-1/c] ]

# Kiepert Hyperbola has equation is bc(b^2-c^2)/x + ca(c^2-a^2)/y + ab(a^2-b^2)/z = 0 where x:y:z are trilinear coordinates
def TriIsOnKiepertHyperbola(a,b,c,t):
    [ta,tb,tc] = t
    return b*c*(b^2-c^2)/ta + c*a*(c^2-a^2)/tb + a*b*(a^2-b^2)/tc

# Kiepert Center is X115 center of Kiepert Hyperbola
def TriKiepertCenter(a,b,c):
    return [(b^2 - c^2)^2/a,(c^2 - a^2)^2/b,(a^2 - b^2)^2/c]

# TriFirstQuarterTurn returns trilinear coordinates of rectangular isosceles triangle on each edge
def TriFirstQuarterTurn(a,b,c):
    [tri_NA,tri_NB,tri_NC] = TriFirstNapoleonTriangle(a,b,c)
    kA = a*tri_NA[0] + b*tri_NA[1] + c*tri_NA[2]
    kB = a*tri_NB[0] + b*tri_NB[1] + c*tri_NB[2]
    kC = a*tri_NC[0] + b*tri_NC[1] + c*tri_NC[2]
    tri_A1 = [          2*tri_NA[0],            2*tri_NA[1]+(sqrt(3)-1)*kA/b  , 2*tri_NA[2]+(sqrt(3)-1)*kA/c ]
    tri_B1 = [ 2*tri_NB[0] +(sqrt(3)-1)*kB/a  ,     2*tri_NB[1],                2*tri_NB[2]+(sqrt(3)-1)*kB/c ]
    tri_C1 = [ 2*tri_NC[0] +(sqrt(3)-1)*kC/a  , 2*tri_NC[1]+(sqrt(3)-1)*kC/b,        2*tri_NC[2]             ]
    return [tri_A1,tri_B1,tri_C1]

# TriSecondQuarterTurn returns trilinear coordinates of rectangular isosceles triangle on each edge
def TriSecondQuarterTurn(a,b,c):
    [tri_NA,tri_NB,tri_NC] = TriSecondNapoleonTriangle(a,b,c)
    kA = a*tri_NA[0] + b*tri_NA[1] + c*tri_NA[2]
    kB = a*tri_NB[0] + b*tri_NB[1] + c*tri_NB[2]
    kC = a*tri_NC[0] + b*tri_NC[1] + c*tri_NC[2]
    tri_A1 = [          2*tri_NA[0],            2*tri_NA[1]+(sqrt(3)-1)*kA/b  , 2*tri_NA[2]+(sqrt(3)-1)*kA/c ]
    tri_B1 = [ 2*tri_NB[0] +(sqrt(3)-1)*kB/a  ,     2*tri_NB[1],                2*tri_NB[2]+(sqrt(3)-1)*kB/c ]
    tri_C1 = [ 2*tri_NC[0] +(sqrt(3)-1)*kC/a  , 2*tri_NC[1]+(sqrt(3)-1)*kC/b,        2*tri_NC[2]             ]
    return [tri_A1,tri_B1,tri_C1]

# TriKiepertTriangle returns Kiepert triangle got when erecting similar isoscele triangles on edges
# g is tangent of half angle phi, where phi is the duplicated (and base) angle for the isoscele triangle
# BE CAREFUL : if Phi = pi/2 then sn_Phi = 1 cs_Phi = 0 and
#   a*(-1) + b*sin(C + pi/2) + c*sin(B+pi/2) = -a + b*cos(C) + c*cos(B) = -a + b*(a^2 + b^2 - c^2)/(2*a*b) + c*(a^2 + c^2 - b^2)/(2*a*c)
#  = (1/2*a)*( - 2*a^2 + a^2 + b^2 - c^2 + a^2 + c^2 - b^2 ) = 0
#  that's why calling TriFirstQuarterTurn
def RT_TriKiepertTriangle(sn_A,cs_A,sn_B,cs_B,sn_C,cs_C,sn_Phi,cs_Phi):
    if sn_Phi != 1:
        sn_A_plus_Phi = sn_A*cs_Phi + cs_A*sn_Phi
        sn_B_plus_Phi = sn_B*cs_Phi + cs_B*sn_Phi
        sn_C_plus_Phi = sn_C*cs_Phi + cs_C*sn_Phi
        tri_A1 = [   -sn_Phi,     sn_C_plus_Phi, sn_B_plus_Phi ]
        tri_B1 = [ sn_C_plus_Phi,    -sn_Phi,    sn_A_plus_Phi ]
        tri_C1 = [ sn_B_plus_Phi, sn_A_plus_Phi,    -sn_Phi    ]
    else:
        #print "specific point"
        [tri_A1,tri_B1,tri_C1] = TriFirstQuarterTurn(a,b,c)
    return [tri_A1,tri_B1,tri_C1]

def TriKiepertTriangle(a,b,c,g):
    sn_Phi = 2*g/(1+g^2) ; cs_Phi = (1-g^2)/(1+g^2)
    [sn_A,sn_B,sn_C] = TriTriangleSines(a,b,c)
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    return RT_TriKiepertTriangle(sn_A,cs_A,sn_B,cs_B,sn_C,cs_C,sn_Phi,cs_Phi)

# Kiepert Hyperbola point is the point where lines between Kiepert Triangle vertex and ABC vertex concour
def RT_TriKiepertHyperbolaPoint(sn_A,cs_A,sn_B,cs_B,sn_C,cs_C,sn_Phi,cs_Phi):
    sn_A_plus_Phi = sn_A*cs_Phi + cs_A*sn_Phi
    sn_B_plus_Phi = sn_B*cs_Phi + cs_B*sn_Phi
    sn_C_plus_Phi = sn_C*cs_Phi + cs_C*sn_Phi
    return [sn_B_plus_Phi*sn_C_plus_Phi,sn_C_plus_Phi*sn_A_plus_Phi,sn_A_plus_Phi*sn_B_plus_Phi]

def TriKiepertHyperbolaPoint(a,b,c,g):
    sn_Phi = 2*g/(1+g^2) ; cs_Phi = (1-g^2)/(1+g^2)
    [sn_A,sn_B,sn_C] = TriTriangleSines(a,b,c)
    [cs_A,cs_B,cs_C] = TriCosines(a,b,c)
    return RT_TriKiepertHyperbolaPoint(sn_A,cs_A,sn_B,cs_B,sn_C,cs_C,sn_Phi,cs_Phi)

# Generalized Euler line for quadrilateral : line thru midpoints of diagonals, returned as two points
def TriQuadEulerLine(a,b,c,t1,t2,t3,t4):
    t5 = TriMidpoint(a,b,c,t1,t3)
    t6 = TriMidpoint(a,b,c,t2,t4)
    return [t5,t6]

# TriCyclicQuadIncenter return incenter of the cyclic quad ABDC : intersection of Euler line and A-bisector
# point D has trilinear coordinates [-1,x,c*x/(a*x - b)]
def TriCyclicQuadIncenter(a,b,c,x):
    t1 =  [2*a*b*c*x^2 - a^2*b*x - a^2*c*x - 2*b^2*c*x - b*c^2*x + c^3*x + a*b^2 + a*b*c,
              a*(2*a*b*x^2 - a^2*x - 2*b^2*x + c^2*x + a*b),
              a*(2*a*b*x^2 - a^2*x - 2*b^2*x + c^2*x + a*b)]
    return t1

def TriCyclicQuadLastVertice(a,b,c,x):
    t1 = [-1,x,c*x/(a*x - b)]
    return t1

# TriOrthogonalQuadLastVertice returns trilinear coordinates of D with a = BC, b = AC, c = AB, d = BD, e = CD
# If h is altitude of D on BC : h^2 + l^2 = d^2 and h^2 + (a - l)^2 = h^2 + l^2 + a^2 - 2*a*l = e^2
# => 2*a*l = a^2 + d^2 - e^2 => l = (a^2 + d^2 - e^2)/(2*a) and h^2 = d^2 - l^2
# OLD :
#def TriOrthogonalQuadLastVertice(a,b,c,S,d,e):
#    qa = a^2 ; qd = d^2 ; qe = e^2
#    ql = (qa + qd - qe)^2/(4*qa) ; l = sqrt(ql)
#    qh = qd - ql ; h = sqrt(qh)
#    ct = l/h ; cu = (a-l)/h
#    t1 = TriConwayFormula_A(a,b,c,S,ct,cu)
#    return t1
def TriOrthogonalQuadLastVertice(a,b,c,S,d,e):
    S1 = TriArea(a,d,e) # Area of triangle BDC complementary with ABC
    t1 = [-2*a*b*c,  c*( (a^2 + b^2 - c^2) + S/S1*(a^2 - d^2 + e^2) ), b*( (a^2 - b^2 + c^2) + S/S1*(a^2 + d^2 - e^2) )]
    return t1

# TriCircularSectorCentroid returns centroid of circle sector radius c, angle 2*alpha
# reference triangle is ABC, B is circle center and BC is angle bisector
# BL = 4*c*sn_alpha/(6*alpha) = 2/3*sn_alpha/alpha
def TriCircularSectorCentroid(a,b,c,alpha,sn_alpha):
    t1 = [0, 3*a*alpha - 2*c*sn_alpha, 2*b*sn_alpha]
    return t1

# Yff central triangle functions

# Trilinear coordinates or various points
def Yff_Coordinates(a,b,c):
    # Yff central triangle coordinates
    tri_1 = [1,0,0]
    tri_2 = [0,1,0]
    tri_3 = [0,0,1]
    # Yff main triangle
    tri_4 = [ -(a + b - c)/c , (a + b + c)*b/((b + c)*c) , (a + b)*(a - b - c)/((b + c)*c) ]
    tri_5 = [ -(a + b + c)*a*b/((a + c)^2*c) , (a + b - c)*b/((a + c)*c) , (a + b)*(a - b + c)*b/((a + c)^2*c) ]
    tri_6 = [ (a - b + c)*b/(a + c)^2 , -(a - b - c)*b/((a + c)*(b + c)) , -(a + b + c)*b*c/((a + c)^2*(b + c)) ]
    # Another triangle thru intersection points
    tri_7 = [ -(a + b + c)*(b + c)/((a + c)*c) , (a + b - c)/c , (a + b)*(a - b + c)/((a + c)*c) ]
    tri_8 = [ -(a + b - c)*b/((a + c)*c) , (a + b + c)*b/((b + c)*c) , (a + b)*(a - b - c)*b/((a + c)*(b + c)*c) ]
    tri_9 = [ -(a + b)*(a - b + c)*(b + c)/((a + c)*a*c) , (a + b)*(a - b - c)/(a*c) , (a + b + c)*(a + b)^2/((a + c)*a*c) ]
    # Yff intersections of edges lines of ABC triangle with main triangle UVW edges
    tri_23_45_79 = [0,-1,(b+a)/c]
    tri_23_56_78 = [0,-1,b/(c+a)]
    tri_13_45_89 = [-1,0,(b+a)/c]
    tri_13_46_78 = [-1,0,a/(c+b)]
    tri_12_56_89 = [b/(a+c),-1,0]
    tri_12_46_79 = [(b+c)/a,-1,0]
    #
    return [tri_1,tri_2,tri_3,tri_4,tri_5,tri_6,tri_7,tri_8,tri_9,tri_23_45_79,tri_23_56_78,tri_13_45_89,tri_13_46_78,tri_12_56_89,tri_12_46_79]

# TriYffCenter returns Yff center of congruence
def TriYffCenter(a,b,c):
    t_Yc = [sqrt(b*c/(-a+b+c)),sqrt(a*c/(a-b+c)),sqrt(a*b/(a+b-c))]
    return t_Yc

# Squared Area of main triangle
def Yff_MainTriangleQuadArea(a,b,c):
    QuadArea = 78767921/339061792*(a^3 - a^2*b - a*b^2 + b^3 - a^2*c - 2*a*b*c - b^2*c - a*c^2 - b*c^2 + c^3)^4*(a + b + c)^2/((a^2 + b^2 - c^2)^2*(a^2 - b^2 + c^2)^2*(a^2 - b^2 - c^2)^2)
    return QuadArea

# Area of main triangle
def Yff_MainTriangleArea(a,b,c):
    Area = sqrt(78767921/339061792)*(a^3 - a^2*b - a*b^2 + b^3 - a^2*c - 2*a*b*c - b^2*c - a*c^2 - b*c^2 + c^3)^2*(a + b + c)/((a^2 + b^2 - c^2)*(a^2 - b^2 + c^2)*(a^2 - b^2 - c^2))
    return Area

# Quadrances of edges lengths main triangle
def Yff_MainTriangleQuadrances(a,b,c):
    q45 =  1/84765448*(21191362*a^4 - 42382724*a^2*b^2 + 21191362*b^4 - 42382724*a^2*c^2 + 42382724*b^2*c^2 + 21191362*c^4 + 78767921*a^2 + 157535842*a*b + 78767921*b^2 + 157535842*a*c + 157535842*b*c + 78767921*c^2)*(a^3 - a^2*b - a*b^2 + b^3 - a^2*c - 2*a*b*c - b^2*c - a*c^2 - b*c^2 + c^3)^2*c^2/((a^2 - b^2 + c^2)^2*(a^2 - b^2 - c^2)^2*b^2)
    q46 =  1/84765448*(21191362*a^4*b^4 - 42382724*a^2*b^6 + 21191362*b^8 + 42382724*a^4*b^2*c^2 + 42382724*a^2*b^4*c^2 - 84765448*b^6*c^2 + 21191362*a^4*c^4 + 42382724*a^2*b^2*c^4 + 127148172*b^4*c^4 - 42382724*a^2*c^6 - 84765448*b^2*c^6 + 21191362*c^8 + 78767921*a^2*b^4 + 157535842*a*b^5 + 78767921*b^6 + 157535842*a*b^4*c + 157535842*b^5*c - 157535842*a^2*b^2*c^2 - 315071684*a*b^3*c^2 - 78767921*b^4*c^2 - 315071684*a*b^2*c^3 - 315071684*b^3*c^3 + 78767921*a^2*c^4 + 157535842*a*b*c^4 - 78767921*b^2*c^4 + 157535842*a*c^5 + 157535842*b*c^5 + 78767921*c^6)*(a^3 - a^2*b - a*b^2 + b^3 - a^2*c - 2*a*b*c - b^2*c - a*c^2 - b*c^2 + c^3)^2/((a^2 + b^2 - c^2)^2*(a^2 - b^2 + c^2)^2*b^2*c^2)
    q56 =  1/84765448*(21191362*a^4 - 42382724*a^2*b^2 + 21191362*b^4 - 42382724*a^2*c^2 + 42382724*b^2*c^2 + 21191362*c^4 + 78767921*a^2 + 157535842*a*b + 78767921*b^2 + 157535842*a*c + 157535842*b*c + 78767921*c^2)*(a^3 - a^2*b - a*b^2 + b^3 - a^2*c - 2*a*b*c - b^2*c - a*c^2 - b*c^2 + c^3)^2*b^2/((a^2 + b^2 - c^2)^2*(a^2 - b^2 - c^2)^2*c^2)
    return [q45,q46,q56]

# Squared Area of central triangle
def Yff_CentralTriangleQuadArea(a,b,c):
    QuadArea = (a+b+c)*(a+b-c)*(a-b+c)*(-a+b+c)/16
    return QuadArea

# Area of triangle p1-p2-p3
def Yff_CentralTriangleArea(a,b,c):
    QuadArea = Yff_CentralTriangleQuadArea(a,b,c)
    Area = QuadArea
    return Area

# Quadrances of edges lengths main triangle
def Yff_CentralTriangleQuadrances(a,b,c):
    q23 = a^2
    q13 = b^2
    q12 = c^2
    return [q12,q13,q23]

# Haberdasher functions

# Japanese theorem for cyclic quadrilateral
#    the theorem proves that the four incenters or the four triangles of a triangulation make a rectangle

# TriJapaneseCyclicQuadrilateral_Vertices returns vertices for a quadrilateral ABDC made of two triangles,
# one is ABC the reference triangle and BDC is the other one ; parameter t giving trilinear coordinates for D
def TriJapaneseCyclicQuadrilateral_Vertices(a,b,c,t):
    # Set vertices of quadrilateral ABCD and incenter
    [tri_1,tri_2,tri_3] = TriMainTriangle(a,b,c)  # A = tri_1 ; B = tri_2 ; C = tri_3
    tri_18 = TriIncenter(a,b,c) ; tri_4 = t   # I = tri_18 ; D = tri_4
    # Compute point intersection of diagonals
    l1_4 = TriLine(a,b,c,tri_1,tri_4) ; l2_3 = TriLine(a,b,c,tri_2,tri_3)
    tri_5 = TriLineIntersect(a,b,c,l1_4,l2_3) ; x5 = (tri_5[1]/tri_5[2])/(c/b) # E = tri_5
    #
     # Angle bisector theorem for triangle BDE
    k2 = TriQuadranceTwoPoints(a,b,c,tri_4,tri_2) / TriQuadranceTwoPoints(a,b,c,tri_4,tri_5) ; k = sqrt(k2)
    x6 = (1/k)*((k + 1)*x5 + 1) ; tri_6 = [0,c*x6,b] # F = tri_6
    #  Angle bisector theorem for triangle EAB
    k2 = TriQuadranceTwoPoints(a,b,c,tri_1,tri_2) / TriQuadranceTwoPoints(a,b,c,tri_1,tri_5) ; k = sqrt(k2)
    x7 = (1/k)*((k + 1)*x5 + 1) ; tri_7 = [0,c*x7,b] # G = tri_7
    # Intersect DF and AG to get incenter H of ABD triangle
    l4_6 = TriLine(a,b,c,tri_4,tri_6) ; l1_7 = TriLine(a,b,c,tri_1,tri_7)
    tri_8 = TriLineIntersect(a,b,c,l4_6,l1_7) # H = tri_8
    #
    # Angle bisector theorem for triangle CAE
    k2 = TriQuadranceTwoPoints(a,b,c,tri_1,tri_3) / TriQuadranceTwoPoints(a,b,c,tri_1,tri_5) ; k = sqrt(k2)
    x9 = k*x5/(k + x5 + 1) ; tri_9 = [0,c*x9,b] # M = tri_9
    # Angle bisector theorem for triangle EDC
    k2 = TriQuadranceTwoPoints(a,b,c,tri_4,tri_3) / TriQuadranceTwoPoints(a,b,c,tri_4,tri_5) ; k = sqrt(k2)
    x10 = k*x5/(k + x5 + 1) ; tri_10 = [0,c*x10,b] # N = tri_10
    # Intersect AM and DN to get incenter of ACD triangle
    l1_9 = TriLine(a,b,c,tri_1,tri_9) ; l4_10 = TriLine(a,b,c,tri_4,tri_10)
    tri_11 = TriLineIntersect(a,b,c,l1_9,l4_10)  # O = tri_11
    #
    # Compute middle K of OH and then P as reflection of I
    tri_12 = TriMidpoint(a,b,c,tri_8,tri_11) # K = tri_12
    tri_13 = TriReflectionPoint(a,b,c,tri_12,tri_18) # P = tri_13
    #
    # Quadrilaterals ABHI, BDPH, DCOP, CAIO are cyclic..finding circumcenters
    tri_14 = TriTriangleCircumcenter(a,b,c,tri_1,tri_2,tri_18) # Q = tri_14
    tri_15 = TriTriangleCircumcenter(a,b,c,tri_2,tri_4,tri_8)  # R = tri_15
    tri_16 = TriTriangleCircumcenter(a,b,c,tri_4,tri_3,tri_13) # T = tri_16
    tri_17 = TriTriangleCircumcenter(a,b,c,tri_3,tri_1,tri_11) # U = tri_17
    #
    vertices = [tri_1,tri_2,tri_3,tri_4,tri_5,tri_6,tri_7,tri_8,tri_9,tri_10,tri_11,tri_12,tri_13,tri_14,tri_15,tri_16,tri_17,tri_18]
    return vertices

# TriJapaneseCyclicQuadrilateral_Incenters returns only incenters, of cyclic quadrilateral ABC + extra vertice D given by t

def TriJapaneseCyclicQuadrilateral_Incenters(a,b,c,t):
    [tri_1,tri_2,tri_3,tri_4,tri_5,tri_6,tri_7,tri_8,tri_9,tri_10,tri_11,tri_12,tri_13,tri_14,tri_15,tri_16,tri_17,tri_18] = TriJapaneseCyclicQuadrilateral_Vertices(a,b,c,t)
    return [tri_18,tri_8,tri_11,tri_13] # I,H,O,P

# TriJapaneseCyclicQuadrilateral_Vertices returns edges for a drawing of the five cyclic quadriletrals icluding one rectangle
def TriJapaneseCyclicQuadrilateral_Edges():
    edge_1  =  [ 1, 2   ]    # AB
    edge_2  =  [ 2, 4   ]    # BD
    edge_3  =  [ 4, 3   ]    # DC
    edge_4  =  [ 3, 1   ]    # CA
    edge_5  =  [ 18, 8  ]    # IH
    edge_6  =  [ 8, 13  ]    # HP
    edge_7  =  [ 13, 11 ]    # PO
    edge_8  =  [ 11, 18 ]    # OI
    edge_9 =   [ 1, 18  ]    # AI
    edge_10 =  [ 2, 8   ]    # BH
    edge_11 =  [ 3, 11  ]    # CO
    edge_12 =  [ 4, 13  ]    # DP
    edges = [edge_1,edge_2,edge_3,edge_4,edge_5,edge_6,edge_7,edge_8,edge_9,edge_10,edge_11,edge_12]
    #
    return edges


# TriHaberdasher_A_Outside_LimitPointsSegment returns true when segment of limits points doesn't exist or
# A is outside that segment
def TriHaberdasher_A_Outside_LimitPointsSegment(a,b,c,S):
    [ctA,ctB,ctC] = TriConwayTriangleCotangents(a,b,c,S)
    outside = False
    if ctB == 0:
        outside = True
    else:
        u = ctA/ctB
        outside = ( (2*u + 1) > 0 )
    return outside

# TriHaberdasher_Parameters returns parameters for Haberdasher building
def TriHaberdasher_Parameters(a,b,c,y):
    r1   = a^2 - b^2 + c^2 - (a^2 - b^2 - 3*c^2)*y
    r2   =  - 3*a^2 - b^2 + c^2  - (a^2 - 5*b^2 + c^2)*y
    r3   = 3 + y
    r4   = (a^2 - 3*b^2 - c^2)*y + a^2 + b^2 - c^2
    r5   = - 1 + y
    r6   = - 5*a^2 - b^2 + c^2  - 4*(a^2 - 2*b^2 + c^2)*y + (a^2 + b^2 - 5*c^2)*y^2
    r7   = a^2 - b^2 - 3*c^2
    r8   = a^2 - 5*b^2 + c^2
    r9   = a^2 - 3*b^2 - c^2
    r10  = a^2 + b^2 - c^2
    r11  = a^2 - b^2
    return [r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11]

# TriHaberdasher_Coordinates returns trilinear coordinates or various points
def TriHaberdasher_Coordinates(a,b,c,y):
    # Get parameters
    [r1,r2,r3,r4,r5,r6,r7,r8,r9,r10,r11] = TriHaberdasher_Parameters(a,b,c,y)
    # Reference triangle
    tri_1 = [1,0,0]  # A
    tri_2 = [0,1,0]  # B
    tri_3 = [0,0,1]  # C
    # Midpoints of edges
    tri_4 = TriMidpoint(a,b,c,tri_1,tri_2)  # D = (A+B)/2
    tri_5 = TriMidpoint(a,b,c,tri_2,tri_3)  # E = (B+C)/2
    # Set cutting point on edge
    tri_6 = [c*y,0,a]  # F
    # Compute first inner cutting point
    l5_6 = TriLine(a,b,c,tri_5,tri_6) # line EF
    l4_7 = [-a*r1,b*r1,c*r2]  # line DG
    tri_7 = TriLineIntersect(a,b,c,l5_6,l4_7) # G
    # Compute second inner cutting point
    tri_8 = TriReflectionPoint(a,b,c,tri_6,tri_1)  # H
    tri_9 = TriMidpoint(a,b,c,tri_8,tri_3)   # J = (C+H)/2
    l9_10 = [-a*r3*r4,b*r6,c*r4*r5] # line JK
    tri_10 = TriLineIntersect(a,b,c,l5_6,l9_10) # K
    # Compute points after rotation
    tri_11 = TriReflectionPoint(a,b,c,tri_6,tri_7) # L
    tri_12 = TriReflectionPoint(a,b,c,tri_6,tri_4) # M
    tri_13 = TriReflectionPoint(a,b,c,tri_9,tri_5) # N
    tri_14 = TriReflectionPoint(a,b,c,tri_9,tri_10) # O
    l11_12 = TriLine(a,b,c,tri_11,tri_12) # line LM
    l13_14 = TriLine(a,b,c,tri_13,tri_14) # line NO
    tri_15 = TriLineIntersect(a,b,c,l11_12,l13_14) # P
    # Compute limiting points on line cuts
    l4_16 = [-a*r7, b*r7, c*r8] # line DQ
    l1_5 = TriLine(a,b,c,tri_1,tri_5) # line AE
    tri_16 = TriLineIntersect(a,b,c,l4_16,l1_5) # Q
    tri_17 = TriMidpoint(a,b,c,tri_1,tri_3) # R = (A+C)/2
    l5_17 = TriLine(a,b,c,tri_5,tri_17) # line ER
    l4_18 = [-a*c,b*c,-r11] # line DT
    tri_18 = TriLineIntersect(a,b,c,l5_17,l4_18) # T
    # Compute projection of E on AC
    l1_3 = TriLine(a,b,c,tri_1,tri_3) # line AC
    l5_19 = [-a*r9, b*r10,-c*r10]
    tri_19 = TriLineIntersect(a,b,c,l1_3,l5_19)
    return [tri_1,tri_2,tri_3,tri_4,tri_5,tri_6,tri_7,tri_8,tri_9,tri_10,tri_11,tri_12,tri_13,tri_14,tri_15,tri_16,tri_17,tri_18,tri_19]

def TriHaberdasher_Limit_Points_Exist(a,b,c,S):
    [ctA,ctB,ctC] = TriConwayTriangleCotangents(a,b,c,S)
    exist = False
    if ctB == 0:
        exist = (ctA  != -1)
    else:
        u = ctA/ctB ; v = ctC/ctB
        #exist = ((u^2 - 2*u*v + v^2 - 4*u - 4*v)*(u + 1)*(v + 1) >= 0)
        if ((u == -1/2) and (v == 3/2)):
            exist = False
        else:
            exist = ((u^2 - 2*u*v + v^2 - 4*u - 4*v) > 0)
    return exist

# TriHaberdasher_Limit_Points returns points for roots of equation
# (2*ctA + ctB)*y^2 - 2*(ctA - ctB + ctC)*y + (ctB + 2*ctC) = 0
# or (2*u + 1)*y^2 - 2*(u + v - 1)*y + (2*v + 1) = 0
def TriHaberdasher_Limit_Points_Parameters(a,b,c,S):
    [ctA,ctB,ctC] = TriConwayTriangleCotangents(a,b,c,S)
    if ctB == 0:
        # equation is (ctA )*y^2 - (ctA + ctC)*y + (ctC)  = (y - 1)*(ctA*y - ctC) = 0
        # or because ctA ctC = 1 : (y - 1)*(ctA*y - 1/ctA) = 0
        if ctA < ctC :
            y1 = 1/ctA^2 ; y2 = 1
        else:
            y1 = 1 ; y2 = 1/ctA^2
    else:
        u = ctA/ctB ; v = ctC/ctB
        if u == -1/2:
            y1 = y2 = (1/2)*(2*v + 1)/(v - 3/2)  # v different from 3/2
        else:
            disc = (u - v)^2 - 4*(u + v) # u^2 - 2*u*v + v^2 - 4*u - 4*v
            m = u + v - 1 ; n = 2*u + 1
            y1 = (m + sqrt(disc))/n
            y2 = (m - sqrt(disc))/n
            if ctB < 0:
                y0 = y1 ; y1 = y2 ; y2 = y0 # swap y1,y2
    return [y1,y2]

def TriHaberdasher_Limit_Points(a,b,c,S):
    [y1,y2] = TriHaberdasher_Limit_Points_Parameters(a,b,c,S)
    return [[c*y1,0,a],[c*y2,0,a]]

# TriHaberdasher_MakeRectangle_Exist checks if
def TriHaberdasher_MakeRectangle_Exist(a,b,c,S,k):
    [ctA,ctB,ctC] = TriConwayTriangleCotangents(a,b,c,S)
    if ctB == 0:
        u = ctA
        if 2*k*u - 4*u^2 - 1 == 0:
            exist = False
        else:
            exist = ( 2*k*ctA^2 - ctA + 2*k >= 0 )
    else:
        u = ctA/ctB ; v = ctC/ctB ; w = sqrt(u + v + u*v)
        disc =  8*k*(u + v)*w  - 4*(u + v + u*v); print "disc = ",disc
        exist = ( disc >= 0 )
    return exist

# TriHaberdasher_MakeRectangle returns cutting point in order to make rectangle with width/height = k
# from equation (2*k*w - 4*u - 1  - v)*y^2 + 2*(2*k*w - 1 + v)*y + 2*k*w - 1 - v = 0  when ctB is not zero
# or (4*ctA^2 - 2*k*ctA + 1)*y^2 - 2*(2*k*ctA + 1)*y - 2*k*ctA + 1  = 0 when ctB = 0
def TriHaberdasher_MakeRectangle_Parameters(a,b,c,S,k):
    [ctA,ctB,ctC] = TriConwayTriangleCotangents(a,b,c,S)
    if ctB == 0:
        u = ctA
        n = 4*u^2 - 2*k*u + 1
        if n == 0:
            y = - (2*u^2)/(4*u^2 + 2)
        else:
            disc = (2*k*u + 1)^2 - (4*u^2 - 2*k*u + 1)*(- 2*k*u + 1)
            m = 2*k*u + 1
            y1 = (m - sqrt(disc))/n
            y2 = (m + sqrt(disc))/n
    else:
        u = ctA/ctB ; v = ctC/ctB ; w = sqrt(u + v + u*v)  # w = 1/ctB^2 because cotA cotB + cotA cotC + cotB cotC = 1
        n = 2*k*w - 4*u - 1 - v
        if n == 0:
            if 2*k*w - 1 + v == 0:
                y = 0.0 # choose any y, because v = -2*u, and any y is root
            else:
                y = -(1/2)*(2*k*w - 1 - v)/(2*k*w - 1 + v)
        else:
            disc =  8*k*(u + v)*w  - 4*(u + v + u*v)  # w^2 has been replaced with (u + v + u*v)
            m = - (2*k*w - 1 + v)
            y1 = (m - sqrt(disc))/n
            y2 = (m + sqrt(disc))/n
    return [y1,y2]

def TriHaberdasher_MakeRectangle_Points(a,b,c,S):
    [y1,y2] = TriHaberdasher_MakeRectangle_Parameters(a,b,c,S,k)
    return [[c*y1,0,a],[c*y2,0,a]]

# TriHaberdasher_Bisector returns y such as F is on the bisector line of AB and angle B bisector (BI)
# we use angle bisector theorem where X intersection of BI with AC, I incenter : AY/XY = AB/BX
# returned value y is given too by bisector theorem but for triangle ABC and cevian BY
def TriHaberdasher_Bisector(a,b,c):
    dBX = sqrt((a + b + c)*(a - b + c)*a*c)/(a + c)
    dAX = b*c/(a + c)
    dCX =  a*b/(a + c)
    dAY = c*dAX/(c + dBX)
    dXY = dAX*dBX/(c + dBX)
    dCY = dCX + dXY
    y = dCY/dAY
    return y

# TriPinwheelSquare1 is compass-and-ruler computation, splitting square FGUB in 25 tiles
# 23 tiles are triangles and 2 tiles are cyclic quadrilaterals
def TriPinwheelSquare1(a,b,c):
    tri_1 = [1,0,0]                                   # F
    tri_2 = [0,1,0]                                   # G
    tri_3 = [0,0,1]                                   # U
    tri_4 = TriMidpoint(a,b,c,tri_2,tri_3)            # N = (G+U)/2
    tri_5 = TriReflectionPoint(a,b,c,tri_4,tri_1)     # B = 2D - A
    tri_6 = TriMidpoint(a,b,c,tri_1,tri_3)            # H = (A+C)/2
    tri_7 = TriMidpoint(a,b,c,tri_5,tri_3)            # P = (B+C)/2
    l5_6 = TriLine(a,b,c,tri_5,tri_6)                 # BH
    l2_7 = TriLine(a,b,c,tri_2,tri_7)                 # GP
    tri_8 = TriLineIntersect(a,b,c,l5_6,l2_7)         # A = BH n GP
    tri_9 = TriMidpoint(a,b,c,tri_5,tri_2)            # D = (B+G)/2
    l9_1 = TriLine(a,b,c,tri_9,tri_1)                 # DF
    l2_6 = TriLine(a,b,c,tri_2,tri_6)                 # GH
    tri_10 = TriLineIntersect(a,b,c,l9_1,l2_6)        # W = DF n GH
    tri_11 = TriMidpoint(a,b,c,tri_8,tri_1)           # J = (A+F)/2
    tri_12 = TriMidpoint(a,b,c,tri_6,tri_3)           # X = (H+U)/2
    tri_13 = TriMidpoint(a,b,c,tri_7,tri_3)           # C = (P+U)/2
    tri_14 = TriMidpoint(a,b,c,tri_7,tri_12)          # R = (P+X)/2
    tri_15 = TriMidpoint(a,b,c,tri_1,tri_6)           # K = (F+H)/2
    tri_16 = TriMidpoint(a,b,c,tri_1,tri_2)           # O = (F+G)/2
    tri_17 = TriMidpoint(a,b,c,tri_1,tri_16)          # Q = (F+O)/2
    tri_18 = TriMidpoint(a,b,c,tri_9,tri_2)           # V = (D+G)/2
    tri_19 = TriMidpoint(a,b,c,tri_18,tri_10)         # Z=S = (V+W)/2
    tri_20 = TriMidpoint(a,b,c,tri_9,tri_18)          # Y = (D+V)/2
    tri_21 = TriMidpoint(a,b,c,tri_19,tri_18)         # T = (Z+V)/2
    tri_22 = TriMidpoint(a,b,c,tri_2,tri_18)          # L = (G+V)/2
    tri_23 = TriMidpoint(a,b,c,tri_5,tri_9)           # E = (B+D)/2
    tri_24 = TriMidpoint(a,b,c,tri_5,tri_23)          # I = (B+E)/2
    tri_25 = TriMidpoint(a,b,c,tri_9,tri_23)          # M = (D+E)/2
    return [tri_1,tri_2,tri_3,tri_4,tri_5,tri_6,tri_7,tri_8,tri_9,tri_10,tri_11,tri_12,tri_13,tri_14,tri_15,tri_16,tri_17,tri_18,tri_19,tri_20,tri_21,tri_22,tri_23,tri_24,tri_25]

# TriPinwheelSquare_Vertices is same than TriPinwheelSquare1, returning barycentric coordinates with respect to reference triangle FGU
def TriPinwheelSquare_Vertices():
    sq2 = sqrt(2)
    tri_1  = [1, 0, 0]                # F
    tri_2  = [0, 1, 0]                # G
    tri_3  = [0, 0, 1]                # U
    tri_4  = [0, 1, 1]                # N
    tri_5  = [-1, sq2, sq2]           # B
    tri_6  = [1, 0, sq2]              # H
    tri_7  = [-1, sq2, 2*sq2]         # P
    tri_8  = [-sq2, 3, 4]             # A
    tri_9  = [-1, 2*sq2, sq2]         # D
    tri_10 = [sq2, 4, 2]              # W
    tri_11 = [3, 3*sq2, 4*sq2]        # J
    tri_12 = [1, 0, 3*sq2]            # X
    tri_13 = [-1, sq2, 4*sq2]         # C
    tri_14 = [-sq2, 4, 14]            # R
    tri_15 = [3*sq2, 0, 2]            # K
    tri_16 = [1, sq2, 0]              # O
    tri_17 = [3*sq2, 2, 0]            # Q
    tri_18 = [-1, 4*sq2, sq2]         # V
    tri_19 = [0, 3, 1]                # Z
    tri_20 = [-3*sq2, 16, 6]          # Y
    tri_21 = [-sq2, 14, 4]            # T
    tri_22 = [-1, 8*sq2, sq2]         # L
    tri_23 = [-3*sq2, 8, 6]           # E
    tri_24 = [-7, 8*sq2, 7*sq2]       # I
    tri_25 = [-5, 8*sq2, 5*sq2]       # M
    vertices = [tri_1,tri_2,tri_3,tri_4,tri_5,tri_6,tri_7,tri_8,tri_9,tri_10,tri_11,tri_12,tri_13,tri_14,tri_15,tri_16,tri_17,tri_18,tri_19,tri_20,tri_21,tri_22,tri_23,tri_24,tri_25]
    return vertices

# TriPinwheelSquare_Edges returns edges, straight line segments betweens vertices, for tiling of square FGUB
def TriPinwheelSquare_Edges():
    edge_1  =  [ 3, 1 ]    # UF
    edge_2  =  [ 1, 2 ]    # FG
    edge_3  =  [ 2, 5 ]    # GB
    edge_4  =  [ 5, 3 ]    # BU
    edge_5  =  [ 5, 6 ]    # BH
    edge_6  =  [ 7, 8 ]    # PA
    edge_7  =  [ 8, 1 ]    # AF
    edge_8  =  [ 8, 9 ]    # AD
    edge_9  =  [ 9, 10 ]   # DW
    edge_10 =  [ 10, 8 ]   # WA
    edge_11 =  [ 10, 11 ]  # WJ
    edge_12 =  [ 7, 12 ]   # PX
    edge_13 =  [ 12, 8 ]   # XA
    edge_14 =  [ 13, 14 ]  # CR
    edge_15 =  [ 14, 8 ]   # RA
    edge_16 =  [ 8, 15 ]   # AK
    edge_17 =  [ 15, 11 ]  # KJ
    edge_18 =  [ 17, 11 ]  # QJ
    edge_19 =  [ 11, 16 ]  # JO
    edge_20 =  [ 16, 10 ]  # OW
    edge_21 =  [ 16, 18 ]  # OV
    edge_22 =  [ 18, 10 ]  # VW
    edge_23 =  [ 9, 19 ]   # DZ
    edge_24 =  [ 19, 16 ]  # ZO
    edge_25 =  [ 20, 21 ]  # YT
    edge_26 =  [ 21, 16 ]  # TO
    edge_27 =  [ 22, 16 ]  # LO
    edge_28 =  [ 8, 23 ]   # AE
    edge_29 =  [ 8, 24 ]   # AM
    edges = [edge_1,edge_2,edge_3,edge_4,edge_5,edge_6,edge_7,edge_8,edge_9,edge_10,edge_11,edge_12,edge_13,edge_14,edge_15,edge_16,edge_17,edge_18,edge_19,edge_20,edge_21,edge_22,edge_23,edge_24,edge_25,edge_26,edge_27,edge_28,edge_29]
    #
    return edges

# TriPinwheelTriangle1 is compass-and-ruler computation, splitting Pinwheel triangle ABC into 25 tiles
# 23 tiles are triangles and 2 tiles are cyclic quadrilaterals
def TriPinwheelTriangle1(a,b,c):
    tri_1  = [1,0,0]   # B
    tri_2  = [0,1,0]   # A
    tri_3  = [0,0,1]   # C
    tri_4  = TriMidpoint(a,b,c,tri_2,tri_1) # D = (A+B)/2
    tri_5  = TriMidpoint(a,b,c,tri_2,tri_3) # F = (A+C)/2
    tri_6  = TriMidpoint(a,b,c,tri_3,tri_5) # I = (C+F)/2
    tri_7  = TriMidpoint(a,b,c,tri_4,tri_5) # P = (D+F)/2
    tri_8  = TriMidpoint(a,b,c,tri_2,tri_5) # J = (A+F)/2
    tri_9  = TriMidpoint(a,b,c,tri_4,tri_6) # W = (D+I)/2
    tri_10 = TriMidpoint(a,b,c,tri_1,tri_4) # Y = (B+D)/2
    l2_9   = TriLine(a,b,c,tri_2,tri_9)   # AW
    l6_10  = TriLine(a,b,c,tri_6,tri_10)  # IY
    tri_11 = TriLineIntersect(a,b,c,l2_9,l6_10) # R = AW n IY
    l5_11  = TriLine(a,b,c,tri_5,tri_11)  # FR
    l1_3   = TriLine(a,b,c,tri_1,tri_3)   # BC
    tri_12 = TriLineIntersect(a,b,c,l5_11,l1_3) # E = FR n BC
    tri_13 = TriMidpoint(a,b,c,tri_1,tri_12) # G = (B+E)/2
    tri_14 = TriMidpoint(a,b,c,tri_3,tri_12) # H = (C+E)/2
    tri_15 = TriMidpoint(a,b,c,tri_5,tri_13) # O = (F+G)/2
    tri_16 = TriMidpoint(a,b,c,tri_5,tri_14) # K = (F+H)/2
    tri_17 = TriMidpoint(a,b,c,tri_5,tri_15) # Q = (F+O)/2
    tri_18 = TriMidpoint(a,b,c,tri_4,tri_13) # M = (D+G)/2
    tri_19 = TriMidpoint(a,b,c,tri_4,tri_18) # N = (D+M)/2
    tri_20 = TriMidpoint(a,b,c,tri_13,tri_18) # L = (G+M)/2
    tri_21 = TriMidpoint(a,b,c,tri_18,tri_7) # T = (M+P)/2
    tri_22 = TriMidpoint(a,b,c,tri_18,tri_21) # S = (M+T)/2
    return [tri_1,tri_2,tri_3,tri_4,tri_5,tri_6,tri_7,tri_8,tri_9,tri_10,tri_11,tri_12,tri_13,tri_14,tri_15,tri_16,tri_17,tri_18,tri_19,tri_20,tri_21,tri_22]


def TriPinwheelTriangle_Vertices():
    sq5 = sqrt(5)
    tri_1  = [ 1 ,0, 0 ]               # B
    tri_2  = [ 0, 1, 0 ]               # A
    tri_3  = [ 0, 0, 1 ]               # C
    tri_4  = [ sq5, 2, 0 ]             # D = (A+B)/2
    tri_5  = [ 0, 1, sq5 ]             # F = (A+C)/2
    tri_6  = [ 0, 1, 3*sq5 ]           # I = (C+F)/2
    tri_7  = [ 5, 4*sq5, 10 ]          # P = (D+F)/2
    tri_8  = [ 0, 3*sq5, 5 ]           # J = (A+F)/2
    tri_9  = [ 5, 3*sq5, 15 ]          # W = (D+I)/2
    tri_10 = [ 3*sq5, 2, 0 ]           # Y = (B+D)/2
    tri_11 = [ 3, sq5, 9 ]             # R = AW n IY
    tri_12 = [ 3, 0, 4]                # E = FR n BC
    tri_13 = [ 2, 0, 1 ]               # G = (B+E)/2
    tri_14 = [ 3, 0, 14]               # H = (C+E)/2
    tri_15 = [ 4, sq5, 7 ]             # O = (F+G)/2
    tri_16 = [ 3, 2*sq5, 24 ]          # K = (F+H)/2
    tri_17 = [ 4, 3*sq5, 17 ]          # Q = (F+O)/2
    tri_18 = [ 13, 2*sq5, 4 ]          # M = (D+G)/2
    tri_19 = [ 23, 6*sq5, 4 ]          # N = (D+M)/2
    tri_20 = [ 29, 2*sq5 , 12 ]        # L = (G+M)/2
    tri_21 = [ 9, 3*sq5, 7 ]           # T = (M+P)/2
    tri_22 = [ 22, 5*sq5, 11 ]         # S = (M+T)/2
    return [tri_1,tri_2,tri_3,tri_4,tri_5,tri_6,tri_7,tri_8,tri_9,tri_10,tri_11,tri_12,tri_13,tri_14,tri_15,tri_16,tri_17,tri_18,tri_19,tri_20,tri_21,tri_22]


# TriPinwheelTriangle_Edges returns edges, straight line segments betweens vertices, for tiling of square FGUB
def TriPinwheelTriangle_Edges():
    edge_1  =  [ 2, 1 ]    # AB
    edge_2  =  [ 1, 3 ]    # BC
    edge_3  =  [ 3, 2 ]    # CA
    edge_4  =  [ 4, 7 ]    # DP
    edge_5  =  [ 7, 2 ]    # PA
    edge_6  =  [ 8, 7 ]    # JP
    edge_7  =  [ 4, 13 ]   # DG
    edge_8  =  [ 13, 5 ]   # GF
    edge_9  =  [ 5, 14 ]   # FH
    edge_10 =  [ 12, 15 ]  # EO
    edge_11 =  [ 15, 4 ]   # OD
    edge_12 =  [ 7, 15 ]   # PO
    edge_13 =  [ 15, 8 ]   # OJ
    edge_14 =  [ 15, 11]   # OR
    edge_15 =  [ 11, 16 ]  # RK
    edge_16 =  [ 16, 6 ]   # KI
    edge_17 =  [ 12, 16 ]  # EK
    edge_18 =  [ 8, 17 ]   # JQ
    edge_19 =  [ 17, 11 ]  # QR
    edge_20 =  [ 11, 12 ]  # RE
    edge_21 =  [ 1, 18 ]   # BM
    edge_22 =  [ 18, 15 ]  # MO
    edge_23 =  [ 1, 19 ]   # BN
    edge_24 =  [ 1, 20 ]   # BL
    edge_25 =  [ 20, 15 ]  # LO
    edge_26 =  [ 18, 22 ]  # MS
    edge_27 =  [ 22, 21 ]  # ST
    edge_28 =  [ 21, 7 ]   # TP
    edge_29 =  [ 19, 22 ]  # NS
    edge_30 =  [ 22, 15 ]  # SO
    edge_31 =  [ 3, 16 ]   # CK
    #
    edges = [edge_1,edge_2,edge_3,edge_4,edge_5,edge_6,edge_7,edge_8,edge_9,edge_10,edge_11,edge_12,edge_13,edge_14,edge_15,edge_16,edge_17,edge_18,edge_19,edge_20,edge_21,edge_22,edge_23,edge_24,edge_25,edge_26,edge_27,edge_28,edge_29,edge_30,edge_31]
    #
    return edges

# TriNewtonLinesSegments returns Newton lines for a reference triangle ABC and a ratio k
def TriNewtonLinesSegments(a,b,c,k):
    # Newton line with respect to C is parallel with AB and goes thru tri_1 and tri_2 points
    tri_1 = [0, c, b*k]
    tri_2 = [c, 0, a*k]
    # Newton line with respect to A is parallel with BC and goes thru tri_3 and tri_4 points
    tri_3 = [c*k, 0, a]
    tri_4 = [b*k, a, 0]
    # Newton line with respect to B is parallel with CA and goes thru tri_5 and tri_6 points
    tri_5 = [b, a*k, 0]
    tri_6 = [0, c*k, b]
    return [tri_1,tri_2,tri_3,tri_4,tri_5,tri_6]

# TriNewtonLinesRatios returns Newton lines ratios for a reference triangle ABC and a point given by trilinear coordinates
def TriNewtonLinesRatios(a,b,c,S,t):
    # Convert trilinear coordinates to barycentric and then areal coordinates
    m = TriBarycentric(a,b,c,tri_L)
    u = TriBarycentric2ArealCoordinates(S,m)
    # Retrieving ratios from the area of sub-triangles
    k1 = u[0] / (u[1] + u[2])
    k2 = u[1] / (u[2] + u[0])
    k3 = u[2] / (u[0] + u[1])
    return [k1,k2,k3]


# Find all roots is general function to find all roots of equation (copy from : https://ask.sagemath.org/question/7823/numerically-find-all-roots-in-an-interval/)
def find_all_roots(f, a, b, eps=0.000001):
    roots = []
    intervals_to_check = [(a,b)]
    while intervals_to_check:
        start, end = intervals_to_check.pop()
        try:
            root = find_root(f, start, end)
        except RuntimeError:
            continue
        if root in roots:
            continue
        if abs(f(root)) < 1:
            roots.append(root)
        intervals_to_check.extend([(start, root-eps), (root+eps, end)])
    roots.sort()
    return roots

# -- next dissection functions are used for dissecting a triangle ABC into five 120 degrees triangles

# TriTriangleDissectionCubic set coefficients for the cubic  in u = a^2/S ; v = b^2/S ; w = c^2/S
# 3*w^2*x^3 + w*(3*(u - v - 2*w) - 4*sqrt(3))*x^2 + (6*v*w + 3*w^2  + 8*sqrt(3)*w - 16)*x + (2*sqrt(3)*(u - v - w) + 8 - 3*v*w) = 0

def TriTriangleDissectionCubic(u,v,w):
    sq3 = sqrt(3).n()
    k3 = 3*w^2
    k2 = w*(3*(u - v - 2*w) - 4*sq3)
    k1 = 6*v*w + 3*w^2  + 8*sq3*w - 16
    k0 = 2*sq3*(u - v - w) + 8 - 3*v*w
    return [k0,k1,k2,k3]

# TriTriangleGuessDissection returns roots x for dissection a triangle in five parts when we know qa = BC^2 (quadrance of the base edge of ABC),
#  qb = CA^2 (quadrance of one another edge ABC)  qd = BD^2 = (1-x)^2 AB^2 (quadrance of one part of the last edge)
def  TriTriangleGuessDissection(qa,qb,qd):
    u = qb/qa ; v = qd/qa
    x = var('x')
    f1(x) = (u^4 - 4*u^3  + 6*u^2  - 4*u + 1)*x^10 + ( - 9*u^4 - 6*u^3*v  + 36*u^3 + 18*u^2*v - 54*u^2 - 18*u*v + 36*u + 6*v - 9)*x^9 + (37*u^4 + 50*u^3*v + 18*u^2*v^2 - 148*u^3 - 146*u^2*v - 36*u*v^2 + 222*u^2 + 142*u*v + 18*v^2 - 148*u - 46*v + 37)*x^8 + (- 92*u^4 - 188*u^3*v - 123*u^2*v^2 - 18*u*v^3 + 368*u^3 + 536*u^2*v + 234*u*v^2 + 18*v^3  - 552*u^2 - 508*u*v - 111*v^2 + 368*u + 160*v  - 92)*x^7 + (154*u^4 + 416*u^3*v + 363*u^2*v^2 + 90*u*v^3 + 9*v^4 - 616*u^3 - 1166*u^2*v - 668*u*v^2 - 90*v^3  + 924*u^2 + 1084*u*v + 309*v^2 - 616*u - 334*v + 154)*x^6 + (- 182*u^4 - 592*u^3*v - 597*u^2*v^2 - 186*u*v^3 - 36*v^4 + 728*u^3 + 1648*u^2*v + 1090*u*v^2 + 204*v^3 - 1092*u^2 - 1520*u*v - 513*v^2 + 728*u + 464*v - 182)*x^5+ (154*u^4 + 556*u^3*v + 591*u^2*v^2 + 206*u*v^3 + 63*v^4 - 616*u^3 - 1558*u^2*v - 1108*u*v^2 - 274*v^3 + 924*u^2 + 1448*u*v + 558*v^2 - 616*u - 446*v + 154)*x^4 + (- 92*u^4 - 340*u^3*v - 357*u^2*v^2 - 138*u*v^3 - 63*v^4 + 368*u^3 + 976*u^2*v + 718*u*v^2 + 234*v^3 - 552*u^2 - 932*u*v - 405*v^2  + 368*u + 296*v - 92)*x^3 + (37*u^4 + 128*u^3*v + 129*u^2*v^2 + 62*u*v^3 + 37*v^4 - 148*u^3 - 386*u^2*v - 292*u*v^2 - 124*v^3 + 222*u^2+ 388*u*v + 189*v^2 - 148*u - 130*v + 37)*x^2+ (- 9*u^4 - 26*u^3*v - 27*u^2*v^2 - 18*u*v^3 - 10*v^4 + 36*u^3 + 86*u^2*v + 70*u*v^2 + 36*v^3 - 54*u^2 - 94*u*v - 51*v^2 + 36*u + 34*v - 9)*x + (u^4 + 2*u^3*v + 3*u^2*v^2 + 2*u*v^3 + v^4 - 4*u^3 - 8*u^2*v - 8*u*v^2 - 4*v^3 + 6*u^2 + 10*u*v + 6*v^2 - 4*u - 4*v  + 1)
    # Searching first if root exist
    try:
        x0 = find_root(f1,0.001,0.999,maxiter=100)
    except RuntimeError:
        x0 = -1
    if x0 == -1:
        return []
    # If found one, searching for all of them
    roots1 = find_all_roots(f1, 0.001, 0.999)  # here we can have extra roots because f1=0 equation coming by squaring area S
    # select only roots satisfying cubic equation
    roots = []
    for x0 in roots1:
        qc = qd/(1-x0)^2 ; S2 = TriSquaredAreaFromQuadrances(qa,qb,qc); S = sqrt(S2)
        [k0,k1,k2,k3] = TriTriangleDissectionCubic(qa/S,qb/S,qc/S)
        f(x) = k3*x^3 + k2*x^2 + k1*x + k0
        if abs(f(x0)) < 0.0000001:
            roots.append(x0)
    return roots

# TriTriangleGuessOneDissection returns root x for dissection a triangle in five parts when we know qa = BC^2 (quadrance of the base edge of ABC),
#  qb = CA^2 (quadrance of one another edge ABC)  qd = BD^2 = (1-x)^2 AB^2 (quadrance of one part of the last edge)
def  TriTriangleGuessOneDissection(qa,qb,qd):
    u = qb/qa ; v = qd/qa
    x = var('x')
    f1(x) = (u^4 - 4*u^3  + 6*u^2  - 4*u + 1)*x^10 + ( - 9*u^4 - 6*u^3*v  + 36*u^3 + 18*u^2*v - 54*u^2 - 18*u*v + 36*u + 6*v - 9)*x^9 + (37*u^4 + 50*u^3*v + 18*u^2*v^2 - 148*u^3 - 146*u^2*v - 36*u*v^2 + 222*u^2 + 142*u*v + 18*v^2 - 148*u - 46*v + 37)*x^8 + (- 92*u^4 - 188*u^3*v - 123*u^2*v^2 - 18*u*v^3 + 368*u^3 + 536*u^2*v + 234*u*v^2 + 18*v^3  - 552*u^2 - 508*u*v - 111*v^2 + 368*u + 160*v  - 92)*x^7 + (154*u^4 + 416*u^3*v + 363*u^2*v^2 + 90*u*v^3 + 9*v^4 - 616*u^3 - 1166*u^2*v - 668*u*v^2 - 90*v^3  + 924*u^2 + 1084*u*v + 309*v^2 - 616*u - 334*v + 154)*x^6 + (- 182*u^4 - 592*u^3*v - 597*u^2*v^2 - 186*u*v^3 - 36*v^4 + 728*u^3 + 1648*u^2*v + 1090*u*v^2 + 204*v^3 - 1092*u^2 - 1520*u*v - 513*v^2 + 728*u + 464*v - 182)*x^5+ (154*u^4 + 556*u^3*v + 591*u^2*v^2 + 206*u*v^3 + 63*v^4 - 616*u^3 - 1558*u^2*v - 1108*u*v^2 - 274*v^3 + 924*u^2 + 1448*u*v + 558*v^2 - 616*u - 446*v + 154)*x^4 + (- 92*u^4 - 340*u^3*v - 357*u^2*v^2 - 138*u*v^3 - 63*v^4 + 368*u^3 + 976*u^2*v + 718*u*v^2 + 234*v^3 - 552*u^2 - 932*u*v - 405*v^2  + 368*u + 296*v - 92)*x^3 + (37*u^4 + 128*u^3*v + 129*u^2*v^2 + 62*u*v^3 + 37*v^4 - 148*u^3 - 386*u^2*v - 292*u*v^2 - 124*v^3 + 222*u^2+ 388*u*v + 189*v^2 - 148*u - 130*v + 37)*x^2+ (- 9*u^4 - 26*u^3*v - 27*u^2*v^2 - 18*u*v^3 - 10*v^4 + 36*u^3 + 86*u^2*v + 70*u*v^2 + 36*v^3 - 54*u^2 - 94*u*v - 51*v^2 + 36*u + 34*v - 9)*x + (u^4 + 2*u^3*v + 3*u^2*v^2 + 2*u*v^3 + v^4 - 4*u^3 - 8*u^2*v - 8*u*v^2 - 4*v^3 + 6*u^2 + 10*u*v + 6*v^2 - 4*u - 4*v  + 1)
    # Searching first if root exist
    try:
        x0 = find_root(f1,0.001,0.999,maxiter=100)
    except RuntimeError:
        x0 = -1
    return x0

def TriTriangleDissectionRootByQuadrances(qa,qb,qc,S):
    u = qa/S ; v = qb/S ; w = qc/S
    [k0,k1,k2,k3] = TriEvalRealList(TriTriangleDissectionCubic(u,v,w))
    x = var('x')
    f(x) = k3*x^3 + k2*x^2 + k1*x + k0  #; print "f = ",f
    #
    # find unique real root between 0,1 :  D = (1-x)A + xB
    if k0 == 0:
        x = -1 # 0 is root
    else:
        if (k0 + k1 + k2 + k3) == 0:
            x = -1 # 1 is root
        else:
            try:
                x = find_root(f,0,1)
            except RuntimeError:
                x = -1 # no root between 0 and 1
    return x

def TriTriangleDissectionRoot(a,b,c,S):
    return TriTriangleDissectionRootByQuadrances(a^2,b^2,c^2,S)


# RT_TriTriangleDissectionFiveParts dissects a triangle in five obtuse (120 degrees) subtriangles, given x
#   tri_1 is point D on AB segment
#   tri_2 is point E, first Fermat point of triangl BCD
#   tri_3 is point F on CE such as DE and AF are parallel
# OLD :
#def RT_TriTriangleDissectionFiveParts(a,b,c,S,x):
#    tri_1 =  [-(x - 1)/a, x/b, 0]
#    tri_2 =  [(3*a^2*x - 3*b^2*x + 3*c^2*x - 4*sqrt(3)*S*x + 3*a^2 + 3*b^2 - 3*c^2 + 4*sqrt(3)*S)*(3*a^2 - 3*b^2 + 3*c^2 + 4*sqrt(3)*S)*(1 - x)/a,
#          (3*a^2*x - 3*b^2*x + 3*c^2*x - 4*sqrt(3)*S*x + 3*a^2 + 3*b^2 - 3*c^2 + 4*sqrt(3)*S)*(3*a^2*x - 3*b^2*x - 3*c^2*x + 4*sqrt(3)*S*x - 3*a^2 + 3*b^2 + 3*c^2 + 4*sqrt(3)*S)/b,
#          (6*c^2*x + 3*a^2 - 3*b^2 - 3*c^2 - 4*sqrt(3)*S)*(3*a^2 - 3*b^2 + 3*c^2 + 4*sqrt(3)*S)*(x - 1)/c]
#    tri_3 =  [(6*c^2*x + 3*a^2 - 3*b^2 - 3*c^2 + 4*sqrt(3)*S)*(3*a^2 - 3*b^2 + 3*c^2 + 4*sqrt(3)*S)*(1 - x)/a,
#          (3*a^2*x - 3*b^2*x - 3*c^2*x + 4*sqrt(3)*S*x - 3*a^2 + 3*b^2 + 3*c^2 + 4*sqrt(3)*S)*(6*c^2*x + 3*a^2 - 3*b^2 - 3*c^2 + 4*sqrt(3)*S)/b,
#          6*(6*c^2*x + 3*a^2 - 3*b^2 - 3*c^2 - 4*sqrt(3)*S)*c*(x - 1)]
#    return [tri_1,tri_2,tri_3]  # trilinear coordinates for D,E,F
# NEW
def RT_TriTriangleDissectionFiveParts(a,b,c,S,x):
    qa = a^2 ; qb = b^2 ; qc = c^2
    u = qa/qc ; v = qb/qc ; w = S/(qc*sqrt(3))
    z = u - v + 1 ; q = u + v - 1 ; r = v - u + 1
    tri_1 =  [(1 - x)/a, x/b, 0]
    # Removing factor 9*qc^2
    tri_2 =  [(x*z + q + 4*(1 - x)*w) * (z + 4*w) * (1 - x) /a,
              (x*z + q + 4*(1 - x)*w) * ((1 - x)*r + 4*(1 + x)*w) /b,
              (2*x - r - 4*w) * (z + 4*w) * (x - 1) /c]
    tri_3 =  [(2*x - r + 4*w) * (z + 4*w) * (1 - x) /a,
              ((1 - x)*r + 4*(1+x)*w) * (2*x - r + 4*w) /b,
              2*(2*x - r - 4*w) * (x - 1) /c]
    return [tri_1,tri_2,tri_3]  # trilinear coordinates for D,E,F


# TriTriangleDissectionFiveParts dissects a triangle in five obtuse (120 degrees) subtriangles, returning trilinear coordinates of :
#    D on AB edge , E internal point adjacent to 3 subtriangles, F internal point adjacent to 4 subtriangles
# F is such as AF is parallel with E given a cubic equation f = 0
def TriTriangleDissectionFiveParts(a,b,c,S):
    x = TriTriangleDissectionRoot(a,b,c,S)
    [tri_1,tri_2,tri_3] = RT_TriTriangleDissectionFiveParts(a,b,c,S,x)
    return [tri_1,tri_2,tri_3]  # trilinear coordinates for D,E,F

# TriTriangleDissectionCircle returns the squared circumradius and other quadrances for the triangle expansion
def TriTriangleDissectionCircle(qa,qb,qc,S,x):
    a = sqrt(qa) ; b = sqrt(qb) ; c = sqrt(qc)
    [tri_1,tri_2,tri_3] = RT_TriTriangleDissectionFiveParts(a,b,c,S,x) # trilinear coordinates for D,E,F
    [tri_4,tri_5,tri_6] = TriMainTriangle(a,b,c) # trilinear coordinates for A,B,C
    q34 = TriQuadranceTwoPoints(a,b,c,tri_3,tri_4) # Quad(A,F)
    q13 = TriQuadranceTwoPoints(a,b,c,tri_1,tri_3) # Quad(D,F)
    q15 = TriQuadranceTwoPoints(a,b,c,tri_1,tri_5) # Quad(B,D)
    q46 = TriQuadranceTwoPoints(a,b,c,tri_4,tri_6) # Quad(A,C)
    q56 = TriQuadranceTwoPoints(a,b,c,tri_5,tri_6) # Quad(B,C)
    qO = RT_CyclicQuadSquaredCircumradius(sqrt(q56),sqrt(q46),sqrt(q13) + sqrt(q34),sqrt(q15))
    return [qO,q15,q46]  #squared circumradius, Quad(B,D) and Quad(A,C)


def TriTriangleDissectionFivePartsQuadrances(qa,qb,qc,S,x):
    u = qa/qc ; v = qb/qc ; w = S/(qc*sqrt(3))
    # adding shortcuts
    q7 =  2/3*(u^2 - 2*u*v + v^2 - 2*u - 2*v + 1)*(u*x - v*x - 12*w*x + 2*x^2 + u + v + 12*w - 3*x + 1)/(2*u^2*x - 4*u*v*x + 2*v^2*x + 8*u*w*x - 8*v*w*x + 16*w*x^2 - 2*u^2 + 4*u*v - 2*v^2 + 8*u*w + 8*v*w - 4*u*x - 4*v*x - 24*w*x + 4*u + 4*v + 8*w + 2*x - 2)^2
    q8 = 4/3*(u - v + 4*w + 1)^2*(-u^2 - v^2 + 2*u*v + 2*u + 2*v - 1)/(u^2 - 2*u*v + v^2 - 2*u - 2*v + 1)
    q9 =  -2/3*(u*x - v*x - 12*w*x + 2*x^2 + u + v + 12*w - 3*x + 1)
    q0 =  -8/9*(4*u^2 - 8*u*v + 4*v^2 + 9*u*x - 9*v*x - 12*w*x + 6*x^2 + u + v + 12*w - 3*x + 1)
    # defining D from  q1 = qBD ; q2 = CD
    q1 = qc*(1 - x)^2
    q2 = qc*(x*(x-1) + u*x + v*(1 - x))
    # defining E from  q3 = qBE ; q4 = CE
    q3 = -qc*q7*(u - v + 4*w + 1)^2*(1 - x)^2
    q4 = -qc*q7*((u - v - 4*w + 1)*x + u + v + 4*w - 1)^2
    # defining F from  q5 = qBF ; q6 = CF
    q5 =  qc*q0*(1 - x)^2/q8
    q6 =  qc*q9*(u - v + 4*w + 2*x - 1)^2/q8
    return [q1,q2,q3,q4,q5,q6]

# TriTriangleDissectionFivePartsValidity returns boolean of validity of dissection according to the fact
# that F is in inside ABC triangle (equivalently on the good side of the line parallel with AB thru C)
# When x = 0 or x = 1, we have a special triangle and dissection is not valid
def TriTriangleDissectionFivePartsValidity(qa,qb,qc,S,x):
    u = qa/qc ; v = qb/qc ; w = S/(qc*sqrt(3))
    if (x < 0) or (x > 1):
        return False
    else:
        if (u^2 - 2*u*v + v^2 - 2*u + v + 1) == 0:
            return False # x = 0 is real root
        else:
            if (u^2 - 2*u*v + v^2 + u - 2*v + 1) == 0:
                return False # x = 1 is real root
            else:
                if (2*x + u - v - 1 + 4*w)*(x - 1 - 4*w) >= 0:
                    return False
                else:
                    return True

# TriCompleteQuad returns a complete quadrilateral wrt ABC given a point t with ta:ta:tb trilinear coordinates
def TriCompleteQuad(a,b,c,t):
    [ta,tb,tc] = t
    t1 = [ta, tb, 0] ; t2 = [0,1,0] ; t3 = [0, tb, tc] ; t4 = [ta, tb, tc]  # t1t2t3t4 is quadrilateral
    t5 = [ta, 2*tb, tc]                  # intersection of diagonals t1t3-t2t4
    t6 = [1,0,0] ; t7 = [0,0,1]          # intersection of extended edges
    t8 = [ta, 0, tc] ; t9 = [ta, 0, -tc] # two extra points in the harmonic line
    return [t1,t2,t3,t4,t5,t6,t7,t8,t9]

# TriQuadLunesCenters returns trilinear coordinates for centers of lunes thru t5 for quadrilateral t1t2t3t4
def TriQuadLunesCenters(a,b,c,t1,t2,t3,t4,t5):
    lbis12 = TriEvalReal(TriSegmentBissectorLine(a,b,c,t1,t2))
    lbis15 = TriEvalReal(TriSegmentBissectorLine(a,b,c,t1,t5))
    t12 = TriEvalReal(TriLineIntersect(a,b,c,lbis12,lbis15))
    lbis34 = TriEvalReal(TriSegmentBissectorLine(a,b,c,t3,t4))
    lbis35 = TriEvalReal(TriSegmentBissectorLine(a,b,c,t3,t5))
    t34 = TriEvalReal(TriLineIntersect(a,b,c,lbis34,lbis35))
    lbis23 = TriEvalReal(TriSegmentBissectorLine(a,b,c,t2,t3))
    t23 = TriEvalReal(TriLineIntersect(a,b,c,lbis23,lbis35))
    lbis41 = TriEvalReal(TriSegmentBissectorLine(a,b,c,t4,t1))
    t41 = TriEvalReal(TriLineIntersect(a,b,c,lbis41,lbis15))
    return [t12,t23,t34,t41]

# TriQuadBiMedians returns bimedians points for a quadrilateral t1t2t3t4
def TriQuadBiMedians(a,b,c,t1,t2,t3,t4):
    t12 = TriMidpoint(a,b,c,t1,t2)
    t34 = TriMidpoint(a,b,c,t3,t4)
    t23 = TriMidpoint(a,b,c,t2,t3)
    t41 = TriMidpoint(a,b,c,t4,t1)
    return [t12,t34,t23,t41]

# TriQuadCentroid returns centroid of a quad lamina as intersection of the two bimedians
def TriQuadCentroid(a,b,c,t1,t2,t3,t4):
    [t12,t34,t23,t41] = TriQuadBiMedians(a,b,c,t1,t2,t3,t4)
    t12 = TriEvalReal(t12); t34 = TriEvalReal(t34); t23 = TriEvalReal(t23); t41 = TriEvalReal(t41)
    lm1 = TriLine(a,b,c,t12,t34) ; lm1 = TriEvalReal(lm1)
    lm2 = TriLine(a,b,c,t23,t41) ; lm2 = TriEvalReal(lm2)
    t5 = TriLineIntersect(a,b,c,lm1,lm2) ; t5 = TriEvalReal(t5)
    return t5

# RT_TriQuadCentroid returns centroid of a quad lamina as intersection of the two bimedians
def RT_TriQuadCentroid(a,b,c,t1,t2,t3,t4):
    [t12,t34,t23,t41] = TriQuadBiMedians(a,b,c,t1,t2,t3,t4)
    t12 = TriFactor(t12); t34 = TriFactor(t34); t23 = TriFactor(t23); t41 = TriFactor(t41)
    lm1 = TriLine(a,b,c,t12,t34) ; lm1 = TriFactor(lm1)
    lm2 = TriLine(a,b,c,t23,t41) ; lm2 = TriFactor(lm2)
    t5 = TriLineIntersect(a,b,c,lm1,lm2) ; t5 = TriFactor(t5)
    return t5

# TriCircularSectorCentroid returns centroid of a circular sector on edge t1-t2, angle gamma
def RT_TriCircularSectorCentroid(a,b,c,S,t1,t2,gamma,cs_gamma,sn_gamma):
    q = (a*cs_gamma - 2*a*sn_gamma/(3*gamma))^2
    [t3,t9] = TriTriangleApexVertex(a,b,c,S,q,t1,t2)
    return t3

def TriCircularSectorCentroid(a,b,c,S,t1,t2,gamma):
    cs_gamma = cos(gamma).n() ; sn_gamma = sin(gamma).n()
    t4 = TriEvalReal(t1); t5 = TriEvalReal(t2)
    t3 = TriEvalReal(RT_TriCircularSectorCentroid(a,b,c,S,t4,t5,gamma.n(),cs_gamma,sn_gamma))
    return t3


# TriHTree returns list of parallelograms in the H-tree recursive decomposition of a parallelogram given by
# three corners t1,t2,t3
# 'edges' parameter is a list of couples (e[0],e[1]) of two trilinear coordinates
# 'depth' parameter is depth of the H-tree
def TriHTree(a,b,c,depth,t1,t2,t3,edges):
    t4 = TriMidpoint(a,b,c,t1,t3) # center of symmetry of parallelogram (intersection of diagonals)
    if depth > 0:
        t5 = TriMidpoint(a,b,c,t2,t3)  # midpoint of one edge of parallelogram
        t6 = TriReflectionPoint(a,b,c,t4,t5) # midpoint on the opposite edge
        t7 = TriReflectionPoint(a,b,c,t4,t2) # last corner of parallelogram t1-t2-t3-t7
        t8 = TriMidpoint(a,b,c,t2,t6) # centroid of 1st sub-parallelogram
        t9 = TriMidpoint(a,b,c,t5,t7) # centroid of 2nd sub-parallelogram
        TriHTree(a,b,c,depth - 1,t6,t1,t2,edges)
        TriHTree(a,b,c,depth - 1,t5,t3,t7,edges)
        edges.append([t8,t9]) # add edge between centroids


# TriConicEquation returns equation of conic in trilinear coordinates given five points t1,t2,t3,t4,t5 in
# trilinear coordinates
# We use ray of conics. Points t1 = [ta1,tb1,tc1] ... t5 = [ta5,tb5,tc5] with t1t2t3t4 a not self-crossing quad
# Equations of lines are then :
# l12 =  [(tb2*tc1 - tb1*tc2), (ta1*tc2 - ta2*tc1), (ta2*tb1 - ta1*tb2)]
# l34 =  [(tb4*tc3 - tb3*tc4), (ta3*tc4 - ta4*tc3), (ta4*tb3 - ta3*tb4)]
# l14 =  [(tb4*tc1 - tb1*tc4), (ta1*tc4 - ta4*tc1), (ta4*tb1 - ta1*tb4)]
# l23 =  [(tb3*tc2 - tb2*tc3), (ta2*tc3 - ta3*tc2), (ta3*tb2 - ta2*tb3)]
# Point t is on conic from the ray given by parameter delta iff :
#   (TriOnLine(t,l12)* TriOnLine(t,l34) + delta*TriOnLine(t,l14)* TriOnLine(t,l23)) = 0
# because we simply multiply equations to get lines equations to get the degenerate conic (l12 = 0,l34 = 0) and to
# get the degenerate conic (l14 = 0,l23 = 0)
# The fifth point t5 helps for getting parameter delta, substituting t=t5 in the equation.
#
# Because equation is very long to write we compute coefficients k2,k1 of equation in alpha,beta,gamma for point t=[alpha,beta,gamma]
#  k2*alpha^2 + k1*alpha*beta + ... + 0 = 0
# in one specific function TriConicEquationParameters and get the other ones by permutation
# Last coefficient for the alpha^0*beta^0*gamma^0 term is 0.

# TriConicEquationParameters returns k1,k2 parameters using derivative()  (double partial derivayes in alpha, beta )

def TriConicEquationParameters(t1,t2,t3,t4,t5):
    [ta1,tb1,tc1] = t1; [ta2,tb2,tc2] = t2; [ta3,tb3,tc3] = t3; [ta4,tb4,tc4] = t4; [ta5,tb5,tc5] = t5
    k2 =  - ta3*ta5*tb2*tb4^2*tb5*tc1^2*tc2*tc3 + ta2*ta5*tb3*tb4^2*tb5*tc1^2*tc2*tc3 + ta3*ta4*tb2*tb4*tb5^2*tc1^2*tc2*tc3
    k2 += - ta2*ta4*tb3*tb4*tb5^2*tc1^2*tc2*tc3 + ta3*ta5*tb1*tb4^2*tb5*tc1*tc2^2*tc3 - ta1*ta5*tb3*tb4^2*tb5*tc1*tc2^2*tc3
    k2 += - ta3*ta4*tb1*tb4*tb5^2*tc1*tc2^2*tc3 + ta1*ta4*tb3*tb4*tb5^2*tc1*tc2^2*tc3 - ta2*ta5*tb1*tb4^2*tb5*tc1*tc2*tc3^2
    k2 += + ta1*ta5*tb2*tb4^2*tb5*tc1*tc2*tc3^2 + ta2*ta4*tb1*tb4*tb5^2*tc1*tc2*tc3^2 - ta1*ta4*tb2*tb4*tb5^2*tc1*tc2*tc3^2
    k2 += + ta4*ta5*tb2*tb3^2*tb5*tc1^2*tc2*tc4 - ta2*ta5*tb3^2*tb4*tb5*tc1^2*tc2*tc4 - ta3*ta4*tb2*tb3*tb5^2*tc1^2*tc2*tc4
    k2 += + ta2*ta3*tb3*tb4*tb5^2*tc1^2*tc2*tc4 - ta4*ta5*tb1*tb3^2*tb5*tc1*tc2^2*tc4 + ta1*ta5*tb3^2*tb4*tb5*tc1*tc2^2*tc4
    k2 += + ta3*ta4*tb1*tb3*tb5^2*tc1*tc2^2*tc4 - ta1*ta3*tb3*tb4*tb5^2*tc1*tc2^2*tc4 - ta4*ta5*tb2^2*tb3*tb5*tc1^2*tc3*tc4
    k2 += + ta3*ta5*tb2^2*tb4*tb5*tc1^2*tc3*tc4 + ta2*ta4*tb2*tb3*tb5^2*tc1^2*tc3*tc4 - ta2*ta3*tb2*tb4*tb5^2*tc1^2*tc3*tc4
    k2 += + ta4*ta5*tb1^2*tb3*tb5*tc2^2*tc3*tc4 - ta3*ta5*tb1^2*tb4*tb5*tc2^2*tc3*tc4 - ta1*ta4*tb1*tb3*tb5^2*tc2^2*tc3*tc4
    k2 += + ta1*ta3*tb1*tb4*tb5^2*tc2^2*tc3*tc4 + ta4*ta5*tb1*tb2^2*tb5*tc1*tc3^2*tc4 - ta1*ta5*tb2^2*tb4*tb5*tc1*tc3^2*tc4
    k2 += - ta2*ta4*tb1*tb2*tb5^2*tc1*tc3^2*tc4 + ta1*ta2*tb2*tb4*tb5^2*tc1*tc3^2*tc4 - ta4*ta5*tb1^2*tb2*tb5*tc2*tc3^2*tc4
    k2 += + ta2*ta5*tb1^2*tb4*tb5*tc2*tc3^2*tc4 + ta1*ta4*tb1*tb2*tb5^2*tc2*tc3^2*tc4 - ta1*ta2*tb1*tb4*tb5^2*tc2*tc3^2*tc4
    k2 += + ta2*ta5*tb1*tb3^2*tb5*tc1*tc2*tc4^2 - ta1*ta5*tb2*tb3^2*tb5*tc1*tc2*tc4^2 - ta2*ta3*tb1*tb3*tb5^2*tc1*tc2*tc4^2
    k2 += + ta1*ta3*tb2*tb3*tb5^2*tc1*tc2*tc4^2 - ta3*ta5*tb1*tb2^2*tb5*tc1*tc3*tc4^2 + ta1*ta5*tb2^2*tb3*tb5*tc1*tc3*tc4^2
    k2 += + ta2*ta3*tb1*tb2*tb5^2*tc1*tc3*tc4^2 - ta1*ta2*tb2*tb3*tb5^2*tc1*tc3*tc4^2 + ta3*ta5*tb1^2*tb2*tb5*tc2*tc3*tc4^2
    k2 += - ta2*ta5*tb1^2*tb3*tb5*tc2*tc3*tc4^2 - ta1*ta3*tb1*tb2*tb5^2*tc2*tc3*tc4^2 + ta1*ta2*tb1*tb3*tb5^2*tc2*tc3*tc4^2
    k2 += - ta4*ta5*tb2*tb3^2*tb4*tc1^2*tc2*tc5 + ta3*ta5*tb2*tb3*tb4^2*tc1^2*tc2*tc5 + ta2*ta4*tb3^2*tb4*tb5*tc1^2*tc2*tc5
    k2 += - ta2*ta3*tb3*tb4^2*tb5*tc1^2*tc2*tc5 + ta4*ta5*tb1*tb3^2*tb4*tc1*tc2^2*tc5 - ta3*ta5*tb1*tb3*tb4^2*tc1*tc2^2*tc5
    k2 += - ta1*ta4*tb3^2*tb4*tb5*tc1*tc2^2*tc5 + ta1*ta3*tb3*tb4^2*tb5*tc1*tc2^2*tc5 + ta4*ta5*tb2^2*tb3*tb4*tc1^2*tc3*tc5
    k2 += - ta2*ta5*tb2*tb3*tb4^2*tc1^2*tc3*tc5 - ta3*ta4*tb2^2*tb4*tb5*tc1^2*tc3*tc5 + ta2*ta3*tb2*tb4^2*tb5*tc1^2*tc3*tc5
    k2 += - ta4*ta5*tb1^2*tb3*tb4*tc2^2*tc3*tc5 + ta1*ta5*tb1*tb3*tb4^2*tc2^2*tc3*tc5 + ta3*ta4*tb1^2*tb4*tb5*tc2^2*tc3*tc5
    k2 += - ta1*ta3*tb1*tb4^2*tb5*tc2^2*tc3*tc5 - ta4*ta5*tb1*tb2^2*tb4*tc1*tc3^2*tc5 + ta2*ta5*tb1*tb2*tb4^2*tc1*tc3^2*tc5
    k2 += + ta1*ta4*tb2^2*tb4*tb5*tc1*tc3^2*tc5 - ta1*ta2*tb2*tb4^2*tb5*tc1*tc3^2*tc5 + ta4*ta5*tb1^2*tb2*tb4*tc2*tc3^2*tc5
    k2 += - ta1*ta5*tb1*tb2*tb4^2*tc2*tc3^2*tc5 - ta2*ta4*tb1^2*tb4*tb5*tc2*tc3^2*tc5 + ta1*ta2*tb1*tb4^2*tb5*tc2*tc3^2*tc5
    k2 += - ta3*ta5*tb2^2*tb3*tb4*tc1^2*tc4*tc5 + ta2*ta5*tb2*tb3^2*tb4*tc1^2*tc4*tc5 + ta3*ta4*tb2^2*tb3*tb5*tc1^2*tc4*tc5
    k2 += - ta2*ta4*tb2*tb3^2*tb5*tc1^2*tc4*tc5 + ta3*ta5*tb1^2*tb3*tb4*tc2^2*tc4*tc5 - ta1*ta5*tb1*tb3^2*tb4*tc2^2*tc4*tc5
    k2 += - ta3*ta4*tb1^2*tb3*tb5*tc2^2*tc4*tc5 + ta1*ta4*tb1*tb3^2*tb5*tc2^2*tc4*tc5 - ta2*ta5*tb1^2*tb2*tb4*tc3^2*tc4*tc5
    k2 += + ta1*ta5*tb1*tb2^2*tb4*tc3^2*tc4*tc5 + ta2*ta4*tb1^2*tb2*tb5*tc3^2*tc4*tc5 - ta1*ta4*tb1*tb2^2*tb5*tc3^2*tc4*tc5
    k2 += + ta3*ta5*tb1*tb2^2*tb3*tc1*tc4^2*tc5 - ta2*ta5*tb1*tb2*tb3^2*tc1*tc4^2*tc5 - ta1*ta3*tb2^2*tb3*tb5*tc1*tc4^2*tc5
    k2 += + ta1*ta2*tb2*tb3^2*tb5*tc1*tc4^2*tc5 - ta3*ta5*tb1^2*tb2*tb3*tc2*tc4^2*tc5 + ta1*ta5*tb1*tb2*tb3^2*tc2*tc4^2*tc5
    k2 += + ta2*ta3*tb1^2*tb3*tb5*tc2*tc4^2*tc5 - ta1*ta2*tb1*tb3^2*tb5*tc2*tc4^2*tc5 + ta2*ta5*tb1^2*tb2*tb3*tc3*tc4^2*tc5
    k2 += - ta1*ta5*tb1*tb2^2*tb3*tc3*tc4^2*tc5 - ta2*ta3*tb1^2*tb2*tb5*tc3*tc4^2*tc5 + ta1*ta3*tb1*tb2^2*tb5*tc3*tc4^2*tc5
    k2 += - ta2*ta4*tb1*tb3^2*tb4*tc1*tc2*tc5^2 + ta1*ta4*tb2*tb3^2*tb4*tc1*tc2*tc5^2 + ta2*ta3*tb1*tb3*tb4^2*tc1*tc2*tc5^2
    k2 += - ta1*ta3*tb2*tb3*tb4^2*tc1*tc2*tc5^2 + ta3*ta4*tb1*tb2^2*tb4*tc1*tc3*tc5^2 - ta1*ta4*tb2^2*tb3*tb4*tc1*tc3*tc5^2
    k2 += - ta2*ta3*tb1*tb2*tb4^2*tc1*tc3*tc5^2 + ta1*ta2*tb2*tb3*tb4^2*tc1*tc3*tc5^2 - ta3*ta4*tb1^2*tb2*tb4*tc2*tc3*tc5^2
    k2 += + ta2*ta4*tb1^2*tb3*tb4*tc2*tc3*tc5^2 + ta1*ta3*tb1*tb2*tb4^2*tc2*tc3*tc5^2 - ta1*ta2*tb1*tb3*tb4^2*tc2*tc3*tc5^2
    k2 += - ta3*ta4*tb1*tb2^2*tb3*tc1*tc4*tc5^2 + ta2*ta4*tb1*tb2*tb3^2*tc1*tc4*tc5^2 + ta1*ta3*tb2^2*tb3*tb4*tc1*tc4*tc5^2
    k2 += - ta1*ta2*tb2*tb3^2*tb4*tc1*tc4*tc5^2 + ta3*ta4*tb1^2*tb2*tb3*tc2*tc4*tc5^2 - ta1*ta4*tb1*tb2*tb3^2*tc2*tc4*tc5^2
    k2 += - ta2*ta3*tb1^2*tb3*tb4*tc2*tc4*tc5^2 + ta1*ta2*tb1*tb3^2*tb4*tc2*tc4*tc5^2 - ta2*ta4*tb1^2*tb2*tb3*tc3*tc4*tc5^2
    k2 += + ta1*ta4*tb1*tb2^2*tb3*tc3*tc4*tc5^2 + ta2*ta3*tb1^2*tb2*tb4*tc3*tc4*tc5^2 - ta1*ta3*tb1*tb2^2*tb4*tc3*tc4*tc5^2
    #
    k1 =    ta3*ta5^2*tb2*tb4^2*tc1^2*tc2*tc3 - ta2*ta5^2*tb3*tb4^2*tc1^2*tc2*tc3 - ta3*ta4^2*tb2*tb5^2*tc1^2*tc2*tc3
    k1 += + ta2*ta4^2*tb3*tb5^2*tc1^2*tc2*tc3 - ta3*ta5^2*tb1*tb4^2*tc1*tc2^2*tc3 + ta1*ta5^2*tb3*tb4^2*tc1*tc2^2*tc3
    k1 += + ta3*ta4^2*tb1*tb5^2*tc1*tc2^2*tc3 - ta1*ta4^2*tb3*tb5^2*tc1*tc2^2*tc3 + ta2*ta5^2*tb1*tb4^2*tc1*tc2*tc3^2
    k1 += - ta1*ta5^2*tb2*tb4^2*tc1*tc2*tc3^2 - ta2*ta4^2*tb1*tb5^2*tc1*tc2*tc3^2 + ta1*ta4^2*tb2*tb5^2*tc1*tc2*tc3^2
    k1 += - ta4*ta5^2*tb2*tb3^2*tc1^2*tc2*tc4 + ta2*ta5^2*tb3^2*tb4*tc1^2*tc2*tc4 + ta3^2*ta4*tb2*tb5^2*tc1^2*tc2*tc4
    k1 += - ta2*ta3^2*tb4*tb5^2*tc1^2*tc2*tc4 + ta4*ta5^2*tb1*tb3^2*tc1*tc2^2*tc4 - ta1*ta5^2*tb3^2*tb4*tc1*tc2^2*tc4
    k1 += - ta3^2*ta4*tb1*tb5^2*tc1*tc2^2*tc4 + ta1*ta3^2*tb4*tb5^2*tc1*tc2^2*tc4 + ta4*ta5^2*tb2^2*tb3*tc1^2*tc3*tc4
    k1 += - ta3*ta5^2*tb2^2*tb4*tc1^2*tc3*tc4 - ta2^2*ta4*tb3*tb5^2*tc1^2*tc3*tc4 + ta2^2*ta3*tb4*tb5^2*tc1^2*tc3*tc4
    k1 += - ta4*ta5^2*tb1^2*tb3*tc2^2*tc3*tc4 + ta3*ta5^2*tb1^2*tb4*tc2^2*tc3*tc4 + ta1^2*ta4*tb3*tb5^2*tc2^2*tc3*tc4
    k1 += - ta1^2*ta3*tb4*tb5^2*tc2^2*tc3*tc4 - ta4*ta5^2*tb1*tb2^2*tc1*tc3^2*tc4 + ta1*ta5^2*tb2^2*tb4*tc1*tc3^2*tc4
    k1 += + ta2^2*ta4*tb1*tb5^2*tc1*tc3^2*tc4 - ta1*ta2^2*tb4*tb5^2*tc1*tc3^2*tc4 + ta4*ta5^2*tb1^2*tb2*tc2*tc3^2*tc4
    k1 += - ta2*ta5^2*tb1^2*tb4*tc2*tc3^2*tc4 - ta1^2*ta4*tb2*tb5^2*tc2*tc3^2*tc4 + ta1^2*ta2*tb4*tb5^2*tc2*tc3^2*tc4
    k1 += - ta2*ta5^2*tb1*tb3^2*tc1*tc2*tc4^2 + ta1*ta5^2*tb2*tb3^2*tc1*tc2*tc4^2 + ta2*ta3^2*tb1*tb5^2*tc1*tc2*tc4^2
    k1 += - ta1*ta3^2*tb2*tb5^2*tc1*tc2*tc4^2 + ta3*ta5^2*tb1*tb2^2*tc1*tc3*tc4^2 - ta1*ta5^2*tb2^2*tb3*tc1*tc3*tc4^2
    k1 += - ta2^2*ta3*tb1*tb5^2*tc1*tc3*tc4^2 + ta1*ta2^2*tb3*tb5^2*tc1*tc3*tc4^2 - ta3*ta5^2*tb1^2*tb2*tc2*tc3*tc4^2
    k1 += + ta2*ta5^2*tb1^2*tb3*tc2*tc3*tc4^2 + ta1^2*ta3*tb2*tb5^2*tc2*tc3*tc4^2 - ta1^2*ta2*tb3*tb5^2*tc2*tc3*tc4^2
    k1 += + ta4^2*ta5*tb2*tb3^2*tc1^2*tc2*tc5 - ta3^2*ta5*tb2*tb4^2*tc1^2*tc2*tc5 - ta2*ta4^2*tb3^2*tb5*tc1^2*tc2*tc5
    k1 += + ta2*ta3^2*tb4^2*tb5*tc1^2*tc2*tc5 - ta4^2*ta5*tb1*tb3^2*tc1*tc2^2*tc5 + ta3^2*ta5*tb1*tb4^2*tc1*tc2^2*tc5
    k1 += + ta1*ta4^2*tb3^2*tb5*tc1*tc2^2*tc5 - ta1*ta3^2*tb4^2*tb5*tc1*tc2^2*tc5 - ta4^2*ta5*tb2^2*tb3*tc1^2*tc3*tc5
    k1 += + ta2^2*ta5*tb3*tb4^2*tc1^2*tc3*tc5 + ta3*ta4^2*tb2^2*tb5*tc1^2*tc3*tc5 - ta2^2*ta3*tb4^2*tb5*tc1^2*tc3*tc5
    k1 += + ta4^2*ta5*tb1^2*tb3*tc2^2*tc3*tc5 - ta1^2*ta5*tb3*tb4^2*tc2^2*tc3*tc5 - ta3*ta4^2*tb1^2*tb5*tc2^2*tc3*tc5
    k1 += + ta1^2*ta3*tb4^2*tb5*tc2^2*tc3*tc5 + ta4^2*ta5*tb1*tb2^2*tc1*tc3^2*tc5 - ta2^2*ta5*tb1*tb4^2*tc1*tc3^2*tc5
    k1 += - ta1*ta4^2*tb2^2*tb5*tc1*tc3^2*tc5 + ta1*ta2^2*tb4^2*tb5*tc1*tc3^2*tc5 - ta4^2*ta5*tb1^2*tb2*tc2*tc3^2*tc5
    k1 += + ta1^2*ta5*tb2*tb4^2*tc2*tc3^2*tc5 + ta2*ta4^2*tb1^2*tb5*tc2*tc3^2*tc5 - ta1^2*ta2*tb4^2*tb5*tc2*tc3^2*tc5
    k1 += + ta3^2*ta5*tb2^2*tb4*tc1^2*tc4*tc5 - ta2^2*ta5*tb3^2*tb4*tc1^2*tc4*tc5 - ta3^2*ta4*tb2^2*tb5*tc1^2*tc4*tc5
    k1 += + ta2^2*ta4*tb3^2*tb5*tc1^2*tc4*tc5 - ta3^2*ta5*tb1^2*tb4*tc2^2*tc4*tc5 + ta1^2*ta5*tb3^2*tb4*tc2^2*tc4*tc5
    k1 += + ta3^2*ta4*tb1^2*tb5*tc2^2*tc4*tc5 - ta1^2*ta4*tb3^2*tb5*tc2^2*tc4*tc5 + ta2^2*ta5*tb1^2*tb4*tc3^2*tc4*tc5
    k1 += - ta1^2*ta5*tb2^2*tb4*tc3^2*tc4*tc5 - ta2^2*ta4*tb1^2*tb5*tc3^2*tc4*tc5 + ta1^2*ta4*tb2^2*tb5*tc3^2*tc4*tc5
    k1 += - ta3^2*ta5*tb1*tb2^2*tc1*tc4^2*tc5 + ta2^2*ta5*tb1*tb3^2*tc1*tc4^2*tc5 + ta1*ta3^2*tb2^2*tb5*tc1*tc4^2*tc5
    k1 += - ta1*ta2^2*tb3^2*tb5*tc1*tc4^2*tc5 + ta3^2*ta5*tb1^2*tb2*tc2*tc4^2*tc5 - ta1^2*ta5*tb2*tb3^2*tc2*tc4^2*tc5
    k1 += - ta2*ta3^2*tb1^2*tb5*tc2*tc4^2*tc5 + ta1^2*ta2*tb3^2*tb5*tc2*tc4^2*tc5 - ta2^2*ta5*tb1^2*tb3*tc3*tc4^2*tc5
    k1 += + ta1^2*ta5*tb2^2*tb3*tc3*tc4^2*tc5 + ta2^2*ta3*tb1^2*tb5*tc3*tc4^2*tc5 - ta1^2*ta3*tb2^2*tb5*tc3*tc4^2*tc5
    k1 += + ta2*ta4^2*tb1*tb3^2*tc1*tc2*tc5^2 - ta1*ta4^2*tb2*tb3^2*tc1*tc2*tc5^2 - ta2*ta3^2*tb1*tb4^2*tc1*tc2*tc5^2
    k1 += + ta1*ta3^2*tb2*tb4^2*tc1*tc2*tc5^2 - ta3*ta4^2*tb1*tb2^2*tc1*tc3*tc5^2 + ta1*ta4^2*tb2^2*tb3*tc1*tc3*tc5^2
    k1 += + ta2^2*ta3*tb1*tb4^2*tc1*tc3*tc5^2 - ta1*ta2^2*tb3*tb4^2*tc1*tc3*tc5^2 + ta3*ta4^2*tb1^2*tb2*tc2*tc3*tc5^2
    k1 += - ta2*ta4^2*tb1^2*tb3*tc2*tc3*tc5^2 - ta1^2*ta3*tb2*tb4^2*tc2*tc3*tc5^2 + ta1^2*ta2*tb3*tb4^2*tc2*tc3*tc5^2
    k1 += + ta3^2*ta4*tb1*tb2^2*tc1*tc4*tc5^2 - ta2^2*ta4*tb1*tb3^2*tc1*tc4*tc5^2 - ta1*ta3^2*tb2^2*tb4*tc1*tc4*tc5^2
    k1 += + ta1*ta2^2*tb3^2*tb4*tc1*tc4*tc5^2 - ta3^2*ta4*tb1^2*tb2*tc2*tc4*tc5^2 + ta1^2*ta4*tb2*tb3^2*tc2*tc4*tc5^2
    k1 += + ta2*ta3^2*tb1^2*tb4*tc2*tc4*tc5^2 - ta1^2*ta2*tb3^2*tb4*tc2*tc4*tc5^2 + ta2^2*ta4*tb1^2*tb3*tc3*tc4*tc5^2
    k1 += - ta1^2*ta4*tb2^2*tb3*tc3*tc4*tc5^2 - ta2^2*ta3*tb1^2*tb4*tc3*tc4*tc5^2 + ta1^2*ta3*tb2^2*tb4*tc3*tc4*tc5^2
    return [k1,k2]

def TriConicEquation(t1,t2,t3,t4,t5):
    [ta1,tb1,tc1] = t1; [ta2,tb2,tc2] = t2; [ta3,tb3,tc3] = t3; [ta4,tb4,tc4] = t4; [ta5,tb5,tc5] = t5
    [kab,kaa] = TriConicEquationParameters([ta1,tb1,tc1],[ta2,tb2,tc2],[ta3,tb3,tc3],[ta4,tb4,tc4],[ta5,tb5,tc5])
    [kbc,kbb] = TriConicEquationParameters([tb1,tc1,ta1],[tb2,tc2,ta2],[tb3,tc3,ta3],[tb4,tc4,ta4],[tb5,tc5,ta5])
    [kca,kcc] = TriConicEquationParameters([tc1,ta1,tb1],[tc2,ta2,tb2],[tc3,ta3,tb3],[tc4,ta4,tb4],[tc5,ta5,tb5])
    # equation :  kaa*alpha^2 + kbb*beta^2 + kcc*gamma^2 + kbc*beta*gamma + kca*gamma*alpha + kab*alpha*beta = 0
    #             for point [alpha,beta,gamma]
    return [kaa,kbb,kcc,kab,kbc,kca]

# TriIntersectConicLine returns two points of intersection between conic and a line
def TriIntersectConicLine(k,l):
    [kaa,kbb,kcc,kab,kbc,kca] = k
    [la,lb,lc] = l
    #
    disc2 =  kbc^2*la^2 - 4*kbb*kcc*la^2 - 2*kbc*kca*la*lb + 4*kab*kcc*la*lb + kca^2*lb^2 - 4*kaa*kcc*lb^2 - 2*kab*kbc*la*lc + 4*kbb*kca*la*lc + 4*kaa*kbc*lb*lc - 2*kab*kca*lb*lc + kab^2*lc^2 - 4*kaa*kbb*lc^2 ; disc = sqrt(disc2)
    #
    u = 2*(kcc*lb^2 - kbc*lb*lc + kbb*lc^2)
    if u != 0:
        v = -(2*kcc*la*lb - kbc*la*lc - kca*lb*lc + kab*lc^2)
        w = kbc*la*lb - kca*lb^2 - 2*kbb*la*lc + kab*lb*lc
        t1 = [u,v - lc*disc, w + lb*disc] ; t2 = [u,v + lc*disc, w - lb*disc]
    else:
        u = 2*(kaa*lc^2 - kca*lc*la + kcc*la^2)
        if u != 0:
            v = -(2*kaa*lb*lc - kca*lb*la - kab*lc*la + kbc*la^2)
            w = kca*lb*lc - kab*lc^2 - 2*kcc*lb*la + kbc*lc*la
            t1 = [w + lc*disc, u, v - la*disc] ; t2 = [w - lc*disc ,u, v + la*disc]
        else:
            u = 2*(kbb*la^2 - kab*la*lb + kaa*lb^2)
            if u != 0:
                v = -(2*kbb*lc*la - kab*lc*lb - kbc*la*lb + kca*lb^2)
                w = kab*lc*la - kbc*la^2 - 2*kaa*lc*lb + kca*la*lb
                t1 = [v - lb*disc, w + la*disc, u] ; t2 = [v + lb*disc, w - la*disc ,u]
            else:
                if la^2 + lb^2 == 0:
                    if kaa != 0:
                        disc2 = kab^2 - 4*kaa*kbb ; disc = sqrt(disc2)
                        t1 = [kab + disc ,-2*kaa,0] ; t2 = [kab - disc ,-2*kaa,0]
                    else:
                        t1 = [1,0,0] ; t2 = [0,1,0]
                else:
                    if la^2 + lc^2 == 0:
                        if kcc != 0:
                            disc2 = kca^2 - 4*kaa*kcc ; disc = sqrt(disc2)
                            t1 = [-2*kcc,0,kca + disc] ; t2 = [-2*kcc, 0, kca - disc]
                        else:
                            t1 = [1,0,0] ; t2 = [0,0,1]
                    else:
                        if lb^2 + lc^2 == 0:
                            if kbb != 0:
                                disc2 = kbc^2 - 4*kbb*kcc ; disc = sqrt(disc2)
                                t1 = [0, kbc + disc  -2*kbb] ; t2 = [0, kbc - disc, -2*kbb]
                            else:
                                t1 = [0,1,0] ; t2 = [0,0,1]
                        else:
                            print "** NOT YET IMPLEMENTED **"
                            t1 = [0,0,0] ; t2 = [0,0,0] # TBD : special case
    #
    return [t1,t2]

# TriConicTangent returns tangent line for a conic at a given point
def TriConicTangent(k,t):
    [kaa,kbb,kcc,kab,kbc,kca] = k
    [alpha,beta,gamma] = t
    #
    la = beta^2*kab^2 - 4*beta^2*kaa*kbb - 4*beta*gamma*kaa*kbc + 2*beta*gamma*kab*kca + gamma^2*kca^2 - 4*gamma^2*kaa*kcc
    if la != 0:
        lb = -alpha*beta*kab^2 + 4*alpha*beta*kaa*kbb + 2*alpha*gamma*kaa*kbc - beta*gamma*kab*kbc - alpha*gamma*kab*kca + 2*beta*gamma*kbb*kca + gamma^2*kbc*kca - 2*gamma^2*kab*kcc
        lc = 2*alpha*beta*kaa*kbc + beta^2*kab*kbc - alpha*beta*kab*kca - 2*beta^2*kbb*kca - beta*gamma*kbc*kca - alpha*gamma*kca^2 + 4*alpha*gamma*kaa*kcc + 2*beta*gamma*kab*kcc
    else:
        lb = kab^2 - 4*kaa*kbb
        lc = -2*kaa*kbc + kab*kca
    #
    return [la,lb,lc]

# TriConicCentroid returns center of the conic
def TriConicCentroid(a,b,c,k):
    [kaa,kbb,kcc,kab,kbc,kca] = k
    u = kaa ; v = kbb ; w = kcc ; f = kbc/2 ; g = kca/2 ; h = kab/2
    t = TriConicCenter(a,b,c,u,v,w,f,g,h)
    return t

# TriScaledSteinerCircumellipse returns parameters for a scaled Steiner Circumellipse of ABC triangle, factor q = k^2
def TriScaledSteinerCircumellipse(a,b,c,q):
    kaa = (q - 1)*a^2 ; kbb = (q - 1)*b^2 ; kcc = (q - 1)*c^2
    kab = (2*q + 1)*a*b ; kbc = (2*q + 1)*b*c ; kca = (2*q + 1)*a*c
    return [kaa,kbb,kcc,kab,kbc,kca]

# TriScaledSteinerCircumellipsePoints returns five point on the scaled circumellipse.
# the five points are images of A,B,C an two other points under a circle inversion, center G of ellipse, and radius k
def TriScaledSteinerCircumellipsePoints(a,b,c,k):
    w = (1 - 2*k)/(1 + k) ; z = (1 - k)/(1 + 2*k)
    t1 =  [1/a, z/b, z/c]
    t2 =  [z/a, 1/b, z/c]
    t3 =  [z/a, z/b, 1/c]
    t4 =  [1/a, w/b, 1/c]
    t5 =  [1/a, 1/b, w/c]
    return [t1,t2,t3,t4,t5]

# TriSteinerCircumellipseAxis returns axis of Steiner Circumellipse of reference triangle ABC
def TriSteinerCircumellipseAxis(a,b,c):
    w = sqrt(a^4 + b^4 + c^4 - a^2*b^2  - a^2*c^2 - b^2*c^2)
    major =  [(a^4*b^2 - 2*a^2*b^4 + a^4*c^2 + b^4*c^2 - 2*a^2*c^4 + b^2*c^4 + w*(a^2*b^2 + a^2*c^2 - 2*b^2*c^2))*a,
              -(2*a^4*b^2 - a^2*b^4 - a^4*c^2 - b^4*c^2 - a^2*c^4 + 2*b^2*c^4 - w*(a^2*b^2 - 2*a^2*c^2 + b^2*c^2))*b,
              (a^4*b^2 + a^2*b^4 - 2*a^4*c^2 - 2*b^4*c^2 + a^2*c^4 + b^2*c^4 + w*(- 2*a^2*b^2 + a^2*c^2 + b^2*c^2))*c]
    minor =  [(a^6 - 3*a^4*b^2 + 2*a^2*b^4 - 3*a^4*c^2 - a^2*b^2*c^2 + b^4*c^2 + 2*a^2*c^4 + b^2*c^4 - a^4*w + a^2*b^2*w + a^2*c^2*w + 2*b^2*c^2*w)*a*(b + c)*(b - c),
              -(2*a^4*b^2 - 3*a^2*b^4 + b^6 + a^4*c^2 - a^2*b^2*c^2 - 3*b^4*c^2 + a^2*c^4 + 2*b^2*c^4 + a^2*b^2*w - b^4*w + 2*a^2*c^2*w + b^2*c^2*w)*(a + c)*(a - c)*b,
              (a^4*b^2 + a^2*b^4 + 2*a^4*c^2 - a^2*b^2*c^2 + 2*b^4*c^2 - 3*a^2*c^4 - 3*b^2*c^4 + c^6 + 2*a^2*b^2*w + a^2*c^2*w + b^2*c^2*w - c^4*w)*(a + b)*(a - b)*c]
    return [major,minor]

# TriRegiomontanusPointOfView returns point M on BC, for best view of a statue on BA such as AD/AB = k where D is base of statue
# distance BM is quadratic mean of BD and BA
def TriRegiomontanusPointOfView(a,b,c,k):
    l = c*sqrt(1 - k)
    t = [0, -c*l + a*c, b*l]
    return t

print "...Trilinear coordinates module loaded"









