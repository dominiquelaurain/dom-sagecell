print "Matrix representation of Conics (MR) module loading ..."

# Conics is represented by quadratic form q = A x^2 + B xy + C y^2 + D x + E y + F = 0
# with can be written as x^T Mq x = 0 where Mq is matrix
#             (A   B/2  D/2)
#     Mq  =   (B/2  C   E/2)
#             (D/2 E/2   F )

# We use here symbolic ring SR for computations.
# Values are corced to real only for ^(1/n) operations

# General function definining conic and its inverse
def MrcConic(q):
    [A,B,C,D,E,F] = q
    mq = matrix(SR,3,3)
    mq[0,0] =  A  ; mq[0,1] = B/2 ; mq[0,2] = D/2
    mq[1,0] = B/2 ; mq[1,1] =  C  ; mq[1,2] = E/2
    mq[2,0] = D/2 ; mq[2,1] = E/2 ; mq[2,2] =  F
    return mq

def MrcNullConic():
    return MrcConic([0,0,0,0,0,0])

def MrcIdentity():
    return MrcConic([1,0,1,0,0,1])

def MrcTrace(mq):
    return mq[0,0] + mq[1,1] + mq[2,2]

def MrcConicParameters(mq):
    A = mq[0,0] ; B = 2*mq[0,1] ; C = mq[1,1] ; D = 2*mq[0,2] ; E = 2*mq[1,2] ; F = mq[2,2]
    return [A,B,C,D,E,F]

def MrcConicEquation(mq,x,y):
    [A,B,C,D,E,F] = MrcConicParameters(mq)
    x,y = var('x,y')
    return A*x^2 + B*x*y + C*y^2 + D*x + E*y + F

# MrcConicSectionFunctions returns functions (A,B,C,D,E,F) of conic section defined by  A x^2 + B xy + C y^2 + D x + E y + F = 0
# it is defined from five points on the conic, three points out of five must not be aligned
#  x1^2 A + x1y1 B + y1^2 C + x1 D + y1 E = -F
#  x2^2 A + x2y2 B + y2^2 C + x2 D + y2 E = -F
#  x3^2 A + x3y3 B + y3^2 C + x3 D + y3 E = -F
#  x4^2 A + x4y4 B + y4^2 C + x4 D + y4 E = -F
#  x5^2 A + x5y5 B + y5^2 C + x5 D + y5 E = -F
def MrcConicSectionFunctions(p1,p2,p3,p4,p5):
    [x1,y1] = p1 ; [x2,y2] = p2 ; [x3,y3] = p3 ; [x4,y4] = p4 ; [x5,y5] = p5
    x,y,A,B,C,D,E,F = var('x,y,A,B,C,D,E,F') ; f(x,y) = x^2*A + x*y*B + y^2*C + x*D + y*E + F
    # OLD res = solve([f(x1,y1),f(x2,y2),f(x3,y3),f(x4,y4),f(x5,y5),F == -1],A,B,C,D,E,F) ... but NOK because we can have F == 0
    res = solve([f(x1,y1),f(x2,y2),f(x3,y3),f(x4,y4),f(x5,y5),A + B + C + D + E + F == 1],A,B,C,D,E,F)
    A = res[0][0].rhs()
    B = res[0][1].rhs()
    C = res[0][2].rhs()
    D = res[0][3].rhs()
    E = res[0][4].rhs()
    F = res[0][5].rhs()
    return [A,B,C,D,E,F]


# MrcFunctionsTruncate truncate to 0 too small real values
def MrcFunctionsTruncate(fn):
    for i in xrange(5):
        if abs(fn[i]) < 0.0000001:
            fn[i] = 0
    return fn

# MrcIsDegenerate returns 0 if the conic is degenerate
def MrcIsDegenerate(mq):
    return mq.determinant()

# MrcMinor returns type of conic from determinant of minor matrix reomiving last row and last column
def MrcType(mq):
    m = matrix(SR,2,2)
    m[0,0] = mq[0,0] ; m[0,1] = mq[0,1]
    m[1,0] = mq[1,0] ; m[1,1] = mq[1,1]
    return m.determinant()

#  MrcIsHyperbola,  MrcIsHyperbola and  MrcIsHyperbola use MrcType to check type of conic
def MrcIsHyperbola(mq):
    return sign(MrcType(mq)) + 1

def MrcIsParabola(mq):
    return sign(MrcType(mq))

def MrcIsEllipse(mq):
    return sign(MrcType(mq)) - 1

# MrcIsCircle check directly if conic is circle
def MrcIsCircle(mq):
    quad = ((mq[0,0]-mq[1,1])^2 + mq[0,1]^2 + mq[1,0]^2)
    return quad

# MrcIsParallelLines check if conic is made of two parallel lines
def MrcIsParallelLines(mq):
    return MrcIsDegenerate(mq)^2 + MrcIsParabola(mq)^2

# MrcIsIntersectingLines check if conic is made of two intersecting lines
def MrcIsIntersectingLines(mq):
    return MrcIsDegenerate(mq)^2 + MrcIsHyperbola(mq)^2

# MrcIsSinglePoint check if conic is made of a single point
def MrcIsSinglePoint(mq):
    return MrcIsDegenerate(mq)^2 + MrcIsEllipse(mq)^2

# MrcConicCenter returns center point of conic, if not parabola (B^2 = 4AC)
def MrcCenter(mq):
    [A,B,C,D,E,F] = MrcConicParameters(mq)
    xc = (B*E-2*C*D)/(4*A*C-B^2)
    yc = (D*B-2*A*E)/(4*A*C-B^2)
    return [xc,yc]

# MrcConicEccentricy returns eccentricity of conics
def MrcEccentricity(mq):
    [A,B,C,D,E,F] = MrcConicParameters(mq)
    nu = -sign(mq.determinant())
    ro = sqrt((A - C)^2 + B^2)
    e2 = 2*ro/(nu*(A + C)  + ro)
    e = sqrt(e2)
    return e

# MrcEllipsisSemiAxisQuadrances returns quadrances (=square of length) for semi-axis
def MrcEllipsisSemiAxisQuadrances(mq):
    [A,B,C,D,E,F] = MrcConicParameters(mq)
    k = sqrt((A - C)^2 + B^2)
    l = A + C
    m = B^2 - 4*A*C
    n = 2*(A*E^2  + C*D^2 + F*B^2 - B*D*E - 4*A*C*F)
    qa = n/(m*(k - l)) ; qb = n/(m*(- k - l))
    return [qa,qb]

# MrcEllipsisSlopeCotangent returns cotangent of 2*phi where phi is angle of semi-axis with x-axis
def MrcEllipsisSlopeCotangent(mq):
    [A,B,C,D,E,F] = MrcConicParameters(mq)
    ct = (A - C)/B
    return ct


# MrcCharacteristicPolynomial returns characteristic polynomial : det(M -xI)
def MrcCharacteristicPolynomial(x,mq):
    mc = matrix(SR,3,3)
    mc[0,0] = mq[0,0] - x ; mc[0,1] = mq[0,1] ; mc[0,2] = mq[0,2]
    mc[1,0] = mq[1,0] ; mc[1,1] = mq[1,1] - x ; mc[1,2] = mq[1,2]
    mc[2,0] = mq[2,0] ; mc[2,1] = mq[2,1] ; mc[2,2] = mq[2,2] - x
    return mc.determinant()

# MrcSkewSymmetricMatrix returns the skew symmetrix matrix associated to a 3D vector
def MrcSkewSymmetricMatrix(v):
    ms = matrix(SR,3,3)
    ms[0,0] = 0     ; ms[0,1] = -v[2] ; ms[0,2] = v[1]
    ms[1,0] = v[2]  ; ms[1,1] =  0    ; ms[1,2] = -v[0]
    ms[2,0] = -v[1] ; ms[2,1] =  v[0] ; ms[2,2] =  0
    return ms

# MrcMaxElement returns indices of element with max absolute value
def MrcMaxElement(M):
    #print "MaxElement M = ",M
    maxvalue = 0
    im = 0 ; jm = 0
    for i in xrange(3):
        for j in xrange(3):
            if abs(M[i][j]) > maxvalue :
                maxvalue = abs(M[i][j])
                im = i ; jm = j
    return [im,jm]

# MrcTruncate reset to 0 smallest elements
def MrcTruncate(M):
    M1 = M
    for i in xrange(3):
        for j in xrange(3):
            M1[i,j] = M1[i,j].n() # eval to real
            if abs(M1[i,j]) < 0.0000001:
                M1[i,j] = 0.0 # truncate to 0.0
    return M1

# MrcGetLines returns lines l,m from conic defined by non symmetrix matrix rank 1 if l and m not parallel
# and rank 2 if l and m parallel
def MrcGetLines(mq1):
    #print "Rank of matrix : ",mq1.rank()
    # Choose non vanishing value as element of maximum absolute value
    [i,j] = MrcMaxElement(mq1)
    l = mq1[i,:] ; m = mq1[:,j].transpose()
    return [l,m]

# MrcDegenerateConic2Lines returns l,m from degenerate conic mq
def MrcDegenerateConic2Lines(mq):
    #print "Rank of matrix mq : ",mq.rank()
    # Compute adjoint matrix
    mqa = mq.adjoint() #; print "Adjoint matrix : \n",mqa
    # Choose max absolute value element
    [i,j] = MrcMaxElement(mqa) ; v = mqa[i,j] ; s = sign(v) ; v = sqrt(s*v)
    # Compute column vector and normalize it (dividing by square root of max absolute value element)
    p1 = mqa[:,j] #; print "vector p1 : ",p1
    p = [s*p1[0,0]/v,s*p1[1,0]/v,s*p1[2,0]/v] #; print "p = ",p
    # Construct a "skew" symmetric matrix S lxm
    Slxm = MrcSkewSymmetricMatrix(p) #; print "Skew Symmetrix matrix : \n",Slxm
    # Add to degenerate conic
    mq1 = mq + Slxm #; print "mq1 matrix = \n",mq1
    [l,m] = MrcGetLines(mq1) #; print "l = ",l, " ; m = ",m
    return [l,m]

# MrcEigenValuesSymmetric returns eigenvalues of symmetrix matrix mq, roots x of det(xI - mq) = 0
def MrcEigenValuesSymmetric(mq):
    p1 = mq[0,1]^2 + mq[0,2]^2 + mq[1,2]^2
    if (p1 == 0):
        # mq is diagonal
        eig1 = mq(0,9)
        eig2 = mq(1,1)
        eig3 = mq(2,2)
    else:
        k = MrcTrace(mq)/3
        p2 = (mq[0,0] - k)^2 + (mq[1,1] - k)^2 + (mq[2,2] - k)^2 + 2 * p1
        p = sqrt(p2 / 6)
        mq1 = (1/p)*(mq - k*MrcIdentity())
        r = mq1.determinant()/2
        # In exact arithmetic for a symmetric matrix  -1 <= r <= 1
        # but computation error can leave it slightly outside this range.
        if (r <= -1):
            phi = pi / 3
        else:
            if (r >= 1):
                phi = 0
            else:
                phi = acos(r)/3
    #
    # the eigenvalues satisfy eig3 <= eig2 <= eig1
    eig1 = k + 2*p*cos(phi)
    eig3 = k + 2*p*cos(phi + (2*pi/3))
    eig2 = 3*k - eig1 - eig3     # since trace(A) = eig1 + eig2 + eig3
    #
    return [eig1,eig2,eig3]

# MrcRealThirdRoot(x) returns real root of y^3 = x, without the weird y = (x)^(1/3) bad behaviour when x < 0
def MrcRealThirdRoot(x):
    if x < 0:
        y = - (-x)^(1/3)
    else:
        y = (x)^(1/3)
    return y

# MrcEigenValues returns eigenvalues of matrix mq, roots x of det(xI - mq) = 0
def MrcEigenValues(mq):
    f(x) = mq.charpoly('x') #; print "f : ",f
    cf = f.list() #; print " cf = ",cf
    a = cf[3] ; b = cf[2] ; c = cf[1] ; d = cf[0] # equation ax^3 + bx^2 + cx + d = 0
    Q = (3*a*c - b^2)/(9*a^2)
    R = (9*a*b*c - 27*a^2*d - 2*b^3)/(54*a^3)
    U = Q^3 + R^2 ; V = 1
    if U < 0:
            U = -U ; V = I
    S = MrcRealThirdRoot(R + V*sqrt(U))
    T = MrcRealThirdRoot(R - V*sqrt(U))
    eig1 = S + T -b/(3*a)
    eig2 = -(S + T)/2 - b/(3*a) + I*(sqrt(3)/2)*(S - T)
    eig3 = -(S + T)/2 - b/(3*a) - I*(sqrt(3)/2)*(S - T)
    #
    return [eig1,eig2,eig3]

# MrcIntersectConics returns intersections points of two conics
def MrcIntersectConics(mq1,mq2):
    # Compute null conic
    mnull = MrcNullConic()
    # Compute matrix ma = mq1 (-mq2)^-1
    ma = mq1*(-mq2).inverse() #; print "mq = mq1 (-mq2)^-1 matrix = \n",ma
    # Eigenvalues are roots of characteristic polynomial det(ma - xI) or root of det(mq1 + x mq2)
    #print "ma = ",ma
    #e = ma.eigenvalues()
    e = MrcEigenValues(ma)
    #print "e = ",e
    #
    # Annoying but sometimes one need to get EXACT values for symbolic computation and sometimes
    # one need next code for APPROXIMATE real values
    for i in xrange(3):
        e[i] = e[i].n()    # don't limit precision too much, to stay near EXACT values
    #print "e = ",e
    #
    # Select last real value
    kr = e[0]
    for k in e:
        #print "k = ",k
        # Check not null conic
        if  ((mq1 + k*mq2).adjoint() == mnull):
            print "null conic for k = \n",k
        else:
            if k.is_real():
                kr = k
    #print "kr = ",kr, " ; type(kr) = ",type(kr)
    # Compute matrix of the degenerate conic (determinant of matrix is zero)
    mq = (mq1 + kr*mq2) #; print "mq matrix = \n",mq
    # Compute lines intersecting
    [l,m] = MrcDegenerateConic2Lines(mq)
    # Intersecting lines with first conic to get points
    x,y = var('x,y')
    lf(x,y) = (l*matrix(SR,3,1,[x,y,1])).determinant() #; print "lf : ",lf
    mf(x,y) = (m*matrix(SR,3,1,[x,y,1])).determinant() #; print "mf : ",mf
    f1(x,y) = (matrix(SR,1,3,[x,y,1])*mq1*matrix(SR,3,1,[x,y,1])).determinant()  #; print "f1 : ",f1
    res = solve([f1 == 0,(lf*mf) == 0],x,y)
    return res


# MrcExactIntersectConics returns intersections points of two conics
def MrcExactIntersectConics(mq1,mq2):
    # Compute null conic
    mnull = MrcNullConic()
    # Compute matrix ma = mq1 (-mq2)^-1
    ma = mq1*(-mq2).inverse() ; print "mq = mq1 (-mq2)^-1 matrix = \n",ma
    # Eigenvalues are roots of characteristic polynomial det(ma - xI) or root of det(mq1 + x mq2)
    #print "ma = ",ma
    #e = ma.eigenvalues()
    e = MrcEigenValues(ma)
    #
    # Select last real value
    kr = e[0]
    for k in e:
        #print "k = ",k
        # Check not null conic
        if  ((mq1 + k*mq2).adjoint() == mnull):
            print "null conic for k = \n",k
        else:
            kr = k
    #print "kr = ",kr, " ; type(kr) = ",type(kr)
    # Compute matrix of the degenerate conic (determinant of matrix is zero)
    mq = (mq1 + kr*mq2) ; print "mq matrix = \n",mq
    # Compute lines intersecting
    [l,m] = MrcDegenerateConic2Lines(mq)
    # Intersecting lines with first conic to get points
    x,y = var('x,y')
    lf(x,y) = (l*matrix(SR,3,1,[x,y,1])).determinant() #; print "lf : ",lf
    mf(x,y) = (m*matrix(SR,3,1,[x,y,1])).determinant() #; print "mf : ",mf
    f1(x,y) = (matrix(SR,1,3,[x,y,1])*mq1*matrix(SR,3,1,[x,y,1])).determinant()  #; print "f1 : ",f1
    res = solve([f1 == 0,(lf*mf) == 0],x,y)
    return res

print "... Matrix representation of Conics module loaded"









