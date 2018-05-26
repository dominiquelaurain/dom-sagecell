print "Arithmetic (AC) module loading ..."

# General usage function

# CalkinWilfNext returns next rational number in the Calkin Wilf sequence
def CalkinWilfNext(q):
    q = 1 / (2*floor(q) - q + 1)
    return q

# CalkinWilfLevelFirst returns first rational number in the Calkin Wilf sequence at depth n
# because that index of that number is 2^n which can be written 10...0  (with n bits one)
# the rational has continued fraction representation [0,n-1,1] and tequals 0 + 1 / (n - 1 + 1/1) = 1/n
def CalkinWilfLevelFirst(n):
    return 1/n

# NormalisedAdditiveFactor returns normalised factors of integer n as bits list [b0,b1,..,bs] of 2-exponents such as
# n = sum( (-1)^i 2^bi) for i from 0 to s
def NormalisedAdditiveFactor(n,addfact):
    m = 0
    while n != 0:
        # dividing by 2 until odd
        k = 0 ; u = 1
        while n % 2 == 0:
            k = k + 1 ; u = 2*u
            n = n / 2
        # add exponent to list
        m = m + k ; addfact.append(m)
        # decrase n by 1 (because odd) and negate
        n = - (n - 1)

# CalkinWilfElementCff returns the continued fraction for m-th term in Calkin-Wilf sequence
# OBSOLETE (too much complicated)
#def CalkinWilfElementCff(m):
#    n = floor(log(m,2)) # n is depth in the tree
#    j = m + 1 - 2^n #  j is index at the depth  (1,2,..)
#    addfact = [] ; NormalisedAdditiveFactor(j,addfact)
#    # computation of the continued fraction from the additive factors of j
#    cf = [] ; b = 0
#    for bi in addfact:
#        cf.append(bi - b)
#        b = bi
#    cf.append(n + 1 - b)
#    return cf

# CalkinWilfElementCff returns the continued fraction for m-th term in Calkin-Wilf sequence
# the continued fraction is deduced from binary representation of m
# it starts with 0 with m is even, then count numbers of consecutive bits with same value
def CalkinWilfElementCff(m):
    m1 = m ; cf = [] ; b = m % 2 ; k = 0
    if b == 0:
        cf.append(0)  # add 0 if even number in order to have odd number of elements for cff
    while m1 > 0:
        bn = m1 % 2
        if (bn != b):
            cf.append(k)
            b = bn ; k = 0
        k = k + 1
        m1 = (m1 - bn) / 2
    cf.append(k)
    return cf

# CalkinWilfElement returns rational m-th term sequence in Calkin-Wilf sequence
def CalkinWilfElement(m):
    cf = CalkinWilfElementCff(m)
    return continued_fraction(cf).value()

# HyperbinaryExpansionsNumber returns number of hyperbinary expansions for integer m
# m = a + 2b + 4c + 8d + 16e + 32f + ... where coefficients can be 0,1,2
# the function returns number of ways writing m as sum of at most two same power of two
def HyperbinaryExpansionsNumber(m):
    q = CalkinWilfElement(m)
    return q.denominator()

# EvenSizeRationalCff adds an element to the list if odd number of elements for continued fraction ql for rational
# Using a trick that  an = (an - 1) + 1/1  so ql = [a0,a1,..,an] = [a0,a1,..,(an - 1),1]
def EvenSizeRationalCff(ql):
    n = len(ql)
    if is_odd(n) and n > 0:
        ql[n-1] = ql[n-1] - 1
        ql.append(1)

# CalkinWilfIndex returns depth n (1,2,..) and index j (1,2,..) for a rational q in the Calkin-Wilf tree
# there are 1 + 2 + 2^(n-1) = 2^n - 1 rational vertices in levels 1,2,..,n-1
def CalkinWilfIndex(q):
    ql = continued_fraction_list(q)
    EvenSizeRationalCff(ql) # we need even number of elements in the continued fraction
    n = sum(ql)
    # deducing j from the fact that rational j-th vertex in the n-th level is given by continued fraction
    #  [bk, (bk-1 - bk),...,(b0 - b1),(n - b0)]
    # it's the inverse of the code in the CalkinWilfElement function : we use the continued fraction to retrieve
    # back the bj additive factors
    j = 0
    for ai in reversed(ql):
        for d in xrange(ai):
            j = 2*j
        j = j + 1
        j = -j
    j = j + 1
    return [n,j]

# SternBrocotMap returns rational element at the same place (depth,index) than q in the Calkin-Wilf tree
def SternBrocotMap(q):
    ql = continued_fraction_list(q)
    EvenSizeRationalCff(ql)
    # Add 1 to first term and substract 1 to last term
    ql[0] += 1 ; ql[len(ql)-1] -= 1
    # Rewrite backwards
    ql = list(reversed(ql))
    # Deduce a rational number and inverse it
    r = 1/continued_fraction(ql).value()
    return r

print "...Arithmetic module loaded"









