n=50
# Choose random prime r of right size
r = random_prime(2^(n+1), lbound=2^n)
h = 2
s = 3
assert(4*h-s>0)
for D in xrange(1,10000):
    vest = isqrt((5*r-1)/D)
    for vdelta in xrange(-10000,10000):
        V = vest + vdelta
        t = D*V^2-5*r+1
        if t%s==0:
            n = h*r
            q = n+t-1
            if 4*q - t^2 >= 0 and is_prime(q):
                print('q:',q,'r:',r,'t:',t,'D:',D,'V:',V)
