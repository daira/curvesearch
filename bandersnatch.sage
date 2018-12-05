def f(x):
    return 2*x^2 + 1

prevx = 1
prediction = 1
u = 0
i = 0
while i < 10:
    if u % 2 == 1:
        x2 = 2*(u^4 + u^2)
    else:
        x2 = (u^4 + 2*u^2)/2

    if squarefree_part(x2) == 1:
        x = isqrt(x2)
        print("* x = %d (pred. %d) = 2*%d, f(x) = %d = 8*%d + 1, x/prevx = %f" % (x, prediction, x//2, f(x), (f(x)-1)//8, 1.0*x/prevx))
        assert squarefree_part(f(x)) == 1
        assert i < 5 or squarefree_part(f(prediction)) == 1
        prediction = x*x//prevx
        prevx = x
        DV2 = 12*x^4 + 16*x^3 + 12*x^2 + 8*x + 3
        assert f(x).divides(DV2)
        R = DV2//f(x)
        D = squarefree_part(R)
        print("%d^2 * %s (D = 10^%f) for %f bits" % (isqrt(f(x)), factor(R), log(D)/log(10), log(36*x^4)/log(2)))
        i += 1

    u += 1

while True:
    x = prediction
    print("# x = %d = 2*%d, f(x) = %d = 8*%d + 1, x/prevx = %f" % (x, x//2, f(x), (f(x)-1)//8, 1.0*x/prevx))
    assert squarefree_part(f(x)) == 1
    prediction = x*x//prevx
    prevx = x
    DV2 = 12*x^4 + 16*x^3 + 12*x^2 + 8*x + 3
    assert f(x).divides(DV2)
    R = DV2//f(x)
#    D = squarefree_part(R)
    print("%d^2 * %s (R = 10^%f) for %f bits" % (isqrt(f(x)), R, log(R)/log(10), log(36*x^4)/log(2)))
