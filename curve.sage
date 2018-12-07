print('s|t')
for s in xrange(1,100):
    for t in xrange(s):
        if (t^2-3*t+3)%s==0:
            print(s,t)

