#
# mutates a sequence according to the Jukes-Cantor model
# X is an numpy.array with entries in {0,1,2,3}
# t is branch length
# mu is mutation rate
#


def mutate(X, t, mu):
    import numpy.random as rand
    L = len(X)
    mutatedSeq = X.copy()

    numMutation = rand.poisson(L * mu * t)
    for i in range(numMutation):
        site = rand.randint(0, L)
        mutatedSeq[site] = rand.randint(0, 4)

    return mutatedSeq
