import itertools

N = 12
init = [[i, i + 1] for i in range(1, N, 2)]
count = 0
def full(N):
    full_set = set([])
    full_set.clear()
    for i in range(1, N + 1):
        for j in range(i + 1, N + 1):
            for k in range(j + 1, N + 1):
                for l in range(k + 1, N + 1):
                    full_set.add((i, j, k, l))
    return full_set


def transp(n, k = 0):
    global count
    init = [[(i + k)%(2*n), (i + 1 + k)%(2*n)] for i in range(1, 2*n, 2)]

    for j in range(3):
        for i in range( n-1):
            init[i][0], init[i+1][0] = init[i+1][0], init[i][0]
            count += 1
            yield init

        init[-1] = list(reversed(init[-1]))

        for i in range(n-2,-1,-1):
            init[i][1], init[i+1][1] = init[i+1][1], init[i][1]
            count += 1
            yield init

        init[0] = list(reversed(init[0]))


def perm1(init, flag=False):
    if flag:
        # k = 1
        # o = -1
        # a = init[-1][k]
        # for i in range(len(init) - 1, 0, -1):
        #     init[i][k] = init[i + o][k]
        # init[0][k] = a
        k = 1
        o = 1
        a = init[0][k]
        for i in range(len(init) - 1):
            init[i][k] = init[i + o][k]
        init[-1][k] = a
    else:
        k = 0
        o = 1
        a = init[0][k]
        for i in range(len(init) - 1):
            init[i][k] = init[i + o][k]
        init[-1][k] = a
    return init

print()
def perm2(init, k=0, flag=False):
    # if flag:
    #     p=1
    init[-1 - k][0], init[-1 - k][1] = init[-1 - k][1], init[-1 - k][0]
    init = perm1(init,flag)
    init[-1 - k][0], init[-1 - k][1] = init[-1 - k][1], init[-1 - k][0]
    return init


def pr(init):
    a = [i[0] for i in init]
    b = [i[1] for i in init]
    print(a)
    print(b, "\n")


def get4prod(init):
    l = len(init)
    s = set([])
    s.clear()
    for i in range(l):
        for j in range(i + 1, l):
            s.add(tuple(sorted(init[i] + init[j])))
    return s


def strategy1(init, init1=None):
    l = len(init)
    s = set([])
    s.clear()
    print(N * (N - 1) * (N - 2) * (N - 3) / 24)
    print((N - 1) * (N // 2 - 1) * (N // 2) // 2)
    flag = False

    s = s.union(get4prod(init))
    for o in range(N//2 + 1):
        for j in range(N-1):
            for i in range(l - 1):
                init = perm2(init, k=i, flag=flag)
                s = s.union(get4prod(init))
                init = perm1(init, flag)
                s = s.union(get4prod(init))
            # for i in range(0,3,2):
            #     init[0][1], init[1][0] = init[1][0], init[0][1]
            # init[0][1], init[1][0] = init[1][0], init[0][1]
            # init = perm1(init)
            print(len(s))

        # flag = not flag
        print('--------')
        if o < N//2:
            init[0][1], init[o][0] = init[o][0], init[0][1]
        # init[-1], init[-2] = init[-2], init[-1]
    # print(full_set.difference(s))

    if init1 is not None:
        init = init1
        for j in range(1):
            for i in range(l - 1):
                init = perm2(init, k=i)
                s = s.union(get4prod(init))
                init = perm1(init)
                s = s.union(get4prod(init))

            print(len(s))
    return s



def noninter():
    s = set()
    for i in itertools.permutations(range(N)):
        s.add(tuple(sorted((tuple(sorted([i[j-1] + 1, i[j] + 1])) for j in range(1, N, 2)))))
    return s
# noninter()



if __name__ == "__main__":
    L_theor = []
    L_exp = []
    C = []
    Qubits = range(4, 31, 2)
    for N in Qubits:
        s = set()
        for k in range(N):
            for pair in transp(N//2, k):
                s = s.union(get4prod(pair))
        C.append(count)
        count = 0
        L_theor.append(N*(N-1)*(N-2)*(N-3)/6/4)
        print(N*(N-1)*(N-2)*(N-3)/6/4)
        L_exp.append(len(s))
        print(len(s))

import matplotlib.pyplot as plt

plt.plot([el//2 for el in Qubits], L_theor, color = "Black", lw = "3",)
plt.plot([el//2 for el in Qubits], L_exp, color = "Red", lw = "3", label="Number of prod 4 gamma")
plt.plot([el//2 for el in Qubits], C, color = "Purple", lw = "3",label="Number of transposition")
plt.ylabel('2 pauli operations',fontdict={'size': 22})
plt.xlabel('number of qubits',fontdict={'size': 22})
plt.grid()
plt.show()