
n = 22
wires = list(range(n))
P  = list(range(n))
dic = {}
# dic = {(i,j) : 0 for i in range(n) for j in range(i+1,n)}
# for el in [(i,j,k) for i in range(n) for j in range(i+1, n) for k in range(j+1,n)]:
#     dic[el] = 0

for el in [(i,j,k, r) for i in range(n) for j in range(i+1, n) for k in range(j+1,n) for r in range(k+1, n)]:
    dic[el] = 0


# def generatePairs(wires):
    
def gen_ad_pair(wires):
    s = set()
    for i in range(n-1):
        s.add(tuple(sorted([wires[i], wires[i+1]])))
    return s

def gen_ad_triples(wires):
    s = set()
    for i in range(n-2):
        s.add(tuple(sorted([wires[i], wires[i + 1], wires[i + 2]])))
    return s

def gen_ad_fours(wires):
    s = set()
    for i in range(n-3):
        s.add(tuple(sorted([wires[i], wires[i + 1], wires[i + 2], wires[i+3]])))
    return s

def add_ad_pairs(wires):
    s = gen_ad_pair(wires)
    for pair in dic:
        if pair in s:
            dic[pair] += 1
            # for _pair in s:
            #     if (_pair[0] != pair[1]) and (pair[0] != _pair[1]):
            #         dic[pair].add(_pair)

    s = gen_ad_triples(wires)
    for toto in dic:
        if toto in s:
            dic[toto] += 1

    s = gen_ad_fours(wires)
    for toto in dic:
        if toto in s:
            dic[toto] += 1

def to_T(P):
    T = []
    i=0
    while (i < n):
        T.append(P[i])
        i += 2

    if i == n+1:
        i = n - 2
    else:
        i = n - 1
    while (i >0):
        T.append(P[i])
        i -= 2 
    return T
T = []
i=0
while (i < n):
    T.append(wires[i])
    i += 2

if i == n+1:
    i = n - 2
else:
    i = n - 1
while (i >0):
    T.append(wires[i])
    i -= 2 

i=0

def adswap(wires, i, j):
    wires[i], wires[j] = wires[j], wires[i]     


def layer_swap_1(wires, init=0, step=0, S=P, through=False):
    new_wires = [i for i in wires]
    i = init
    if not through:
        while (i < n - 1):
            adswap(new_wires, S[i], S[i + 1])
            i += 2 + step
    else:
        while (i < n + j - 2):
            adswap(new_wires, S[i%n], S[(i + 1)%n])
            i += 2 + step
    return new_wires


def layer_swap(wires, num_neig=4, j=0, S=P, step=0, inter=False, displ=0, inv=False):

    size = len(wires)
    new_wires = [None]*len(wires)
    for i in range(size):
        new_wires[i] = wires[i]

    l = num_neig

    if inv:
        if not inter:
            for i in range(0, (size - l*(j%2))//l//2):
                for k in range(l):
                    adswap(new_wires, S[2*l*i + k + l*(j%2)], S[2*l*i + l + k + l*(j%2)])
                for k in range(l//2):
                    adswap(new_wires, S[2*l*i + l + k  + l*(j%2)], S[2*l*i + l + l - k - 1 + l*(j%2)])
        else:
            for i in range(0, (size)//l//2):
                for k in range(l):
                    adswap(new_wires, S[(2*l*i + k + l*(j%2)) % size ], S[(2*l*i + l + l - k - 1 + l*(j%2)) % size])
                for k in range(l//2):
                    adswap(new_wires, S[2*l*i + l + k  + l*(j%2)], S[2*l*i + l + l - k - 1 + l*(j%2)])
    else:
        if not inter:
            for i in range(0, (size - l*(j%2))//l//2):
                for k in range(l):
                    adswap(new_wires, S[2*l*i + k + l*(j%2)], S[2*l*i + l + k + l*(j%2)])
        else:
            for i in range(0, (size)//l//2):
                for k in range(l):
                    adswap(new_wires, S[(2*l*i + k + l*(j%2)) % size ], S[(2*l*i + l + k + l*(j%2)) % size])
    return new_wires

def PSN(wires):
    for i in range(n//2):
        add_ad_pairs(wires)
        wires = layer_swap(wires, 1, 0)
        add_ad_pairs(wires)
        wires = layer_swap(wires, 1, 1)
    return list(reversed(wires))


print(wires)
# print(to_T(wires))

# for i in range(n//2):
T = to_T(P)
copy_wires = [i for i in wires]

# wires = PSN(wires)
for j in range(0,n):
    for i in range(n//2):
        wires = PSN(wires)
        wires = layer_swap_1(wires, init=0, step=0, S=T, through=False)
        wires = layer_swap_1(wires, init=1, step=0, S=T, through=False)
        # wires = layer_swap_1(wires, j % 3, step=0, S=to_T(T))
        # if all(copy_wires[i]== wires[i] for i in range(n)):
        #     print(j)
    wires = layer_swap_1(wires, 0 % 3, step=1, S=T, through=True)
    # wires = layer_swap_1(wires, (j+1) % 3, step=1, S=T)
    # wires = layer_swap_1(wires, (j+2) % 3, step=1, S=T)
i = 0
j = 0
for pair in dic:
    if dic[pair] == 0:
        i += 1
        print(pair, " : ", dic[pair])
    j += dic[pair]
print(i)
print(j)
print(len(dic))
