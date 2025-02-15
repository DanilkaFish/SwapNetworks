
n = 20
wires = list(range(n))
P  = list(range(n))

dic = {(i,j) : 0 for i in range(n) for j in range(i+1,n)}
for el in [(i,j,k) for i in range(n) for j in range(i+1, n) for k in range(j+1,n)]:
    dic[el] = 0

# for el in [(i,j,k, r) for i in range(n) for j in range(i+1, n) for k in range(j+1,n) for r in range(k+1, n)]:
#     dic[el] = 0

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

def layer_two_swap(wires, S=T):
    new_wires = [i for i in wires]
    i = 0
    while (i<n-1):
        adswap(new_wires, S[i], S[((i + 1) % n)])
        i += 2

    i = 1
    while (i<n-1):
        adswap(new_wires, S[i], S[((i + 1) % n)])
        i += 2
    print(new_wires)
    return new_wires

    
def layer_via_one(wires, init1, init2, S=T):
    i = 0
    # init2 = init2 % (n//2)
    # init1 = init1 % (n//2)
    # print("before", wires)
    while (i < n):
        if init1*2 <= i <= init1*2 + 1:
            pass
        else:
            adswap(wires, S[i % n], S[(i + 1) % n])
        i += 2
    wires = PSN(wires)
    i = 1
    while (i < n):
        if (init2*2  + 1<= i <= init2*2 + 2) :
            pass
        else:
            adswap(wires, S[i % n], S[(i + 1) % n])
        i += 2
    print("after", wires)
    return wires

def layer_swap(wires, num_neig=4, j=0, S=P, inter=False, displ=0, inv=False):

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

def layer_pp_swap(wires, i, j, S=T):
    i = i % 2
    if i == 0:
        i = j
        while (i < n - 3 - j):
            adswap(wires, S[i % n], S[(i + 2) % n])
            adswap(wires, S[(i + 1) % n], S[(i + 3) % n])
            i += 4
    else:
        i = 2 + j
        while (i < n - 3 - j):
            adswap(wires, S[i % n], S[(i + 2) % n])
            adswap(wires, S[(i + 1) % n], S[(i + 3) % n])
            i += 4
        adswap(wires, S[(-1) % n - j], S[(-2) % n - j])
        
    return wires 

def PSN(wires):
    for i in range(n//2):
        add_ad_pairs(wires)
        wires = layer_swap(wires, 1, 0)
        wires = layer_swap(wires, 1, 1)
    return list(reversed(wires))


print(wires)
print(to_T(wires))

# for i in range(n//2):
copy_wires = [i for i in wires]
for j in range(0, n//2):
# for j in range(n//4 - 1):
#     wires = layer_two_swap(wires)
#     wires = PSN(wires)
    wires = PSN(wires)
    # wires = layer_via_one(wires, j,  (i ) % (n//2))
    # wires = layer_pp_swap(wires, j, j//(n//2))
    wires = layer_swap(wires, 1, 0, S=T, inter=False)
    wires = PSN(wires)
    wires = layer_swap(wires, 1, 1, S=T, inter=False, inv=True)
    # wires = copy_wires
    # # print(to_T(wires))
    # wires = layer_swap(wires, 1, 0, S=T)
    # wires = layer_swap(wires, 1, 1, S=T)
    # print(to_T(wires))

    # wires = PSN(wires)
    # wires = layer_swap(wires, 1, 0, S=T)
    # wires = layer_swap(wires, 1, 1, S=T)
    # wires = layer_swap(wires, 1, 1)
    # wires = list(range(n))
    # wires = layer_two_swap(wires)
    # wires = layer_via_one(wires, (i) % (n//2),  (i) % (n//2))
    # wires = layer_swap(wires, 1, 0)
    # wires = layer_swap(wires, 1, 1)
    # wires = layer_via_one(wires, n,  j % (n//2))
    # wires = layer_via_one(wires, (n//2 - 1), (n//2 - 1))
print(wires)
i = 0
for pair in dic:
    if dic[pair] != 0:
        i += 1
        print(pair, " : ", dic[pair])
print(i)