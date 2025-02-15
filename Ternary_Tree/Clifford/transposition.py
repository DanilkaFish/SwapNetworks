
n = 14
wires = list(range(n))
n = len(wires)
P = []
P = [i for i in wires]
T = []
F = []
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
while (i < n):
    F.append(T[i])
    i += 2

if i == n+1:
    i = n - 2
else:
    i = n - 1
while (i >0):
    F.append(T[i])
    i -= 2 


dic = {(i,j) : 0 for i in range(n) for j in range(i+1,n)}
# for el in [(i,j,k) for i in range(n) for j in range(i+1, n) for k in range(j+1,n)]:
    # dic[el] = 0

for el in [(i,j,k, r) for i in range(n) for j in range(i+1, n) for k in range(j+1,n) for r in range(k+1, n)]:
    dic[el] = 0

def adswap(wires,i, j):
    wires[i], wires[j] = wires[j], wires[i] 

def layer_two_swap(wires, S, k=0):

    size = len(wires)
    new_wires = [None]*len(wires)
    for i in range(size):
        new_wires[i] = wires[i]

    # print(new_wires)
    for i in range(0, (size - k%2)//2):
        adswap(new_wires, S[2*i + k%2], S[2*i + 1 + k%2])
    # print(new_wires)
    # for i in range(1, (size-1)//2 ):
    #     adswap(new_wires, S[1+2*i], S[2+2*i])
    return new_wires

def layer_swap(wires, num_neig=4, j=0):

    size = len(wires)
    new_wires = [None]*len(wires)
    for i in range(size):
        new_wires[i] = wires[i]

    # print(new_wires)
    l = num_neig
    for i in range(0, (size - l*(j%2))//l//2):
        for k in range(l):
            adswap(new_wires, 2*l*i + k + l*(j%2), 2*l*i + l + k + l*(j%2))
    
    # print(new_wires)
    # for i in range(0, (size - 3)//6):
    #     for k in range(3):
    #         adswap(new_wires, 6*i + k + 3, 6*i + 6 + k)

    return new_wires

def layer_three_swap_(wires):
    size = len(wires)
    new_wires = [None]*len(wires)
    for i in range(size):
        new_wires[i] = wires[i]

    for i in range(0, size//3):
        adswap(new_wires, 3*i, 3*i + 1)

    for i in range(0, (size-1)//3):
        adswap(new_wires, 3*i + 1, 3*i + 2)

    for i in range(0, (size-2)//3):
        adswap(new_wires, 3*i + 2, 3*i + 3)
    return new_wires

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


def PSN(wires):
    for i in range(n//2):
        add_ad_pairs(wires)
        wires = layer_swap(wires, 1, 0)
        # add_ad_pairs(wires)
        wires = layer_swap(wires, 1, 1)
    return list(reversed(wires))

def TSN(wires):
    # print("TSN : ", wires)
    for k in range(n//2):
        wires = PSN(wires)
        wires = layer_two_swap(wires, T, 0)
        # wires = PSN(wires)
        wires = layer_two_swap(wires, T, 1)
        # wires = PSN(wires)
    for k in range(n//2):
        wires = layer_two_swap(wires, T, 0)
        wires = layer_two_swap(wires, T, 1)
    return wires

# print(dic)
print("TSN : ", wires)
# for k in range(n//2):

#     wires = TSN(wires)
#     wires = layer_two_swap(wires, F, 0)
#     wires = TSN(wires)
#     wires = layer_two_swap(wires, F, 1)
# print(wires)

print("TSN : ", wires)
# wires = TSN(wires)
for k in range(n//2):

    wires = TSN(wires)
    wires = layer_two_swap(wires, F, 0)
    # wires = TSN(wires)
    wires = layer_two_swap(wires, F, 1)
print(wires)

for pair in dic:
    if dic[pair] == 0:
        print(pair, " : ", dic[pair])

# print(T)
# print(P)