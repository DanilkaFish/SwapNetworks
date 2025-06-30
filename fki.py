from copy import deepcopy

def next_number(reversed_number_array: list):
    flag = True
    for i in range(len(reversed_number_array) - 1):
        if (reversed_number_array[i] < reversed_number_array[i + 1] - 1):
            reversed_number_array[i] += 1
            flag = False
    if flag:
        if reversed_number_array[-1] < 9:
            reversed_number_array[-1] += 1
        else:
            reversed_number_array.append(0)
            k = 0
            for i in range(len(reversed_number_array)):
                reversed_number_array[i] = k 
                k += 1
    
N = int(input())
n_arr = [0]
for i in range(N):
    next_number(n_arr)

res = 0
prod = 1
for num in n_arr:
    res += num * prod
    prod *= 10
print(res)
    
    
