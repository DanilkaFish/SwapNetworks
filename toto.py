import numpy as np 
def _gen_mask(self):
    size = 1 << self.num_qubits
    ne = self.num_electrons
    indices = np.array([(i,j) for i in range(size) for j in range(size) if i.bit_count() != ne or j.bit_count() != ne])
    flat_indices = indices[:, 0] * size + indices[:, 1]
    return flat_indices

if __name__ == "__main__":
    size = 1 << 2
    i = 8
    r = i.bit_count()
    ones = np.ones((size,size))
    indices = np.array([(i,j) for i in range(size) for j in range(size) if i.bit_count() != 1 or j.bit_count() != 1])
    flat_indices = indices[:, 0] * size + indices[:, 1]

    print(ones)
    # for index in index_array:
    #     print(index, ones[index])
    #     ones[*index] = 0
    ones.flat[flat_indices] = 0
    print(ones)