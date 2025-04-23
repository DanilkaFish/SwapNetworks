from sympy import sqf


class fermions:
    def __init__(self, n: int):
        self.n = n
        self.list = list(range(2*n))
        self.fours = {} 
        for i in range(1):
            for j in range(i + 1, 2*self.n):
                for k in range(j + 1, 2*self.n):
                    for l in range(k + 1, 2*self.n):
                        self.fours[(i,j,k,l)] = 0
                        
    def get_left(self, i : int):
        return self.list[i]

    def get_right(self, i : int):
        return self.list[- i + 2 * self.n - 1]
    
    def move(self, i):
        self.list[i % (self.n * 2)], self.list[(i + 1) % (self.n*2)] = \
            self.list[(i + 1) % (self.n*2)], self.list[i % (self.n * 2)]

        
    def generate_pairs(self):
        pairs = [(self.get_left(i), self.get_right(i)) for i in range(self.n)]
        return pairs 
    
    def generate_fours(self):
        pairs = self.generate_pairs()

        for i in range(self.n):
            for j in range(i + 1, self.n):
                four = tuple(sorted(pairs[i] + pairs[j]))
                if four in self.fours:
                    self.fours[four] = self.fours.get(four, 0) + 1
        return self.fours
    
    def __str__(self):
        s = ""
        pairs = self.generate_pairs()
        for i in range(self.n):
            s += str(pairs[i]) + "\n"
        return s
    
if __name__ == "__main__":
    n = 6
    ferms = fermions(n)
    
    i = 0
    for _ in range(n + 1):
    # for _ in range(n*2*n - 2*n - n + 1 ):
        print(ferms)
        ferms.move(i)
        i += 1
        if (i  + 1) % n == 0:
            print(ferms)
            ferms.move(i)
            i += 1
        fours = ferms.generate_fours()
    print(ferms)
    # for key, num in fours.items():
    #     if num != 0:
    #         print(key, " : ", num)