from __future__ import annotations
from copy import deepcopy
"""
3
!..
#.?
##.


3
!#.
##.
..?
"""


class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def __eq__(self, point):
        return self.x == point.x and self.y == point.y
    def __add__(self, point: Point):
        return Point(self.x + point.x, self.y + point.y)
    
    def __repr__(self):
        return str(self.x) + " " + str(self.y)
    
directions = [
    Point(-1, 0),  
    Point(1, 0),   
    Point(0, -1),  
    Point(0, 1),   
    Point(-1, -1), 
    Point(-1, 1),  
    Point(1, -1),  
    Point(1, 1)    
]

class Labirinth:
    def __init__(self, A, init_point: Point, fin_point: Point, N: int):
        self.N = N 
        self.init_labirith = deepcopy(A)
        self.labirinth = deepcopy(A)# # -- препятсвие, i - величина пути, изначально -1
        self.fin_point = fin_point
        self.init_point = init_point
        self.current_points = [init_point]

    def propagate(self):
        while(len(self.current_points)):
            new_current_points = []
            for point in self.current_points:
                for dir in directions:
                    if (self[dir + point] != '#'):
                        if (self[dir + point] > self[point] + 1 or self[dir + point] == -1):
                            new_current_points.append(dir + point)
                            self[dir + point] = self[point] + 1
            self.current_points = new_current_points

    def get_reversed_path(self):
        if self[self.fin_point] in {"#", -1} :
            return "NO"
        point = self.fin_point
        while (point != self.init_point):
            for dir in directions:
                if (self[dir + point] != '#') and (self[dir + point] < self[point]):
                    self.init_labirith[point.x][point.y] = "+"
                    point = dir + point
                    break
        s = "YES\n"
        for i in range(self.N):
            for j in range(self.N):
                if self.init_labirith[i][j] == "#":
                    s += "#"
                if self.init_labirith[i][j] == -1:
                    s += "."
                if self.init_labirith[i][j] == "+":
                    s += "+"
                if self.init_labirith[i][j] == "0":
                    s += "!"
            s += "\n"
        return s
        


    def __getitem__(self, point: Point):
        if (0<= point.x < self.N) and (0 <= point.y < self.N):
            return self.labirinth[point.x][point.y]
        return '#'

    def __setitem__(self, point: Point, value: int):
        self.labirinth[point.x][point.y] = value

def read():
    N = int(input())
    A = []
    for i in range(N):
        s = input()
        A.append([])
        for j, char in enumerate(s):
            if char == ".":
                A[-1].append(-1)
            if char == "!":
                A[-1].append(0)
                init_point = Point(i,j)
            if char == "?":
                A[-1].append(-1)
                fin_point = Point(i,j)
            if char == "#":
                A[-1].append("#")
    return Labirinth(A, init_point, fin_point, N)
L = read()
L.propagate()
print(L.get_reversed_path())