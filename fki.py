
def fill_distances(matrix):
    n = len(matrix)
    m = len(matrix[0])
    
    output = [[-1 for _ in range(m)] for _ in range(n)]
    
    filled = []
    
    for i in range(n):
        for j in range(m):
            if matrix[i][j] == 1:
                output[i][j] = 0
                filled.append((i, j))
    
    dirs = [(-1, 0), (1, 0), (0, -1), (0, 1),
                  (-1, -1), (-1, 1), (1, -1), (1, 1)]
    
    while filled:
        x, y = filled.pop(0)
        for dx, dy in dirs:
            nx = x + dx
            ny = y + dy
            if 0 <= nx < n and 0 <= ny < m and output[nx][ny] == -1:
                output[nx][ny] = output[x][y] + 1
                filled.append((nx, ny))
    return output

N, M = [int(el) for el in input().split()]

matrix = [[int(el) for el in input().split()] for i in range(N)]

result = fill_distances(matrix)
for row in result:
    print(' '.join(map(str, row)))