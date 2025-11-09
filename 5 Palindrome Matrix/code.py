nm = input().split()
n = int(nm[0])
m = int(nm[1])
# Authors: Morteza Farrokhnejad, Ali Farrokhnejad
grid = []
for i in range(n):
    nums = list(input())
    for j in range(m):
        if nums[j] == ".":
            grid.append(None)
        else:
            grid.append(int(nums[j]))

par = [i for i in range(n*m)]
rank = [0 for i in range(n*m)]

for i in range(n):
    j = 0
    while j < m:
        if grid[i * m + j] is None:
            j += 1
            continue
        start = j
        while j < m and grid[i * m + j] is not None:
            j += 1
        end = j - 1
        for k in range(int((end - start + 1) / 2)):
            a = i * m + start + k
            b = i * m + end - k
            
            root_a = a
            while par[root_a] != root_a:
                root_a = par[root_a]
            cur = a
            while par[cur] != root_a:
                nxt = par[cur]
                par[cur] = root_a
                cur = nxt
            
            root_b = b
            while par[root_b] != root_b:
                root_b = par[root_b]
            cur = b
            while par[cur] != root_b:
                nxt = par[cur]
                par[cur] = root_b
                cur = nxt
            if root_a != root_b:
                if rank[root_a] < rank[root_b]:
                    root_a, root_b = root_b, root_a
                par[root_b] = root_a
                if rank[root_a] == rank[root_b]:
                    rank[root_a] += 1

for j in range(m):
    i = 0
    while i < n:
        if grid[i*m + j] is None:
            i += 1
            continue
        start = i
        while i < n and grid[i * m + j] is not None:
            i += 1
        end = i - 1
        for k in range(int((end - start + 1) / 2)):
            a = (start + k) * m + j
            b = (end - k) * m + j 
            
            root_a = a
            while par[root_a] != root_a:
                root_a = par[root_a]
            cur = a
            while par[cur] != root_a:
                nxt = par[cur]
                par[cur] = root_a
                cur = nxt
            
            root_b = b
            while par[root_b] != root_b:
                root_b = par[root_b]
            cur = b
            while par[cur] != root_b:
                nxt = par[cur]
                par[cur] = root_b
                cur = nxt
            if root_a != root_b:
                if rank[root_a] < rank[root_b]:
                    root_a, root_b = root_b, root_a
                par[root_b] = root_a
                if rank[root_a] == rank[root_b]:
                    rank[root_a] += 1

sgmnts = {}
for i in range(n):
    for j in range(m):
        if grid[i * m + j] is not None:
            ind = i * m + j
            root = ind
            while par[root] != root:
                root = par[root]
            cur = ind
            while par[cur] != root:
                nxt = par[cur]
                par[cur] = root
                cur = nxt
            if root not in sgmnts:
                sgmnts[root] = []
            sgmnts[root].append((i, j, grid[i * m + j]))

res = ['.' if x is None else str(x) for x in grid]
for root, cells in sgmnts.items():
    mat = [c for i, j, c in cells]
    fin_cost = float('inf')
    digit = -1
    
    for d in range(10):
        cost = sum(abs(c - d) for c in mat)
        if cost < fin_cost or (cost == fin_cost and d < digit):
            fin_cost = cost
            digit = d
    for i, j, c in cells:
        res[i * m + j] = str(digit)

for i in range(n):
    print("".join(res[i * m : (i + 1) * m]))