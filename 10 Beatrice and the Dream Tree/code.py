from collections import deque, defaultdict

def uniqueNodes(N,adjacent):
    if N ==0:
        return 0
        
    parent =[0] * (N + 1)
    children =[[] for _ in range(N + 1)]
    visitedN =[False] * (N + 1)
    q = deque([1])
    visitedN[1] =True
    while q:
        u = q.popleft()
        for v in adjacent[u]:
            if not visitedN[v]:
                visitedN[v] = True
                parent[v]= u
                children[u].append(v)
                q.append(v)
    
    pending =[0] * (N + 1)
    for i in range(1, N + 1):
        pending[i]=len(children[i])
    processQ = deque()
    for i in range(1,N+ 1):
        if pending[i] == 0:
            processQ.append(i)
    subtreeId =[0] * (N + 1)
    idMap={}
    counter= 0
    while processQ:
        u = processQ.popleft()
        childIDs= [subtreeId[v] for v in children[u]]
        childIDs.sort()
        key= tuple(childIDs)
        if key not in idMap:
            counter +=1
            idMap[key] =counter
        subtreeId[u] =idMap[key]
        
        if u == 1:
            continue
        p = parent[u]
        pending[p] -=1
        if pending[p] == 0:
            processQ.append(p)
    
    unique =[False] * (N + 1)
    unique[1] =True
    q =deque([1])
    
    while q:
        u= q.popleft()
        if not unique[u]:
            continue
        groups = defaultdict(int)
        for v in children[u]:
            groups[subtreeId[v]] += 1
        for v in children[u]:
            if groups[subtreeId[v]] == 1:
                unique[v] =True
                q.append(v)
    
    count= sum(1 for i in range(1, N + 1) if unique[i])
    return count

while True:
    try:
        line = input().strip()
        N= int(line)
        adjacent =[[] for _ in range(N + 1)]
        for _ in range(N - 1):
            line= input().strip()
            u, v = map(int, line.split())
            adjacent[u].append(v)
            adjacent[v].append(u)
        result = uniqueNodes(N, adjacent)
        print(result)
    except EOFError:
        break
    except ValueError:
        break

    # Authors: Morteza Farrokhenjad, Ali Farrokhnejad