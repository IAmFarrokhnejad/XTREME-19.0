import sys
sys.setrecursionlimit(1000000)
input = sys.stdin.readline
# Author: Morteza Farrokhnejad, Ali Farrokhnejad

n = int(input())
distances = list(map(int, input().split()))
length = 2 * n

for d in distances:
    if d < 0 or d >= length:
        print("No")
        sys.exit(0)

distances.sort()

def backtrack(p, result, openS, used_mask):
    if p == length:
        return used_mask == (1 << n) - 1 and not openS

    # Try '('
    result[p] = '('
    openS.append(p)
    if backtrack(p + 1, result, openS, used_mask):
        return True
    openS.pop()

    # Try ')'
    if openS:
        openP = openS.pop()
        distance = p - openP - 1

       
        for i in range(n):
            if not (used_mask >> i) & 1 and distances[i] == distance:
                new_mask = used_mask | (1 << i)
                result[p] = ')'
                if backtrack(p + 1, result, openS, new_mask):
                    return True
                break

        openS.append(openP)

    return False


result = [''] * length
if backtrack(0, result, [], 0):
    print("Yes")
    print(''.join(result))
else:
    print("No")