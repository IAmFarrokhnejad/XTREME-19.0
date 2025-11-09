import sys

# Author: Morteza Farrokhnejad, Ali Farrokhnejad


input = sys.stdin.read
data = input().split()
index = 0
T = int(data[index])
index += 1
for _ in range(T):
    N = int(data[index])
    K = int(data[index + 1])
    S = data[index + 2]
    index += 3
    a = [1 if c == 'S' else 0 for c in S]
    M = max(0, N - K + 1)
    prefix_xor = [0]
    count = 0
    impossible = False
    for j in range(M):
        left = max(0, j - K + 1)
        sum_prev = prefix_xor[j] ^ prefix_xor[left]
        x = a[j] ^ sum_prev
        if x == 1:
            count += 1
        new_prefix = prefix_xor[j] ^ x
        prefix_xor.append(new_prefix)
    if not impossible:
        for j in range(M, N):
            left = max(0, j - K + 1)
            right = M - 1
            if left > right:
                current_xor = 0
            else:
                current_xor = prefix_xor[right + 1] ^ prefix_xor[left]
            if current_xor != a[j]:
                impossible = True
                break
    if impossible:
        print(-1)
    else:
        print(count)