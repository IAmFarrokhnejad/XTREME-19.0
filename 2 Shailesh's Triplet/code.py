import math

# Author: Morteza Farrokhnejad, Ali Farrokhnejad
cases = int(input())


for i in range(cases):
    n = int(input())
    
    if n % 2 == 1 or n&(n-1) == 0:
        print(-1)
        continue
    else:
        a = int(n - (1 << (n.bit_length() - 1))) 
        b = n // 2
        c = n // 2 + (1 << (n.bit_length() - 1)) 
        print(a, b, int(c))


        