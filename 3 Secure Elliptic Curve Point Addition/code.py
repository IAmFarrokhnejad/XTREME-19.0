import math
# Authors: Morteza Farrokhnejad, Ali Farrokhnejad
# https://stackoverflow.com/questions/4798654/modular-multiplicative-inverse-function-in-python
def mod_inverse(x, p):
    return pow(x, -1, p)

cases = int(input())

for i in range(cases):
    params = list(map(int, input().split()))
    a = params[0]
    b = params[1]
    p = params[2]
    x1 = params[3]
    y1 = params[4]
    x2 = params[5]
    y2 = params[6]
    
    if x1 == x2:
        if y1 != y2:
            print("POINT_AT_INFINITY")
            continue
        elif x1 == x2 and y1 == y2:
            if y1 == 0:
                print("POINT_AT_INFINITY")
                continue
            else: 
                slope = ((3 * x1 * x1 + a) * mod_inverse((2 * y1) % p, p)) % p
    else:
        slope = ((y2 - y1) * mod_inverse((x2-x1) % p, p)) % p
        
    x3 = (slope * slope - x1 - x2) % p
    r3 = (slope * x3 + (y1 - slope * x1)) % p
    y3 = (-r3) % p
    
    x3 = (x3 + p) % p
    y3 = (y3 + p) % p
    
    
    print(f"{x3} {y3}")