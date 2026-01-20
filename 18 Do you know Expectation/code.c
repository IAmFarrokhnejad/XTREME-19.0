#include <stdio.h>
#include <math.h>

// Author: Morteza Farrokhnejad, Ali Farrokhnejad

int main() {
    int n, k;
    const int maximum = 1024;
    scanf("%d %d", &n, &k);
    
    int a[500000];
    for (int i = 0; i < n; i++) {
        scanf("%d", &a[i]);
    }
    double dp[maximum];
    double newDp[maximum];
    for (int i = 0; i < maximum; i++) {
        dp[i] = 0.0;
    }
    
    
    dp[0] = 1.0;
    for (int i = 0; i < n; i++) {
        int num = a[i];
        for (int j = 0; j < maximum; j++) 
        {
            newDp[j] = 0.0;
        }
        
        for (int xorValue = 0; xorValue < maximum; xorValue++) {
            if (dp[xorValue] > 0) {
                newDp[xorValue] += dp[xorValue] * 0.5;
                int newXor = xorValue ^ num;
                newDp[newXor] += dp[xorValue] * 0.5;
            }
        }
        
        for (int j = 0; j < maximum; j++) {
            dp[j] = newDp[j];
        }
    }
    
    double expected = 0.0;
    for (int xorValue = 0; xorValue < maximum; xorValue++) {
        if (dp[xorValue] > 0) {
            double power = pow((double)xorValue, (double)k);
            expected += power * dp[xorValue];
        }
    }
    printf("%.2f\n", expected);
    return 0;
}