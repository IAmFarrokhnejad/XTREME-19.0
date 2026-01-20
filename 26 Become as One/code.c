#include <stdio.h>

int main() {
    
    
    int N, Q;
    scanf("%d %d", &N, &Q);
    int A[N];
    for (int i = 0; i < N; i++) 
    {
        scanf("%d", &A[i]);
        
    }
    
    long long pfx[N + 1];
    pfx[0] = 0;
    for (int i = 1; i <= N; i++) {
        pfx[i] = pfx[i - 1] + (1LL << A[i - 1]);
    }
    for (int q = 0; q < Q; q++) {
        int L, R;
        scanf("%d %d", &L, &R);
        long long s = pfx[R] - pfx[L - 1];
        if ((s & (s - 1)) == 0) {
            printf("Yes\n");
        } else {
            printf("No\n");
        }}
    return 0;
}

// Author: Morteza Farrokhnejad, Ali Farrokhnejad