
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
// Author: Morteza Farrokhnejad, Ali Farrokhnejad

int main() {
    int N, k;
    scanf("%d %d", &N, &k);
    
    int* cost =(int*)calloc(N, sizeof(int));
    int* cs = (int*)malloc(N * sizeof(int));
    int cCount =0;
    
    for (int i =1; i < N; i++) {
        scanf("%d", &cost[i]);
        
        bool already_exists = false;
        for (int j= 0; j < cCount; j++) {
            if (cs[j] == cost[i]) {
                already_exists = true;
                break;
            }
        }
        if (!already_exists) {
            cs[cCount++]= cost[i];
        }
    }
    
    const int MAX_S = 5000;
    bool* possible = (bool*)calloc(MAX_S + 1, sizeof(bool));
    possible[0] = true;
    
    for (int i = 0; i < cCount; i++) {
        int c = cs[i];
        for (int j = c; j <= MAX_S; j++) {
            possible[j] = possible[j] || possible[j - c];
        }
    }
    int Q;
    scanf("%d", &Q);
    
    for (int q = 0; q < Q; q++) {
        int A, B;
        scanf("%d %d", &A, &B);
        
        int dist = cost[A] + cost[B];
        
        if (k < dist) {
            printf("No\n");
            continue;
        }
        
        int r = k - dist;
        if (r % 2 != 0) {
            printf("No\n");
            continue;
        }
        
        int s = r / 2;
        
        if (possible[s]) {
            printf("Yes\n");} else {
            printf("No\n");
        }}
    
    free(cs);
    free(cost);
    free(possible);
    
    return 0;
}