#include <bits/stdc++.h>

// Fenwick: https://cp-algorithms.com/data_structures/fenwick.html
// Largest sum sub: https://www.geeksforgeeks.org/dsa/largest-sum-subarray-least-k-numbers/
// K-th largest sum: https://www.geeksforgeeks.org/dsa/k-th-largest-sum-contiguous-subarray/

// Author: Morteza Farrokhnejad, Ali Farrokhnejad
using namespace std;

#define MAXLIMIT 100010

long long S[25];

int A[MAXLIMIT], P[MAXLIMIT], B[MAXLIMIT];

long long fen_count[MAXLIMIT], fen_sumj[MAXLIMIT];

int N, K;

void update(long long fen[], int idx, long long val) {
    while (idx <= N + 1) {
        fen[idx] += val;
        idx += idx & -idx;
    }
}

long long query(long long fen[], int idx) {
    long long res = 0;
    while (idx > 0) {
        res += fen[idx];
        idx -= idx & -idx;
    }
    return res;
}

int main() {
    scanf("%d %d", &N, &K);
    for (int i = 1; i <= N; i++) {
        scanf("%d", &A[i]);
    }
    memset(S, 0, sizeof(S));
    for (int t = 1; t <= 20; t++) {
        int M = 1 << t;
        for (int i = 1; i <= N; i++) {
            B[i] = (A[i] >= M ? 1 : 0);
        }
        P[0] = 0;
        for (int i = 1; i <= N; i++) {
            P[i] = P[i - 1] + B[i];
        }
        memset(fen_count, 0, sizeof(fen_count));
        memset(fen_sumj, 0, sizeof(fen_sumj));
        update(fen_count, 1, 1);
        update(fen_sumj, 1, 0);
        for (int i = 1; i <= N; i++) {
            for (int l = 1; l <= K; l++) {
                int thresh = P[i] - l;
                if (thresh >= 0) {
                    long long cnt = query(fen_count, thresh + 1);
                    long long sj = query(fen_sumj, thresh + 1);
                    long long contrib = (long long)i * cnt - sj;
                    S[l] += contrib;
                }
            }
            int key = P[i] + 1;
            update(fen_count, key, 1);
            update(fen_sumj, key, (long long)i);
        }
    }
    for (int l = 1; l <= K; l++) {
        if (l > 1) printf(" ");
        printf("%lld", S[l]);
    }
    printf("\n");
    return 0;
}