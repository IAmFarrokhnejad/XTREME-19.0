#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <stdint.h>


// Authors: Morteza Farrokhnejad, Ali Farrokhnejad


#define NEG_INF_LL (LLONG_MIN / 4)

typedef struct {
    long long *m;
    long long *b;
    int sz;
    int cap;
    int ptr;
} Hull;

static void hull_init(Hull *h, int cap) {
    h->m = (long long*)malloc(sizeof(long long) * cap);
    h->b = (long long*)malloc(sizeof(long long) * cap);
    h->sz = 0;
    h->cap = cap;
    h->ptr = 0;
}

static void hull_free(Hull *h) {
    free(h->m);
    free(h->b);
    h->m = h->b = NULL;
    h->sz = h->cap = h->ptr = 0;
}

static int hull_bad(Hull *h, long long m1, long long b1, long long m2, long long b2, long long m3, long long b3) {
    __int128 left = (__int128)(b2 - b1) * (__int128)(m2 - m3);
    __int128 right = (__int128)(b3 - b2) * (__int128)(m1 - m2);
    return left >= right;
}

static void hull_add(Hull *h, long long m, long long b) {
    if (h->sz >= h->cap) {
        int newcap = h->cap * 2;
        if (newcap < h->cap + 4) newcap = h->cap + 4;
        h->m = (long long*)realloc(h->m, sizeof(long long) * newcap);
        h->b = (long long*)realloc(h->b, sizeof(long long) * newcap);
        h->cap = newcap;
    }
    while (h->sz >= 2) {
        long long m1 = h->m[h->sz - 2];
        long long b1 = h->b[h->sz - 2];
        long long m2 = h->m[h->sz - 1];
        long long b2 = h->b[h->sz - 1];
        if (hull_bad(h, m1, b1, m2, b2, m, b)) {
            h->sz--;
            if (h->ptr > h->sz - 1) h->ptr = h->sz - 1;
        } else break;
    }
    h->m[h->sz] = m;
    h->b[h->sz] = b;
    h->sz++;
}

static long long hull_query(Hull *h, long long x) {
    if (h->sz == 0) return NEG_INF_LL;
    while (h->ptr + 1 < h->sz) {
        __int128 v1 = (__int128)h->m[h->ptr] * x + (__int128)h->b[h->ptr];
        __int128 v2 = (__int128)h->m[h->ptr + 1] * x + (__int128)h->b[h->ptr + 1];
        if (v2 >= v1) h->ptr++;
        else break;
    }
    __int128 res = (__int128)h->m[h->ptr] * x + (__int128)h->b[h->ptr];
    if (res > LLONG_MAX) return LLONG_MAX;
    if (res < LLONG_MIN) return LLONG_MIN;
    return (long long)res;
}

int main() {
    int N, M;
    if (scanf("%d", &N) != 1) return 0;
    char *alpha = (char*)malloc((N + 5));
    if (scanf("%s", alpha) != 1) return 0;
    if (scanf("%d", &M) != 1) return 0;
    char *beta = (char*)malloc((M + 5));
    if (scanf("%s", beta) != 1) return 0;
    long long match, mismatch, gap;
    if (scanf("%lld %lld %lld", &match, &mismatch, &gap) != 3) return 0;

    int foo = N;
    int bar = M;
    int W = bar + 1;
    long long *dp = (long long*)malloc(sizeof(long long) * (foo + 1) * (bar + 1));
    long long *holeX = (long long*)malloc(sizeof(long long) * (foo + 1) * (bar + 1));
    long long *holeY = (long long*)malloc(sizeof(long long) * (foo + 1) * (bar + 1));
    if (!dp || !holeX || !holeY) {
        return 0;
    }

    for (int i = 0; i <= foo; ++i) {
        for (int j = 0; j <= bar; ++j) {
            int idx = i * W + j;
            dp[idx] = NEG_INF_LL;
            holeX[idx] = NEG_INF_LL;
            holeY[idx] = NEG_INF_LL;
        }
    }

    dp[0 * W + 0] = 0;
    for (int i = 1; i <= foo; ++i) {
        __int128 t = i;
        __int128 cost128 = t * t * (__int128)gap;
        long long cost = (long long)cost128;
        dp[i * W + 0] = cost;
        holeX[i * W + 0] = cost;
    }
    for (int j = 1; j <= bar; ++j) {
        __int128 t = j;
        __int128 cost128 = t * t * (__int128)gap;
        long long cost = (long long)cost128;
        dp[0 * W + j] = cost;
        holeY[0 * W + j] = cost;
    }

    Hull *colHull = (Hull*)malloc(sizeof(Hull) * (bar + 1));
    for (int j = 0; j <= bar; ++j) {
        hull_init(&colHull[j], 4);
    }

    for (int i = 1; i <= foo; ++i) {
        int t = i - 1;
        for (int j = 0; j <= bar; ++j) {
            int idx = t * W + j;
            long long best = dp[idx];
            if (holeX[idx] > best) best = holeX[idx];
            __int128 mm = (__int128)(-2) * (__int128)gap * (__int128)t;
            __int128 bb = (__int128)best + (__int128)gap * (__int128)t * (__int128)t;
            long long mll = (long long)mm;
            long long bll = (long long)bb;
            hull_add(&colHull[j], mll, bll);
        }

        Hull rowHull;
        hull_init(&rowHull, 4);

        for (int j = 1; j <= bar; ++j) {
            int tcol = j - 1;
            int idx_t = i * W + tcol;
            long long bestRow = dp[idx_t];
            if (holeY[idx_t] > bestRow) bestRow = holeY[idx_t];
            __int128 mmr = (__int128)(-2) * (__int128)gap * (__int128)tcol;
            __int128 bbr = (__int128)bestRow + (__int128)gap * (__int128)tcol * (__int128)tcol;
            hull_add(&rowHull, (long long)mmr, (long long)bbr);

            int idx_diag = (i - 1) * W + (j - 1);
            long long bestPrev = dp[idx_diag];
            if (holeX[idx_diag] > bestPrev) bestPrev = holeX[idx_diag];
            if (holeY[idx_diag] > bestPrev) bestPrev = holeY[idx_diag];
            long long matchscore = (alpha[i - 1] == beta[j - 1]) ? match : mismatch;
            dp[i * W + j] = bestPrev + matchscore;

            long long qcol = hull_query(&colHull[j], i);
            if (qcol == NEG_INF_LL) {
                holeX[i * W + j] = NEG_INF_LL;
            } else {
                __int128 add = (__int128)gap * (__int128)i * (__int128)i;
                __int128 res = (__int128)qcol + add;
                if (res > LLONG_MAX) holeX[i * W + j] = LLONG_MAX;
                else if (res < LLONG_MIN) holeX[i * W + j] = LLONG_MIN;
                else holeX[i * W + j] = (long long)res;
            }

            long long qrow = hull_query(&rowHull, j);
            if (qrow == NEG_INF_LL) {
                holeY[i * W + j] = NEG_INF_LL;
            } else {
                __int128 add = (__int128)gap * (__int128)j * (__int128)j;
                __int128 res = (__int128)qrow + add;
                if (res > LLONG_MAX) holeY[i * W + j] = LLONG_MAX;
                else if (res < LLONG_MIN) holeY[i * W + j] = LLONG_MIN;
                else holeY[i * W + j] = (long long)res;
            }
        }

        hull_free(&rowHull);
    }

    long long finalScore = dp[foo * W + bar];
    if (holeX[foo * W + bar] > finalScore) finalScore = holeX[foo * W + bar];
    if (holeY[foo * W + bar] > finalScore) finalScore = holeY[foo * W + bar];

    printf("%lld\n", finalScore);


    for (int j = 0; j <= bar; ++j) hull_free(&colHull[j]);
    free(colHull);
    free(dp);
    free(holeX);
    free(holeY);
    free(alpha);
    free(beta);

    return 0;
}