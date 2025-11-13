#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <stdint.h>

#define N 1000000
#define MAXP 80000 
#define TH 100      
#define MAXQ 100000  

typedef long long ll;

static int bit[N + 2];
static inline void bit_add(int idx, int val) {
    while (idx <= N) {
        bit[idx] += val;
        idx += idx & -idx;
    }
}
static inline int bit_sum(int idx) {
    int s = 0;
    while (idx > 0) {
        s += bit[idx];
        idx -= idx & -idx;
    }
    return s;
}

static int primes[MAXP];
static int pc = 0;
static bool is_prime_arr[N + 1];

static void sieve_primes(void) {
    memset(is_prime_arr, 1, sizeof(is_prime_arr));
    is_prime_arr[0] = is_prime_arr[1] = 0;
    for (int i = 2; (long long)i * i <= N; ++i) {
        if (is_prime_arr[i]) {
            for (int j = i * i; j <= N; j += i) is_prime_arr[j] = 0;
        }
    }
    pc = 0;
    for (int i = 2; i <= N; ++i) if (is_prime_arr[i]) primes[pc++] = i;
}


static int *storage_vals = NULL;
static int storage_cap = 0;
static int storage_len = 0;
static int s_start[MAXP];
static int s_len[MAXP];

static void ensure_storage(int need) {
    if (storage_len + need <= storage_cap) return;
    int newcap = storage_cap ? storage_cap * 2 : 1 << 20;
    while (newcap < storage_len + need) newcap <<= 1;
    storage_vals = (int*)realloc(storage_vals, newcap * sizeof(int));
    storage_cap = newcap;
}

static inline int lower_bound_int(const int *a, int n, int x) {
    int l = 0, r = n;
    while (l < r) {
        int m = (l + r) >> 1;
        if (a[m] < x) l = m + 1;
        else r = m;
    }
    return l;
}
static inline int upper_bound_int(const int *a, int n, int x) {
    int l = 0, r = n;
    while (l < r) {
        int m = (l + r) >> 1;
        if (a[m] <= x) l = m + 1;
        else r = m;
    }
    return l;
}

typedef struct { int L, R, idx; } Query;


int main(void) {
    ios:;
    sieve_primes();


    storage_cap = 0; storage_len = 0; storage_vals = NULL;
    for (int i = 0; i < pc; ++i) {
        s_start[i] = storage_len;
        s_len[i] = 0;
        int p = primes[i];
        for (int j = 0; j < pc; ++j) {
            long long prod = (long long)p * primes[j];
            if (prod > N) break;
            ensure_storage(1);
            storage_vals[storage_len++] = (int)prod;
            s_len[i]++;
        }
    }

    int heavy_idx[MAXP], light_idx[MAXP];
    int hc = 0, lc = 0;
    for (int i = 0; i < pc; ++i) {
        if (s_len[i] > TH) heavy_idx[hc++] = i;
        else light_idx[lc++] = i;
    }

    int T;
    if (scanf("%d", &T) != 1) return 0;
    Query queries[MAXQ];
    int qcount = T;
    for (int i = 0; i < T; ++i) {
        int L, R;
        scanf("%d %d", &L, &R);
        if (L < 1) L = 1;
        if (R > N) R = N;
        queries[i].L = L;
        queries[i].R = R;
        queries[i].idx = i;
    }
    int *U = (int*)malloc((2 * T + 5) * sizeof(int));
    int Ulen = 0;
    for (int i = 0; i < T; ++i) {
        U[Ulen++] = queries[i].R;
        int Lm1 = queries[i].L - 1;
        if (Lm1 < 0) Lm1 = 0;
        U[Ulen++] = Lm1;
    }
    // unique sort U
    qsort(U, Ulen, sizeof(int), (int(*)(const void*,const void*)) (int(*)(const int*, const int*)) strcmp); 
    free(U);
    U = (int*)malloc((2 * T + 5) * sizeof(int));
    Ulen = 0;
    for (int i = 0; i < T; ++i) {
        U[Ulen++] = queries[i].R;
        int Lm1 = queries[i].L - 1;
        if (Lm1 < 0) Lm1 = 0;
        U[Ulen++] = Lm1;
    }
    // integer comparator:
    int cmp_int(const void *a, const void *b) {
        int va = *(const int*)a, vb = *(const int*)b;
        if (va < vb) return -1;
        if (va > vb) return 1;
        return 0;
    }
    qsort(U, Ulen, sizeof(int), cmp_int);
    int w = 0;
    for (int i = 0; i < Ulen; ++i) if (i == 0 || U[i] != U[i-1]) U[w++] = U[i];
    Ulen = w;

    int *Rq_idx = (int*)malloc(T * sizeof(int));
    int *Lm1_idx = (int*)malloc(T * sizeof(int));
    for (int i = 0; i < T; ++i) {
        int R = queries[i].R;
        int Lm1 = queries[i].L - 1; if (Lm1 < 0) Lm1 = 0;
        Rq_idx[i] = lower_bound_int(U, Ulen, R);
        Lm1_idx[i] = lower_bound_int(U, Ulen, Lm1);
    }

    uint32_t **heavy_pref = NULL;
    if (hc > 0) {
        heavy_pref = (uint32_t**)malloc(hc * sizeof(uint32_t*));
        for (int h = 0; h < hc; ++h) {
            heavy_pref[h] = (uint32_t*)malloc((size_t)Ulen * sizeof(uint32_t));
            int pi = heavy_idx[h];
            int *vals = &storage_vals[s_start[pi]];
            int vlen = s_len[pi];
            int ptr = 0;
            uint32_t cnt = 0;
            for (int k = 0; k < Ulen; ++k) {
                int bound = U[k];
                while (ptr < vlen && vals[ptr] <= bound) { cnt++; ptr++; }
                heavy_pref[h][k] = cnt;
            }
        }
    }
    int *cntA = (int*)calloc((N + 1), sizeof(int));
    long long total_pairs = 0;

    for (int ii = 0; ii < lc; ++ii) {
        int pi = light_idx[ii];
        int *vals = &storage_vals[s_start[pi]];
        int vlen = s_len[pi];
        for (int i = 0; i < vlen; ++i) {
            int a = vals[i];
            for (int j = i + 1; j < vlen; ++j) {
                (void)vals[j];
                cntA[a]++;
                total_pairs++;
            }
        }
    }

    int *offset = (int*)malloc((N + 2) * sizeof(int));
    offset[0] = 0;
    for (int x = 1; x <= N + 1; ++x) offset[x] = offset[x-1] + cntA[x-1];
    int LP = offset[N+1-1];
    if (LP != (int)total_pairs) {
        LP = 0;
        for (int x = 0; x <= N; ++x) LP += cntA[x];
    }
    int *flatBs = (int*)malloc(((size_t)LP) * sizeof(int));
    for (int x = 0; x <= N; ++x) cntA[x] = 0;

    for (int ii = 0; ii < lc; ++ii) {
        int pi = light_idx[ii];
        int *vals = &storage_vals[s_start[pi]];
        int vlen = s_len[pi];
        for (int i = 0; i < vlen; ++i) {
            int a = vals[i];
            for (int j = i + 1; j < vlen; ++j) {
                int b = vals[j];
                int pos = offset[a] + cntA[a];
                flatBs[pos] = b;
                cntA[a]++;
            }
        }
    }

    Query *qs_sorted = (Query*)malloc(T * sizeof(Query));
    memcpy(qs_sorted, queries, T * sizeof(Query));
    int cmpQ(const void *pa, const void *pb) {
        const Query *A = pa, *B = pb;
        return B->L - A->L; // descending
    }
    qsort(qs_sorted, T, sizeof(Query), cmpQ);

    memset(bit, 0, sizeof(bit));
    ll *answers = (ll*)malloc(T * sizeof(ll));
    int curA = N + 1;
    for (int qi = 0; qi < T; ++qi) {
        int L = qs_sorted[qi].L;
        int R = qs_sorted[qi].R;
        if (L < 1) L = 1;
        if (R > N) R = N;
        for (int a = curA - 1; a >= L; --a) {
            int cnt = cntA[a];
            if (cnt) {
                int start = offset[a];
                for (int t = 0; t < cnt; ++t) {
                    int b = flatBs[start + t];
                    bit_add(b, 1);
                }
            }
        }
        curA = L;
        int light_contrib = bit_sum(R);

        int orig_idx = qs_sorted[qi].idx;
        int Rid = Rq_idx[orig_idx];
        int Lmid = Lm1_idx[orig_idx];
        ll heavy_contrib = 0;
        for (int h = 0; h < hc; ++h) {
            uint32_t cntR = heavy_pref[h][Rid];
            uint32_t cntL = (Lmid >= 0) ? heavy_pref[h][Lmid] : 0;
            uint32_t cnt = cntR - cntL;
            if (cnt > 1) heavy_contrib += (ll)cnt * (cnt - 1) / 2;
        }
        answers[orig_idx] = (ll)light_contrib + heavy_contrib;
    }

    for (int i = 0; i < T; ++i) {
        printf("%lld\n", answers[i]);
    }

    /* cleanup */
    free(storage_vals);
    free(U);
    free(Rq_idx);
    free(Lm1_idx);
    for (int h = 0; h < hc; ++h) free(heavy_pref[h]);
    free(heavy_pref);
    free(cntA);
    free(offset);
    free(flatBs);
    free(qs_sorted);
    free(answers);

    return 0;
}

// Authors: Morteza Farrokhnejad, Ali Farrokhnejad