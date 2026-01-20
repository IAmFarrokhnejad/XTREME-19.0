#include <bits/stdc++.h>
using namespace std;
// Author: Morteza Farrokhnejad, Ali Farrokhnejad

const int MAXX = 1 << 21;
const int N = 7;
const int E = 21;

long long dp[MAXX];
bool is_valid[MAXX];

int edge_id[N+1][N+1];
pair<int, int> edges[E];

int parent[N+1];
int rankk[N+1];

int find(int x) {
    if (parent[x] != x) parent[x] = find(parent[x]);
    return parent[x];
}

void union_sets(int x, int y) {
    x = find(x);
    y = find(y);
    if (x != y) {
        if (rankk[x] < rankk[y]) swap(x, y);
        parent[y] = x;
        if (rankk[x] == rankk[y]) rankk[x]++;
    }
}

void precompute() {
    int idx = 0;
    memset(edge_id, -1, sizeof(edge_id));
    for (int i = 1; i <= N; ++i) {
        for (int j = i + 1; j <= N; ++j) {
            edges[idx] = {i, j};
            edge_id[i][j] = idx;
            edge_id[j][i] = idx;
            idx++;
        }
    }

    is_valid[0] = false;
    for (int mask = 1; mask < MAXX; ++mask) {
        // Init UF
        for (int i = 1; i <= N; ++i) {
            parent[i] = i;
            rankk[i] = 0;
        }
        int deg[N+1] = {0};
        for (int bit = 0; bit < E; ++bit) {
            if (mask & (1 << bit)) {
                int a = edges[bit].first;
                int b = edges[bit].second;
                deg[a]++;
                deg[b]++;
                union_sets(a, b);
            }
        }
        int odd = 0;
        bool has_root[N+1] = {false};
        for (int v = 1; v <= N; ++v) {
            if (deg[v] > 0) {
                has_root[find(v)] = true;
            }
            odd += (deg[v] & 1);
        }
        int comps = 0;
        for (int i = 1; i <= N; ++i) {
            if (has_root[i]) comps++;
        }
        is_valid[mask] = (comps == 1 && (odd == 0 || odd == 2));
    }

    memset(dp, 0, sizeof(dp));
    for (int mask = 0; mask < MAXX; ++mask) {
        dp[mask] = is_valid[mask] ? 1LL : 0LL;
    }

    for (int bit = 0; bit < E; ++bit) {
        for (int mask = 0; mask < MAXX; ++mask) {
            if (mask & (1 << bit)) {
                dp[mask] += dp[mask ^ (1 << bit)];
            }
        }
    }
}

int main() {
    precompute();
    int T;
    scanf("%d", &T);
    for (int t = 0; t < T; ++t) {
        int M;
        scanf("%d", &M);
        int S = 0;
        for (int i = 0; i < M; ++i) {
            int x, y;
            scanf("%d %d", &x, &y);
            if (x > y) swap(x, y);
            int id = edge_id[x][y];
            S |= (1 << id);
        }
        printf("%lld\n", dp[S]);
    }
    return 0;
}