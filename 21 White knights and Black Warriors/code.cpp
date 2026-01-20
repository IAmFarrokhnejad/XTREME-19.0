// Segment Tree: https://www.geeksforgeeks.org/dsa/segment-tree-data-structure/

#include <bits/stdc++.h>
using namespace std;

const int INF = 1e9;
const int MAXN = 1000005;

int N, Q;
vector<int> adj[MAXN];
int color[MAXN];
int distW[MAXN];

int parent[MAXN], depth[MAXN], heavy[MAXN], head[MAXN], pos[MAXN], sz[MAXN];
int curPos = 0;
int baseArray[MAXN];

struct SegTree {
    int n;
    vector<int> tree;

    SegTree(int n) : n(n), tree(4 * n, INF) {}

    void build(vector<int>& vals, int idx, int l, int r) {
        if (l == r) {
            tree[idx] = vals[l];
            return;
        }
        int mid = (l + r) / 2;
        build(vals, 2 * idx, l, mid);
        build(vals, 2 * idx + 1, mid + 1, r);
        tree[idx] = min(tree[2 * idx], tree[2 * idx + 1]);
    }

    int query(int idx, int l, int r, int ql, int qr) {
        if (qr < l || ql > r) return INF;
        if (ql <= l && r <= qr) return tree[idx];
        int mid = (l + r) / 2;
        return min(query(2 * idx, l, mid, ql, qr),
                   query(2 * idx + 1, mid + 1, r, ql, qr));
    }

    int query(int l, int r) { return query(1, 0, n - 1, l, r); }
};

void bfs_from_whites() {
    fill(distW, distW + N + 1, INF);
    queue<int> q;
    for (int i = 1; i <= N; i++) {
        if (color[i] == 1) {
            distW[i] = 0;
            q.push(i);
        }
    }
    while (!q.empty()) {
        int u = q.front(); q.pop();
        for (int v : adj[u]) {
            if (distW[v] > distW[u] + 1) {
                distW[v] = distW[u] + 1;
                q.push(v);
            }
        }
    }
}

int dfs(int u, int p) {
    sz[u] = 1;
    parent[u] = p;
    depth[u] = (p == -1 ? 0 : depth[p] + 1);
    int maxSub = 0;
    heavy[u] = -1;

    for (int v : adj[u]) {
        if (v == p) continue;
        int sub = dfs(v, u);
        sz[u] += sub;
        if (sub > maxSub) {
            maxSub = sub;
            heavy[u] = v;
        }
    }
    return sz[u];
}

void decompose(int u, int h) {
    head[u] = h;
    pos[u] = curPos++;
    baseArray[pos[u]] = distW[u];

    if (heavy[u] != -1)
        decompose(heavy[u], h);

    for (int v : adj[u]) {
        if (v != parent[u] && v != heavy[u])
            decompose(v, v);
    }
}

int queryPath(int u, int v, SegTree& seg) {
    int res = INF;
    while (head[u] != head[v]) {
        if (depth[head[u]] < depth[head[v]]) swap(u, v);
        res = min(res, seg.query(pos[head[u]], pos[u]));
        u = parent[head[u]];
    }
    if (depth[u] > depth[v]) swap(u, v);
    res = min(res, seg.query(pos[u], pos[v]));
    return res;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    cin >> N >> Q;

    for (int i = 1; i <= N; i++) cin >> color[i];

    for (int i = 0; i < N - 1; i++) {
        int a, b;
        cin >> a >> b;
        adj[a].push_back(b);
        adj[b].push_back(a);
    }

    bfs_from_whites();

    dfs(1, -1);
    decompose(1, 1);

    vector<int> vals(curPos);
    for (int i = 0; i < curPos; i++) vals[i] = baseArray[i];
    SegTree seg(curPos);
    seg.build(vals, 1, 0, curPos - 1);

    while (Q--) {
        int u, v;
        cin >> u >> v;
        cout << queryPath(u, v, seg) << '\n';
    }

    return 0;
}
// Author: Morteza Farrokhnejad, Ali Farrokhnejad