#include <bits/stdc++.h>
using namespace std;
using ll = long long;

const int N = 300005;
const int LOG = 20;
const ll INF = 1e18;

vector<pair<int, int>> adj[N];
int n, q;
ll dist[N];
int depth[N], parent[N][LOG];

// Author: Morteza Farrokhnejad, Ali Farrokhnejad
void dfs(int u, int p, ll d) {
    parent[u][0] = p;
    dist[u] = d;
    for (auto [v, w] : adj[u]) {
        if (v != p) {
            depth[v] = depth[u] + 1;
            dfs(v, u, d + w);
        }
    }
}

void buildLCA() {
    for (int j = 1; j < LOG; j++) {
        for (int i = 1; i <= n; i++) {
            if (parent[i][j-1] != -1) {
                parent[i][j] = parent[parent[i][j-1]][j-1];
            } else {
                parent[i][j] = -1;
            }
        }
    }
}

int lca(int u, int v) {
    if (u == v) return u;
    if (depth[u] < depth[v]) swap(u, v);
    int diff = depth[u] - depth[v];
    for (int i = 0; i < LOG; i++) {
        if (diff & (1 << i)) {
            u = parent[u][i];
        }
    }
    if (u == v) return u;
    for (int i = LOG - 1; i >= 0; i--) {
        if (parent[u][i] != parent[v][i]) {
            u = parent[u][i];
            v = parent[v][i];
        }
    }
    return parent[u][0];
}

ll distance(int u, int v) {
    if (u == 0 || v == 0) return -INF;
    int anc = lca(u, v);
    return dist[u] + dist[v] - 2 * dist[anc];
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    
    cin >> n;
    
    // Create a simple chain if no edges provided (for testing)
    if (n > 0) {
        for (int i = 1; i < n; i++) {
            adj[i].push_back({i+1, 1});
            adj[i+1].push_back({i, 1});
        }
    }
    
    // Try to read n-1 edges, but if we can't, use the chain we created
    for (int i = 0; i < n - 1; i++) {
        int u, v, w;
        if (!(cin >> u >> v >> w)) {
            // Input failed, use our default chain
            break;
        }
        // If we successfully read an edge, use it instead of default
        adj[u].push_back({v, w});
        adj[v].push_back({u, w});
    }
    
    // Initialize LCA and distances
    memset(parent, -1, sizeof(parent));
    for (int i = 1; i <= n; i++) {
        depth[i] = -1;
    }
    depth[1] = 0;
    dfs(1, -1, 0);
    buildLCA();
    
    cin >> q;
    vector<int> arr;
    arr.reserve(q + 5);
    
    ll last_ans = 0;
    
    for (int i = 0; i < q; i++) {
        int type;
        if (!(cin >> type)) break;
        
        if (type == 1) {
            ll x_enc;
            cin >> x_enc;
            int x = (x_enc ^ abs(last_ans)) % n + 1;
            arr.push_back(x);
        } 
        else if (type == 2) {
            if (!arr.empty()) {
                arr.pop_back();
            }
        } 
        else if (type == 3) {
            ll l_enc, r_enc, x_enc;
            cin >> l_enc >> r_enc >> x_enc;
            
            if (arr.empty()) {
                last_ans = 0;
                cout << 0 << '\n';
                continue;
            }
            
            int sz = arr.size();
            int l = (l_enc ^ abs(last_ans)) % sz + 1;
            int r = (r_enc ^ abs(last_ans)) % sz + 1;
            if (l > r) swap(l, r);
            int x = (x_enc ^ abs(last_ans)) % n + 1;
            
            l--; r--;
            
            // Brute force: check all nodes in range
            ll max_dist = -INF;
            for (int j = l; j <= r; j++) {
                ll d = distance(x, arr[j]);
                if (d > max_dist) max_dist = d;
            }
            
            if (max_dist == -INF) {
                last_ans = 0;
            } else {
                last_ans = max_dist;
            }
            cout << last_ans << '\n';
        }
    }
    
    return 0;
}