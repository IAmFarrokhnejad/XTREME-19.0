#include <bits/stdc++.h>
using namespace std;
//Author: Morteza Farrokhnejad, Ali Farrokhnejad

const int INF = 2000000000;
const int MAXN = 600010;

int tree[3][4 * MAXN];
string s;
int n;

void build(int d, int node, int l, int r) {
    if (l == r) {
        tree[d][node] = (s[l - 1] - '0' == d ? l : INF);
        return;
    }
    int m = (l + r) / 2;
    build(d, node * 2, l, m);
    build(d, node * 2 + 1, m + 1, r);
    tree[d][node] = min(tree[d][node * 2], tree[d][node * 2 + 1]);
}

void update(int d, int node, int l, int r, int pos) {
    if (l == r) {
        tree[d][node] = INF;
        return;
    }
    int m = (l + r) / 2;
    if (pos <= m)
        update(d, node * 2, l, m, pos);
    else
        update(d, node * 2 + 1, m + 1, r, pos);
    tree[d][node] = min(tree[d][node * 2], tree[d][node * 2 + 1]);
}

int query(int d, int node, int l, int r, int ql, int qr) {
    if (ql > r || qr < l) return INF;
    if (ql <= l && r <= qr) return tree[d][node];
    int m = (l + r) / 2;
    int p1 = query(d, node * 2, l, m, ql, qr);
    int p2 = query(d, node * 2 + 1, m + 1, r, ql, qr);
    return min(p1, p2);
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);
    
    int T;
    cin >> T;
    while (T--) {
        cin >> s;
        n = s.length();

        // Build segment trees for digits 0, 1, 2
        for (int d = 0; d < 3; ++d) {
            build(d, 1, 1, n);
        }

        int k = n / 6;
        int need[6] = {1, 1, 2, 0, 1, 2};

        for (int p = 0; p < k; ++p) {
            int prev = 0;
            for (int j = 0; j < 6; ++j) {
                int d = need[j];
                int pos = query(d, 1, 1, n, prev + 1, n);
                if (j) cout << " ";
                cout << pos;
                update(d, 1, 1, n, pos);
                prev = pos;
            }
            cout << "\n";
        }
    }
    return 0;
}