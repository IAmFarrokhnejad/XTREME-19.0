// https://en.wikipedia.org/wiki/Composition_of_functions
// https://mathworld.wolfram.com/MobiusTransformation.html
// FFT/NTT from codeforces (used in bitonic)
// Authors: Morteza Farrokhnejad, Ali Farrokhnejad

#include <bits/stdc++.h>
using namespace std;
using ll = long long;

const ll MOD = 998244353;
const ll G = 3;

ll modpow(ll a, ll e) {
    a %= MOD; if (a < 0) a += MOD;
    ll r = 1;
    while (e) {
        if (e & 1) r = (r * a) % MOD;
        a = (a * a) % MOD;
        e >>= 1;
    }
    return r;
}
ll invmod(ll x) { x %= MOD; if (x < 0) x += MOD; return modpow(x, MOD-2); }

struct NTT {
    vector<int> rev;
    vector<ll> roots;
    NTT() { roots.resize(1); roots[0] = 1; }

    void ensure(int n) {
        int logn = __builtin_ctz(n);
        if ((int)rev.size() == n) return;
        rev.assign(n, 0);
        for (int i = 1; i < n; ++i)
            rev[i] = (rev[i>>1] >> 1) | ((i&1) << (logn-1));
        if ((int)roots.size() <= logn) {
            int old = (int)roots.size();
            roots.resize(logn+1);
            for (int k = old; k <= logn; ++k)
                roots[k] = modpow(G, (MOD-1) >> k);
        }
    }

    void ntt(vector<ll>& a, bool invert) {
        int n = (int)a.size();
        ensure(n);
        for (int i = 0; i < n; ++i) if (i < rev[i]) swap(a[i], a[rev[i]]);
        for (int len = 1, lvl = 1; len < n; len <<= 1, ++lvl) {
            ll wlen = roots[lvl];
            if (invert) wlen = invmod(wlen);
            for (int i = 0; i < n; i += (len << 1)) {
                ll w = 1;
                for (int j = 0; j < len; ++j) {
                    ll u = a[i+j];
                    ll v = a[i+j+len] * w % MOD;
                    a[i+j] = u + v; if (a[i+j] >= MOD) a[i+j] -= MOD;
                    a[i+j+len] = u - v; if (a[i+j+len] < 0) a[i+j+len] += MOD;
                    w = w * wlen % MOD;
                }
            }
        }
        if (invert) {
            ll invn = invmod(n);
            for (int i = 0; i < n; ++i) a[i] = a[i] * invn % MOD;
        }
    }

    vector<ll> multiply(const vector<ll>& A, const vector<ll>& B) {
        if (A.empty() || B.empty()) return {};
        int need = (int)A.size() + (int)B.size() - 1;
        int n = 1; while (n < need) n <<= 1;
        vector<ll> fa(n,0), fb(n,0);
        for (size_t i = 0; i < A.size(); ++i) fa[i] = A[i];
        for (size_t i = 0; i < B.size(); ++i) fb[i] = B[i];
        ntt(fa, false); ntt(fb, false);
        for (int i = 0; i < n; ++i) fa[i] = fa[i] * fb[i] % MOD;
        ntt(fa, true);
        fa.resize(need);
        while (!fa.empty() && fa.back() == 0) fa.pop_back();
        return fa;
    }
} ntt_engine;

// helpers
vector<ll> trunc(const vector<ll>& A, size_t n) {
    vector<ll> R(min(A.size(), n));
    for (size_t i = 0; i < R.size(); ++i) R[i] = A[i];
    return R;
}

vector<ll> poly_inv(const vector<ll>& A, int n) {
    vector<ll> R(1, invmod(A[0]));
    int cur = 1;
    while (cur < n) {
        cur <<= 1;
        vector<ll> Aac = trunc(A, cur);
        vector<ll> R2 = ntt_engine.multiply(R, R);
        R2 = ntt_engine.multiply(R2, Aac);
        R2.resize(cur);
        if ((int)R.size() < cur) R.resize(cur, 0);
        for (int i = 0; i < cur; ++i) {
            ll val = ((2LL * (i < (int)R.size() ? R[i] : 0) % MOD) - (i < (int)R2.size() ? R2[i] : 0)) % MOD;
            if (val < 0) val += MOD;
            R[i] = val;
        }
    }
    R.resize(n);
    while (!R.empty() && R.back() == 0) R.pop_back();
    return R;
}

pair<vector<ll>, vector<ll>> poly_divmod(const vector<ll>& A, const vector<ll>& B) {
    int nA = (int)A.size(), nB = (int)B.size();
    if (nA < nB) return {vector<ll>{0}, A};
    int n = nA - nB + 1;
    vector<ll> Ar(n), Br(nB);
    for (int i = 0; i < n; ++i) Ar[i] = A[nA - 1 - i];
    for (int i = 0; i < nB; ++i) Br[i] = B[nB - 1 - i];
    vector<ll> invBr = poly_inv(Br, n);
    vector<ll> Qr = ntt_engine.multiply(Ar, invBr);
    Qr.resize(n);
    vector<ll> Q(n);
    for (int i = 0; i < n; ++i) Q[i] = Qr[n - 1 - i];
    vector<ll> BQ = ntt_engine.multiply(B, Q);
    vector<ll> R = A;
    if (R.size() < BQ.size()) R.resize(BQ.size(), 0);
    for (size_t i = 0; i < BQ.size(); ++i) {
        R[i] = (R[i] - BQ[i]) % MOD;
        if (R[i] < 0) R[i] += MOD;
    }
    while (!R.empty() && R.back() == 0) R.pop_back();
    while (!Q.empty() && Q.back() == 0) Q.pop_back();
    return {Q, R};
}

vector<vector<ll>> build_prod_tree(const vector<ll>& X) {
    int m = (int)X.size();
    if (m == 0) return {};
    vector<vector<ll>> tree(4*m);
    function<void(int,int,int)> build = [&](int v, int l, int r) {
        if (l == r) {
            tree[v] = vector<ll>{(MOD - X[l]) % MOD, 1};
            return;
        }
        int mid = (l + r) >> 1;
        build(v<<1, l, mid);
        build(v<<1|1, mid+1, r);
        tree[v] = ntt_engine.multiply(tree[v<<1], tree[v<<1|1]);
    };
    build(1, 0, m-1);
    return tree;
}

vector<ll> multipoint_eval(const vector<ll>& P, const vector<ll>& X) {
    int q = (int)X.size();
    vector<ll> res(q, 0);
    if (q == 0) return res;
    auto tree = build_prod_tree(X);
    int nTree = (int)tree.size();
    vector<vector<ll>> rem(nTree);
    function<void(int,int,int,const vector<ll>&)> solve = [&](int v, int l, int r, const vector<ll>& poly) {
        if (poly.empty()) rem[v].clear();
        else rem[v] = poly_divmod(poly, tree[v]).second;
        if (l == r) return;
        int mid = (l + r) >> 1;
        solve(v<<1, l, mid, rem[v]);
        solve(v<<1|1, mid+1, r, rem[v]);
    };
    solve(1, 0, q-1, P);
    function<void(int,int,int)> collect = [&](int v, int l, int r) {
        if (l == r) {
            res[l] = rem[v].empty() ? 0 : rem[v][0] % MOD;
            return;
        }
        int mid = (l + r) >> 1;
        collect(v<<1, l, mid);
        collect(v<<1|1, mid+1, r);
    };
    collect(1, 0, q-1);
    return res;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int N, M;
    if (!(cin >> N >> M)) return 0;
    vector<ll> S(N);
    for (int i = 0; i < N; ++i) {
        cin >> S[i];
        S[i] %= MOD; if (S[i] < 0) S[i] += MOD;
    }

    function<vector<ll>(int,int)> build_prod = [&](int l, int r)->vector<ll> {
        if (l == r) return vector<ll>{ (MOD - S[l]) % MOD, 1 };
        int m = (l + r) >> 1;
        auto L = build_prod(l, m);
        auto R = build_prod(m+1, r);
        return ntt_engine.multiply(L, R);
    };

    vector<ll> P;
    if (N > 0) P = build_prod(0, N-1);
    vector<ll> Pder(max(0, (int)P.size() - 1));
    for (int i = 1; i < (int)P.size(); ++i) Pder[i-1] = P[i] * i % MOD;

    vector<int> op_type(M);
    vector<ll> op_val(M);
    for (int i = 0; i < M; ++i) {
        cin >> op_type[i];
        if (op_type[i] == 1) cin >> op_val[i];
        else op_val[i] = 0;
    }

    vector<array<ll,4>> ABCD(M);
    ll a = 1, b = 0, c = 0, d = 1;
    for (int i = 0; i < M; ++i) {
        if (op_type[i] == 1) {
            ll X = op_val[i] % MOD; if (X < 0) X += MOD;
            ll na = (a + X * c) % MOD;
            ll nb = (b + X * d) % MOD;
            a = na; b = nb;
        } else {
            ll oa = a, ob = b, oc = c, od = d;
            a = oc; b = od; c = oa; d = ob;
        }
        ABCD[i][0] = a; ABCD[i][1] = b; ABCD[i][2] = c; ABCD[i][3] = d;
    }

    vector<ll> uniq;
    uniq.reserve(M);
    unordered_map<ll,int> map_idx;
    for (int i = 0; i < M; ++i) {
        a = ABCD[i][0]; b = ABCD[i][1]; c = ABCD[i][2]; d = ABCD[i][3];
        ll cm = (c % MOD + MOD) % MOD;
        if (cm == 0) continue;
        ll invc = invmod(cm);
        ll t = ( (d % MOD + MOD) % MOD ) * invc % MOD;
        ll x = (MOD - t) % MOD;
        if (map_idx.find(x) == map_idx.end()) {
            int id = (int)uniq.size();
            map_idx[x] = id;
            uniq.push_back(x);
        }
    }

    vector<ll> valP, valPd;
    if (!uniq.empty()) {
        valP = multipoint_eval(P, uniq);
        valPd = multipoint_eval(Pder, uniq);
    }

    ll sumS = 0;
    for (int i = 0; i < N; ++i) { sumS += S[i]; if (sumS >= MOD) sumS -= MOD; }

    for (int i = 0; i < M; ++i) {
        a = ABCD[i][0]; b = ABCD[i][1]; c = ABCD[i][2]; d = ABCD[i][3];
        a = (a % MOD + MOD) % MOD; b = (b % MOD + MOD) % MOD;
        c = (c % MOD + MOD) % MOD; d = (d % MOD + MOD) % MOD;
        ll ans = 0;
        if (c == 0) {
            ll invd = invmod(d);
            ans = ((a * sumS + b * (ll)N) % MOD) * invd % MOD;
        } else {
            ll invc = invmod(c);
            ll t = d * invc % MOD;
            ll x = (MOD - t) % MOD;
            int id = map_idx[x];
            ll P_x = valP[id] % MOD; if (P_x < 0) P_x += MOD;
            ll Pd_x = valPd[id] % MOD; if (Pd_x < 0) Pd_x += MOD;
            ll F = (MOD - (Pd_x * invmod(P_x) % MOD)) % MOD;
            ll S1 = invc * F % MOD;
            ll term1 = (a * (ll)N) % MOD * invc % MOD;
            ll coef = (b - (a * d % MOD) * invc % MOD) % MOD; if (coef < 0) coef += MOD;
            ans = (term1 + coef * S1) % MOD;
        }
        if (ans < 0) ans += MOD;
        cout << ans << '\n';
    }
    return 0;
}