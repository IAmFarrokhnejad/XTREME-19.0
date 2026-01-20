#include <bits/stdc++.h>
using namespace std;
// Author: Morteza Farrokhnejad, Ali Farrokhnejad
const long long MOD = 998244353;

// Copied from before
long long modpow(long long base, long long exp, long long mod = MOD) {
    long long result = 1;
    base %= mod;
    while (exp > 0) {
        if (exp & 1) result = result * base % mod;
        base = base * base % mod;
        exp >>= 1;
    }
    return result;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int test_count;
    cin >> test_count;

    while (test_count--) {
        long long n;
        cin >> n;

        if (n % 2 == 0) {
            cout << 0 << '\n';
        } else {
            long long inv_n = modpow(n, MOD - 2);
            cout << inv_n << '\n';
        }
    }

    return 0;
}
