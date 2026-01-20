#include <bits/stdc++.h>

using namespace std;
// Author: Morteza Farrokhnejad, Ali Farrokhnejad

typedef long long ll;

const ll MOD = 1000000007LL;

const int MAX_BITS = 64;

ll dp_cache[MAX_BITS][2][2][2][2];

int x_bits[MAX_BITS], y_bits[MAX_BITS], z_bits[MAX_BITS], s_bits[MAX_BITS];

ll recursive_dp(int position, int tx, int ty, int tz, int ts) {
    if (position == MAX_BITS) return 1LL;
    ll &result = dp_cache[position][tx][ty][tz][ts];
    if (result != -1LL) return result;
    result = 0LL;
    int max_bx = tx ? x_bits[position] : 1;
    int max_by = ty ? y_bits[position] : 1;
    int max_bz = tz ? z_bits[position] : 1;
    for (int bit_x = 0; bit_x <= max_bx; bit_x++) {
        for (int bit_y = 0; bit_y <= max_by; bit_y++) {
            for (int bit_z = 0; bit_z <= max_bz; bit_z++) {
                int bit_xor = bit_x ^ bit_y ^ bit_z;
                if (ts && bit_xor > s_bits[position]) continue;
                int nxt_tx = tx && (bit_x == max_bx);
                int nxt_ty = ty && (bit_y == max_by);
                int nxt_tz = tz && (bit_z == max_bz);
                int nxt_ts = ts && (bit_xor == s_bits[position]);
                result = (result + recursive_dp(position + 1, nxt_tx, nxt_ty, nxt_tz, nxt_ts)) % MOD;
            }
        }
    }
    return result;
}

ll compute_prefix(ll max_x, ll max_y, ll max_z, ll max_n) {
    if (max_x < 0 || max_y < 0 || max_z < 0 || max_n < 0) return 0LL;
    memset(dp_cache, -1, sizeof(dp_cache));
    for (int i = 0; i < MAX_BITS; i++) {
        x_bits[i] = (max_x >> (MAX_BITS - 1 - i)) & 1;
        y_bits[i] = (max_y >> (MAX_BITS - 1 - i)) & 1;
        z_bits[i] = (max_z >> (MAX_BITS - 1 - i)) & 1;
        s_bits[i] = (max_n >> (MAX_BITS - 1 - i)) & 1;
    }
    return recursive_dp(0, 1, 1, 1, 1);
}

ll compute_range_count(ll low_b, ll high_b, ll low_c, ll high_c, ll low_d, ll high_d, ll max_n) {
    if (max_n < 0) return 0LL;
    vector<tuple<ll, ll, ll, int>> terms = {
        {high_b, high_c, high_d, 1},
        {high_b, high_c, low_d - 1, -1},
        {high_b, low_c - 1, high_d, -1},
        {low_b - 1, high_c, high_d, -1},
        {high_b, low_c - 1, low_d - 1, 1},
        {low_b - 1, high_c, low_d - 1, 1},
        {low_b - 1, low_c - 1, high_d, 1},
        {low_b - 1, low_c - 1, low_d - 1, -1}
    };
    ll answer = 0LL;
    for (auto& term : terms) {
        ll ux, uy, uz;
        int sign;
        tie(ux, uy, uz, sign) = term;
        ll value = compute_prefix(ux, uy, uz, max_n);
        ll contrib = (sign == 1) ? (value % MOD) : ((MOD - (value % MOD)) % MOD);
        answer = (answer + contrib) % MOD;
    }
    return answer;
}

ll query_processor(ll low_a, ll high_a, ll low_b, ll high_b, ll low_c, ll high_c, ll low_d, ll high_d) {
    ll count_a = high_a - low_a + 1;
    ll count_b = high_b - low_b + 1;
    ll count_c = high_c - low_c + 1;
    ll count_d = high_d - low_d + 1;
    ll total = count_a % MOD;
    total = (total * (count_b % MOD)) % MOD;
    total = (total * (count_c % MOD)) % MOD;
    total = (total * (count_d % MOD)) % MOD;
    ll upper_count = compute_range_count(low_b, high_b, low_c, high_c, low_d, high_d, high_a);
    ll lower_count = (low_a > 0) ? compute_range_count(low_b, high_b, low_c, high_c, low_d, high_d, low_a - 1) : 0LL;
    ll lose = (upper_count - lower_count + 2LL * MOD) % MOD;
    ll win = (total - lose + MOD) % MOD;
    return win;
}

int main() {
    ios::sync_with_stdio(false);
    cin.tie(NULL);
    int num_queries;
    cin >> num_queries;
    for (int i = 0; i < num_queries; i++) {
        ll low_a, high_a, low_b, high_b, low_c, high_c, low_d, high_d;
        cin >> low_a >> high_a >> low_b >> high_b >> low_c >> high_c >> low_d >> high_d;
        cout << query_processor(low_a, high_a, low_b, high_b, low_c, high_c, low_d, high_d) << '\n';
    }
    return 0;
}