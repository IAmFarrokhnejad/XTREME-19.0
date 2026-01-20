#include <bits/stdc++.h>
using namespace std;

int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    int N;
    cin >> N;
    vector<int> P(N + 1), pos(N + 1);

    for (int i = 0; i < N; ++i) {
        cin >> P[i];
        pos[P[i]] = i;
    }

    long long ans = LLONG_MAX;

    for (int dir = 0; dir < 2; ++dir) {
        vector<int> d(N);

        for (int x = 1; x <= N; ++x) {
            if (dir == 0) {
                d[x - 1] = (pos[x] - (x - 1)) % N;}
            else{
                d[x - 1] = (pos[x] + (x - 1)) % N;
}
            if (d[x - 1] < 0){
                d[x - 1] += N;}
        }

        sort(d.begin(), d.end());

        int max_gap = 0;
        for (int i = 0; i < N - 1; ++i)
            max_gap = max(max_gap, d[i + 1] - d[i]);
        max_gap = max(max_gap, N + d[0] - d[N - 1]);

        long long res = (N - max_gap + 1) / 2;
        ans = min(ans, res);
    }

    cout << ans << "\n";
    return 0;
}
// Author: Morteza Farrokhnejad, Ali Farrokhnejad