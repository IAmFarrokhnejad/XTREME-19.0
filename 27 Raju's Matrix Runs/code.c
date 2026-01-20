#include <stdio.h>
#include <stdlib.h>

#define MAXN 1000
#define MAXM 1000

int N, M;
int grid[MAXN][MAXM];
int dp[MAXN][MAXM];

int dirs[4][2] = { {-1,0}, {1,0}, {0,-1}, {0,1} };

int max(int a, int b) {
    return a > b ? a : b;
}

// Recursive DFS with memoization
int longest_from(int i, int j) {
    if (dp[i][j] != -1)
        return dp[i][j];

    dp[i][j] = 1; // start with current cell

    for (int k = 0; k < 4; k++) {
        int x = i + dirs[k][0];
        int y = j + dirs[k][1];
        if (x >= 0 && x < N && y >= 0 && y < M && grid[x][y] > grid[i][j]) {
            dp[i][j] = max(dp[i][j], 1 + longest_from(x, y));
        }
    }

    return dp[i][j];
}

int main() {
    if (scanf("%d %d", &N, &M) != 2)
        return 1;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            scanf("%d", &grid[i][j]);
            dp[i][j] = -1; // initialize dp
        }
    }

    int max_length = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            max_length = max(max_length, longest_from(i, j));
        }
    }

    printf("%d\n", max_length);
    return 0;
}
// Author: Morteza Farrokhnejad, Ali Farrokhnejad