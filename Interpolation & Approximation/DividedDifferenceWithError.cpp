// Implemented by 2207008

#include <bits/stdc++.h>
using namespace std;

int main() {
    ifstream in("DividedDifferenceWithError_input.txt");
    ofstream out("DividedDifferenceWithError_output.txt");
    
    if (!in || !out) {
        out << "Error opening input/output file!" << endl;
        return 1;
    }

    int n;
    in >> n; // Number of data points (excluding extra point for error)

    vector<double> x(n+1), y(n+1);
    for (int i = 0; i < n; i++)
        in >> x[i] >> y[i]; // Input first n data points (x, y)
    in >> x[n] >> y[n];   // extra point for error term

    double X;
    in >> X; // Value of x to interpolate

    vector<vector<double>> d(n+1, vector<double>(n+1, 0));

    for (int i = 0; i <= n; i++)
        d[i][0] = y[i];

    for (int j = 1; j <= n; j++) {
        for (int i = 0; i <= n - j; i++) {
            d[i][j] = (d[i+1][j-1] - d[i][j-1]) / (x[i+j] - x[i]);
        }
    }

    double ans = d[0][0];
    double p = 1.0;

    for (int j = 1; j < n; j++) {
        p *= (X - x[j-1]);
        ans += p * d[0][j];
    }

    out << fixed << setprecision(6);
    out << "Newton's Divided Difference Interpolation Method for Calculating Error: " << endl;
    out << "Interpolated value : " << ans << endl;

    // Error estimation
    double errP = 1.0;
    for (int i = 0; i < n; i++)
        errP *= (X - x[i]);
    double error = fabs(errP * d[0][n]);
    out << "Estimated error using the extra point: " << error << endl;
    
    return 0;
}
