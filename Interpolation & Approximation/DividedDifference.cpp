// Implemented by 2207008

#include <bits/stdc++.h>
using namespace std;

int main()
{
    ifstream in("DividedDifference_input.txt");
    ofstream out("DividedDifference_output.txt");

    if (!in || !out)
    {
        out << "Error opening input/output file!" << endl;
        return 1;
    }

    int n;
    in >> n; // Number of data points
    vector<double> x(n), y(n);
    for (int i = 0; i < n; i++)
        in >> x[i] >> y[i]; // Input data points (x, y)

    vector<vector<double>> d(n, vector<double>(n, 0));
    for (int i = 0; i < n; i++)
        d[i][0] = y[i];

    for (int j = 1; j < n; j++)
    {
        for (int i = 0; i < n - j; i++)
        {
            d[i][j] = (d[i + 1][j - 1] - d[i][j - 1]) / (x[i + j] - x[i]);
        }
    }

    double X;
    in >> X; // Value of x to interpolate

    double ans = d[0][0];
    double p = 1.0;

    for (int j = 1; j < n; j++)
    {
        p *= (X - x[j - 1]);
        ans += p * d[0][j];
    }

    out << "Interpolated value: " << ans << endl;
    in.close();
    out.close();
    return 0;
}
