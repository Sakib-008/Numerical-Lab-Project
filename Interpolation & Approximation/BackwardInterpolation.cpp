// Implemented by 2207008

#include <bits/stdc++.h>
using namespace std;

void printDifferenceTable(const vector<vector<double>> &d, int n, ofstream &out)
{
    out << "\nBackward Difference Table:" << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j <= i; j++)
        {
            out << setw(10) << d[i][j] << " ";
        }
        out << endl;
    }
}

int main()
{
    ifstream in("BackwardInterpolation_input.txt");
    ofstream out("BackwardInterpolation_output.txt");

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

    double h = x[1] - x[0];

    for (int i = 1; i < n - 1; i++)
    {
        if (fabs((x[i + 1] - x[i]) - h) > 1e-9)
        {
            out << "Error: x values are not equally spaced." << endl;
            in.close();
            out.close();
            return 1;
        }
    }

    double X;
    in >> X; // Value of x to interpolate
    
    vector<vector<double>> d(n, vector<double>(n, 0));
    for (int i = 0; i < n; i++)
        d[i][0] = y[i];
    for (int j = 1; j < n; j++)
    {
        for (int i = n - 1; i >= j; i--)
        {
            d[i][j] = d[i][j - 1] - d[i - 1][j - 1];
        }
    }
    double u = (X - x[n - 1]) / h;
    double p = 1.0;
    double f = 1.0;
    double ans = d[n - 1][0];
    for (int j = 1; j < n; j++)
    {
        p *= (u + (j - 1));
        f *= j;
        ans += (p * d[n - 1][j]) / f;
    }
    out << "Newton's Backward Interpolation Method : " << endl;
    printDifferenceTable(d, n, out);
    out << fixed << setprecision(3);
    out << "\nInterpolated answer = " << ans << endl;
    in.close();
    out.close();
    return 0;
}