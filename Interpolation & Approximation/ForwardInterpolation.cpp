// Implemented by 2207008

#include <bits/stdc++.h>
using namespace std;

void printDifferenceTable(const vector<vector<double>> &d, int n, ofstream &out)
{
    out << "\nDifference Table:" << endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n - i; j++)
        {
            out << setw(10) << d[i][j] << " ";
        }
        out << endl;
    }
}

int main()
{
    ifstream in("ForwardInterpolation_input.txt");
    ofstream out("ForwardInterpolation_output.txt");

    if (!in || !out)
    {
        cout << "Error opening input/output file!" << endl;
        return 1;
    }

    int n; // Number of data points
    in >> n;
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

    double X; // Value of x to interpolate
    in >> X;
    vector<vector<double>> d(n, vector<double>(n, 0));
    for (int i = 0; i < n; i++)
        d[i][0] = y[i];
    for (int j = 1; j < n; j++)
    {
        for (int i = 0; i < n - j; i++)
        {
            d[i][j] = d[i + 1][j - 1] - d[i][j - 1];
        }
    }
    double u = (X - x[0]) / h;
    double p = 1.0;
    double f = 1.0;
    double ans = d[0][0];
    for (int j = 1; j < n; j++)
    {
        p *= (u - (j - 1));
        f *= j;
        ans += (d[0][j] * p) / f;
    }
    out << "Newton's Forward Interpolation Method : " << endl;
    printDifferenceTable(d, n, out);
    out << fixed << setprecision(3);
    out << "\nInterpolated answer : " << ans << endl;
    in.close();
    out.close();
    return 0;
}
