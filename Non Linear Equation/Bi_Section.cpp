// Bi_Section implemented by 2207008

#include <bits/stdc++.h>
using namespace std;

double f(double x, vector<double> &v)
{
    double result = 0;
    int n = v.size();
    for (int i = 0; i < n; i++)
    {
        int deg = n - i - 1;
        result += v[i] * pow(x, deg);
    }
    return result;
}

int main()
{
    ifstream in("BiSection_input.txt");
    ofstream out("BiSection_output.txt");

    if (!in || !out)
    {
        cout << "Error opening input/output file!" << endl;
        return 1;
    }

    int n;
    in >> n; // Degree of the polynomial equation
    double stepSize, e;
    in >> stepSize >> e; // Step size and tolerance
    vector<double> eqn(n + 1);
    for (int i = 0; i <= n; i++)
    {
        in >> eqn[i]; // Coefficients of the polynomial
    }

    out << "Bisection method : " << endl;

    out << "The inputted equation is : ";

    for (int i = 0; i <= n; i++)
    {
        if (eqn[i] == 0)
            continue;
        if (i != 0 && eqn[i] > 0)
            out << " + ";
        if (eqn[i] < 0)
            out << " - ";
        out << abs(eqn[i]);
        int deg = n - i;
        if (deg > 0)
            out << "x";
        if (deg > 1)
            out << "^" << deg;
        if (i == n)
            out << " = 0" << endl;
    }

    double xmax = sqrt((eqn[1] / eqn[0]) * (eqn[1] / eqn[0]) - 2 * (eqn[2] / eqn[0]));
    double scanStart = -xmax, scanEnd = xmax;
    vector<pair<double, double>> intervals;

    double prev = scanStart;
    for (double cur = scanStart + stepSize; cur <= scanEnd; cur += stepSize)
    {
        if (f(prev, eqn) * f(cur, eqn) < 0)
        {
            intervals.push_back({prev, cur});
        }
        prev = cur;
    }

    if (intervals.size() == 0)
    {
        out << "No initial guess found in this range" << endl;
        in.close();
        out.close();
        return -1;
    }

    out << fixed << setprecision(3);
    for (int i = 0; i < intervals.size(); i++)
    {
        double x1 = intervals[i].first;
        double x2 = intervals[i].second;
        double x0, f0, f1, f2;
        int step = 0;
        out << "Bracket for root " << i + 1 << " : [" << x1 << ", " << x2 << "]" << endl;
        do
        {
            f1 = f(x1, eqn);
            f2 = f(x2, eqn);
            x0 = (x1 + x2) / 2;
            f0 = f(x0, eqn);
            if (f1 * f0 < 0)
                x2 = x0;
            else
                x1 = x0;
            step++;
        } while (fabs(f0) > e && step < 1000);

        out << "Iteration needed : " << step << endl;
        out << "Root " << i + 1 << " : " << x0 << endl;
        out << endl;
    }
    in.close();
    out.close();
    return 0;
}
