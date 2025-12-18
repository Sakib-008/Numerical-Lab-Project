// RK_Method implemented by 2207008

#include <bits/stdc++.h>
using namespace std;

double f(double x, double y)
{
    return 2 * x + 1;
}

int main()
{
    ifstream in("RK_input.txt");
    ofstream out("RK_output.txt");

    if (!in || !out)
    {
        cout << "Error opening input/output file!" << endl;
        return 1;
    }

    double x0, y0, xn, yn, h;
    in >> x0 >> y0; // Initial values x0 and y0
    in >> xn;       // Final values of x
    in >> h;        // Step size

    if (h <= 0)
    {
        cout << "Step size h must be positive!" << endl;
        return 1;
    }

    if (xn <= x0)
    {
        cout << "xn must be greater than x0!" << endl;
        return 1;
    }

    int steps = floor((xn - x0) / h);

    if (steps == 0)
    {
        cout << "Step size is too large for the given range!" << endl;
        return 1;
    }

    out << fixed << setprecision(3);
    yn = y0;

    out << "RK4 Method: " << endl;
    for (int i = 1; i <= steps; i++)
    {
        double k1 = h * f(x0, y0);
        double k2 = h * f(x0 + h / 2, y0 + k1 / 2);
        double k3 = h * f(x0 + h / 2, y0 + k2 / 2);
        double k4 = h * f(x0 + h, y0 + k3);
        yn = y0 + (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        x0 = x0 + h;
        y0 = yn;
        out << "Step " << i << ": x = " << x0 << "  y = " << y0 << endl;
    }

    out << "\nFinal Result: yn = " << yn << endl;
    in.close();
    out.close();

    return 0;
}