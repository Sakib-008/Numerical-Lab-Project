#include <bits/stdc++.h>

using namespace std;

// The linear equation format is: y=a+bx
void Linear(ofstream &out, vector<double> x, vector<double> y, int newX)
{
  int n = x.size();
  double Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
  for (int i = 0; i < n; i++)
  {
    Sx += x[i];
    Sy += y[i];
    Sxx += x[i] * x[i];
    Sxy += x[i] * y[i];
  }
  double b = (n * Sxy - Sx * Sy) / (n * Sxx - Sx * Sx);

  double a = (Sy - b * Sx) / n;

  out << "The Linear Equation is: " << endl;
  out << "y = " << a << " + " << b << "x" << endl;

  double newY = a + b * newX;

  out << "For x = " << newX << ", the predicted value of y is: ";
  out << newY << endl;
}

int main(void)
{
  ifstream in("LinearEquationInput.txt");
  ofstream out("LinearEquationOutput.txt");

  int n;
  in >> n;
  vector<double> x(n), y(n);

  for (int i = 0; i < n; i++)
  {
    in >> x[i] >> y[i];
  }

  out << setprecision(2) << fixed;

  out << "Data points entered:" << endl;
  for (int i = 0; i < n; i++)
  {
    out << "(" << x[i] << ", " << y[i] << ")\n";
  }

  int newX;
  in >> newX;
  Linear(out, x, y, newX);

  in.close();
  out.close();
  return 0;
}