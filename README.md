# Table of Contents

- [Solution of Linear Equations](#solution-of-linear-equations)

  - [Gauss Elimination Method](#gauss-elimination-method)
    - [Theory](#gauss-elimination-theory)
    - [Code](#gauss-elimination-code)
    - [Input](#gauss-elimination-input)
    - [Output](#gauss-elimination-output)
  - [Gauss Jordan Elimination Method](#gauss-jordan-elimination-method)
    - [Theory](#gauss-jordan-theory)
    - [Code](#gauss-jordan-code)
    - [Input](#gauss-jordan-input)
    - [Output](#gauss-jordan-output)
  - [LU Decomposition Method](#lu-decomposition-method)
    - [Theory](#lu-decomposition-theory)
    - [Code](#lu-decomposition-code)
    - [Input](#lu-decomposition-input)
    - [Output](#lu-decomposition-output)
  - [Matrix Inversion Method](#matrix-inversion-method)
    - [Theory](#matrix-inversion-theory)
    - [Code](#matrix-inversion-code)
    - [Input](#matrix-inversion-input)
    - [Output](#matrix-inversion-output)

- [Solution of Non-Linear Equations](#solution-of-non-linear-equations)

  - [Bisection Method](#bisection-method)
    - [Theory](#bisection-theory)
    - [Code](#bisection-code)
    - [Input](#bisection-input)
    - [Output](#bisection-output)
  - [False Position Method](#false-position-method)
    - [Theory](#false-position-theory)
    - [Code](#false-position-code)
    - [Input](#false-position-input)
    - [Output](#false-position-output)
  - [Secant Method](#secant-method)
    - [Theory](#secant-theory)
    - [Code](#secant-code)
    - [Input](#secant-input)
    - [Output](#secant-output)
  - [Newton-Raphson Method](#newton-raphson-method)
    - [Theory](#newton-raphson-theory)
    - [Code](#newton-raphson-code)
    - [Input](#newton-raphson-input)
    - [Output](#newton-raphson-output)

- [Solution of Differential Equations](#solution-of-differential-equations)

  - [Runge-Kutta Method](#runge-kutta-method)
    - [Theory](#runge-kutta-theory)
    - [Code](#runge-kutta-code)
    - [Input](#runge-kutta-input)
    - [Output](#runge-kutta-output)

- [Numerical Integration](#numerical-integration)

  - [Simpson's one-third rule](#simpsons-one-third-rule)
    - [Theory](#simpsons-one-third-rule-theory)
    - [Code](#simpsons-one-third-rule-code)
    - [Input](#simpsons-one-third-rule-input)
    - [Output](#simpsons-one-third-rule-output)
  - [Simpson's three-eighths rule](#simpsons-three-eighths-rule)
    - [Theory](#simpsons-three-eighths-rule-theory)
    - [Code](#simpsons-three-eighths-rule-code)
    - [Input](#simpsons-three-eighths-rule-input)
    - [Output](#simpsons-three-eighths-rule-output)

- [Numerical Differentiation](#numerical-differentiation)

  - [Differentiation](#differentiation)
    - [Theory](#differentiation-theory)
    - [Code](#differentiation-code)
    - [Input](#differentiation-input)
    - [Output](#differentiation-output)

- [Curve Fitting Regression](#curve-fitting-regression)

  - [Linear Equation](#linear-equation)
    - [Theory](#linear-equation-theory)
    - [Code](#linear-equation-code)
    - [Input](#linear-equation-input)
    - [Output](#linear-equation-output)
  - [Polynomial Equation](#polynomial-equation)
    - [Theory](#polynomial-equation-theory)
    - [Code](#polynomial-equation-code)
    - [Input](#polynomial-equation-input)
    - [Output](#polynomial-equation-output)
  - [Transcendental Equation](#transcendental-equation)
    - [Theory](#transcendental-equation-theory)
    - [Code](#transcendental-equation-code)
    - [Input](#transcendental-equation-input)
    - [Output](#transcendental-equation-output)

- [Interpolation and Approximation](#interpolation-and-approximation)

  - [Newton's Forward Interpolation](#newtons-forward-interpolation)
    - [Theory](#newtons-forward-interpolation-theory)
    - [Code](#newtons-forward-interpolation-code)
    - [Input](#newtons-forward-interpolation-input)
    - [Output](#newtons-forward-interpolation-output)
  - [Newton's Backward Interpolation](#newtons-backward-interpolation)
    - [Theory](#newtons-backward-interpolation-theory)
    - [Code](#newtons-backward-interpolation-code)
    - [Input](#newtons-backward-interpolation-input)
    - [Output](#newtons-backward-interpolation-output)
  - [Divided Difference Interpolation](#divided-difference-interpolation)
    - [Theory](#divided-difference-interpolation-theory)
    - [Code](#divided-difference-interpolation-code)
    - [Input](#divided-difference-interpolation-input)
    - [Output](#divided-difference-interpolation-output)
  - [Divided Difference Interpolation with Error](#divided-difference-interpolation-with-error)
    - [Theory](#divided-difference-interpolation-with-error-theory)
    - [Code](#divided-difference-interpolation-with-error-code)
    - [Input](#divided-difference-interpolation-with-error-input)
    - [Output](#divided-difference-interpolation-with-error-output)

---

### Solution of Linear Equations

### Gauss Elimination Method

#### Gauss Elimination Theory

Gauss Elimination is a method to solve a system of linear equations by converting the matrix into Row Echelon Form using elementary row operations. In this method, the given system of equations is written in augmented matrix form [A | H]. Then, forward elimination is performed to the matrix which converts the matrix into an upper triangular (Row Echelon Form) matrix. In upper triangular matrix, all elements below the main diagonal are zero. After that, the values of the variables are computed using back substitution starting from the last equation in Row Echelon Form matrix.

Algorithm:

1. The system of linear equations is written in augmented matrix form.
2. The first element of the first row is chosen as the pivot element.
3. All elements below the pivot element are converted to zero using elementary row operations.
4. Moving to the next row and next column, the forward elimination process is performed until the matrix converts into an upper triangular (row echelon form) matrix.
5. Lastly, back substitution is used to find the solution of the linear equations.

#### Gauss Elimination Code

```cpp
// Gauss Elimination Method implemented by 2207009
#include <bits/stdc++.h>
using namespace std;

void printMatrix(ofstream &out, vector<vector<double>> mat)
{
  int n = mat.size();
  int m = mat[0].size();
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
      out << mat[i][j] << " ";
    out << endl;
  }
}

void gaussElimination(ofstream &out, vector<vector<double>> A, int n)
{
  int flag = 0;

  for (int i = 0; i < n; i++)
  {
    if (A[i][i] == 0)
    {
      for (int r = i + 1; r < n; r++)
      {
        if (A[r][i] != 0)
        {
          swap(A[i], A[r]);
          break;
        }
      }
      if (A[i][i] == 0)
      {
        flag = 1;
        break;
      }
    }

    for (int j = i + 1; j < n; j++)
    {
      double ratio = A[j][i] / A[i][i];
      for (int k = 0; k < n + 1; k++)
      {
        A[j][k] = A[j][k] - ratio * A[i][k];
      }
    }
  }

  if (flag == 1)
  {
    bool noSolution = false;
    bool infiniteSolution = false;

    for (int i = 0; i < n; i++)
    {
      bool allZero = true;
      for (int j = 0; j < n; j++)
      {
        if (A[i][j] > 0)
        {
          allZero = false;
          break;
        }
      }

      if (allZero && A[i][n] != 0)
        noSolution = true;
      if (allZero && A[i][n] == 0)
        infiniteSolution = true;
    }
    if (noSolution)
    {
      out << "No solution" << endl;
    }
    else if (infiniteSolution)
    {
      out << "Infinite solutions" << endl;
    }
  }
  else
  {
    out << "The Row Echelon Form is: " << endl;
    printMatrix(out, A);
    vector<double> X(n);

    for (int i = n - 1; i >= 0; i--)
    {
      double sum = A[i][n];
      for (int j = i + 1; j < n; j++)
      {
        sum -= A[i][j] * X[j];
      }
      X[i] = sum / A[i][i];
    }

    out << "The unique solution is: " << endl;
    for (int i = 0; i < n; i++)
      out << "x" << i + 1 << " = " << X[i] << (i == n - 1 ? "" : " , ");
    out << endl;
  }
}

int main(void)
{
  ifstream in("GaussEliminationInput.txt");
  ofstream out("GaussEliminationOutput.txt");

  int T;
  in >> T;

  while (T--)
  {
    int n;
    in >> n;

    vector<vector<double>> A(n, vector<double>(n + 1));

    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
        in >> A[i][j];
      in >> A[i][n];
    }

    out << setprecision(3) << fixed;

    out << "The system of equations are: " << endl;
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        if (j > 0)
          out << (A[i][j] < 0 ? " - " : " + ");

        out << abs(A[i][j]) << "x" << j + 1;
      }
      out << " = " << A[i][n] << endl;
    }
    gaussElimination(out, A, n);
  }

  in.close();
  out.close();

  return 0;
}
```

#### Gauss Elimination Input

```
3
5
2 1 -1 3 2 9
1 3 2 -1 1 8
3 2 4 1 -2 20
2 1 3 2 1 17
1 -1 2 3 4 15
2
1 1 2
2 2 4
2
1 1 2
2 2 5
```

##### Input Format

```
The input is taken from a file named GaussEliminationInput.txt.

The first line of input contains an integer T - the number of test cases.

For each test case:

The first line contains an integer n - the number of equations.

The next n lines each contain n + 1 real numbers (augmented matrix).
```

#### Gauss Elimination Output

```
The system of equations are:
2.000x1 + 1.000x2 - 1.000x3 + 3.000x4 + 2.000x5 = 9.000
1.000x1 + 3.000x2 + 2.000x3 - 1.000x4 + 1.000x5 = 8.000
3.000x1 + 2.000x2 + 4.000x3 + 1.000x4 - 2.000x5 = 20.000
2.000x1 + 1.000x2 + 3.000x3 + 2.000x4 + 1.000x5 = 17.000
1.000x1 - 1.000x2 + 2.000x3 + 3.000x4 + 4.000x5 = 15.000
The Row Echelon Form is:
2.000 1.000 -1.000 3.000 2.000 9.000
0.000 2.500 2.500 -2.500 0.000 3.500
0.000 0.000 5.000 -3.000 -5.000 5.800
0.000 0.000 0.000 1.400 3.000 3.360
0.000 0.000 0.000 0.000 1.857 2.200
The unique solution is:
x1 = 5.154 , x2 = -1.000 , x3 = 2.262 , x4 = -0.138 , x5 = 1.185
The system of equations are:
1.000x1 + 1.000x2 = 2.000
2.000x1 + 2.000x2 = 4.000
Infinite solutions
The system of equations are:
1.000x1 + 1.000x2 = 2.000
2.000x1 + 2.000x2 = 5.000
No solution
```

##### Output Format

```
The output is written to a file named GaussEliminationOutput.txt.

For each test case, the program prints:

The system of linear equations.

Nature of the solution: (Unique solution, No solution, Infinite solutions)

If the solution is unique, then print The Row Echelon Form Matrix.

The solution vector.

else If there is no solution, then print No solution.

else If there is infinite solutions, then print Infinite solutions.

All floating-point values are printed with 3 decimal places.
```

---

### Gauss Jordan Elimination Method

#### Gauss Jordan Theory

The Gauss Jordan Elimination Method is used to solve a system of linear equations by using forward and backward elimination. In this method, the given system of equations is written in augmented matrix form [A | H]. Then, the matrix is converted into Reduced Row Echelon Form. That is, the diagonal elements of the matrix are made 1 and all other elements in each pivot column are made zero. After that, the solution is obtained from the constant vector part of the Reduced Row Echelon Form matrix.

Algorithm:

1. The system of linear equations is written in augmented matrix form.
2. A pivot element is selected and made equal to 1 by dividing the entire row by the pivot element.
3. All other elements in the pivot column made zero by using forward and backward elimination.
4. Move to the next pivot position and repeat the process until the coefficient part of the matrix is converted into an identity matrix.
5. The solution is obtained from the constant vector part of the matrix.

#### Gauss Jordan Code

```cpp
// Gauss Jordan Elimination Method implemented by 2207009
#include <bits/stdc++.h>
using namespace std;

void printMat(ofstream &out, vector<vector<double>> mat)
{
  int n = mat.size();
  int m = mat[0].size();
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < m; j++)
      out << mat[i][j] << " ";
    out << endl;
  }
}

void gaussJordanElimination(ofstream &out, vector<vector<double>> A, int n)
{
  int flag = 0;
  for (int i = 0; i < n; i++)
  {
    if (A[i][i] == 0)
    {
      for (int r = i + 1; r < n; r++)
      {
        if (A[r][i] != 0)
        {
          swap(A[i], A[r]);
          break;
        }
      }
      if (A[i][i] == 0)
      {
        flag = 1;
        break;
      }
    }

    double temp = A[i][i];
    for (int k = 0; k < n + 1; k++)
      A[i][k] = A[i][k] / temp;

    for (int j = i + 1; j < n; j++)
    {
      double ratio = A[j][i] / A[i][i];
      for (int k = 0; k < n + 1; k++)
      {
        A[j][k] = A[j][k] - ratio * A[i][k];
      }
    }

    for (int j = i - 1; j >= 0; j--)
    {
      double ratio = A[j][i] / A[i][i];
      for (int k = 0; k < n + 1; k++)
      {
        A[j][k] = A[j][k] - ratio * A[i][k];
      }
    }
  }

  if (flag == 1)
  {
    bool noSolution = false;
    bool infiniteSolution = false;

    for (int i = 0; i < n; i++)
    {
      bool allZero = true;
      for (int j = 0; j < n; j++)
      {
        if (A[i][j] > 0)
        {
          allZero = false;
          break;
        }
      }
      if (allZero && A[i][n] != 0)
        noSolution = true;
      if (allZero && A[i][n] == 0)
        infiniteSolution = true;
    }

    if (noSolution)
    {
      out << "No solution" << endl;
    }
    else if (infiniteSolution)
    {
      out << "Infinite solutions" << endl;
    }
  }
  else
  {
    out << "The Reduced Row Echelon Form is: " << endl;
    printMat(out, A);
    out << "The unique solution is: " << endl;
    for (int i = 0; i < n; i++)
      out << "x" << i + 1 << " = " << A[i][n] << (i == n - 1 ? "" : " , ");
    out << endl;
  }
}

int main(void)
{
  ifstream in("GaussJordanEliminationInput.txt");
  ofstream out("GaussJordanEliminationOutput.txt");

  int T;
  in >> T;

  while (T--)
  {
    int n;
    in >> n;
    vector<vector<double>> A(n, vector<double>(n + 1));

    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n + 1; j++)
        in >> A[i][j];
    }

    out << setprecision(3) << fixed;

    out << "The system of equations are: " << endl;
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        if (j > 0)
        {
          out << (A[i][j] < 0 ? " - " : " + ");
        }
        out << abs(A[i][j]) << "x" << j + 1;
      }
      out << " = " << A[i][n] << endl;
    }
    gaussJordanElimination(out, A, n);
  }

  in.close();
  out.close();

  return 0;
}
```

#### Gauss Jordan Input

```
3
5
2 1 -1 3 2 9
1 3 2 -1 1 8
3 2 4 1 -2 20
2 1 3 2 1 17
1 -1 2 3 4 15
2
1 1 2
2 2 4
2
1 1 2
2 2 5
```

##### Input Format

```
The input is taken from a file named GaussJordanEliminationInput.txt.

The first line of input contains an integer T - the number of test cases.

For each test case:

The first line contains an integer n - the number of equations.

The next n lines each contain n + 1 real numbers (augmented matrix).
```

#### Gauss Jordan Output

```
The system of equations are:
2.000x1 + 1.000x2 - 1.000x3 + 3.000x4 + 2.000x5 = 9.000
1.000x1 + 3.000x2 + 2.000x3 - 1.000x4 + 1.000x5 = 8.000
3.000x1 + 2.000x2 + 4.000x3 + 1.000x4 - 2.000x5 = 20.000
2.000x1 + 1.000x2 + 3.000x3 + 2.000x4 + 1.000x5 = 17.000
1.000x1 - 1.000x2 + 2.000x3 + 3.000x4 + 4.000x5 = 15.000
The Reduced Row Echelon Form is:
1.000 0.000 0.000 0.000 0.000 5.154
0.000 1.000 0.000 0.000 0.000 -1.000
0.000 0.000 1.000 0.000 0.000 2.262
0.000 0.000 0.000 1.000 0.000 -0.138
0.000 0.000 0.000 0.000 1.000 1.185
The unique solution is:
x1 = 5.154 , x2 = -1.000 , x3 = 2.262 , x4 = -0.138 , x5 = 1.185
The system of equations are:
1.000x1 + 1.000x2 = 2.000
2.000x1 + 2.000x2 = 4.000
Infinite solutions
The system of equations are:
1.000x1 + 1.000x2 = 2.000
2.000x1 + 2.000x2 = 5.000
No solution
```

##### Output Format

```
The output is written to a file named GaussJordanEliminationOutput.txt.

For each test case, the program prints:

The system of linear equations.

Nature of the solution: (Unique solution, No solution, Infinite solutions)

If the solution is unique, then print the Reduced Row Echelon Form Matrix.

The solution vector.

else If there is no solution, then print No solution.

else If there is infinite solutions, then print Infinite solutions.

All floating-point values are printed with 3 decimal places.
```

---

### LU Decomposition Method

Implemented by 2207008

#### LU Decomposition Theory

LU factorrization is a matrix decomposition technique used to solve a system of linear equations efficientlly. This method is also called Cholesky method, Doolittle's method, Crout's method. In this method, a square matrix A is decomposed into the product of two triangluar matrices. A = LU. Here L is a lower triangular matrix whose all diagonal elements are 1 and U is an upper triangular matrix. It simplifies the system AX = B by converting it into two simpler systems. AX = B --> LUX = B --> LY = B where (UX = Y). Then solves LY = B by forward substitution and UX = Y by backward substitution. This method is useful when solving system with the same coefficient matrix A but the constant vector B changes. <br>
Calculation steps:<br>

1. First row of U : u<sub>11</sub>, u<sub>12</sub>, u<sub>13</sub>, ...<br>
2. First column of L : l<sub>21</sub>, l<sub>31</sub>, ...<br>
3. Second row of U : u<sub>22</sub>, u<sub>23</sub>, ...<br>
4. Second column of L : l<sub>32</sub>, ...<br>
5. Third row of U : u<sub>33</sub>, ...

Solution Classification:<br>

1. After decomposition if U[i][i] is zero and Y[i] is non-zero then system has no solution.
2. If both are zero, the system has infinite solution.
3. Otherwise, the system has a unique solution.

#### LU Decomposition Code

```cpp
// LU Decomposition implemented by 2207008

#include <bits/stdc++.h>
using namespace std;

int main()
{
    ifstream in("LU_input.txt");
    ofstream out("LU_output.txt");

    if (!in || !out)
    {
        cout << "Error opening input/output file!" << endl;
        return 1;
    }

    out << fixed << setprecision(3);
    out << "LU Factorization method : " << endl;

    int t; // Number of test cases
    in >> t;
    for (int cs = 1; cs <= t; cs++)
    {
        int n;
        in >> n; // Number of equations.
        vector<vector<double>> A(n, vector<double>(n)), L(n, vector<double>(n, 0)), U(n, vector<double>(n, 0));
        vector<double> B(n), X(n), Y(n);

        // Entering the augmented matrix.
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                in >> A[i][j];
            }
            in >> B[i];
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = i; j < n; j++)
            {
                double sum = 0;
                for (int k = 0; k < i; k++)
                {
                    sum += L[i][k] * U[k][j];
                }
                U[i][j] = A[i][j] - sum;
            }
            for (int j = i; j < n; j++)
            {
                if (i == j)
                    L[i][j] = 1;
                else
                {
                    double sum = 0;
                    for (int k = 0; k < i; k++)
                    {
                        sum += L[j][k] * U[k][i];
                    }
                    L[j][i] = (A[j][i] - sum) / U[i][i];
                }
            }
        }
        out << "Output for case " << cs << " : " << endl;
        out << "\nLower Triangular Matrix (L) : " << endl;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                out << L[i][j] << " ";
            }
            out << endl;
        }
        out << "\nUpper Triangular Matrix (U) : " << endl;
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                out << U[i][j] << " ";
            }
            out << endl;
        }
        for (int i = 0; i < n; i++)
        {
            double sum = 0;
            for (int j = 0; j < i; j++)
            {
                sum += L[i][j] * Y[j];
            }
            Y[i] = B[i] - sum;
        }
        bool soln = true;
        for (int i = 0; i < n; i++)
        {
            if (fabs(U[i][i]) < 1e-9 && fabs(Y[i]) > 1e-9)
            {
                out << "\nNo Solution!" << endl;
                out << "\n------------------------------\n"
                    << endl;
                soln = false;
                break;
            }
            else if (fabs(U[i][i]) < 1e-9 && fabs(Y[i]) < 1e-9)
            {
                out << "\nInfinite Solution!" << endl;
                out << "\n------------------------------\n"
                    << endl;
                soln = false;
                break;
            }
        }
        if (soln == false)
            continue;
        for (int i = n - 1; i >= 0; i--)
        {
            double sum = 0;
            for (int j = i + 1; j < n; j++)
            {
                sum += U[i][j] * X[j];
            }
            X[i] = (Y[i] - sum) / U[i][i];
        }
        out << "\nSystem has unique solution!" << endl;
        out << "\nSolution : " << endl;
        for (int i = 0; i < n; i++)
        {
            out << "x" << i + 1 << " = " << X[i] << endl;
        }
        out << "\n------------------------------\n"
            << endl;
    }
    in.close();
    out.close();
    return 0;
}
```

#### LU Decomposition Input

```
3
5
2 1 -1 3 2 9
1 3 2 -1 1 8
3 2 4 1 -2 20
2 1 3 2 1 17
1 -1 2 3 4 15
2
1 1 2
2 2 4
2
1 1 2
2 2 5
```

##### Input Format

```
The input is taken from a file named LU_input.txt.

The first line contains an integer t — the number of test cases.

For each test case:

The first line contains an integer n — the number of equations (size of the system).

The next n lines each contain n + 1 real numbers (augmented matrix).
```

#### LU Decomposition Output

```
LU Factorization method :
Output for case 1 :

Lower Triangular Matrix (L) :
1.000 0.000 0.000 0.000 0.000
0.500 1.000 0.000 0.000 0.000
1.500 0.200 1.000 0.000 0.000
1.000 0.000 0.800 1.000 0.000
0.500 -0.600 0.800 1.714 1.000

Upper Triangular Matrix (U) :
2.000 1.000 -1.000 3.000 2.000
0.000 2.500 2.500 -2.500 0.000
0.000 0.000 5.000 -3.000 -5.000
0.000 0.000 0.000 1.400 3.000
0.000 0.000 0.000 0.000 1.857

System has unique solution!

Solution :
x1 = 5.154
x2 = -1.000
x3 = 2.262
x4 = -0.138
x5 = 1.185

------------------------------

Output for case 2 :

Lower Triangular Matrix (L) :
1.000 0.000
2.000 1.000

Upper Triangular Matrix (U) :
1.000 1.000
0.000 0.000

Infinite Solution!

------------------------------

Output for case 3 :

Lower Triangular Matrix (L) :
1.000 0.000
2.000 1.000

Upper Triangular Matrix (U) :
1.000 1.000
0.000 0.000

No Solution!

------------------------------


```

##### Output Format

```
The output is written to a file named LU_output.txt.

For each test case, the program prints:

The case number.

Lower Triangular Matrix L.

Upper Triangular Matrix U.

Nature of the solution: (No solution, Infinite soltuion, unique solution)

If the solution is unique, then print the solution vector.

A separator line after each test case.

All floating-point values are printed with 3 decimal places.
```

---

### Matrix Inversion Method

#### Matrix Inversion Theory

The Matrix Inversion Method is used to solve a system of linear equations when the coefficient matrix is square and invertible. In this method, the inverse of the coefficient matrix is calculated and then multiplied with the constant matrix to determine the solution.

Algorithm:

1. The coefficients of the system of linear equations are written into matrix form.
2. Check if the matrix is square and its determinant is non-zero.
3. Determine the cofactor matrix of the coefficient matrix.
4. Perform transpose operation to convert the cofactor matrix into adjoint matrix.
5. Divide all the elements of the adjoint matrix by the determinant of the coefficient matrix to convert the adjoint matrix into inverse matrix.
6. Multiply the inverse matrix with constant matrix to determine the solution of the system of equations.

#### Matrix Inversion Code

```cpp
// Matrix Inversion Method implemented by 2207009
#include <bits/stdc++.h>
using namespace std;

vector<vector<double>> cofactorMatrix(vector<vector<double>> A, int p, int q)
{
  int n = A.size();
  vector<vector<double>> temp;
  temp.reserve(n - 1);
  for (int i = 0; i < n; i++)
  {
    if (i == p)
      continue;
    vector<double> row;
    row.reserve(n - 1);
    for (int j = 0; j < n; j++)
    {
      if (j == q)
        continue;
      row.push_back(A[i][j]);
    }
    temp.push_back(row);
  }
  return temp;
}

double determinantMatrix(vector<vector<double>> A)
{
  int n = A.size();
  if (n == 1)
    return A[0][0];
  if (n == 2)
    return A[0][0] * A[1][1] - A[0][1] * A[1][0];

  double detA = 0.0;

  for (int j = 0; j < n; j++)
  {
    vector<vector<double>> cofactor = cofactorMatrix(A, 0, j);
    double value = determinantMatrix(cofactor);
    if (j % 2 == 1)
      value = -value;
    detA += value * A[0][j];
  }
  return detA;
}

vector<vector<double>> transposeMatrix(vector<vector<double>> A)
{
  int n = A.size();
  vector<vector<double>> T(n, vector<double>(n));
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
      T[j][i] = A[i][j];
  }
  return T;
}

vector<vector<double>> adjointMatrix(vector<vector<double>> A)
{
  int n = A.size();
  if (n == 1)
  {
    A[0][0] = 1;
    return A;
  };
  vector<vector<double>> adjoint(n, vector<double>(n));
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
    {
      vector<vector<double>> cofactor = cofactorMatrix(A, i, j);
      double value = determinantMatrix(cofactor);
      if ((i + j) % 2 == 1)
        value = -value;
      adjoint[i][j] = value;
    }
  }
  return transposeMatrix(adjoint);
}

vector<vector<double>> inverseMatrix(vector<vector<double>> A)
{
  double detA = determinantMatrix(A);
  if (detA == 0)
    throw runtime_error("Singular Matrix");
  int n = A.size();
  vector<vector<double>> adjoint = adjointMatrix(A);
  vector<vector<double>> inverse(n, vector<double>(n));
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
      inverse[i][j] = adjoint[i][j] / detA;
  }
  return inverse;
}

vector<double> calculateSolution(vector<vector<double>> invA, vector<double> B)
{
  int n = B.size();
  vector<double> X(n, 0);
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
      X[i] += invA[i][j] * B[j];
  }
  return X;
}

void printMat(ofstream &out, vector<vector<double>> A)
{
  int n = A.size();
  int m = A[0].size();
  for (int i = 0; i < n; i++)
  {
    for (int j = 0; j < n; j++)
      out << A[i][j] << " ";
    out << endl;
  }
}

int main(void)
{
  ifstream in("MatrixInversionInput.txt");
  ofstream out("MatrixInversionOutput.txt");

  int T;
  in >> T;

  while (T--)
  {

    int n;
    in >> n;

    vector<vector<double>> A(n, vector<double>(n));
    vector<double> B(n);

    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
        in >> A[i][j];
      in >> B[i];
    }

    out << setprecision(3) << fixed;

    out << "The system of equations are: " << endl;
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
      {
        if (j > 0)
        {
          out << (A[i][j] < 0 ? " - " : " + ");
        }
        out << abs(A[i][j]) << "x" << j + 1;
      }
      out << " = " << B[i] << endl;
    }

    double detA = determinantMatrix(A);

    if (detA > 0)
    {
      vector<vector<double>> transposeMat = adjointMatrix(A);
      out << "The transpose matrix is: " << endl;
      printMat(out, transposeMat);
      vector<vector<double>> invA = inverseMatrix(A);
      out << "The final inverse matrix is: " << endl;
      printMat(out, invA);
      vector<double> X = calculateSolution(invA, B);
      out << "The unique solution is: " << endl;
      for (int i = 0; i < n; i++)
        out << "x" << i + 1 << " = " << X[i] << (i == n - 1 ? "" : " , ");
      out << endl;
    }
    else
    {
      for (int i = 0; i < n; i++)
      {
        if (A[i][i] == 0)
        {
          for (int r = i + 1; r < n; r++)
          {
            if (A[r][i] != 0)
            {
              swap(A[i], A[r]);
              swap(B[i], B[r]);
              break;
            }
          }
          if (A[i][i] == 0)
            break;
        }

        for (int j = i + 1; j < n; j++)
        {
          double ratio = A[j][i] / A[i][i];
          for (int k = 0; k < n; k++)
            A[j][k] = A[j][k] - ratio * A[i][k];
          B[j] = B[j] - ratio * B[i];
        }
      }

      bool noSolution = false;
      bool infiniteSolution = false;
      for (int i = 0; i < n; i++)
      {
        bool allZero = true;
        for (int j = 0; j < n; j++)
        {
          if (A[i][j] > 0)
          {
            allZero = false;
            break;
          }
        }

        if (allZero && B[i] != 0)
          noSolution = true;
        if (allZero && B[i] == 0)
          infiniteSolution = true;
      }

      if (noSolution)
        out << "No solution" << endl;
      else if (infiniteSolution)
        out << "Infinite solutions" << endl;
    }
  }
  in.close();
  out.close();
  return 0;
}
```

#### Matrix Inversion Input

```
3
5
2 1 -1 3 2 9
1 3 2 -1 1 8
3 2 4 1 -2 20
2 1 3 2 1 17
1 -1 2 3 4 15
2
1 1 2
2 2 4
2
1 1 2
2 2 5
```

##### Input Format

```
The input is taken from a file named MatrixInversionInput.txt.

The first line of input contains an integer T - the number of test cases.

For each test case:

The first line contains an integer n - the number of equations.

The next n lines each contain n + 1 real numbers (augmented matrix).
```

#### Matrix Inversion Output

```
The system of equations are:
2.000x1 + 1.000x2 - 1.000x3 + 3.000x4 + 2.000x5 = 9.000
1.000x1 + 3.000x2 + 2.000x3 - 1.000x4 + 1.000x5 = 8.000
3.000x1 + 2.000x2 + 4.000x3 + 1.000x4 - 2.000x5 = 20.000
2.000x1 + 1.000x2 + 3.000x3 + 2.000x4 + 1.000x5 = 17.000
1.000x1 - 1.000x2 + 2.000x3 + 3.000x4 + 4.000x5 = 15.000
The transpose matrix is:
25.000 25.000 125.000 -245.000 105.000
-0.000 0.000 -65.000 130.000 -65.000
-16.000 -3.000 -15.000 45.000 -10.000
-3.000 -29.000 -80.000 175.000 -75.000
4.000 17.000 20.000 -60.000 35.000
The final inverse matrix is:
0.385 0.385 1.923 -3.769 1.615
-0.000 0.000 -1.000 2.000 -1.000
-0.246 -0.046 -0.231 0.692 -0.154
-0.046 -0.446 -1.231 2.692 -1.154
0.062 0.262 0.308 -0.923 0.538
The unique solution is:
x1 = 5.154 , x2 = -1.000 , x3 = 2.262 , x4 = -0.138 , x5 = 1.185
The system of equations are:
1.000x1 + 1.000x2 = 2.000
2.000x1 + 2.000x2 = 4.000
Infinite solutions
The system of equations are:
1.000x1 + 1.000x2 = 2.000
2.000x1 + 2.000x2 = 5.000
No solution
```

##### Output Format

```
The output is written to a file named MatrixInversionOutput.txt.

For each test case, the program prints:

The system of linear equations.

The Transpose Matrix.

The Inverse Matrix.

Nature of the solution: (Unique solution, No solution, Infinite solutions)

If the solution is unique, then print the solution vector.

else If there is no solution, then print No solution.

else If there is infinite solutions, then print Infinite solutions.

All floating-point values are printed with 3 decimal places.
```

---

### Solution of Non-Linear Equations

### Bisection Method

Implemented by 2207008

#### Bisection Theory

The Bisection Method is a numerical method used to approximate a real root of a nonlinear equation of the form: f(x) = 0. It is also called Binary chopping or half-interval method. It is a braketing method which states that if a continuous function f(x) gives opposite signs at two points a and b, then at least one real root exists in the interval (a,b).<br>

Algorithm:

1. Let x1 = a and x2 = b. Define x0 as the midpoint between a and b. x0 = (x1 + x2) / 2.
2. Evaluate f(x0), f(x1), f(x2).
3. If f(x0) = 0, then the root is x0.
4. Select the subinterval where the sign change occurs based on the following conditions:<br>
   If f(x0) \* f(x1) < 0, then root is between x0 and x1. So take [x0, x1] as the new subinterval.<br>
   If f(x0) \* f(x2) < 0, then root is between x0 and x2. So take [x0, x2] as the new subinterval.
5. Repeat the step 3 while fabs(f(x0)) > E. (Solution is close enough to zero).

#### Bisection Code

```cpp
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
```

#### Bisection Input

```
4
0.5 0.0001
1 0 -5 0 4
```

##### Input Format

```
The input is read from a file named BiSection_input.txt.

The first line contains an integer n (Degree of the polynomial equation).

The second line contains two real numbers: stepSize (step size for scanning the interval) and e (tolerance or error limit).

The third line contains n + 1 real numbers : Coefficients of the polynomial in descending order of degree.
```

#### Bisection Output

```
Bisection method :
The inputted equation is : 1x^4 - 5x^2 + 4 = 0
Bracket for root 1 : [-2.162, -1.662]
Iteration needed : 15
Root 1 : -2.000

Bracket for root 2 : [-1.162, -0.662]
Iteration needed : 13
Root 2 : -1.000

Bracket for root 3 : [0.838, 1.338]
Iteration needed : 13
Root 3 : 1.000

Bracket for root 4 : [1.838, 2.338]
Iteration needed : 15
Root 4 : 2.000


```

##### Output Format:

```
The output is written to a file named BiSection_output.txt.

The name of method is printed.

The given polynomial equation is printed at first.

For each detected root, the program prints:

The bracketing interval where a root exists

Number of iterations required to approximate the root.

The approximated root value

A blank line after each root

All floating-point values are printed with 3 decimal places.
```

---

### False Position Method

Implemented by 2207007

#### False Position Theory

The False Position Method(regular falsi) is a numerical method that is used to find the real root of a nonlinear equation `f(x) = 0`.
Like the Bisection Method, it is a bracketing method, which needs two points a and b, such that `f(a) · f(b) < 0` (different signs). Instead of finding the midpoint, it uses a straight line between `(a, f(a))` and `(b, f(b))` and the point of intersection with the x-axis.
Algorithm:

1. Select initial points a and b such that `f(a) · f(b) < 0`
2. Calculate c:

   `c = b - f(b) · (b - a) / (f(b) - f(a))`

3. Calculate f(c)
4. Check for convergence: if `|f(c)| < Tolerance(epsilon)`, then c is the root
5. Update the bracket based on the sign of f(c):

If `f(a) · f(c) > 0`, then a = c
Else, b = c

Repeat steps 2-5 until converged tolerance or the maximum number of iterations is reached.

#### False Position Code

```cpp
#include <bits/stdc++.h>
using namespace std;

double f(double x, vector<double> &v) {
    double result = 0;
    int n = v.size();
    for (int i = 0; i < n; i++) {
        int deg = n - i - 1;
        result += v[i] * pow(x, deg);
    }
    return result;
}

// False Position for one bracket
pair<double, int> falsePosition(double a, double b, vector<double> &eqn, double tol, int maxIt = 1000) {
    double fa = f(a, eqn);
    double fb = f(b, eqn);
    if (fa * fb >= 0) return {0, 0};

    double c = a;
    int steps = 0;
    while (steps < maxIt) {
        c = b - fb * (b - a) / (fb - fa);
        double fc = f(c, eqn);

        steps++; // increment iteration

        if (fabs(fc) < tol) break;

        if (fa * fc < 0) {
            b = c;
            fb = fc;
        } else {
            a = c;
            fa = fc;
        }
    }
    return {c, steps};
}

int main() {
    ifstream in("F_input.txt");
    ofstream out("F_output.txt");

    if (!in || !out) return 1;

    int n;
    in >> n;
    double stepSize, tol;
    in >> stepSize >> tol;

    vector<double> eqn(n + 1);
    for (int i = 0; i <= n; i++) in >> eqn[i];

    out << "False Position Method\n";
    out << "Equation: ";
    for (int i = 0; i <= n; i++) {
        if (eqn[i] == 0) continue;
        if (i != 0 && eqn[i] > 0) out << " + ";
        if (eqn[i] < 0) out << " - ";
        out << abs(eqn[i]);
        int deg = n - i;
        if (deg > 0) out << "x";
        if (deg > 1) out << "^" << deg;
    }
    out << " = 0\n\n";

    // Fixed scanning range
    double scanStart = -10, scanEnd = 10;
    vector<pair<double, double>> intervals;

    double prev = scanStart;
    for (double cur = scanStart + stepSize; cur <= scanEnd; cur += stepSize) {
        if (f(prev, eqn) * f(cur, eqn) < 0)
            intervals.push_back({prev, cur});
        prev = cur;
    }

    if (intervals.empty()) {
        out << "No brackets found in this range\n";
        return 0;
    }

    out << left << setw(10) << "Root"
        << setw(12) << "Iterations"
        << setw(20) << "Bracket" << endl;
    out << "----------------------------------------\n";

    for (int i = 0; i < intervals.size(); i++) {
        auto res = falsePosition(intervals[i].first, intervals[i].second, eqn, tol);
        out << setw(10) << fixed << setprecision(6) << res.first
            << setw(12) << res.second
            << "[" << intervals[i].first << ", " << intervals[i].second << "]" << endl;
    }

    cout << "Output written to FB_output.txt\n";
    return 0;
}

```

#### False Position Input

```txt
4
0.1 0.0001
1 10 -5 0 4
```

##### Input Format:

```
The input is read from a file named F_input.txt.

The inputs are taken as

deg  -> Degree of the polynomial equation

a0 a1 a2 ... an -> Coefficients of the polynomial in descending order

step -> Step size for scanning the interval to find initial brackets

eps -> Convergence tolerance for the root
```

#### False Position Output

```txt
False Position Method
Equation: 1x^4 + 10x^3 - 5x^2 + 4 = 0

Root      Iterations  Bracket
----------------------------------------
-0.610497 4           [-0.700000, -0.600000]

```

##### Output Format

```
The output is written to a file named F_output.txt.

The polynomial equation is printed.

In a table Root, Iterations, Bracket are printed for each detected root:
Root      Iter   Bracket
------------------------------
x1        iter1  [l1, r1]
x2        iter2  [l2, r2]

```

---

### Secant Method

Implemented by 2207007

#### Secant Theory

The Secant Method is a numerical technique for finding roots of equations f(x) = 0. It's similar to Newton-Raphson but doesn't require the derivative. Instead, it approximates the derivative using two previous points.

The method draws a secant line through two points on the curve and finds where it crosses the x-axis. This intersection becomes the next approximation.

Algorithm:

1. Choose two initial guesses x₀ and x₁
2. Calculate the next approximation:
   `x_{n+1} = x_n - f(x_n) · (x_n - x_{n-1}) / (f(x_n) - f(x_{n-1}))`
3. Evaluate `f(x_{n+1})`
4. Check convergence: if `|x_{n+1} - x_n| < epsilon`, then `x_{n+1}` is the root
5. Update: `x_{n-1} = x_n` and `x_n = x_{n+1}`
6. Repeat steps 2-5 until convergence or maximum iterations reached

The Secant Method converges slower than Newton-Raphson but faster than bisection. Its main advantage is that it doesn't need the derivative, making it useful when `f'(x)` is difficult to compute or doesn't exist.

#### Secant Code

```cpp
#include <bits/stdc++.h>
using namespace std;

double f(double x) {
    return x*x - 4.0;
}

double df(double x) {
    return 2.0*x;
}
double secant(double x_n, double x_n1, double epsilon, int max_iteration) {
    double x_n2;
    for (int i = 0; i < max_iteration; i++) {
        x_n2 = x_n1 - ((x_n1 - x_n) / (f(x_n1) - f(x_n))) * f(x_n1);
        if (abs(x_n2 - x_n1) < epsilon) {
            return x_n2;
        }
        x_n = x_n1;
        x_n1 = x_n2;
    }
    return NAN;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    double range_start, range_end, step, epsilon;
    int max_iteration;
    fin >> range_start >> range_end >> step >> epsilon >> max_iteration;

    vector<double> raphson_points;
    vector<pair<double, double>> secant_pairs;

    for (double i = range_start; i < range_end; i += step) {
        if (f(i) * f(i + step) < 0) {
            raphson_points.push_back(i + step / 2);
            secant_pairs.push_back({i, i + step});
        }
    }
    if (raphson_points.size() == 0) {
        raphson_points.push_back(0.5);
    }

    for (auto pr : secant_pairs) {
        double root = secant(pr.first, pr.second, epsilon, max_iteration);
        fout << "Root using Secant: " << root << endl;
    }

    fin.close();
    fout.close();
    return 0;
}

```

#### Secant Input

```
-5 5 0.1 1e-7 100000
```

##### Input Format

```
The input is read from a file named input.txt. It contains:

range_start → Start of the interval to scan for roots

range_end → End of the interval

step → Step size for scanning the interval to detect brackets

epsilon → Tolerance for the root

max_iteration → Maximum number of iterations allowed
```

#### Secant Output

```
Root using Secant: -2
Root using Secant: 2
```

##### Output Format

```
The output is written to a file named output.txt.

For each detected interval where the function changes sign, the program prints the root calculated using the Secant Method:

Root using Secant: value
```

---

### Newton-Raphson Method

Implemented by 2207007

#### Newton-Raphson Theory

The Newton-Raphson Method is a numerical technique for finding roots of equations `f(x) = 0`. It uses the derivative of the function to improve initial guess, converging faster to the solution.

The method draws a tangent line at the current point and finds where it crosses the x-axis. This intersection becomes the next approximation.

Algorithm:

1. Choose an initial guess x₀
2. Calculate the next approximation:
   `x_{n+1} = x_n - f(x_n) / f'(x_n)`
3. Evaluate `f(x_{n+1})`
4. Check convergence: if `|x_{n+1} - x_n| < epsilon`, then `x_{n+1}` is the root
5. Update: `x_n = x_{n+1}`
6. Repeat steps 2-5 until convergence or maximum iterations reached

The method converges rapidly when the initial guess is close to the actual root. However, it requires the derivative `f'(x)` and may fail if the derivative is zero or the starting point is too far from the root.

#### Newton-Raphson Code

```cpp
#include <bits/stdc++.h>
using namespace std;

double f(double x) {
    return x*x - 4.0;
}

double df(double x) {
    return 2.0*x;
}

double NewtonRaphson(double x_n, double epsilon, int max_iteration) {
    double x_n1;
    for (int i = 0; i < max_iteration; i++) {
        x_n1 = x_n - (f(x_n) / df(x_n));
        if (abs(x_n1 - x_n) < epsilon) {
            return x_n1;
        }
        x_n = x_n1;
    }
    return NAN;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    double range_start, range_end, step, epsilon;
    int max_iteration;
    fin >> range_start >> range_end >> step >> epsilon >> max_iteration;

    vector<double> raphson_points;
    vector<pair<double, double>> secant_pairs;

    for (double i = range_start; i < range_end; i += step) {
        if (f(i) * f(i + step) < 0) {
            raphson_points.push_back(i + step / 2);
            secant_pairs.push_back({i, i + step});
        }
    }
    if (raphson_points.size() == 0) {
        raphson_points.push_back(0.5);
    }

    for (double pt : raphson_points) {
        double root = NewtonRaphson(pt, epsilon, max_iteration);
        fout << "Root using Newton Raphson: " << root << endl;
    }

    fin.close();
    fout.close();
    return 0;
}
```

#### Newton-Raphson Input

```
-5 5 0.1 1e-7 100000
```

##### Input Format

```
The input is read from a file named input.txt. It contains:

range_start -> Start of the interval to scan for roots
range_end -> End of the interval
step -> Step size for scanning the interval
epsilon -> Tolerance for the root
max_iteration -> Maximum number of iterations allowed
```

#### Newton-Raphson Output

```
Root using Newton Raphson: -2
Root using Newton Raphson: 2
```

##### Output Format

```
The output is written to a file named output.txt.

For each initial guess detected within the interval, the program prints the root calculated using the Newton-Raphson method:

Root using Newton Raphson: value
```

---

### Solution of Differential Equations

### Runge-Kutta Method

Implemented by 2207008

#### Runge-Kutta Theory

The Runge-Kutta 4th Order (RK4) method is a numerical technique to approximate the solution of a first-order ordinary differential equation (ODE) which is of the the form: dy/dx = f(x, y), y(x<sub>0</sub>) = y<sub>0</sub>
. It provies a good balance between accuracy and computational efficiency.<br>

Algorithm:

1. Given a initial value of x --> x<sub>0</sub> and y --> y<sub>0</sub>.
2. Given a step size h, the value of y at the next point is computed using four intermediate slopes: <br>
   k<sub>1</sub> = h _ f(x<sub>n</sub>, y<sub>n</sub>)<br>
   k<sub>2</sub> = h _ f(x<sub>n</sub> + h/2, y<sub>n</sub> + k<sub>1</sub>/2)<br>
   k<sub>3</sub> = h _ f(x<sub>n</sub> + h/2, y<sub>n</sub> + k<sub>2</sub>/2)<br>
   k<sub>4</sub> = h _ f(x<sub>n</sub> + h, y<sub>n</sub> + k<sub>3</sub>)

3. The next value of y is then calculated as:<br>
   y<sub>n+1</sub> = y<sub>n</sub> + (k<sub>1</sub> + 2 _ k<sub>2</sub> + 2 _ k<sub>3</sub> + k<sub>4</sub>) / 6

#### Runge-Kutta Code

```cpp
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
        in.close();
        out.close();
        return 1;
    }

    double x0, y0, xn, yn, h;
    in >> x0 >> y0; // Initial values x0 and y0
    in >> xn;       // Final values of x
    in >> h;        // Step size

    if (h <= 0)
    {
        cout << "Step size h must be positive!" << endl;
        in.close();
        out.close();
        return 1;
    }

    if (xn <= x0)
    {
        cout << "xn must be greater than x0!" << endl;
        in.close();
        out.close();
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

```

#### Runge-Kutta Input

```
0 1
0.2
0.1
```

##### Input Format

```
The input is read from a file named RK_input.txt.

The first line contains two real numbers: x0 (initial value of the independent variable) and y0 (initial value of the dependent variable)

The second line contains a real number: xn (final value of x)

The third line contains a real number: h (step size)
```

#### Runge-Kutta Output

```
RK4 Method:
Step 1: x = 0.100  y = 1.110
Step 2: x = 0.200  y = 1.240

Final Result: yn = 1.240

```

##### Output Format

```
The output is written to a file named RK_output.txt.

At first the method name is printed.

For each step, the updated values of x and y are displayed.

Finally, the computed value of y at x = xn is printed.

All numerical values are printed with 3 decimal places.
```

---

### Numerical Integration

### Simpson's one-third rule

#### Simpson's one-third rule Theory

Simpson's 1/3 Rule is a numerical integration method used to approximate the definite integral of a function. It works by fitting parabolas (quadratic curves) through sets of three points and calculating the area under these parabolas.

The method divides the interval `[a, b]` into an even number of subintervals and approximates the curve with parabolic segments, providing better accuracy than simpler methods like the Trapezoidal Rule.

Algorithm:

1. Divide the interval `[a, b]` into n equal subintervals (n must be even)
2. Calculate the step size:
   `h = (b - a) / n`
3. Generate points:
   `x_i = a + i · h` for `i = 0, 1, 2, ..., n`
4. Evaluate the function at each point:
   `f(x_0), f(x_1), f(x_2), ..., f(x_n)`
5. Apply Simpson's 1/3 formula:
   `I ≈ (h/3) · [f(x_0) + 4·Σf(x_odd) + 2·Σf(x_even) + f(x_n)]`

   Where:

   - `f(x_0)` and `f(x_n)` are the first and last terms
   - `Σf(x_odd)` is the sum of function values at odd indices
   - `Σf(x_even)` is the sum of function values at even indices (excluding endpoints)

Simpson's 1/3 Rule provides excellent accuracy for smooth functions and requires an even number of intervals to work properly.

#### Simpson's one-third rule Code

```cpp
#include <bits/stdc++.h>
using namespace std;

#define vd vector<double>

ofstream fout;

double func(double x)
{
    return sqrt(x);
}

void simpsonOneThird(double a, double b, double h)
{
    fout<<"Simpson's 1/3:\n";
    int n = (b - a) / h;
    vd x(n + 1), fx(n + 1);

    for (int i = 0; i <= n; i++)
    {
        x[i] = a + i * h;
        fx[i] = func(x[i]);
    }

    double res = 0;
    for (int i = 0; i <= n; i++)
    {
        if (i == 0 || i == n)
            res += fx[i];
        else if (i % 2 != 0)
            res += 4 * fx[i];
        else
            res += 2 * fx[i];
    }

    fout << res * (h / 3);
}

int main()
{
    // input format
    //a b h
    ifstream fin("input.txt");
    fout.open("output.txt");

    double a, b, n;
    fin >> a >> b >> n;
    double h = (b-a)/n;

    simpsonOneThird(a, b, h);
    cout<<"Output in output.txt";
    fin.close();
    fout.close();
    return 0;
}

```

#### Simpson's one-third rule Input

```
0 4.5 90
```

##### Input Format

```
a -> lower limit
b -> upper limit
n -> total steps(slices)
```

#### Simpson's one-third rule Output

```
Simpson's 1/3:
Result = 6.22222
```

##### Output Format

```
The output is printed rounded with 2 decimal numbers after the floating point
```

### Simpson's three-eighths rule rule

#### Simpson's three-eighths rule Theory

Simpson's 3/8 Rule Theory

Simpson's 3/8 Rule is a numerical integration method used to approximate the definite integral of a function. It works by fitting cubic curves (third-degree polynomials) through sets of four points and calculating the area under these curves.

The method divides the interval `[a, b]` into subintervals where the number of intervals must be a multiple of 3, providing higher accuracy than Simpson's 1/3 Rule for certain functions.

Algorithm:

1. Divide the interval `[a, b]` into n equal subintervals (n must be a multiple of 3)

2. Calculate the step size:
   `h = (b - a) / n`

3. Generate points:
   `x_i = a + i · h` for `i = 0, 1, 2, ..., n`

4. Evaluate the function at each point:
   `f(x_0), f(x_1), f(x_2), ..., f(x_n)`

5. Apply Simpson's 3/8 formula:
   `I ≈ (3h/8) · [f(x_0) + 3·Σf(x_non-multiple-of-3) + 2·Σf(x_multiple-of-3) + f(x_n)]`

   Where:

   - `f(x_0)` and `f(x_n)` are the first and last terms
   - `Σf(x_non-multiple-of-3)` is the sum of function values at indices not divisible by 3
   - `Σf(x_multiple-of-3)` is the sum of function values at indices divisible by 3 (excluding endpoints)

Simpson's 3/8 Rule is particularly useful when the number of intervals is a multiple of 3 and provides better accuracy for functions with higher-order derivatives.

#### Simpson's three-eighths rule Code

```cpp
#include <bits/stdc++.h>
using namespace std;

#define vd vector<double>

ofstream fout;

double func(double x)
{
    return sqrt(x);
}


void simpsonThreeEight(double a, double b, double h)
{
    fout<<"Simpson's 3/8:\n";
    int n = round((b - a) / h);
    if (n % 3 != 0)
    {
        fout << "Error: n must be multiple of 3 for Simpson 3/8 rule\n";
        return;
    }

    h = (b - a) / n;

    vd x(n + 1), fx(n + 1);
    for (int i = 0; i <= n; i++)
    {
        x[i] = a + i * h;
        fx[i] = func(x[i]);
    }

    double res = fx[0] + fx[n];
    for (int i = 1; i < n; i++)
    {
        if (i % 3 == 0)
            res += 2 * fx[i];
        else
            res += 3 * fx[i];
    }

    fout << res * (3 * h / 8);
}

int main()
{
    // input format
    //a b h
    ifstream fin("input.txt");
    fout.open("output.txt");

    double a, b, n;
    fin >> a >> b >> n;
    double h = (b-a)/n;

    simpsonThreeEight(a, b, h);
    cout<<"Output in output.txt";
    fin.close();
    fout.close();
    return 0;
}

```

#### Simpson's three-eighths rule Input

```
0 4.5 90
```

##### Input Format

```
a -> lower limit
b -> upper limit
n -> total steps(slices)
```

#### Simpson's three-eighths rule Output

```
Simpson's 3/8:
6.36
```

##### Output Format

```
The output is printed rounded with 2 decimal numbers after the floating point
```

---

### Numerical Differentiation

### Differentiation

#### Differentiation Theory

Forward Differentiation is a numerical method used to approximate the derivative of a function at a given point using forward difference formulas. It constructs a difference table from function values and applies formulas to estimate first and second derivatives.

The method uses equally spaced points and builds a table of forward differences, then applies Newton's forward difference formulas to calculate derivatives.

Algorithm:

1. Divide the interval `[a, b]` into n equal subintervals
2. Calculate the step size:
   `h = (b - a) / n`
3. Generate points and evaluate the function:
   `x_i = a + i · h` and `y_i = f(x_i)` for `i = 0, 1, 2, ..., n`
4. Construct the forward difference table:
   - `Δ⁰y_i = y_i`
   - `Δʲy_i = Δʲ⁻¹y_{i+1} - Δʲ⁻¹y_i` for `j = 1, 2, ..., n`
5. Calculate the normalized distance from the first point:
   `u = (p - x_0) / h`
6. Apply Newton's forward difference formula for first derivative:
   `f'(p) = [Δy_0 + ((2u-1)/2)·Δ²y_0 + ((3u²-6u+2)/6)·Δ³y_0] / h`
7. Apply formula for second derivative:
   `f''(p) = [Δ²y_0 + (u-1)·Δ³y_0] / h²`

Forward differentiation works best when the point p is near the beginning of the interval. The accuracy depends on the step size h and the smoothness of the function.

#### Differentiation Code

```cpp
#include <bits/stdc++.h>
using namespace std;

double func(double x) {
    return x*x*x - 2*x*x + x;
}

void forwardDifferentiation(double a, double b, int n, double p, ofstream &fout) {
    double h = (b - a) / n;

    vector<double> X(n+1), Y(n+1);
    for (int i = 0; i <= n; i++) {
        X[i] = a + i*h;
        Y[i] = func(X[i]);
    }

    // Forward difference table
    vector<vector<double>> d(n+1, vector<double>(n+1, 0));
    for (int i = 0; i <= n; i++)
        d[i][0] = Y[i];

    for (int j = 1; j <= n; j++)
        for (int i = 0; i <= n - j; i++)
            d[i][j] = d[i+1][j-1] - d[i][j-1];

    fout << "Forward Difference Table:\n";
    for (int i = 0; i <= n; i++) {
        for (int j = 0; j <= n-i; j++)
            fout << d[i][j] << "\t";
        fout << "\n";
    }

    double u = (p - X[0]) / h;

    double fp = 0, fpp = 0;

    if (n >= 2) {
        fp = (d[0][1] + ((2*u - 1)/2.0)*d[0][2] + ((3*u*u - 6*u + 2)/6.0)*d[0][3]) / h;
        fpp = (d[0][2] + (u - 1)*d[0][3]) / (h*h);
    } else if (n == 1) {
        fp = d[0][1] / h;
        fpp = 0;
    }

    fout << "\nf'(p) = " << fp << endl;
    fout << "f''(p) = " << fpp << endl;
}

int main() {
    ifstream fin("input.txt");
    ofstream fout("output.txt");

    double a, b, p;
    int n;

    fin >> a >> b >> n >> p;

    forwardDifferentiation(a, b, n, p, fout);

    fin.close();
    fout.close();
    return 0;
}

```

#### Differentiation Input

```
0 2 4 1.5
```

##### Input Format

```
a -> starting point of the interval
b -> ending point of the interval
n -> number of subintervals (determines how many points to use)
p -> the specific point where you want to calculate the derivative
```

#### Differentiation Output

```
Forward Difference Table:
0	0.125	-0.25	0.75	0
0.125	-0.125	0.5	0.75
0	0.375	1.25
0.375	1.625
2

f'(p) = 1.75
f''(p) = 5

```

##### Output Format

```
The forward Difference table is printed and then the first and second order derivative are calculated
```

---

### Curve Fitting Regression

### Linear Equation

#### Linear Equation Theory

Linear equation can be solved using curve fitting which is also known as linear regression. It is a numerical technique to find the best fit straight line that represents the relationship between two variables based on a given set of data points. The relationship between the variables is expressed in the form y=a+bx. Where a is the intercept and 𝑏 is the slope of the line. The values of a and b are determined using the least squares principle. When the best fit equation is obtained, it can be used to predict the value of the dependent variable (y) for any given value of the independent variable (x).

​
Formula:

let, n is the number of data points. Then,

b = n∑xy−(∑x)(∑y)​ / n∑x^2−(∑x)^2

a = ∑y−b∑x / n

The best fit linear equation is:

y = a + bx

#### Linear Equation Code

```cpp
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
```

#### Linear Equation Input

```
7
1 3
2 4
3 4
4 5
5 8
6 9
7 10
8
```

##### Input Format

```
The input is taken from a file named LinearEquationInput.txt.

The first line of input contains an integer n - the number of data points.

The next n lines each contain two real numbers x and y - the coordinates of the data points.

The last line of input contains an integer newX - for which the predicted value of y is to be calculated.
```

#### Linear Equation Output

```
Data points entered:
(1.00, 3.00)
(2.00, 4.00)
(3.00, 4.00)
(4.00, 5.00)
(5.00, 8.00)
(6.00, 9.00)
(7.00, 10.00)
The Linear Equation is:
y = 1.14 + 1.25x
For x = 8, the predicted value of y is: 11.14
```

##### Output Format

```
The output is written to a file named LinearEquationOutput.txt.

All the entered data points are printed.

Then the equation of the best fit straight line is printed in the form
y = a + bx.

Finally, the predicted value of y corresponding to x = newX is printed.
```

### Polynomial Equation

#### Polynomial Equation Theory

Polynomial equation can be solved using curve fitting. It is a numerical technique to find the best fit polynomial curve that represents the relationship between two variables based on a given set of data points. In this method, the relationship between the variables is assumed to be a polynomial of degree m and is expressed in the form: y=a0​+a1​x^1+a2​x^2+⋯+am​x^m. Here, a0,a1,a2,...,am are the coefficients of the polynomial. These coefficients are determined using the least squares principle. After the best fit polynomial equation is obtained, it can be used to predict the value of the dependent variable (y) for any given value of the independent variable (x).

Formula:

For each coefficient `ai` (i = 0, 1, ..., m), the normal equation is:

Σ(y_k \* x_k^i) = a0Σ(x_k^i) + a1Σ(x_k^(i+1)) + a2Σ(x_k^(i+2)) + ... + amΣ(x_k^(i+m))

The normal equations are written in matrix form.

The best fit polynomial equation is:

y = a0 + a1x + a2x^2 + ... + amx^m

Algorithm:

1. Create a set of equations using the least squares method based on the given data points.

2. Convert the equations in matrix form.

3. Solve the matrix using Gaussian elimination to find the coefficients of the polynomial.

4. Determine the best fit polynomial equation using the calculated coefficients.

5. Substitute the given value into the polynomial equation to calculate the predicted result.

#### Polynomial Equation Code

```cpp
#include <bits/stdc++.h>

using namespace std;

void printPolynomial(ofstream &out, vector<double> coeff, int m)
{
  out << "y = ";
  bool firstTerm = true;

  for (int i = 0; i < m; i++)
  {
    if (coeff[i] == 0)
      continue;

    if (!firstTerm)
      out << " + ";

    if (i == 0)
    {
      out << coeff[i];
    }
    else if (i == 1)
    {
      if (coeff[i] == 1)
        out << "x";
      else
        out << coeff[i] << "x";
    }
    else
    {
      if (coeff[i] == 1)
        out << "x^" << i;
      else
        out << coeff[i] << "x^" << i;
    }
    firstTerm = false;
  }
  out << endl;
}

void Polynomial(ofstream &out, vector<double> x, vector<double> y, int m, double newX, int n)
{
  int size = m + 1;
  vector<vector<double>> A(size, vector<double>(size, 0));
  vector<double> B(size, 0);

  for (int i = 0; i < size; i++)
  {
    for (int j = 0; j < size; j++)
    {
      for (int k = 0; k < n; k++)
      {
        A[i][j] += pow(x[k], i + j);
      }
    }
    for (int k = 0; k < n; k++)
    {
      B[i] += y[k] * pow(x[k], i);
    }
  }

  vector<double> coeff(size, 0);
  for (int i = 0; i < size; i++)
  {
    int maxi = i;
    for (int k = i + 1; k < size; k++)
    {
      if (fabs(A[k][i]) > fabs(A[maxi][i]))
      {
        maxi = k;
      }
    }
    if (A[maxi][i] == 0)
    {
      out << "Matrix is singular." << endl;
      return;
    }
    if (maxi != i)
    {
      swap(A[i], A[maxi]);
      swap(B[i], B[maxi]);
    }

    for (int j = i + 1; j < size; j++)
    {
      double factor = A[j][i] / A[i][i];
      for (int k = i; k < size; k++)
      {
        A[j][k] -= factor * A[i][k];
      }
      B[j] -= factor * B[i];
    }
  }
  for (int i = size - 1; i >= 0; i--)
  {
    coeff[i] = B[i];
    for (int j = i + 1; j < size; j++)
    {
      coeff[i] -= A[i][j] * coeff[j];
    }
    coeff[i] /= A[i][i];
  }

  printPolynomial(out, coeff, m);

  double newY = 0.0;
  for (int i = 0; i < size; i++)
  {
    newY += coeff[i] * pow(newX, i);
  }
  out << "Predicted value at x = " << newX << " is y = " << newY << endl;
}

int main(void)
{
  ifstream in("PolynomialEquationInput.txt");
  ofstream out("PolynomialEquationOutput.txt");

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

  int m;
  in >> m;

  double newX;
  in >> newX;

  Polynomial(out, x, y, m, newX, n);

  in.close();
  out.close();
  return 0;
}
```

#### Polynomial Equation Input

```
5
1 6
2 11
3 18
4 27
5 38
3
6
```

##### Input Format

```
The input is taken from a file named PolynomialEquationInput.txt.

The first line of input contains an integer n - the number of data points.

The next n lines each contain two real numbers x and y - the coordinates of the data points.

The next line contains an integer m - the degree of the polynomial.

The last line contains a real number newX, for which the predicted value of y is to be calculated.
```

#### Polynomial Equation Output

```
Data points entered:
(1.00, 6.00)
(2.00, 11.00)
(3.00, 18.00)
(4.00, 27.00)
(5.00, 38.00)
y = 3.00 + 2.00x + 1.00x^2
For x = 6.00, predicted value of y is: 51.00
```

##### Output Format

```
The output is written to a file named PolynomialEquationOutput.txt.

All the entered data points are printed.

Then the best fit polynomial equation is printed.

Finally, the predicted value of y for x = newX is printed.
```

### Transcendental Equation

#### Transcendental Equation Theory

Transcendental equation can be solved using curve fitting which is also known as power regression. It is a numerical technique used to find the best fit curve when the relationship between two variables is non-linear and follows a power law. The relationship between the variables is expressed in the form y=ax^b, where a and b are constants. To apply the least squares principle, the equation is first converted into a linear form by taking logarithm on both sides, which gives lny=lna+blnx. Assuming lny as Y and lnx as X. Then, the equation becomes a straight line equation. The values of a and b are then determined using the least squares method. After that, the best fit curve equation is obtained. Then, it can be used to predict the value of the dependent variable (y) for any given value of the independent variable (x).

Formula:

Let, n is the number of data points.

Given equation is, y=ax^b

taking logarithm on both sides:

lny=lna+blnx

Let,

Y=lny, X=lnx, A=lna

Now,
​

b = n∑XY−(∑X)(∑Y) / n∑X^2 - (∑X)^2

A = ∑Y−b∑X / n

a = e^A

The best fit transcendental equation is:

y = ax^b

#### Transcendental Equation Code

```cpp
#include <bits/stdc++.h>

using namespace std;

// The transcendental equation format is: y=ax^b
// lny = lna + blnx (Taking ln on both sides)
// Y = A + bX
// Here, Y = lny, X = lnx, A = lna
void Transcendental(ofstream &out, vector<double> x, vector<double> y, int newX)
{
  int n = x.size();
  double Sx = 0, Sy = 0, Sxx = 0, Sxy = 0;
  for (int i = 0; i < n; i++)
  {
    Sx += log(x[i]);
    Sy += log(y[i]);
    Sxx += log(x[i]) * log(x[i]);
    Sxy += log(x[i]) * log(y[i]);
  }
  double b = (n * Sxy - Sx * Sy) / (n * Sxx - Sx * Sx);

  double A = (Sy - b * Sx) / n;
  double a = exp(A);

  out << "The Transcendental Equation is: " << endl;
  out << "y = " << a << "x^" << b << endl;

  double newY = a * pow(newX, b);
  out << "For x = " << newX << ", the predicted value of y is: ";
  out << newY << endl;
}

int main(void)
{
  ifstream in("TranscendentalEquationInput.txt");
  ofstream out("TranscendentalEquationOutput.txt");

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
  Transcendental(out, x, y, newX);

  in.close();
  out.close();
  return 0;
}
```

#### Transcendental Equation Input

```
5
1 50
2 80
3 96
4 120
5 145
6
```

##### Input Format

```
The input is taken from a file named TranscendentalEquationInput.txt.

The first line of input contains an integer n - the number of data points.

The next n lines each contain two real numbers x and y - the coordinates of the data points.

The last line of input contains an integer newX - for which the predicted value of y is to be calculated.
```

#### Transcendental Equation Output

```
Data points entered:
(1.00, 50.00)
(2.00, 80.00)
(3.00, 96.00)
(4.00, 120.00)
(5.00, 145.00)
The Transcendental Equation is:
y = 49.88x^0.64
For x = 6, the predicted value of y is: 157.62
```

##### Output Format

```
The output is written to a file named TranscendentalEquationOutput.txt.

All the entered data points are printed.

Then the best fit transcendental equation is printed in the form
y=ax^b.

Finally, the predicted value of y corresponding to x = newX is printed.
```

---

### Interpolation and Approximation

### Newton's Forward Interpolation

Implemented by 2207008

#### Newton's Forward Interpolation Theory

Newton’s Forward Interpolation is used to approximate the value of a function at a point X using a set of equally spaced data points. This method is used when the interpolation point is near the begining of the given data points. The method constructs a forward difference table from the given data points. It is calculated by repeatedly subtracting each element from the one below it in the previous column, starting with the y-values in the first column. Then the method uses the table and computes the interpolated value with the formula : <br>
y = y<sub>0</sub> + u Δy<sub>0</sub> + (u(u-1)/2!) Δ<sup>2</sup>y<sub>0</sub> + (u(u-1)(u-2)/3!) Δ<sup>3</sup>y<sub>0</sub> + ... + (u(u-1)(u-2)...(u-n+1)/n!) Δ<sup>n</sup>y<sub>0</sub><br>
Here u = (X - x<sub>0</sub>)/h and h = x<sub>1</sub> - x<sub>0</sub>

Algorithm:

1. Given the number of data data point and the data points (x, y).
2. Check if the x-values are equally spaced. If not, then the interpolation is not possible with this method.
3. Read the data point X to interpolate on
4. Compute the forward difference table.
5. Compute u and h.
6. Use the formula above to calculate the interpolated value. <br>

#### Newton's Forward Interpolation Code

```cpp
// Implemented by 2207008

#include <bits/stdc++.h>
using namespace std;

void printDifferenceTable(const vector<vector<double>> &d, int n, ofstream &out)
{
    out << "\nForward Difference Table:" << endl;
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

```

#### Newton's Forward Interpolation Input

```
4
3 180
5 150
7 120
9 90
4
```

##### Input Format

```
The input is read from a file named ForwardInterpolation_input.txt.

The first line contains an integer n (Number of data points).

The next n lines each contain two real numbers: xi yi (the data points), where xi are equally spaced.

The last line contains a real number X, the value at which interpolation is needed.

```

#### Newton's Forward Interpolation Output

```
Newton's Forward Interpolation Method :

Forward Difference Table:
       180        -30          0          0
       150        -30          0
       120        -30
        90

Interpolated answer : 165.000

```

##### Output Format

```
The output is written to a file named ForwardInterpolation_output.txt.

The name of the method is printed first.

Then the forward difference table is printed.

Finally, the interpolated value at X is printed with 3 decimal places.

```

### Newton's Backward Interpolation

Implemented by 2207008

#### Newton's Backward Interpolation Theory

Newton’s Backward Interpolation is used to approximate the value of a function at a point X using a set of equally spaced data points. This method is used when the interpolation point is near the end of the given data points. The method constructs a backward difference table from the given data points. It is calculated by repeatedly subtracting each element in the previous column from the element above it, starting with the y-values in the first column. Then the method uses the table and computes the interpolated value with the formula:<br>
y = y<sub>n</sub> + v ∇y<sub>n</sub> + (v(v+1)/2!) ∇<sup>2</sup>y<sub>n</sub> + (v(v+1)(v+2)/3!) ∇<sup>3</sup>y<sub>n</sub> + ... + (v(v+1)(v+2)...(v+n-1)/n!) ∇<sup>n</sup>y<sub>n</sub><br>
Here v = (X - x<sub>n</sub>)/h and h = x<sub>1</sub> - x<sub>0</sub>

Algorithm:

1. Given the number of data points and the data points (x, y).
2. Check if the x-values are equally spaced. If not, then interpolation is not possible with this method.
3. Read the data point X to interpolate on.
4. Compute the backward difference table.
5. Compute u and h.
6. Use the formula above to calculate the interpolated value.

#### Newton's Backward Interpolation Code

```cpp
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
```

#### Newton's Backward Interpolation Input

```
5
24 28.06
28 30.19
32 32.75
36 34.94
40 40
33
```

##### Input Format

```
The input is read from a file named BackwardInterpolation_input.txt.

The first line contains an integer n (Number of data points).

The next n lines each contain two real numbers: xi yi (the data points), where xi are equally spaced.

The last line contains a real number X, the value at which interpolation is needed.

```

#### Newton's Backward Interpolation Output

```
Newton's Backward Interpolation Method :

Difference Table:
     28.06
     30.19       2.13
     32.75       2.56       0.43
     34.94       2.19      -0.37       -0.8
        40       5.06       2.87       3.24       4.04

Interpolated answer = 33.275

```

##### Output Format

```
The output is written to a file named BackwardInterpolation_output.txt.

The name of the method is printed first.

Then the backward difference table is printed.

Finally, the interpolated value at X is printed with 3 decimal places.

```

### Divided Difference Interpolation

#### Divided Difference Interpolation Theory

Newton’s Divided Difference Interpolation is used to approximate the value of a function at a point X using a set of data points, which may be unequally spaced. The method constructs a divided difference table from the given data points. It is computed in the following method.<br>
f[x<sub>i</sub>, x<sub>j</sub>] = (f(x<sub>i</sub>) - f(x<sub>j</sub>)) / (x<sub>i</sub> - x<sub>j</sub>)<br>
f[x<sub>i</sub>, x<sub>j</sub>, x<sub>k</sub>] = (f[x<sub>i</sub>, x<sub>j</sub>] - f[x<sub>j</sub>, x<sub>k</sub>]) / (x<sub>i</sub> - x<sub>k</sub>)<br>
f[x<sub>i</sub>, x<sub>j</sub>, x<sub>k</sub>, x<sub>l</sub>] = (f[x<sub>i</sub>, x<sub>j</sub>, x<sub>k</sub>] - f[x<sub>j</sub>, x<sub>k</sub>, x<sub>l</sub>]) / (x<sub>i</sub> - x<sub>l</sub>)<br>

Then the interpolated value is calculated using the formula:<br>
f(x<sub>n</sub>) = f(x<sub>0</sub>) + (x - x<sub>0</sub>) f[x<sub>1</sub>, x<sub>0</sub>] + (x - x<sub>0</sub>)(x - x<sub>1</sub>) f[x<sub>2</sub>, x<sub>1</sub>, x<sub>0</sub>] + (x - x<sub>0</sub>)(x - x<sub>1</sub>)(x - x<sub>2</sub>) f[x<sub>3</sub>, x<sub>2</sub>, x<sub>1</sub>, x<sub>0</sub>] + ... + (x - x<sub>0</sub>)(x - x<sub>1</sub>)...(x - x<sub>n-1</sub>) f[x<sub>n</sub>, x<sub>n-1</sub>, ..., x<sub>0</sub>]<br>

Algorithm:

1. Read the number of data points and the data points (x, y).
2. Read the value X to interpolate.
3. Initialize the divided difference table with y-values in the first column.
4. Compute higher-order divided differences using the formula above.
5. Use the divided difference formula to compute the interpolated value at X.

#### Divided Difference Interpolation Code

```cpp
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
    out << fixed << setprecision(6);
    out << "Newton's Divided Difference Interpolation Method : " << endl;
    out << "Interpolated value: " << ans << endl;
    in.close();
    out.close();
    return 0;
}

```

#### Divided Difference Interpolation Input

```
4
1 0
4 1.386294
6 1.79175
5 1.609438
2
```

##### Input Format

```
The input is read from a file named DividedDifference_input.txt.

The first line contains an integer n (Number of data points).

The next n lines each contain two real numbers: xi yi (the data points), which may be unequally spaced.

The last line contains a real number X, the value at which interpolation is needed.

```

#### Divided Difference Interpolation Output

```
Newton's Divided Difference Interpolation Method :
Interpolated value: 0.628762

```

##### Output Format

```
The output is written to a file named DividedDifference_output.txt.

The name of the method is printed first.

Finally, the interpolated value at X is printed with 6 decimal places.

```

### Divided Difference Interpolation with Error

Implemented by 2207008

#### Divided Difference Interpolation with Error Theory

Newton’s Divided Difference Interpolation with error calculation is used to approximate the error. If an additional data point f(x<sub>n+1</sub>) is given then the truncation error can be approximated as : <br>
R<sub>n</sub> ≅ f[x<sub>n+1</sub>, x<sub>n</sub>, x<sub>n-1</sub>, ..., x<sub>1</sub>, x<sub>0</sub>] (x - x<sub>0</sub>)(x - x<sub>1</sub>) ... (x - x<sub>n</sub>)<br>

Algorithm:

1. Read the number of data points n and the first n data points (x, y).
2. Read the extra point for error estimation.
3. Read the value X to interpolate.
4. Initialize the divided difference table with y-values in the first column.
5. Compute higher-order divided differences as before.
6. Compute the interpolated value using the Newton divided difference formula.
7. Compute the estimated error using the extra point with the above formula given above.

#### Divided Difference Interpolation with Error Code

```cpp
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

```

#### Divided Difference Interpolation with Error Input

```
3
1 0
4 1.386294
6 1.79175
5 1.609438
2
```

##### Input Format

```
The input is read from a file named DividedDifferenceWithError_input.txt.

The first line contains an integer n (Number of data points excluding extra point for error estimation).

The next n lines each contain two real numbers: xi yi (the data points).

The (n+1)-th line contains an extra point (x, y) used for error estimation.

The last line contains a real number X, the value at which interpolation is needed.

```

#### Divided Difference Interpolation with Error Output

```
Newton's Divided Difference Interpolation Method for Calculating Error:
Interpolated value : 0.565846
Estimated error using the extra point: 0.062916

```

##### Output Format

```
The output is written to a file named DividedDifferenceWithError_output.txt.

The name of the method is printed first.

Then the interpolated value at X is printed with 6 decimal places.

Finally, the estimated error using the extra point is printed.

```
