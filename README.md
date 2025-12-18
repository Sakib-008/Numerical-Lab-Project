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

```python
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

```python
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

LU factorrization is a matrix decomposition technique used to solve a system of linear equations efficientlly. This method is also called Cholesky method, Doolittle's method, Crout's method. In this method, a square matrix A is decomposed into the product of two triangluar matrices. A = LU. Here L is a lower triangular matrix whose all diagonal elements are 1 and U is an upper triangular matrix. It simplifies the system AX = B by converting it into two simpler systems. AX = B --> LUX = B --> LY = B where (UX = Y). Then solves LY = B by forward substitution and UX = Y. This method is useful when solving system with the same coefficient matrix A but the constant vector B changes. <br>
Calculation steps:<br>

1. First row of U : u11, u12, u13...
2. First column of L : l21, l31...
3. Second row of U : u22, u23...
4. Second column of L : l32...
5. Third row of U : u33... <br>

Solution Classification:<br>

1. After decomposition if U[i][i] is zero and Y[i] is non-zero then system has no solution.
2. If both are zero, the system has infinite solution.
3. Otherwise, the system has a unique solution.

#### LU Decomposition Code

```python
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
    out<<"LU Factorization method : "<<endl;

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
                        out << "\n------------------------------\n" << endl;
                soln = false;
                break;
            }
            else if (fabs(U[i][i]) < 1e-9 && fabs(Y[i]) < 1e-9)
            {
                out << "\nInfinite Solution!" << endl;
                        out << "\n------------------------------\n" << endl;
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
        out << "\n------------------------------\n" << endl;
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

```python
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
Algorithm:<br>

1. Let x1 = a and x2 = b. Define x0 as the midpoint between a and b. x0 = (x1 + x2) / 2.
2. Evaluate f(x0), f(x1), f(x2).
3. If f(x0) = 0, then the root is x0.
4. Select the subinterval where the sign change occurs based on the following conditions:<br>
   If f(x0) _ f(x1) < 0, then root is between x0 and x1. So take [x0, x1] as the new subinterval.<br>
   If f(x0) _ f(x2) < 0, then root is between x0 and x2. So take [x0, x2] as the new subinterval.
5. Repeat the step 3 while fabs(f(x0)) > E. (Solution is close enough to zero).

#### Bisection Code

```python
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

#### False Position Theory

[Add your theory content here]

#### False Position Code

```python
# Add your code here
```

#### False Position Input

```
[Add your input format here]
```

#### False Position Output

```
[Add your output format here]
```

---

### Secant Method

#### Secant Theory

[Add your theory content here]

#### Secant Code

```python
# Add your code here
```

#### Secant Input

```
[Add your input format here]
```

#### Secant Output

```
[Add your output format here]
```

---

### Newton-Raphson Method

#### Newton-Raphson Theory

[Add your theory content here]

#### Newton-Raphson Code

```python
# Add your code here
```

#### Newton-Raphson Input

```
[Add your input format here]
```

#### Newton-Raphson Output

```
[Add your output format here]
```

---

### Solution of Differential Equations

### Runge-Kutta Method

Implemented by 2207008

#### Runge-Kutta Theory

The Runge-Kutta 4th Order (RK4) method is a numerical technique to approximate the solution of a first-order ordinary differential equation (ODE) which is of the the form: dy/dx = f(x, y), y(x0) = y0. It provies a good balance between accuracy and computational efficiency.<br>

Algorithm:

1. Given a initial value of x --> x0 and y --> y0.
2. Given a step size h, the value of y at the next point is computed using four intermediate slopes: <br>
   k1 = h * f(xn, yn)<br>
   k2 = h * f(xn + h/2. yn + k1/2)<br>
   k3 = h * f(xn + h/2, yn + k2/2)<br>
   k4 = h * f(xn + h, y0 + k3)
3. The next value of y is then calculated as:<br>
   yn+1 = yn + (k1 + 2 * k2 + 2 * k3 + k4) / 6

#### Runge-Kutta Code

```python
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
