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

Gauss Elimination is a method to solve a system of linear equations by converting the matrix into Row Echelon Form using row operations. It solves the linear equations by using back substitution in the Row Echelon Form matrix.

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

---

### Gauss Jordan Elimination Method

#### Gauss Jordan Theory

Gauss Jordan Elimination is a method to solve a system of linear equations by converting the matrix into Reduced Row Echelon Form using row operations. The constant vector of the Reduced Row Echelon Form matrix is the solution of linear equations.

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

---

### LU Decomposition Method

#### LU Decomposition Theory

[Add your theory content here]

#### LU Decomposition Code

```python
# Add your code here
```

#### LU Decomposition Input

```
[Add your input format here]
```

#### LU Decomposition Output

```
[Add your output format here]
```

---

### Matrix Inversion Method

#### Matrix Inversion Theory

Matrix Inversion Method solves the system of linear equations (AX = B) by finding the inverse matrix (A^-1) of the corresponding co-efficient matrix (A) of the system of linear equations. Then it determines the solution by multiplying the inverse matrix with the constant vector (X = A^-1 x B).

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

---

### Solution of Non-Linear Equations

### Bisection Method

#### Bisection Theory

[Add your theory content here]

#### Bisection Code

```python
# Add your code here
```

#### Bisection Input

```
[Add your input format here]
```

#### Bisection Output

```
[Add your output format here]
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

#### Runge-Kutta Theory

[Add your theory content here]

#### Runge-Kutta Code

```python
# Add your code here
```

#### Runge-Kutta Input

```
[Add your input format here]
```

#### Runge-Kutta Output

```
[Add your output format here]
```

---
