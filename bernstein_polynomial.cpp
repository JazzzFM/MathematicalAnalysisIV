# include <cmath>
# include <cstdlib>
# include <cstring>
# include <ctime>
# include <iomanip>
# include <iostream>
# include "bernstein_polynomial.hpp"
# include <string>

using namespace std;


double *bernstein_matrix(int n){
//  Objetivo:
//
// BERNSTEIN_MATRIX devuelve la matriz de Bernstein.
//
//  Discusión:
//
// La matriz de Bernstein de orden N es una matriz A de NxN que se puede utilizar para
// transformar un vector de coeficientes de base de potencia C que representan un polinomio
// P(X) a un vector de coeficiente de base de Bernstein correspondiente B:
//
// B = A * C
//
// Los N vectores de base de potencia se ordenan como (1,X,X^2,...X^(N-1)) y los N
// Vectores de base de Bernstein como ((1-X)^(N-1), X*(1_X)^(N-2),...,X^(N-1)).
//
// Para N = 5, la matriz tiene la forma:
//
// 1 -4   6  -4  1
// 0  4 -12  12 -4
// 0  0   6 -12  6
// 0  0   0   4 -4
// 0  0   0   0  1
//
// y los números en cada columna representan los coeficientes en la potencia
// desarrollo en serie de un polinomio de Bernstein, de modo que
//
// B(5,4) = - 4x^4 + 12x^3 - 12x^2 + 4x
//
// Parámetros:
//	Entrada, int N, el orden de la matriz.
//
// Salida:
// 	doble BERNSTEIN_MATRIX[N*N], la matriz de Bernstein.
//

  double *a;
  int i;
  int j;

  a = new double[n*n];
 
  for(j = 0; j < n; j++){
      
      for(i = 0; i <= j; i++){
          a[i+j*n] = r8_mop(j-i)*r8_choose(n-1-i, j-i)*r8_choose(n-1,i);
   }
      for(i = j+1; i < n; i++){
         
	 a[i+j*n] = 0.0;
    }
  }
  
  return a;
}

double bernstein_matrix_determinant(int n){
//  Objetivo:
//	BERNSTEIN_MATRIX_DETERMINANT devuelve el determinante de la matriz de BERNSTEIN.
//
// Parámetros:
// 	Entrada, int N, el orden de la matriz.
//
// Salida: 
// 	doble BERNSTEIN_MATRIX_DETERMINANT, el determinante.

  int i;
  double value;

  value = 1.0;
  for(i = 0; i < n; i++){
	value = value*r8_choose(n-1, i);
  }

  return value;
}

double *bernstein_matrix_inverse(int n){
//
//  Objetivo:
//  	BERNSTEIN_MATRIX_INVERSE devuelve la matriz de Bernstein inversa.
//
//  Discusión:
//
// La matriz inversa de Bernstein de orden N es una matriz A NxN que puede
// ser usado para transformar un vector de coeficientes de base de Bernstein B
// representando un polinomio P(X) a una base de potencia correspondiente
// coeficiente vector C:
//
// C = A * B
//
// Los N vectores de base de potencia se ordenan como (1,X,X^2,...X^(N-1)) y los N
// Vectores de base de Bernstein como ((1-X)^(N-1), X*(1_X)^(N-2),...,X^(N-1)).
//
// Para N = 5, la matriz tiene la forma:
//
// 1  1   1   1   1
// 0  1/4 1/2 3/4 1
// 0  0   1/6 1/2 1
// 0  0   0   1/4 1
// 0  0   0   0   1
//
//
// Parámetros:
// 	Entrada, int N, el orden de la matriz.
//
// Salida:
// 	doble BERNSTEIN_MATRIX_INVERSE[N*N], la matriz inversa de Bernstein.
//
  
  double *a;
  int i;
  int j;

  a = new double[n*n];

  for(j = 0; j < n; j++){
    
	  for(i = 0; i <= j; i++){
              
              a[i+j*n] = r8_choose(j, i)/r8_choose(n-1, i);
           }
    
           for(i = j + 1; i < n; i++){
              a[i+j*n] = 0.0;
           }
  }

  return a;
}

double *bernstein_poly_01(int n, double x){
//  Objetivo:
//	BERNSTEIN_POLY_01 evalúa los polinomios de Bernstein basados en [0,1].
//
//  Discusión:
// 	 Se supone que los polinomios de Bernstein se basan en [0,1].
//
// La fórmula es:
// 	B(N,I)(X) = [N!/(I!*(N-I)!)] * (1-X)^(N-I) * X^I
//
// Primeros valores:
//
//    B(0,0)(X) = 1
//
//    B(1,0)(X) =      1-X
//    B(1,1)(X) =                X
//
//    B(2,0)(X) =     (1-X)^2
//    B(2,1)(X) = 2 * (1-X)    * X
//    B(2,2)(X) =                X^2
//
//    B(3,0)(X) =     (1-X)^3
//    B(3,1)(X) = 3 * (1-X)^2 * X
//    B(3,2)(X) = 3 * (1-X)   * X^2
//    B(3,3)(X) =               X^3
//
//    B(4,0)(X) =     (1-X)^4
//    B(4,1)(X) = 4 * (1-X)^3 * X
//    B(4,2)(X) = 6 * (1-X)^2 * X^2
//    B(4,3)(X) = 4 * (1-X)   * X^3
//    B(4,4)(X) =               X^4
//
// Valores especiales:
//
// B(N,I)(X) tiene un valor máximo único en X = I/N.
//
// B(N,I)(X) tiene un I-fold cero en 0 y N-I fold zero en 1.
//
// B(N,I)(1/2) = C(N,K) / 2^N
//
// Para X y N fijos, los polinomios suman 1:
//
// Suma ( 0 <= I <= N ) B(N,I)(X) = 1
//
// Parámetros:
//
// Entrada:
// 	- int N, el grado de los polinomios de Bernstein
//    	para ser utilizado. Para cualquier N, hay un conjunto de N+1 polinomios de Bernstein,
// 	cada uno de grado N, que forman una base para polinomios en [0,1].
//
// 	- doble X, el punto de evaluación.
//
// Salida:
// 	- double BERNSTEIN_POLY[N+1], los valores de N+1 Polinomios de Bernstein en X.


  double *bern;
  int i;
  int j;

  bern = new double[n+1];

  if(n == 0){
    
	  bern[0] = 1.0;
  }else if(0 < n){

    bern[0] = 1.0 - x;
    bern[1] = x;
 
    for(i = 2; i <= n; i++){

      bern[i] = x * bern[i-1];
      
      for(j = i - 1; 1 <= j; j--){
        bern[j] = x*bern[j-1] + (1.0-x)*bern[j];
      
      }
      
      bern[0] = (1.0 - x)*bern[0];
    }
  }
  return bern;
}

double *bernstein_poly_01_matrix(int m, int n, double x[]){
//  Objetivo:
//	- BERNSTEIN_POLY_01_MATRIX evalúa los polinomios de Bernstein basados en [0,1].
//
//  Discusión:
//	Se supone que los polinomios de Bernstein se basan en [0,1].
//	La fórmula es:
//		B(N,I)(X) = [N!/(I!*(N-I)!)] * (1-X)^(N-I) * X^I
//
//    Primeros valores:
//
//    B(0,0)(X) = 1
//
//    B(1,0)(X) =      1-X
//    B(1,1)(X) =                X
//
//    B(2,0)(X) =     (1-X)^2
//    B(2,1)(X) = 2 * (1-X)    * X
//    B(2,2)(X) =                X^2
//
//    B(3,0)(X) =     (1-X)^3
//    B(3,1)(X) = 3 * (1-X)^2 * X
//    B(3,2)(X) = 3 * (1-X)   * X^2
//    B(3,3)(X) =               X^3
//
//    B(4,0)(X) =     (1-X)^4
//    B(4,1)(X) = 4 * (1-X)^3 * X
//    B(4,2)(X) = 6 * (1-X)^2 * X^2
//    B(4,3)(X) = 4 * (1-X)   * X^3
//    B(4,4)(X) =               X^4
//
// Valores especiales:
//
// B(N,I)(X) tiene un valor máximo único en X = I/N.
//
// B(N,I)(1/2) = C(N,K) / 2^N
//
// Para X y N fijos, los polinomios suman 1:
//
// Suma( 0 <= I <= N ) B(N,I)(X) = 1
//
// Parámetros:
//
// 	- Entrada: int M, el número de puntos de evaluación.
//	- Entrada: int N, el grado de los polinomios de Bernstein
//    		para ser utilizado. Para cualquier N, hay un conjunto de N+1 polinomios de Bernstein,
// 		cada uno de grado N, que forman una base para polinomios en [0,1].
// 	- Entrada, doble X[M], los puntos de evaluación.
//
// Salida:
// 	- doble BERNSTEIN_POLY_01_MATRIX[M*(N+1)], los valores de N+1 Polinomios de Bernstein en los puntos de evaluación.

  double *b;
  int i;
  int j;
  int k;

  b = new double[m*(n+1)];

  for(i = 0; i < m; i++){
      
      if(n == 0){
          b[i+0*m] = 1.0;
    }

    else if(0 < n){
      
      b[i+0*m] = 1.0 - x[i];
      b[i+1*m] = x[i];
 
      for(j = 2; j <= n; j++){
        
          b[i+j*m] = x[i] * b[i+(j-1)*m];
        
	  for(k = j - 1; 1 <= k; k--){
          
              b[i+k*m] = x[i]*b[i+(k-1)*m] + (1.0 - x[i])*b[i+k*m];
          }

        b[i+0*m] = (1.0 - x[i])*b[i+0*m];
      
      }
    }
  }
  return b;
}


void bernstein_poly_01_values(int &n_data, int &n, int &k, double &x, double &b){
//  Objetivo:
//	BERNSTEIN_POLY_01_VALUES devuelve algunos valores de los polinomios de Bernstein.
//
//  Discusión:
//	Se supone que los polinomios de Bernstein se basan en [0,1].
//	La fórmula para los polinomios de Bernstein es
//		B(N,I)(X) = [N!/(I!(N-I)!)] * (1-X)^(N-I) * X^I
//
// En Mathematica, la función puede ser evaluada por: Binomial[n,i] * (1-x)^(n-i) * x^i
//
// Primeros valores:
//
//    B(0,0)(X) = 1
//
//    B(1,0)(X) =      1-X
//    B(1,1)(X) =                X
//
//    B(2,0)(X) =     (1-X)^2
//    B(2,1)(X) = 2 * (1-X)    * X
//    B(2,2)(X) =                X^2
//
//    B(3,0)(X) =     (1-X)^3
//    B(3,1)(X) = 3 * (1-X)^2  * X
//    B(3,2)(X) = 3 * (1-X)    * X^2
//    B(3,3)(X) =                X^3
//
//    B(4,0)(X) =     (1-X)^4
//    B(4,1)(X) = 4 * (1-X)^3  * X
//    B(4,2)(X) = 6 * (1-X)^2  * X^2
//    B(4,3)(X) = 4 * (1-X)    * X^3
//    B(4,4)(X) =                X^4
// 
// Valores especiales:
//
// B(N,I)(X) tiene un valor máximo único en X = I/N.
//
// B(N,I)(1/2) = C(N,K) / 2^N
//
// Para X y N fijos, los polinomios suman 1:
//
// Sum( 0 <= I <= N ) B(N,I)(X) = 1
//
// Parámetros:
//
// Entrada/salida, int &N_DATA. El usuario establece N_DATA en 0 antes de que
//    	primera llamada. En cada llamada, la rutina incrementa N_DATA en 1 y
// 	devuelve los datos correspondientes; cuando no hay más datos, el
// 	el valor de salida de N_DATA volverá a ser 0.
//
// Salida, int &N, el grado del polinomio.
//
// Salida, int &K, el índice del polinomio.
//
// Salida, doble &X, el argumento del polinomio.
//
// Salida, double &B, el valor del polinomio B(N,K)(X).
//

# define N_MAX 15

  static double b_vec[N_MAX] = {
     0.1000000000000000E+01,
     0.7500000000000000E+00,
     0.2500000000000000E+00,
     0.5625000000000000E+00,
     0.3750000000000000E+00,
     0.6250000000000000E-01,
     0.4218750000000000E+00,
     0.4218750000000000E+00,
     0.1406250000000000E+00,
     0.1562500000000000E-01,
     0.3164062500000000E+00,
     0.4218750000000000E+00,
     0.2109375000000000E+00,
     0.4687500000000000E-01,
     0.3906250000000000E-02 };

  static int k_vec[N_MAX] = {
    0,
    0, 1,
    0, 1, 2,
    0, 1, 2, 3,
    0, 1, 2, 3, 4 };

  static int n_vec[N_MAX] = {
    0,
    1, 1,
    2, 2, 2,
    3, 3, 3, 3,
    4, 4, 4, 4, 4 };

  static double x_vec[N_MAX] = {
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00,
     0.25E+00 };

  if (n_data < 0){
      n_data = 0;
  }

  n_data = n_data + 1;

  if(N_MAX < n_data){

    n_data = 0;
    n = 0;
    k = 0;
    x = 0.0;
    b = 0.0;

  }
  else{

    n = n_vec[n_data-1];
    k = k_vec[n_data-1];
    x = x_vec[n_data-1];
    b = b_vec[n_data-1];
  }

  return;
# undef N_MAX
}

double *bernstein_poly_ab(int n, double a, double b, double x){
//  Objetivo:
//
// BERNSTEIN_POLY_AB evalúa en X los polinomios de Bernstein basados en [A,B].
//
//  Discusión:
//
// La fórmula es:
//
// BERNA(N,I)(X) = [N!/(I!*(N-I)!)] * (B-X)^(N-I) * (X-A)^I / (B-A)^N
//
// Primeros valores
//
//    B(0,0)(X) =   1
//
//    B(1,0)(X) = (      B-X                ) / (B-A)
//    B(1,1)(X) = (                 X-A     ) / (B-A)
//
//    B(2,0)(X) = (     (B-X)^2             ) / (B-A)^2
//    B(2,1)(X) = ( 2 * (B-X)    * (X-A)    ) / (B-A)^2
//    B(2,2)(X) = (                (X-A)^2  ) / (B-A)^2
//
//    B(3,0)(X) = (     (B-X)^3             ) / (B-A)^3
//    B(3,1)(X) = ( 3 * (B-X)^2  * (X-A)    ) / (B-A)^3
//    B(3,2)(X) = ( 3 * (B-X)    * (X-A)^2  ) / (B-A)^3
//    B(3,3)(X) = (                (X-A)^3  ) / (B-A)^3
//
//    B(4,0)(X) = (     (B-X)^4             ) / (B-A)^4
//    B(4,1)(X) = ( 4 * (B-X)^3  * (X-A)    ) / (B-A)^4
//    B(4,2)(X) = ( 6 * (B-X)^2  * (X-A)^2  ) / (B-A)^4
//    B(4,3)(X) = ( 4 * (B-X)    * (X-A)^3  ) / (B-A)^4
//    B(4,4)(X) = (                (X-A)^4  ) / (B-A)^4
// Parámetros:
//
// 	- Entrada, int N, el grado de los polinomios de Bernstein
//    	para ser utilizado. Para cualquier N, hay un conjunto de N+1 polinomios de Bernstein,
// 	cada uno de grado N, que forman una base para polinomios en [A,B].
//
// 	- Ingrese, double A, B, los puntos finales del intervalo en el que
// 	los polinomios se van a basar. A y B no deben ser iguales.
//
// 	- Entrada, doble X, el punto en el que los polinomios van a ser evaluados.
//
// 	- Salida, doble BERNSTEIN_POLY_AB[N+1], los valores de N+1 Polinomios de Bernstein en X.

  double *bern;
  int i;
  int j;

  if(b == a){
    
    cerr << "\n";
    cerr << "BERNSTEIN_POLYNOMIAL_AB - Fatal error!\n";
    cerr << "  A = B = " << a << "\n";
    exit ( 1 );
  }

  bern = new double[n+1];

  if(n == 0){
      bern[0] = 1.0;

   }else if(0 < n){
    
    bern[0] = (b - x)/(b - a);
    bern[1] = (x - a)/(b - a);
 
    for(i = 2; i <= n; i++){
      
      bern[i] = (x - a)*bern[i-1]/(b - a);
      
      for(j = i - 1; 1 <= j; j--){
      
          bern[j] = ((b - x)*bern[j] + (x - a)*bern[j-1])/(b- a);
      
	}
      
      bern[0] = ( b - x ) * bern[0] / ( b - a );
    }
  }

  return bern;
}


double *bernstein_poly_ab_approx(int n, double a, double b, double ydata[], int nval, double xval[]){
//  Objetivo:
//	- BERNSTEIN_POLY_AB_APPROX: Aproximación de Bernstein a F(X) en [A,B].
//
// Fórmula:
//	 BPAB(F)(X) = sum(0 <= I <= N ) F(X(I)) * B_BASE(I,X)
//
//    dónde
//	- X(I) = ( ( N - I ) * A + I * B ) / N
// 	- B_BASE(I,X) es el valor del I-ésimo polinomio de base de Bernstein en X.
//
//  Discusión:
//	El polinomio de Bernstein BPAB(F) para F(X) sobre [A,B] es una aproximación,
// 	no es un interpolador; en otras palabras, no se garantiza que su valor sea igual
// 	el de F en cualquier punto en particular. Sin embargo, para un intervalo fijo
// 	[A,B], si dejamos que N aumente, el polinomio de Bernstein converge
// 	uniformemente a F en todas partes en [A,B], siempre que F sea continua.
// 	Incluso si F no es continua, pero está acotada, el polinomio converge
// 	puntualmente a F(X) en todos los puntos de continuidad. Por otro lado,
// 	la convergencia es bastante lenta comparada con otras interpolaciones
//	 y esquemas de aproximación.
//  Parámetros:
//
// 	Entrada, int N, el grado del polinomio de Bernstein
//    	para ser utilizado. N debe ser al menos 0.
//
// 	Ingrese, double A, B, los puntos finales del intervalo en el que
// 	se basa el aproximado. A y B no deben ser iguales.
//
// 	Entrada, doble YDATA[N+1], los valores de datos en N+1 por igual
// 	puntos espaciados en [A,B]. Si N = 0, entonces el punto de evaluación debería
// 	ser 0.5 * (A + B). De lo contrario, el punto de evaluación debería ser
// 	( (N-I)*A + I*B ) / N ).
//
// 	Ingrese, int NVAL, el número de puntos en los que el
// 	se va a evaluar el aproximado.
//
// 	Entrada, double XVAL[NVAL], el punto en el que Bernstein
// 	se va a evaluar la aproximante del polinomio. Las entradas de XVAL no
// 	tiene que estar en el intervalo [A,B].
//
// 	Salida, doble BPAB_APPROX[NVAL], los valores de Bernstein
// 	polinomio aproximado para F, basado en [A,B], evaluado en XVAL.

  double *bvec;
  int i;
  double *yval;

  yval = new double[nval];

  for(i = 0; i < nval; i++){
    //
    //  Evaluate the Bernstein basis polynomials at XVAL.
     bvec = bernstein_poly_ab ( n, a, b, xval[i] );
   
    //
    //  Now compute the sum of YDATA(I) * BVEC(I).

    yval[i] = r8vec_dot_product(n + 1, ydata, bvec);

    delete [] bvec;
  }

  return yval;
}

int i4_max(int i1, int i2){
//
//  Purpose:
//
//    I4_MAX returns the maximum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, are two integers to be compared.
//
//    Output, int I4_MAX, the larger of I1 and I2.
//
  
  int value;

  if(i2 < i1){
    
      value = i1;
  }
  else{
    
      value = i2;
  }

  return value;
}

int i4_min(int i1, int i2){
//
//  Purpose:
//
//    I4_MIN returns the minimum of two I4's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    13 October 1998
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I1, I2, two integers to be compared.
//
//    Output, int I4_MIN, the smaller of I1 and I2.
//

  int value;

  if(i1 < i2){
    
      value = i1;
  }
  else{
      value = i2;
  }

  return value;
}


double r8_choose(int n, int k){
//
//  Purpose:
//
//    R8_CHOOSE computes the binomial coefficient C(N,K) as an R8.
//
//  Discussion:
//
//    The value is calculated in such a way as to avoid overflow and
//    roundoff.  The calculation is done in R8 arithmetic.
//
//    The formula used is:
//
//      C(N,K) = N! / ( K! * (N-K)! )
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    ML Wolfson, HV Wright,
//    Algorithm 160:
//    Combinatorial of M Things Taken N at a Time,
//    Communications of the ACM,
//    Volume 6, Number 4, April 1963, page 161.
//
//  Parameters:
//
//    Input, int N, K, the values of N and K.
//
//    Output, double R8_CHOOSE, the number of combinations of N
//    things taken K at a time.
//
  int i;
  int mn;
  int mx;
  double value;

  mn = i4_min ( k, n - k );

  if(mn < 0){
  
      value = 0.0;
  }
  else if(mn == 0){
      value = 1.0;
  }
  else{
    
    mx = i4_max(k, n - k);
    value = (double)(mx + 1);

    for(i = 2; i <= mn; i++){

      value = (value*(double)(mx + i))/(double)i;
    }
  }
  return value;
}

double r8_max(double x, double y){
//
//  Purpose:
//
//    R8_MAX returns the maximum of two R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, Y, the quantities to compare.
//
//    Output, double R8_MAX, the maximum of X and Y.
//

  double value;

  if(y < x){

    value = x;
  
  }else{

    value = y;
  }

  return value;
}

double r8_mop(int i){
//
//  Purpose:
//
//    R8_MOP returns the I-th power of -1 as an R8 value.
//
//  Discussion:
//
//    An R8 is an double value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    16 November 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int I, the power of -1.
//
//    Output, double R8_MOP, the I-th power of -1.
//
  
  double value;

  if ((i % 2) == 0){
    
      value = 1.0;
  }
  else{
    
      value = -1.0;
  }

  return value;
}

double r8_uniform_01(int &seed){
//
//  Purpose:
//
//    R8_UNIFORM_01 returns a unit pseudorandom R8.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = ( 16807 * seed ) mod ( 2^31 - 1 )
//      u = seed / ( 2^31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      R8_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 April 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Second Edition,
//    Springer, 1987,
//    ISBN: 0387964673,
//    LC: QA76.9.C65.B73.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, December 1986, pages 362-376.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation,
//    edited by Jerry Banks,
//    Wiley, 1998,
//    ISBN: 0471134031,
//    LC: T57.62.H37.
//
//    Peter Lewis, Allen Goodman, James Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, Number 2, 1969, pages 136-143.
//
//  Parameters:
//
//    Input/output, int &SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double R8_UNIFORM_01, a new pseudorandom variate, 
//    strictly between 0 and 1.
//

  const int i4_huge = 2147483647;
  int k;
  double r;

  if (seed == 0){
    cerr << "\n";
    cerr << "R8_UNIFORM_01 - Fatal error!\n";
    cerr << "  Input value of SEED = 0.\n";
    exit ( 1 );
  }

  k = seed/127773;

  seed = 16807*(seed-k*127773)-k*2836;

  if (seed < 0){

    seed = seed + i4_huge;
  }

  r = (double)(seed)*4.656612875E-10;

  return r;
}

double r8mat_is_identity(int n, double a[]){
//
//  Purpose:
//
//    R8MAT_IS_IDENTITY determines if an R8MAT is the identity.
//
//  Discussion:
//
//    An R8MAT is a matrix of real ( kind = 8 ) values.
//
//    The routine returns the Frobenius norm of A - I.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 July 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the matrix.
//
//    Output, double R8MAT_IS_IDENTITY, the Frobenius norm
//    of the difference matrix A - I, which would be exactly zero
//    if A were the identity matrix.
//

  double error_frobenius;
  int i;
  int j;
  double t;

  error_frobenius = 0.0;

  for(i = 0; i < n; i++){
    for(j = 0; j < n; j++){

      if(i == j){

        t = a[i+j*n] - 1.0;
      
      }
      else{
        
	t = a[i+j*n];
      }
      error_frobenius = error_frobenius + t*t;
    }
  }

  error_frobenius = sqrt ( error_frobenius );

  return error_frobenius;
}

double *r8mat_mm_new(int n1, int n2, int n3, double a[], double b[]){
//
//  Purpose:
//
//    R8MAT_MM_NEW multiplies two matrices.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as the function value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    18 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N1, N2, N3, the order of the matrices.
//
//    Input, double A[N1*N2], double B[N2*N3], the matrices to multiply.
//
//    Output, double R8MAT_MM[N1*N3], the product matrix C = A * B.
//

  double *c;
  int i;
  int j;
  int k;

  c = new double[n1*n3];

  for(i = 0; i < n1; i ++){
    for(j = 0; j < n3; j++){

      c[i+j*n1] = 0.0;

      for(k = 0; k < n2; k++){
          
	   c[i+j*n1] = c[i+j*n1] + a[i+k*n1] * b[k+j*n2];
      }
    }
  }

  return c;
}

double *r8mat_mv_new(int m, int n, double a[], double x[]){
//
//  Purpose:
//
//    R8MAT_MV_NEW multiplies a matrix times a vector.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    For this routine, the result is returned as the function value.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    11 April 2007
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns of the matrix.
//
//    Input, double A[M,N], the M by N matrix.
//
//    Input, double X[N], the vector to be multiplied by A.
//
//    Output, double R8MAT_MV_NEW[M], the product A*X.
//

  int i;
  int j;
  double *y;

  y = new double[m];

  for(i = 0; i < m; i++){

    y[i] = 0.0;

    for(j = 0; j < n; j++){

      y[i] = y[i] + a[i+j*m] * x[j];
    }
  }

  return y;
}

double r8mat_norm_fro(int m, int n, double a[]){
//
//  Purpose:
//
//    R8MAT_NORM_FRO returns the Frobenius norm of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    The Frobenius norm is defined as
//
//      R8MAT_NORM_FRO = sqrt (
//        sum ( 1 <= I <= M ) sum ( 1 <= j <= N ) A(I,J)**2 )
//    The matrix Frobenius norm is not derived from a vector norm, but
//    is compatible with the vector L2 norm, so that:
//
//      r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, double A[M*N], the matrix whose Frobenius
//    norm is desired.
//
//    Output, double R8MAT_NORM_FRO, the Frobenius norm of A.
//

  int i;
  int j;
  double value;

  value = 0.0;
  for(j = 0; j < n; j++){
      for( i = 0; i < m; i++){

      value = value + pow(a[i+j*m], 2 );
    }
  }

  value = sqrt(value);

  return value;
}

void r8mat_print(int m, int n, double a[], string title){
//
//  Purpose:
//
//    R8MAT_PRINT prints an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//    Entry A(I,J) is stored as A[I+J*M]
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    10 September 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows in A.
//
//    Input, int N, the number of columns in A.
//
//    Input, double A[M*N], the M by N matrix.
//
//    Input, string TITLE, a title.
//
  
    r8mat_print_some(m, n, a, 1, 1, m, n, title);

  return;
}

void r8mat_print_some(int m, int n, double a[], int ilo, int jlo, int ihi, int jhi, string title){
//
//  Purpose:
//
//    R8MAT_PRINT_SOME prints some of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    20 August 2010
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, the number of rows of the matrix.
//    M must be positive.
//
//    Input, int N, the number of columns of the matrix.
//    N must be positive.
//
//    Input, double A[M*N], the matrix.
//
//    Input, int ILO, JLO, IHI, JHI, designate the first row and
//    column, and the last row and column to be printed.
//
//    Input, string TITLE, a title.
//

# define INCX 5

  int i;
  int i2hi;
  int i2lo;
  int j;
  int j2hi;
  int j2lo;

  cout << "\n";
  cout << title << "\n";

  if( m <= 0 || n <= 0 ){
    
    cout << "\n";
    cout << "  (None)\n";
    return;
  }

    //  Print the columns of the matrix, in strips of 5.

  for(j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX){
    
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min(j2hi, n);
    j2hi = i4_min(j2hi, jhi);

    cout << "\n";
   //  For each column J in the current range...
   //  Write the header.

    cout << "  Col:    ";
    for(j = j2lo; j <= j2hi; j++){

      cout << setw(7) << j - 1 << "       ";
    }

    cout << "\n";
    cout << "  Row\n";
    cout << "\n";
    
    //  Determine the range of the rows in this strip.

    i2lo = i4_max(ilo, 1);
    i2hi = i4_min(ihi, m);

    for(i = i2lo; i <= i2hi; i++){

    //  Print out (up to) 5 entries in row I, that lie in the current strip.

      cout << setw(5) << i - 1 << ": ";

      for(j = j2lo; j <= j2hi; j++){
        cout << setw(12) << a[i-1+(j-1)*m] << "  ";
      }
      cout << "\n";
    }
  }

  return;
# undef INCX
}

double *r8mat_zeros_new(int m, int n){
//
//  Purpose:
//
//    R8MAT_ZEROS_NEW returns a new zeroed R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 October 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int M, N, the number of rows and columns.
//
//    Output, double R8MAT_ZEROS_NEW[M*N], the new zeroed matrix.
//

  double *a;
  int i;
  int j;

  a = new double[m*n];

  for(j = 0; j < n; j++){
      for(i = 0; i < m; i++){

      a[i+j*m] = 0.0;
    }
  }
  return a;
}

double r8vec_dot_product(int n, double a1[], double a2[]){
//
//  Purpose:
//
//    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vectors.
//
//    Input, double A1[N], A2[N], the two vectors to be considered.
//
//    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
//

  int i;
  double value;

  value = 0.0;
  for(i = 0; i < n; i++){
    
     value = value + a1[i] * a2[i];
  }

  return value;
}

double *r8vec_linspace_new(int n, double a_first, double a_last){
//
//  Purpose:
//
//    R8VEC_LINSPACE_NEW creates a vector of linearly spaced values.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    29 March 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A_FIRST, A_LAST, the first and last entries.
//
//    Output, double R8VEC_LINSPACE_NEW[N], a vector of linearly spaced data.
//

  double *a;
  int i;

  a = new double[n];

  if(n == 1){
  
      a[0] = ( a_first + a_last ) / 2.0;
  }
  else{
      for(i = 0; i < n; i++){

      a[i] = ((double)(n-1-i)*a_first + (double)(i)*a_last)/(double)(n-1);
    }
  }
  return a;
}

double r8vec_sum(int n, double a[]){
//
//  Purpose:
//
//    R8VEC_SUM returns the sum of an R8VEC.
//
//  Discussion:
//
//    An R8VEC is a vector of R8's.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    15 October 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the number of entries in the vector.
//
//    Input, double A[N], the vector.
//
//    Output, double R8VEC_SUM, the sum of the vector.
//

  int i;
  double value;

  value = 0.0;
  
  for(i = 0; i < n; i++){

    value = value + a[i];
  }

  return value;
}

void timestamp(){
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    08 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct std::tm *tm_ptr;
  size_t len;
  std::time_t now;

  now = std::time(NULL);
  tm_ptr = std::localtime(&now);

  len = std::strftime(
		time_buffer, 
		TIME_SIZE, 
		"%d %B %Y %I:%M:%S %p", 
		tm_ptr );

  std::cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}
