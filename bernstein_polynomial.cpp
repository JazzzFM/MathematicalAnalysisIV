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
// Parámetros:
//
// Entrada, int I1, I2, son dos enteros a comparar.
//
// Salida, int I4_MAX, el mayor de I1 e I2.

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
//  Objetivo:
//	I4_MIN devuelve el mínimo de dos I4.
//
// Parámetros:
//	 Entrada, int I1, I2, dos enteros a comparar.
//	 Salida, int I4_MIN, el menor de I1 e I2.

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
//  Objetivo:
//	 R8_CHOOSE calcula el coeficiente binomial C(N,K) como un R8.
//
//  Discusión:
//
// 	El valor se calcula de forma que se evite el desbordamiento y
//    	redondear. El cálculo se realiza en aritmética R8.
//
// 	La fórmula utilizada es:
//
// 	C(N,K) = N! / ( K! * (N-K)! )
// Parámetros:
//
// Ingrese, int N, K, los valores de N y K.
//
// Salida:
// 	 -double R8_CHOOSE, el número de combinaciones de N
// 		cosas tomadas K a la vez.
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
//  Objetivo:
//
// R8_MAX devuelve el máximo de dos R8.
// Parámetros:
//	- Entrada, doble X, Y, las cantidades a comparar.
//
// Salida:
// double R8_MAX, el máximo de X e Y.

  double value;

  if(y < x){

    value = x;
  
  }else{

    value = y;
  }

  return value;
}

double r8_mop(int i){
//  Objetivo:
//	- R8_MOP devuelve la potencia I-ésima de -1 como un valor R8.
//
//  Discusión:
//	- Un R8 es un valor doble.
// Parámetros:
// 	Entrada, int I, la potencia de -1.
// 	Salida, doble R8_MOP, la i-ésima potencia de -1.
  
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
// Objetivo:
//
// R8_UNIFORM_01 devuelve una unidad pseudoaleatoria R8.
//
//  Discusión:
//
// Esta rutina implementa la recursividad
//
// semilla = ( 16807 * semilla ) mod ( 2^31 - 1 )
// u = semilla / ( 2^31 - 1 )
//
// La aritmética de enteros nunca requiere más de 32 bits,
// incluyendo un bit de signo.
//
// Si la semilla inicial es 12345, entonces los tres primeros cálculos son
// Entrada Salida R8_UNIFORM_01
// sed sed
//
// 12345 207482415 0,096616
// 207482415 1790989824 0,833995
// 1790989824 2035175616 0.947702
//
// Parámetros:
//
// Entrada/salida, int &SEED, el valor "sed". Normalmente, esto
// el valor no debe ser 0. En la salida, SEED se ha actualizado.
//
// Salida, doble R8_UNIFORM_01, una nueva variable pseudoaleatoria,
// estrictamente entre 0 y 1.

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
//  Objetivo:
//	R8MAT_IS_IDENTITY determina si un R8MAT es la identidad.
//
//  Discusión:
//	Un R8MAT es una matriz de valores reales (tipo = 8).
//	La rutina devuelve la norma de Frobenius de A - I.
//
// Parámetros:
//	Entrada, int N, el orden de la matriz.
//	Entrada, doble A[N*N], la matriz.
//
//	Salida, doble R8MAT_IS_IDENTITY, la norma Frobenius
// 	de la matriz de diferencias A - I, que sería exactamente cero
// 	si A fuera la matriz identidad.

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
//  Objetivo:
//	R8MAT_MM_NEW multiplica dos matrices.
//
//  Discusión:
//	Un R8MAT es una matriz doblemente dimensionada de valores R8, almacenada como un vector
//	en orden de columnas principales.
//
// 	Para esta rutina, el resultado se devuelve como el valor de la función.
//
// Parámetros:
//	Entrada, int N1, N2, N3, el orden de las matrices.
//	Entrada, double A[N1*N2], double B[N2*N3], las matrices a multiplicar.
//	Salida, doble R8MAT_MM[N1*N3], la matriz del producto C = A * B.


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
//  Objetivo:
//	R8MAT_MV_NEW multiplica una matriz por un vector.
//
//  Discusión:
//	Un R8MAT es una matriz doblemente dimensionada de valores R8, almacenada como un vector
// 	en orden de columnas principales.
//
//	Para esta rutina, el resultado se devuelve como el valor de la función.
//
// Parámetros:
//	Ingrese, int M, N, el número de filas y columnas de la matriz.
//	Entrada, doble A[M,N], la matriz M por N.
//	Entrada, doble X[N], el vector a multiplicar por A.
//	Salida, doble R8MAT_MV_NEW[M], el producto A*X.

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
//  Objetivo:
//	R8MAT_NORM_FRO devuelve la norma Frobenius de un R8MAT.
//
//  Discusión:
//	Un R8MAT es una matriz doblemente dimensionada de valores R8, almacenada como un vector
// 	en orden de columnas principales.
//
// La norma de Frobenius se define como
//
// R8MAT_NORM_FRO = raíz cuadrada(suma(1 <= I <= M ) * suma (1 <= j <= N) A(I,J)**2 )
// La norma matricial de Frobenius no se deriva de una norma vectorial, sino
// es compatible con la norma del vector L2, por lo que:
//
// r8vec_norm_l2 ( A * x ) <= r8mat_norm_fro ( A ) * r8vec_norm_l2 ( x ).
//
// Parámetros:
//	- Ingrese, int M, el número de filas en A.
//	- Ingrese, int N, el número de columnas en A.
//	- Entrada, doble A[M*N], la matriz cuyo Frobenius se desea la norma.
//	- Salida, doble R8MAT_NORM_FRO, la norma Frobenius de A.

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
//  Objetivo:
//	- R8MAT_PRINT imprime un R8MAT.
//
//  Discusión:
// 	- Un R8MAT es una matriz doblemente dimensionada de valores R8, almacenada como un vector
// 	en orden de columnas principales.
//
//	 La entrada A(I,J) se almacena como A[I+J*M]
// Parámetros:
// 	- Ingrese, int M, el número de filas en A.
// 	- Ingrese, int N, el número de columnas en A.
//	- Entrada, doble A[M*N], la matriz M por N.
//	- Entrada, cadena TITLE, un título.
  
    r8mat_print_some(m, n, a, 1, 1, m, n, title);

  return;
}

void r8mat_print_some(int m, int n, double a[], int ilo, int jlo, int ihi, int jhi, string title){
//  Objetivo:
//	- R8MAT_PRINT_SOME imprime algo de un R8MAT.
//
//  Discusión:
//	- Un R8MAT es una matriz doblemente dimensionada de valores R8, almacenada como un vector en orden de columnas principales.
//
// Parámetros:
//	- Ingrese, int M, el número de filas de la matriz. M debe ser positivo.
// 	- Entrada, int N, el número de columnas de la matriz. N debe ser positivo.
// 	- Entrada, doble A[M*N], la matriz.
//	- Ingrese, int ILO, JLO, IHI, JHI, designe la primera fila y columna, y la última fila y columna a imprimir.
//	- Entrada, cadena TITLE, un título.


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


  for(j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX){
    
    j2hi = j2lo + INCX - 1;
    j2hi = i4_min(j2hi, n);
    j2hi = i4_min(j2hi, jhi);

    cout << "\n";
	// Para cada columna J en el rango actual...
    	// Escribe el encabezado.

    cout << "  Col:    ";
    for(j = j2lo; j <= j2hi; j++){

      cout << setw(7) << j - 1 << "       ";
    }

    cout << "\n";
    cout << "  Row\n";
    cout << "\n";
    

    i2lo = i4_max(ilo, 1);
    i2hi = i4_min(ihi, m);

    for(i = i2lo; i <= i2hi; i++){
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
//  Objetivo:
// 	R8MAT_ZEROS_NEW devuelve un nuevo R8MAT puesto a cero.
//  Discusión:
//	Un R8MAT es una matriz doblemente dimensionada de valores R8, almacenada como un vector
// 	en orden de columnas principales.
// Parámetros:
//	- Ingrese, int M, N, el número de filas y columnas.
//	- Salida, doble R8MAT_ZEROS_NEW[M*N], la nueva matriz puesta a cero.

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
//  Objetivo:
//	R8VEC_DOT_PRODUCT calcula el producto escalar de un par de R8VEC.
//  Discusión:
// 	Un R8VEC es un vector de R8.
// Parámetros:
//	Ingrese, int N, el número de entradas en los vectores.
//	Entrada, doble A1[N], A2[N], los dos vectores a considerar.
//	Salida, doble R8VEC_DOT_PRODUCT, el producto escalar de los vectores.

  int i;
  double value;

  value = 0.0;
  for(i = 0; i < n; i++){
    
     value = value + a1[i] * a2[i];
  }

  return value;
}

double *r8vec_linspace_new(int n, double a_first, double a_last){
//  Objetivo:
//	 R8VEC_LINSPACE_NEW crea un vector de valores espaciados linealmente.
//  Discusión:
//	Un R8VEC es un vector de R8.
// Parámetros:
//	Entrada, int N, el número de entradas en el vector.
// 	Entrada, doble A_FIRST, A_LAST, la primera y la última entrada.
// 	Salida, doble R8VEC_LINSPACE_NEW[N], un vector de datos espaciados linealmente.

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
//  Objetivo:
// 	R8VEC_SUM devuelve la suma de un R8VEC.
//  Discusión:
// 	Un R8VEC es un vector de R8.
// Parámetros:
//	Entrada, int N, el número de entradas en el vector.
//	Entrada, doble A[N], el vector.
// 	Salida, doble R8VEC_SUM, la suma del vector.

  int i;
  double value;

  value = 0.0;
  
  for(i = 0; i < n; i++){

    value = value + a[i];
  }

  return value;
}

void timestamp(){
//  Objetivo:
//	TIMESTAMP imprime la fecha YMDHMS actual como una marca de tiempo.
//  Ejemplo:
// 	31 de mayo de 2001 09:45:54 a. m.
// Parámetros:
//	Ninguno

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
