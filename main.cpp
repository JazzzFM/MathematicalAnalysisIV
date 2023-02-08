# include <cmath>
# include <cstdlib>
# include <iomanip>
# include <iostream>
# include "bernstein_polynomial.hpp"
# include <string>
#include <chrono>
#include <thread>

using namespace std;

int main();
void bernstein_matrix_test();
void bernstein_matrix_test2();
void bernstein_matrix_determinant_test();
void bernstein_matrix_inverse_test();
void bernstein_poly_01_test();
void bernstein_poly_01_test2();
void bernstein_poly_01_matrix_test();
void bernstein_poly_ab_test();
void bernstein_poly_ab_approx_test();
void bernstein_to_legendre_test();
void bernstein_to_power_test();
void bernstein_vandermonde_test();

string var;

int main ( ){
//  Objetivo: MAIN es el programa principal de BERNSTEIN_POLYNOMIAL_PRB.
  
	timestamp();
	cout << "\n";
  	cout << "BERNSTEIN_POLYNOMIAL_PRB\n";
  	cout << "  Test de la biblioteca BERNSTEIN_POLYNOMIAL.\n";
  
  	std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  	bernstein_matrix_test();
	std::this_thread::sleep_for(std::chrono::milliseconds(1000));

  	bernstein_matrix_determinant_test();
  	std::this_thread::sleep_for(std::chrono::milliseconds(2000));

  	std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  	bernstein_matrix_inverse_test();
  	std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  
  	std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  	bernstein_poly_01_test();
  	std::this_thread::sleep_for(std::chrono::milliseconds(100));
  
  	std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  	bernstein_poly_01_test2();
  	std::this_thread::sleep_for(std::chrono::milliseconds(2000));
  
  	std::this_thread::sleep_for(std::chrono::milliseconds(1000));
  	bernstein_poly_01_matrix_test();

  	std::this_thread::sleep_for(std::chrono::milliseconds(1000));  
  	bernstein_poly_ab_test(); 
  	std::this_thread::sleep_for(std::chrono::milliseconds(2000));

  	bernstein_poly_ab_approx_test();
  	td::this_thread::sleep_for(std::chrono::milliseconds(1000));
  
  	bernstein_to_legendre_test();
  	std::this_thread::sleep_for(std::chrono::milliseconds(3000));

 	bernstein_to_power_test();
	std::this_thread::sleep_for(std::chrono::milliseconds(3000));
  
	bernstein_vandermonde_test();
 	std::this_thread::sleep_for(std::chrono::milliseconds(3000));
  
  	cout << "\n";
  	cout << "BERNSTEIN_POLYNOMIAL_PRB\n";
  	cout << " Fin de ejecución.\n";
  	cout << " \n";
  	timestamp();

  return 0;
}

void bernstein_matrix_test(){
//  Propósito: BERNSTEIN_MATRIX_TEST prueba BERNSTEIN_MATRIX.

  double *a;
  int n;

  cout << "\n";
  cout << "BERNSTEIN_MATRIX_TEST\n";
   cout << " BERNSTEIN_MATRIX devuelve una matriz A que transforma a\n";
   cout << " vector de coeficiente polinomial de la base de potencia a\n";
   cout << " la base de Bernstein.\n";
  n = 5;
  a = bernstein_matrix (n);
  r8mat_print(n, n, a, " Matriz Bernstein A de dimensión  5:");

  delete [] a;

  return;
}

void bernstein_matrix_test2(){
// Propósito: BERNSTEIN_MATRIX_TEST2 prueba BERNSTEIN_MATRIX.
//  Discusión:
// 	Aquí usamos la matriz de Bernstein para describir 
// 	un polinomio de Bernstein en términos de los monomios clásicos.

  double *a;
  double *ax;
  int i;
  int k;
  int n;
  double *x;
    
   cout << "\n";
   cout << "PRUEBA06\n";
   cout << " BERNSTEIN_MATRIX devuelve una matriz A que\n";
   cout << " transforma un vector coeficiente polinomial\n";
   cout << " de la base de Bernstein a la base de potencia.\n";
   cout << " Podemos usar esto para obtener valores explícitos de\n";
   cout << " Coeficientes polinómicos de Bernstein de cuarto grado como\n";
   cout << "\n";
   cout << " b(4,K)(X) = C4 * x^4\n";
   cout << " + C3 * x^3\n";
   cout << " + C2 * x^2\n";
   cout << " + C1 * x\n";
   cout << " + C0 * 1\n";

  n = 5;
  cout << "\n";
  cout << "     K       C4           C3            C2";
  cout << "            C1             C0\n";
  cout << "\n";

  a = bernstein_matrix(n);
  x = new double[n];

  for(k = 0; k < n; k++){
  	
	for(i = 0; i < n; i++){
      	    x[i] = 0.0;
         }
    
    x[k] = 1.0;
    ax = r8mat_mv_new(n, n, a, x );

    cout << "  " << setw(4) << k << "  ";
    for(i = 0; i < n; i++){

      cout << "  " << setw(14) << ax[i];
    }
    
	cout << "\n";
  }

  delete [] a;
  delete [] ax;
  delete [] x;

  return;
}

void bernstein_matrix_determinant_test(){
//  Propósito:  
//  	BERNSTEIN_MATRIX_DETERMINANT_TEST prueba BERNSTEIN_MATRIX_DETERMINANT.


  double *a;
  double a_norm_frobenius;
  double d1;
  int n;
   

   cout << "\n";
   cout << "BERNSTEIN_MATRIX_DETERMINANT_TEST\n";
   cout << " BERNSTEIN_MATRIX_DETERMINANT calcula el determinante de\n";
   cout << " la matriz de Bernstein.\n";
   cout << "\n";
   cout << " N ||A|| det(A)\n";
   cout << " calculado\n";
   cout << "\n";

  for(n = 5; n <= 15; n++){

    a = bernstein_matrix(n);
    
    a_norm_frobenius = r8mat_norm_fro(n, n, a);

    d1 = bernstein_matrix_determinant ( n );

    cout << "  " << setw(4) << n
         << "  " << setw(14) << a_norm_frobenius
         << "  " << setw(14) << d1 << "\n";

    delete [] a;
  }

  return;
}

void bernstein_matrix_inverse_test(){
//  Propósito:
//    BERNSTEIN_MATRIX_INVERSE_TEST prueba BERNSTEIN_MATRIX_INVERSE.

  double *a;
  double a_norm_frobenius;
  double *b;
  double b_norm_frobenius;
  double *c;
  double error_norm_frobenius;
  int n;

   cout << "\n";
   cout << "BERNSTEIN_MATRIX_INVERSE_TEST\n";
   cout << " BERNSTEIN_MATRIX_INVERSE calcula el inverso de\n";
   cout << " Matriz de Bernstein A.\n";
   cout << "\n";
   cout << " N ||A|| ||inv(A)|| ||I-A*inv(A)||\n";
   cout << "\n";
  
 for(n = 5; n <= 15; n++){
    
    a = bernstein_matrix(n);
    
    a_norm_frobenius = r8mat_norm_fro(n, n, a);

    b = bernstein_matrix_inverse(n);
    
    b_norm_frobenius = r8mat_norm_fro(n, n, b);

    c = r8mat_mm_new(n, n, n, a, b);

    error_norm_frobenius = r8mat_is_identity(n, c);

    cout << "  " << setw(4) << n
         << "  " << setw(14) << a_norm_frobenius
         << "  " << setw(14) << b_norm_frobenius
         << "  " << setw(14) << error_norm_frobenius << "\n";

    delete [] a;
    delete [] b;
    delete [] c;
  }
  return;
}

void bernstein_poly_01_test(){
//  Propósito:  BERNSTEIN_POLY_01_TEST prueba BERNSTEIN_POLY_01.

  double b;
  double *bvec;
  int k;
  int n;
  int n_data;
  double x;

   cout << "\n";
   cout << "BERNSTEIN_POLY_01_TEST:\n";
   cout << " BERNSTEIN_POLY_01 evalúa los polinomios de Bernstein\n";
   cout << " basado en el intervalo [0,1].\n";
   cout << "\n";
   cout << " N K X Exacto BP01(N,K)(X)\n";
   cout << "\n";

  n_data = 0;

  while(true){
    
    bernstein_poly_01_values(n_data, n, k, x, b);

    if (n_data == 0){
      break;
    }

    bvec = bernstein_poly_01(n, x);

    cout << "  " << setw(4) << n
         << "  " << setw(4) << k
         << "  " << setw(7) << x
         << "  " << setw(14) << b
         << "  " << setw(14) << bvec[k] << "\n";

    delete [] bvec;
  }

  return;
}

void bernstein_poly_01_test2(){
// Propósito: 
// 	BERNSTEIN_POLY_01_TEST2 prueba BERNSTEIN_POLY_01.
// Discusión:
// 	 Aquí probamos la propiedad Partition-of-Unity.
 
  double *bvec;
  int n;
  int n_data;
  int seed;
  double x;

   cout << "\n";
   cout << "BERNSTEIN_POLY_01_TEST2:\n";
   cout << " BERNSTEIN_POLY_01 evalúa los polinomios de Bernstein\n";
   cout << " basado en el intervalo [0,1].\n";
   cout << "\n";
   cout << " Aquí probamos la partición de la propiedad: \n";
   cout << "\n";
   cout << " N X Suma (0 <= K <= N) BP01(N,K)(X)\n";
   cout << "\n";

  seed = 123456789;

  for(n = 0; n <= 10; n++){
    
    x = r8_uniform_01 ( seed );

    bvec = bernstein_poly_01 ( n, x );

    cout << "  " << setw(4) << n
         << "  " << setw(7) << x
         << "  " << setw(14) << r8vec_sum ( n + 1, bvec ) << "\n";

    delete [] bvec;
  }
  return;
}

void bernstein_poly_01_matrix_test(){
// Propósito: 
// 	BERNSTEIN_POLY_01_MATRIX_TEST prueba BERNSTEIN_POLY_01_MATRIX.

  double *b;
  int m;
  int n;
  double *x;

   cout << "\n";
   cout << "BERNSTEIN_POLY_01_MATRIX_TEST\n";
   cout << " BERNSTEIN_POLY_01_MATRIX recibe M valores de datos X,\n";
   cout << " y un grado N, y devuelve una matriz B Mx(N+1) tal que\n";
   cout << " B(i,j) es el j-ésimo polinomio de Bernstein evaluado en el.\n";
   cout << "i-ésimo valor de datos.\n";

  m = 5;
  x = r8vec_linspace_new(m, 0.0, 1.0);
  n = 1;
  b = bernstein_poly_01_matrix(m, n, x);
  r8mat_print(m, n + 1, b, "  B(5,1+1):");
  delete [] b;
  delete [] x;

  m = 5;
  x = r8vec_linspace_new(m, 0.0, 1.0);
  n = 4;
  b = bernstein_poly_01_matrix(m, n, x);
  r8mat_print(m, n + 1, b, "  B(5,4+1):");
  delete [] b;
  delete [] x;

  m = 10;
  x = r8vec_linspace_new(m, 0.0, 1.0);
  n = 4;
  b = bernstein_poly_01_matrix(m, n, x);
  r8mat_print (m, n + 1, b, "  B(10,4+1):");
  delete [] b;
  delete [] x;

  m = 3;
  x = r8vec_linspace_new(m, 0.0, 1.0);
  n = 5;
  b = bernstein_poly_01_matrix(m, n, x);
  r8mat_print (m, n + 1, b, "  B(3,5+1):");
  delete [] b;
  delete [] x;

  return;
}

void bernstein_poly_ab_test(){
//  Propósito: 
//  	BERNSTEIN_POLY_AB_TEST prueba BERNSTEIN_POLY_AB.

  double a;
  double b;
  double *bern;
  int k;
  int n = 10;
  double x;
  

  cout << "\n";
  cout << "BERNSTEIN_POLY_AB_TEST\n";
  cout << " BERNSTEIN_POLY_AB evalúa los polinomios de Bernstein sobre un\n";
  cout << " intervalo arbitrario [A,B].\n";
  cout << "\n";
  cout << " Aquí, demostramos que \n";
  cout << " BPAB(N,K,A1,B1)(X1) = BPAB(N,K,A2,B2)(X2)\n";
  cout << "siempre que solo eso\n";
  cout << " (X1-A1)/(B1-A1) = (X2-A2)/(B2-A2).\n";

  x = 0.3;
  a = 0.0;
  b = 1.0;

  bern = bernstein_poly_ab(n, a, b, x);
 
  cout << "\n";
  cout << "     N     K     A        B        X       BPAB(N,K,A,B)(X)\n";
  cout << "\n";

  for(k = 0; k <= n; k++){
    
  cout << "  " << setw(4) << n
         << "  " << setw(4) << k
         << "  " << setw(7) << a
         << "  " << setw(7) << b
         << "  " << setw(7) << x
         << "  " << setw(14) << bern[k] << "\n";
  }

  delete [] bern;
 
  x = 1.3;
  a = 1.0;
  b = 2.0;
  bern = bernstein_poly_ab(n, a, b, x);
 
  cout << "\n";
  cout << "     N     K     A        B        X       BPAB(N,K,A,B)(X)\n";
  cout << "\n"; 

  for(k = 0; k <= n; k++){
    
     cout << "  " << setw(4) << n
         << "  " << setw(4) << k
         << "  " << setw(7) << a
         << "  " << setw(7) << b
         << "  " << setw(7) << x
         << "  " << setw(14) << bern[k] << "\n";
  }

  delete [] bern;

  x = 2.6;
  a = 2.0;
  b = 4.0;
  
  bern = bernstein_poly_ab(n, a, b, x);
 
  cout << "\n";
  cout << "     N     K     A        B        X       BPAB(N,K,A,B)(X)\n";
  cout << "\n";
 
  for (k = 0; k <= n; k++){
    
      cout << "  " << setw(4) << n
         << "  " << setw(4) << k
         << "  " << setw(7) << a
         << "  " << setw(7) << b
         << "  " << setw(7) << x
         << "  " << setw(14) << bern[k] << "\n";
  }

  delete [] bern;

  return;
}

void bernstein_poly_ab_approx_test(){
//  Propósito: 
//  	BERNSTEIN_POLY_AB_APPROX_TEST prueba BERNSTEIN_POLY_AB_APPROX.

  double a;
  double b;
  double error_max;
  int i;
  int maxdata = 20;
  int ndata;
  int nsample;
  int nval = 501;
  double *xdata;
  double *xval;
  double *ydata;
  double *yval;

  cout << "\n";
  cout << "BERNSTEIN_POLY_AB_APPROX_TEST\n";
  cout << "BERNSTEIN_POLY_AB_APPROX evalúa el polinomio de Bernstein\n";
  cout << " aproximado a una funcion F(X).\n";
  
  a = 1.0;
  b = 3.0;

  cout << "\n";
  cout << "     N      Max Error\n";
  cout << "\n";

  for(ndata = 0; ndata <= maxdata; ndata++){
    xdata = new double[ndata+1];
    ydata = new double[ndata+1];
    
    for(i = 0; i <= ndata; i++){
      
      if(ndata == 0){
        
          xdata[i] = 0.5*(a+b);
      }else{
        
           xdata[i] = ((double)(ndata - i)*a + 
			(double)(i)*b)/(double)(ndata);
      }

      ydata[i] = sin(xdata[i]);
    }
    
    xval = r8vec_linspace_new(nval, a, b);

    error_max = 0.0;

    yval = bernstein_poly_ab_approx(ndata, a, b, ydata, nval, xval);

    error_max = 0.0;
    
    for(i = 0; i < nval; i++){
      error_max = r8_max ( error_max, fabs ( yval[i] - sin ( xval[i] ) ) );
    
    }
    
    cout << "  " << setw(4) << ndata
         << "  " << setw(14) << error_max << "\n";

    delete [] xdata;
    delete [] xval;
    delete [] ydata;
    delete [] yval;
  }
  return;
}

void bernstein_to_legendre_test(){
//  Propósito:  
//  	BERNSTEIN_TO_LEGENDRE_TEST prueba BERNSTEIN_TO_LEGENDRE

  double *a;
  double *b;
  double *c;
  double e;
  int n = 5;
  
  cout << "\n";
  cout << "BERNSTEIN_TO_LEGENDRE_TEST:\n";
  cout << " BERNSTEIN_TO_LEGENDRE devuelve la matriz A que mapea\n";
  cout << " coeficientes polinómicos de Bernstein a la forma de Legendre.\n";

  a = bernstein_to_legendre(n);
  r8mat_print(n + 1, n + 1, a, "  A = bernstein_to_legendre(5):" );

  b = legendre_to_bernstein(n);
  r8mat_print(n + 1, n + 1, b, "  B = legendre_to_bernstein(5):" );

  c = r8mat_mm_new(n + 1, n + 1, n + 1, a, b);
  e = r8mat_is_identity(n + 1, c);

  cout << "\n";
  cout << "  ||A*B-I|| = " << e << "\n";

  delete [] a;
  delete [] b;
  delete [] c;

  return;
}

void bernstein_to_power_test(){
//  Propósito:
//  	 BERNSTEIN_TO_POWER_TEST prueba BERNSTEIN_TO_POWER.
  
  double *a;
  double *b;
  double *c;
  double e;
  int n = 5;

  cout << "\n";
  cout << "BERNSTEIN_TO_POWER_TEST:\n";
  cout << " BERNSTEIN_TO_POWER devuelve la matriz A que mapea\n";
  cout << " coeficientes polinómicos de Bernstein a la forma de Power.\n";

  a = bernstein_to_power(n);
  r8mat_print(n + 1, n + 1, a, "  A = bernstein_to_power(5):");

  b = power_to_bernstein(n);
  r8mat_print(n + 1, n + 1, b, "  B = power_to_bernstein(5):");

  c = r8mat_mm_new(n + 1, n + 1, n + 1, a, b);
  e = r8mat_is_identity(n + 1, c);
  
  cout << "\n";
  cout << "  ||A*B-I|| = " << e << "\n";

  delete [] a;
  delete [] b;
  delete [] c;

  return;
}

void bernstein_vandermonde_test(){
//  Propósito:  
//  	BERNSTEIN_VANDERMONDE_TEST prueba BERNSTEIN_VANDERMOND

  double *a;
  int n;

cout << "\n";
   cout << "BERNSTEIN_VANDERMONDE_TEST\n";
   cout << " BERNSTEIN_VANDERMONDE devuelve una matriz NxN cuya entrada (I,J)\n";
   cout << " es el valor del J-ésimo polinomio de Bernstein de grado N-1\n";
   cout << " evaluado en el i-ésimo punto igualmente espaciado en [0,1].\n";
  n = 8;
  a = bernstein_vandermonde(n);
  r8mat_print(n, n, a, "------> Bernstein Vandermonde(8): ");

  delete [] a;

  return;
}
