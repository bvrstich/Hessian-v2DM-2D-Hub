/**
 * @mainpage 
 * This is an implementation of the dual only, potential reduction interior point method
 * for optimizing the second order density matrix using the P, Q, G and T N-representability conditions.
 * Compiling can be done with the options PQ, PQG, PQGT1, PQGT2 and PQGT (for all conditions active) with logical consequences for the program.
 * @author Brecht Verstichel, Ward Poelmans
 * @date 22-02-2010
 */
#include <iostream>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;

//includes all important headers and defines which conditions are
//going to be used:
#include "include.h"

/**
 * In the main the actual program is run.\n 
 * We start from the unity density matrix normed on the particle number and minimize the 
 * ojective function:\n\n
 * Tr (Gamma H) - t * ln(det P(Gamma)) \n\n
 * Once the minimum is found the parameter t is reduced and a new search is initiated,
 * this goes on until convergence is reached.\n
 * The potential is minimized using the Newton-Raphson method and the resulting linear system
 * is solved via the linear conjugate gradient method.
 */
int main(void) {

   srand(time(NULL));

   cout.precision(10);

   const int L = 6;//dimension of the lattice
   const int N = 36;//nr of particles

   Tools::init(L,N);

   Hamiltonian::init();

   TPM::init();
   PHM::init();
   DPM::init();
   PPHM::init();

   PPHM_ns::init();

   TPTPM::init();

   Gradient::init();

   PPHM pphm;
   pphm.fill_Random();

   int L2 = Tools::gL2();
   int L4 = L2*L2;
   int L6 = L4*L2;
   int L8 = L6*L2;
   int L10 = L8*L2;

   double *******ppharray0 = new double ****** [L2];

   for(int K = 0;K < L2;++K){

      ppharray0[K] = new double ***** [2];

      for(int S_ab = 0;S_ab < 2;++S_ab){

         ppharray0[K][S_ab] = new double **** [L2];

         for(int a = 0;a < L2;++a){

            ppharray0[K][S_ab][a] = new double *** [L2];

            for(int b = 0;b < L2;++b){

               ppharray0[K][S_ab][a][b] = new double ** [2];

               for(int S_de = 0;S_de < 2;++S_de){

                  ppharray0[K][S_ab][a][b][S_de] = new double * [L2];

               }

            }

         }

      }

   }

   ppharray0[0][0][0][0][0][0] = new double [4*L10];

   for(int K = 1;K < L2;++K)
      ppharray0[K][0][0][0][0][0] =  ppharray0[K - 1][0][0][0][0][0] + 4*L8;

   for(int K = 0;K < L2;++K)
      for(int S_ab = 1;S_ab < 2;++S_ab)
         ppharray0[K][S_ab][0][0][0][0] =  ppharray0[K][S_ab - 1][0][0][0][0] + 2*L8;

   for(int K = 0;K < L2;++K)
      for(int S_ab = 0;S_ab < 2;++S_ab)
         for(int a = 1;a < L2;++a)
            ppharray0[K][S_ab][a][0][0][0] =  ppharray0[K][S_ab][a - 1][0][0][0] + 2*L6;

   for(int K = 0;K < L2;++K)
      for(int S_ab = 0;S_ab < 2;++S_ab)
         for(int a = 0;a < L2;++a)
            for(int b = 1;b < L2;++b)
               ppharray0[K][S_ab][a][b][0][0] =  ppharray0[K][S_ab][a][b - 1][0][0] + 2*L4;

   for(int K = 0;K < L2;++K)
      for(int S_ab = 0;S_ab < 2;++S_ab)
         for(int a = 0;a < L2;++a)
            for(int b = 0;b < L2;++b)
               for(int S_de = 1;S_de < 2;++S_de)
                  ppharray0[K][S_ab][a][b][S_de][0] =  ppharray0[K][S_ab][a][b][S_de - 1][0] + L4;

   for(int K = 0;K < L2;++K)
      for(int S_ab = 0;S_ab < 2;++S_ab)
         for(int a = 0;a < L2;++a)
            for(int b = 0;b < L2;++b)
               for(int S_de = 0;S_de < 2;++S_de)
                  for(int d = 1;d < L2;++d)
                     ppharray0[K][S_ab][a][b][S_de][d] =  ppharray0[K][S_ab][a][b][S_de][d - 1] + L2;

   double *****ppharray1 = new double **** [L2];

   for(int K = 0;K < L2;++K){

      ppharray1[K] = new double *** [L2];

      for(int a = 0;a < L2;++a){

         ppharray1[K][a] = new double ** [L2];

         for(int b = 0;b < L2;++b){

            ppharray1[K][a][b] = new double * [L2];

         }

      }

   }

   ppharray1[0][0][0][0] = new double [L10];

   for(int K = 1;K < L2;++K)
      ppharray1[K][0][0][0] = ppharray1[K - 1][0][0][0] + L8;

   for(int K = 0;K < L2;++K)
      for(int a = 1;a < L2;++a)
         ppharray1[K][a][0][0] = ppharray1[K][a - 1][0][0] + L6;

   for(int K = 0;K < L2;++K)
      for(int a = 0;a < L2;++a)
         for(int b = 1;b < L2;++b)
            ppharray1[K][a][b][0] = ppharray1[K][a][b - 1][0] + L4;

   for(int K = 0;K < L2;++K)
      for(int a = 0;a < L2;++a)
         for(int b = 0;b < L2;++b)
            for(int d = 1;d < L2;++d)
               ppharray1[K][a][b][d] = ppharray1[K][a][b][d - 1] + L2;

   pphm.convert(ppharray0,ppharray1);

   TPTPM tpmm;
   tpmm.dpt2_pph(ppharray0,ppharray1);

   delete [] ppharray0[0][0][0][0][0][0];

   for(int K = 0;K < L2;++K){

      for(int S_ab = 0;S_ab < 2;++S_ab){

         for(int a = 0;a < L2;++a){

            for(int b = 0;b < L2;++b){

               for(int S_de = 0;S_de < 2;++S_de){

                  delete [] ppharray0[K][S_ab][a][b][S_de];

               }

               delete [] ppharray0[K][S_ab][a][b];

            }

            delete [] ppharray0[K][S_ab][a];

         }

         delete [] ppharray0[K][S_ab];

      }

      delete [] ppharray0[K];

   }

   delete [] ppharray0;

   delete [] ppharray1[0][0][0][0];

   for(int K = 0;K < L2;++K){

      for(int a = 0;a < L2;++a){

         for(int b = 0;b < L2;++b){

            delete [] ppharray1[K][a][b];

         }

         delete [] ppharray1[K][a];

      }

      delete [] ppharray1[K];

   }

   delete [] ppharray1;

/*
   int L2 = Tools::gL2();
   int L4 = L2*L2;
   int L6 = L4*L2;
   int L8 = L6*L2;

   double **ppharray = new double * [2*L2];

   for(int B = 0;B < L2;++B)//S = 1/2
      ppharray[B] = new double [4*L8];

   for(int B = L2;B < 2*L2;++B)//S = 3/2
      ppharray[B] = new double [L8];

   pphm.convert(ppharray);

   TPTPM tpmm;
   tpmm.dpt2_pph(ppharray);

   //remove the array
   for(int B = 0;B < 2*L2;++B)
      delete [] ppharray[B];

   delete [] ppharray;
*/
   /*
      Newton newton;

   //hamiltoniaan
   TPM ham;
   ham.hubbard(1.0);

   TPM rdm;
   rdm.unit();

   double t = 1.0;
   double tolerance = 1.0e-5;

   int tot_iter = 0;

   //outer iteration: scaling of the potential barrier
   while(t > 1.0e-12){

   cout << t << "\t" << rdm.trace() << "\t" << rdm.ddot(ham) << "\t";

   int nr_newton_iter = 0;

   double convergence = 1.0;

   //inner iteration: 
   //Newton's method for finding the minimum of the current potential
   while(convergence > tolerance){

   ++nr_newton_iter;

   SUP P;

   P.fill(rdm);

   P.invert();

   //fill the Newton object with the correct information, and solve for Delta
   newton.construct(t,ham,P);

   //dit wordt de stap:
   TPM delta;
   delta.convert(newton.gGradient());

   //line search
   double a = delta.line_search(t,P,ham);

   //rdm += a*delta;
   rdm.daxpy(a,delta);

   convergence = a*a*delta.ddot(delta);

   }

   cout << nr_newton_iter << endl;

   t /= 2.0;

   //what is the tolerance for the newton method?
   tolerance = 1.0e-5*t;

   if(tolerance < 1.0e-12)
   tolerance = 1.0e-12;

   tot_iter += nr_newton_iter;

   }

   cout << endl;

   cout << "Final Energy:\t" << ham.ddot(rdm) << endl;
   cout << endl;
   cout << "Final Spin:\t" << rdm.S_2() << endl;

   cout << endl;
   cout << "Total nr of Newton steps = " << tot_iter << endl;
   */
      Gradient::clear();

   TPTPM::clear();

   PPHM_ns::clear();

   PPHM::clear();
   DPM::clear();
   PHM::clear();
   TPM::clear();

   Hamiltonian::clear();

   Tools::clear();

   return 0;

}
