#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::vector;
using std::cout;
using std::endl;
using std::ios;

#include "include.h"

/**
 * standard constructor:
 */
SPSPM::SPSPM() : Matrix(Tools::gL2()) { }

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix tpmm_c
 * @param spmm_c object that will be copied into this.
 */
SPSPM::SPSPM(const SPSPM &spmm_c) : Matrix(spmm_c){ }

/**
 * destructor
 */
SPSPM::~SPSPM(){ }

ostream &operator<<(ostream &output,const SPSPM &spmm_p){

   for(int a = 0;a < Tools::gL2();++a)
      for(int e = 0;e < Tools::gL2();++e)
         output << a << "\t" << e << "\t" << spmm_p(a,e) << endl;

   return output;

}

/**
 * construct a SPSPM by doubly-tracing out the direct product of two TPM's
 */
void SPSPM::dpt2(double scale,const TPM &Q){

   for(int a = 0;a < Tools::gL2();++a)
      for(int e = a;e < Tools::gL2();++e){

         (*this)(a,e) = 0.0;

         //first S = 0
         for(int k = 0;k < Tools::gL2();++k){

            int l = Hamiltonian::gadjoint_sum(a,k,e);

            (*this)(a,e) += (Q(0,a,k,e,l) * Q(0,a,k,e,l) )/ ( TPM::gnorm(a,k) * TPM::gnorm(a,k) * TPM::gnorm(e,l) * TPM::gnorm(e,l) );

         }

         double ward = 0.0;

         //then S = 1
         for(int k = 0;k < Tools::gL2();++k){

            int l = Hamiltonian::gadjoint_sum(a,k,e);

            ward += Q(1,a,k,e,l) * Q(1,a,k,e,l);

         }

         (*this)(a,e) += 3.0 * ward;

         (*this)(a,e) *= scale;


      }

   this->symmetrize();

}

/**
 * construct the doubly-traced direct product of two PHM matrices
 */
void SPSPM::dpt2(double scale,double **pharray){

   int L2 = Tools::gL2();

   for(int a = 0;a < Tools::gL2();++a)
      for(int e = a;e < Tools::gL2();++e){

         (*this)(a,e) = 0.0;

         //first S = 0
         for(int k = 0;k < Tools::gL2();++k){

            int K = Hamiltonian::add(a,k);

            (*this)(a,e) += pharray[K][a + e*L2] * pharray[K][a + e*L2];

         }

         double ward = 0.0;

         //then S = 1
         for(int k = 0;k < Tools::gL2();++k){

            int K = Hamiltonian::add(a,k);

            ward += pharray[K + L2][a + e*L2] * pharray[K + L2][a + e*L2];

         }

         (*this)(a,e) += 3.0 * ward;

         (*this)(a,e) *= scale;


      }

   this->symmetrize();

}

/**
 * construct a SPSPM by quadruple tracing the direct product of two DPM's, input the array
 */
void SPSPM::dpt4(double scale,double **dparray){

   int L2 = Tools::gL2();
   int L4 = L2*L2;
   int L6 = L2*L4;
   int L8 = L2*L6;

   for(int a = 0;a < L2;++a)
      for(int e = a;e < L2;++e){

         (*this)(a,e) = 0.0;

         double ward = 0.0;

         //first S = 1/2
         for(int l = 0;l < L2;++l)
            for(int k = 0;k < L2;++k){

               int K = Hamiltonian::add(a,l,k);

               for(int S_al = 0;S_al < 2;++S_al)
                  for(int S_en = 0;S_en < 2;++S_en)
                     for(int n = 0;n < L2;++n){

                        ward += dparray[K][a + l*L2 + e*L4 + n*L6 + S_al*L8 + 2*S_en*L8] * dparray[K][a + l*L2 + e*L4 + n*L6 + S_al*L8 + 2*S_en*L8]

                           / ( TPM::gnorm(a,l) * TPM::gnorm(a,l) * TPM::gnorm(e,n) * TPM::gnorm(e,n) );


                     }

            }

         (*this)(a,e) = 2.0 * ward;

         ward = 0.0;

         //then S = 3/2
         for(int l = 0;l < L2;++l)
            for(int k = 0;k < L2;++k){

               int K = Hamiltonian::add(a,l,k);

               for(int n = 0;n < L2;++n)
                  ward += dparray[K + L2][a + l*L2 + e*L4 + n*L6] * dparray[K + L2][a + l*L2 + e*L4 + n*L6];

            }

         (*this)(a,e) += 4.0 * ward;

         //scale
         (*this)(a,e) *= 0.5 * scale;

      }

   this->symmetrize();

}

/**
 * construct a SPSPM by quadruple skew-tracing the direct product of two PPHM's, input the array
 */
void SPSPM::dpw4(double scale,double **ppharray){

   int L2 = Tools::gL2();
   int L4 = L2*L2;
   int L6 = L4*L2;
   int L8 = L6*L2;

   *this = 0.0;

   //first S = 1/2
   for(int S_kl = 0;S_kl < 2;++S_kl){

      for(int k = 0;k < L2;++k)
         for(int l = k + S_kl;l < L2;++l)
            for(int a = 0;a < L2;++a){

               int K_pph = Hamiltonian::add(k,l,a);

               for(int S_mn = 0;S_mn < 2;++S_mn){

                  for(int m = 0;m < L2;++m)
                     for(int n = m + S_mn;n < L2;++n){

                        int e = Hamiltonian::adjoint(K_pph,m,n);

                        (*this)(a,e) += 2.0 * scale * ppharray[K_pph][k + l*L2 + m*L4 + n*L6 + S_kl*L8 + 2*S_mn*L8] 
                        
                           * ppharray[K_pph][k + l*L2 + m*L4 + n*L6 + S_kl*L8 + 2*S_mn*L8];

                     }

               }

            }

   }

   //then S = 3/2
   for(int k = 0;k < L2;++k)
      for(int l = k + 1;l < L2;++l)
         for(int a = 0;a < L2;++a){

            int K_pph = Hamiltonian::add(k,l,a);

            for(int m = 0;m < L2;++m)
               for(int n = m + 1;n < L2;++n){

                  int e = Hamiltonian::adjoint(K_pph,m,n);

                  (*this)(a,e) += 4.0 * scale * ppharray[K_pph + L2][k + l*L2 + m*L4 + n*L6] * ppharray[K_pph + L2][k + l*L2 + m*L4 + n*L6];

               }

         }

}
