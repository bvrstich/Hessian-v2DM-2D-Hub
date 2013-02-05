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
TPSPM::TPSPM() : RecMat(TPTPM::gn(),Tools::gL2()) { }

/**
 * copy constructor: constructs RecMat object
 * @param tpspm_c object that will be copied into this.
 */
TPSPM::TPSPM(const TPSPM &tpspm_c) : RecMat(tpspm_c){ }

/**
 * destructor
 */
TPSPM::~TPSPM(){ }

/**
 * construct a TPSPM by tracing the direct product of two TPM's
 */
void TPSPM::dpt(double scale,const TPM &Q){

  int B,I,J;

  int S;

  int K,a,b,c,d;

   for(int i = 0;i < TPTPM::gn();++i){

      B = TPTPM::gtpmm2t(i,0);

      I = TPTPM::gtpmm2t(i,1);
      J = TPTPM::gtpmm2t(i,2);

      S = TPM::gblock_char(B,0);
      K = TPM::gblock_char(B,1);

      a = TPM::gt2s(B,I,0);
      b = TPM::gt2s(B,I,1);
      c = TPM::gt2s(B,J,0);
      d = TPM::gt2s(B,J,1);

      for(int k = 0;k < Tools::gL2();++k){

         int l = Hamiltonian::gadjoint(K,k);

         if(k == l)
            (*this)(i,k) = 4.0 * scale * Q(S,a,b,k,l) * Q(S,c,d,k,l);
         else
            (*this)(i,k) = 2.0 * scale * Q(S,a,b,k,l) * Q(S,c,d,k,l);

      }
   }

}

/**
 * construct the singly-traced antisymmetrized direct product of two PHM matrices
 */
void TPSPM::dpt(double scale,double **pharray){

   int L2 = Tools::gL2();

   int B,I,J;

   int S;

   int sign;

   int a,b,c,d;
   int c_,d_;

   for(int i = 0;i < TPTPM::gn();++i){

      B = TPTPM::gtpmm2t(i,0);

      I = TPTPM::gtpmm2t(i,1);
      J = TPTPM::gtpmm2t(i,2);

      S = TPM::gblock_char(B,0);

      sign = 1 - 2*S;

      a = TPM::gt2s(B,I,0);
      b = TPM::gt2s(B,I,1);
      c = TPM::gt2s(B,J,0);
      d = TPM::gt2s(B,J,1);

      c_ = Hamiltonian::gbar(c);
      d_ = Hamiltonian::gbar(d);

      for(int k = 0;k < Tools::gL2();++k){

         (*this)(i,k) = 0.0;

         for(int Z = 0;Z < 2;++Z){

            //(a,d,c,b)
            int P = Hamiltonian::add(a,d_);

            double ward = pharray[P + Z*L2][a + k*L2] * pharray[P + Z*L2][c + k*L2];

            //(b,d,c,a)
            P = Hamiltonian::add(b,d_);

            ward += sign * pharray[P + Z*L2][b + k*L2] * pharray[P + Z*L2][c + k*L2];

            //(a,c,d,b)
            P = Hamiltonian::add(a,c_);

            ward += sign * pharray[P + Z*L2][a + k*L2] * pharray[P + Z*L2][d + k*L2];

            //(b,c,d,a)
            P = Hamiltonian::add(b,c_);

            ward += pharray[P + Z*L2][b + k*L2] * pharray[P + Z*L2][d + k*L2];

            (*this)(i,k) += 2.0 * (2.0*Z + 1.0) * Tools::g6j(0,0,Z,S) * ward;

         }

         (*this)(i,k) *= scale;

      }

   }

}

/**
 * construct a TPSPM object by triple-tracing the direct product of two DPM objects, here the DPM has already been transformed to a double **array
 */
void TPSPM::dpt3(double scale,double **dparray){

   int L2 = Tools::gL2();
   int L4 = L2*L2;
   int L6 = L2*L4;
   int L8 = L2*L6;

   int B,I,J;

   int S;

   int a,b,c,d;

   for(int i = 0;i < TPTPM::gn();++i){

      B = TPTPM::gtpmm2t(i,0);

      I = TPTPM::gtpmm2t(i,1);
      J = TPTPM::gtpmm2t(i,2);

      S = TPM::gblock_char(B,0);

      a = TPM::gt2s(B,I,0);
      b = TPM::gt2s(B,I,1);
      c = TPM::gt2s(B,J,0);
      d = TPM::gt2s(B,J,1);

      for(int e = 0;e < L2;++e){

         double ward = 0.0;

         //first S = 1/2
         for(int k = 0;k < L2;++k){//abk and cdk

            int K = Hamiltonian::add(a,b,k);

            //first S_ep = 0
            for(int p = 0;p < e;++p)
               ward += dparray[K][a + b*L2 + e*L4 + p*L6 + S*L8] * dparray[K][c + d*L2 + e*L4 + p*L6 + S*L8];

            ward += 2.0 * dparray[K][a + b*L2 + e*L4 + e*L6 + S*L8] * dparray[K][c + d*L2 + e*L4 + e*L6 + S*L8];

            for(int p = e + 1;p < L2;++p)
               ward += dparray[K][a + b*L2 + e*L4 + p*L6 + S*L8] * dparray[K][c + d*L2 + e*L4 + p*L6 + S*L8];


            //then S_ep = 1
            for(int p = 0;p < e;++p)
               ward += dparray[K][a + b*L2 + e*L4 + p*L6 + S*L8 + 2*L8] * dparray[K][c + d*L2 + e*L4 + p*L6 + S*L8 + 2*L8];

            for(int p = e + 1;p < L2;++p)
               ward += dparray[K][a + b*L2 + e*L4 + p*L6 + S*L8 + 2*L8] * dparray[K][c + d*L2 + e*L4 + p*L6 + S*L8 + 2*L8];

         }

         (*this)(i,e) = 2.0/(2*S + 1.0) * ward;

         //then S = 3/2, only when
         if(S == 1){

            ward = 0.0;

            for(int k = 0;k < L2;++k){//abk and cdk

               int K = Hamiltonian::add(a,b,k);

               //only S_ep = 1 term possible
               for(int p = 0;p < e;++p)
                  ward += dparray[K + L2][a + b*L2 + e*L4 + p*L6] * dparray[K + L2][c + d*L2 + e*L4 + p*L6];

               for(int p = e + 1;p < L2;++p)
                  ward += dparray[K + L2][a + b*L2 + e*L4 + p*L6] * dparray[K + L2][c + d*L2 + e*L4 + p*L6];

            }

            (*this)(i,e) += 4.0/3.0 * ward;

         }

         (*this)(i,e) *= scale;

      }

   }

}

/**
 * construct a TPSPM object by once tracing, and twice skew-tracing the direct product of two PPHM objects, here the PPHM has already been transformed to a double **array
 */
void TPSPM::dptw2(double scale,double **ppharray){

   int L2 = Tools::gL2();
   int L4 = L2*L2;
   int L6 = L4*L2;
   int L8 = L6*L2;

   int B,I,J;

   int S;

   int a,b,c,d;

   for(int i = 0;i < TPTPM::gn();++i){

      B = TPTPM::gtpmm2t(i,0);

      I = TPTPM::gtpmm2t(i,1);
      J = TPTPM::gtpmm2t(i,2);

      S = TPM::gblock_char(B,0);

      a = TPM::gt2s(B,I,0);
      b = TPM::gt2s(B,I,1);
      c = TPM::gt2s(B,J,0);
      d = TPM::gt2s(B,J,1);

      for(int e = 0;e < L2;++e){

         double ward = 0.0;

         //first S = 1/2
         for(int S_kl = 0;S_kl < 2;++S_kl)
            for(int k = 0;k < L2;++k)
               for(int l = k + S_kl;l < L2;++l){

                  int K_pph = Hamiltonian::add(k,l,e);

                  ward += ppharray[K_pph][a + b*L2 + k*L4 + l*L6 + S*L8 + 2*S_kl*L8] * ppharray[K_pph][c + d*L2 + k*L4 + l*L6 + S*L8 + 2*S_kl*L8];

               }

         (*this)(i,e) = 2.0/(2.0*S + 1.0) * ward;

         if(S == 1){

            ward = 0.0;

            for(int k = 0;k < L2;++k)
               for(int l = k + 1;l < L2;++l){

                  int K_pph = Hamiltonian::add(k,l,e);

                  //first S_kl = 0
                  ward += ppharray[K_pph + L2][a + b*L2 + k*L4 + l*L6] * ppharray[K_pph + L2][c + d*L2 + k*L4 + l*L6];

               }

            (*this)(i,e) += 4.0/3.0 * ward;

         }

         (*this)(i,e) *= scale;

      }
   }

}

/**
 * construct a TPSPM object by triple skew-tracing the direct product of two PPHM objects, here the PPHM has already been transformed to a double **array
 */
void TPSPM::dpw3(double scale,double **ppharray){

   int L2 = Tools::gL2();
   int L4 = L2*L2;
   int L6 = L4*L2;
   int L8 = L6*L2;

   int B,I_i,J_i;

   int S;

   int a,b,c,d;
   int a_,b_;

   double ward;
   int k,K_pph;

   int sign;

   for(int i = 0;i < TPTPM::gn();++i){

      B = TPTPM::gtpmm2t(i,0);

      I_i = TPTPM::gtpmm2t(i,1);
      J_i = TPTPM::gtpmm2t(i,2);

      S = TPM::gblock_char(B,0);

      sign = 1 - 2*S;

      a = TPM::gt2s(B,I_i,0);
      b = TPM::gt2s(B,I_i,1);
      c = TPM::gt2s(B,J_i,0);
      d = TPM::gt2s(B,J_i,1);

      a_ = Hamiltonian::gbar(a);
      b_ = Hamiltonian::gbar(b);

      for(int e = 0;e < L2;++e){

         int e_ = Hamiltonian::gbar(e);

         (*this)(i,e) = 0.0;

         //S'' = 1/2
         for(int J = 0;J < 2;++J){

            ward = 0.0;

            for(int S_mn = 0;S_mn < 2;++S_mn)
               for(int m = 0;m < L2;++m)
                  for(int n = m + S_mn;n < L2;++n){

                     K_pph = Hamiltonian::add(m,n,e_);

                     //1)
                     k = Hamiltonian::adjoint(K_pph,d,a_);

                     ward += ppharray[K_pph][m + n*L2 + k*L4 + d*L6 + S_mn*L8 + 2*J*L8] * ppharray[K_pph][m + n*L2 + k*L4 + b*L6 + S_mn*L8 + 2*J*L8];

                     //2)
                     k = Hamiltonian::adjoint(K_pph,d,b_);

                     ward += sign * ppharray[K_pph][m + n*L2 + k*L4 + d*L6 + S_mn*L8 + 2*J*L8] * ppharray[K_pph][m + n*L2 + k*L4 + a*L6 + S_mn*L8 + 2*J*L8];

                     //3)
                     k = Hamiltonian::adjoint(K_pph,c,a_);

                     ward += sign * ppharray[K_pph][m + n*L2 + k*L4 + c*L6 + S_mn*L8 + 2*J*L8] * ppharray[K_pph][m + n*L2 + k*L4 + b*L6 + S_mn*L8 + 2*J*L8];

                     //4)
                     k = Hamiltonian::adjoint(K_pph,c,b_);

                     ward += ppharray[K_pph][m + n*L2 + k*L4 + c*L6 + S_mn*L8 + 2*J*L8] * ppharray[K_pph][m + n*L2 + k*L4 + a*L6 + S_mn*L8 + 2*J*L8];

                  }

            (*this)(i,e) += 2.0 * (2.0*J + 1.0) * Tools::g6j(0,0,S,J) * ward;

         }

         ward = 0.0;

         //S'' = 3/2
         for(int m = 0;m < L2;++m)
            for(int n = m + 1;n < L2;++n){

               K_pph = Hamiltonian::add(m,n,e_);

               //1)
               k = Hamiltonian::adjoint(K_pph,d,a_);

               ward += ppharray[K_pph + L2][m + n*L2 + k*L4 + d*L6] * ppharray[K_pph + L2][m + n*L2 + k*L4 + b*L6];

               //2)
               k = Hamiltonian::adjoint(K_pph,d,b_);

               ward += sign * ppharray[K_pph + L2][m + n*L2 + k*L4 + d*L6] * ppharray[K_pph + L2][m + n*L2 + k*L4 + a*L6];

               //3)
               k = Hamiltonian::adjoint(K_pph,c,a_);

               ward += sign * ppharray[K_pph + L2][m + n*L2 + k*L4 + c*L6] * ppharray[K_pph + L2][m + n*L2 + k*L4 + b*L6];

               //4)
               k = Hamiltonian::adjoint(K_pph,c,b_);

               ward += ppharray[K_pph + L2][m + n*L2 + k*L4 + c*L6] * ppharray[K_pph + L2][m + n*L2 + k*L4 + a*L6];

            }

         (*this)(i,e) += 4.0 * 3.0 * Tools::g6j(0,0,S,1) * ward;

         (*this)(i,e) *= scale;

      }
   }

}
