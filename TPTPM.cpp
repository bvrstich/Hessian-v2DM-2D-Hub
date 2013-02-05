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

int ***TPTPM::t2tpmm;
vector< vector<int> > TPTPM::tpmm2t;

/**
 * initialize the static lists
 */
void TPTPM::init(){

   t2tpmm = new int ** [Tools::gM()];

   for(int B = 0;B < Tools::gM();++B){

      t2tpmm[B] = new int * [TPM::gdim(B)];

      for(int i = 0;i < TPM::gdim(B);++i)
      t2tpmm[B][i] = new int [TPM::gdim(B)];

   }

   vector<int> v(3);

   int tpmm = 0;

   for(int B = 0;B < Tools::gM();++B){

      for(int i = 0;i < TPM::gdim(B);++i)
         for(int j = i;j < TPM::gdim(B);++j){

            v[0] = B;
            v[1] = i;
            v[2] = j;

            tpmm2t.push_back(v);

            t2tpmm[B][i][j] = tpmm;
            t2tpmm[B][j][i] = tpmm;

            ++tpmm;

         }

   }

}

/**
 * deallocate the static lists
 */
void TPTPM::clear(){

   for(int B = 0;B < Tools::gM();++B){

      for(int i = 0;i < TPM::gdim(B);++i)
         delete [] t2tpmm[B][i];

      delete [] t2tpmm[B];

   }

   delete [] t2tpmm;

}

/**
 * standard constructor:
 */
TPTPM::TPTPM() : Matrix(tpmm2t.size()) { }

/**
 * copy constructor: constructs Matrix object of dimension M*(M - 1)/2 and fills it with the content of matrix tpmm_c
 * @param tpmm_c object that will be copied into this.
 */
TPTPM::TPTPM(const TPTPM &tpmm_c) : Matrix(tpmm_c){ }

/**
 * destructor
 */
TPTPM::~TPTPM(){ }

/**
 * access the elements of the matrix in tp mode
 * @param B block index of the first two indices
 * @param I first tp index that forms the tpmm row index i together with J
 * @param J second tp index that forms the tpmm row index i together with I
 * @param B_ block index of the second two indices
 * @param K first tp index that forms the tpmm column index j together with L
 * @param L second tp index that forms the tpmm column index j together with K
 * @return the number on place TPTPM(i,j)
 */
double TPTPM::operator()(int B,int I,int J,int B_,int K,int L) const{

   int i = t2tpmm[B][I][J];
   int j = t2tpmm[B_][K][L];

   return (*this)(i,j);

}

ostream &operator<<(ostream &output,const TPTPM &tpmm_p){

   int B,I,J,B_,K,L;

   int a,b,c,d;
   int e,z,t,h;

   for(int i = 0;i < TPTPM::gn();++i){

      B = tpmm_p.tpmm2t[i][0];
      I = tpmm_p.tpmm2t[i][1];
      J = tpmm_p.tpmm2t[i][2];

      a = TPM::gt2s(B,I,0);
      b = TPM::gt2s(B,I,1);

      c = TPM::gt2s(B,J,0);
      d = TPM::gt2s(B,J,1);

      for(int j = i;j < TPTPM::gn();++j){

         B_ = tpmm_p.tpmm2t[j][0]; 
         K = tpmm_p.tpmm2t[j][1];
         L = tpmm_p.tpmm2t[j][2];

         e = TPM::gt2s(B_,K,0);
         z = TPM::gt2s(B_,K,1);

         t = TPM::gt2s(B_,L,0);
         h = TPM::gt2s(B_,L,1);

         output << i << "\t" << j << "\t|\t(" << B << ")\t" << I << "\t" << J << "\t(" << B_ << ")\t" << K << "\t" << L << "\t|\t" << 

            "(" << a << "," << b << "," << c << "," << d << ")\t(" << e << "," << z << "," << t << "," << h << ")\t|\t" << tpmm_p(i,j) << endl;

      }

   }

   return output;

}

/**
 * @return the dimension of a TPTPM matrix
 */
int TPTPM::gn(){

   return tpmm2t.size();

}

/**
 * access to the lists from outside the class
 */
int TPTPM::gt2tpmm(int B,int I,int J){

   return t2tpmm[B][I][J];

}

/**
 * access to the lists from outside the class
 * @param option == 0 return B, == 1 return a, == 2 return b
 */
int TPTPM::gtpmm2t(int i,int option){

   return tpmm2t[i][option];

}

/**
 * construct the antisymmetrized direct product of two PHM matrices
 */
void TPTPM::dp(double **pharray){

   int L2 = Tools::gL2();

   int B,B_;

   int sign,sign_;

   int a,b,c,d;
   int e,z,t,h;

   int d_;
   int t_,h_;

   int I_i,J_i,K_i,L_i;

   int S,S_;

   for(int i = 0;i < gn();++i){

      B = tpmm2t[i][0];

      S = TPM::gblock_char(B,0);

      sign = 1 - 2*S;

      I_i = tpmm2t[i][1];
      J_i = tpmm2t[i][2];

      a = TPM::gt2s(B,I_i,0);
      b = TPM::gt2s(B,I_i,1);
      c = TPM::gt2s(B,J_i,0);
      d = TPM::gt2s(B,J_i,1);

      d_ = Hamiltonian::gbar(d); 

      for(int j = i;j < gn();++j){

         B_ = tpmm2t[j][0];

         S_ = TPM::gblock_char(B_,0);

         sign_ = 1 - 2*S_;

         K_i = tpmm2t[j][1];
         L_i = tpmm2t[j][2];

         e = TPM::gt2s(B_,K_i,0);
         z = TPM::gt2s(B_,K_i,1);
         t = TPM::gt2s(B_,L_i,0);
         h = TPM::gt2s(B_,L_i,1);

         t_ = Hamiltonian::gbar(t); 
         h_ = Hamiltonian::gbar(h); 

         (*this)(i,j) = 0.0;

         int P = Hamiltonian::add(a,d_);
         int P_ = Hamiltonian::gbar(P);

         //(a,d,c,b)_(e,h,t,z) and (b,c,d,a)_(z,t,h,e)
         if(P == Hamiltonian::add(e,h_)){

            for(int Z = 0;Z < 2;++Z){

               (*this)(i,j) += (2*Z + 1.0) * Tools::g6j(0,0,Z,S) * Tools::g6j(0,0,Z,S_) * ( pharray[P + Z*L2][a + e*L2] * pharray[P + Z*L2][c + t*L2]

                     + pharray[P + Z*L2][a + t*L2]* pharray[P + Z*L2][c + e*L2] + pharray[P_ + Z*L2][b + z*L2] * pharray[P_ + Z*L2][d + h*L2]

                     + pharray[P_ + Z*L2][b + h*L2] * pharray[P_ + Z*L2][d + z*L2] ); 

            }

         }

         //(a,d,c,b)_(z,h,t,e) and (b,c,d,a)_(e,t,h,z)
         if(P == Hamiltonian::add(z,h_)){

            for(int Z = 0;Z < 2;++Z){

               (*this)(i,j) += sign_ * (2*Z + 1.0) * Tools::g6j(0,0,Z,S) * Tools::g6j(0,0,Z,S_) * ( pharray[P + Z*L2][a + z*L2] * pharray[P + Z*L2][c + t*L2]

                     + pharray[P + Z*L2][a + t*L2]* pharray[P + Z*L2][c + z*L2] + pharray[P_ + Z*L2][b + e*L2] * pharray[P_ + Z*L2][d + h*L2]

                     + pharray[P_ + Z*L2][b + h*L2] * pharray[P_ + Z*L2][d + e*L2] ); 

            }

         }

         //(a,d,c,b)_(e,t,h,z) and (b,c,d,a)_(z,h,t,e)
         if(P == Hamiltonian::add(e,t_) ){

            for(int Z = 0;Z < 2;++Z){

               (*this)(i,j) += sign_ * (2*Z + 1.0) * Tools::g6j(0,0,Z,S) * Tools::g6j(0,0,Z,S_) * ( pharray[P + Z*L2][a + e*L2] * pharray[P + Z*L2][c + h*L2]

                     + pharray[P + Z*L2][a + h*L2]* pharray[P + Z*L2][c + e*L2]  + pharray[P_ + Z*L2][b + z*L2] * pharray[P_ + Z*L2][d + t*L2]

                     + pharray[P_ + Z*L2][b + t*L2] * pharray[P_ + Z*L2][d + z*L2] ); 

            }

         }

         //(a,d,c,b)_(z,t,h,e) and (b,c,d,a)_(e,h,t,z)
         if(P == Hamiltonian::add(z,t_)){

            for(int Z = 0;Z < 2;++Z){

               (*this)(i,j) += (2*Z + 1.0) * Tools::g6j(0,0,Z,S) * Tools::g6j(0,0,Z,S_) * ( pharray[P + Z*L2][a + z*L2] * pharray[P + Z*L2][c + h*L2]

                     + pharray[P + Z*L2][a + h*L2]* pharray[P + Z*L2][c + z*L2] + pharray[P_ + Z*L2][b + e*L2] * pharray[P_ + Z*L2][d + t*L2]

                     + pharray[P_ + Z*L2][b + t*L2] * pharray[P_ + Z*L2][d + e*L2] ); 

            }

         }

         P = Hamiltonian::add(b,d_);
         P_ = Hamiltonian::gbar(P);

         //(b,d,c,a)_(e,h,t,z) and (a,c,d,b)_(z,t,h,e)
         if(P == Hamiltonian::add(e,h_)){

            for(int Z = 0;Z < 2;++Z){

               (*this)(i,j) += sign * (2*Z + 1.0) * Tools::g6j(0,0,Z,S) * Tools::g6j(0,0,Z,S_) * ( pharray[P + Z*L2][b + e*L2] * pharray[P + Z*L2][c + t*L2]

                     + pharray[P + Z*L2][b + t*L2]* pharray[P + Z*L2][c + e*L2] + pharray[P_ + Z*L2][a + z*L2] * pharray[P_ + Z*L2][d + h*L2]

                     + pharray[P_ + Z*L2][a + h*L2] * pharray[P_ + Z*L2][d + z*L2] ); 

            }

         }

         //(b,d,c,a)_(z,h,t,e) and (a,c,d,b)_(e,t,h,z)
         if(P == Hamiltonian::add(z,h_)){

            for(int Z = 0;Z < 2;++Z){

               (*this)(i,j) += sign * sign_ * (2*Z + 1.0) * Tools::g6j(0,0,Z,S) * Tools::g6j(0,0,Z,S_) * ( pharray[P + Z*L2][b + z*L2] * pharray[P + Z*L2][c + t*L2]

                     + pharray[P + Z*L2][b + t*L2]* pharray[P + Z*L2][c + z*L2] + pharray[P_ + Z*L2][a + e*L2] * pharray[P_ + Z*L2][d + h*L2]

                     + pharray[P_ + Z*L2][a + h*L2] * pharray[P_ + Z*L2][d + e*L2] ); 

            }

         }

         //(b,d,c,a)_(e,t,h,z) and (a,c,d,b)_(z,h,t,e)
         if(P == Hamiltonian::add(e,t_)){

            for(int Z = 0;Z < 2;++Z){

               (*this)(i,j) += sign * sign_ * (2*Z + 1.0) * Tools::g6j(0,0,Z,S) * Tools::g6j(0,0,Z,S_) * ( pharray[P + Z*L2][b + e*L2] * pharray[P + Z*L2][c + h*L2]

                     + pharray[P + Z*L2][b + h*L2]* pharray[P + Z*L2][c + e*L2] + pharray[P_ + Z*L2][a + z*L2] * pharray[P_ + Z*L2][d + t*L2]

                     + pharray[P_ + Z*L2][a + t*L2] * pharray[P_ + Z*L2][d + z*L2] ); 

            }

         }

         //(b,d,c,a)_(z,t,h,e) and (a,c,d,b)_(e,h,t,z)
         if(P == Hamiltonian::add(z,t_)){

            for(int Z = 0;Z < 2;++Z){

               (*this)(i,j) += sign * (2*Z + 1.0) * Tools::g6j(0,0,Z,S) * Tools::g6j(0,0,Z,S_) * ( pharray[P + Z*L2][b + z*L2] * pharray[P + Z*L2][c + h*L2]

                     + pharray[P + Z*L2][b + h*L2]* pharray[P + Z*L2][c + z*L2] + pharray[P_ + Z*L2][a + e*L2] * pharray[P_ + Z*L2][d + t*L2]

                     + pharray[P_ + Z*L2][a + t*L2] * pharray[P_ + Z*L2][d + e*L2] ); 

            }

         }

      }
   }
}

/**
 * construct a TPTPM object by double tracing the direct product of two DPM's, already transformed to a double **array
 */
void TPTPM::dpt2(double **dparray){

   int L2 = Tools::gL2();
   int L4 = L2*L2;
   int L6 = L2*L4;
   int L8 = L2*L6;

   int B,B_;

   int a,b,c,d;
   int e,z,t,h;

   int I_i,J_i,K_i,L_i;

   int S,S_;

   for(int i = 0;i < gn();++i){

      B = tpmm2t[i][0];

      S = TPM::gblock_char(B,0);

      I_i = tpmm2t[i][1];
      J_i = tpmm2t[i][2];

      a = TPM::gt2s(B,I_i,0);
      b = TPM::gt2s(B,I_i,1);
      c = TPM::gt2s(B,J_i,0);
      d = TPM::gt2s(B,J_i,1);

      for(int j = i;j < gn();++j){

         B_ = tpmm2t[j][0];

         S_ = TPM::gblock_char(B_,0);

         K_i = tpmm2t[j][1];
         L_i = tpmm2t[j][2];

         e = TPM::gt2s(B_,K_i,0);
         z = TPM::gt2s(B_,K_i,1);
         t = TPM::gt2s(B_,L_i,0);
         h = TPM::gt2s(B_,L_i,1);

         double ward = 0.0;

         //first S = 1/2 part
         for(int k = 0;k < L2;++k){

            int K_dp = Hamiltonian::add(a,b,k);

            ward += dparray[K_dp][a + b*L2+ e*L4 + z*L6 + S*L8 + 2*S_*L8] * dparray[K_dp][c + d*L2+ t*L4 + h*L6 + S*L8 + 2*S_*L8]

               + dparray[K_dp][a + b*L2+ t*L4 + h*L6 + S*L8 + 2*S_*L8] * dparray[K_dp][c + d*L2+ e*L4 + z*L6 + S*L8 + 2*S_*L8];

         }

         (*this)(i,j) = 2.0 / ( (2*S + 1.0)*(2*S_ + 1.0) ) * ward;

         if(S == 1 && S_ == 1){//only then contribution from S = 3/2 part

            ward = 0.0;

            for(int k = 0;k < L2;++k){

               int K_dp = Hamiltonian::add(a,b,k);

               ward += dparray[K_dp + L2][a + b*L2+ e*L4 + z*L6] * dparray[K_dp + L2][c + d*L2+ t*L4 + h*L6]

                  + dparray[K_dp + L2][a + b*L2+ t*L4 + h*L6] * dparray[K_dp + L2][c + d*L2+ e*L4 + z*L6];

            }

            (*this)(i,j) += 4.0 / 9.0  * ward;

         }

      }
   }

}

/**
 * construct a TPTPM by double tracing the direct product of two PPHM matrices, already translated to 'array' for for faster access
 */
void TPTPM::dpt2_pph(double **ppharray){

   int L2 = Tools::gL2();
   int L4 = L2*L2;
   int L6 = L4*L2;
   int L8 = L6*L2;

   int B,B_;

   int a,b,c,d;
   int e,z,t,h;

   int I_i,J_i,K_i,L_i;

   int S,S_;

   for(int i = 0;i < gn();++i){

      B = tpmm2t[i][0];

      S = TPM::gblock_char(B,0);

      I_i = tpmm2t[i][1];
      J_i = tpmm2t[i][2];

      a = TPM::gt2s(B,I_i,0);
      b = TPM::gt2s(B,I_i,1);
      c = TPM::gt2s(B,J_i,0);
      d = TPM::gt2s(B,J_i,1);

      for(int j = i;j < gn();++j){

         B_ = tpmm2t[j][0];

         S_ = TPM::gblock_char(B_,0);

         K_i = tpmm2t[j][1];
         L_i = tpmm2t[j][2];

         e = TPM::gt2s(B_,K_i,0);
         z = TPM::gt2s(B_,K_i,1);
         t = TPM::gt2s(B_,L_i,0);
         h = TPM::gt2s(B_,L_i,1);

         double ward = 0.0;

         //first S = 1/2 part
         for(int k = 0;k < L2;++k){

            int K_pph = Hamiltonian::add(a,b,k);

            ward += ppharray[K_pph][a + b*L2 + e*L4 + z*L6 + S*L8 + 2*S_*L8] * ppharray[K_pph][c + d*L2 + t*L4 + h*L6 + S*L8 + 2*S_*L8]

               + ppharray[K_pph][a + b*L2 + t*L4 + h*L6 + S*L8 + 2*S_*L8] * ppharray[K_pph][c + d*L2 + e*L4 + z*L6 + S*L8 + 2*S_*L8];

         }

         (*this)(i,j) = 2.0 / ( (2*S + 1.0)*(2*S_ + 1.0) ) * ward;

         if(S == 1 && S_ == 1){//only then contribution from S = 3/2 part

            ward = 0.0;

            for(int k = 0;k < L2;++k){

               int K_pph = Hamiltonian::add(a,b,k);

               ward += ppharray[K_pph + L2][a + b*L2 + e*L4 + z*L6] * ppharray[K_pph + L2][c + d*L2 + t*L4 + h*L6]

                  + ppharray[K_pph + L2][a + b*L2 + t*L4 + h*L6] * ppharray[K_pph + L2][c + d*L2 + e*L4 + z*L6];

            }

            (*this)(i,j) += 4.0 / 9.0  * ward;

         }


      }
   }

}

/**
 * construct a TPTPM by once tracing and once skew-tracing the direct product of two PPHM matrices, already translated to 'array' for for faster access
 */
void TPTPM::dptw(double **ppharray){

   int L2 = Tools::gL2();
   int L4 = L2*L2;
   int L6 = L4*L2;
   int L8 = L6*L2;

   int B,B_;

   int a,b,c,d;
   int e,z,t,h;

   int a_,b_;

   int I_i,J_i,K_i,L_i;

   int S,S_;

   int sign;

   double ward;

   for(int i = 0;i < gn();++i){

      B = tpmm2t[i][0];

      S = TPM::gblock_char(B,0);

      sign = 1 - 2*S;

      I_i = tpmm2t[i][1];
      J_i = tpmm2t[i][2];

      a = TPM::gt2s(B,I_i,0);
      b = TPM::gt2s(B,I_i,1);
      c = TPM::gt2s(B,J_i,0);
      d = TPM::gt2s(B,J_i,1);

      a_ = Hamiltonian::gbar(a);
      b_ = Hamiltonian::gbar(b);

      for(int j = 0;j < gn();++j){

         B_ = tpmm2t[j][0];

         S_ = TPM::gblock_char(B_,0);

         K_i = tpmm2t[j][1];
         L_i = tpmm2t[j][2];

         e = TPM::gt2s(B_,K_i,0);
         z = TPM::gt2s(B_,K_i,1);
         t = TPM::gt2s(B_,L_i,0);
         h = TPM::gt2s(B_,L_i,1);

         (*this)(i,j) = 0.0;

         //S'' = 1/2 first
         for(int J = 0;J < 2;++J){

            ward = 0.0;

            for(int k = 0;k < L2;++k){

               //(da;bc)
               int K_pph = Hamiltonian::add(k,d,a_);

               ward += ppharray[K_pph][e + z*L2 + k*L4 + d*L6 + S_*L8 + 2*J*L8] *  ppharray[K_pph][t + h*L2 + k*L4 + b*L6 + S_*L8 + 2*J*L8]

                  + ppharray[K_pph][e + z*L2 + k*L4 + b*L6 + S_*L8 + 2*J*L8] *  ppharray[K_pph][t + h*L2 + k*L4 + d*L6 + S_*L8 + 2*J*L8];

               K_pph = Hamiltonian::add(k,d,b_);

               ward += sign * ( ppharray[K_pph][e + z*L2 + k*L4 + d*L6 + S_*L8 + 2*J*L8] *  ppharray[K_pph][t + h*L2 + k*L4 + a*L6 + S_*L8 + 2*J*L8]

                  + ppharray[K_pph][e + z*L2 + k*L4 + a*L6 + S_*L8 + 2*J*L8] *  ppharray[K_pph][t + h*L2 + k*L4 + d*L6 + S_*L8 + 2*J*L8] );

               K_pph = Hamiltonian::add(k,c,a_);

               ward += sign * ( ppharray[K_pph][e + z*L2 + k*L4 + c*L6 + S_*L8 + 2*J*L8] *  ppharray[K_pph][t + h*L2 + k*L4 + b*L6 + S_*L8 + 2*J*L8]

                  + ppharray[K_pph][e + z*L2 + k*L4 + b*L6 + S_*L8 + 2*J*L8] *  ppharray[K_pph][t + h*L2 + k*L4 + c*L6 + S_*L8 + 2*J*L8] );

               K_pph = Hamiltonian::add(k,c,b_);

               ward += ppharray[K_pph][e + z*L2 + k*L4 + c*L6 + S_*L8 + 2*J*L8] *  ppharray[K_pph][t + h*L2 + k*L4 + a*L6 + S_*L8 + 2*J*L8]

                  + ppharray[K_pph][e + z*L2 + k*L4 + a*L6 + S_*L8 + 2*J*L8] *  ppharray[K_pph][t + h*L2 + k*L4 + c*L6 + S_*L8 + 2*J*L8];

            }

            (*this)(i,j) += ward * (2*J + 1.0) * Tools::g6j(0,0,S,J);

         }

         (*this)(i,j) *= 2.0 / (2*S_ + 1.0);

         //then S'' = 3/2, only
         if(S_ == 1){

            ward = 0.0;

            for(int k = 0;k < L2;++k){

               //(da;bc)
               int K_pph = Hamiltonian::add(k,d,a_);

               ward += ppharray[K_pph + L2][e + z*L2 + k*L4 + d*L6] *  ppharray[K_pph + L2][t + h*L2 + k*L4 + b*L6]

                  + ppharray[K_pph + L2][e + z*L2 + k*L4 + b*L6] *  ppharray[K_pph + L2][t + h*L2 + k*L4 + d*L6];

               K_pph = Hamiltonian::add(k,d,b_);

               ward += sign * ( ppharray[K_pph + L2][e + z*L2 + k*L4 + d*L6] *  ppharray[K_pph + L2][t + h*L2 + k*L4 + a*L6]

                  + ppharray[K_pph + L2][e + z*L2 + k*L4 + a*L6] *  ppharray[K_pph + L2][t + h*L2 + k*L4 + d*L6] );

               K_pph = Hamiltonian::add(k,c,a_);

               ward += sign * ( ppharray[K_pph + L2][e + z*L2 + k*L4 + c*L6] *  ppharray[K_pph + L2][t + h*L2 + k*L4 + b*L6]

                  + ppharray[K_pph + L2][e + z*L2 + k*L4 + b*L6] *  ppharray[K_pph + L2][t + h*L2 + k*L4 + c*L6] );

               K_pph = Hamiltonian::add(k,c,b_);

               ward += ppharray[K_pph + L2][e + z*L2 + k*L4 + c*L6] *  ppharray[K_pph + L2][t + h*L2 + k*L4 + a*L6]

                  + ppharray[K_pph + L2][e + z*L2 + k*L4 + a*L6] *  ppharray[K_pph + L2][t + h*L2 + k*L4 + c*L6];

            }

            (*this)(i,j) +=  4.0 * ward * Tools::g6j(0,0,S,1);

         }

      }
   }

}


/**
 * construct a TPTPM by twice skew-tracing the direct product of two PPHM matrices, already translated to 'array' for for faster access
 */
void TPTPM::dpw2(double **ppharray){

   int L2 = Tools::gL2();
   int L4 = L2*L2;
   int L6 = L4*L2;
   int L8 = L6*L2;

   int B,B_;

   int a,b,c,d;
   int e,z,t,h;

   int a_,b_;
   int e_,z_;

   int I_i,J_i,K_i,L_i;

   int S,S_;

   int sign,sign_;

   double ward;

   int K_pph,m;

   for(int i = 0;i < gn();++i){

      B = tpmm2t[i][0];

      S = TPM::gblock_char(B,0);

      sign = 1 - 2*S;

      I_i = tpmm2t[i][1];
      J_i = tpmm2t[i][2];

      a = TPM::gt2s(B,I_i,0);
      b = TPM::gt2s(B,I_i,1);
      c = TPM::gt2s(B,J_i,0);
      d = TPM::gt2s(B,J_i,1);

      a_ = Hamiltonian::gbar(a);
      b_ = Hamiltonian::gbar(b);

      for(int j = i;j < gn();++j){

         B_ = tpmm2t[j][0];

         S_ = TPM::gblock_char(B_,0);

         sign_ = 1 - 2*S_;

         K_i = tpmm2t[j][1];
         L_i = tpmm2t[j][2];

         e = TPM::gt2s(B_,K_i,0);
         z = TPM::gt2s(B_,K_i,1);
         t = TPM::gt2s(B_,L_i,0);
         h = TPM::gt2s(B_,L_i,1);

         e_ = Hamiltonian::gbar(e);
         z_ = Hamiltonian::gbar(z);

         (*this)(i,j) = 0.0;

         //first S'' = 1/2
         for(int J = 0;J < 2;++J)
            for(int J_ = 0;J_ < 2;++J_){

               ward = 0.0;

               for(int k = 0;k < L2;++k){

                  K_pph = Hamiltonian::add(k,d,a_);

                  //1) da;bc he;zt
                  m = Hamiltonian::adjoint(K_pph,h,e_);

                  ward += ppharray[K_pph][k + d*L2 + m*L4 + h*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + b*L2 + m*L4 + z*L6 + J*L8 + 2*J_*L8]

                     + ppharray[K_pph][k + d*L2 + m*L4 + z*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + b*L2 + m*L4 + h*L6 + J*L8 + 2*J_*L8];

                  //2) da;bc hz;et
                  m = Hamiltonian::adjoint(K_pph,h,z_);

                  ward += sign_ * ( ppharray[K_pph][k + d*L2 + m*L4 + h*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + b*L2 + m*L4 + e*L6 + J*L8 + 2*J_*L8]

                     + ppharray[K_pph][k + d*L2 + m*L4 + e*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + b*L2 + m*L4 + h*L6 + J*L8 + 2*J_*L8] );

                  //3) da;bc te;zh
                  m = Hamiltonian::adjoint(K_pph,t,e_);

                  ward += sign_ * (ppharray[K_pph][k + d*L2 + m*L4 + t*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + b*L2 + m*L4 + z*L6 + J*L8 + 2*J_*L8]

                     + ppharray[K_pph][k + d*L2 + m*L4 + z*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + b*L2 + m*L4 + t*L6 + J*L8 + 2*J_*L8] );

                  //4) da;bc tz;eh
                  m = Hamiltonian::adjoint(K_pph,t,z_);

                  ward += ppharray[K_pph][k + d*L2 + m*L4 + t*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + b*L2 + m*L4 + e*L6 + J*L8 + 2*J_*L8]

                     + ppharray[K_pph][k + d*L2 + m*L4 + e*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + b*L2 + m*L4 + t*L6 + J*L8 + 2*J_*L8];

                  K_pph = Hamiltonian::add(k,d,b_);

                  //5) db;ac he;zt
                  m = Hamiltonian::adjoint(K_pph,h,e_);

                  ward += sign * ( ppharray[K_pph][k + d*L2 + m*L4 + h*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + a*L2 + m*L4 + z*L6 + J*L8 + 2*J_*L8]

                     + ppharray[K_pph][k + d*L2 + m*L4 + z*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + a*L2 + m*L4 + h*L6 + J*L8 + 2*J_*L8] );

                  //6) db;ac hz;et
                  m = Hamiltonian::adjoint(K_pph,h,z_);

                  ward += sign * sign_ * ( ppharray[K_pph][k + d*L2 + m*L4 + h*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + a*L2 + m*L4 + e*L6 + J*L8 + 2*J_*L8]

                     + ppharray[K_pph][k + d*L2 + m*L4 + e*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + a*L2 + m*L4 + h*L6 + J*L8 + 2*J_*L8] );

                  //7) db;ac te;zh
                  m = Hamiltonian::adjoint(K_pph,t,e_);

                  ward += sign * sign_ * ( ppharray[K_pph][k + d*L2 + m*L4 + t*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + a*L2 + m*L4 + z*L6 + J*L8 + 2*J_*L8]

                     + ppharray[K_pph][k + d*L2 + m*L4 + z*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + a*L2 + m*L4 + t*L6 + J*L8 + 2*J_*L8] );

                  //8) db;ac tz;eh
                  m = Hamiltonian::adjoint(K_pph,t,z_);

                  ward += sign * ( ppharray[K_pph][k + d*L2 + m*L4 + t*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + a*L2 + m*L4 + e*L6 + J*L8 + 2*J_*L8]

                     + ppharray[K_pph][k + d*L2 + m*L4 + e*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + a*L2 + m*L4 + t*L6 + J*L8 + 2*J_*L8] );

                  K_pph = Hamiltonian::add(k,c,a_);

                  //9) ca;bd he;zt
                  m = Hamiltonian::adjoint(K_pph,h,e_);

                  ward += sign * ( ppharray[K_pph][k + c*L2 + m*L4 + h*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + b*L2 + m*L4 + z*L6 + J*L8 + 2*J_*L8]

                     + ppharray[K_pph][k + c*L2 + m*L4 + z*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + b*L2 + m*L4 + h*L6 + J*L8 + 2*J_*L8] );

                  //10) ca;bd hz;et
                  m = Hamiltonian::adjoint(K_pph,h,z_);

                  ward += sign * sign_ * ( ppharray[K_pph][k + c*L2 + m*L4 + h*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + b*L2 + m*L4 + e*L6 + J*L8 + 2*J_*L8]

                     + ppharray[K_pph][k + c*L2 + m*L4 + e*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + b*L2 + m*L4 + h*L6 + J*L8 + 2*J_*L8] );

                  //11) ca;bd te;zh
                  m = Hamiltonian::adjoint(K_pph,t,e_);

                  ward += sign * sign_ * ( ppharray[K_pph][k + c*L2 + m*L4 + t*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + b*L2 + m*L4 + z*L6 + J*L8 + 2*J_*L8]

                     + ppharray[K_pph][k + c*L2 + m*L4 + z*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + b*L2 + m*L4 + t*L6 + J*L8 + 2*J_*L8] );

                  //12) ca;bd tz;eh
                  m = Hamiltonian::adjoint(K_pph,t,z_);

                  ward += sign * ( ppharray[K_pph][k + c*L2 + m*L4 + t*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + b*L2 + m*L4 + e*L6 + J*L8 + 2*J_*L8]

                     + ppharray[K_pph][k + c*L2 + m*L4 + e*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + b*L2 + m*L4 + t*L6 + J*L8 + 2*J_*L8] );

                  K_pph = Hamiltonian::add(k,c,b_);

                  //13) cb;ad he;zt
                  m = Hamiltonian::adjoint(K_pph,h,e_);

                  ward += ppharray[K_pph][k + c*L2 + m*L4 + h*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + a*L2 + m*L4 + z*L6 + J*L8 + 2*J_*L8]

                     + ppharray[K_pph][k + c*L2 + m*L4 + z*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + a*L2 + m*L4 + h*L6 + J*L8 + 2*J_*L8];

                  //14) cb;ad hz;et
                  m = Hamiltonian::adjoint(K_pph,h,z_);

                  ward += sign_ * ( ppharray[K_pph][k + c*L2 + m*L4 + h*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + a*L2 + m*L4 + e*L6 + J*L8 + 2*J_*L8]

                     + ppharray[K_pph][k + c*L2 + m*L4 + e*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + a*L2 + m*L4 + h*L6 + J*L8 + 2*J_*L8] );

                  //15) cb;ad te;zh
                  m = Hamiltonian::adjoint(K_pph,t,e_);

                  ward += sign_ * ( ppharray[K_pph][k + c*L2 + m*L4 + t*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + a*L2 + m*L4 + z*L6 + J*L8 + 2*J_*L8]

                     + ppharray[K_pph][k + c*L2 + m*L4 + z*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + a*L2 + m*L4 + t*L6 + J*L8 + 2*J_*L8] );

                  //16) cb;ad tz;eh
                  m = Hamiltonian::adjoint(K_pph,t,z_);

                  ward += ppharray[K_pph][k + c*L2 + m*L4 + t*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + a*L2 + m*L4 + e*L6 + J*L8 + 2*J_*L8]

                     + ppharray[K_pph][k + c*L2 + m*L4 + e*L6 + J*L8 + 2*J_*L8] * ppharray[K_pph][k + a*L2 + m*L4 + t*L6 + J*L8 + 2*J_*L8];

               }//end of loop over k
               
               (*this)(i,j) += (2*J + 1.0) * (2*J_ + 1.0) * Tools::g6j(0,0,S,J) * Tools::g6j(0,0,S_,J_) * ward;

            }//end of loops over J,J_

         (*this)(i,j) *= 2.0;

         //S = 3/2
         ward = 0.0;

         for(int k = 0;k < L2;++k){

            K_pph = Hamiltonian::add(k,d,a_);

            //1) da;bc he;zt
            m = Hamiltonian::adjoint(K_pph,h,e_);

            ward += ppharray[K_pph + L2][k + d*L2 + m*L4 + h*L6] * ppharray[K_pph + L2][k + b*L2 + m*L4 + z*L6]

               + ppharray[K_pph + L2][k + d*L2 + m*L4 + z*L6] * ppharray[K_pph + L2][k + b*L2 + m*L4 + h*L6];

            //2) da;bc hz;et
            m = Hamiltonian::adjoint(K_pph,h,z_);

            ward += sign_ * ( ppharray[K_pph + L2][k + d*L2 + m*L4 + h*L6] * ppharray[K_pph + L2][k + b*L2 + m*L4 + e*L6]

               + ppharray[K_pph + L2][k + d*L2 + m*L4 + e*L6] * ppharray[K_pph + L2][k + b*L2 + m*L4 + h*L6] );

            //3) da;bc te;zh
            m = Hamiltonian::adjoint(K_pph,t,e_);

            ward += sign_ * ( ppharray[K_pph + L2][k + d*L2 + m*L4 + t*L6] * ppharray[K_pph + L2][k + b*L2 + m*L4 + z*L6]

               + ppharray[K_pph + L2][k + d*L2 + m*L4 + z*L6] * ppharray[K_pph + L2][k + b*L2 + m*L4 + t*L6] );

            //4) da;bc tz;eh
            m = Hamiltonian::adjoint(K_pph,t,z_);

            ward += ppharray[K_pph + L2][k + d*L2 + m*L4 + t*L6] * ppharray[K_pph + L2][k + b*L2 + m*L4 + e*L6]

               + ppharray[K_pph + L2][k + d*L2 + m*L4 + e*L6] * ppharray[K_pph + L2][k + b*L2 + m*L4 + t*L6];

            K_pph = Hamiltonian::add(k,d,b_);

            //5) db;ac he;zt
            m = Hamiltonian::adjoint(K_pph,h,e_);

            ward += sign * ( ppharray[K_pph + L2][k + d*L2 + m*L4 + h*L6] * ppharray[K_pph + L2][k + a*L2 + m*L4 + z*L6]

               + ppharray[K_pph + L2][k + d*L2 + m*L4 + z*L6] * ppharray[K_pph + L2][k + a*L2 + m*L4 + h*L6] );

            //6) db;ac hz;et
            m = Hamiltonian::adjoint(K_pph,h,z_);

            ward += sign * sign_ * ( ppharray[K_pph + L2][k + d*L2 + m*L4 + h*L6] * ppharray[K_pph + L2][k + a*L2 + m*L4 + e*L6]

               + ppharray[K_pph + L2][k + d*L2 + m*L4 + e*L6] * ppharray[K_pph + L2][k + a*L2 + m*L4 + h*L6] );

            //7) db;ac te;zh
            m = Hamiltonian::adjoint(K_pph,t,e_);

            ward += sign * sign_ * ( ppharray[K_pph + L2][k + d*L2 + m*L4 + t*L6] * ppharray[K_pph + L2][k + a*L2 + m*L4 + z*L6]

               + ppharray[K_pph + L2][k + d*L2 + m*L4 + z*L6] * ppharray[K_pph + L2][k + a*L2 + m*L4 + t*L6] );

            //8) db;ac tz;eh
            m = Hamiltonian::adjoint(K_pph,t,z_);

            ward += sign * ( ppharray[K_pph + L2][k + d*L2 + m*L4 + t*L6] * ppharray[K_pph + L2][k + a*L2 + m*L4 + e*L6]

               + ppharray[K_pph + L2][k + d*L2 + m*L4 + e*L6] * ppharray[K_pph + L2][k + a*L2 + m*L4 + t*L6] );

            K_pph = Hamiltonian::add(k,c,a_);

            //9) ca;bd he;zt
            m = Hamiltonian::adjoint(K_pph,h,e_);

            ward += sign * ( ppharray[K_pph + L2][k + c*L2 + m*L4 + h*L6] * ppharray[K_pph + L2][k + b*L2 + m*L4 + z*L6]

               + ppharray[K_pph + L2][k + c*L2 + m*L4 + z*L6] * ppharray[K_pph + L2][k + b*L2 + m*L4 + h*L6] );

            //10) ca;bd hz;et
            m = Hamiltonian::adjoint(K_pph,h,z_);

            ward += sign * sign_ * ( ppharray[K_pph + L2][k + c*L2 + m*L4 + h*L6] * ppharray[K_pph + L2][k + b*L2 + m*L4 + e*L6]

               + ppharray[K_pph + L2][k + c*L2 + m*L4 + e*L6] * ppharray[K_pph + L2][k + b*L2 + m*L4 + h*L6] );

            //11) ca;bd te;zh
            m = Hamiltonian::adjoint(K_pph,t,e_);

            ward += sign * sign_ * ( ppharray[K_pph + L2][k + c*L2 + m*L4 + t*L6] * ppharray[K_pph + L2][k + b*L2 + m*L4 + z*L6]

               + ppharray[K_pph + L2][k + c*L2 + m*L4 + z*L6] * ppharray[K_pph + L2][k + b*L2 + m*L4 + t*L6] );

            //12) ca;bd tz;eh
            m = Hamiltonian::adjoint(K_pph,t,z_);

            ward += sign * ( ppharray[K_pph + L2][k + c*L2 + m*L4 + t*L6] * ppharray[K_pph + L2][k + b*L2 + m*L4 + e*L6]

               + ppharray[K_pph + L2][k + c*L2 + m*L4 + e*L6] * ppharray[K_pph + L2][k + b*L2 + m*L4 + t*L6] );

            K_pph = Hamiltonian::add(k,c,b_);

            //13) cb;ad he;zt
            m = Hamiltonian::adjoint(K_pph,h,e_);

            ward += ppharray[K_pph + L2][k + c*L2 + m*L4 + h*L6] * ppharray[K_pph + L2][k + a*L2 + m*L4 + z*L6]

               + ppharray[K_pph + L2][k + c*L2 + m*L4 + z*L6] * ppharray[K_pph + L2][k + a*L2 + m*L4 + h*L6];

            //14) cb;ad hz;et
            m = Hamiltonian::adjoint(K_pph,h,z_);

            ward += sign_ * ( ppharray[K_pph + L2][k + c*L2 + m*L4 + h*L6] * ppharray[K_pph + L2][k + a*L2 + m*L4 + e*L6]

               + ppharray[K_pph + L2][k + c*L2 + m*L4 + e*L6] * ppharray[K_pph + L2][k + a*L2 + m*L4 + h*L6] );

            //15) cb;ad te;zh
            m = Hamiltonian::adjoint(K_pph,t,e_);

            ward += sign_ * ( ppharray[K_pph + L2][k + c*L2 + m*L4 + t*L6] * ppharray[K_pph + L2][k + a*L2 + m*L4 + z*L6]

               + ppharray[K_pph + L2][k + c*L2 + m*L4 + z*L6] * ppharray[K_pph + L2][k + a*L2 + m*L4 + t*L6] );

            //16) cb;ad tz;eh
            m = Hamiltonian::adjoint(K_pph,t,z_);

            ward += ppharray[K_pph + L2][k + c*L2 + m*L4 + t*L6] * ppharray[K_pph + L2][k + a*L2 + m*L4 + e*L6]

               + ppharray[K_pph + L2][k + c*L2 + m*L4 + e*L6] * ppharray[K_pph + L2][k + a*L2 + m*L4 + t*L6];

         }//end of loop over k

         (*this)(i,j) += 4.0 * 9.0 * Tools::g6j(0,0,S,1) * Tools::g6j(0,0,S_,1) * ward;

      }
   }

}
