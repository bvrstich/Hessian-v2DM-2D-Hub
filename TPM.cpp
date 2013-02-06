#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::cout;
using std::endl;
using std::ios;

#include "include.h"

vector< vector<int> > *TPM::t2s;
int ***TPM::s2t;

int **TPM::block_char;

double **TPM::norm;

double TPM::Sa;
double TPM::Sb;
double TPM::Sc;

/**
 * static function that initializes the static variables
 */
void TPM::init(){

   //allocate stuff
   t2s = new vector< vector<int> > [Tools::gM()];

   s2t = new int ** [Tools::gM()];

   for(int B = 0;B < Tools::gM();++B){

      s2t[B] = new int * [Tools::gL2()];

      for(int a = 0;a < Tools::gL2();++a)
         s2t[B][a] = new int [Tools::gL2()];

   }

   block_char = new int * [Tools::gM()];

   for(int B = 0;B < Tools::gM();++B)
      block_char[B] = new int [2];

   norm = new double * [Tools::gL2()];

   for(int a = 0;a < Tools::gL2();++a)
      norm[a] = new double [Tools::gL2()];

   vector<int> v(2);

   int block = 0;

   //tp index
   int t;

   //loop over the K_x K_y blocks
   for(int K = 0;K < Tools::gL2();++K){

      t = 0;

      //S = 0
      block_char[block][0] = 0;
      block_char[block][1] = K;

      for(int a = 0;a < Tools::gL2();++a)
         for(int b = a;b < Tools::gL2();++b){

            if( ( (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0)) % Tools::gL() == Hamiltonian::ga_xy(K,0) ) 

                  && ( (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1)) % Tools::gL() == Hamiltonian::ga_xy(K,1) ) ){

               v[0] = a;
               v[1] = b;

               t2s[block].push_back(v);

               s2t[block][a][b] = t;
               s2t[block][b][a] = t;

               if(a == b)
                  norm[a][b] = 1.0/(std::sqrt(2.0));
               else
                  norm[a][b] = 1.0;

               norm[b][a] = norm[a][b];

               ++t;

            }

         }

      t = 0;

      //S = 1
      block_char[Tools::gL2() + block][0] = 1;
      block_char[Tools::gL2() + block][1] = K;

      for(int a = 0;a < Tools::gL2();++a)
         for(int b = a + 1;b < Tools::gL2();++b){

            if( ( (Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0)) % Tools::gL() == Hamiltonian::ga_xy(K,0) ) 

                  && ( (Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1)) % Tools::gL() == Hamiltonian::ga_xy(K,1) ) ){

               v[0] = a;
               v[1] = b;

               t2s[Tools::gL2() + block].push_back(v);

               s2t[Tools::gL2() + block][a][b] = t;
               s2t[Tools::gL2() + block][b][a] = t;

               ++t;

            }

         }

      ++block;

   }

  //only for overlapmatrix coefficients!
   int M = Tools::gM();
   int N = Tools::gN();

   Sa = 1.0;
   Sb = 0.0;
   Sc = 0.0;

#ifdef __Q_CON
   Sa += 1.0;
   Sb += (4.0*N*N + 2.0*N - 4.0*N*M + M*M - M)/(N*N*(N - 1.0)*(N - 1.0));
   Sc += (2.0*N - M)/((N - 1.0)*(N - 1.0));
#endif

#ifdef __G_CON
   Sa += 4.0;
   Sc += (2.0*N - M - 2.0)/((N - 1.0)*(N - 1.0));
#endif

#ifdef __T1_CON
   Sa += M - 4.0;
   Sb += (M*M*M - 6.0*M*M*N -3.0*M*M + 12.0*M*N*N + 12.0*M*N + 2.0*M - 18.0*N*N - 6.0*N*N*N)/( 3.0*N*N*(N - 1.0)*(N - 1.0) );
   Sc -= (M*M + 2.0*N*N - 4.0*M*N - M + 8.0*N - 4.0)/( 2.0*(N - 1.0)*(N - 1.0) );
#endif

#ifdef __T2_CON
   Sa += 5.0*M - 8.0;
   Sb += 2.0/(N - 1.0);
   Sc += (2.0*N*N + (M - 2.0)*(4.0*N - 3.0) - M*M)/(2.0*(N - 1.0)*(N - 1.0));
#endif

}

/**
 * static function that deallocates the static variables
 */
void TPM::clear(){

   delete [] t2s;

   for(int B = 0;B < Tools::gM();++B){

      for(int a = 0;a < Tools::gL2();++a)
         delete [] s2t[B][a];

      delete [] s2t[B];

      delete [] block_char[B];

   }

   delete [] s2t;

   delete [] block_char;

   for(int a = 0;a < Tools::gL2();++a)
      delete [] norm[a];

   delete [] norm;

}

/**
 * standard constructor for a spinsymmetrical, translationally invariant tp matrix on a 2D lattice: 
 * constructs BlockMatrix object with 2 * L^2 blocks, L^2 for S = 0 and L^2 for S = 1,
 */
TPM::TPM() : BlockMatrix(Tools::gM()) {

   //set the dimension of the blocks
   for(int B = 0;B < Tools::gL2();++B)//S = 0
      setMatrixDim(B,t2s[B].size(),1);

   for(int B = Tools::gL2();B < Tools::gM();++B)//S = 1
      setMatrixDim(B,t2s[B].size(),3);

}

/**
 * copy constructor for a spinsymmetrical, translationally invariant tp matrix on a 2D lattice:
 * constructs BlockMatrix object with 2 * L^2 blocks, L^2 for S = 0 and L^2 for S = 1,
 * @param tpm_c The TPM object to be copied into (*this)
 */
TPM::TPM(const TPM &tpm_c) : BlockMatrix(tpm_c){ }

/**
 * destructor
 */
TPM::~TPM(){ }

ostream &operator<<(ostream &output,const TPM &tpm_p){

   int S,K;

   for(int B = 0;B < tpm_p.gnr();++B){

      S = tpm_p.block_char[B][0];
      K = tpm_p.block_char[B][1];

      output << "S =\t" << S << "\tK_x =\t" << Hamiltonian::ga_xy(K,0) << "\tK_y =\t" << Hamiltonian::ga_xy(K,1) <<

         "\tdimension =\t" << tpm_p.gdim(B) << "\tdegeneracy =\t" << tpm_p.gdeg(B) << std::endl;

      output << std::endl;

      for(int i = 0;i < tpm_p.gdim(B);++i)
         for(int j = 0;j < tpm_p.gdim(B);++j){

            output << i << "\t" << j << "\t|\t" << tpm_p.t2s[B][i][0] << "\t" << tpm_p.t2s[B][i][1]

               << "\t" << tpm_p.t2s[B][j][0] << "\t" << tpm_p.t2s[B][j][1] << "\t" << tpm_p(B,i,j) << endl;

         }

      output << std::endl;

   }

   return output;

}

/**
 * access the elements of the the blocks in sp mode, the symmetry or antisymmetry of the blocks is automatically accounted for:\n\n
 * Antisymmetrical for S = 1, symmetrical in the sp orbitals for S = 0\n\n
 * @param S The tp spin quantumnumber
 * @param a first sp index that forms the tp row index i of block B(S,K_x,K_y), together with b
 * @param b second sp index that forms the tp row index i of block B(S,K_x,K_y), together with a
 * @param c first sp index that forms the tp column index j of block B(S,K_x,K_y), together with d
 * @param d second sp index that forms the tp column index j of block B(S,K_x,K_y), together with c
 * @return the number on place TPM(B,i,j) with the right phase.
 */
double TPM::operator()(int S,int a,int b,int c,int d) const{

   int K_x = ( Hamiltonian::ga_xy(a,0) + Hamiltonian::ga_xy(b,0) )%Tools::gL();
   int K_y = ( Hamiltonian::ga_xy(a,1) + Hamiltonian::ga_xy(b,1) )%Tools::gL();

   //check if momentum is ok
   if( K_x != ( Hamiltonian::ga_xy(c,0) + Hamiltonian::ga_xy(d,0) )%Tools::gL() )
      return 0;

   if( K_y !=  ( Hamiltonian::ga_xy(c,1) + Hamiltonian::ga_xy(d,1) )%Tools::gL() )
      return 0;

   int K = Hamiltonian::gxy_a(K_x,K_y);
   int B = K + S*Tools::gL2();

   if(S == 0){

      int i = s2t[B][a][b];
      int j = s2t[B][c][d];

      return (*this)(B,i,j);

   }
   else{//S = 1

      if( (a == b) || (c == d) )
         return 0;
      else{

         int i = s2t[B][a][b];
         int j = s2t[B][c][d];

         int phase = 1;

         if(a > b)
            phase *= -1;
         if(c > d)
            phase *= -1;

         return phase * (*this)(B,i,j);

      }

   }

}

/**
 * construct the spinsymmetrical hubbard hamiltonian in momentum space with on site repulsion U
 * @param U onsite repulsion term
 */
void TPM::hubbard(double U){

   int a,b,c,d;//sp indices

   double ward = 1.0/(Tools::gN() - 1.0);

   for(int B = 0;B < gnr();++B){

      for(int i = 0;i < gdim(B);++i){

         a = t2s[B][i][0];
         b = t2s[B][i][1];

         for(int j = i;j < gdim(B);++j){

            c = t2s[B][j][0];
            d = t2s[B][j][1];

            //init
            (*this)(B,i,j) = 0;

            //hopping (kinetic energy):
            if(i == j){

               (*this)(B,i,i) = -2.0 * ward * ( cos( 2.0 * Hamiltonian::ga_xy(a,0) * 3.141592653589793238462 / (double) Tools::gL())  

                     + cos( 2.0 * Hamiltonian::ga_xy(a,1) * 3.141592653589793238462 / (double) Tools::gL()) 

                     + cos( 2.0 * Hamiltonian::ga_xy(b,0) * 3.141592653589793238462 / (double) Tools::gL()) 

                     + cos( 2.0 * Hamiltonian::ga_xy(b,1) * 3.141592653589793238462 / (double) Tools::gL()) );

            }

            //on-site repulsion
            if(B < Tools::gL()*Tools::gL()){//only spin zero

               double hard = 2.0*U / (double) (Tools::gL()*Tools::gL());

               if(a == b)
                  hard /= std::sqrt(2.0);

               if(c == d)
                  hard /= std::sqrt(2.0);

               (*this)(B,i,j) += hard;

            }

         }
      }

   }

   this->symmetrize();

}

/**
 * initialize this onto the unitmatrix with trace N*(N - 1)/2
 */
void TPM::unit(){

   double ward = Tools::gN()*(Tools::gN() - 1.0)/(Tools::gM()*(Tools::gM() - 1.0));

   for(int B = 0;B < gnr();++B){

      for(int i = 0;i < gdim(B);++i){

         (*this)(B,i,i) = ward;

         for(int j = i + 1;j < gdim(B);++j)
            (*this)(B,i,j) = (*this)(B,j,i) = 0.0;

      }
   }

}

/**
 * @return The expectation value of the total spin for the TPM.
 */
double TPM::S_2() const{

   double ward = 0.0;

   int S;

   for(int B = 0;B < gnr();++B){

      S = block_char[B][0];

      if(S == 0){

         for(int i = 0;i < gdim(B);++i)
            ward += -1.5 * (Tools::gN() - 2.0)/(Tools::gN() - 1.0) * (*this)(B,i,i);

      }
      else{

         for(int i = 0;i < this->gdim(B);++i)
            ward += 3.0 * ( -1.5 * (Tools::gN() - 2.0)/(Tools::gN() - 1.0) + 2.0 ) * (*this)(B,i,i);

      }

   }

   return ward;

}

/**
 * access to the lists from outside of the class
 */
int TPM::gt2s(int B,int i,int option){

   return t2s[B][i][option];

}

/**
 * access to the lists from outside of the class
 */
int TPM::gs2t(int B,int a,int b){

   return s2t[B][a][b];

}

/**
 * access to the lists from outside of the class
 */
int TPM::gblock_char(int B,int option){

   return block_char[B][option];

}

/**
 * access to the lists from outside of the class
 */
int TPM::gdim(int B){

   return t2s[B].size();

}

/**
 * convert a 2DM object from vector to matrix/TPM form
 * @param grad input Gradient object
 */
void TPM::convert(const Gradient &grad){

   int tpmm;
   int S;

   for(int B = 0;B < Tools::gM();++B){

      S = block_char[B][0];

      for(int i = 0;i < gdim(B);++i)
         for(int j = i;j < gdim(B);++j){

            tpmm = TPTPM::gt2tpmm(B,i,j);

            (*this)(B,i,j) = grad[tpmm]/ ( 2.0 * Gradient::gnorm(tpmm) * (2.0 * S + 1.0) );

         }

   }

   this->symmetrize();

}

/**
 * The spincoupled Q map
 * @param option = 1, regular Q map , = -1 inverse Q map
 * @param tpm_d the TPM of which the Q map is taken and saved in this.
 */
void TPM::Q(int option,const TPM &tpm_d){

   double a = 1;
   double b = 1.0/(Tools::gN()*(Tools::gN() - 1.0));
   double c = 1.0/(Tools::gN() - 1.0);

   this->Q(option,a,b,c,tpm_d);

}

/**
 * The spincoupled Q-like map: see primal-dual.pdf for more info (form: Q^S(A,B,C)(TPM) )
 * @param option = 1, regular Q-like map , = -1 inverse Q-like map
 * @param A factor in front of the two particle piece of the map
 * @param B factor in front of the no particle piece of the map
 * @param C factor in front of the single particle piece of the map
 * @param tpm_d the TPM of which the Q-like map is taken and saved in this.
 */
void TPM::Q(int option,double A,double B,double C,const TPM &tpm_d){

   //for inverse
   if(option == -1){

      B = (B*A + B*C*Tools::gM() - 2.0*C*C)/( A * (C*(Tools::gM() - 2.0) -  A) * ( A + B*Tools::gM()*(Tools::gM() - 1.0) - 2.0*C*(Tools::gM() - 1.0) ) );
      C = C/(A*(C*(Tools::gM() - 2.0) - A));
      A = 1.0/A;

   }

   SPM spm;
   spm.bar(C,tpm_d);

   //de trace*2 omdat mijn definitie van trace in berekeningen over alle (alpha,beta) loopt
   double ward = B*tpm_d.trace()*2.0;

   int a,b;

   for(int B = 0;B < gnr();++B){

      for(int i = 0;i < gdim(B);++i){

         a = t2s[B][i][0];
         b = t2s[B][i][1];

         //tp part is only nondiagonal part
         for(int j = i;j < gdim(B);++j)
            (*this)(B,i,j) = A * tpm_d(B,i,j);

         (*this)(B,i,i) += ward - spm[a] - spm[b];

      }

   }

   this->symmetrize();

}

/**
 * access to the TPM norms from outside the class
 */
double TPM::gnorm(int a,int b){

   return norm[a][b];

}

/**
 * The G down map, maps a PHM object onto a TPM object using the G map.
 * @param phm input PHM
 */
void TPM::G(const PHM &phm){

   SPM spm;
   spm.bar(1.0/(Tools::gN() - 1.0),phm);

   int a,b,c,d;

   //the conjugated indices
   int a_,b_,c_,d_;

   int S;

   int sign;

   for(int B = 0;B < gnr();++B){

      S = block_char[B][0];

      sign = 1 - 2*S;

      for(int i = 0;i < gdim(B);++i){

         a = t2s[B][i][0];
         b = t2s[B][i][1];

         a_ = Hamiltonian::gbar(a);
         b_ = Hamiltonian::gbar(b);

         //tp part is only nondiagonal part
         for(int j = i;j < gdim(B);++j){

            c = t2s[B][j][0];
            d = t2s[B][j][1];

            c_ = Hamiltonian::gbar(c);
            d_ = Hamiltonian::gbar(d);

            (*this)(B,i,j) = 0.0;

            //four ph terms:
            for(int Z = 0;Z < 2;++Z)
               (*this)(B,i,j) -= (2.0*Z + 1.0) * Tools::g6j(0,0,S,Z) * ( phm(Z,a,d_,c,b_) + phm(Z,b,c_ ,d, a_) 

                     + sign * ( phm(Z,b,d_ ,c, a_ ) + phm(Z,a,c_,d,b_) ) );

            //norm:
            if(a == b)
               (*this)(B,i,j) /= std::sqrt(2.0);

            if(c == d)
               (*this)(B,i,j) /= std::sqrt(2.0);

         }

         (*this)(B,i,i) += spm[a] + spm[b];

      }

   }

   this->symmetrize();

}

/** 
 * The T1-down map that maps a DPM on TPM. This is just a Q-like map using the TPM::bar (dpm) as input.
 * @param dpm the input DPM matrix
 */
void TPM::T(const DPM &dpm){

   TPM tpm;
   tpm.bar(1.0,dpm);

   double a = 1;
   double b = 1.0/(3.0*Tools::gN()*(Tools::gN() - 1.0));
   double c = 0.5/(Tools::gN() - 1.0);

   this->Q(1,a,b,c,tpm);

}

/**
 * Construct a spincoupled, translationally invariant TPM matrix out of a spincoupled, translationally invariant DPM matrix.
 * For the definition and derivation see symmetry.pdf
 * @param scale scalefactor for the traced DPM
 * @param dpm input DPM
 */
void TPM::bar(double scale,const DPM &dpm){

   int a,b,c,d;

   double ward;

   //first the S = 0 part, easiest:
   for(int B = 0;B < Tools::gL2();++B){

      for(int i = 0;i < gdim(B);++i){

         a = t2s[B][i][0];
         b = t2s[B][i][1];

         for(int j = i;j < gdim(B);++j){

            c = t2s[B][j][0];
            d = t2s[B][j][1];

            (*this)(B,i,j) = 0.0;

            //only total S = 1/2 can remain because cannot couple to S = 3/2 with intermediate S = 0
            for(int e = 0;e < Tools::gL()*Tools::gL();++e)
               (*this)(B,i,j) += 2.0 * dpm(0,0,a,b,e,0,c,d,e);

            (*this)(B,i,j) *= scale;

         }
      }

   }

   //then the S = 1 part:
   for(int B = Tools::gL2();B < Tools::gM();++B){

      for(int i = 0;i < gdim(B);++i){

         a = t2s[B][i][0];
         b = t2s[B][i][1];

         for(int j = i;j < gdim(B);++j){

            c = t2s[B][j][0];
            d = t2s[B][j][1];

            (*this)(B,i,j) = 0.0;

            for(int Z = 0;Z < 2;++Z){//loop over the dpm blocks: S = 1/2 and 3/2 = Z + 1/2

               ward = 0.0;

               for(int e = 0;e < Tools::gL()*Tools::gL();++e)
                  ward += dpm(Z,1,a,b,e,1,c,d,e);

               ward *= (2 * (Z + 0.5) + 1.0)/3.0;

               (*this)(B,i,j) += ward;

            }

            (*this)(B,i,j) *= scale;

         }
      }

   }

   this->symmetrize();

}


/**
 * The spincoupled T2-down map that maps a PPHM on a TPM object.
 * @param pphm Input PPHM object
 */
void TPM::T(const PPHM &pphm){

   //first make the bar tpm
   TPM tpm;
   tpm.bar(1.0,pphm);

   //then make the bar phm
   PHM phm;
   phm.bar(1.0,pphm);

   //also make the bar spm with the correct scale factor
   SPM spm;
   spm.bar(0.5/(Tools::gN() - 1.0),pphm);

   int a,b,c,d;
   int a_,b_,c_,d_;
   int sign;

   double norm;

   int S;

   for(int B = 0;B < gnr();++B){//loop over the blocks

      S = block_char[B][0];

      sign = 1 - 2*S;

      for(int i = 0;i < gdim(B);++i){

         a = t2s[B][i][0];
         b = t2s[B][i][1];

         //and for access to the hole elements:
         a_ = Hamiltonian::gbar(a);
         b_ = Hamiltonian::gbar(b);

         for(int j = i;j < gdim(B);++j){

            c = t2s[B][j][0];
            d = t2s[B][j][1];

            c_ = Hamiltonian::gbar(c);
            d_ = Hamiltonian::gbar(d);

            //determine the norm for the basisset
            norm = 1.0;

            if(S == 0){

               if(a == b)
                  norm /= std::sqrt(2.0);

               if(c == d)
                  norm /= std::sqrt(2.0);

            }

            //first the tp part
            (*this)(B,i,j) = tpm(B,i,j);

            //sp part is diagonal for translationaly invariance
            if(i == j)
               (*this)(B,i,j) += spm[a_] + spm[b_];

            for(int Z = 0;Z < 2;++Z){

               (*this)(B,i,j) -= norm * (2.0 * Z + 1.0) * Tools::g6j(0,0,S,Z) * ( phm(Z,d,a_,b,c_) + sign * phm(Z,d,b_,a,c_) 
               
                     + sign * phm(Z,c,a_,b,d_) +  phm(Z,c,b_,a,d_) );

            }

         }
      }

   }

   this->symmetrize();

}

/**
 * The bar function that maps a PPHM object onto a TPM object by tracing away the last pair of incdices of the PPHM
 * @param pphm Input PPHM object
 */
void TPM::bar(double scale,const PPHM &pphm){

   int a,b,c,d;
   int Z;

   double ward;

   for(int B = 0;B < gnr();++B){//loop over the tp blocks

      Z = block_char[B][0];//spin of the TPM - block

      for(int i = 0;i < gdim(B);++i){

         a = t2s[B][i][0];
         b = t2s[B][i][1];

         for(int j = i;j < gdim(B);++j){

            c = t2s[B][j][0];
            d = t2s[B][j][1];

            (*this)(B,i,j) = 0.0;

            for(int S = 0;S < 2;++S){//loop over three particle spin: 1/2 and 3/2

               ward = (2.0*(S + 0.5) + 1.0)/(2.0*Z + 1.0);

               for(int e = 0;e < Tools::gL()*Tools::gL();++e)
                  (*this)(B,i,j) += ward * pphm(S,Z,a,b,e,Z,c,d,e);

            }

            (*this)(B,i,j) *= scale;

         }

      }

   }

   this->symmetrize();

}

/**
 * Collaps a SUP matrix S onto a TPM matrix like this:\n\n
 * sum_i Tr (S u^i)f^i = this
 * @param option = 0, project onto full symmetric matrix space, = 1 project onto traceless symmetric matrix space
 * @param S input SUP
 */
void TPM::collaps(int option,const SUP &S){

   *this = S.gI();

   TPM hulp;

   hulp.Q(1,S.gQ());

   *this += hulp;

#ifdef __G_CON
   hulp.G(S.gG());

   *this += hulp;
#endif

#ifdef __T1_CON
   hulp.T(S.gT1());

   *this += hulp;
#endif

#ifdef __T2_CON
   hulp.T(S.gT2());

   *this += hulp;
#endif

   if(option == 1)
      this->proj_Tr();

}

/**
 * orthogonal projection onto the space of traceless matrices
 */
void TPM::proj_Tr(){

   double ward = (2.0 * this->trace())/(Tools::gM()*(Tools::gM() - 1.0));

   for(int B = 0;B < gnr();++B)
      for(int i = 0;i < gdim(B);++i)
         (*this)(B,i,i) -= ward;

}

/**
 * ( Overlapmatrix of the U-basis ) - map, maps a TPM onto a different TPM, this map is actually a Q-like map
 * for which the paramaters a,b and c are calculated in primal_dual.pdf. Since it is a Q-like map the inverse
 * can be taken as well.
 * @param option = 1 direct overlapmatrix-map is used , = -1 inverse overlapmatrix map is used
 * @param tpm_d the input TPM
 */
void TPM::S(int option,const TPM &tpm_d){

   this->Q(option,Sa,Sb,Sc,tpm_d);

}
