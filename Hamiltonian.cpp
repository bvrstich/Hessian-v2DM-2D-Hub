#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::endl;
using std::ios;

#include "include.h"

int **Hamiltonian::xy_a;
int **Hamiltonian::a_xy;

/**
 * function that allocates and constructs the lists.
 */
void Hamiltonian::init(){

   //allocate
   xy_a = new int * [Tools::gL()];

   for(int x = 0;x < Tools::gL();++x)
      xy_a[x] = new int [Tools::gL()];

   a_xy = new int * [Tools::gL2()];

   for(int i = 0;i < Tools::gL2();++i)
      a_xy[i] = new int [2];

   //construct:
   int a = 0;

   for(int x = 0;x < Tools::gL();++x)
      for(int y = 0;y < Tools::gL();++y){

         a_xy[a][0] = x;
         a_xy[a][1] = y;

         xy_a[x][y] = a;

         a++;

      }

}

/**
 * deallocates the lists.
 */
void Hamiltonian::clear(){

   //delete xy_a
   for(int x = 0;x < Tools::gL();++x)
      delete [] xy_a[x];

   delete [] xy_a;

   //delete a_xy
   for(int i = 0;i < Tools::gL2();++i)
      delete [] a_xy[i];

   delete [] a_xy;

}

/**
 * print the list
 */
void Hamiltonian::print(){

   for(int i = 0;i < Tools::gL2();++i)
      std::cout << i << "\t" << a_xy[i][0]<< "\t" << a_xy[i][1] << std::endl;

}

/**
 * access to the a_xy list from outside this class
 * @param a the sp-index
 * @param option can be 0 or 1
 * @return k_x if option == 0 , k_y if option == 1
 */
int Hamiltonian::ga_xy(int a,int option){

   return a_xy[a][option];

}

/**
 * access to the list xy_a from outside this class
 * @param k_x the x momentum
 * @param k_y the y momentum
 * @return the sp index corresponding to k_x and k_y
 */
int Hamiltonian::gxy_a(int k_x,int k_y){
   
   return xy_a[k_x][k_y];

}

/**
 * transform the "particle momentum" sp-index to the "hole momentum" sp-index
 * @param a the input sp-index
 */
int Hamiltonian::bar(int a){

   return xy_a[(-a_xy[a][0] + Tools::gL())%Tools::gL()][(-a_xy[a][1] + Tools::gL())%Tools::gL()];

}

/**
 * find the sp index which combines with sp-index k to sum up to sp index K.
 */
int Hamiltonian::adjoint(int K,int k){

  return xy_a[(a_xy[K][0] - a_xy[k][0] + Tools::gL())%Tools::gL()][(a_xy[K][1] - a_xy[k][1] + Tools::gL())%Tools::gL()];
}

/**
 * find the sp index which combines with sp-index e to sum up to sp index formed by (a + k).
 */
int Hamiltonian::adjoint_sum(int a,int k,int e){

  return xy_a[(a_xy[a][0] + a_xy[k][0] - a_xy[e][0] + Tools::gL())%Tools::gL()][(a_xy[a][1] + a_xy[k][1] - a_xy[e][1] + Tools::gL())%Tools::gL()];

}

/**
 * add two sp indices together to form a new sp index (i.e. add the separate x and y momenta and recombine)
 */
int Hamiltonian::add(int a,int b){

   return xy_a[(a_xy[a][0] + a_xy[b][0])%Tools::gL()][(a_xy[a][1] + a_xy[b][1])%Tools::gL()];

}

/**
 * find the sp index which combines with sp-indices a and b to sum up to sp index K.
 */
int Hamiltonian::adjoint(int K,int a,int b){

  return xy_a[(a_xy[K][0] - a_xy[a][0] - a_xy[b][0] + 2*Tools::gL())%Tools::gL()][(a_xy[K][1] - a_xy[a][1] - a_xy[b][1] + 2*Tools::gL())%Tools::gL()];

}

/**
 * add three sp indices together to form a new sp index (i.e. add the separate x and y momenta and recombine)
 */
int Hamiltonian::add(int a,int b,int c){

   return xy_a[(a_xy[a][0] + a_xy[b][0] + a_xy[c][0])%Tools::gL()][(a_xy[a][1] + a_xy[b][1] + a_xy[c][1])%Tools::gL()];

}
