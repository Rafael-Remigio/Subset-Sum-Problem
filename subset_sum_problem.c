//
// AED, November 2021
//
// Solution of the first practical assignement (subset sum problem)
//
// Place your student numbers and names here
//

#if __STDC_VERSION__ < 199901L
# error "This code must must be compiled in c99 mode or later (-std=c99)" // to handle the unsigned long long data type
#endif
#ifndef STUDENT_H_FILE
# define STUDENT_H_FILE "000000.h"
#endif


//
// include files
//

#include <stdio.h>
#include <stdlib.h>
#include "elapsed_time.h"
#include STUDENT_H_FILE

//Sorts

  //Bubble Sort
    void bubbleSort(integer_t arr[], int n)
    {
        // used just for the graphs and to show how effective heapSort is

      int i, j;
      int swapped;
      for (i = 0; i < n-1; i++)
      {
        swapped = 0;
        for (j = 0; j < n-i-1; j++)
        {
            if (arr[j] > arr[j+1]) 
            {
              swap(&arr[j], &arr[j+1]);
              swapped = 1;
            }
        }
    
        if (swapped == 0)
            break;
      }
    }
  //

  //Heap Sort
    void swap(integer_t *a, integer_t *b) {
        integer_t temp = *a;
        *a = *b;
        *b = temp;
    }
      
    void heapify(integer_t arr[], int n, int i) {
        // Find largest among root, left child and right child
        int largest = i;
        int left = 2 * i + 1;
        int right = 2 * i + 2;
      
        if (left < n && arr[left] > arr[largest])
          largest = left;
      
        if (right < n && arr[right] > arr[largest])
          largest = right;
      
        // Swap and continue heapifying if root is not largest
        if (largest != i) {
          swap(&arr[i], &arr[largest]);
          heapify(arr, n, largest);
        }
    }
      

    void heapSort(integer_t arr[], int n) {
        // Build max heap
        for (int i = n / 2 - 1; i >= 0; i--)
          heapify(arr, n, i);
      
        // Heap sort
        for (int i = n - 1; i >= 0; i--) {
          swap(&arr[0], &arr[i]);
      
          // Heapify root element to get highest element at root again
          heapify(arr, i, 0);
        }
    }
  //
// 

//Algorithms

  // Brute Force Iterative Algorithm 
    int Bf_Iter(int n, integer_t *p, integer_t desired_sum){

      int comb = 0;
      integer_t test_sum;

      for(comb;comb<(1<<n);comb++){   // changes the combination
        
        test_sum =0;

        for(int bit=0; bit<n ;bit++){
          
          int mask = (1<<bit);
          if((comb & mask)!=0){test_sum += p[bit];}   
        }   

        if(test_sum == desired_sum){break;} // tests if the sum is equal to the desired sum 
      }

      return comb;
        
    }
  //

  // Brute Force Recursive Algorithm  
    int Bf_recur( unsigned int n,unsigned int m,integer_t *p,double sum, int comb,integer_t desired_sum){
        
      if(m == n){ // if is in last step of recursion
        
        if (sum == desired_sum){
          return comb;    // return combination different from 0 if found
                
        }
        else return 0;      // else returns 0
      }

        
      int result = Bf_recur(n,m + 1u,p,sum , comb ,desired_sum);  
      if (result == 0){         // if combination is 0 it means we must take another path          
    
        return Bf_recur(n,m + 1u,p,sum + p[m] ,comb+ pow(2,m),desired_sum);      
      }   

      return result;
    }
  //

  // Smart Brute Force Recursive Algorithm  
    integer_t Bf_recur_smart( unsigned int n,int m,integer_t *p,integer_t sum, integer_t comb,integer_t desired_sum)
    {  // this function run trougth the array from bigger to smaller and if the current_sum is bigger then the desired_sum we exit, saving lots of time

      if (sum == desired_sum){    // found the solucion so we return the combination
        return comb;
      }

      if(sum > desired_sum){ // if current sum > desired_sum this branch does not contain the solucion
        return 0;
      }

      if(m == -1) {               // if we reached the end of the array
        return 0 ;
      }

      integer_t power = (1ull);
      integer_t result = Bf_recur_smart(n,m-1,p,sum + p[m], comb + (power<<m),desired_sum);  
        
      if (result == 0)  {             // if combination is 0 it means we must take another path         
        return Bf_recur_smart(n,m-1,p,sum  ,comb,desired_sum);      
      }   
      return result;
    }
  //
  
  // Meet in the middle Algorithm
    //Create Arrays needed to the mitm
      void calcsubarray(integer_t a[], integer_t x[], int n, int c){

        integer_t s;
        for (int i=0; i<(1<<n); i++)
        {
          s = 0;
          for (int j=0; j<n; j++){
            if (i & (1<<j)){
              s += a[j+c];
            }    
          }

          if(s >= 0){
            x[i] = s;  
          }
        }
      } 
    //Algorithm
      int mitm(int n, integer_t *p, integer_t desired_sum){
    
        // pega o tamanho
        int size_X = 1<<(n/2);
        int size_Y = 1<<(n-n/2);
        
        // arranja espaço para as somas
        integer_t *X = malloc(size_X*sizeof(integer_t));
        integer_t *Y = malloc(size_Y*sizeof(integer_t)); 

        // enche os arrays com as somas respetivas
        calcsubarray(p, X, n/2, 0);
        calcsubarray(p, Y, n-n/2, n/2);
        
        // Sorta os arrays 
        heapSort(X, size_X);
        heapSort(Y, size_Y);
        
        // loopa comparando e tal (como o stor explicou)

        int i= 0;
        int j= size_Y - 1;
        while(i< size_X && j >= 0){
            integer_t s = X[i]+Y[j]; 
            if(s == desired_sum ){ 
                free(X);
                free(Y);
                return 1;   
            }else if(s < desired_sum){
                i++;
            }else{
                j--;
            } 
            
            
        }

        return 0;
      }
    //
  //

  // Faster Meet in the middle Algorithm (The difference being in the fact that this one sort the arrays while creating them)

    //Create Arrays needed to the mitm
      void faster_calcsubarray(integer_t a[], integer_t x[], int n, int c){

        for (int i=0; i<n; i++)
        {

          integer_t i1 = (1<<n)-(1<<i);
          integer_t j1 = i1;
          integer_t k1 = (1<<n)-(2<<i);

          while(i1 < (1<<n)){

            if(x[i1] <= x[j1] + a[i+c]){
              x[k1++] = x[i1++];
            }else{
              x[k1++] = x[j1++] + a[i+c];   
            }
          }
          while(j1 < (1<<n)){
            x[j1++] += a[i+c];
          }
        }

      } 
    //Algorithm
      int faster_mitm(int n, integer_t *p, integer_t desired_sum){
      
          // pega o tamanho
          int size_X = 1<<(n/2);
          int size_Y = 1<<(n-n/2);
          
          // arranja espaço para as somas
          integer_t *X = malloc(size_X*sizeof(integer_t));
          integer_t *Y = malloc(size_Y*sizeof(integer_t)); 

          for(int i=0;i<size_X;i++){
              X[i]=0;
          }  
          for(int i=0;i<size_Y;i++){
              Y[i]=0;
          } 
              
          // enche os arrays com as somas respetivas
          faster_calcsubarray(p, X, n/2, 0);
          faster_calcsubarray(p, Y, n-n/2, n/2); 
          
        

          // loopa comparando e tal (como o stor explicou)

          int i= 0;
          int j= size_Y - 1;
          while(i< size_X && j >= 0){
              integer_t s = X[i]+Y[j];
              
              if(s == desired_sum){ 
                  
                  free(X);
                  free(Y);
                  
                  return 1;   
              }else if(s < desired_sum){
                  i++;
              }else{
                  j--;
              } 
              
              
          }
        
          return 0;
      }
    //
  //

//


//
// main program
//

int main(void)
{
  fprintf(stderr,"Program configuration:\n");
  fprintf(stderr,"  min_n ....... %d\n",min_n);
  fprintf(stderr,"  max_n ....... %d\n",max_n);
  fprintf(stderr,"  n_sums ...... %d\n",n_sums);
  fprintf(stderr,"  n_problems .. %d\n",n_problems);
  fprintf(stderr,"  integer_t ... %d bits\n",8 * (int)sizeof(integer_t));
  //
  // place your code here

  


  //
  return 0;
}