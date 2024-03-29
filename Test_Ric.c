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
# define STUDENT_H_FILE "000000_extra.h"
#endif


//
// include files
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "elapsed_time.h"
#include STUDENT_H_FILE

//Sorts

  void swap(integer_t *a, integer_t *b) {
    integer_t temp = *a;
    *a = *b;
    *b = temp;
  }

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
     
  //
// 

//Algorithms
/*
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
    //Create Array(x[]) composed by all the possible sums of a given Array(a[])
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
    //Algorithm itself
      int mitm(int n, integer_t *p, integer_t desired_sum){
    
        // Get sub-arrays sizes
        int size_X = 1<<(n/2);
        int size_Y = 1<<(n-n/2);
        
        // Allocate space for them
        integer_t *X = malloc(size_X*sizeof(integer_t));
        integer_t *Y = malloc(size_Y*sizeof(integer_t)); 

        // Create the sub arrays
        calcsubarray(p, X, n/2, 0);
        calcsubarray(p, Y, n-n/2, n/2);
        
        // Sort them
        heapSort(X, size_X);
        heapSort(Y, size_Y);
        
         
        int i= 0;
        int j= size_Y - 1;
        while(i< size_X && j >= 0){

          integer_t s = X[i]+Y[j]; 
          if(s == desired_sum ){ 
        
            free(X);    //freeing the space of the arrays
            free(Y);
            return 1;   // return 1 if its found
          }else if(s < desired_sum){
            i++;
          }else{
            j--;
          } 
        }

        // return 0 if is not found
        return 0;
      }
    //
  //

  // Faster Meet in the middle Algorithm (The difference being in the fact that this one sort the arrays while creating them)

    //Create a Sorted Array(x[]) composed by all the possible sums of a given Array(a[])
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
    //Algorithm itself
      int faster_mitm(int n, integer_t *p, integer_t desired_sum){
      
        // Get sub-arrays sizes
        int size_X = 1<<(n/2);
        int size_Y = 1<<(n-n/2);
          
        // Allocate space for them
        integer_t *X = malloc(size_X*sizeof(integer_t));
        integer_t *Y = malloc(size_Y*sizeof(integer_t)); 

        for(int i=0;i<size_X;i++){
          X[i]=0;
        }  
        for(int i=0;i<size_Y;i++){
          Y[i]=0;
        } 
              
        // Create the sub arrays and Sort them
        faster_calcsubarray(p, X, n/2, 0);
        faster_calcsubarray(p, Y, n-n/2, n/2); 
          
         
        int i= 0;
        int j= size_Y - 1;
        while(i< size_X && j >= 0){
          integer_t s = X[i]+Y[j];
              
          if(s == desired_sum){ 
                  
            free(X);    //freeing the space of the arrays
            free(Y);                 
            return 1;   // return 1 if its found
          }else if(s < desired_sum){
            i++;
          }else{
            j--;
          }  
        }

        // return 0 if is not found
        return 0;
      }
    //
  //
  
  // Schroeppel and Shamir technique

    //Create Array(x[][2]) composed by all the possible sums of a given Array(a[]) and the cominations
      void calcsubarray2d(integer_t a[], integer_t (*x)[2], int n, int c){

        integer_t s;
        
        for (int i=0; i<(1<<n); i++)
        {
          s = 0;
          integer_t cmb = 0;
          
          for (int j=0; j<n; j++){
            if (i & (1<<j)){
              s += a[j+c];
              cmb += pow(2, j+c);
            }    
          }

          if(s >= 0){
            x[i][0] = s;
            x[i][1] = cmb;  
          }
        }
      } 
    //

    // Heapify for Min Heap
      void min_heapify(int (*arr)[2], int n, integer_t (*A)[2] , integer_t (*B)[2] ,int i){

        int smallest = i; // Initialize largest as root
        int l = 2 * i + 1; // down left = 2*i + 1
        int r = 2 * i + 2; // down right = 2*i + 2
    
        // If left child is larger than root
        if (l < n && A[arr[l][0]][0] + B[arr[l][1]][0]  <= A[arr[smallest][0]][0] + B[arr[smallest][1]][0])
          smallest = l;
    
        // If right child is larger than largest so far
        if (r < n && A[arr[r][0]][0] + B[arr[r][1]][0]  <= A[arr[smallest][0]][0] + B[arr[smallest][1]][0])
          smallest = r;
    
        // If largest is not root
        if (smallest != i) {

          int temp_array_0 = arr[smallest][0]; int temp_array_1 = arr[smallest][1];
          arr[smallest][0] =arr[i][0];
          arr[smallest][1] =arr[i][1];
          arr[i][0] = temp_array_0;
          arr[i][1] = temp_array_1;

          // Recursively heapify the affected sub-tree
          min_heapify(arr, n, A, B ,smallest);
        }
      }
    //  
    
    // Create Min Heap
      void generateMinHeap(int (*minheap)[2] , integer_t (*A)[2],integer_t (*B)[2],int size_a,int size_b){

        // Generate tree from array indexes
        int heap_iter = 0;
        for(int i = 0; i < size_a; i++){
          for (int j = 0; j < size_b; j++){
            minheap[heap_iter][0] = i;
            minheap[heap_iter][1] = j;
            heap_iter+=1;
          }
        }
        
        // Perform reverse level order traversal
        // from last non-leaf node and heapify
        int startIdx = ( (size_a*size_b) / 2) - 1;

        // each node
        for (int i = startIdx; i >= 0; i--) {
          min_heapify(minheap, size_a*size_b , A, B ,i);
        }

      }
    //

    // Heapify for Max Heap
      void max_heapify(int (*arr)[2], int n, integer_t (*A)[2] , integer_t (*B)[2] ,int i){

        int largest = i; // Initialize largest as root
        int l = 2 * i + 1; // down left = 2*i + 1
        int r = 2 * i + 2; // down right = 2*i + 2
  
        // If left child is larger than root
        if (l < n && A[arr[l][0]][0] + B[arr[l][1]][0]  >= A[arr[largest][0]][0] + B[arr[largest][1]][0])
          largest = l;
  
        // If right child is larger than largest so far
        if (r < n && A[arr[r][0]][0] + B[arr[r][1]][0]  >= A[arr[largest][0]][0] + B[arr[largest][1]][0])
          largest = r;
  
        // If largest is not root
        if (largest != i) {
          int temp_array_0 = arr[largest][0]; int temp_array_1 = arr[largest][1];
          arr[largest][0] =arr[i][0];
          arr[largest][1] =arr[i][1];
          arr[i][0] = temp_array_0;
          arr[i][1] = temp_array_1;

          // Recursively heapify the affected sub-tree
          max_heapify(arr, n, A, B ,largest);
        }
      }
    //

    // Create Max Heap
      void generateMaxHeap(int (*maxheap)[2] , integer_t (*C)[2],integer_t (*D)[2],int size_c,int size_d){

        // Generate tree from array indexes
        int heap_iter = 0;
        for(int i = 0; i < size_c; i++){
          for (int j = 0; j < size_d; j++){
            maxheap[heap_iter][0] = i;
            maxheap[heap_iter][1] = j;
            heap_iter+=1;
          }
        }

        // Perform reverce level order traversal
        // from last non-leaf node and heapify
        int startIdx = ( (size_c*size_d) / 2) - 1;
        // each node
        for (int i = startIdx; i >= 0; i--) {
          max_heapify(maxheap, size_c*size_d , C, D ,i);
        }
      }
    //

    //Algorithm itself
      integer_t SS(int n, integer_t *p, integer_t desired_sum){

        // pega o tamanho
        int a = (n/2)/2;
        int b = ((n/2) - (n/2)/2);
        int c = (n - n/2)/2;
        int d = ((n - n/2) - (n - n/2)/2);

        // arranja espaço para as somas
        integer_t (*A)[2] = malloc(sizeof(integer_t)*(1<<a) * 2);
        integer_t (*B)[2] = malloc(sizeof(integer_t)*(1<<b) * 2);
        integer_t (*C)[2] = malloc(sizeof(integer_t)*(1<<c) * 2);
        integer_t (*D)[2] = malloc(sizeof(integer_t)*(1<<d) * 2);

        calcsubarray2d(p, A, a, 0);
        calcsubarray2d(p, B, b, a);
        calcsubarray2d(p, C, c, a + b);
        calcsubarray2d(p, D, d, a + b + c);

        heapSort2d(A,(1<<a));
        heapSort2d(B,(1<<b));
        heapSort2d(C,(1<<c));
        heapSort2d(D,(1<<d));

        int (*minheap)[2] = malloc(sizeof(int)* 2 * ((1<<a) * (1<<b)));

        generateMinHeap(minheap, A, B,(1<<a),(1<<b));

        int (*maxheap)[2] = malloc(sizeof(int)* 2 * ((1<<c) * (1<<d)));

        generateMaxHeap(maxheap, C, D,(1<<c),(1<<d));

        integer_t soma;
        integer_t soma2;
        for (int i = 0;i<((1<<c) * (1<<d));i++){
          if (C[maxheap[i][0]][0] + D[maxheap[i][1]][0] > desired_sum){
            continue;
          }

          for (int j = 0;j<((1<<a) * (1<<b));j++){
            soma = A[minheap[j][0]][0] + B[minheap[j][1]][0] + C[maxheap[i][0]][0] + D[maxheap[i][1]][0]; 
            soma2 = A[minheap[j][0]][1] + B[minheap[j][1]][1] + C[maxheap[i][0]][1] + D[maxheap[i][1]][1]; 
            if (soma == desired_sum){
              free(A);
              free(B);
              free(C);
              free(D);
              free(minheap);
              free(maxheap); 
              return soma2;
            }
          }
        }
          
           
        free(A);
        free(B);
        free(C);
        free(D);
        free(minheap);
        free(maxheap);
        return 0;
      }
    //

  //
  */
  // Schroeppel and Shamir technique by Teles
     
    integer_t currentMinHeapSize = 0;
    integer_t currentMaxHeapSize = 0;
    integer_t minHeap[1 << 23];
    integer_t maxHeap[1 << 23];
    integer_t indicesMin[1 << 23][2];
    integer_t indicesMax[1 << 23][2];
    integer_t iMin, iMax, jMin, jMax;
    // This works, it swaps two elements of an array
    void swapMin(integer_t i, integer_t j){
    
        integer_t temp = minHeap[i];
        minHeap[i] = minHeap[j];
        minHeap[j] = temp;
    
        temp = indicesMin[i][0];
        indicesMin[i][0] = indicesMin[j][0];
        indicesMin[j][0] = temp;
    
        temp = indicesMin[i][1];
        indicesMin[i][1] = indicesMin[j][1];
        indicesMin[j][1] = temp;
    }
    // This works, it swaps two elements of an array
    void swapMax(integer_t i, integer_t j){
    
        integer_t temp = maxHeap[i];
        maxHeap[i] = maxHeap[j];
        maxHeap[j] = temp;
    
        temp = indicesMax[i][0];
        indicesMax[i][0] = indicesMax[j][0];
        indicesMax[j][0] = temp;
    
        temp = indicesMax[i][1];
        indicesMax[i][1] = indicesMax[j][1];
        indicesMax[j][1] = temp;
    }
     
    // Returns left child node index
    integer_t leftChild(integer_t i){
    
        return i * 2 + 1;
    }
    // Returns right child node index
    integer_t rightChild(integer_t i){
    
        return i * 2 + 2;
    }
     
     
    void heapify(integer_t i){
    
        // If this node isn't a leaf and is greater than any of its children
        if (!(i + 1 > (currentMinHeapSize / 2) && i < currentMinHeapSize))
        {
            if (minHeap[i] > minHeap[i * 2 + 1] || minHeap[i] > minHeap[i * 2 + 2])
            {
                // Swap with the left child and heapify him
                if (minHeap[i * 2 + 1] < minHeap[i * 2 + 2])
                {
                    swapMin(i, i * 2 + 1);
                    heapify(i * 2 + 1);
                }
                else
                {
                    swapMin(i, i * 2 + 2);
                    heapify(i * 2 + 2);
                }
            }
        }
    }

    void insertHeap(integer_t number, integer_t i, integer_t j){
    
        minHeap[currentMinHeapSize] = number;
        indicesMin[currentMinHeapSize][0] = i;
        indicesMin[currentMinHeapSize][1] = j;
        // Index of where we are at
        integer_t index = currentMinHeapSize;
        currentMinHeapSize++;
    
        // If the child is bigger than the parent, swap, and then go to the parent's index
        while (minHeap[index] < minHeap[(index - 1) / 2])
        {
            swapMin(index, (index - 1) / 2);
            index = (index - 1) / 2;
        }
    }

    integer_t pop(){
    
        integer_t pop = minHeap[0];
        minHeap[0] = minHeap[currentMinHeapSize - 1];
        indicesMin[0][0] = indicesMin[currentMinHeapSize - 1][0];
        indicesMin[0][1] = indicesMin[currentMinHeapSize - 1][1];
        //printf("What happens at pop() when currentSize = %d?\n", *currentSize);
        //printArraysMultidimensional(*currentSize, indicesArray);
    
        //printArraysInt(*currentSize, heap);
    
        currentMinHeapSize--;
        heapify(0);
    
        return pop;
    }

    void maxHeapify(integer_t i){
    
        // If this node isn't a leaf and is smaller than any of its children
        if (!(i + 1 > (currentMaxHeapSize / 2) && i < currentMaxHeapSize))
        {
            if (maxHeap[i] < maxHeap[i * 2 + 1] || maxHeap[i] < maxHeap[i * 2 + 2])
            {
                // Swap with the left child and heapify him
                if (maxHeap[i * 2 + 1] > maxHeap[i * 2 + 2])
                {
                    swapMax(i, i * 2 + 1);
                    maxHeapify(i * 2 + 1);
                }
                else
                {
                    swapMax(i, i * 2 + 2);
                    maxHeapify(i * 2 + 2);
                }
            }
        }
    }

    integer_t maxPop(){
    
        integer_t pop = maxHeap[0];
        maxHeap[0] = maxHeap[currentMaxHeapSize - 1];
        indicesMax[0][0] = indicesMax[currentMaxHeapSize - 1][0];
        indicesMax[0][1] = indicesMax[currentMaxHeapSize - 1][1];
    
        currentMaxHeapSize--;
        maxHeapify(0);
    
        return pop;
    }
     
    void insertMaxHeap(integer_t number, integer_t i, integer_t j){
    
        maxHeap[currentMaxHeapSize] = number;
        indicesMax[currentMaxHeapSize][0] = i;
        indicesMax[currentMaxHeapSize][1] = j;
        // Index of where we are at
        integer_t index = currentMaxHeapSize;
        currentMaxHeapSize++;
    
        // If the child is bigger than the parent, swap, and then go to the parent's index
        while (maxHeap[index] > maxHeap[(index - 1) / 2] && index > 0)
        {
            //printf("Inside while loop\n");
            swapMax(index, (index - 1) / 2);
            index = (index - 1) / 2;
        }
    }
     
    integer_t *mergerSort(int nSize, integer_t *pArray, integer_t *a){
      
      integer_t i1, j1, k1;
      a[(1 << nSize) - 2] = 0;
      a[(1 << nSize) - 1] = pArray[0];
      for (int i = 1; i < nSize; i++)
      {
          i1 = (1 << nSize) - (1 << i);
          j1 = i1;
          k1 = (1 << nSize) - (2 << i);
          while (i1 < (1 << nSize) && j1 < (1 << nSize))
          {
              if (a[i1] <= a[j1] + pArray[i])
              {
                  a[k1] = a[i1];
                  k1++;
                  i1++;
              }
              else
              {
                  a[k1] = a[j1] + pArray[i];
                  k1++;
                  j1++;
              }
          }
          while (j1 < (1 << nSize))
          {
              a[j1] += pArray[i];
              j1++;
          }
      }

      return a;
    }
    
    //Algorithm itself
      int SS_T(int n, integer_t *p, integer_t desired_sum){
 
        int L1Size;
        if (n % 4 != 0)
            L1Size = n / 4 + 1;
        else
            L1Size = n / 4;
        int L2Size = n / 4 + (n % 4) / 2;
        int R2Size = n / 4 + (n % 4) / 3;
        int R1Size = n / 4;
    
        // Create arrays with correct sizes
        // This is working correctly
        integer_t L1[L1Size];
        integer_t L2[L2Size];
        integer_t R2[R2Size];
        integer_t R1[R1Size];
    
        for (int i = 0; i < L1Size; i++)
            L1[i] = p[i];
        for (int i = L1Size; i < L1Size + L2Size; i++)
            L2[i - L1Size] = p[i];
        for (int i = L1Size + L2Size; i < L1Size + L2Size + R2Size; i++)
            R2[i - (L1Size + L2Size)] = p[i];
        for (int i = L1Size + L2Size + R2Size; i < L1Size + L2Size + R2Size + R1Size; i++)
            R1[i - (L1Size + L2Size + R2Size)] = p[i];
    
        integer_t fancyL1Size = 1 << L1Size;
        integer_t fancyL2Size = 1 << L2Size;
        integer_t fancyR2Size = 1 << R2Size;
        integer_t fancyR1Size = 1 << R1Size;
    
        integer_t fancyL1[fancyL1Size];
        integer_t fancyL2[fancyL2Size];
        integer_t fancyR2[fancyR2Size];
        integer_t fancyR1[fancyR1Size];
    
        mergerSort(L1Size, L1, fancyL1);
        mergerSort(L2Size, L2, fancyL2);
        mergerSort(R2Size, R2, fancyR2);
        mergerSort(R1Size, R1, fancyR1);

         

           
        for (int k = 0; k < fancyL1Size; k++)
          insertHeap(fancyL1[k] /* + fancyL2[0] which is always 0*/, k, 0);
        // Do the same, but j = fancyR1Size-1
        for (int k = 0; k < fancyR2Size; k++)
          insertMaxHeap(fancyR2[k] + fancyR1[fancyR1Size - 1], k, fancyR1Size - 1);
    
        integer_t K = pow(2,n);
        integer_t minTop, maxTop;
 
        for (integer_t count = 0; count < K; count++){
      
          minTop = minHeap[0];
          maxTop = maxHeap[0];
  
          iMin = indicesMin[0][0], jMin = indicesMin[0][1];
          iMax = indicesMax[0][0], jMax = indicesMax[0][1];
  
          integer_t sum = minTop + maxTop;
  
          if (sum == desired_sum){
             
            return 1;
          }
          else if (sum < desired_sum){
          
            pop();
            if (jMin + 1 < fancyL2Size)
                insertHeap(fancyL1[iMin] + fancyL2[jMin + 1], iMin, jMin + 1);
          }
          else
          {
            maxPop();
            if (jMax > 0)
              insertMaxHeap(fancyR2[iMax] + fancyR1[jMax - 1], iMax, jMax - 1);
          }
        }
        
        return 0;
      }
    //
     
  //
//

//Converter (converts a decimal number to binary, used to get the commbinations)
  char *Converter(int n,integer_t x, char *sol){
    for(int bit=0; bit<n ;bit++){
        integer_t mask = (1ull<<bit);   
        if((x & mask)!=0){ sol[bit]='1';}else{ sol[bit]='0';}
    }
    sol[n]='\0';

    return sol;
  }
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
  
 
  /* Setting up file */
  FILE *fp_1 = NULL;
  FILE *fp_2 = NULL;
  FILE *fp_3 = NULL;
  FILE *fp_4 = NULL;
  FILE *fp_5 = NULL;
  FILE *fp_6 = NULL;
  FILE *fp_7 = NULL;
  FILE *fp_8 = NULL;
  FILE *fp_9 = NULL;
  FILE *fp_10 = NULL;
  FILE *fp_11 = NULL;
  FILE *fp_12 = NULL;
  FILE *fp_13 = NULL;
  FILE *fp_14 = NULL;
  fp_1 = fopen("data_1.log", "a");
  fp_2 = fopen("data_1_max.log", "a");
  fp_3 = fopen("data_2.log", "a");
  fp_4 = fopen("data_2_max.log", "a");
  fp_5 = fopen("data_3.log", "a");
  fp_6 = fopen("data_3_max.log", "a");
  fp_7 = fopen("data_4.log", "a");
  fp_8 = fopen("data_4_max.log", "a");
  fp_9 = fopen("data_5.log", "a");
  fp_10 = fopen("data_5_max.log", "a");
  fp_11 = fopen("data_6.log", "a");
  fp_12 = fopen("data_6_max.log", "a");
  fp_13 = fopen("data_7.log", "a");
  fp_14 = fopen("data_7_max.log", "a");
    
 

  // Loop for n's
  for(int i = 0;i <n_problems;i++){
  
    printf("--------------------------- \n");
                
    // get n and p
    int n = all_subset_sum_problems[i].n;  
    printf("n =  %i\n\n",n);
    integer_t *p = all_subset_sum_problems[i].p;     

    // Set up to get the combinations 
    char comb_bin[n+1];

    // Setting up time variables    
    double dt_bf_i = 0;  
    double dt_bf_i_max = 0;  
    double dt_bf_r = 0;
    double dt_bf_r_max = 0;
    double dt_bf_i_s = 0;
    double dt_bf_i_s_max = 0;
    double dt_mitm = 0;     
    double dt_mitm_max = 0;
    double dt_f_mitm = 0;     
    double dt_f_mitm_max = 0;
    double dt_ss = 0;     
    double dt_ss_max = 0;

    // Loop for sum's
    for(int j = 0;j < n_sums;j++){
           
      integer_t sum = all_subset_sum_problems[i].sums[j];
      double tmp_dt;

      /*
      // Iterative
      tmp_dt = cpu_time();   
      int comb = Bf_Iter(n, p, sum);
      tmp_dt = cpu_time() - tmp_dt;
      if(tmp_dt > dt_bf_i_max){
        dt_bf_i_max = tmp_dt;
      }
      dt_bf_i += tmp_dt;


      // Recursive
      tmp_dt = cpu_time();   
      int comb_rec= Bf_recur(n,0, p,0,0,sum);
      tmp_dt = cpu_time() - tmp_dt;        
      if(tmp_dt > dt_bf_r_max){
        dt_bf_r_max = tmp_dt;
      }
      dt_bf_r += tmp_dt;


      // Recursive Smart
      tmp_dt = cpu_time();   
      integer_t comb_smart= Bf_recur_smart(n,n-1, p,0,0,sum);
      tmp_dt = cpu_time() - tmp_dt;
      if(tmp_dt > dt_bf_i_s_max){
        dt_bf_i_s_max = tmp_dt;
      }
      dt_bf_i_s += tmp_dt;
            
            
      // Meet in the middle
      tmp_dt = cpu_time();   
      int x= mitm(n, p, sum);
      tmp_dt = cpu_time() - tmp_dt;
      if(tmp_dt > dt_mitm_max){
        dt_mitm_max = tmp_dt;
      }
      dt_mitm += tmp_dt;

            
      // Fast Meet in the middle
      tmp_dt = cpu_time();   
      int y= faster_mitm(n, p, sum);
      tmp_dt = cpu_time() - tmp_dt;
      if(tmp_dt > dt_f_mitm_max){
        dt_f_mitm_max = tmp_dt;
      }
      dt_f_mitm += tmp_dt;
        
      */
      // Schroeppel and Shamir technique

      tmp_dt = cpu_time();    
      int z= SS_T(n, p, sum);
      tmp_dt = cpu_time() - tmp_dt;
      if(tmp_dt > dt_ss_max){
        dt_ss_max = tmp_dt;
      }
      dt_ss += tmp_dt;
      
 
      // print results
      printf("-------------------------------------------------\n");/*
      printf("Brute force                           %d,  %lld || %i -> %s  \n", j ,sum,comb,   Converter(n, comb, comb_bin));
      printf("Brute force recursiva                 %d,  %lld || %i -> %s  \n", j ,sum,comb_rec,   Converter(n, comb_rec, comb_bin));
      printf("Brute force recur smart               %d,  %lld || %lld -> %s  \n", j ,sum,comb_smart,   Converter(n, comb_smart, comb_bin));  
      printf("Meet in the middle                    %d,  %lld || %i  \n", j ,sum, x);
      printf("Faster meet in the middle             %d,  %lld || %i  \n", j ,sum, y);*/
      printf("Schroeppel and Shamir technique       %d|| %d \n", j, z); 
            
    }

    // store times (commented beacause we already took the data needed)
    /*
    fprintf(fp_1,"%i %f \n",n, dt_bf_i/20);
    fprintf(fp_2,"%i %f \n",n, dt_bf_i_max);
    fprintf(fp_3,"%i %f \n",n, dt_bf_r/20);
    fprintf(fp_4,"%i %f \n",n, dt_bf_r_max);
    fprintf(fp_5,"%i %f \n",n, dt_bf_i_s/20);
    fprintf(fp_6,"%i %f \n",n, dt_bf_i_s_max);
    fprintf(fp_7,"%i %f \n",n, dt_mitm/20);
    fprintf(fp_8,"%i %f \n",n, dt_mitm_max);
    fprintf(fp_9,"%i %f \n",n, dt_f_mitm/20);
    fprintf(fp_10,"%i %f \n",n, dt_f_mitm_max);
    fprintf(fp_11,"%i %f \n",n, dt_ss/20);
    fprintf(fp_12,"%i %f \n",n, dt_ss_max);
    fprintf(fp_13,"%i %f \n",n, dt_ss/20);
    fprintf(fp_14,"%i %f \n",n, dt_ss_max);*/

  }     


   
  return 0;
}