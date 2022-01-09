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

    void heapify2d(integer_t (*arr)[2], int n, int i) {
        // Find largest among root, left child and right child
        int largest = i;
        int left = 2 * i + 1;
        int right = 2 * i + 2;
      
        if (left < n && arr[left][0] > arr[largest][0])
          largest = left;
      
        if (right < n && arr[right][0] > arr[largest][0])
          largest = right;
      
        // Swap and continue heapifying if root is not largest
        if (largest != i) {
          swap(&arr[i][0], &arr[largest][0]);
          swap(&arr[i][1], &arr[largest][1]);
          heapify2d(arr, n, largest);
        }
    }
      

    void heapSort2d(integer_t (*arr)[2], int n) {
        // Build max heap
        for (int i = n / 2 - 1; i >= 0; i--)
          heapify2d(arr, n, i);
      
        // Heap sort
        for (int i = n - 1; i >= 0; i--) {
          swap(&arr[0][0], &arr[i][0]);
          swap(&arr[0][1], &arr[i][1]);
      
          // Heapify root element to get highest element at root again
          heapify2d(arr, i, 0);
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
        
        /* Go through the array X from start to end, and through Y form end to start, testing if we get the desired sum from putting them together */
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
          
        /* Go through the array X from start to end, and through Y form end to start, testing if we get the desired sum from putting them together */
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

  // Schroeppel and Shamir technique by Teles
     
    
     
    integer_t iMaxIndices[1<<23];
    integer_t jMaxIndices[1<<23];
    
    integer_t iMinIndices[1<<23];
    integer_t jMinIndices[1<<23];
    

    integer_t minHeap[1<<23];
    integer_t minSize = 0;
    
    integer_t maxHeap[1<<23];
    integer_t maxSize = 0;

    integer_t iMin, iMax, jMin, jMax;



    void MinHeapify(integer_t i){
     
      int l = 2 * i + 1;
      int r = 2 * i + 2;
     
      // If this node isn't a leaf and is smaller than any of its children
        if (!(i + 1 > (minSize / 2) && i < minSize))
        {
            if (minHeap[i] > minHeap[l] || minHeap[i] > minHeap[r])
            {
                // Swap with the left child and heapify him
                if (minHeap[l] < minHeap[r])
                {
                    swap(&minHeap[i], &minHeap[l]);
    
                    swap(&iMinIndices[i], &iMinIndices[l]);
                    swap(&jMinIndices[i], &jMinIndices[l]);
                
                    MinHeapify(l);
                }
                else
                {
                    swap(&minHeap[i], &minHeap[r]);
    
                    swap(&iMinIndices[i], &iMinIndices[r]);
                    swap(&jMinIndices[i], &jMinIndices[r]);
                
                    MinHeapify(r);
                }
            }
        }
    }
    void InsertMinHeap(integer_t newNum, integer_t i, integer_t j){
       
      minHeap[minSize] = newNum;
      iMinIndices[minSize]=i;
      jMinIndices[minSize]=j;
      integer_t pos = minSize;
      minSize++;
       
      while (minHeap[pos] < minHeap[(pos - 1) / 2]){
        
        swap(&minHeap[pos], &minHeap[(pos - 1) / 2]);
        swap(&iMinIndices[pos], &iMinIndices[(pos - 1) / 2]);
        swap(&jMinIndices[pos], &jMinIndices[(pos - 1) / 2]);
        pos = (pos - 1) / 2;
      }

       
    }
    void MinPop(){
        // replace first node by last and delete last
        minHeap[0] = minHeap[minSize - 1];
        iMinIndices[0] = iMinIndices[minSize - 1];
        jMinIndices[0] = jMinIndices[minSize - 1];
        minSize--;
    
        MinHeapify(0);
    }
    
    
    void MaxHeapify(integer_t i){
     
      integer_t l = 2 * i + 1;
      integer_t r = 2 * i + 2;
      // If this node isn't a leaf and is smaller than any of its children
        if (!(i + 1 > (maxSize / 2) && i < maxSize))
        {
            if (maxHeap[i] < maxHeap[l] || maxHeap[i] < maxHeap[r])
            {
                // Swap with the left child and heapify him
                if (maxHeap[l] > maxHeap[r])
                {
                    swap(&maxHeap[i], &maxHeap[l]);
    
                    swap(&iMaxIndices[i], &iMaxIndices[l]);
                    swap(&jMaxIndices[i], &jMaxIndices[l]);
                
                    MaxHeapify(l);
                }
                else
                {
                    swap(&maxHeap[i], &maxHeap[r]);
    
                    swap(&iMaxIndices[i], &iMaxIndices[r]);
                    swap(&jMaxIndices[i], &jMaxIndices[r]);
                
                    MaxHeapify(r);
                }
            }
        }
       
    }
    void InsertMaxHeap(integer_t newNum, integer_t i, integer_t j){
    
       
      maxHeap[maxSize] = newNum;
      iMaxIndices[maxSize]=i;
      jMaxIndices[maxSize]=j;
      integer_t pos = maxSize;
      maxSize++;
      while (maxHeap[pos] > maxHeap[(pos - 1) / 2] && pos> 0){
        
        swap(&maxHeap[pos], &maxHeap[(i - 1) / 2]);
        swap(&iMaxIndices[pos], &iMaxIndices[(pos - 1) / 2]);
        swap(&jMaxIndices[pos], &jMaxIndices[(pos - 1) / 2]);
        pos = (pos - 1) / 2;
      }
       
    }
    void MaxPop(){
      // replace first node by last and delete last
      maxHeap[0] = maxHeap[maxSize - 1];
      iMaxIndices[0] = iMaxIndices[maxSize - 1];
      jMaxIndices[0] = jMaxIndices[maxSize - 1];
      maxSize--;
    
      
      MaxHeapify(0);
      
    }

    //Algorithm itself
      int SS_T(int n, integer_t *p, integer_t desired_sum){
 
        int aSize = (n/2)/2;
        int bSize = ((n/2) - (n/2)/2);
        int cSize = (n - n/2)/2;
        int dSize = ((n - n/2) - (n - n/2)/2);
        
      
        int aCombSize = 1<<aSize;
        int bCombSize = 1<<bSize;
        int cCombSize = 1<<cSize;
        int dCombSize = 1<<dSize;
      
        // arranja espaço para as somas
        integer_t *A = malloc(sizeof(integer_t)*aCombSize);
        integer_t *B = malloc(sizeof(integer_t)*bCombSize);
        integer_t *C = malloc(sizeof(integer_t)*cCombSize);
        integer_t *D = malloc(sizeof(integer_t)*dCombSize);
      
        faster_calcsubarray(p, A, aSize, 0);
        faster_calcsubarray(p, B, bSize, aSize);
        faster_calcsubarray(p, C, cSize, aSize + bSize);
        faster_calcsubarray(p, D, dSize, aSize + bSize + cSize);
        
        for (int k = 0; k < bCombSize; k++){
          InsertMinHeap(B[k], k, 0);
        }

        for (int k = 0; k < aCombSize; k++){
          InsertMaxHeap(A[k] + C[cCombSize - 1], k, cCombSize - 1);
        }
         
        integer_t K = pow(2,n);
        integer_t min, max;
      
        for (integer_t i = 0; i < K; i++){
        
          max = maxHeap[0];
          min = minHeap[0];
      
          iMin = iMinIndices[0];
          jMin = jMinIndices[0];
          
          iMax = iMaxIndices[0];
          jMax = jMaxIndices[0];
      
          integer_t sum = min + max;

          
          if (sum == desired_sum){
            return 1;
          }
          else if (sum < desired_sum){
          
            MinPop();
            if (jMin + 1 < dCombSize){ 
               
              InsertMinHeap(B[iMin] + D[jMin+1], iMin, jMin+1);
            }
          }
          else{
          
            MaxPop();
            if (jMax > 0){ 
              InsertMaxHeap(A[iMax] + C[jMax-1], iMax, jMax-1);
            }
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
  for(int i = 0;i < 1;i++){
  
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
    for(int j = 0;j < 1;j++){
           
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
      printf("Schroeppel and Shamir technique       %d,  %lld || %lld \n", j ,sum,z); 
            
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