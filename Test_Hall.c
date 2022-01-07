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
  /*
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
  */
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
    // é
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
    //
  //

  // Schroeppel and Shamir technique by Teles
     
    void swapInt(int *a, int *b) {
      integer_t temp = *a;
      *a = *b;
      *b = temp;
    } 

    int iMinIndices[1<<17];
    int jMinIndices[1<<17];
    int iMaxIndices[1<<17];
    int jMaxIndices[1<<17];
    
    // Halll System
      typedef struct Heap Heap;
      struct Heap{

        // type = 0 (Heap) || type = 1 (MaxHeap)
        char type;
        //
        integer_t *arr;
        // Current Size of the Heap
        integer_t size;
        // Maximum capacity of the heap
        integer_t capacity;
      };

      integer_t parent(integer_t i){

        // Get the index of the parent
        return (i - 1) / 2;
      }

      integer_t left_child(integer_t i){

        return (2 * i + 1);
      }

      integer_t right_child(integer_t i){

        return (2 * i + 2);
      }

      integer_t get_min(Heap *heap){

        // Only for MinHeaps
        // Return the root node element,
        // since that's the minimum
        return heap->arr[0];
      }

      integer_t get_max(Heap *heap){

        // Only for MaxHeaps
        // Return the root node element,
        // since that's the maximum
        return heap->arr[0];
      }

      Heap *init_heap(integer_t capacity){

        Heap *heap = (Heap *)calloc(1ull, sizeof(Heap));
        heap->arr = (integer_t *)calloc(capacity, sizeof(integer_t));
        heap->capacity = capacity;
        heap->size = 0;
        return heap;
      }

      Heap *insert_heap(Heap *heap, integer_t element, char type, int i, int j){

        if (type == 0){
        
            // Inserts an element to the min heap
            // We first add it to the bottom (last level)
            // of the tree, and keep swapping with it's parent
            // if it is lesser than it. We keep doing that until
            // we reach the root node. So, we will have inserted the
            // element in it's proper position to preserve the min heap property
            if (heap->size == heap->capacity)
            {
                fprintf(stderr, "Cannot insert %llu. Heap is already full!\n", element);
                return heap;
            }
 
            // We can add it. Increase the size and add it to the end
            heap->size++;
            heap->arr[heap->size - 1ull] = element;
            iMinIndices[heap->size]=i;
            jMinIndices[heap->size]=j;

            // Keep swapping until we reach the root
            integer_t curr = heap->size - 1ull;
            // As long as you aren't in the root node, and while the
            // parent of the last element is greater than it
            while (curr > 0 && heap->arr[parent(curr)] > heap->arr[curr])
            {
                // Swap
                integer_t temp = heap->arr[parent(curr)];
                heap->arr[parent(curr)] = heap->arr[curr];
                heap->arr[curr] = temp;
                // Update the current index of element
                curr = parent(curr);
            }
            return heap;
        }
        else if (type == 1){
        
            // Inserts an element to the max heap
            // We first add it to the bottom (last level)
            // of the tree, and keep swapping with it's parent
            // if it is grater than it. We keep doing that until
            // we reach the root node. So, we will have inserted the
            // element in it's proper position to preserve the max heap property
            if (heap->size == heap->capacity)
            {
                fprintf(stderr, "Cannot insert %llu. Heap is already full!\n", element);
                return heap;
            }
           

            // We can add it. Increase the size and add it to the end
            heap->size++;
            heap->arr[heap->size - 1ull] = element;
            iMaxIndices[heap->size]=i;
            jMaxIndices[heap->size]=j;

            // Keep swapping until we reach the root
            integer_t curr = heap->size - 1ull;
            // As long as you aren't in the root node, and while the
            // parent of the last element is lesser than it
            while (curr > 0 && heap->arr[parent(curr)] < heap->arr[curr])
            {
                // Swap
                integer_t temp = heap->arr[parent(curr)];
                heap->arr[parent(curr)] = heap->arr[curr];
                heap->arr[curr] = temp;
                // Update the current index of element
                curr = parent(curr);
            }
            return heap;
        }
        return heap;
      }

      Heap *heapify(Heap *heap, integer_t index, char type){

        switch (type)
        {
        case 0:
            // Rearranges the heap as to maintain
            // the min-heap property
            if (heap->size <= 1ull)
                return heap;

            integer_t left = left_child(index);
            integer_t right = right_child(index);

            // Variable to get the smallest element of the subtree
            // of an element an index
            integer_t smallest = index;

            // If the left child is smaller than this element, it is
            // the smallest
            if (left < heap->size && heap->arr[left] < heap->arr[index])
                smallest = left;

            // Similarly for the right, but we are updating the smallest element
            // so that it will definitely give the least element of the subtree
            if (right < heap->size && heap->arr[right] < heap->arr[smallest])
                smallest = right;

            // Now if the current element is not the smallest,
            // swap with the current element. The min heap property
            // is now satisfied for this subtree. We now need to
            // recursively keep doing this until we reach the root node,
            // the point at which there will be no change!
            if (smallest != index)
            {
                integer_t temp = heap->arr[index];
                heap->arr[index] = heap->arr[smallest];
                heap->arr[smallest] = temp;

                swapInt(&iMinIndices[index], &iMinIndices[smallest]);
                swapInt(&jMinIndices[index], &jMinIndices[smallest]);

                heap = heapify(heap, smallest, 0);
            }

            return heap;
            break;
        case 1:
            // Rearranges the heap as to maintain
            // the max-heap property
            if (heap->size <= 1ull)
                return heap;

            left = left_child(index);
            right = right_child(index);

            // Variable to get the greatest element of the subtree
            // of an element an index
            integer_t greatest = index;

            // If the left child is greatest than this element, it is
            // the greatest
            if (left < heap->size && heap->arr[left] > heap->arr[index])
                greatest = left;

            // Similarly for the right, but we are updating the greatest element
            // so that it will definitely give the least element of the subtree
            if (right < heap->size && heap->arr[right] > heap->arr[greatest])
                greatest = right;

            // Now if the current element is not the greatest,
            // swap with the current element. The max heap property
            // is now satisfied for this subtree. We now need to
            // recursively keep doing this until we reach the root node,
            // the point at which there will be no change!
            if (greatest != index)
            {
                integer_t temp = heap->arr[index];
                heap->arr[index] = heap->arr[greatest];
                heap->arr[greatest] = temp;

                swapInt(&iMaxIndices[index], &iMaxIndices[greatest]);
                swapInt(&jMaxIndices[index], &jMaxIndices[greatest]);

                heap = heapify(heap, greatest, 1);
            }

            return heap;
            break;
        }
        return heap;
      }     

      Heap *delete_minimum(Heap *heap){

        // Deletes the minimum element, at the root
        if (!heap || heap->size == 0)
            return heap;

        integer_t size = heap->size;
        integer_t last_element = heap->arr[size - 1];

        // Update root value with the last element
        heap->arr[0] = last_element;
        swapInt(&iMinIndices[0], &iMinIndices[size - 1]);
        swapInt(&jMinIndices[0], &jMinIndices[size - 1]);

        // Now remove the last element, by decreasing the size
        heap->size--;
        size--;

        // We need to call heapify(), to maintain the min-heap
        // property
        heap = heapify(heap, 0, 0);
        return heap;
      }

      Heap *delete_maximum(Heap *heap){

        // Deletes the maximum element, at the root
        if (!heap || heap->size == 0)
            return heap;

        integer_t size = heap->size;
        integer_t last_element = heap->arr[size - 1];

        // Update root value with the last element
        heap->arr[0] = last_element;
        swapInt(&iMaxIndices[0], &iMaxIndices[size - 1]);
        swapInt(&jMaxIndices[0], &jMaxIndices[size - 1]);

        // Now remove the last element, by decreasing the size
        heap->size--;
        size--;

        // We need to call heapify(), to maintain the max-heap
        // property
        heap = heapify(heap, 0, 1);
        return heap;
      }

      Heap *delete_element(Heap *heap, integer_t index, char type){

        switch (type)
        {
        case 0:
            // Deletes an element, indexed by index
            // Ensure that it's lesser than the current root
            heap->arr[index] = get_min(heap) - 1;

            // Now keep swapping, until we update the tree
            integer_t curr = index;
            while (curr > 0 && heap->arr[parent(curr)] > heap->arr[curr])
            {
                integer_t temp = heap->arr[parent(curr)];
                heap->arr[parent(curr)] = heap->arr[curr];
                heap->arr[curr] = temp;
                curr = parent(curr);
            }

            // Now simply delete the minimum element
            heap = delete_minimum(heap);
            return heap;
            break;
        case 1:
            // Deletes an element, indexed by index
            // Ensure that it's lesser than the current root
            heap->arr[index] = get_min(heap) - 1;

            // Now keep swapping, until we update the tree
            curr = index;
            while (curr > 0 && heap->arr[parent(curr)] < heap->arr[curr])
            {
                integer_t temp = heap->arr[parent(curr)];
                heap->arr[parent(curr)] = heap->arr[curr];
                heap->arr[curr] = temp;
                curr = parent(curr);
            }

            // Now simply delete the maximum element
            heap = delete_maximum(heap);
            return heap;
            break;
            break;
        }
        return heap;
      }
    
    //

    // Setting Heaps  
       
    //

    //Algorithm itself
      int SS_T(int n, integer_t *p, integer_t desired_sum){
        
        integer_t capacity = 1<<17;

        Heap *MinHeap = (Heap *)calloc(1ull, sizeof(Heap));
        MinHeap->arr = (integer_t *)calloc(capacity, sizeof(integer_t));
        MinHeap->capacity = capacity;
        MinHeap->size = 1;
                
        Heap *MaxHeap = (Heap *)calloc(1ull, sizeof(Heap));
        MaxHeap->arr = (integer_t *)calloc(capacity, sizeof(integer_t));
        MaxHeap->capacity = capacity;
        MaxHeap->size = 1; 


 
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
            
        for (int k = 0; k < aCombSize; k++){
          insert_heap(MinHeap, A[k], 0, k, 0); 
        }

        for (int k = 0; k < (cCombSize); k++){
          insert_heap(MaxHeap, C[k] + D[dCombSize - 1], 1, k, dCombSize - 1);  
        }
         
        integer_t K = pow(2,n);
        integer_t min, max;
      
        for (integer_t i = 0; i < K; i++){
        
          max = MaxHeap->arr[0];
          min = MinHeap->arr[0];
      
          int iMin = iMinIndices[0]; int jMin = jMinIndices[0];
          int iMax = iMaxIndices[0]; int jMax = jMaxIndices[0];
      
          integer_t sum = min + max;

          printf("%lld + %lld = %lld \n", max, min, sum);

          printf("----------\n");
          for(int i=0;i<MaxHeap->size;i++){
            printf("%lld, ", MaxHeap->arr[i]);
          }
          printf("\n-------\n");

          if (sum == desired_sum)
            return 1;
          else if (sum < desired_sum){
          
            printf("<\n");
            delete_minimum(MinHeap);
            if (jMin + 1 < bCombSize){ 
              insert_heap(MinHeap, A[iMin] + B[jMin + 1], 0, iMin, jMin + 1);  
            }
          }
          else{
             
            printf(">\n");
            delete_maximum(MaxHeap);
            if (jMax > 0){ 
              insert_heap(MaxHeap, C[iMax] + D[jMax - 1], 1, iMax, jMax - 1); 
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
  
  
 

  // Loop for n's
  for(int i = 0;i < 1;i++){
  
    printf("--------------------------- \n");
                
    // get n and p
    int n = all_subset_sum_problems[i].n;  
    printf("n =  %i\n\n",n);
    integer_t *p = all_subset_sum_problems[i].p;     

   

 

    // Loop for sum's
    for(int j = 0;j < 1;j++){
           
      integer_t sum = all_subset_sum_problems[i].sums[j];
     

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
 
      int z= SS_T(n, p, sum);
       
      
 
      // print results
      printf("-------------------------------------------------\n");/*
      printf("Brute force                           %d,  %lld || %i -> %s  \n", j ,sum,comb,   Converter(n, comb, comb_bin));
      printf("Brute force recursiva                 %d,  %lld || %i -> %s  \n", j ,sum,comb_rec,   Converter(n, comb_rec, comb_bin));
      printf("Brute force recur smart               %d,  %lld || %lld -> %s  \n", j ,sum,comb_smart,   Converter(n, comb_smart, comb_bin));  
      printf("Meet in the middle                    %d,  %lld || %i  \n", j ,sum, x);
      printf("Faster meet in the middle             %d,  %lld || %i  \n", j ,sum, y);*/
      printf("Schroeppel and Shamir technique       %d,  %lld || %d \n", j ,sum,z); 
            
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