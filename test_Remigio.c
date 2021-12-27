#include <stdio.h> 
#include <stdlib.h>

#include "elapsed_time.h"
#include "000000.h"
#include <string.h>
#include <math.h>
#define min_n       10
#define max_n       57
#define n_sums      20
#define n_problems  (max_n - min_n + 1)

int Bf_Iter(int n, integer_t *p, integer_t desired_sum){

    int comb = 0;
    integer_t test_sum;

    for(comb;comb<(1<<n);comb++){
        test_sum =0;

        for(int bit=0; bit<n ;bit++){
            int mask = (1<<bit);
            if((comb & mask)!=0){test_sum += p[bit];}
        }   

        if(test_sum == desired_sum){break;}
    }

    return comb;
    
}
 
int Bf_recur( unsigned int n,unsigned int m,integer_t *p,double sum, int comb,integer_t desired_sum){
    
    if(m == n)
    {
        if (sum == desired_sum){
            return comb;
            
        }
        else return 0;
    }

    
    int stuff = Bf_recur(n,m + 1u,p,sum , comb ,desired_sum);  
    if (stuff == 0)  {                 
 
        return Bf_recur(n,m + 1u,p,sum + p[m] ,comb+ pow(2,m),desired_sum);      
    }   
    return stuff;
}

integer_t Bf_recur_smart( unsigned int n,int m,integer_t *p,integer_t sum, integer_t comb,integer_t desired_sum)
{  

    if (sum == desired_sum){
        return comb;
    }
    if(sum > desired_sum){
        return 0;
    }
    if(m == -1)
    { 
        
        return 0 ;
    }
    integer_t power = (1ull);
    integer_t stuff = Bf_recur_smart(n,m-1,p,sum + p[m], comb + (power<<m),desired_sum);  
    if (stuff == 0)  {                 
        return Bf_recur_smart(n,m-1,p,sum  ,comb,desired_sum);      
    }   
    return stuff;
}

  

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

void swap2(int *a[2], int *b[2]) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

 
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
  
// Main function to do heap sort
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

void bubbleSort(integer_t arr[], int n)
{
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



int mitm(int n, integer_t *p, integer_t desired_sum){
 
    // pega o tamanho
    int size_X = 1<<(n/2);
    int size_Y = 1<<(n-n/2);
    printf("size x = %i and size y = %i  ",size_X,size_Y);
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


void min_heapify(int (*arr)[2], int n, integer_t A[] , integer_t B[] ,int i)
{
    int smallest = i; // Initialize largest as root
    int l = 2 * i + 1; // down left = 2*i + 1
    int r = 2 * i + 2; // down right = 2*i + 2
 
    // If left child is larger than root
    if (l < n && A[arr[l][0]] + B[arr[l][1]]  <= A[arr[smallest][0]] + B[arr[smallest][1]])
        smallest = l;
 
    // If right child is larger than largest so far
    if (r < n && A[arr[r][0]] + B[arr[r][1]]  <= A[arr[smallest][0]] + B[arr[smallest][1]])
        smallest = r;
 
    // If largest is not root
    if (smallest != i) {
        swap2(&arr[i],&arr[smallest]);
 
        // Recursively heapify the affected sub-tree
        min_heapify(arr, n, A, B ,smallest);
    }

}
void generateMinHeap(int (*minheap)[2] , integer_t A[],integer_t B[],int size_a,int size_b){

        // Generate tree from array indexes
        int heap_iter = 0;
        for(int i = 0; i < size_a; i++){
            for (int j = 0; j < size_b; j++){
                minheap[heap_iter][0] = i;
                minheap[heap_iter][1] = j;
                heap_iter+=1;
            }
        }

        for (int i = 0;i< size_a*size_b;i++){
            printf("\nminheap[%i] = [%i,%i]\t=  %llu",i,minheap[i][0],minheap[i][1], A[minheap[i][0]] + B[minheap[i][1]]);
        }
        
        
        // Perform reverse level order traversal
        // from last non-leaf node and heapify
        int startIdx = ( (size_a*size_b) / 2) - 1;
        printf("\n---%i--------------",startIdx);
        // each node
        for (int i = startIdx; i >= 0; i--) {
            min_heapify(minheap, size_a*size_b , A, B ,i);
        }
        for (int i = 0;i< size_a*size_b;i++){
            printf("\nminheap[%i] = [%i,%i]\t=  %llu",i,minheap[i][0],minheap[i][1], A[minheap[i][0]] + B[minheap[i][1]]);
        }
}

int SS(int n, integer_t *p, integer_t desired_sum){

    // pega o tamanho
    int a = (n/2)/2;
    int b = ((n/2) - (n/2)/2);
    int c = (n - n/2)/2;
    int d = ((n - n/2) - (n - n/2)/2);

    // arranja espaço para as somas
        integer_t *A = malloc(sizeof(integer_t)*(1<<a));
        integer_t *B = malloc(sizeof(integer_t)*(1<<b));
        integer_t *C = malloc(sizeof(integer_t)*(1<<c));
        integer_t *D = malloc(sizeof(integer_t)*(1<<d));

        calcsubarray(p, A, a, 0);
        calcsubarray(p, B, b, a);
        calcsubarray(p, C, c, a + b);
        calcsubarray(p, D, d, a + b + c);

        heapSort(A,(1<<a));
        heapSort(B,(1<<b));
        heapSort(C,(1<<c));
        heapSort(D,(1<<d));

        int (*minheap)[2] = malloc(sizeof(int)* 2 * ((1<<a) * (1<<b)));

        generateMinHeap(minheap, A, B,(1<<a),(1<<b));

        for(int i = 0;i < (1<<a);i++){
            for(int j = 0;j < (1<<b);j++){
                for(int h = 0;h < (1<<c);h++){
                    for(int g = 0;g < (1<<d);g++){
                        integer_t soma  = A[i] + B[j] + C[h]+D[g];
                        if (desired_sum==soma){
                            printf("\nsum is %llu \n",soma);
                            printf("\nindexes are : %i -- %i -- %i -- %i\n",i,j,h,g);
                                    free(A);
                                    free(B);
                                    free(C);
                                    free(D);
                                    free(minheap);
                            return 1;
                        }
                    }   
                }
            }
        }
        
        
        printf("\n not found \n");
        free(A);
        free(B);
        free(C);
        free(D);
        free(minheap);
        return 0;
}



char *Converter(int n,int x, char *sol){
    for(int bit=0; bit<n ;bit++){
        int mask = (1<<bit);   
        if((x & mask)!=0){ sol[bit]='1';}else{ sol[bit]='0';}
    }
    sol[n]='\0';

    return sol;
}

int main(void){


 
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
	
	remove("data.log");    
    /* Open for the first time the file provided as argument */
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
    
 
    
    printf("\n");
    printf("Program configuration:\n");
    printf("  min_n ....... %d\n",min_n);
    printf("  max_n ....... %d\n",max_n);
    printf("  n_sums ...... %d\n",n_sums);
    printf("  n_problems .. %d\n",n_problems);
    printf("  integer_t ... %d bits\n",8 * (int)sizeof(integer_t));
    printf("\n");
     

   // start looping for n's
    for(int i = 0;i < 2;i++)
    {
        printf("--------------------------- \n");
                
        // get n
        int n = all_subset_sum_problems[i].n;  
        printf("n =  %i\n\n",n);
        
        // get p and sums
        integer_t *p = all_subset_sum_problems[i].p;     

        char comb_bin[n+1];

        
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

        // loop for sums
        for(int j = 0;j < n_sums;j++)
        {   
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
            */
            
            // Meet in the middle
/*             tmp_dt = cpu_time();   
            int x= mitm(n, p, sum);
            tmp_dt = cpu_time() - tmp_dt;

            if(tmp_dt > dt_mitm_max){
                dt_mitm_max = tmp_dt;
            }
            dt_mitm += tmp_dt;
 */
            int y = SS(n, p, sum);
            /*

            tmp_dt = cpu_time();   
            int y= faster_mitm(n, p, sum);
            tmp_dt = cpu_time() - tmp_dt;

            if(tmp_dt > dt_f_mitm_max){
                dt_f_mitm_max = tmp_dt;
            }
            dt_f_mitm += tmp_dt;*/
 
            // print results
            //printf("-------------------------------------------------\n");
            //printf("Brute force             %d,  %lld || %i -> %s  \n", j ,sum,comb,   Converter(n, comb, comb_bin));
            //printf("Brute force recursiva   %d,  %lld || %i -> %s  \n", j ,sum,comb_rec,   Converter(n, comb_rec, comb_bin));
            //printf("Brute force recur smart %d,  %lld || %lld -> %s  \n", j ,sum,comb_smart,   Converter(n, comb_smart, comb_bin));
            //printf("Brute force             %d,  %lld || %i -> %s \n", j ,sum, comb, Converter(n, comb, comb_bin));
/*             printf("Meet in the middle      %d,  %lld || %i  \n", j ,sum, x); */
            //printf("Faster meet in the middle      %d,  %lld || %i  \n", j ,sum, y);
            
        }

        // store times 
        /*
        fprintf(fp_1,"%i %f \n",n, dt_bf_i/20);
        fprintf(fp_2,"%i %f \n",n, dt_bf_i_max);
        fprintf(fp_3,"%i %f \n",n, dt_bf_r/20);
        fprintf(fp_4,"%i %f \n",n, dt_bf_r_max);
        fprintf(fp_5,"%i %f \n",n, dt_bf_i_s/20);
        fprintf(fp_6,"%i %f \n",n, dt_bf_i_s_max);*/
/*         fprintf(fp_7,"%i %f \n",n, dt_mitm/20);
        fprintf(fp_8,"%i %f \n",n, dt_mitm_max); *//*
        fprintf(fp_9,"%i %f \n",n, dt_f_mitm/20);
        fprintf(fp_10,"%i %f \n",n, dt_f_mitm_max);*/

    }   
    
}
