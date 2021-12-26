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
     
    // Sorta os arrays (Acho que isto pode ser o problema)
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

    // enche os arrays com as somas respetivas
    faster_calcsubarray(p, X, n/2, 0);
    faster_calcsubarray(p, Y, n-n/2, n/2); 
     
    printf("\n");
    printf("\n");
    for(int i=0;i<size_X;i++){
        printf("%lld,", X[i]);
    } 
    printf("\n------------\n");
    for(int i=0;i<size_Y;i++){
        printf("%lld,", Y[i]);
    }
    printf("\n");
    printf("\n");

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
    
 
    
    printf("\n");
    printf("Program configuration:\n");
    printf("  min_n ....... %d\n",min_n);
    printf("  max_n ....... %d\n",max_n);
    printf("  n_sums ...... %d\n",n_sums);
    printf("  n_problems .. %d\n",n_problems);
    printf("  integer_t ... %d bits\n",8 * (int)sizeof(integer_t));
    printf("\n");
     

   // start looping for n's
    for(int i = 0;i < 1;i++)
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
            

            // Meet in the middle
            tmp_dt = cpu_time();   
            int x= mitm(n, p, sum);
            tmp_dt = cpu_time() - tmp_dt;

            if(tmp_dt > dt_mitm_max){
                dt_mitm_max = tmp_dt;
            }
            dt_mitm += tmp_dt;
            */


            int y= faster_mitm(n, p, sum);
 
            // print results
            //printf("-------------------------------------------------\n");
            //printf("Brute force             %d,  %lld || %i -> %s  \n", j ,sum,comb,   Converter(n, comb, comb_bin));
            //printf("Brute force recursiva   %d,  %lld || %i -> %s  \n", j ,sum,comb_rec,   Converter(n, comb_rec, comb_bin));
            //printf("Brute force recur smart %d,  %lld || %lld -> %s  \n", j ,sum,comb_smart,   Converter(n, comb_smart, comb_bin));
            //printf("Brute force             %d,  %lld || %i -> %s \n", j ,sum, comb, Converter(n, comb, comb_bin));
            //printf("Meet in the middle      %d,  %lld || %i  \n", j ,sum, x);
            printf("Faster meet in the middle      %d,  %lld || %i  \n", j ,sum, y);
            
        }

        // store times 
        /*
        fprintf(fp_1,"%i %f \n",n, dt_bf_i/20);
        fprintf(fp_2,"%i %f \n",n, dt_bf_i_max);
        fprintf(fp_3,"%i %f \n",n, dt_bf_r/20);
        fprintf(fp_4,"%i %f \n",n, dt_bf_r_max);
        fprintf(fp_5,"%i %f \n",n, dt_bf_i_s/20);
        fprintf(fp_6,"%i %f \n",n, dt_bf_i_s_max);
        fprintf(fp_7,"%i %f \n",n, dt_mitm/20);
        fprintf(fp_8,"%i %f \n",n, dt_mitm_max);*/

    }   
    
}
