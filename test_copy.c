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
 

int Bf_recur( unsigned int n,unsigned int m,integer_t *p,double sum, int comb,integer_t desired_sum)
{   
    
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

int Bf_recur_smart( unsigned int n,int m,integer_t *p,int sum, int comb,integer_t desired_sum)
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

    int stuff = Bf_recur_smart(n,m-1,p,sum + p[m], comb + pow(2,m),desired_sum);  
    if (stuff == 0)  {                 
        return Bf_recur_smart(n,m-1,p,sum  ,comb,desired_sum);      
    }   
    return stuff;
}

 

void merge(integer_t arr[], integer_t a[], int l, int m, int r)
{
    integer_t i, j, k;
    integer_t n1 = m - l + 1;
    integer_t n2 = r - m;
  
    /* create temp arrays */
    integer_t L[n1], R[n2];
    integer_t L_2[n1], R_2[n2];
  
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++){
        L[i] = arr[l + i];
        L_2[i] = a[l + i];}
    for (j = 0; j < n2; j++){
        R[j] = arr[m + 1 + j];
        R_2[j] = a[m + 1 + j];}
  
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            a[k] =  L_2[i];
            i++;
        }
        else {
            arr[k] = R[j];
            a[k] = R_2[i];
            j++;
        }
        k++;
    }
  
    /* Copy the remaining elements of L[], if there
    are any */
    while (i < n1) {
        arr[k] = L[i];
        a[k] = L_2[i];
        i++;
        k++;
    }
  
    /* Copy the remaining elements of R[], if there
    are any */
    while (j < n2) {
        arr[k] = R[j];
        a[k] = R_2[j];
        j++;
        k++;
    }
}
  
 
void mergeSort(integer_t arr[],integer_t a[], int l, int r)
{
    if (l < r) {
        // Same as (l+r)/2, but avoids overflow for
        // large l and h
        integer_t m = l + (r - l) / 2;
  
        // Sort first and second halves
        mergeSort(arr,a, l, m);
        mergeSort(arr,a, m + 1, r);
  
        merge(arr,a, l, m, r);
    }
}










void calcsubarray(integer_t a[], integer_t x[], integer_t b[], int n, int c)
{
    for (integer_t i=0; i<(1<<n); i++)
    {
        integer_t s = 0;
        for (integer_t j=0; j<n; j++){
            if (i & (1<<j)){
                s += a[j+c];
                b[i] += pow(2, j);
            }    
                 
        }
        if(s > 0){
          x[i] = s;  
        }

        b[i] = b[i]<<c;

        

        
    }
} 


// insertion sort q roubei so prar 
void Sort(integer_t arr[], integer_t n)
{
    integer_t i, key, j;
    for (i = 1; i < n; i++) {
        key = arr[i];
        j = i - 1;
  
        while (j >= 0 && arr[j] > key) {
            arr[j + 1] = arr[j];
            j = j - 1;
        }
        arr[j + 1] = key;
    }
}



int mitm(int n, integer_t *p, integer_t desired_sum){
 
    // pega o tamanho
    integer_t size_X = 1<<(n/2);
    integer_t size_Y = 1<<(n-n/2);


    // arranja espaço para as somas
    integer_t *X = malloc(size_X*sizeof(integer_t));
    integer_t *Y = malloc(size_Y*sizeof(integer_t));

    integer_t *a = malloc(size_X*sizeof(integer_t));
    integer_t *b = malloc(size_X*sizeof(integer_t));

    
     

    /*
    
        integer_t a[n/2];
        memcpy(a, p, (n/2) * sizeof(integer_t));
        integer_t b[(n+1)/2];
        memcpy(b, &p[(n/2)], ((n+1)/2) * sizeof(integer_t));
    */

 
    // enche os arrays com as somas respetivas
    calcsubarray(p, X, a, n/2, 0);
    calcsubarray(p, Y, b,  n-n/2, n/2);

   

 
    // Sorta os arrays (Acho que isto pode ser o problema)
    mergeSort(X,a, 0, size_X-1);
    mergeSort(Y,b, 0, size_Y-1);


     
 

   
    integer_t max;

    
    // loopa comparando e tal (como o stor explicou)
    for(integer_t i=0; i< size_X;){
        for(integer_t j=size_Y -1; j>=0;){
            
            if(X[i]+Y[j] == desired_sum ){ 

                printf("\n %lld + %lld = %lld\n", X[i], Y[j], desired_sum);
                return a[i] + b[j];   
            }else if(X[i]+Y[j] < desired_sum){
                i++;
            }else if(X[i]+Y[j] > desired_sum){
                j--;
            } 
           max = X[i]+Y[j];
        }
        return max; //retorna a ultima soma se deu merda
    }

    // liberta o espaco
    free(X);
    free(Y);

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

int main(void)
{           
    /* Setting up file */
    FILE *fp_1 = NULL;
    FILE *fp_2 = NULL;
    FILE *fp_3 = NULL;
    FILE *fp_4 = NULL;
	
	remove("data.log");    
    /* Open for the first time the file provided as argument */
    fp_1 = fopen("data_1.log", "a");
    fp_2 = fopen("data_2.log", "a");
    fp_3 = fopen("data_3.log", "a");
    fp_4 = fopen("data_4.log", "a");
    integer_t pi[5] = {1,3,5,7,8};
    printf("\n\nResposta = %d",Bf_recur_smart(5,4,pi ,0,0,9));
    
    printf("\n");
    printf("Program configuration:\n");
    printf("  min_n ....... %d\n",min_n);
    printf("  max_n ....... %d\n",max_n);
    printf("  n_sums ...... %d\n",n_sums);
    printf("  n_problems .. %d\n",n_problems);
    printf("  integer_t ... %d bits\n",8 * (int)sizeof(integer_t));
    printf("\n");
     

   // start looping for n's
    for(int i = 0;i < n_problems;i++)
    {
        printf("--------------------------- \n");
                
        // get n
        int n = all_subset_sum_problems[i].n;  
        printf("n =  %i\n\n",n);
        
        // get p and sums
        integer_t *p = all_subset_sum_problems[i].p;    
        integer_t *sums = all_subset_sum_problems[i].sums;    

        //create space to store results
        char comb_bin[n+1];

        double dt_bf_i = 0;  
        double dt_bf_r = 0;
        double dt_bf_i_s = 0;
        double dt_mitm = 0;     
        // loop for sums
        for(int j = 0;j < n_sums;j++)
        {   
            integer_t sum = all_subset_sum_problems[i].sums[j];

            // run function and take time 
            double tmp_dt = cpu_time();   
            int comb = Bf_Iter(n, p, sum);
            tmp_dt = cpu_time() - tmp_dt;
            dt_bf_i += tmp_dt;

            tmp_dt = cpu_time();   
            //int comb_rec= Bf_recur(n,0, p,0,0,sum);
            tmp_dt = cpu_time() - tmp_dt;
            dt_bf_r += tmp_dt;

            tmp_dt = cpu_time();   
            //int comb_smart= Bf_recur_smart(n,n-1, p,0,0,sum);
            tmp_dt = cpu_time() - tmp_dt;
            dt_bf_i_s += tmp_dt;

            tmp_dt = cpu_time();   
           
            int x= mitm(n, p, sum);
            tmp_dt = cpu_time() - tmp_dt;
            dt_mitm += tmp_dt;
 
            // print results
            printf("-------------------------------------------------\n");
            printf("Brute force             %d,  %lld || %i -> %s  \n", j ,sum,comb,   Converter(n, comb, comb_bin));
            //printf("Brute force recursiva   %d,  %lld || %i -> %s  \n", j ,sum,comb_rec,   Converter(n, comb_rec, comb_bin));
            //printf("Brute force recur smart %d,  %lld || %i -> %s  \n", j ,sum,comb_smart,   Converter(n, comb_smart, comb_bin));
            //printf("Brute force             %d,  %lld || %i -> %s \n", j ,sum, comb, Converter(n, comb, comb_bin));
            printf("Meet in the middle      %d,  %lld || %i -> %s \n", j ,sum, x, Converter(n, x, comb_bin));
            
        }

        // store times 
        //fprintf(fp_1,"%i %f \n",n, dt_bf_i);
        //fprintf(fp_2,"%i %f \n",n, dt_bf_r);
        //fprintf(fp_3,"%i %f \n",n, dt_bf_i_s);
        //fprintf(fp_4,"%i %f \n",n, dt_mitm);

    }   


    /*
     for(int i = 0;i < 15;i++)
    {
        printf("--------------------------- \n");
 
        // get n
        int n = all_subset_sum_problems[i].n;    
        char comb_bin[n+1];
        printf("n =  %i\n\n",n);
        // get p and sums
        integer_t *p = all_subset_sum_problems[i].p;    
        integer_t *sums = all_subset_sum_problems[i].sums; 
        
        for(int j = 0;j < n_sums;j++){   
            integer_t sum = all_subset_sum_problems[i].sums[j];
            int combination = print_all_sums_recursive(n,0, p,0,0,sum);
            printf("\nsoma = %llu --- solution = %s ", sum ,Converter(n, combination, comb_bin));

        }

    } 
    return 0;

    */
    
}
