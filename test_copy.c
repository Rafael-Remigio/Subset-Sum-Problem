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

integer_t Bf_Iter(int n, integer_t *p, integer_t desired_sum){

    integer_t comb = 0;
    integer_t test_sum;

    for(comb;comb<(1<<n);comb++){
        test_sum =0;

        for(int bit=0; bit<n ;bit++){
            integer_t mask = (1<<bit);
            if((comb & mask)!=0){test_sum += p[bit];}
        }   

        if(test_sum == desired_sum){break;}
    }

    return comb;
    
}
 

integer_t Bf_recur( unsigned int n,unsigned int m,integer_t *p,double sum, integer_t comb,integer_t desired_sum)
{   
    
    if(m == n)
    {
        if (sum == desired_sum){
            return comb;
            
        }
        else return 0;
    }

    
    integer_t stuff = Bf_recur(n,m + 1u,p,sum , comb ,desired_sum);  
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

 

void swap(integer_t *a, integer_t *b) {
    integer_t temp = *a;
    *a = *b;
    *b = temp;
}
  
void heapify(integer_t arr[],integer_t a[], int n, integer_t i) {
    // Find largest among root, left child and right child
    integer_t largest = i;
    integer_t left = 2 * i + 1;
    integer_t right = 2 * i + 2;
  
    if (left < n && arr[left] > arr[largest])
      largest = left;
  
    if (right < n && arr[right] > arr[largest])
      largest = right;
  
    // Swap and continue heapifying if root is not largest
    if (largest != i) {
      swap(&arr[i], &arr[largest]);
      swap(&a[i], &a[largest]);
      heapify(arr, a, n, largest);
    }
}
  
// Main function to do heap sort
void heapSort(integer_t arr[],integer_t a[], int n) {
    // Build max heap
    for (int i = n / 2 - 1; i >= 0; i--)
      heapify(arr, a, n, i);
  
    // Heap sort
    for (int i = n - 1; i >= 0; i--) {
      swap(&arr[0], &arr[i]);
      swap(&a[0], &a[i]);
  
      // Heapify root element to get highest element at root again
      heapify(arr, a, i, 0);
    }
}





void calcsubarray(integer_t a[], integer_t x[], integer_t b[], int n, int c)
{   
    // para cada soma 
    for (integer_t i=0; i<(1<<n); i++)
    {
        integer_t s = 0;
        integer_t cmb = 0;
        // para cada p[]
        for (integer_t j=0; j<n; j++){
            // se i tiver alguma coisa em comum com 2^j
            if (i & (1<<j)){
                s += a[j+c];
                cmb += pow(2, j+c);
            }    
                 
        }
        if(s >= 0){
          x[i] = s;  
          b[i] = cmb;  
        }

        

    }
} 
 

integer_t mitm(int n, integer_t *p, integer_t desired_sum){
 
    // pega o tamanho
    integer_t size_X = 1<<(n/2);
    integer_t size_Y = 1<<(n-n/2);


    // arranja espaço para as somas
    integer_t *X = malloc(size_X*sizeof(integer_t));
    integer_t *Y = malloc(size_Y*sizeof(integer_t));

    integer_t *a = malloc(size_X*sizeof(integer_t));
    integer_t *b = malloc(size_Y*sizeof(integer_t));

    
     

    /*
    
        integer_t a[n/2];
        memcpy(a, p, (n/2) * sizeof(integer_t));
        integer_t b[(n+1)/2];
        memcpy(b, &p[(n/2)], ((n+1)/2) * sizeof(integer_t));
    */

 
    // enche os arrays com as somas respetivas
    calcsubarray(p, X, a, n/2, 0);
    calcsubarray(p, Y, b,  n-n/2, n/2);

   

 
    // Sorta os arrays  
    heapSort(X, a, size_X);
    heapSort(Y, b, size_Y);

    /*
    printf("\n");
    for(int i=0; i<size_X;i++){
        printf("%lld, ", X[i]);
        printf("%lld, ", a[i]);
    }
    printf("\n");
    */

    // loopa comparando e tal (como o stor explicou)
    for(integer_t i=0; i< size_X;){
        for(integer_t j=size_Y -1; j>=0;){
            
            if(X[i]+Y[j] == desired_sum ){ 
                //printf("\n %lld + %lld = %lld\n", X[i], Y[j], desired_sum);
                // liberta o espaco
                free(X);
                free(Y);
                integer_t r = a[i] + b[j];
                free(a);
                free(b);
                return r;   
            }else if(X[i]+Y[j] < desired_sum){
                i++;
            }else if(X[i]+Y[j] > desired_sum){
                j--;
            } 
           
        }
        
    }

    

    return 0;
}



char *Converter(int n,integer_t x, char *sol){
    for(int bit=0; bit<n ;bit++){
        integer_t mask = (1ull<<bit);   
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
    
    printf("\n");
    printf("Program configuration:\n");
    printf("  min_n ....... %d\n",min_n);
    printf("  max_n ....... %d\n",max_n);
    printf("  n_sums ...... %d\n",n_sums);
    printf("  n_problems .. %d\n",n_problems);
    printf("  integer_t ... %d bits\n",8 * (int)sizeof(integer_t));
    printf("\n");
     

   // start looping for n's
    for(int i = 20;i < n_problems;i++)
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
            //integer_t comb = Bf_Iter(n, p, sum);
            tmp_dt = cpu_time() - tmp_dt;
            dt_bf_i += tmp_dt;

            tmp_dt = cpu_time();   
            //integer_t comb_rec= Bf_recur(n,0, p,0,0,sum);
            tmp_dt = cpu_time() - tmp_dt;
            dt_bf_r += tmp_dt;

            tmp_dt = cpu_time();   
            integer_t comb_smart= Bf_recur_smart(n,n-1, p,0,0,sum);
            tmp_dt = cpu_time() - tmp_dt;
            dt_bf_i_s += tmp_dt;

            tmp_dt = cpu_time();   
           
            integer_t x= mitm(n, p, sum);
            tmp_dt = cpu_time() - tmp_dt;
            dt_mitm += tmp_dt;
 
            // print results
            printf("-------------------------------------------------\n");
            //printf("Brute force             %d,  %lld || %lld -> %s  \n", j ,sum,comb,   Converter(n, comb, comb_bin));
            //printf("Brute force recursiva   %d,  %lld || %lld -> %s  \n", j ,sum,comb_rec,   Converter(n, comb_rec, comb_bin));
            printf("Brute force recur smart %d,  %lld || %llu -> %s  \n", j ,sum,comb_smart,   Converter(n, comb_smart, comb_bin)); 
            printf("Meet in the middle      %d,  %lld || %llu -> %s \n", j ,sum, x, Converter(n, x, comb_bin));
            
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
