#include <stdio.h> 
#include "elapsed_time.h"
#include "000000.h"

#define min_n       10
#define max_n       57
#define n_sums      20
#define n_problems  (max_n - min_n + 1)

int Bf_Iter(int n, integer_t *p, integer_t *sums, integer_t desired_sum){
    
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
    FILE *fp = NULL;
	
	remove("data.log");    
    /* Open for the first time the file provided as argument */
    fp = fopen("data.log", "a");


    
    printf("\n");
    printf("Program configuration:\n");
    printf("  min_n ....... %d\n",min_n);
    printf("  max_n ....... %d\n",max_n);
    printf("  n_sums ...... %d\n",n_sums);
    printf("  n_problems .. %d\n",n_problems);
    printf("  integer_t ... %d bits\n",8 * (int)sizeof(integer_t));
    printf("\n");
     

    // start looping for n's
    for(int i = 0;i < 15;i++)
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

        
        double dt = 0;   
        // loop for sums
        for(int j = 0;j < n_sums;j++)
        {   
            integer_t sum = all_subset_sum_problems[i].sums[j];

            // run function and take time 
            double tmp_dt = cpu_time();   
            int comb = Bf_Iter(n, p, sums, sum);
            tmp_dt = cpu_time() - tmp_dt;
            dt += tmp_dt;

            // print results
            printf("%d,  %lld -> %s\n", j ,sum, Converter(n, comb, comb_bin));
        }

        // store times 
        fprintf(fp,"%i %f \n",n, dt);

    }

    return 0;

    
    
}
