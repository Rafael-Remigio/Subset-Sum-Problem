#include <stdio.h> 
#include "000000.h"

#define min_n       10
#define max_n       57
#define n_sums      20
#define n_problems  (max_n - min_n + 1)

char *BruteForce(int n, integer_t *p, integer_t *sums, integer_t desired_sum, char *comb_bin){
    
    int comb = 0;
    integer_t test_sum;

    for(comb=0;comb<(1<<n);comb++){

        test_sum =0;

        for(int bit=0; bit<n ;bit++){
            int mask = (1<<bit);
            if((comb & mask)!=0){test_sum += p[bit];}
        }   

        if(test_sum == desired_sum){break;}

    }

    
    for(int bit=0; bit<n ;bit++){
        int mask = (1<<bit);   
        if((comb & mask)!=0){ comb_bin[bit]='1';}else{ comb_bin[bit]='0';}


    }
    comb_bin[n]='\0';

    return comb_bin;
}


int main(void)
{   
    
    printf("\n");
    printf("Program configuration:\n");
    printf("  min_n ....... %d\n",min_n);
    printf("  max_n ....... %d\n",max_n);
    printf("  n_sums ...... %d\n",n_sums);
    printf("  n_problems .. %d\n",n_problems);
    printf("  integer_t ... %d bits\n",8 * (int)sizeof(integer_t));
    printf("\n");
    //
    // for each n
    //
    for(int i = 0;i < 1;i++)
    {
        int n = all_subset_sum_problems[i].n; // the value of n
        integer_t *p = all_subset_sum_problems[i].p;    
        integer_t *sums = all_subset_sum_problems[i].sums;    
        char comb_bin[n];

        for(int j = 0;j < n_sums;j++)
        {
            printf("\n");
            printf("%lld -> %s\n",all_subset_sum_problems[i].sums[j], BruteForce(n, p, sums, all_subset_sum_problems[i].sums[j], comb_bin));
        }
    }

    return 0;

    
    
}
