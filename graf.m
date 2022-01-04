clear;
close all;
clc;

A_1 = load("data_1.log");
A_2 = load("data_1_max.log");
A_3 = load("data_2.log");
A_4 = load("data_2_max.log");
A_5 = load("data_3.log");
A_6 = load("data_3_max.log");
A_7 = load("data_4.log");
A_8 = load("data_4_max.log");
A_9 = load("data_5.log");
A_10 = load("data_5_max.log");
A_11 = load("data_6.log");
A_12 = load("data_6_max.log");
A_13 = load("data_7.log");
A_14 = load("data_7_max.log");
 


 %se coluda 2==1  %coluna
n_V1_1 = A_1(:, 1);
f_V1_1 = A_1(:, 2); 

n_V1_2 = A_2(:, 1);
f_V1_2 = A_2(:, 2); 

n_V1_3 = A_3(:, 1);
f_V1_3 = A_3(:, 2);

n_V1_4 = A_4(:, 1);
f_V1_4 = A_4(:, 2);

n_V1_5 = A_5(:, 1);
f_V1_5 = A_5(:, 2); 

n_V1_6 = A_6(:, 1);
f_V1_6 = A_6(:, 2); 

n_V1_7 = A_7(:, 1);
f_V1_7 = A_7(:, 2); 

n_V1_8 = A_8(:, 1);
f_V1_8 = A_8(:, 2);

n_V1_9 = A_9(:, 1);
f_V1_9 = A_9(:, 2); 

n_V1_10 = A_10(:, 1);
f_V1_10 = A_10(:, 2);


n_V1_11 = A_11(:, 1);
f_V1_11 = A_11(:, 2); 

n_V1_12 = A_12(:, 1);
f_V1_12 = A_12(:, 2);

n_V1_13 = A_13(:, 1);
f_V1_13 = A_13(:, 2);

n_V1_14 = A_14(:, 1);
f_V1_14 = A_14(:, 2);

    
  
%semilogy(n_V1_7, f_V1_7, n_V1_8, f_V1_8, n_V1_7, f_V1_7./10, n_V1_8, f_V1_8./10);
 
%semilogy(n_V1_1, f_V1_1, n_V1_2, f_V1_2, n_V1_3, f_V1_3, n_V1_4, f_V1_4, n_V1_5, f_V1_5, n_V1_6, f_V1_6, n_V1_7, f_V1_7, n_V1_8, f_V1_8, n_V1_9, f_V1_9, n_V1_10, f_V1_10, n_V1_11, f_V1_11, n_V1_12, f_V1_12, n_V1_13, f_V1_13, n_V1_14, f_V1_14);
semilogy(n_V1_1, f_V1_1, n_V1_3, f_V1_3, n_V1_5, f_V1_5, n_V1_7, f_V1_7, n_V1_9, f_V1_9, n_V1_11, f_V1_11, n_V1_13, f_V1_13 );

%semilogy( n_V1_7, f_V1_7, n_V1_8, f_V1_8,n_V1_9, f_V1_9, n_V1_10, f_V1_10, n_V1_11, f_V1_11, n_V1_12, f_V1_12, n_V1_13, f_V1_13, n_V1_14, f_V1_14);

%legend('Meet in the meedle','worst Meet in the meedle','Fast Meet in the meedle','worst Fast Meet in the meedle','Schroeppel and Shamir technique','worst Schroeppel and Shamir technique','Schroeppel and Shamir technique Teles','worst Schroeppel and Shamir technique Teles');
%legend('Iterative', 'worst Iterative','Recursive','worst Recursive','Recursive Smart','worst Recursive Smart','Meet in the meedle','worst Meet in the meedle','Fast Meet in the meedle','worst Fast Meet in the meedle','Schroeppel and Shamir technique','worst Schroeppel and Shamir technique','Schroeppel and Shamir technique Teles','worst Schroeppel and Shamir technique Teles');
legend('Iterative',  'Recursive', 'Recursive Smart', 'Meet in the meedle' ,'Fast Meet in the meedle' ,'Schroeppel and Shamir technique', 'Schroeppel and Shamir technique Teles' );
grid on


%% 
figure(1);    
semilogy(n_V1_1, f_V1_1, n_V1_2, f_V1_2); 
xlim([0,60]);
figure(2);   
semilogy(n_V1_3, f_V1_3, n_V1_4, f_V1_4);
xlim([0,60]);
figure(3);     
semilogy(n_V1_5, f_V1_5,n_V1_6, f_V1_6); 
xlim([0,60]);
figure(4);     
semilogy(n_V1_7, f_V1_7,n_V1_8, f_V1_8); 
xlim([0,60]);
figure(5);     
semilogy(n_V1_9, f_V1_9,n_V1_10, f_V1_10); 
xlim([0,60]);
figure(6);    
semilogy(n_V1_11, f_V1_11, n_V1_12, f_V1_12); 
xlim([0,60]);  


 

