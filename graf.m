clear;
clc;

A_1 = load("data_1.log");
A_2 = load("data_2.log");
A_3 = load("data_3.log");


 %se coluda 2==1  %coluna
n_V1_1 = A_1(:, 1);
f_V1_1 = A_1(:, 2); 

n_V1_2 = A_2(:, 1);
f_V1_2 = A_2(:, 2); 

n_V1_3 = A_3(:, 1);
f_V1_3 = A_3(:, 2); 


plot(n_V1_1, log(f_V1_1), n_V1_2, log(f_V1_2), n_V1_3, log(f_V1_3));

xlim([0,60]);
ylim([-23,9]);