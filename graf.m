clear;
close all;
clc;

A_1 = load("data_1.log");
A_2 = load("data_2.log");
A_3 = load("data_3.log");
A_4 = load("data_4.log");
A_5 = load("data_4_max.log");


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


semilogy(n_V1_1, f_V1_1, n_V1_2, f_V1_2, n_V1_3, f_V1_3, n_V1_4, f_V1_4, n_V1_5, f_V1_5);
grid on
xlim([0,60]); 