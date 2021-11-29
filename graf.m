clear;
clc;

A = load("data.log");


 %se coluda 2==1  %coluna
n_V1 = A(:, 1);
f_V1 = A(:, 2); 


plot(n_V1, f_V1);
xlim([0,60]);
ylim([0,60]);