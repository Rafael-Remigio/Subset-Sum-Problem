import math
import pathlib
import os.path

file1 = open('\Trabalho-AED\Solutions_Teless', 'r')
Lines = file1.readlines()
 

# Strips the newline character
for line in Lines:
    Linha = line.split(" ")
    if Linha[0] != "-------------------------------------------------\n":
        print(Linha[-1])
