% This is meant to test the function minksum, which is to take minkowski
% summation of two sets. The input is two sets, output is the unique minkowski sum
% of these two sets.
clear;
clc;
A=[1 2;3 4];
B=[2 4;5 3];
C=minksum(A,B) % the result should be [3 6;6 5;5 8;8 7]
A=[1 2;3 4];
B=[3 4;1 2];
C=minksum(A,B) % the result should be [2 4;4 6;6 8]