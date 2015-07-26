% This is meant to test the function intersection£¬ which is to find the
% intersection of two lines. The input is four points, each two defining
% one line. The output is the intersection of these two lines.
clear;
clc;
A=[0 0];
B=[1 1];
C=[0 1];
D=[1 2];
[X,Y]=intersec(A,B,C,D) % two parallel lines should have no intersection, but in this special case, we set it to be the original point
A=[0 1];
B=[1 1];
C=[0 2];
D=[2 0];
[X,Y]=intersec(A,B,C,D) % the result should be [1,1]
A=[0 -2];
B=[-2 0];
C=[0 0];
D=[-1 -1];
[X,Y]=intersec(A,B,C,D) % the result should be [-1,-1]
