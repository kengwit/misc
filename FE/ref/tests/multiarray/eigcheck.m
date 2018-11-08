clear all;
close all;
clc;

format long;

A = [0 1 2;
     1 2 3;
     2 3 4];
[V,D]=eig(A);

Arecon = D(1,1)*V(:,1)*V(:,1)'+D(2,2)*V(:,2)*V(:,2)'+D(3,3)*V(:,3)*V(:,3)';

V
D
Arecon

