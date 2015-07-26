% This is to find the optimal threshold for decentralized sequential hypothesis 
% testing using fminsearchbnd
tic 
options1=optimset('Display','Iter','MaxFunEvals',inf,'MaxIter',inf, 'TolX', 0.001, 'TolFun', 1);
[bound, fval, exitflag, output]=fminsearchbnd(@asymp_decentrl,[0.01 0.98 0.02 0.98],[0 0 0 0],[1 1 1 1],options1)
toc