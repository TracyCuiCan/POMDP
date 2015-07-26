% This is to find the optimal threshold for centralized sequential hypothesis 
% testing using fminsearchbnd
tic 
options=optimset('Display','iter','MaxFunEvals',inf,'MaxIter',inf, 'TolX', 0.001, 'TolFun', 0.001);
[bound, fval, exitflag, output]=fminsearchbnd(@asymp_centrl,[0.001 0.999],[0 0],[1 1],options)
toc