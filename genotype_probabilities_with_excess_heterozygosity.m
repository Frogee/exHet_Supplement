% MATLAB M-file to derive the general solution for the probability of 
%    a marker in a genotype class given a selfing population, F_{t}, for t
%    generations. This is supplemental information for Truong & McCormick et 
%    al (2014) where we incorporate a heterozygosity zygotic viability term.
%    See "Resolution of genetic map expansion caused by excess 
%    heterozygosity in plant recombinant inbred populations" for more 
%    information.
%
% Written by Sandra Truong at Texas A&M University, July 2014.
% Provided as is, without warranty, without guarantee.
% *     This program is free software; you can redistribute it and/or
% *     modify it under the terms of the GNU General Public License,
% *     version 3, as published by the Free Software Foundation.
% *
% *     This program is distributed in the hope that it will be useful,
% *     but without any warranty; without even the implied warranty of
% *     merchantability or fitness for a particular purpose.  See the GNU
% *     General Public License, version 3, for more details.
% *
% *     A copy of the GNU General Public License, version 3, is available
% *     at http://www.r-project.org/Licenses/GPL-3
% ##############################################################################

% assign variables:
% r is the recombination frequency, and t is the generation interval
syms r t;
%
% h is amount of heterozygosity maintained in each generation and can be
% parameterized given generation t. That is if H is the amount of
% heterozygosity in an F_{t} population, then h^(t-1)=H
syms h;
%
% u is the viability of Aa to AA and aa 
syms u;
u = (2*h)/(2-2*h);
%
% d is a parameter necessary to weigh to u appropriately
syms d;
d = 2*((1-r)^2)+8*u*r*(1-r)+ 2*(r^2)+ 2*(u^2)*(((1-r)^2)+(r^2));

% Transition probability matrix for 5 classes of genotypes
T = [
1,0,(1-h^1)/2, (2*((1-r^1)^2))/d^1, (2*(r^2))/d^1;
0,1,(1-h^1)/2, (2*(r^2))/d^1, (2*((1-r^1)^2))/d^1;
0,0, (h^1), (8*u^1*r^1*(1-r^1)^1)/d^1, (8*u^1*r^1*(1-r^1)^1)/d^1;
0,0,0, (2*u^2*(1-r^1)^2)/d^1, (2*u^2*r^2)/d^1;
0,0,0, (2*u^2*r^2)/d^1, (2*u^2*(1-r^1)^2)/d^1];

% Take eigenvalues of Transition probability matrix to set up system of
% equations to find the general solution given generation t for all 5
% classes
eigT = eig(T);

% qit=[
%   p(class 1 in generation t);
%   p(class 2 in generation t);
%   p(class 3 in generation t);
%   p(class 4 in generation t);
%   p(class 5 in generation t)];
% Initialize probability of class given generation t. For an F_{1} (t=1)
% from the initial mating of homozygous parents (ie AABB x aabb), all
% individuals in the F_{1} are of class 4 (ie AaBb in coupling (AB/ab))
qi1=[0;0;0;1;0];
qi2 = T*qi1; 
qi3 = T*qi2; 
qi4 = T*qi3;
                                                                 
% bclass = [
%   p(class in F1);
%   p(class in F2);
%   p(class in F3);
%   p(class in F4)];
% Set up the frequences directly in F_{1}, F_{2}, F_{3}, and F_{4} for each
% class
b1=[qi1(1,1);qi2(1,1);qi3(1,1);qi4(1,1)];
b2=[qi1(2,1);qi2(2,1);qi3(2,1);qi4(2,1)];
b3=[qi1(3,1);qi2(3,1);qi3(3,1);qi4(3,1)];
b4=[qi1(4,1);qi2(4,1);qi3(4,1);qi4(4,1)];
b5=[qi1(5,1);qi2(5,1);qi3(5,1);qi4(5,1)];
% Set up the 4 linear equations (for each generation t=1,2,3,4)
A=[
eigT(1,1)^1 eigT(2,1)^1 eigT(3,1)^1 eigT(4,1)^1;
eigT(1,1)^2 eigT(2,1)^2 eigT(3,1)^2 eigT(4,1)^2;
eigT(1,1)^3 eigT(2,1)^3 eigT(3,1)^3 eigT(4,1)^3;
eigT(1,1)^4 eigT(2,1)^4 eigT(3,1)^4 eigT(4,1)^4];

% We now have a system of 4 linear equations with 4 unknowns for each class
% A*[coefficients of general solution]=bclass
x1=linsolve(A,b1);
x2=linsolve(A,b2);
x3=linsolve(A,b3);
x4=linsolve(A,b4);
x5=linsolve(A,b5);

q_it=[eigT(1,1)^t eigT(2,1)^t eigT(3,1)^t eigT(4,1)^t];

% pclass is the probability of class i (where i=1,2,3,4,5) given
% heterozygosity maintained h, recombination r, and generation t
p1=q_it*x1;
p2=q_it*x2;
p3=q_it*x3;
p4=q_it*x4;
p5=q_it*x5;

% The assignments below are to obtain the chance of ultimate absorption in
% to the absorbing states (class 1 and class 2 genotypes) from the
% transient states (class 3, class 4, and class 5 genotypes). More
% specifically, since all F_{1} individuals are theoretically of class 4
% transient state, then R(I-Q)^{-1}[1,2] and R(I-Q)^{-1}[2,2] are
% probabilities of ultimate aborption into class 1 and 2, respectively. 
%
% (I-Q)
iminusq = [
(1-h^1), -(8*u^1*r^1*(1-r^1)^1)/d^1, -(8*u^1*r^1*(1-r^1)^1)/d^1;
0, (d^1-2*u^2*(1-r^1)^2)/d^1, -(2*u^2*r^2)/d^1;
0, -(2*u^2*r^2)/d^1, (d^1-2*u^2*(1-r^1)^2)/d^1];
%
% (I-Q)^{-1} 
inviminusq=inv(iminusq);
%
% R 
rm = [
(1-h^1)/2, (2*((1-r^1)^2))/d^1 , (2*r^2)/d^1;
(1-h^1)/2, (2*r^2)/d^1 , (2*((1-r^1)^2))/d^1];
%
% R(I-Q)^{-1} 
rminviminusq=rm*inviminusq;
% This ends the calculation of the chance of ultimate absorption into the
% absorbing states. One can check that for h=1/2, R(I-Q)^{-1}[1,2] and 
% R(I-Q)^{-1}[2,2] are 1/(1+2r) and 2r/(1+2r), respectively, as found by
% Bulmer in his book Mathematical Theory of Quantitative Genetics (1980). 
