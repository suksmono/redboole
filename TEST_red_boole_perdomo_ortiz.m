%------------------------------------------------------------------------
% redBoole: Boolean reduction of k-body spin interaction to 2-body
%  for implementation of quantum computing on D-Wave
%  Created by: suksmono@{STEI-ITB,MDR Inc.}
%  Example: Case taken from Protein Foldings's of Perdomo-Ortiz paper, 
%   "Finding low-energy ...", Sci.Rep.
%------------------------------------------------------------------------
clear all; clc;
addpath functions; % s2q, q2s, g2bdy ... are located here
%% define symbolic boolean variables
syms q1 q2 q3 q4; 
%% The problem: Eq.(S12) in q-domain  
S12= 4 -3*q1 + 4*q2 -4*q1*q2 -q3 ...
     + q1*q3 -2*q2*q3 +4*q4 -2*q1*q4 ...
     - 8*q2*q4 + 5*q1*q2*q4 -2*q3*q4 ... 
     + 5*q2*q3*q4 - q1*q2*q3*q4;

%% reduce 3- and 4-body to 2-body, do substition,compensate with H2body 
% q1*q2 ->q5; q3*q4 ->q6 ; H2body(qi,qj,qk;dij)
% from the paper: d12=6, and d34=4
syms q5 q6 d12 d34;
d12=6; d34=4;
%
S12_2bdy= simplify( ...
            subs(S12,{q1*q2, q3*q4},{q5,q6}) ...
            + H2body(q1,q2,q5,d12) ...
            + H2body(q3,q4,q6,d34) ...
        ) ;
%% intermediate result
% S12_2bdy = 14*q5 + 10*q6 - 3*q1 + 4*q2 - q3 + 4*q4 - q5*q6 - 12*q5*q1 
%           - 12*q5*q2 + 5*q6*q2 + 5*Q5*q4 - 8*q6*q3 - 8*q6*q4 + 6*q1*q2 
%           + q1*q3 - 2*q1*q4 - 2*q2*q3 - 8*q2*q4 + 4*q3*q4 + 4

%% THE FINAL RESULT
syms s1 s2 s3 s4 s5 s6;
S12_s = simplify( ... 
                 subs( S12_2bdy,...
                       {q1,q2,q3,q4,q5,q6},...
                       {s2q(s1),s2q(s2),s2q(s3),s2q(s4),s2q(s5),s2q(s6)}...
                 )...
                );
disp(sprintf('The 2-body interaction is:\n %s',char(S12_s)));
% S12_s = 10+ (13*s1)/4 + (3*s2)/4 + (7*s3)/4 + s4/4 - 2*s5 - 2*s6 
%         + (3*s1*s2)/2 + (s1*s3)/4 - (s1*s4)/2 - (s2*s3)/2 - 3*s1*s5 - 2*s2*s4 
%         - 3*s2*s5 + s3*s4 + (5*s2*s6)/4 - 2*s3*s6 + (5*s4*s5)/4 - 2*s4*s6 - (s5*s6)/4
%% which is consistent with the paper !!!
 