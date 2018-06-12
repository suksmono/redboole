%% use symbolic computation for Boolean reduction
% redBoole: Boolean reduction of k-body spin interaction to 2-body
% Example: Case taken from Perdomo-Ortiz paper, 
%   "Finding low-energy ...", Sci.Rep.
% Eq.(S12)= 4 -3*q1 + 4*q2 -4*q1*q2 -q3
%           + q1*q3 -2*q2*q3 +4*q4 -2*q1*q4
%           -8*q2*q4 + 5*q1*q2*q4 -2*q3*q4 
%           + 5*q2*q3*q4 - q1*q2*q3*q4
syms q1 q2 q3 q4 s(x) q(x); % define boolean variables and functions
s(x)=1-2*x;
%% problem in q-domeain  
S12=  4 -3*q1 + 4*q2 -4*q1*q2 -q3 ...
     + q1*q3 -2*q2*q3 +4*q4 -2*q1*q4 ...
     - 8*q2*q4 + 5*q1*q2*q4 -2*q3*q4 ... 
     + 5*q2*q3*q4 - q1*q2*q3*q4;

%% reduce 3- and 4-body to 2-body, do substition  
% q1*q2 ->Q5; q3*q4 ->Q6 ;% use capital Q instead of q_tilde
% H_5(q1,q2,Q5;d5)=d5*(3*Q5+q1*q2-2*q1*Q5-2*q2*q5)
% H_6(q3,q4,Q5;d6)=d6*(3*Q6+q3*q4-2*q3*Q6-2*q4*Q6)
syms H_5 H_6 Q5 Q6 d12 d34
%% tr d12=6, and d34=4, from the paper
d12=6;
d34=4;
H_12=d12*(3*Q5+q1*q2-2*q1*Q5-2*q2*Q5);
H_34=d34*(3*Q6+q3*q4-2*q3*Q6-2*q4*Q6);
%
S12_2bdy = simplify( 4 -3*q1 + 4*q2 -4*Q5 -q3 ...
      + q1*q3 -2*q2*q3 +4*q4 -2*q1*q4 ...
      - 8*q2*q4 + 5*Q5*q4 -2*Q6 ... 
      + 5*q2*Q6 - Q5*Q6 + ...
      + H_12 + H_34);
%% Result
% S12_2bdy = 14*Q5 + 10*Q6 - 3*q1 + 4*q2 - q3 + 4*q4 - Q5*Q6 - 12*Q5*q1 
%           - 12*Q5*q2 + 5*Q6*q2 + 5*Q5*q4 - 8*Q6*q3 - 8*Q6*q4 + 6*q1*q2 
%           + q1*q3 - 2*q1*q4 - 2*q2*q3 - 8*q2*q4 + 4*q3*q4 + 4
  
%% transform to s
syms s1 s2 s3 s4 s5 s6;
q(x)=(1/2)*(1-x); % ransform formula
q1=q(s1);
q2=q(s2);
q3=q(s3);
q4=q(s4);
Q5=q(s5);
Q6=q(s6);
%%
S12_s = simplify(14*Q5 + 10*Q6 - 3*q1 + 4*q2 - q3 + 4*q4 - Q5*Q6 - 12*Q5*q1 ...
          - 12*Q5*q2 + 5*Q6*q2 + 5*Q5*q4 - 8*Q6*q3 - 8*Q6*q4 + 6*q1*q2 ...
          + q1*q3 - 2*q1*q4 - 2*q2*q3 - 8*q2*q4 + 4*q3*q4 + 4);
%% show result
S12_s,
% THE RESULT
% S12_s = 10+ (13*s1)/4 + (3*s2)/4 + (7*s3)/4 + s4/4 - 2*s5 - 2*s6 
%         + (3*s1*s2)/2 + (s1*s3)/4 - (s1*s4)/2 - (s2*s3)/2 - 3*s1*s5 - 2*s2*s4 
%         - 3*s2*s5 + s3*s4 + (5*s2*s6)/4 - 2*s3*s6 + (5*s4*s5)/4 - 2*s4*s6 - (s5*s6)/4
%% CORRECT