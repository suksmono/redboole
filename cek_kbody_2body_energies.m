%% Check consistencies of q<->s transform in
% converting k-body to 2-body conversion
addpath functions;
%% define variables

syms q1 q2 q3 q4;
syms s1 s2 s3 s4;
%% dummy hamiltonian
H=q1*q2*q3;
%% ------------------------------------------------------
% should be simplified into H2bdy=q4*q3+H(q1,q2;q4,d4);
% for d4=8 => H =d4(3*q4+q1*q2-2*q1*q4-2*q2*q4)
% H2bdy=q3*q4 + 24*q4 + 8q1*q2 - 16*q1*q4 - 16*q2*q4
% ------------------------------------------------------
d4=8;
H2b= simplify( ...
            subs(H,{q1*q2},{q4}) ...
            + H2body(q1,q2,q4,d4) ...
        ) ;
disp(sprintf('The 2-body interaction H2b:\n %s',char(H2b)));
% all possible states and related energy
% [q4 q3 q2 q1]
vbin4=[ 0 0 0 0;
        0 0 0 1;
        0 0 1 0;
        0 0 1 1;
        0 1 0 0;
        0 1 0 1;
        0 1 1 0;
        0 1 1 1;
        1 0 0 0;
        1 0 0 1;
        1 0 1 0;
        1 0 1 1;
        1 1 0 0;
        1 1 0 1;
        1 1 1 0;
        1 1 1 1 ...
                ];
 % LIST all of energy
 for m=1:2^4
     v1=vbin4(m,:);
     q1=v1(4); q2=v1(3); q3=v1(2); q4=v1(1);
     E3= double(and(and(q1,q2),q3));
     E_ppr = 1+double(q1)-double(q2)+double(q3) ...
             + double(q4)-double(and(and(q1,q2),q3)) ...
             + double(and(and(and(q3,q4),q2),q1));
     %% 24*q4 + 8*q1*q2 - 16*q1*q4 - 16*q2*q4 + q3*q4
     E4 = 24*double(q4)+ 8*double(and(q1,q2))-...
          16*double(and(q1,q4))- 16*double(and(q2,q4))+ double(and(q3,q4));
     disp(sprintf('%2i>> v1= %i %i %i %i, E3=%g | E4=%g | Eppr=%g',m, v1,E3, E4, E_ppr));
 end
 
