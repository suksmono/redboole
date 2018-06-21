%% simulation of Perdomo's paper
% A.Perdomo et al., "Construction of model Hamiltonian ...," 
% PRA, 78, 012320(2008)
% converting k-body to 2-body conversion
%% Conclusion: delta >max(H)
%NB: (1) possibly H should have unique minimum
addpath functions;
%% define variables

syms q1 q2 q3 q4 q5 q6;
syms s1 s2 s3 s4 s5 s6;
%% dummy hamiltonian
H=1+q1-q2+q3+q4-q1*q2*q3+q1*q2*q3*q4;
%% delta should be larger than max(H), in this case >4
d5=5;
d6=6;
H2b= simplify( ...
            subs(H,{q1*q2, q3*q4},{q5, q6}) ...
            + H2body(q1,q2,q5,d5) ...
            + H2body(q3,q4,q6,d6) ...
        ) ;
disp(sprintf('The 2-body interaction H2b:\n %s',char(H2b)));
%% simplified into
%   q1 - q2 + q3 + q4 + 15*q5 + 18*q6 
% + 5*q1*q2 - 10*q1*q5 - 10*q2*q5 
% + 6*q3*q4 - q3*q5 - 12*q3*q6 
% - 12*q4*q6 + q5*q6 + 1

% all possible states and related energy
% [q4 q3 q2 q1]
vbin4=[ ...
        0 0 0 0 0 0;
        0 0 0 0 0 1;
        0 0 0 0 1 0;
        0 0 0 0 1 1;
        0 0 0 1 0 0;
        0 0 0 1 0 1;
        0 0 0 1 1 0;
        0 0 0 1 1 1;
        0 0 1 0 0 0;
        0 0 1 0 0 1;
        0 0 1 0 1 0;
        0 0 1 0 1 1;
        0 0 1 1 0 0;
        0 0 1 1 0 1;
        0 0 1 1 1 0;
        0 0 1 1 1 1;

        0 1 0 0 0 0;
        0 1 0 0 0 1;
        0 1 0 0 1 0;
        0 1 0 0 1 1;
        0 1 0 1 0 0;
        0 1 0 1 0 1;
        0 1 0 1 1 0;
        0 1 0 1 1 1;
        0 1 1 0 0 0;
        0 1 1 0 0 1;
        0 1 1 0 1 0;
        0 1 1 0 1 1;
        0 1 1 1 0 0;
        0 1 1 1 0 1;
        0 1 1 1 1 0;
        0 1 1 1 1 1;
        
        1 0 0 0 0 0;
        1 0 0 0 0 1;
        1 0 0 0 1 0;
        1 0 0 0 1 1;
        1 0 0 1 0 0;
        1 0 0 1 0 1;
        1 0 0 1 1 0;
        1 0 0 1 1 1;
        1 0 1 0 0 0;
        1 0 1 0 0 1;
        1 0 1 0 1 0;
        1 0 1 0 1 1;
        1 0 1 1 0 0;
        1 0 1 1 0 1;
        1 0 1 1 1 0;
        1 0 1 1 1 1;

        1 1 0 0 0 0;
        1 1 0 0 0 1;
        1 1 0 0 1 0;
        1 1 0 0 1 1;
        1 1 0 1 0 0;
        1 1 0 1 0 1;
        1 1 0 1 1 0;
        1 1 0 1 1 1;
        1 1 1 0 0 0;
        1 1 1 0 0 1;
        1 1 1 0 1 0;
        1 1 1 0 1 1;
        1 1 1 1 0 0;
        1 1 1 1 0 1;
        1 1 1 1 1 0;
        1 1 1 1 1 1;
];
 % LIST all of energy
 for m=1:2^6
     v1=vbin4(m,:);
     q1=v1(6); q2=v1(5); q3=v1(4); q4=v1(3); q5=v1(2); q6=v1(1);
     E_4bdy = 1+double(q1)-double(q2 )+double(q3) ...
             + double(q4)-double(and(and(q1,q2),q3)) ...
             + double(and(and(and(q3,q4),q2),q1));
    %% 2-body form:
    %   q1 - q2 + q3 + q4 + 15*q5 + 18*q6 
    % + 5*q1*q2 - 10*q1*q5 - 10*q2*q5 
    % + 6*q3*q4 - q3*q5 - 12*q3*q6 
    % - 12*q4*q6 + q5*q6 + 1
     E_2bdy= double(q1)- double(q2)+ double(q3)+ double(q4)+ 15*double(q5) + 18*double(q6) ...
            + 5*double(and(q1,q2)) -10*double(and(q1,q5))-10*double(and(q2,q5))...
            + 6* double(and(q3,q4))- double(and(q3,q5)) -12*double(and(q3,q6)) ...
            -12*double(and(q4,q6))+ double(and(q5,q6)) ...
            + 1;
     disp(sprintf('%2i>> v1= %i%i %i%i%i%i, E4b=%g | E2b=%g',m, v1,E_4bdy, E_2bdy));
 end
 %% RESULTS
%  1>> v1= 00 0000, E4b=1 | E2b=1
%  2>> v1= 00 0001, E4b=2 | E2b=2
%%
%  3>> v1= 00 0010, E4b=0 | E2b=0
%%
%  4>> v1= 00 0011, E4b=1 | E2b=6
%  5>> v1= 00 0100, E4b=2 | E2b=2
%  6>> v1= 00 0101, E4b=3 | E2b=3
%  7>> v1= 00 0110, E4b=1 | E2b=1
%  8>> v1= 00 0111, E4b=1 | E2b=7
%  9>> v1= 00 1000, E4b=2 | E2b=2
% 10>> v1= 00 1001, E4b=3 | E2b=3
% 11>> v1= 00 1010, E4b=1 | E2b=1
% 12>> v1= 00 1011, E4b=2 | E2b=7
% 13>> v1= 00 1100, E4b=3 | E2b=9
% 14>> v1= 00 1101, E4b=4 | E2b=10
% 15>> v1= 00 1110, E4b=2 | E2b=8
% 16>> v1= 00 1111, E4b=3 | E2b=14
% 17>> v1= 01 0000, E4b=1 | E2b=16
% 18>> v1= 01 0001, E4b=2 | E2b=7
% 19>> v1= 01 0010, E4b=0 | E2b=5
% 20>> v1= 01 0011, E4b=1 | E2b=1
% 21>> v1= 01 0100, E4b=2 | E2b=16
% 22>> v1= 01 0101, E4b=3 | E2b=7
% 23>> v1= 01 0110, E4b=1 | E2b=5
% 24>> v1= 01 0111, E4b=1 | E2b=1
% 25>> v1= 01 1000, E4b=2 | E2b=17
% 26>> v1= 01 1001, E4b=3 | E2b=8
% 27>> v1= 01 1010, E4b=1 | E2b=6
% 28>> v1= 01 1011, E4b=2 | E2b=2
% 29>> v1= 01 1100, E4b=3 | E2b=23
% 30>> v1= 01 1101, E4b=4 | E2b=14
% 31>> v1= 01 1110, E4b=2 | E2b=12
% 32>> v1= 01 1111, E4b=3 | E2b=8
% 33>> v1= 10 0000, E4b=1 | E2b=19
% 34>> v1= 10 0001, E4b=2 | E2b=20
% 35>> v1= 10 0010, E4b=0 | E2b=18
% 36>> v1= 10 0011, E4b=1 | E2b=24
% 37>> v1= 10 0100, E4b=2 | E2b=8
% 38>> v1= 10 0101, E4b=3 | E2b=9
% 39>> v1= 10 0110, E4b=1 | E2b=7
% 40>> v1= 10 0111, E4b=1 | E2b=13
% 41>> v1= 10 1000, E4b=2 | E2b=8
% 42>> v1= 10 1001, E4b=3 | E2b=9
% 43>> v1= 10 1010, E4b=1 | E2b=7
% 44>> v1= 10 1011, E4b=2 | E2b=13
% 45>> v1= 10 1100, E4b=3 | E2b=3
% 46>> v1= 10 1101, E4b=4 | E2b=4
% 47>> v1= 10 1110, E4b=2 | E2b=2
% 48>> v1= 10 1111, E4b=3 | E2b=8
% 49>> v1= 11 0000, E4b=1 | E2b=35
% 50>> v1= 11 0001, E4b=2 | E2b=26
% 51>> v1= 11 0010, E4b=0 | E2b=24
% 52>> v1= 11 0011, E4b=1 | E2b=20
% 53>> v1= 11 0100, E4b=2 | E2b=23
% 54>> v1= 11 0101, E4b=3 | E2b=14
% 55>> v1= 11 0110, E4b=1 | E2b=12
% 56>> v1= 11 0111, E4b=1 | E2b=8
% 57>> v1= 11 1000, E4b=2 | E2b=24
% 58>> v1= 11 1001, E4b=3 | E2b=15
% 59>> v1= 11 1010, E4b=1 | E2b=13
% 60>> v1= 11 1011, E4b=2 | E2b=9
% 61>> v1= 11 1100, E4b=3 | E2b=18
% 62>> v1= 11 1101, E4b=4 | E2b=9
% 63>> v1= 11 1110, E4b=2 | E2b=7
% 64>> v1= 11 1111, E4b=3 | E2b=3
 %% Check consistencies of q<->s transform in

