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
%% construct s-domain 2-body hamiltonian
H2b_sdom = simplify( ... 
                 subs( H2b,...
                       {q1,q2,q3,q4,q5,q6},...
                       {q2s(s1),q2s(s2),q2s(s3),q2s(s4),q2s(s5),q2s(s6)}...
                 )...
                );
%% Normalize
H2b_sdom=H2b_sdom/max(abs(coeffs(H2b_sdom)));
disp(sprintf('The 2-body interaction in s-domain:\n %s',char(H2b_sdom)));
% result
% H2b_sdom =  (3*s1)/41 + (7*s2)/41 + (5*s3)/41 + (4*s4)/41 - (10*s5)/41 - (13*s6)/41 
%           + (5*s1*s2)/41 - (10*s1*s5)/41 - (10*s2*s5)/41 + (6*s3*s4)/41 - (s3*s5)/41 
%           - (12*s3*s6)/41 - (12*s4*s6)/41 + (s5*s6)/41 + 1

% all possible states and related energy
% [q6 q5 q4 q3 q2 q1]
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
     vs1=vq2s(q1); vs2=vq2s(q2); vs3=vq2s(q3); 
     vs4=vq2s(q4); vs5=vq2s(q5); vs6=vq2s(q6);
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
    % H2b_sdom =  (3*s1)/41 + (7*s2)/41 + (5*s3)/41 + (4*s4)/41 - (10*s5)/41 - (13*s6)/41 
    %           + (5*s1*s2)/41 - (10*s1*s5)/41 - (10*s2*s5)/41 + (6*s3*s4)/41 - (s3*s5)/41 
    %           - (12*s3*s6)/41 - (12*s4*s6)/41 + (s5*s6)/41 + 1
    %     sCf=coeffs(H2b_sdom);
    %% indices are not correct
    %     E_2bdy_sdom = sCf(1)*vs1 + sCf(2)*vs2 + sCf(3)*vs3 + sCf(4)*vs4 + sCf(5)*vs5 + sCf(6)*vs6 ...
    %         + sCf(7)*vs1*vs2 + sCf(8)*vs1*vs5 + sCf(9)*vs2*vs5 + sCf(10)*vs3*vs4 + sCf(11)*vs3*vs5 ...
    %         + sCf(12)*vs3*vs6 + sCf(13)*vs4*vs6 + sCf(14)*vs5*vs6 + 1;

    E_2bdy_sdom = (3/41)*vs1 + (7/41)*vs2 + (5/41)*vs3 + (4/41)*vs4 + -(10/41)*vs5 -(13/41)*vs6 ...
        + (5/41)*vs1*vs2 + (-10/41)*vs1*vs5 + (-10/41)*vs2*vs5 + (6/41)*vs3*vs4 + (-1/41)*vs3*vs5 ...
        + (-12/41)*vs3*vs6 + (-12/41)*vs4*vs6 + (1/41)*vs5*vs6 + 1;

    disp(sprintf('%2i>> v1= %i%i %i%i%i%i, E4b=%2.3g| E2b_q=%2.3g| E2b_s= %2.3g',m, v1,E_4bdy, E_2bdy, double(E_2bdy_sdom)) );
 end
%% RESULTS
%  1>> v1= 00 0000, E4b= 1| E2b_q= 1| E2b_s= 0.0976
%  2>> v1= 00 0001, E4b= 2| E2b_q= 2| E2b_s= 0.195
%  3>> v1= 00 0010, E4b= 0| E2b_q= 0| E2b_s=  0
%  4>> v1= 00 0011, E4b= 1| E2b_q= 6| E2b_s= 0.585
%  5>> v1= 00 0100, E4b= 2| E2b_q= 2| E2b_s= 0.195
%  6>> v1= 00 0101, E4b= 3| E2b_q= 3| E2b_s= 0.293
%  7>> v1= 00 0110, E4b= 1| E2b_q= 1| E2b_s= 0.0976
%  8>> v1= 00 0111, E4b= 1| E2b_q= 7| E2b_s= 0.683
%  9>> v1= 00 1000, E4b= 2| E2b_q= 2| E2b_s= 0.195
% 10>> v1= 00 1001, E4b= 3| E2b_q= 3| E2b_s= 0.293
% 11>> v1= 00 1010, E4b= 1| E2b_q= 1| E2b_s= 0.0976
% 12>> v1= 00 1011, E4b= 2| E2b_q= 7| E2b_s= 0.683
% 13>> v1= 00 1100, E4b= 3| E2b_q= 9| E2b_s= 0.878
% 14>> v1= 00 1101, E4b= 4| E2b_q=10| E2b_s= 0.976
% 15>> v1= 00 1110, E4b= 2| E2b_q= 8| E2b_s= 0.78
% 16>> v1= 00 1111, E4b= 3| E2b_q=14| E2b_s= 1.37
% 17>> v1= 01 0000, E4b= 1| E2b_q=16| E2b_s= 1.56
% 18>> v1= 01 0001, E4b= 2| E2b_q= 7| E2b_s= 0.683
% 19>> v1= 01 0010, E4b= 0| E2b_q= 5| E2b_s= 0.488
% 20>> v1= 01 0011, E4b= 1| E2b_q= 1| E2b_s= 0.0976
% 21>> v1= 01 0100, E4b= 2| E2b_q=16| E2b_s= 1.56
% 22>> v1= 01 0101, E4b= 3| E2b_q= 7| E2b_s= 0.683
% 23>> v1= 01 0110, E4b= 1| E2b_q= 5| E2b_s= 0.488
% 24>> v1= 01 0111, E4b= 1| E2b_q= 1| E2b_s= 0.0976
% 25>> v1= 01 1000, E4b= 2| E2b_q=17| E2b_s= 1.66
% 26>> v1= 01 1001, E4b= 3| E2b_q= 8| E2b_s= 0.78
% 27>> v1= 01 1010, E4b= 1| E2b_q= 6| E2b_s= 0.585
% 28>> v1= 01 1011, E4b= 2| E2b_q= 2| E2b_s= 0.195
% 29>> v1= 01 1100, E4b= 3| E2b_q=23| E2b_s= 2.24
% 30>> v1= 01 1101, E4b= 4| E2b_q=14| E2b_s= 1.37
% 31>> v1= 01 1110, E4b= 2| E2b_q=12| E2b_s= 1.17
% 32>> v1= 01 1111, E4b= 3| E2b_q= 8| E2b_s= 0.78
% 33>> v1= 10 0000, E4b= 1| E2b_q=19| E2b_s= 1.85
% 34>> v1= 10 0001, E4b= 2| E2b_q=20| E2b_s= 1.95
% 35>> v1= 10 0010, E4b= 0| E2b_q=18| E2b_s= 1.76
% 36>> v1= 10 0011, E4b= 1| E2b_q=24| E2b_s= 2.34
% 37>> v1= 10 0100, E4b= 2| E2b_q= 8| E2b_s= 0.78
% 38>> v1= 10 0101, E4b= 3| E2b_q= 9| E2b_s= 0.878
% 39>> v1= 10 0110, E4b= 1| E2b_q= 7| E2b_s= 0.683
% 40>> v1= 10 0111, E4b= 1| E2b_q=13| E2b_s= 1.27
% 41>> v1= 10 1000, E4b= 2| E2b_q= 8| E2b_s= 0.78
% 42>> v1= 10 1001, E4b= 3| E2b_q= 9| E2b_s= 0.878
% 43>> v1= 10 1010, E4b= 1| E2b_q= 7| E2b_s= 0.683
% 44>> v1= 10 1011, E4b= 2| E2b_q=13| E2b_s= 1.27
% 45>> v1= 10 1100, E4b= 3| E2b_q= 3| E2b_s= 0.293
% 46>> v1= 10 1101, E4b= 4| E2b_q= 4| E2b_s= 0.39
% 47>> v1= 10 1110, E4b= 2| E2b_q= 2| E2b_s= 0.195
% 48>> v1= 10 1111, E4b= 3| E2b_q= 8| E2b_s= 0.78
% 49>> v1= 11 0000, E4b= 1| E2b_q=35| E2b_s= 3.41
% 50>> v1= 11 0001, E4b= 2| E2b_q=26| E2b_s= 2.54
% 51>> v1= 11 0010, E4b= 0| E2b_q=24| E2b_s= 2.34
% 52>> v1= 11 0011, E4b= 1| E2b_q=20| E2b_s= 1.95
% 53>> v1= 11 0100, E4b= 2| E2b_q=23| E2b_s= 2.24
% 54>> v1= 11 0101, E4b= 3| E2b_q=14| E2b_s= 1.37
% 55>> v1= 11 0110, E4b= 1| E2b_q=12| E2b_s= 1.17
% 56>> v1= 11 0111, E4b= 1| E2b_q= 8| E2b_s= 0.78
% 57>> v1= 11 1000, E4b= 2| E2b_q=24| E2b_s= 2.34
% 58>> v1= 11 1001, E4b= 3| E2b_q=15| E2b_s= 1.46
% 59>> v1= 11 1010, E4b= 1| E2b_q=13| E2b_s= 1.27
% 60>> v1= 11 1011, E4b= 2| E2b_q= 9| E2b_s= 0.878
% 61>> v1= 11 1100, E4b= 3| E2b_q=18| E2b_s= 1.76
% 62>> v1= 11 1101, E4b= 4| E2b_q= 9| E2b_s= 0.878
% 63>> v1= 11 1110, E4b= 2| E2b_q= 7| E2b_s= 0.683
% 64>> v1= 11 1111, E4b= 3| E2b_q= 3| E2b_s= 0.293

