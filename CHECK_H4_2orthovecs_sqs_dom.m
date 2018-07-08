% CHECK_H4_sqs_dom.m
%% use symbolic computation for Boolean reduction
% s={-1,1}; q={0,1}
% keep in mind that boolean reduction works on q-domain
%            problem               ancilla
% case: H2=[ q11 q12 q13 q14;       | q15
%            q21 q22 q23 q24;       | q25
%            q31 q32 q33 q34;       | q35
%            q41 q42 q43 q44]       | q45

% off-diagonal term is dij
% ----------------------------------------------------------------------
% 1. Formulate Hamiltonian: Hks
% 2. Transform to q-domain: Hks-> Hkq
% 3. Boolean reduction: Hkq->H2q
% 4. Transform-back to s-domain: H2q->H2s
% 5. Extract Ising params and apply SA
% ----------------------------------------------------------------------
clear all;
addpath functions;
syms q11 q12 q13 q14; % define boolean variables 
syms q21 q22 q23 q24; % define boolean variables 
syms q31 q32 q33 q34; % define boolean variables 
syms q41 q42 q43 q44; % define boolean variables

% syms q17 q18 q19 q20;
%% 1. Formulate Hamiltonian Hks
% ==> d12=s1*s3 + s2*s4 ;
syms s11 s12 s13 s14;
syms s21 s22 s23 s24;
syms s31 s32 s33 s34;
syms s41 s42 s43 s44;
%% store as matrix of symbols
q=[ q11 q12 q13 q14; ...
    q21 q22 q23 q24; ...
    q31 q32 q33 q34; ...
    q41 q42 q43 q44 ...
    ];

s=[ s11 s12 s13 s14; ...
    s21 s22 s23 s24; ...
    s31 s32 s33 s34; ...
    s41 s42 s43 s44 ...
    ];

    
%% define substitute of qi*qj-> qk
%all qi*qj, for <v1,v2>: 
% q1*q5->q_{M*M+1}, q2*q6->q_{M*M+2}, q3*q7->_{M*M+3}, and q4*q8->q_{M*M+4}:
% in general: qi*q_{i+M*(c-1)}; where c is rcolumn number, M: H-order
%  to be replaced with q_{M*M+M*(c-1)+r}, r: row number
d12=transpose(s(:,1))*s(:,2); 
disp(sprintf('d12 in s-domain: \n%s',char(d12)));
% > s11*s12 + s21*s22 + s31*s32 + s41*s42
d12_2=expand(d12^2);
%% change all si^2 to 1
Hks=simplify( ...
                 subs(d12_2, ...
                 {s11^2, s21^2, s31^2, s41^2, s12^2, s22^2, s32^2, s42^2}, ...
                 {1, 1, 1, 1, 1, 1, 1, 1} ...
                 ) ...
             );            
disp(sprintf('Hks: \n%s',char(Hks)));
%% Hks =  2*s11*s12*s21*s22 + 2*s11*s12*s31*s32 + 2*s11*s12*s41*s42 
%       + 2*s21*s22*s31*s32 + 2*s21*s22*s41*s42 + 2*s31*s32*s41*s42 + 4 
%  max(Hks)=2*6+2=16

 
%% 2. Transform to q-domain: Hks-> Hkq
Hkq=expand(simplify( ...
                subs(Hks, ...
                {s11, s21, s31, s41, s12, s22, s32, s42}, ...
                {s2q(q11), s2q(q21), s2q(q31),s2q(q41), ...
                 s2q(q12), s2q(q22) ,s2q(q32),s2q(q42)}...
                )...
            ));
disp(sprintf('Hkq: \n%s',char(Hkq)));
%%
% ------------------------------------------------------------------------
% Hkq = 
% - 12*q12 - 12*q21 - 12*q22 - 12*q31 - 12*q32 - 12*q41 - 12*q42 - 12*q11 
% +24*q11*q12+ 8*q11*q21 + 8*q11*q22 + 8*q12*q21 + 8*q12*q22 + 8*q11*q31 
% + 8*q11*q32 + 8*q12*q31 + 24*q21*q22 + 8*q12*q32 + 8*q11*q41 + 8*q21*q31 
% + 8*q11*q42 + 8*q12*q41 + 8*q21*q32 + 8*q22*q31 + 8*q12*q42 + 8*q22*q32 
% + 8*q21*q41 + 8*q21*q42 + 8*q22*q41 + 24*q31*q32 + 8*q22*q42 + 8*q31*q41 
% + 8*q31*q42 + 8*q32*q41 + 8*q32*q42 + 24*q41*q42
%%
% - 16*q11*q12*q21 - 16*q11*q12*q22 - 16*q11*q12*q31 - 16*q11*q21*q22 
% - 16*q11*q12*q32 - 16*q12*q21*q22 - 16*q11*q12*q41 - 16*q11*q12*q42 
% - 16*q11*q31*q32 - 16*q21*q22*q31 - 16*q12*q31*q32 - 16*q21*q22*q32 
% - 16*q21*q22*q41 - 16*q21*q31*q32 - 16*q21*q22*q42 - 16*q22*q31*q32 
% - 16*q11*q41*q42 - 16*q12*q41*q42 - 16*q21*q41*q42 - 16*q31*q32*q41 
% - 16*q22*q41*q42 - 16*q31*q32*q42 - 16*q31*q41*q42 - 16*q32*q41*q42 
%%
% + 32*q11*q12*q21*q22 + 32*q11*q12*q31*q32 + 32*q11*q12*q41*q42 
% + 32*q21*q22*q31*q32 + 32*q21*q22*q41*q42 + 32*q31*q32*q41*q42 + 16

%  inter-vector products 
% ------------------------------------------------------------------------

%%  3. Boolean reduction: Hkq->H2q
% choose delta=2*Hmax = 2*16=32
% q11*q12->q17, q21*q22->q18 , q31*q32->q19, q41*q42->q20; 
% [ 11 .... 14| 15]
% [ 21 .... 24| 25]  
% substitute of col1*col2
syms q15 q25 q35 q45;
qs12=[q15 q25 q35 q45]; %
syms s15 s25 s35 s45;
ss12=[s15 s25 s35 s45];

%
d=32;
H2q = simplify( ...
      subs( Hkq, ...
         {q(1,1)*q(1,2), q(2,1)*q(2,2), q(3,1)*q(3,2), q(4,1)*q(4,2)}, ...
         {qs12(1), qs12(2), qs12(3), qs12(4)}...
                )...
            + H2body(q(1,1),q(1,2),qs12(1),d) ...
            + H2body(q(2,1),q(2,2),qs12(2),d) ...
            + H2body(q(3,1),q(3,2),qs12(3),d) ...
            + H2body(q(4,1),q(4,2),qs12(4),d) ...
            );
disp(sprintf('H2q: \n%s',char(expand(H2q))));

%% RESULTS
% H2q: 
%  120*q15 - 12*q12 - 12*q11 - 12*q21 - 12*q22 + 120*q25 - 12*q31 
% - 12*q32 + 120*q35 - 12*q41 - 12*q42 + 120*q45 
% + 32*q11*q12 - 64*q11*q15 - 64*q12*q15 + 8*q11*q21 + 8*q11*q22 
% + 8*q12*q21 + 8*q12*q22 - 16*q11*q25 - 16*q15*q21 - 16*q12*q25 
% - 16*q15*q22 + 32*q15*q25 + 8*q11*q31 + 8*q11*q32 + 8*q12*q31 
% + 32*q21*q22 + 8*q12*q32 - 16*q11*q35 - 16*q15*q31 - 64*q21*q25
% - 16*q12*q35 - 16*q15*q32 - 64*q22*q25 + 32*q15*q35 + 8*q11*q41 
% + 8*q21*q31 + 8*q11*q42 + 8*q12*q41 + 8*q21*q32 + 8*q22*q31 
% + 8*q12*q42 + 8*q22*q32 - 16*q11*q45 - 16*q15*q41 - 16*q21*q35 
% - 16*q25*q31 - 16*q12*q45 - 16*q15*q42 - 16*q22*q35 - 16*q25*q32 
% + 32*q15*q45 + 32*q25*q35 + 8*q21*q41 + 8*q21*q42 + 8*q22*q41 
% + 32*q31*q32 + 8*q22*q42 - 16*q21*q45 - 16*q25*q41 - 64*q31*q35 
% - 16*q22*q45 - 16*q25*q42 - 64*q32*q35 + 32*q25*q45 + 8*q31*q41 
% + 8*q31*q42 + 8*q32*q41 + 8*q32*q42 - 16*q31*q45 - 16*q35*q41 
% - 16*q32*q45 - 16*q35*q42 + 32*q35*q45 + 32*q41*q42 - 64*q41*q45 - 64*q42*q45 + 16


%% 4. Transform-back to s-domain: H2q->H2s
H2s = simplify( ...
            subs(H2q, ...
                { q11,q21,q31,q41, ... 
                  q12,q22,q32,q42,...
                  q15,q25,q35,q45}, ...
                { q2s(s11), q2s(s21),q2s(s31),q2s(s41), ...
                  q2s(s12), q2s(s22),q2s(s32),q2s(s42), ...
                  q2s(s15), q2s(s25),q2s(s35),q2s(s45) ...
                  }...
                )...
            );
disp(sprintf('H2s 2-body in s-domain: \n%s',...
    char(expand(H2s))));
% H2s 2-body in s-domain: 
%   14*s11 + 14*s12 - 28*s15 + 14*s21 + 14*s22 - 28*s25 + 14*s31 
% + 14*s32 - 28*s35 + 14*s41 + 14*s42 - 28*s45 
% + 8*s11*s12 - 16*s11*s15 - 16*s12*s15 + 2*s11*s21 + 2*s11*s22 
% + 2*s12*s21 + 2*s12*s22 - 4*s11*s25 - 4*s15*s21 - 4*s12*s25
% - 4*s15*s22 + 8*s15*s25 + 2*s11*s31 + 2*s11*s32 + 2*s12*s31 
% + 8*s21*s22 + 2*s12*s32 - 4*s11*s35 - 4*s15*s31 - 16*s21*s25 
% - 4*s12*s35 - 4*s15*s32 - 16*s22*s25 + 8*s15*s35 + 2*s11*s41 
% + 2*s21*s31 + 2*s11*s42 + 2*s12*s41 + 2*s21*s32 + 2*s22*s31 
% + 2*s12*s42 + 2*s22*s32 - 4*s11*s45 - 4*s15*s41 - 4*s21*s35 
% - 4*s25*s31 - 4*s12*s45 - 4*s15*s42 - 4*s22*s35 - 4*s25*s32 
% + 8*s15*s45 + 8*s25*s35 + 2*s21*s41 + 2*s21*s42 + 2*s22*s41 
% + 8*s31*s32 + 2*s22*s42 - 4*s21*s45 - 4*s25*s41 - 16*s31*s35 
% - 4*s22*s45 - 4*s25*s42 - 16*s32*s35 + 8*s25*s45 + 2*s31*s41 
% + 2*s31*s42 + 2*s32*s41 + 2*s32*s42 - 4*s31*s45 - 4*s35*s41 
% - 4*s32*s45 - 4*s35*s42 + 8*s35*s45 + 8*s41*s42 - 16*s41*s45 
% - 16*s42*s45 + 112


% ------------------------------------------------------------
%  EXTRACT THE RESULTS FROM SYMBOLIC SOLUTION
% ------------------------------------------------------------
% now, we are working on orthogonalization of 2-binary vectors
% q11 ... q41, q12 .... q42, q15, ..., q45 =? 4*3=12 qubits 
%     [cx,tx]=coeffs(H2s,[s1 s2 s3 s4 s5 s6]);
[cx,tx]=coeffs(H2s,[ s11 s21 s31 s41 ...
                     s12 s22 s32 s42 ...
                     s15 s25 s35 s45 ...
                     ]);
% number of  qubits: q11, ..., q45 => 3*4=12
% [ s11*s21, s11*s31, s11*s41, s11*s12, s11*s22, s11*s32, s11*s42, s11*s15,
%   s11*s25, s11*s35, s11*s45, s11, s21*s31, s21*s41, s12*s21, s21*s22,
%   s21*s32, s21*s42, s15*s21, s21*s25, s21*s35, s21*s45, s21, s31*s41, 
%   s12*s31, s22*s31, s31*s32, s31*s42, s15*s31, s25*s31, s31*s35, s31*s45, 
%   s31, s12*s41, s22*s41, s32*s41, s41*s42, s15*s41, s25*s41, s35*s41, 
%   s41*s45, s41, s12*s22, s12*s32, s12*s42, s12*s15, s12*s25, s12*s35, 
%   s12*s45, s12, s22*s32, s22*s42, s15*s22, s22*s25, s22*s35, s22*s45, 
%   s22, s32*s42, s15*s32, s25*s32, s32*s35, s32*s45, s32, s15*s42, s25*s42, 
%   s35*s42, s42*s45, s42, s15*s25, s15*s35, s15*s45, s15, s25*s35, s25*s45, 
%   s25, s35*s45, s35, s45, 1]

%% create lookup tbl for product of s_ij*s_mn
sx=[transpose(s(:,1)) transpose(s(:,2)) ss12];
tblook=transpose(sx)*sx;
NQ=12;
Jij=zeros(NQ,NQ); % init Jij: coupling coefficients
hi=zeros(NQ,1);
b=0; % bias
%% fill in b
NTERM=length(cx);
for m=1:NTERM
    if (tx(m)==1) 
        b=cx(m);
    end
end;
%% fill in hi
for m=1:NQ
    tsym=sx(m);
    for nn=1:NTERM
        if(tx(nn)==tsym)
            hi(m)=cx(nn);
        end;
    end;
end;
%% fill in Jij, we only consider upper diagonal of Jij
for m=1:NQ
    for nn=m+1:NQ
        % get quadratic terms s_ij*s_mn
        tquad = tblook(m,nn);
        % now get the cofficient \
            for kkk =1:NTERM
                if(tx(kkk)==tquad)
                    Jij(m,nn)=cx(kkk);
                end
            end;
    end;
end;


% if 0==1
% ------------------------------------------------------------
% SEND THE RESULTS (Jij, hi, b) TO SIMULATOR (S.A.)
% ------------------------------------------------------------
ITR1=100;%2*10*1000; %set 1 for random init
vStart=0;
%% create a vector of [(N-1)x(N-1),1] spin with random {-,+] values
vSpin0=vRandSpin(NQ);

% Etresh= 0;%-2*(4*K-1)*(2*K-1);%N/4; % target energy
% annealing schedule
maxIter=1*10*100;
Ptresh=pMCSchedule(maxIter,0.5);
%-- Ptresh is the schedule
k=0;
vSpin=vSpin0;
% % init energy setting
Ecurr = vEnergyBx(b, hi, Jij,vSpin);
% NFLIP=N-NO;
Etresh=0;
while k<maxIter && Ecurr>Etresh
    k=k+1;
    Eprev=Ecurr;
    % try Monte Carlo
    % --flip a row randomly
    RD=randperm(NQ);
    % copy current spin to a template
    tvSpin=vSpin;
    tvSpin(RD(1)) = -1*tvSpin(RD(1));
    Ecurr=vEnergyBx(b, hi, Jij,tvSpin);
    if (Ecurr<Eprev)||rand>Ptresh(k) %accept if lower
        vSpin=tvSpin;
    end;
    Esys=vEnergyBx(b, hi, Jij,vSpin);
    Ecurr=Esys;
    vE(k)=Esys;
    disp(sprintf('LQ%d: Iter->%d, Tresh=%1.4f, E=%4.2f',...
        NQ,k,Ptresh(k),double(Esys)));
%     disp(sprintf('H%d: Iter->%d',N,k));
end;
% ------------------------------------------------------------
% DISPLAY THE RESULTS
% ------------------------------------------------------------
v1=vSpin(1:4); v2=vSpin(5:8);
disp(sprintf('s11= %2g, s21= %2g, s31= %2g, s41= %2g',v1));
disp(sprintf('s12= %2g, s22= %2g, s32= %2g, s42= %2g',v2));
disp(sprintf('<v1,v2> = %2g', v1'*v2));
% H2=reshape(vSpin(1:4),2,2);
% H2*H2',

if length(vE)>1
    imagesc(1);
    plot(vE);
end;
% end
