%% use symbolic computation for Boolean reduction
% s={-1,1}; q={0,1}
% keep in mind that boolean reduction works on q-domain
% case: H2=[q1 q2; q3 q4], off-diagonal term is d12
% ----------------------------------------------------------------------
% 1. Formulate Hamiltonian: Hks
% 2. Transform to q-domain: Hks-> Hkq
% 3. Boolean reduction: Hkq->H2q
% 4. Transform-back to s-domain: H2q->H2s
% 5. Extract Ising params and apply SA
% ----------------------------------------------------------------------
clear all;
addpath functions;
syms q1 q2 q3 q4 q5 q6; % define boolean variables and functions
%% 1. Formulate Hamiltonian Hks
% ==> d12=s1*s3 + s2*s4 ;
syms s1 s2 s3 s4 s5 s6;
d12=s1*s3 + s2*s4 ; 
d12_2=expand(d12^2);
%% change all si^2 to 1
Hks=simplify( ...
                 subs(d12_2, ...
                 {s1^2, s2^2, s3^2, s4^2}, ...
                 {1, 1, 1, 1} ...
                 ) ...
             );            
disp(sprintf('d12^2 in s-domain: \n%s',char(Hks)));
%% 2. Transform to q-domain: Hks-> Hkq
Hkq=expand(simplify( ...
                subs(Hks, ...
                {s1, s2, s3, s4}, ...
                {s2q(q1), s2q(q2),s2q(q3),s2q(q4)}...
                )...
            ));
disp(sprintf('Hkq k-body in q-domain: \n%s',...
    char(Hkq)));
%%
% Hkq =  8*q1*q2 - 4*q2 - 4*q3 - 4*q4 - 4*q1 
%   + 8*q1*q3 + 8*q1*q4 + 8*q2*q3 + 8*q2*q4 + 8*q3*q4 
%   - 16*q1*q2*q3 - 16*q1*q2*q4 - 16*q1*q3*q4 - 16*q2*q3*q4 
%   + 32*q1*q2*q3*q4 + 4

%%  3. Boolean reduction: Hkq->H2q
% q1*q2->q5, q3*q4->q6; choose delta=2*Hmax = 8
d5=8; d6=8;
H2q = simplify( ...
            subs(Hkq, ...
                {q1*q2, q3*q4}, ...
                {q5, q6}...
                )...
            + H2body(q1,q2,q5,d5) ...
            + H2body(q3,q4,q6,d6) ...
            );
disp(sprintf('H2q 2-body in q-domain: \n%s',...
    char(expand(H2q))));
%% RESULTS
% 32*q5 - 4*q2 - 4*q3 - 4*q4 - 4*q1 + 32*q6 
% + 8*q1*q2 + 8*q1*q3 + 8*q1*q4 + 8*q2*q3 - 16*q1*q5 + 8*q2*q4 
% - 16*q1*q6 - 16*q2*q5 + 8*q3*q4 - 16*q2*q6 - 16*q3*q5 
% - 16*q3*q6 - 16*q4*q5 - 16*q4*q6 + 32*q5*q6 + 4
%% 4. Transform-back to s-domain: H2q->H2s
H2s = simplify( ...
            subs(H2q, ...
                {q1,q2,q3,q4,q5,q6}, ...
                {q2s(s1), q2s(s2),q2s(s3),q2s(s4),q2s(s5),q2s(s6)}...
                )...
            );
disp(sprintf('H2s 2-body in s-domain: \n%s',...
    char(expand(H2s))));
%H2s 2-body in s-domain: 
%   4*s1 + 4*s2 + 4*s3 + 4*s4 - 8*s5 - 8*s6 
% + 2*s1*s2 + 2*s1*s3 + 2*s1*s4 + 2*s2*s3 - 4*s1*s5 
% + 2*s2*s4 - 4*s1*s6 - 4*s2*s5 + 2*s3*s4 - 4*s2*s6 
% - 4*s3*s5 - 4*s3*s6 - 4*s4*s5 - 4*s4*s6 + 8*s5*s6 + 16

% ------------------------------------------------------------
% EXTRACT THE RESULTS FROM SYMBOLIC SOLUTION
% ------------------------------------------------------------
[cx,tx]=coeffs(H2s,[s1 s2 s3 s4 s5 s6]);
% number of  qubits: q1, q2, q3, q4, q5, q6 => 6
NQ=6;
Jij=zeros(NQ,NQ); % init Jij: coupling coefficients
hi=zeros(NQ,1);
b=0; % bias

NT=length(cx);
for m=1:NT
    tTerm=tx(m);
    cTerm=cx(m);
    %% strange case
    if double(tTerm == s1)==1
        hi(1)=cTerm;
    end;
    if double(tTerm == s2)==1
        hi(2)=cTerm;
    end;
    if double(tTerm == s3)==1
        hi(3)=cTerm;
    end;
    if double(tTerm == s4)==1
        hi(4)=cTerm;
    end;
    if double(tTerm == s5)==1
        hi(5)=cTerm;
    end;
    if double(tTerm == s6)==1
        hi(6)=cTerm;
    end;
    %%
    if double(tTerm == s1*s2)==1
        Jij(1,2)=cTerm;
    end;
    if double(tTerm == s1*s3)==1
        Jij(1,3)=cTerm;
    end;
    if double(tTerm == s1*s4)==1
        Jij(1,4)=cTerm;
    end;
    if double(tTerm == s1*s5)==1
        Jij(1,5)=cTerm;
    end;
    if double(tTerm == s1*s6)==1
        Jij(1,6)=cTerm;
    end;
    %
    if double(tTerm == s2*s3)==1
        Jij(2,3)=cTerm;
    end;
    if double(tTerm == s2*s4)==1
        Jij(2,4)=cTerm;
    end;
    if double(tTerm == s2*s5)==1
        Jij(2,5)=cTerm;
    end;
    if double(tTerm == s2*s6)==1
        Jij(2,6)=cTerm;
    end;
    %
    if double(tTerm == s3*s4)==1
        Jij(3,4)=cTerm;
    end;
    if double(tTerm == s3*s5)==1
        Jij(3,5)=cTerm;
    end;
    if double(tTerm == s3*s6)==1
        Jij(3,6)=cTerm;
    end;
    %
    if double(tTerm == s4*s5)==1
        Jij(4,5)=cTerm;
    end;
    if double(tTerm == s4*s6)==1
        Jij(4,6)=cTerm;
    end;
    %
    if double(tTerm == s5*s6)==1
        Jij(5,6)=cTerm;
    end;
    %
    if double(tTerm == 1)==1
        b=cTerm;
    end;
end;

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
disp(sprintf('\n s1=%2g, s2=%2g, s3=%2g, s4=%2g, s5=%2g, s6=%2g',vSpin));
H2=reshape(vSpin(1:4),2,2);
H2*H2',

if length(vE)>1
    imagesc(1);
    plot(vE);
end;

