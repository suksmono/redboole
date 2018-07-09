% CHECK_H2_sqs_new.m
% solve H2 with quantum ising machine
% use symbolic computation for Boolean reduction
% s={-1,1}; q={0,1}
% keep in mind that boolean reduction works on q-domain
% case: H4=
% |<-problems=M ->| <---------- ancillas=C(M,2)=M(M-1)/2 --------------->|
%[q11 q12 ||q13(1,2)
% q21 q22 ||q23(1,2)

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
% define h-order
M=2;
%
syms q11 q12 q13;  
syms q21 q22 q23; 
% store as matrix of symbols
q=[ q11 q12 q13; ...
    q21 q22 q23  ...
    ];

% ------------------------------------------------------------------------
% 1. Formulate Hamiltonian Hks
% ==> d12=s1*s3 + s2*s4 ;
% ------------------------------------------------------------------------
syms s11 s12 s13;
syms s21 s22 s23;

s=[ s11 s12 s13; ...
    s21 s22 s23  ...
    ];
% ------------------------------------------------------------------------
% calculation of d12, symboliccally 
% ------------------------------------------------------------------------
d12=transpose(s(:,1))*s(:,2); 
d12_2=expand(d12^2);

% ------------------------------------------------------------------------
%change all si^2 to 1
% ------------------------------------------------------------------------
Hks_12=simplify( ...
                 subs(d12_2, ...
                 { s11^2, s12^2, ...
                   s21^2, s22^2 ...
                   }, ...
                  { 1, 1, 1, 1} ...
                 ) ...
             );            
%
Hks=Hks_12;
disp(sprintf('Hks: \n%s',char(Hks)));
 
% ------------------------------------------------------------------------
% 2. Transform to q-domain: Hks-> Hkq
% ------------------------------------------------------------------------
Hkq=expand(simplify( ...
                subs(Hks, ...
                { s11, s12, ...
                  s21, s22  ...
                  }, ...
                { s2q(q11), s2q(q12), ...
                  s2q(q21), s2q(q22)  ...
                  }...
                )...
            ));
disp(sprintf('Hkq: \n%s',char(Hkq)));
% ------------------------------------------------------------------------

% ------------------------------------------------------------------------
% 3. Boolean reduction: Hkq->H2q
% ------------------------------------------------------------------------
%
d=2*4;
H2q=simplify( ...
     subs( Hkq, ...
          { q(1,1)*q(1,2), q(2,1)*q(2,2) ...
            }, ...
          { q(1,3), q(2,3) ...
            } ...
        ) ...
         + H2body(q(1,1),q(1,2),q(1,3),d) ... %(1,2) v
         + H2body(q(2,1),q(2,2),q(2,3),d) ...
    );
% DISPLAY RESULTS
disp(sprintf('H2q: \n%s',char(expand(H2q))));

% ------------------------------------------------------------------------
% 4. Transform-back to s-domain: H2q->H2s
% ------------------------------------------------------------------------
H2s = simplify( ...
            subs(H2q, ...
                { q11, q12, q13, ... 
                  q21, q22, q23  ...
                  }, ...
                { q2s(s11), q2s(s12),q2s(s13), ...
                  q2s(s21), q2s(s22),q2s(s23)  ...
                  }...
                )...
            );
disp(sprintf('H2s: \n%s', char(expand(H2s))));
% ------------------------------------------------------------
% H2s = 
%   4*s11 + 4*s12 - 8*s13 + 4*s21 + 4*s22 - 8*s23 + 2*s11*s12 
% - 4*s11*s13 - 4*s12*s13 + 2*s11*s21 + 2*s11*s22 + 2*s12*s21 
% - 4*s11*s23 + 2*s12*s22 - 4*s13*s21 - 4*s12*s23 - 4*s13*s22 
% + 8*s13*s23 + 2*s21*s22 - 4*s21*s23 - 4*s22*s23 + 16
% ------------------------------------------------------------
 
% ------------------------------------------------------------
%  EXTRACT THE RESULTS FROM SYMBOLIC SOLUTION
% ------------------------------------------------------------
% now, we are working on orthogonalization of 2-binary vectors
% q11 ... q41, q12 .... q42, q15, ..., q45 =? 4*3=12 qubits 
[cx,tx]=coeffs(H2s,[ s11 s12 s13 ...
                     s21 s22 s23 ...
                     ]);
%                      s110 s210 s310 s410 ...
% number of  qubits: q11, ..., q45 => 7*4=28
NQ=2*3;

%% create lookup tbl for product of s_ij*s_mn
sx=[transpose(s(:,1)) transpose(s(:,2)) transpose(s(:,3)) ...
    ];
%%
tblook=transpose(sx)*sx;
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
disp(sprintf('calculating hi ...'));
for m=1:NQ
    tsym=sx(m);
    for nn=1:NTERM
        if(tx(nn)==tsym)
            hi(m)=cx(nn);
        end;
    end;
end;
%% fill in Jij, we only consider upper diagonal of Jij
disp(sprintf('calculating Jij ...'));
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
maxIter=1*10*1000;
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
    disp(sprintf('LQ%d: Iter->%d of %d, Tresh=%1.4f, E=%4.2f',...
        NQ,k,maxIter, Ptresh(k),double(Esys)));
end;
% ------------------------------------------------------------
% DISPLAY THE RESULTS
% ------------------------------------------------------------
v1=vSpin(1:M); v2=vSpin(M+1:2*M); 
disp(sprintf('\ns11= %2g, s21= %2g',v1));
disp(sprintf('s12= %2g, s22= %2g',v2));
disp(sprintf('<v1,v2> = %2g', v1'*v2));

% plot annealing curve
if length(vE)>1
    imagesc(1);
    plot(vE);
end;

