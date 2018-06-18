function val_q= s2q(val_s)
% transform s={-1,+1} -> q={1,0}
% usage: (1) define symbolic variable x
%             syms val_s;
%        (2) fq=s2q(val_s); 
val_q=(1/2)*(1-val_s);
