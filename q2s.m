function val_s= q2s(val_q)
% transform s={0,1} -> q={1,-1}
% usage: (1) define symbolic variable x
%             syms val_s;
%        (2) fq=1-2*val_q; 
val_s=1-2*val_q;
