function H2sub= H2body(x1,x2,y,d12)
% H_k(qi,qj,qk;dk)=dk*(3*qk+qi*qj-2*qi*qk-2*qj*qk)
H2sub=d12*(3*y+x1*x2-2*x1*y-2*x2*y);
