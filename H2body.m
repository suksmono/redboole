function H2sub= H2body(x1,x2,y,d12)
% H_5(q1,q2,Q5;d5)=d5*(3*Q5+q1*q2-2*q1*Q5-2*q2*q5)
H2sub=d12*(3*y+x1*x2-2*x1*y-2*x2*y);
