function[circ_wake,A1,A2] = circ_wake_shift(circ_wake,circ_i,i,N,M)

if i < 3
    A1 = 0;
else
    A1 = sparse(N*(i-1),N*(i-2));
    A1(1:N*(i-2),1:N*(i-2)) = eye(N*(i-2));
end

A2 = sparse(N*(i-1),N*M);
A2(N*(i-1)-N+1:N*(i-1),N*M-N+1:N*M) = eye(N,N);

circ_wake = A1*circ_wake + A2*circ_i;