function[circ,circ_local,w_ind] = solve_system(C1,C1_drag,C2,C2_drag,LL,i,circ_wake,N,M)

if i == 1
    circ = C1\LL;
    w_ind = C1_drag*circ;
else
    circ = C1\(LL-C2*circ_wake);
    w_ind = C1_drag*circ + C2_drag*circ_wake;
end

foo = reshape(circ,N,M)';
foo2 = foo;
foo2(2:M,:) = foo(2:M,:)-foo(1:M-1,:);
circ_local = reshape(foo2',M*N,1);