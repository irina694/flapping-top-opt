function[X_bar_wake,Y_bar_wake,Z_bar_wake,B1,B2] = wake_shift(Xr_i,Yr_i,Zr_i,Xwake_i1,Ywake_i1,Zwake_i1,i,N,M)

B1 = sparse((N+1)*i,(N+1)*(M+1));
B1((N+1)*i-N:(N+1)*i,(N+1)*(M+1)-N:(N+1)*(M+1)) = eye(N+1,N+1);

if i == 1
    B2 = 0;
else
    B2 = sparse((N+1)*i,(N+1)*(i-1));
    B2(1:(N+1)*(i-1),1:(N+1)*(i-1)) = eye((N+1)*(i-1));
end

X_bar_wake = B1*Xr_i + B2*Xwake_i1;
Y_bar_wake = B1*Yr_i + B2*Ywake_i1;
Z_bar_wake = B1*Zr_i + B2*Zwake_i1;