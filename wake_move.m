function[Xwake,Ywake,Zwake,E] = wake_move(X_bar_wake,Y_bar_wake,Z_bar_wake,circ,circ_wake,D1x,D1y,D1z,D2x,D2y,D2z,t_step,N,i,N_w)

if i-1 < N_w
    wake_start = 1;
else
    wake_start = i-N_w;
end

del_X = D1x*circ*t_step + D2x*circ_wake*t_step;
del_Y = D1y*circ*t_step + D2y*circ_wake*t_step;
del_Z = D1z*circ*t_step + D2z*circ_wake*t_step;

E = sparse((N+1)*i,(N+1)*(i-1));
E(1:(N+1)*(i-1),1:(N+1)*(i-1)) = eye((N+1)*(i-1));

Xwake = X_bar_wake + E*del_X;
Ywake = Y_bar_wake + E*del_Y;
Zwake = Z_bar_wake + E*del_Z;
