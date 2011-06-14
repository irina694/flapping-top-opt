function[CL,CT,power,pressure] = compute_forces_new(circ,circ_local,w_ind,t_step,x,y,z,Xc,Yc,Zc,vec,rho,b,c,V_ref,i,M,N)

if i == 1
    U_kin = zeros(M*N,1);
    V_kin = zeros(M*N,1);
    W_kin = zeros(M*N,1);
else
    U_kin = -(Xc(:,i)-Xc(:,i-1))/t_step;
    V_kin = -(Yc(:,i)-Yc(:,i-1))/t_step;
    W_kin = -(Zc(:,i)-Zc(:,i-1))/t_step;
end

A = sparse(N*M,N*M);
for q = 1:N*M
    A(q,q) = 1;
end
for q = 1:N*M-N
    A(q+N,q) = 1;
end
A = A/2/t_step;

if i == 1
    dfdt = A*circ(:,i);
else
    dfdt = A*circ(:,i) - A*circ(:,i-1);
end

x = reshape(x(:,i),N+1,M+1)';
y = reshape(y(:,i),N+1,M+1)';
z = reshape(z(:,i),N+1,M+1)';

chord = reshape(sqrt((.5*(x(2:M+1,1:N)+x(2:M+1,2:N+1))-.5*(x(1:M,1:N)+x(1:M,2:N+1))).^2+(.5*(y(2:M+1,1:N)+y(2:M+1,2:N+1))-.5*(y(1:M,1:N)+y(1:M,2:N+1))).^2+(.5*(z(2:M+1,1:N)+z(2:M+1,2:N+1))-.5*(z(1:M,1:N)+z(1:M,2:N+1))).^2)',M*N,1);
span = reshape(sqrt((.5*(x(1:M,2:N+1)+x(2:M+1,2:N+1))-.5*(x(1:M,1:N)+x(2:M+1,1:N))).^2+(.5*(y(1:M,2:N+1)+y(2:M+1,2:N+1))-.5*(y(1:M,1:N)+y(2:M+1,1:N))).^2+(.5*(z(1:M,2:N+1)+z(2:M+1,2:N+1))-.5*(z(1:M,1:N)+z(2:M+1,1:N))).^2)',M*N,1);

dL = rho*span.*(sqrt(U_kin.^2+V_kin.^2+W_kin.^2).*circ_local(:,i)+chord.*dfdt).*cos(vec.alpha*pi/180);
dD = rho*span.*(-w_ind(:,i).*circ_local(:,i)+chord.*dfdt.*sin(vec.alpha*pi/180));

dL_vec = [dL.*vec.n_lift_x,dL.*vec.n_lift_y,dL.*vec.n_lift_z];
dD_vec = [dD.*vec.n_drag_x,dD.*vec.n_drag_y,dD.*vec.n_drag_z];

F_X = dL_vec(:,1)+dD_vec(:,1);
F_Z = dL_vec(:,3)+dD_vec(:,3);

pressure = dL./chord./span;

power = sum(2*pressure.*(U_kin.*vec.outward_x+V_kin.*vec.outward_y+W_kin.*vec.outward_z).*chord.*span);
    
CT = -sum(sum(F_X*2))/(b*c*.5*rho*V_ref^2);
CL = sum(sum(F_Z*2))/(b*c*.5*rho*V_ref^2);

