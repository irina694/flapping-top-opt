function[LL] = calc_LL_new(Xr,Yr,Zr,Xc,Yc,Zc,t_step,M,N,i)

if i == 1
    U_kin = zeros(M*N,1); V_kin = zeros(M*N,1); W_kin = zeros(M*N,1);
else
    U_kin = -(Xc(:,i)-Xc(:,i-1))/t_step;
    V_kin = -(Yc(:,i)-Yc(:,i-1))/t_step;
    W_kin = -(Zc(:,i)-Zc(:,i-1))/t_step;
end

Xr = reshape(Xr(:,i),N+1,M+1)';
Yr = reshape(Yr(:,i),N+1,M+1)';
Zr = reshape(Zr(:,i),N+1,M+1)';

v1_x = Xr(1:M,1:N)-Xr(2:M+1,2:N+1); v1_y = Yr(1:M,1:N)-Yr(2:M+1,2:N+1); v1_z = Zr(1:M,1:N)-Zr(2:M+1,2:N+1);
v2_x = Xr(2:M+1,1:N)-Xr(1:M,2:N+1); v2_y = Yr(2:M+1,1:N)-Yr(1:M,2:N+1); v2_z = Zr(2:M+1,1:N)-Zr(1:M,2:N+1);
outward_x = (v1_y.*v2_z-v2_y.*v1_z)'; outward_y = (v2_x.*v1_z-v1_x.*v2_z)'; outward_z = (v1_x.*v2_y-v2_x.*v1_y)';
outward = [outward_x(:),outward_y(:),outward_z(:)]./repmat(sqrt(outward_x(:).^2+outward_y(:).^2+outward_z(:).^2),1,3);

LL = -outward(:,1).*U_kin - outward(:,2).*V_kin - outward(:,3).*W_kin;
