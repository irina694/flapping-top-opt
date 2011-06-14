function[vec] = calc_panel_vectors_new(Xr,Yr,Zr,Xc,Yc,Zc,t,t_step,M,N,i)

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
outward_x = outward(:,1); outward_y = outward(:,2); outward_z = outward(:,3);

tau_x = (.5*(Xr(2:M+1,1:N)+Xr(2:M+1,2:N+1))-.5*(Xr(1:M,1:N)+Xr(1:M,2:N+1)))';
tau_y = (.5*(Yr(2:M+1,1:N)+Yr(2:M+1,2:N+1))-.5*(Yr(1:M,1:N)+Yr(1:M,2:N+1)))';
tau_z = (.5*(Zr(2:M+1,1:N)+Zr(2:M+1,2:N+1))-.5*(Zr(1:M,1:N)+Zr(1:M,2:N+1)))';
tau = [tau_x(:),tau_y(:),tau_z(:)]./repmat(sqrt(tau_x(:).^2+tau_y(:).^2+tau_z(:).^2),1,3);

V_h = U_kin.*tau(:,1)+V_kin.*tau(:,2)+W_kin.*tau(:,3);
V_v = U_kin.*outward(:,1)+V_kin.*outward(:,2)+W_kin.*outward(:,3);
alpha = atan2(V_v,V_h)*180/pi;

j1 = [tau(:,1).^3.*cos(alpha*pi/180)+tau(:,1).*outward(:,1).^2.*cos(alpha*pi/180)+tau(:,2).^2.*tau(:,1).*cos(alpha*pi/180)+tau(:,2).^2.*outward(:,1).*sin(alpha*pi/180)-tau(:,2).*outward(:,2).*tau(:,1).*sin(alpha*pi/180)+tau(:,2).*outward(:,2).*outward(:,1).*cos(alpha*pi/180)+tau(:,3).^2.*tau(:,1).*cos(alpha*pi/180)+tau(:,3).^2.*outward(:,1).*sin(alpha*pi/180)-tau(:,3).*outward(:,3).*tau(:,1).*sin(alpha*pi/180)+tau(:,3).*outward(:,3).*outward(:,1).*cos(alpha*pi/180),...
    tau(:,2).*tau(:,1).^2.*cos(alpha*pi/180)+outward(:,2).*tau(:,1).^2.*sin(alpha*pi/180)-tau(:,1).*tau(:,2).*outward(:,1).*sin(alpha*pi/180)+tau(:,1).*outward(:,2).*outward(:,1).*cos(alpha*pi/180)+tau(:,2).^3.*cos(alpha*pi/180)+tau(:,2).*outward(:,2).^2.*cos(alpha*pi/180)+tau(:,3).^2.*tau(:,2).*cos(alpha*pi/180)+tau(:,3).^2.*outward(:,2).*sin(alpha*pi/180)-tau(:,3).*outward(:,3).*tau(:,2).*sin(alpha*pi/180)+tau(:,3).*outward(:,3).*outward(:,2).*cos(alpha*pi/180),...
    tau(:,3).*tau(:,1).^2.*cos(alpha*pi/180)+outward(:,3).*tau(:,1).^2.*sin(alpha*pi/180)-tau(:,1).*tau(:,3).*outward(:,1).*sin(alpha*pi/180)+tau(:,1).*outward(:,3).*outward(:,1).*cos(alpha*pi/180)+tau(:,3).*tau(:,2).^2.*cos(alpha*pi/180)+outward(:,3).*tau(:,2).^2.*sin(alpha*pi/180)-tau(:,2).*tau(:,3).*outward(:,2).*sin(alpha*pi/180)+tau(:,2).*outward(:,3).*outward(:,2).*cos(alpha*pi/180)+tau(:,3).^3.*cos(alpha*pi/180)+tau(:,3).*outward(:,3).^2.*cos(alpha*pi/180)];

j2 = [outward(:,1).*tau(:,1).^2.*cos(alpha*pi/180)+outward(:,1).^3.*cos(alpha*pi/180)+outward(:,2).*tau(:,2).*tau(:,1).*cos(alpha*pi/180)+outward(:,2).*tau(:,2).*outward(:,1).*sin(alpha*pi/180)-outward(:,2).^2.*tau(:,1).*sin(alpha*pi/180)+outward(:,2).^2.*outward(:,1).*cos(alpha*pi/180)+outward(:,3).*tau(:,3).*tau(:,1).*cos(alpha*pi/180)+outward(:,3).*tau(:,3).*outward(:,1).*sin(alpha*pi/180)-outward(:,3).^2.*tau(:,1).*sin(alpha*pi/180)+outward(:,3).^2.*outward(:,1).*cos(alpha*pi/180),...
    outward(:,1).*tau(:,2).*tau(:,1).*cos(alpha*pi/180)+outward(:,1).*outward(:,2).*tau(:,1).*sin(alpha*pi/180)-tau(:,2).*outward(:,1).^2.*sin(alpha*pi/180)+outward(:,2).*outward(:,1).^2.*cos(alpha*pi/180)+outward(:,2).*tau(:,2).^2.*cos(alpha*pi/180)+outward(:,2).^3.*cos(alpha*pi/180)+outward(:,3).*tau(:,3).*tau(:,2).*cos(alpha*pi/180)+outward(:,3).*tau(:,3).*outward(:,2).*sin(alpha*pi/180)-outward(:,3).^2.*tau(:,2).*sin(alpha*pi/180)+outward(:,3).^2.*outward(:,2).*cos(alpha*pi/180),...
    outward(:,1).*tau(:,3).*tau(:,1).*cos(alpha*pi/180)+outward(:,1).*outward(:,3).*tau(:,1).*sin(alpha*pi/180)-tau(:,3).*outward(:,1).^2.*sin(alpha*pi/180)+outward(:,3).*outward(:,1).^2.*cos(alpha*pi/180)+outward(:,2).*tau(:,3).*tau(:,2).*cos(alpha*pi/180)+outward(:,2).*outward(:,3).*tau(:,2).*sin(alpha*pi/180)-tau(:,3).*outward(:,2).^2.*sin(alpha*pi/180)+outward(:,3).*outward(:,2).^2.*cos(alpha*pi/180)+outward(:,3).*tau(:,3).^2.*cos(alpha*pi/180)+outward(:,3).^3.*cos(alpha*pi/180)];

n_lift_x = j2(:,1);
n_lift_y = j2(:,2);
n_lift_z = j2(:,3);
n_drag_x = j1(:,1);
n_drag_y = j1(:,2);
n_drag_z = j1(:,3);

vec.outward_x = outward_x;
vec.outward_y = outward_y;
vec.outward_z = outward_z;
vec.alpha = alpha;
vec.n_lift_x = n_lift_x;
vec.n_lift_y = n_lift_y;
vec.n_lift_z = n_lift_z;
vec.n_drag_x = n_drag_x;
vec.n_drag_y = n_drag_y;
vec.n_drag_z = n_drag_z;
