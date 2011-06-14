function[Xr,Yr,Zr,Xc,Yc,Zc,Q] = calc_grids(x,y,z,roll,pitch,yaw,path_X,path_Y,path_Z,t,M,N,i)

Qr_xx = sparse((N+1)*(M+1),(N+1)*(M+1));
Qr_yy = sparse((N+1)*(M+1),(N+1)*(M+1));
Qr_zz = sparse((N+1)*(M+1),(N+1)*(M+1));

Qc_xx = sparse(N*M,(N+1)*(M+1));
Qc_yy = sparse(N*M,(N+1)*(M+1));
Qc_zz = sparse(N*M,(N+1)*(M+1));

for n = 1:N+1
   for m = 1:M  
      
      Qr_xx((m-1)*(N+1)+n,(m+1-1)*(N+1)+n) = 1/4;
      Qr_xx((m-1)*(N+1)+n,(m-1)*(N+1)+n) = 3/4;
      Qr_yy((m-1)*(N+1)+n,(m-1)*(N+1)+n) = 1;
      Qr_zz((m-1)*(N+1)+n,(m+1-1)*(N+1)+n) = 1/4;
      Qr_zz((m-1)*(N+1)+n,(m-1)*(N+1)+n) = 3/4;
      
   end
   
   Qr_xx((m+1-1)*(N+1)+n,(m+1-1)*(N+1)+n) = 5/4;
   Qr_xx((m+1-1)*(N+1)+n,(m-1)*(N+1)+n) = -1/4;
   Qr_yy((m+1-1)*(N+1)+n,(m-1)*(N+1)+n) = 1;
   Qr_zz((m+1-1)*(N+1)+n,(m+1-1)*(N+1)+n) = 1;
   
end

for n = 1:N
   for m = 1:M      
      
      Qc_xx((m-1)*N+n,(m-1)*(N+1)+n) = 1/4;
      Qc_xx((m-1)*N+n,(m+1-1)*(N+1)+n) = 1/4;
      Qc_xx((m-1)*N+n,(m-1)*(N+1)+n+1) = 1/4;
      Qc_xx((m-1)*N+n,(m+1-1)*(N+1)+n+1) = 1/4;
      Qc_yy((m-1)*N+n,(m-1)*(N+1)+n) = 1/4;
      Qc_yy((m-1)*N+n,(m+1-1)*(N+1)+n) = 1/4;
      Qc_yy((m-1)*N+n,(m-1)*(N+1)+n+1) = 1/4;
      Qc_yy((m-1)*N+n,(m+1-1)*(N+1)+n+1) = 1/4;
      Qc_zz((m-1)*N+n,(m-1)*(N+1)+n) = 1/4;
      Qc_zz((m-1)*N+n,(m+1-1)*(N+1)+n) = 1/4;
      Qc_zz((m-1)*N+n,(m-1)*(N+1)+n+1) = 1/4;
      Qc_zz((m-1)*N+n,(m+1-1)*(N+1)+n+1) = 1/4;  
      
   end
end

Qc_xx = Qc_xx*Qr_xx;
Qc_yy = Qc_yy*Qr_yy;
Qc_zz = Qc_zz*Qr_zz;

T_roll = [1,0,0;0,cos(roll*pi/180),-sin(roll*pi/180);0,sin(roll*pi/180),cos(roll*pi/180)];
T_pitch = [cos(pitch*pi/180),0,sin(pitch*pi/180);0,1,0;-sin(pitch*pi/180),0,cos(pitch*pi/180)];
T_yaw = [cos(yaw*pi/180),-sin(yaw*pi/180),0;sin(yaw*pi/180),cos(yaw*pi/180),0;0,0,1];
T = T_yaw*T_pitch*T_roll;  

Qr_xx = Qr_xx*T(1,1);
Qr_xy = Qr_yy*T(1,2);
Qr_xz = Qr_zz*T(1,3);
Qr_yx = Qr_xx*T(2,1)/T(1,1);
Qr_yy = Qr_yy*T(2,2);
Qr_yz = Qr_zz*T(2,3);
Qr_zx = Qr_xx*T(3,1)/T(1,1);
Qr_zy = Qr_yy*T(3,2)/T(2,2);
Qr_zz = Qr_zz*T(3,3);

Qc_xx = Qc_xx*T(1,1);
Qc_xy = Qc_yy*T(1,2);
Qc_xz = Qc_zz*T(1,3);
Qc_yx = Qc_xx*T(2,1)/T(1,1);
Qc_yy = Qc_yy*T(2,2);
Qc_yz = Qc_zz*T(2,3);
Qc_zx = Qc_xx*T(3,1)/T(1,1);
Qc_zy = Qc_yy*T(3,2)/T(2,2);
Qc_zz = Qc_zz*T(3,3);

% for m = 1:M+1
%    for n = 1:N+1
%       
%       dummy = reshape(Qr_yy*y(:,1),N+1,M+1)'/T(2,2);
%       if dummy(m,n) == 0
%          Qr_yx((m-1)*(N+1)+n,:) = 0; % prevents root of cambered wing from crossing to negative Y..
%          Qr_yy((m-1)*(N+1)+n,:) = 0;
%          Qr_yz((m-1)*(N+1)+n,:) = 0;
%       end
%       
%    end
% end

Xr(:,1) = Qr_xx*x+Qr_xy*y+Qr_xz*z+path_X;
Yr(:,1) = Qr_yx*x+Qr_yy*y+Qr_yz*z+path_Y;
Zr(:,1) = Qr_zx*x+Qr_zy*y+Qr_zz*z+path_Z;
Xc(:,1) = Qc_xx*x+Qc_xy*y+Qc_xz*z+path_X;
Yc(:,1) = Qc_yx*x+Qc_yy*y+Qc_yz*z+path_Y;
Zc(:,1) = Qc_zx*x+Qc_zy*y+Qc_zz*z+path_Z;

Q.Qr = [Qr_xx,Qr_xy,Qr_xz;Qr_yx,Qr_yy,Qr_yz;Qr_zx,Qr_zy,Qr_zz];
Q.Qc = [Qc_xx,Qc_xy,Qc_xz;Qc_yx,Qc_yy,Qc_yz;Qc_zx,Qc_zy,Qc_zz];
