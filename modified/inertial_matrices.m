function[T,ang_vel,ang_acc,Xr_dot_dot,X,Y,Z,omega_body] = inertial_matrices(t,omega,roll,pitch,yaw,flap_amp,plunge_amp,path_X,path_Y,path_Z,Xo,Yo,Zo,N_nodes)

T = zeros(3,3,length(t)); 
ang_vel = zeros(3,3,length(t)); 
ang_acc = zeros(3,3,length(t));
Xr_dot_dot = zeros(length(t),3);
omega_body = zeros(3,length(t));
for i = 1:length(t)

    T_roll = [1,0,0;0,cos(roll(i)*pi/180),-sin(roll(i)*pi/180);0,sin(roll(i)*pi/180),cos(roll(i)*pi/180)];
    T_pitch = [cos(pitch(i)*pi/180),0,sin(pitch(i)*pi/180);0,1,0;-sin(pitch(i)*pi/180),0,cos(pitch(i)*pi/180)];
    T_yaw = [cos(yaw(i)*pi/180),-sin(yaw(i)*pi/180),0;sin(yaw(i)*pi/180),cos(yaw(i)*pi/180),0;0,0,1];
    T(:,:,i) = T_yaw*T_pitch*T_roll;

    ang_vel_roll = -(flap_amp*pi/180)*sin(omega*t(i)+2*pi)*omega*[0,0,0;0,0,-1;0,1,0];
    ang_vel_pitch = [0,0,0;0,0,0;0,0,0];
    ang_vel_yaw = [0,0,0;0,0,0;0,0,0];
    ang_vel(:,:,i) = ang_vel_yaw+T_yaw*ang_vel_pitch*T_yaw'+T_yaw*T_pitch*ang_vel_roll*T_pitch'*T_yaw';

    omega_body(:,i) = [-(flap_amp*pi/180)*sin(omega*t(i)+2*pi)*omega;0;0];
    
    ang_acc_roll = (-flap_amp*pi/180)*(omega^2)*cos(omega*t(i)+2*pi)*[0,0,0;0,0,-1;0,1,0];
    ang_acc_pitch = [0,0,0;0,0,0;0,0,0];
    ang_acc_yaw = [0,0,0;0,0,0;0,0,0];
    ang_acc(:,:,i) = ang_acc_yaw + ang_vel_yaw*T_yaw*ang_vel_pitch*T_yaw' + T_yaw*ang_acc_pitch*T_yaw' + T_yaw*ang_vel_pitch*T_yaw'*ang_vel_yaw' + ...
        ang_vel_yaw*T_yaw*T_pitch*ang_vel_roll*T_pitch'*T_yaw' + T_yaw*ang_vel_pitch*T_pitch*ang_vel_roll*T_pitch'*T_yaw' + T_yaw*T_pitch*ang_acc_roll*T_pitch'*T_yaw' + ...
        T_yaw*T_pitch*ang_vel_roll*T_pitch'*ang_vel_pitch'*T_yaw' + T_yaw*T_pitch*ang_vel_roll*T_pitch'*T_yaw'*ang_vel_yaw';

    Xr_dot_dot(i,1) = 0;
    Xr_dot_dot(i,2) = 0;
    Xr_dot_dot(i,3) = -plunge_amp*omega*omega*cos(omega*t(i)+2*pi);

    dummy = T(:,:,i)*[Xo';Yo';Zo']+repmat([path_X(i);path_Y(i);path_Z(i)],1,N_nodes);
    X(:,i) = dummy(1,:)';
    Y(:,i) = dummy(2,:)';
    Z(:,i) = dummy(3,:)';

end