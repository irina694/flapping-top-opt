% evaluate_fitness:
%   Input: xx is a vector of size 1-by-Nvariables 
%   Output: CL_ave, CT_ave, CP_ave, eff, mass, peak_disp

% Sets up wing geometry and kinematics parameters, number of time steps,
% and error tolerance. Computes and displays the wing geometry at each time
% step.

function [CL_ave,CT_ave,CP_ave,eff,mass,peak_disp] = evaluate_fitness(xx)

% retrieve stored variables 
load globalVars.mat globalVars 
Nlevels = globalVars.Gvars.Nlevels;

% error flag
error = 0;

try
    lXsmap = Gvars.NinitEdges;
    Xsmap = xx(1:lXsmap);

    [startMap,startDir] = extractsSmap(Xsmap);
    map = initMapStruct(startMap,startDir,input);

    [map.Prules.letters,map.Prules.direction,map.Prules.thickRule,...
        map.Prules.elastRule,map.Prules.pressRule] = extractRules(xx(lXsmap+1:end-Gvars.nB));

    map.cells.vertices.type = 'FFFF';

    map = mapLS(map,Nlevels);

    connectivity=removeVertices(map);
    connectivity.tip = xx(end)/100;

catch
    Info = 'Map L-system error';
    error = 1;
end

%% Wing geometry

M = 13;  % # of chordwise panels
N = 26;  % # of spanwise panels

c = .05;  % chord
b = .2;  % span

cutoff = 1E-9;  % cutoff radius where the Biot-Savart does not act

%% Wing Structure:

E_modulus = 40E9;
poisson = 0.3;
rho_wing = 1400;
damping = 20;
t_min = 0.00005;
t_max = 0.0005;

%% External flow conditions:

V_ref = 5;  % U velocity, and reference for some terms
rho = 1.225;

%% Specify kinematics:

omega = 50;  % frequency, rad/s

flap_amp = 30;
plunge_amp = 0;
pitch_angle = 10;

%% Time steps:

N_per = 50;
N_cycles = 3;
N_w = 10;  % number of deforming wake rows

N_steps = N_per*N_cycles;
t_step = 2*pi/omega/N_per;
t = 0:t_step:t_step*N_steps;
size(t)
N_Q_begin = max(N_per*(N_cycles-1)+1,1);

%% Specify kinematics:

path_X = -V_ref*t;
path_Y = 0*t;
path_Z = plunge_amp*cos(omega*t);
roll = flap_amp*cos(omega*t);
pitch = 0*t+pitch_angle;
yaw = 0*t;

%% Generalized alpha parameters

p_inf = .8;

alp_m = (2*p_inf-1)/(p_inf+1);
alp_f = p_inf/(p_inf+1);
bet = .25*(1-alp_m+alp_f)^2;
gam = .5-alp_m+alp_f;

%% Error tolerances

aero_tol = -3;

%% Compute wing grid:

[Xo,Yo,Zo,N_nodes,trielements,N_elements,GlobalDOF,DOF,asmb,bc,fdof] = wing_grid(c,b,M,N);

%% transform the bioTom topology into SIMP topology

connectivity.vertices.coords = [connectivity.vertices.coords(:,2),connectivity.vertices.coords(:,1)];
connectivity.vertices.coords(:,1) = c*(connectivity.vertices.coords(:,1)-1)/-2;
connectivity.vertices.coords(:,2) = b*(connectivity.vertices.coords(:,2)+1)/4;

[X] = bio2simp(connectivity,Xo,Yo,N_elements,trielements);
thickness = 0*X; 
thickness(X==0) = t_min;
thickness(X==1) = t_max;

mass = rho_wing*b*c/N_elements*sum(thickness);

%% Compute inertial matrices:

 [T,ang_vel,ang_acc,Xr_dot_dot,X,Y,Z,omega_body] = inertial_matrices(t,omega,roll,pitch,yaw,flap_amp,plunge_amp,path_X,path_Y,path_Z,Xo,Yo,Zo,N_nodes);
 
%% Initialize terms

Xwake = 0;
Ywake = 0;
Zwake = 0;
circ_wake = 0;

answer = zeros(N_nodes,6);

%% Mass and damping matrix

[C_mat,M_mat] = CM_find_new(Xo,Yo,Zo,trielements,N_elements,asmb,GlobalDOF,fdof,rho_wing,thickness,damping);

%% Stiffness matrix

K = K_find_new(Xo,Yo,Zo,trielements,N_elements,asmb,GlobalDOF,fdof,E_modulus,poisson,thickness);

%% VLM to FEA interpolation matrix

Q_interp1 = kron(speye(M*N),ones(2,1)); % kron computes Kronecker tensor product

%% Start time steps

figure
set(gcf,'position',[178         215        1287         883])

% for i = 1:length(t)
% 
%     %% Compute inertial forces
%     
%     [F_iner,F_iner_full] = F_find_new(Xo,Yo,Zo,trielements,N_elements,asmb,GlobalDOF,fdof,rho_wing,thickness,T(:,:,i),ang_vel(:,:,i),ang_acc(:,:,i),Xr_dot_dot(i,:));
% 
%     %% Update wake circulations
% 
%     if i > 1
%         [circ_wake,A1,A2] = circ_wake_shift(circ_wake,circ(:,i-1),i,N,M);
%     end
%     
%     %% Aeroelastic loop:
% 
%     aero_err = 0; qq = 0; log = [0,0,0];
%     while aero_err > aero_tol || qq < 3
% 
%         qq = qq + 1;
%         
%         %% Update wing shape in local coordinates
% 
%         x(:,i) = Xo + answer(:,1);
%         y(:,i) = Yo + answer(:,2);
%         z(:,i) = Zo + answer(:,3);
% 
%         %% Compute wing shape in global coordinates
% 
%         [Xr(:,i),Yr(:,i),Zr(:,i),Xc(:,i),Yc(:,i),Zc(:,i),Q(i)] = calc_grids(x(:,i),y(:,i),z(:,i),roll(i),pitch(i),yaw(i),path_X(i),path_Y(i),path_Z(i),t,M,N,i);
% 
%          %% Compute source terms L
% 
%         [LL] = calc_LL_new(Xr,Yr,Zr,Xc,Yc,Zc,t_step,M,N,i);
% 
%         %% Compute panel vectors
% 
%         [vec] = calc_panel_vectors_new(Xr,Yr,Zr,Xc,Yc,Zc,t,t_step,M,N,i);
%         
%         %% Add on to wake
% 
%         [X_bar_wake,Y_bar_wake,Z_bar_wake,B1,B2] = wake_shift(Xr(:,i),Yr(:,i),Zr(:,i),Xwake,Ywake,Zwake,i,N,M);
% 
%         %% Wing on wing influence matrix
% 
%         [C1,C1_drag] = calc_C1_new(Xr(:,i),Yr(:,i),Zr(:,i),Xc(:,i),Yc(:,i),Zc(:,i),M,N,vec);
% 
%         %% Wake on wing influence matrix
% 
%         if i > 1
%             [C2,C2_drag] = calc_C2_new(X_bar_wake,Y_bar_wake,Z_bar_wake,Xc(:,i),Yc(:,i),Zc(:,i),M,N,vec,i);
%         else
%             C2 = 0; C2_drag = 0;
%         end
% 
%         %% Solve system of equations
% 
%         [circ(:,i),circ_local(:,i),w_ind(:,i)] = solve_system(C1,C1_drag,C2,C2_drag,LL,i,circ_wake,N,M);
% 
%         %% Compute forces
% 
%         if qq > 1
%             pressure_old = pressure(:,i);
%         end        
% 
%         [CL(i),CT(i),aero_power(i),pressure(:,i)] = compute_forces_new(circ,circ_local,w_ind,t_step,x,y,z,Xc,Yc,Zc,vec,rho,b,c,V_ref,i,M,N);
%         
%         if qq > 1
%             aero_err = log10(norm(pressure(:,i)-pressure_old));
%         end
%         
%         pressure_FEA = Q_interp1*pressure(:,i);
% 
%         if i == 1
% 
%             D_n1 = zeros(GlobalDOF-length(bc),1);
%             Dp_n1 = zeros(GlobalDOF-length(bc),1);
% 
%             F_aero = Faero_find_new(Xo,Yo,Zo,trielements,N_elements,asmb,GlobalDOF,fdof,pressure_FEA);
%             F = F_iner + F_aero;
% 
%             Dpp_n1 = M_mat\(F - K*D_n1 - C_mat*Dp_n1);
% 
%         end
%         
%         if i > 1
% 
%             if qq == 1
%                 F_old = F; D_n = D_n1; Dp_n = Dp_n1; Dpp_n = Dpp_n1;
%             end
%             
%             F_aero = Faero_find_new(Xo,Yo,Zo,trielements,N_elements,asmb,GlobalDOF,fdof,pressure_FEA);
%             F = F_iner + F_aero;
% 
%             Jac = ((1-alp_m)/bet/t_step/t_step)*M_mat + (gam*(1-alp_f)/bet/t_step)*C_mat + (1-alp_f)*K;
%             Res = (1-alp_f)*F + alp_f*F_old - alp_f*K*D_n - alp_f*C_mat*Dp_n + ...
%                 (1-alp_f)*C_mat*((gam/bet/t_step)*D_n + (gam/bet-1)*Dp_n + t_step*(gam/2/bet-1)*Dpp_n) - ...
%                 alp_m*M_mat*Dpp_n + (1-alp_m)*M_mat*((1/bet/t_step/t_step)*(D_n + t_step*Dp_n) + (1/2/bet-1)*Dpp_n);
%             
%             D_n1 = Jac\Res;
%             Dp_n1 = (gam/bet/t_step)*(D_n1-D_n) - (gam/bet-1)*Dp_n - t_step*(gam/2/bet-1)*Dpp_n;
%             Dpp_n1 = (1/bet/t_step/t_step)*(D_n1 - D_n - t_step*Dp_n) - (1/2/bet-1)*Dpp_n;
%             
%         end
%         
%         answer = zeros(GlobalDOF,1); answer(fdof) = D_n1; answer = reshape(answer,6,N_nodes)';
% 
%         disp(['step: ' sprintf('%4i',i) '  iter: ' sprintf('%4i',qq)  ...
%             '  error ' sprintf('%10.4f',aero_err) ...
%             '  CL: ' sprintf('%6.5f',CL(i)) ...
%             '  tip: ' sprintf('%6.5f',answer(N_nodes,3))])
% 
%         log(qq,:) = [qq,CL(i),aero_err];
% 
%         if qq > 200
%             'Aeroelastic Loop Failed!'
%             return
%         end
% 
%     end
%     
%     answer = zeros(GlobalDOF,1); answer(fdof) = D_n1; answer = reshape(answer,6,N_nodes)';
%     tip(i,:) = [answer(N_nodes,1),answer(N_nodes,2),answer(N_nodes,3)];
% 
%     dKE(i) = 2*Dp_n1'*M_mat*Dpp_n1;
%     dPE(i) = 2*Dp_n1'*K*D_n1;
%     
%     foo = reshape(F_iner_full,6,N_nodes)';
%     rigid_inertial(i) = -sum((foo(:,4) + foo(:,3).*Yo - foo(:,2).*Zo)*omega_body(1,i) + ...
%         (foo(:,5) - foo(:,3).*Xo + foo(:,1).*Zo)*omega_body(2,i) + ...
%         (foo(:,6) + foo(:,2).*Xo - foo(:,1).*Yo)*omega_body(3,i));
%     
%     power = aero_power + dKE + rigid_inertial + dPE;
%     
%     if i > 1
% 
%         %% Wing on wake rollup
% 
%         [D1x,D1y,D1z] = calc_D1_new(X_bar_wake,Y_bar_wake,Z_bar_wake,Xr(:,i),Yr(:,i),Zr(:,i),M,N,N_w,i,cutoff);
%         
%         %% Wake on wake rollup
% 
%         [D2x,D2y,D2z] = calc_D2_new(X_bar_wake,Y_bar_wake,Z_bar_wake,N,N_w,i,cutoff);
%         
%         %% Move wake
% 
%         [Xwake,Ywake,Zwake,E] = wake_move(X_bar_wake,Y_bar_wake,Z_bar_wake,circ(:,i),circ_wake,D1x,D1y,D1z,D2x,D2y,D2z,t_step,N,i,N_w);
% 
%     else
% 
%         Xwake = X_bar_wake; Ywake = Y_bar_wake; Zwake = Z_bar_wake;
% 
%     end
%     
%     out = [CL',CT',power'];
%     
%     %% Plotting:
% 
%     x_temp = reshape(Xr(:,i),N+1,M+1)'; y_temp = reshape(Yr(:,i),N+1,M+1)'; z_temp = reshape(Zr(:,i),N+1,M+1)';
%     count = 1;
%     for m = 1:M
%         for n = 1:N
%             fooX(:,count) = [x_temp(m,n);x_temp(m,n+1);x_temp(m+1,n+1);x_temp(m+1,n)]; fooY(:,count) = [y_temp(m,n);y_temp(m,n+1);y_temp(m+1,n+1);y_temp(m+1,n)]; fooZ(:,count) = [z_temp(m,n);z_temp(m,n+1);z_temp(m+1,n+1);z_temp(m+1,n)];
%             count = count + 1;
%         end
%     end
% 
%     x_temp = reshape(Xwake,N+1,i)'; y_temp = reshape(Ywake,N+1,i)'; z_temp = reshape(Zwake,N+1,i)';
%     if i > 1
%         count = 1;
%         for q = 1:i-1
%             for n = 1:N
%                 fooX_wake(:,count) = [x_temp(q,n);x_temp(q,n+1);x_temp(q+1,n+1);x_temp(q+1,n)]; fooY_wake(:,count) = [y_temp(q,n);y_temp(q,n+1);y_temp(q+1,n+1);y_temp(q+1,n)]; fooZ_wake(:,count) = [z_temp(q,n);z_temp(q,n+1);z_temp(q+1,n+1);z_temp(q+1,n)];
%                 count = count + 1;
%             end
%         end
% 
%     end
% 
%     if i <= 1
%         subplot(2,2,1),patch(fooX-path_X(i),fooY,fooZ,circ(:,i)')
%         title('circulation')
%     else
%         subplot(2,2,1),patch(fooX-path_X(i),fooY,fooZ,circ(:,i)')
%         hold on
%         patch(fooX_wake-path_X(i),fooY_wake,fooZ_wake,circ_wake')
%         title('circulation')
%     end
% 
%     axis equal, axis([0 c*10 0 (b/2)*4/3 -2*c 2*c]), grid on, view([-30 30]), colorbar
% 
%     x_temp = reshape(x(:,i),N+1,M+1)'; y_temp = reshape(y(:,i),N+1,M+1)'; z_temp = reshape(z(:,i),N+1,M+1)';
%     count = 1;
%     for m = 1:M
%         for n = 1:N
%             fooX(:,count) = [x_temp(m,n);x_temp(m,n+1);x_temp(m+1,n+1);x_temp(m+1,n)]; fooY(:,count) = [y_temp(m,n);y_temp(m,n+1);y_temp(m+1,n+1);y_temp(m+1,n)]; fooZ(:,count) = [z_temp(m,n);z_temp(m,n+1);z_temp(m+1,n+1);z_temp(m+1,n)];
%             count = count + 1;
%         end
%     end
% 
%     subplot(2,2,2),surf(reshape(Xo,N+1,M+1)',reshape(Yo,N+1,M+1)',reshape(Zo,N+1,M+1)','facecolor','none')
%     hold on
%     patch(fooX,fooY,fooZ,pressure(:,i)','edgecolor','k')
%     axis equal, axis([0 c 0 b/2 -c c]), axis off, grid on, view([-30 30]), colorbar
%     title('pressure')
% 
%     subplot(4,3,7),plot(t(1:i)/(2*pi/omega),CL(1:i),'k.-',t(1:i)/(2*pi/omega),CT(1:i),'r.-')
%     title('C_L, C_T')
%     subplot(4,3,10),plot(t(1:i)/(2*pi/omega),tip(:,1),'k.-',t(1:i)/(2*pi/omega),tip(:,2),'r.-',t(1:i)/(2*pi/omega),tip(:,3),'b.-')
%     title('tip displacement')    
%     subplot(2,3,5),plot(t(1:i)/(2*pi/omega),aero_power,'r-',t(1:i)/(2*pi/omega),dKE,'b-',t(1:i)/(2*pi/omega),rigid_inertial,'g-',t(1:i)/(2*pi/omega),dPE,'m-',t(1:i)/(2*pi/omega),power,'k-')
%     lgd = legend('aero','dKE','iner','dSE','total');
%     set(lgd,'location','NorthWest','fontsize',8,'box','off')
%     title('power')
%     subplot(2,3,6),drawConnectivity(connectivity);
%     xlabel x
%     ylabel y
%     hold on
%     patch(Xo(trielements'),Yo(trielements'),Zo(trielements')-.1,thickness','edgecolor','none')
%     colorbar
%     title thickness
% 
%     pause(0.01)
% 
%     if i < length(t), clf
%     end
% 
% end
% 
% %% Compute time-integrated objective function
% 
% weight(N_Q_begin,1) = .5*(t(N_Q_begin+1)-t(N_Q_begin))/(t(N_steps+1)-t(N_Q_begin));
% for i = N_Q_begin+1:N_steps
%     weight(i,1) = .5*(t(i+1)-t(i-1))/(t(N_steps+1)-t(N_Q_begin));
% end
% weight(N_steps+1,1) = .5*(t(N_steps+1)-t(N_steps))/(t(N_steps+1)-t(N_Q_begin));
% 
% g_obj = [CL;CT;power]*weight;
% 
% eff = .5*g_obj(2)*rho*V_ref^2*b*c*V_ref/g_obj(3);
% 
% CL_ave = g_obj(1);
% CT_ave = g_obj(2);
% CP_ave = g_obj(3)/(.5*rho*V_ref^2*b*c*V_ref);
% 
% peak_disp = max(abs(tip(:,3)));
