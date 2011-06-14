function[Faero] = Faero_find_new(X,Y,Z,trielements,N_elements,asmb,GlobalDOF,fdof,pressure)

Faero = zeros(GlobalDOF,1);

for i = 1:N_elements

    %% Undeformed element nodes in global system:

    Xg1 = [X(trielements(i,1)); Y(trielements(i,1)); Z(trielements(i,1))];
    Xg2 = [X(trielements(i,2)); Y(trielements(i,2)); Z(trielements(i,2))];
    Xg3 = [X(trielements(i,3)); Y(trielements(i,3)); Z(trielements(i,3))];

    %% Transformation matrix for updated configuration:
    
    v1 = (Xg2-Xg1)/sqrt((Xg2(1)-Xg1(1))^2 + (Xg2(2)-Xg1(2))^2 + (Xg2(3)-Xg1(3))^2);
    v13 = (Xg3-Xg1)/sqrt((Xg3(1)-Xg1(1))^2 + (Xg3(2)-Xg1(2))^2 + (Xg3(3)-Xg1(3))^2);
    v3 = [v1(2)*v13(3)-v1(3)*v13(2);-v1(1)*v13(3)+v1(3)*v13(1);v1(1)*v13(2)-v1(2)*v13(1)];
    v3 = v3/sqrt(v3(1)^2 + v3(2)^2 + v3(3)^2);
    v2 = [v3(2)*v1(3)-v3(3)*v1(2);-v3(1)*v1(3)+v3(3)*v1(1);v3(1)*v1(2)-v3(2)*v1(1)];
    v2 = v2/sqrt(v2(1)^2 + v2(2)^2 + v2(3)^2);
    Ro = [v1,v2,v3];
    E = zeros(18,18); E(1:3,1:3) = Ro; E(4:6,4:6) = Ro; E(7:9,7:9) = Ro; E(10:12,10:12) = Ro; E(13:15,13:15) = Ro; E(16:18,16:18) = Ro;
    
    %% Local system information:

    xg1 = Ro'*(Xg1-Xg1);
    xg2 = Ro'*(Xg2-Xg1);
    xg3 = Ro'*(Xg3-Xg1);

    x13 = xg1(1)-xg3(1); x21 = xg2(1)-xg1(1);
    y13 = xg1(2)-xg3(2); y21 = xg2(2)-xg1(2);
    element_area = .5*(y21*x13 - x21*y13);
   
    %% Pressure forces:
    
    Fel = E*(pressure(i)*element_area/3)*[0;0;1;0;0;0;0;0;1;0;0;0;0;0;1;0;0;0];
    
    Faero(asmb(i,:),1) = Faero(asmb(i,:),1) + Fel;
    
end

Faero = Faero(fdof,1);

