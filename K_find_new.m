function[K] = K_find_new(X,Y,Z,trielements,N_elements,asmb,GlobalDOF,fdof,E_modulus,poisson,thickness)

K = zeros(18*18*N_elements,1);
ID = zeros(18*18*N_elements,2);

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

    x12 = xg1(1)-xg2(1); x13 = xg1(1)-xg3(1); x23 = xg2(1)-xg3(1);
    x21 = xg2(1)-xg1(1); x31 = xg3(1)-xg1(1); x32 = xg3(1)-xg2(1);
    y12 = xg1(2)-xg2(2); y13 = xg1(2)-xg3(2); y23 = xg2(2)-xg3(2);
    y21 = xg2(2)-xg1(2); y31 = xg3(2)-xg1(2); y32 = xg3(2)-xg2(2);
    element_area = .5*(y21*x13 - x21*y13);
    l12 = sqrt(x12^2+y12^2); l23 = sqrt(x23^2+y23^2); l31 = sqrt(x31^2+y31^2);

    %% matrices:

    E_constituitive = (E_modulus*thickness(i)/(1-poisson^2))*[1,poisson,0;poisson,1,0;0,0,(1-poisson)/2];
    D_constituitive = (E_modulus*thickness(i)^3/(12*(1-poisson^2)))*[1,poisson,0;poisson,1,0;0,0,(1-poisson)/2];

    P4 = -6*x23/(l23^2);  P5 = -6*x31/(l31^2);  P6 = -6*x12/(l12^2);
    q4 = 3*x23*y23/(l23^2);  q5 = 3*x31*y31/(l31^2);  q6 = 3*x12*y12/(l12^2);
    r4 = 3*y23^2/(l23^2);  r5 = 3*y31^2/(l31^2);  r6 = 3*y12^2/(l12^2);
    t4 = -6*y23/(l23^2);  t5 = -6*y31/(l31^2);  t6 = -6*y12/(l12^2);

    e_gauss = [0.5;0.5;0.0];
    n_gauss = [0.0;0.5;0.5];
    weight = [1/3,1/3,1/3];

    K_DKT = zeros(9,9); K_LST = zeros(9,9);

    for j = 1:length(weight)
        
        e = e_gauss(j);
        n = n_gauss(j);

        %% LST element:

        B_LST = [-y12-y31,0,-1/2*(y12+y31)*(e*y12-n*y31),y31,0,1/2*e*y12^2-1/2*y12*y31+1/2*e*y12*y23+e*y12*y31+1/2*n*y12*y31+1/2*n*y23*y31,y12,0,(y12*y31)/2-(n*y31^2)/2-(e*y12*y23)/2-(e*y12*y31)/2-n*y12*y31-(n*y23*y31)/2;
                 0,-x13-x21,-1/2*(x13+x21)*(e*x21-n*x13),0,x13,1/2*e*x21^2-1/2*x13*x21+e*x13*x21+1/2*e*x21*x32+1/2*n*x13*x21+1/2*n*x13*x32,0,x21,(x13*x21)/2-(n*x13^2)/2-(e*x13*x21)/2-(e*x21*x32)/2-n*x13*x21-(n*x13*x32)/2;
                 -x13-x21,-y12-y31,1/2*n*x13*y12-e*x21*y12-1/2*e*x21*y31-1/2*e*x13*y12+n*x13*y31+1/2*n*x21*y31,x13,y31,1/8*y12*(4*e*x21+x13*(8*e+4*n-4))+1/8*x21*(4*e*y12+y31*(8*e+4*n-4))+1/8*y23*(4*e*x21+4*n*x13)+1/8*x32*(4*e*y12+4*n*y31),x21,y12,-1/8*y31*(4*n*x13+x21*(4*e+8*n-4))-1/8*x13*(4*n*y31+y12*(4*e+8*n-4))-1/8*y23*(4*e*x21+4*n*x13)-1/8*x32*(4*e*y12+4*n*y31)]/2/element_area;
        K_LST = weight(j)*element_area*B_LST'*E_constituitive*B_LST + K_LST;

        %% DKT element:

        H1 = [P6*(1-2*e)+(P5-P6)*n, q6*(1-2*e)-(q5+q6)*n, -4+6*(e+n)+r6*(1-2*e)-n*(r5+r6), -P6*(1-2*e)+n*(P4+P6), q6*(1-2*e)-n*(q6-q4), -2+6*e+r6*(1-2*e)+n*(r4-r6), -n*(P5+P4), n*(q4-q5), -n*(r5-r4)];
        H2 = [t6*(1-2*e)+(t5-t6)*n, 1+r6*(1-2*e)-(r5+r6)*n, -q6*(1-2*e)+n*(q5+q6), -t6*(1-2*e)+n*(t4+t6), -1+r6*(1-2*e)+n*(r4-r6), -q6*(1-2*e)-n*(q4-q6), -n*(t4+t5), n*(r4-r5), -n*(q4-q5)];
        H3 = [-P5*(1-2*n)-e*(P6-P5), q5*(1-2*n)-e*(q5+q6), -4+6*(e+n)+r5*(1-2*n)-e*(r5+r6), e*(P4+P6), e*(q4-q6), -e*(r6-r4), P5*(1-2*n)-e*(P4+P5), q5*(1-2*n)+e*(q4-q5), -2+6*n+r5*(1-2*n)+e*(r4-r5)];
        H4 = [-t5*(1-2*n)-e*(t6-t5), 1+r5*(1-2*n)-e*(r5+r6), -q5*(1-2*n)+e*(q5+q6), e*(t4+t6), e*(r4-r6), -e*(q4-q6), t5*(1-2*n)-e*(t4+t5), -1+r5*(1-2*n)+e*(r4-r5), -q5*(1-2*n)-e*(q4-q5)];
        Bb = (1/(2*element_area))*[y31*H1+y12*H3;-x31*H2-x12*H4;-x31*H1-x12*H3+y31*H2+y12*H4];
        K_DKT = weight(j)*element_area*Bb'*D_constituitive*Bb + K_DKT;
        
    end

    %% Combine LST and DKT:

    K_l = zeros(18,18);
    K_l([1,2,6,7,8,12,13,14,18],[1,2,6,7,8,12,13,14,18]) = K_LST;
    K_l([3,4,5,9,10,11,15,16,17],[3,4,5,9,10,11,15,16,17]) = K_DKT;

    %% Final terms:
    
    Kel = E*K_l*E';

    %% Assemble:

    temp = asmb(i,:)'; temp = temp(:,[ones(18,1)]); temp2 = temp';
    ID((i-1)*18*18+1:i*18*18,:) = [temp(:),temp2(:)];
    
    K((i-1)*18*18+1:i*18*18,:) = Kel(:);
    
end

K = sparse(ID(:,1),ID(:,2),K,GlobalDOF,GlobalDOF);

K = K(fdof,fdof);
