function[C,M] = CM_find_new(X,Y,Z,trielements,N_elements,asmb,GlobalDOF,fdof,rho,thickness,damping)

alpha_b = 3/2;

e_gauss = [1/3;.47014;0.05971;.47014;.10128;.79742;.10128];
n_gauss = [1/3;.47014;.47014;.05971;.10128;.10128;.79742];
weight = [.225;.13239;.13239;.13239;.12593;.12593;.12593];

M = zeros(18*18*N_elements,3);

for i = 1:N_elements

    Xg1 = [X(trielements(i,1)); Y(trielements(i,1)); Z(trielements(i,1))];
    Xg2 = [X(trielements(i,2)); Y(trielements(i,2)); Z(trielements(i,2))];
    Xg3 = [X(trielements(i,3)); Y(trielements(i,3)); Z(trielements(i,3))];

    v1 = (Xg2-Xg1)/sqrt((Xg2(1)-Xg1(1))^2 + (Xg2(2)-Xg1(2))^2 + (Xg2(3)-Xg1(3))^2);
    v13 = (Xg3-Xg1)/sqrt((Xg3(1)-Xg1(1))^2 + (Xg3(2)-Xg1(2))^2 + (Xg3(3)-Xg1(3))^2);
    v3 = [v1(2)*v13(3)-v1(3)*v13(2);-v1(1)*v13(3)+v1(3)*v13(1);v1(1)*v13(2)-v1(2)*v13(1)];
    v3 = v3/sqrt(v3(1)^2 + v3(2)^2 + v3(3)^2);
    v2 = [v3(2)*v1(3)-v3(3)*v1(2);-v3(1)*v1(3)+v3(3)*v1(1);v3(1)*v1(2)-v3(2)*v1(1)];
    v2 = v2/sqrt(v2(1)^2 + v2(2)^2 + v2(3)^2);
    Ro = [v1,v2,v3];
    E = zeros(18,18); E(1:3,1:3) = Ro; E(4:6,4:6) = Ro; E(7:9,7:9) = Ro; E(10:12,10:12) = Ro; E(13:15,13:15) = Ro; E(16:18,16:18) = Ro;
    
    ro1 = Ro'*(Xg1-Xg1); ro2 = Ro'*(Xg2-Xg1); ro3 = Ro'*(Xg3-Xg1);

    x1 = ro1(1); y1 = ro1(2);
    x2 = ro2(1); y2 = ro2(2);
    x3 = ro3(1); y3 = ro3(2);

    x12 = x1-x2; x13 = x1-x3; x23 = x2-x3;
    x21 = x2-x1; x31 = x3-x1; x32 = x3-x2;
    y12 = y1-y2; y13 = y1-y3; y23 = y2-y3;
    y21 = y2-y1; y31 = y3-y1; y32 = y3-y2;
    element_area = .5*(y21*x13 - x21*y13);

    Me = zeros(18,18);

    for j = 1:7

        e = e_gauss(j);
        n = n_gauss(j);

        N = [1-e-n,0,0,0,0,(alpha_b*(1-e-n)/2)*(y12*e-y31*n),e,0,0,0,0,(alpha_b*e/2)*(y23*n-y12*(1-e-n)),n,0,0,0,0,(alpha_b*n/2)*(y31*(1-e-n)-y23*e);
            0,1-e-n,0,0,0,(alpha_b*(1-e-n)/2)*(x21*e-x13*n),0,e,0,0,0,(alpha_b*e/2)*(x32*n-x21*(1-e-n)),0,n,0,0,0,(alpha_b*n/2)*(x13*(1-e-n)-x32*e);
            0,0,(1-e-n)^2*(3-2*(1-e-n))+2*(1-e-n)*e*n,-(1-e-n)^2*(y12*e+y13*n)-.5*(y12+y13)*(1-e-n)*e*n,(1-e-n)^2*(x12*e+x13*n)+.5*(x12+x13)*(1-e-n)*e*n,0,0,0,e^2*(3-2*e)+2*(1-e-n)*e*n,-(e)^2*(y23*n+y21*(1-e-n))-.5*(y23+y21)*(1-e-n)*e*n,(e)^2*(x23*n+x21*(1-e-n))+.5*(x23+x21)*(1-e-n)*e*n,0,0,0,n^2*(3-2*n)+2*(1-e-n)*e*n,-(n)^2*(y31*(1-e-n)+y32*e)-.5*(y31+y32)*(1-e-n)*e*n,(n)^2*(x31*(1-e-n)+x32*e)+.5*(x31+x32)*(1-e-n)*e*n,0];

        Me = Me + rho*thickness(i)*element_area*weight(j)*N'*N;

    end

    Me = E*Me*E';
    
    temp = asmb(i,:)'; temp = temp(:,[ones(18,1)]); temp2 = temp';
    M((i-1)*18*18+1:i*18*18,:) = [temp(:),temp2(:),Me(:)];
      
end

M = sparse(M(:,1),M(:,2),M(:,3),GlobalDOF,GlobalDOF);
M = M(fdof,fdof);

C = damping*M;
