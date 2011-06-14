function[Xo,Yo,Zo,N_nodes,trielements,N_elements,GlobalDOF,DOF,asmb,bc,fdof] = wing_grid(c,b,M,N)

[Xo,Yo] = meshgrid(0:c/M:c, 0:b/2/N:b/2);
Xo = reshape(Xo,(M+1)*(N+1),1);
Yo = reshape(Yo,(M+1)*(N+1),1);
Zo = Xo*0;

count = 1;
for i = 1:M
    for j = 1:N
        trielements(count,:) = [i, j, i+1, j+1, i, j+1];
        trielements(count + 1,:) = [i, j, i+1, j, i+1, j+1];
        count = count + 2;
    end
end

for i = 1:length(trielements)
    trielements(i,1) = trielements(i,2) + (N+1)*trielements(i,1)-N-1;
    trielements(i,2) = trielements(i,4) + (N+1)*trielements(i,3)-N-1;
    trielements(i,3) = trielements(i,6) + (N+1)*trielements(i,5)-N-1;
end
trielements(:,4:6) = [];

N_nodes = length(Xo);
N_elements = length(trielements);
GlobalDOF = N_nodes*6;

DOF = reshape(1:GlobalDOF,6,N_nodes)';
for i = 1:N_elements
    asmb(i,:) = [DOF(trielements(i,1),:),DOF(trielements(i,2),:),DOF(trielements(i,3),:)];
end

count = 1;
for i = 1:N_nodes
    if Yo(i) == 0
        bc(count:count+5,:) = DOF(i,1:6)';
        count = count + 6;
    end
end

fdof = 1:GlobalDOF;
fdof(bc(:,1)) = [];