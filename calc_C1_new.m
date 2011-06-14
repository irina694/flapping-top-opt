function[C1,C1_drag] = calc_C1_new(Xr,Yr,Zr,Xc,Yc,Zc,M,N,vec)

Xr = reshape(Xr,N+1,M+1)';
Yr = reshape(Yr,N+1,M+1)';
Zr = reshape(Zr,N+1,M+1)';

Xc = reshape(Xc,N,M)';
Yc = reshape(Yc,N,M)';
Zc = reshape(Zc,N,M)';

logX(:,1) = reshape((reshape(Xc',M*N,1)*ones(1,M*N))',M*M*N*N,1);
logY(:,1) = reshape((reshape(Yc',M*N,1)*ones(1,M*N))',M*M*N*N,1);
logZ(:,1) = reshape((reshape(Zc',M*N,1)*ones(1,M*N))',M*M*N*N,1);
logX(:,2) = reshape((reshape(Xr(1:M,1:N)',M*N,1)*ones(1,M*N)),M*M*N*N,1);
logY(:,2) = reshape((reshape(Yr(1:M,1:N)',M*N,1)*ones(1,M*N)),M*M*N*N,1);
logZ(:,2) = reshape((reshape(Zr(1:M,1:N)',M*N,1)*ones(1,M*N)),M*M*N*N,1);
logX(:,3) = reshape((reshape(Xr(1:M,2:N+1)',M*N,1)*ones(1,M*N)),M*M*N*N,1);
logY(:,3) = reshape((reshape(Yr(1:M,2:N+1)',M*N,1)*ones(1,M*N)),M*M*N*N,1);
logZ(:,3) = reshape((reshape(Zr(1:M,2:N+1)',M*N,1)*ones(1,M*N)),M*M*N*N,1);
logX(:,4) = reshape((reshape(Xr(2:M+1,2:N+1)',M*N,1)*ones(1,M*N)),M*M*N*N,1);
logY(:,4) = reshape((reshape(Yr(2:M+1,2:N+1)',M*N,1)*ones(1,M*N)),M*M*N*N,1);
logZ(:,4) = reshape((reshape(Zr(2:M+1,2:N+1)',M*N,1)*ones(1,M*N)),M*M*N*N,1);
logX(:,5) = reshape((reshape(Xr(2:M+1,1:N)',M*N,1)*ones(1,M*N)),M*M*N*N,1);
logY(:,5) = reshape((reshape(Yr(2:M+1,1:N)',M*N,1)*ones(1,M*N)),M*M*N*N,1);
logZ(:,5) = reshape((reshape(Zr(2:M+1,1:N)',M*N,1)*ones(1,M*N)),M*M*N*N,1);

TE_ID = reshape(1:M*N*M*N,M*N,M*N); TE_ID = reshape(TE_ID(M*N-N+1:M*N,:),M*N*N,1);

%1
ro = [logX(:,3)-logX(:,2),logY(:,3)-logY(:,2),logZ(:,3)-logZ(:,2)];
r1 = [logX(:,1)-logX(:,2),logY(:,1)-logY(:,2),logZ(:,1)-logZ(:,2)];
r2 = [logX(:,1)-logX(:,3),logY(:,1)-logY(:,3),logZ(:,1)-logZ(:,3)];
aa = [r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2),r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3),r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)];
bb = (r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2)).^2+(r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3)).^2+(r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)).^2;
cc = [r1./repmat(sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2),1,3)-r2./repmat(sqrt(r2(:,1).^2+r2(:,2).^2+r2(:,3).^2),1,3)];
UVW = (aa./bb(:,[1,1,1])).*repmat((ro(:,1).*cc(:,1)+ro(:,2).*cc(:,2)+ro(:,3).*cc(:,3)),1,3)/4/pi;

%2
ro = [logX(:,4)-logX(:,3),logY(:,4)-logY(:,3),logZ(:,4)-logZ(:,3)];
r1 = [logX(:,1)-logX(:,3),logY(:,1)-logY(:,3),logZ(:,1)-logZ(:,3)];
r2 = [logX(:,1)-logX(:,4),logY(:,1)-logY(:,4),logZ(:,1)-logZ(:,4)];
aa = [r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2),r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3),r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)];
bb = (r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2)).^2+(r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3)).^2+(r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)).^2;
cc = [r1./repmat(sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2),1,3)-r2./repmat(sqrt(r2(:,1).^2+r2(:,2).^2+r2(:,3).^2),1,3)];
temp = (aa./bb(:,[1,1,1])).*repmat((ro(:,1).*cc(:,1)+ro(:,2).*cc(:,2)+ro(:,3).*cc(:,3)),1,3)/4/pi;
UVW = UVW + temp;
UVW_stream = temp;

%3
ro = [logX(:,5)-logX(:,4),logY(:,5)-logY(:,4),logZ(:,5)-logZ(:,4)];
r1 = [logX(:,1)-logX(:,4),logY(:,1)-logY(:,4),logZ(:,1)-logZ(:,4)];
r2 = [logX(:,1)-logX(:,5),logY(:,1)-logY(:,5),logZ(:,1)-logZ(:,5)];
aa = [r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2),r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3),r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)];
bb = (r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2)).^2+(r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3)).^2+(r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)).^2;
cc = [r1./repmat(sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2),1,3)-r2./repmat(sqrt(r2(:,1).^2+r2(:,2).^2+r2(:,3).^2),1,3)];
temp = (aa./bb(:,[1,1,1])).*repmat((ro(:,1).*cc(:,1)+ro(:,2).*cc(:,2)+ro(:,3).*cc(:,3)),1,3)/4/pi;
UVW = UVW + temp;
UVW_stream(TE_ID,:) = UVW_stream(TE_ID,:) + temp(TE_ID,:); 

%4
ro = [logX(:,2)-logX(:,5),logY(:,2)-logY(:,5),logZ(:,2)-logZ(:,5)];
r1 = [logX(:,1)-logX(:,5),logY(:,1)-logY(:,5),logZ(:,1)-logZ(:,5)];
r2 = [logX(:,1)-logX(:,2),logY(:,1)-logY(:,2),logZ(:,1)-logZ(:,2)];
aa = [r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2),r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3),r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)];
bb = (r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2)).^2+(r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3)).^2+(r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)).^2;
cc = [r1./repmat(sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2),1,3)-r2./repmat(sqrt(r2(:,1).^2+r2(:,2).^2+r2(:,3).^2),1,3)];
temp = (aa./bb(:,[1,1,1])).*repmat((ro(:,1).*cc(:,1)+ro(:,2).*cc(:,2)+ro(:,3).*cc(:,3)),1,3)/4/pi;
UVW = UVW + temp;
UVW_stream = UVW_stream + temp;

%1r
ro = [logX(:,3)-logX(:,2),logY(:,3)-logY(:,2),logZ(:,3)-logZ(:,2)];
r1 = [logX(:,1)-logX(:,2),-logY(:,1)-logY(:,2),logZ(:,1)-logZ(:,2)];
r2 = [logX(:,1)-logX(:,3),-logY(:,1)-logY(:,3),logZ(:,1)-logZ(:,3)];
aa = [r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2),r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3),r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)];
bb = (r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2)).^2+(r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3)).^2+(r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)).^2;
cc = [r1./repmat(sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2),1,3)-r2./repmat(sqrt(r2(:,1).^2+r2(:,2).^2+r2(:,3).^2),1,3)];
temp = (aa./bb(:,[1,1,1])).*repmat((ro(:,1).*cc(:,1)+ro(:,2).*cc(:,2)+ro(:,3).*cc(:,3)),1,3)/4/pi;
temp(:,2) = -temp(:,2);
UVW = UVW + temp;

%2r
ro = [logX(:,4)-logX(:,3),logY(:,4)-logY(:,3),logZ(:,4)-logZ(:,3)];
r1 = [logX(:,1)-logX(:,3),-logY(:,1)-logY(:,3),logZ(:,1)-logZ(:,3)];
r2 = [logX(:,1)-logX(:,4),-logY(:,1)-logY(:,4),logZ(:,1)-logZ(:,4)];
aa = [r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2),r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3),r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)];
bb = (r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2)).^2+(r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3)).^2+(r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)).^2;
cc = [r1./repmat(sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2),1,3)-r2./repmat(sqrt(r2(:,1).^2+r2(:,2).^2+r2(:,3).^2),1,3)];
temp = (aa./bb(:,[1,1,1])).*repmat((ro(:,1).*cc(:,1)+ro(:,2).*cc(:,2)+ro(:,3).*cc(:,3)),1,3)/4/pi;
temp(:,2) = -temp(:,2);
UVW = UVW + temp;
UVW_stream = UVW_stream + temp;

%3r
ro = [logX(:,5)-logX(:,4),logY(:,5)-logY(:,4),logZ(:,5)-logZ(:,4)];
r1 = [logX(:,1)-logX(:,4),-logY(:,1)-logY(:,4),logZ(:,1)-logZ(:,4)];
r2 = [logX(:,1)-logX(:,5),-logY(:,1)-logY(:,5),logZ(:,1)-logZ(:,5)];
aa = [r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2),r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3),r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)];
bb = (r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2)).^2+(r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3)).^2+(r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)).^2;
cc = [r1./repmat(sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2),1,3)-r2./repmat(sqrt(r2(:,1).^2+r2(:,2).^2+r2(:,3).^2),1,3)];
temp = (aa./bb(:,[1,1,1])).*repmat((ro(:,1).*cc(:,1)+ro(:,2).*cc(:,2)+ro(:,3).*cc(:,3)),1,3)/4/pi;
temp(:,2) = -temp(:,2);
UVW = UVW + temp;
UVW_stream(TE_ID,:) = UVW_stream(TE_ID,:) + temp(TE_ID,:);

%4r
ro = [logX(:,2)-logX(:,5),logY(:,2)-logY(:,5),logZ(:,2)-logZ(:,5)];
r1 = [logX(:,1)-logX(:,5),-logY(:,1)-logY(:,5),logZ(:,1)-logZ(:,5)];
r2 = [logX(:,1)-logX(:,2),-logY(:,1)-logY(:,2),logZ(:,1)-logZ(:,2)];
aa = [r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2),r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3),r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)];
bb = (r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2)).^2+(r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3)).^2+(r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)).^2;
cc = [r1./repmat(sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2),1,3)-r2./repmat(sqrt(r2(:,1).^2+r2(:,2).^2+r2(:,3).^2),1,3)];
temp = (aa./bb(:,[1,1,1])).*repmat((ro(:,1).*cc(:,1)+ro(:,2).*cc(:,2)+ro(:,3).*cc(:,3)),1,3)/4/pi;
temp(:,2) = -temp(:,2);
UVW = UVW + temp;
UVW_stream = UVW_stream + temp;

outward(:,1) = reshape((vec.outward_x*ones(1,M*N))',M*M*N*N,1);
outward(:,2) = reshape((vec.outward_y*ones(1,M*N))',M*M*N*N,1);
outward(:,3) = reshape((vec.outward_z*ones(1,M*N))',M*M*N*N,1);

C1 = reshape(UVW(:,1).*outward(:,1)+UVW(:,2).*outward(:,2)+UVW(:,3).*outward(:,3),M*N,M*N)';

n_lift(:,1) = reshape((vec.n_lift_x*ones(1,M*N))',M*M*N*N,1);
n_lift(:,2) = reshape((vec.n_lift_y*ones(1,M*N))',M*M*N*N,1);
n_lift(:,3) = reshape((vec.n_lift_z*ones(1,M*N))',M*M*N*N,1);

C1_drag = reshape(UVW_stream(:,1).*n_lift(:,1)+UVW_stream(:,2).*n_lift(:,2)+UVW_stream(:,3).*n_lift(:,3),M*N,M*N)';

