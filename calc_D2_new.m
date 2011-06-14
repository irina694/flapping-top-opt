function[D2x,D2y,D2z] = calc_D2_new(X_bar_wake,Y_bar_wake,Z_bar_wake,N,N_w,i,cutoff)

X_bar_wake = reshape(X_bar_wake,N+1,i)';
Y_bar_wake = reshape(Y_bar_wake,N+1,i)';
Z_bar_wake = reshape(Z_bar_wake,N+1,i)';

if i-1 < N_w
    wake_start = 1;
else
    wake_start = i-N_w;
end

logX(:,1) = reshape((reshape(X_bar_wake(wake_start:i-1,:)',(i-wake_start)*(N+1),1)*ones(1,(i-1)*N))',(i-1)*N*(i-wake_start)*(N+1),1);
logY(:,1) = reshape((reshape(Y_bar_wake(wake_start:i-1,:)',(i-wake_start)*(N+1),1)*ones(1,(i-1)*N))',(i-1)*N*(i-wake_start)*(N+1),1);
logZ(:,1) = reshape((reshape(Z_bar_wake(wake_start:i-1,:)',(i-wake_start)*(N+1),1)*ones(1,(i-1)*N))',(i-1)*N*(i-wake_start)*(N+1),1);
logX(:,2) = reshape((reshape(X_bar_wake(2:i,1:N)',(i-1)*N,1)*ones(1,(i-wake_start)*(N+1))),(i-1)*N*(i-wake_start)*(N+1),1);
logY(:,2) = reshape((reshape(Y_bar_wake(2:i,1:N)',(i-1)*N,1)*ones(1,(i-wake_start)*(N+1))),(i-1)*N*(i-wake_start)*(N+1),1);
logZ(:,2) = reshape((reshape(Z_bar_wake(2:i,1:N)',(i-1)*N,1)*ones(1,(i-wake_start)*(N+1))),(i-1)*N*(i-wake_start)*(N+1),1);
logX(:,3) = reshape((reshape(X_bar_wake(2:i,2:N+1)',(i-1)*N,1)*ones(1,(i-wake_start)*(N+1))),(i-1)*N*(i-wake_start)*(N+1),1);
logY(:,3) = reshape((reshape(Y_bar_wake(2:i,2:N+1)',(i-1)*N,1)*ones(1,(i-wake_start)*(N+1))),(i-1)*N*(i-wake_start)*(N+1),1);
logZ(:,3) = reshape((reshape(Z_bar_wake(2:i,2:N+1)',(i-1)*N,1)*ones(1,(i-wake_start)*(N+1))),(i-1)*N*(i-wake_start)*(N+1),1);
logX(:,4) = reshape((reshape(X_bar_wake(1:i-1,2:N+1)',(i-1)*N,1)*ones(1,(i-wake_start)*(N+1))),(i-1)*N*(i-wake_start)*(N+1),1);
logY(:,4) = reshape((reshape(Y_bar_wake(1:i-1,2:N+1)',(i-1)*N,1)*ones(1,(i-wake_start)*(N+1))),(i-1)*N*(i-wake_start)*(N+1),1);
logZ(:,4) = reshape((reshape(Z_bar_wake(1:i-1,2:N+1)',(i-1)*N,1)*ones(1,(i-wake_start)*(N+1))),(i-1)*N*(i-wake_start)*(N+1),1);
logX(:,5) = reshape((reshape(X_bar_wake(1:i-1,1:N)',(i-1)*N,1)*ones(1,(i-wake_start)*(N+1))),(i-1)*N*(i-wake_start)*(N+1),1);
logY(:,5) = reshape((reshape(Y_bar_wake(1:i-1,1:N)',(i-1)*N,1)*ones(1,(i-wake_start)*(N+1))),(i-1)*N*(i-wake_start)*(N+1),1);
logZ(:,5) = reshape((reshape(Z_bar_wake(1:i-1,1:N)',(i-1)*N,1)*ones(1,(i-wake_start)*(N+1))),(i-1)*N*(i-wake_start)*(N+1),1);

%1
ro = [logX(:,3)-logX(:,2),logY(:,3)-logY(:,2),logZ(:,3)-logZ(:,2)];
r1 = [logX(:,1)-logX(:,2),logY(:,1)-logY(:,2),logZ(:,1)-logZ(:,2)];
r2 = [logX(:,1)-logX(:,3),logY(:,1)-logY(:,3),logZ(:,1)-logZ(:,3)];
aa = [r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2),r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3),r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)];
bb = (r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2)).^2+(r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3)).^2+(r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)).^2;
cc = [r1./repmat(sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2),1,3)-r2./repmat(sqrt(r2(:,1).^2+r2(:,2).^2+r2(:,3).^2),1,3)];
UVW = (aa./bb(:,[1,1,1])).*repmat((ro(:,1).*cc(:,1)+ro(:,2).*cc(:,2)+ro(:,3).*cc(:,3)),1,3)/4/pi;
UVW = UVW.*(bb(:,[1,1,1])>cutoff);
UVW(isnan(UVW)) = 0;

%2
ro = [logX(:,4)-logX(:,3),logY(:,4)-logY(:,3),logZ(:,4)-logZ(:,3)];
r1 = [logX(:,1)-logX(:,3),logY(:,1)-logY(:,3),logZ(:,1)-logZ(:,3)];
r2 = [logX(:,1)-logX(:,4),logY(:,1)-logY(:,4),logZ(:,1)-logZ(:,4)];
aa = [r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2),r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3),r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)];
bb = (r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2)).^2+(r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3)).^2+(r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)).^2;
cc = [r1./repmat(sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2),1,3)-r2./repmat(sqrt(r2(:,1).^2+r2(:,2).^2+r2(:,3).^2),1,3)];
temp = (aa./bb(:,[1,1,1])).*repmat((ro(:,1).*cc(:,1)+ro(:,2).*cc(:,2)+ro(:,3).*cc(:,3)),1,3)/4/pi;
temp(isnan(temp)) = 0;
UVW = UVW + temp.*(bb(:,[1,1,1])>cutoff);

%3
ro = [logX(:,5)-logX(:,4),logY(:,5)-logY(:,4),logZ(:,5)-logZ(:,4)];
r1 = [logX(:,1)-logX(:,4),logY(:,1)-logY(:,4),logZ(:,1)-logZ(:,4)];
r2 = [logX(:,1)-logX(:,5),logY(:,1)-logY(:,5),logZ(:,1)-logZ(:,5)];
aa = [r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2),r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3),r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)];
bb = (r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2)).^2+(r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3)).^2+(r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)).^2;
cc = [r1./repmat(sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2),1,3)-r2./repmat(sqrt(r2(:,1).^2+r2(:,2).^2+r2(:,3).^2),1,3)];
temp = (aa./bb(:,[1,1,1])).*repmat((ro(:,1).*cc(:,1)+ro(:,2).*cc(:,2)+ro(:,3).*cc(:,3)),1,3)/4/pi;
temp(isnan(temp)) = 0;
UVW = UVW + temp.*(bb(:,[1,1,1])>cutoff);

%4
ro = [logX(:,2)-logX(:,5),logY(:,2)-logY(:,5),logZ(:,2)-logZ(:,5)];
r1 = [logX(:,1)-logX(:,5),logY(:,1)-logY(:,5),logZ(:,1)-logZ(:,5)];
r2 = [logX(:,1)-logX(:,2),logY(:,1)-logY(:,2),logZ(:,1)-logZ(:,2)];
aa = [r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2),r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3),r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)];
bb = (r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2)).^2+(r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3)).^2+(r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)).^2;
cc = [r1./repmat(sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2),1,3)-r2./repmat(sqrt(r2(:,1).^2+r2(:,2).^2+r2(:,3).^2),1,3)];
temp = (aa./bb(:,[1,1,1])).*repmat((ro(:,1).*cc(:,1)+ro(:,2).*cc(:,2)+ro(:,3).*cc(:,3)),1,3)/4/pi;
temp(isnan(temp)) = 0;
UVW = UVW + temp.*(bb(:,[1,1,1])>cutoff);

%1r
ro = [logX(:,3)-logX(:,2),logY(:,3)-logY(:,2),logZ(:,3)-logZ(:,2)];
r1 = [logX(:,1)-logX(:,2),-logY(:,1)-logY(:,2),logZ(:,1)-logZ(:,2)];
r2 = [logX(:,1)-logX(:,3),-logY(:,1)-logY(:,3),logZ(:,1)-logZ(:,3)];
aa = [r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2),r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3),r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)];
bb = (r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2)).^2+(r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3)).^2+(r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)).^2;
cc = [r1./repmat(sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2),1,3)-r2./repmat(sqrt(r2(:,1).^2+r2(:,2).^2+r2(:,3).^2),1,3)];
temp = (aa./bb(:,[1,1,1])).*repmat((ro(:,1).*cc(:,1)+ro(:,2).*cc(:,2)+ro(:,3).*cc(:,3)),1,3)/4/pi;
temp(:,2) = -temp(:,2);
temp(isnan(temp)) = 0;
UVW = UVW + temp.*(bb(:,[1,1,1])>cutoff);

%2r
ro = [logX(:,4)-logX(:,3),logY(:,4)-logY(:,3),logZ(:,4)-logZ(:,3)];
r1 = [logX(:,1)-logX(:,3),-logY(:,1)-logY(:,3),logZ(:,1)-logZ(:,3)];
r2 = [logX(:,1)-logX(:,4),-logY(:,1)-logY(:,4),logZ(:,1)-logZ(:,4)];
aa = [r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2),r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3),r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)];
bb = (r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2)).^2+(r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3)).^2+(r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)).^2;
cc = [r1./repmat(sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2),1,3)-r2./repmat(sqrt(r2(:,1).^2+r2(:,2).^2+r2(:,3).^2),1,3)];
temp = (aa./bb(:,[1,1,1])).*repmat((ro(:,1).*cc(:,1)+ro(:,2).*cc(:,2)+ro(:,3).*cc(:,3)),1,3)/4/pi;
temp(:,2) = -temp(:,2);
temp(isnan(temp)) = 0;
UVW = UVW + temp.*(bb(:,[1,1,1])>cutoff);

%3r
ro = [logX(:,5)-logX(:,4),logY(:,5)-logY(:,4),logZ(:,5)-logZ(:,4)];
r1 = [logX(:,1)-logX(:,4),-logY(:,1)-logY(:,4),logZ(:,1)-logZ(:,4)];
r2 = [logX(:,1)-logX(:,5),-logY(:,1)-logY(:,5),logZ(:,1)-logZ(:,5)];
aa = [r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2),r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3),r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)];
bb = (r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2)).^2+(r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3)).^2+(r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)).^2;
cc = [r1./repmat(sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2),1,3)-r2./repmat(sqrt(r2(:,1).^2+r2(:,2).^2+r2(:,3).^2),1,3)];
temp = (aa./bb(:,[1,1,1])).*repmat((ro(:,1).*cc(:,1)+ro(:,2).*cc(:,2)+ro(:,3).*cc(:,3)),1,3)/4/pi;
temp(:,2) = -temp(:,2);
temp(isnan(temp)) = 0;
UVW = UVW + temp.*(bb(:,[1,1,1])>cutoff);

%4r
ro = [logX(:,2)-logX(:,5),logY(:,2)-logY(:,5),logZ(:,2)-logZ(:,5)];
r1 = [logX(:,1)-logX(:,5),-logY(:,1)-logY(:,5),logZ(:,1)-logZ(:,5)];
r2 = [logX(:,1)-logX(:,2),-logY(:,1)-logY(:,2),logZ(:,1)-logZ(:,2)];
aa = [r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2),r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3),r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)];
bb = (r1(:,2).*r2(:,3)-r1(:,3).*r2(:,2)).^2+(r1(:,3).*r2(:,1)-r1(:,1).*r2(:,3)).^2+(r1(:,1).*r2(:,2)-r1(:,2).*r2(:,1)).^2;
cc = [r1./repmat(sqrt(r1(:,1).^2+r1(:,2).^2+r1(:,3).^2),1,3)-r2./repmat(sqrt(r2(:,1).^2+r2(:,2).^2+r2(:,3).^2),1,3)];
temp = (aa./bb(:,[1,1,1])).*repmat((ro(:,1).*cc(:,1)+ro(:,2).*cc(:,2)+ro(:,3).*cc(:,3)),1,3)/4/pi;
temp(:,2) = -temp(:,2);
temp(isnan(temp)) = 0;
UVW = UVW + temp.*(bb(:,[1,1,1])>cutoff);

D2x = reshape(UVW(:,1),(i-1)*N,(i-wake_start)*(N+1))';
D2y = reshape(UVW(:,2),(i-1)*N,(i-wake_start)*(N+1))';
D2z = reshape(UVW(:,3),(i-1)*N,(i-wake_start)*(N+1))';

if wake_start > 1
    D2x = [sparse((wake_start-1)*(N+1),(i-1)*N);D2x];
    D2y = [sparse((wake_start-1)*(N+1),(i-1)*N);D2y];
    D2z = [sparse((wake_start-1)*(N+1),(i-1)*N);D2z];
end