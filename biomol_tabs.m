close all;

clear;  d=3; A=1;
saveresults = false;

% parameters for ChebTuck
m = [129,129,129];
tol = 1e-7;

nlist = [256,512,1024,2048];
nrepeat = 1;
errors = zeros(4,3);
runtimes = zeros(4,3);

% load data
for ii = 1:4
n1 = nlist(ii);Pn = 500;
filename = ['data/Data_256-2048/n' num2str(n1) '_Pn' num2str(Pn) '.mat'];
load(filename,"LRed","Hunif")
h = 2/n1;
xcol = -1+h/2:h:1-h/2;

% get the CP tensor
xi = LRed.LAM1C./Hunif^3;
Rl = length(xi);
U = {LRed.CU1, LRed.CU2, LRed.CU3};
% compute the middle slice of the tensor
ns1 = floor(n1/2);
F = CP_get_subtensor(xi,U,{1:n1,1:n1,ns1:ns1+1});
X_slice = F(:,:,1);
% uniform grid
[xx, yy, zz] = ndgrid(xcol,xcol,xcol(ns1:ns1+1));

% % Alg 2
% if n1 == 256
% F = CP_get_subtensor(xi,U,{1:n1,1:n1,1:n1});
% tic
% f = ChebTuck(F,m,[],tol,[]);
% runtimes(ii,1) = toc;
% % compute the ChebTuck error
% X_cheb1 = f(xx,yy,zz);
% X_ChebTuck = X_cheb1(:,:,1)';
% errors(ii,1) = norm(X_slice(:) - X_ChebTuck(:),'inf');
% end

% Alg 3
tic
f = ChebTuck({xi,U},m,[],tol,[],'F');
runtimes(ii,2) = toc;
% compute the ChebTuck error
X_cheb1 = f(xx,yy,zz);
X_ChebTuck = X_cheb1(:,:,1);
errors(ii,2) = norm(X_slice(:) - X_ChebTuck(:),'inf');

% Alg 4
tic
f = ChebTuck({xi,U},m,[],tol);
runtimes(ii,3) = toc;
% compute the ChebTuck error
X_cheb1 = f(xx,yy,zz);
X_ChebTuck = X_cheb1(:,:,1);
errors(ii,3) = norm(X_slice(:) - X_ChebTuck(:),'inf');
end

