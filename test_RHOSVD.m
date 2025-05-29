close all;

clear;  d=3; A=1;
saveresults = false;

n1 = 256;Pn = 500;
filename = ['data\Data_256-2048\n' num2str(n1) '_Pn' num2str(Pn) '.mat'];

load(filename)

xi = LRed.LAM1C./Hunif^3;
Rl = length(xi);
U = {LRed.CU1, LRed.CU2, LRed.CU3};
grids = {xcol,ycol,zcol};

Full1 = CP_get_subtensor(xi,U,{1:n1,1:n1,1:n1});

[core, V] = hosvd_s(Full1, 1e-7);
Full2 = tucker2full(core,V);

norm(Full1-Full2,'fro')
r = size(core);
[core_r1,V_r1] = rhosvd(xi,U,1e-7,r);
Full31 = tucker2full(core_r1,V_r1);
norm(Full1 - Full31,'fro')

[U1U, U1S, U1V] = svds(U{1}, r(1));
[U2U, U2S, U2V] = svds(U{2}, r(2));
[U3U, U3S, U3V] = svds(U{3}, r(3));

% create a diagonal tensor with elements of xi of size Rl x Rl x Rl
xi_diag = zeros(Rl,Rl,Rl);
xi_diag(1:Rl^2+Rl+1:end) = xi;
core_r = tucker2full(xi_diag,{U1S*U1V',U2S*U2V',U3S*U3V'});
V_r = {U1U, U2U, U3U};

Full3 = tucker2full(core_r,V_r);
norm(Full1-Full3,'fro')
