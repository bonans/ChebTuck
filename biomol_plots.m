close all;

clear;  d=3; A=1;
saveresults = true;

n1 = 256;Pn = 500;
filename = ['data/Data_256-2048/n' num2str(n1) '_Pn' num2str(Pn) '.mat'];

load(filename,"LRed","Hunif")

% get the CP tensor
xi = LRed.LAM1C./Hunif^3;
Rl = length(xi);
U = {LRed.CU1, LRed.CU2, LRed.CU3};

h = 2/n1;
xcol = -1+h/2:h:1-h/2;
figure(1)
plot(xcol,LRed.CU1(:,1:200)/Hunif,'LineWidth',2)
grid on;
set(gca,'fontsize',20);
xticks([-1 -0.5 0 0.5 1])
axis tight
set(gcf,'position',[100,100,900,500])
if saveresults
saveas(gcf,'figures/Pn_cp_vectors.png')
end

% compute the middle slice of the tensor
ns1 = floor(n1/2);
F = CP_get_subtensor(xi,U,{1:n1,1:n1,ns1:ns1+1});
X_slice = F(:,:,1);

% parameters for ChebTuck
m = [129,129,129];
tol = 1e-7;

% compute ChebTuck
[f,ff] = ChebTuck({xi,U},m,[],tol);

% compute the ChebTuck error
[xx, yy, zz] = ndgrid(xcol,xcol,xcol(ns1:ns1+1));
X_cheb1 = f(xx,yy,zz);
X_ChebTuck = X_cheb1(:,:,1);

figure(2)
 mesh(xcol,xcol,X_slice)
 %mesh(A2C2)
    set(gca,'fontsize',24);
    view(3);
    grid on;
    light;
    lighting phong;
    material dull
    camlight('left');
    shading interp;
axis tight
if saveresults && n1 == 2048
    saveas(gcf,['figures/Pn_original.png'])
end

% compute Chebyshev interpolation error
X_ChebInterp = zeros(n1,n1);
for i = 1:n1
    X_ChebInterp(i,:) = FCP_eval(xi,ff,[xcol(i)*ones(n1,1),xcol',xcol(ns1)*ones(n1,1)],A);
end

figure(3)
 mesh(xcol,xcol,X_ChebInterp -X_ChebTuck)
 %mesh(A2C2)
    set(gca,'fontsize',24);
    view(3);
    grid on;
    light;
    lighting phong;
    material dull
    camlight('left');
    shading interp;
axis tight
if saveresults
    saveas(gcf,['figures/Pn_compress_err_n' num2str(n1) '.png'])
end

figure(4)
 mesh(xcol,xcol,X_slice -X_ChebTuck)
 %mesh(A2C2)
    set(gca,'fontsize',24);
    view(3);
    grid on;
    light;
    lighting phong;
    material dull
    camlight('left');
    shading interp;
axis tight
if saveresults
    saveas(gcf,['figures/Pn_total_err_n' num2str(n1) '.png'])
    %saveas(gcf,['figures/Pn_total_err_n' num2str(n1) '_nn.png'])
    %saveas(gcf,['figures/Pn_total_err_n' num2str(n1) '_lin.png'])
end