close all;

clear;  d=3; A=1;
saveresults = true;

filename = ['data/lattice/Lattice_', 'n960_ParticlesNumber_2304_Rank_50.mat'];

load(filename)
grids = {xcol,ycol,zcol};

% get the CP tensor
Ut1 = Ldip.Ut1./(2/Lat.x(end));
Ut2 = Ldip.Ut2./(2/Lat.y(end));
Ut3 = Ldip.Ut3./(1/Lat.z(end));
Rp = size(ASC.U1k,2);
Rn = size(AvC.U1k,2);
Rt = Rp + Rn;
LamUt = [ones(Rp,1);-ones(Rn,1)];

% separate the short and long-range parts
num_short = 11;
inds = [1:Rp-num_short,Rp+1:Rt-num_short];
Ut = {Ut1(:,inds),Ut2(:,inds),Ut3(:,inds)};
LamUt = LamUt(inds);
Ut_CP = CP_normalize({LamUt,Ut});

% compute the middle slice of the tensor
ns1=n1/2;
F = CP_get_subtensor(LamUt,Ut,{1:n1,1:n1,ns1:ns1+1});
X_slice = F(:,:,1);

% parameters for ChebTuck
m = [129,129,129];
tol = 1e-7;

% compute ChebTuck
[f,ff] = ChebTuck(Ut_CP,m,[],tol);

% compute the ChebTuck error
[xx, yy, zz] = ndgrid(xcol,xcol,xcol(ns1:ns1+1));
X_cheb1 = f(xx,yy,zz);
X_ChebTuck = X_cheb1(:,:,1);


figure(1)
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
if saveresults
    saveas(gcf,['figures/La_original_Rs' num2str(num_short) '.png'])
end

% compute Chebyshev interpolation error
X_ChebInterp = zeros(n1,n1);
for i = 1:n1
    X_ChebInterp(i,:) = FCP_eval(Ut_CP{1},ff,[xcol(i)*ones(n1,1),xcol',xcol(ns1)*ones(n1,1)],A);
end

figure(2)
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
    saveas(gcf,['figures/La_compress_err_Rs' num2str(num_short) '.png'])
end

figure(3)
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
    saveas(gcf,['figures/La_total_err_Rs' num2str(num_short) '.png'])
end
