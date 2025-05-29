close all;

clear;  d=3; A=1;
saveresults = true;
n1 = 256;
filename = ['data/Newton_CP_-1_1/Newt_canon_' num2str(n1) '.mat'];
load(filename)

% do the plot of the canonical vectors of the full Newton kernel
h1 = 2/n1;
xcol = -A +0.5*h1:h1:A -0.5*h1;
figure(1)
plot(xcol,U1r./h1,'LineWidth',1)
grid on;
set(gca,'fontsize',20);
xticks([-1 -0.5 0 0.5 1])
axis tight
set(gcf,'position',[100,100,900,500])
if saveresults
saveas(gcf,'figures/cp_vectors.png')
end

% do the 1d Cheb error plots
n1 = 8192;
filename = ['data/Newton_CP_-1_1/Newt_canon_' num2str(n1) '.mat'];
load(filename)

h1 = 2/n1;
xcol = -A +0.5*h1:h1:A -0.5*h1;
xcol = xcol';
tol = 1e-7;
U1r = U1r./h1;
R = size(U1r,2);
xi = ones(R,1);
U = {U1r, U1r, U1r};

cheb_degs = 9:10:200;%2.^(4:16) + 1;
N_degs = length(cheb_degs);
errors = zeros(N_degs, R);
for jj = 1:N_degs
    for k = 1:R
        f = chebfun(@(x) interp1(xcol,U1r(:,k),x,'spline','extrap'),cheb_degs(jj));
        Freal = f(xcol);
        errors(jj,k) = norm(Freal - U1r(:,k),'inf')/norm(Freal,'inf');
    end
end

figure(2)
semilogy(cheb_degs,errors(:,1:15),'LineWidth',2)
grid on;
set(gca,'fontsize',22);
xlabel('polynomial degree')
ylabel('relative error')
ylim([1e-13,1e1])
yticks([10^(-20),10^(-16),10^(-12),10^(-8),10^(-4),1])
xticks([40,80,120,160,200])
if saveresults
saveas(gcf,'figures/newton_err_vs_deg_l.png')
end

figure(3)
semilogy(cheb_degs,errors(:,16:end),'LineWidth',2)
grid on;
set(gca,'fontsize',22);
xlabel('polynomial degree')
ylim([1e-13,1e1])
yticks([10^(-20),10^(-16),10^(-12),10^(-8),10^(-4),1])
xticks([40,80,120,160,200])
if saveresults
saveas(gcf,'figures/newton_err_vs_deg_s.png')
end

figure(4)
semilogy(errors(13,:),'LineWidth',2)
grid on;
set(gca,'fontsize',22);
xlabel('k-th canonical vector')
ylim([1e-13,1e1])
yticks([10^(-20),10^(-16),10^(-12),10^(-8),10^(-4),1])
xticks([10,20,30,40])
if saveresults
saveas(gcf,'figures/newton_err_vs_k.png')
end

figure(5)
f1 = chebfun(@(x) interp1(xcol,U1r(:,10),x,'spline','extrap'),129);
coef_f1 = chebcoeffs(f1);
coef_f1(abs(coef_f1)<eps/100) = 0;
subplot(1,2,1)
semilogy(1:2:129,abs(coef_f1(1:2:end)),'LineWidth',2)
grid on;
set(gca,'fontsize',22);
ylabel('Chebyshev coefficients')
ylim([1e-17,1e0])
yticks([10^(-16),10^(-12),10^(-8),10^(-4),1])
xticks([0,40,80,120])

f2 = chebfun(@(x) interp1(xcol,U1r(:,30),x,'spline','extrap'),129);
coef_f2 = chebcoeffs(f2);
coef_f2(abs(coef_f2)<eps/100) = 0;
subplot(1,2,2)
plot(1:2:129,abs(coef_f2(1:2:end)),'LineWidth',2)
grid on;
set(gca,'fontsize',22);
ylim([0.04,0.12])
xticks([0,40,80,120])

set(gcf, 'Position', [100, 100, 1000, 400]);
if saveresults
saveas(gcf,'figures/ls_coefficients.png')
end



