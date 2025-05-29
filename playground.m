close all;
clear;

nn = 5:10:1000;
mm = 5:10:1000;
svals_P = zeros(length(nn),length(mm));
for ii = 1:length(mm)
    ii
    for jj = 1:length(nn)
    m = mm(ii); n = nn(jj);
    P = uni2cheb([-1,1],n,m);
    %F = Vals2ChebCoeffsMat(m);
    svals_P(ii,jj) = norm(P);
    end
end


mm = 5:10:1000;
svals = zeros(length(mm),1);
for ii = 1:length(mm)
    m = mm(ii);
    F = Vals2ChebCoeffsMat(m);
    svals(ii) = norm(F);
end

loglog(mm,mm.^-1,'r'), hold on, loglog(mm,svals,'.')

clear;  d=3; A=1;
saveresults = false;

n1 = 2048;Pn = 500;
filename = ['data\Data_256-2048\n' num2str(n1) '_Pn' num2str(Pn) '.mat'];

load(filename)

% get the CP tensor
xi = LRed.LAM1C./Hunif^3;
Rl = length(xi);
U = {LRed.CU1, LRed.CU2, LRed.CU3};

m = [129,129,129];
tol = 1e-7;

xx = (-1+1/n1):2/n1:(1-1/n1);


f = @(x) interp1(xx,U{1}(:,100),x,'spline','extrap');
fa = chebfun(f,[-1,1],'splitting','on');

f = @(x) interp1(xx,U{1}(:,300),x,'spline','extrap');
fb = chebfun(f,[-1,1],'splitting','on');

f = @(x) interp1(xx,3*U{1}(:,100) + -8*U{1}(:,300),x,'spline','extrap');
fc = chebfun(f,[-1,1],'splitting','on');

aa = chebtech2.vals2coeffs(interp1(xx,U{1}(:,147),chebpts(125,[-1,1]),'spline','extrap'));
bb = chebtech2.vals2coeffs(interp1(xx,U{1}(:,289),chebpts(125,[-1,1]),'spline','extrap'));
cc = chebtech2.vals2coeffs(interp1(xx,U{1}(:,147) + U{1}(:,289),chebpts(125,[-1,1]),'spline','extrap'));


mm = round(2.^(0:0.3:log2(length(fa))+0.3));
ee = 0*mm;
for jj=1:length(mm)
    n = mm(jj); fn=chebfun(f,n); ee(jj) = norm(fa-fn,inf);
end
loglog(mm,mm.^-3,'r'), hold on, loglog(mm,ee,'.')