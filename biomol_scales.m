close all;

clear;  d=3; A=1;
saveresults = true;
Plist = 100:100:600;
mlist = [9,17,33,65,129,257,513,1025];
num_P = length(Plist);
num_m = length(mlist);

svals = zeros(num_P,129);
TuckerRanks = zeros(num_P, num_m);
tol = 1e-7;
label_names = cell(1,num_P);
for ii = 1:num_P
    n1 = 2048;Pn = Plist(ii);
    label_names{ii} = ['N = ' num2str(Pn)];
    for jj = 5:5%1:num_m
        m = mlist(jj) * ones(1,d);
        filename = ['data/Data_256-2048/n' num2str(n1) '_Pn' num2str(Pn) '.mat'];

        load(filename,"LRed","Hunif")

        % get the CP tensor
        xi = LRed.LAM1C./Hunif^3;
        Rl = length(xi);
        U = {LRed.CU1, LRed.CU2, LRed.CU3};

        % compute ChebTuck
        [f,ff] = ChebTuck({xi,U},m,[],tol);
        TuckerRanks(ii,jj) = max(size(f.core));
        if m == 129
            svals(ii,:) = svd(ff{3});
        end
    end
end


figure(1)
semilogy(svals','linewidth',2)
ylabel('singular values')
yticks([10^(-20),10^(-16),10^(-12),10^(-8),10^(-4),1])
set(gca,'fontsize',20);
axis tight
grid on
set(gcf,'position',[100,100,600,500])
if saveresults
saveas(gcf,'figures/svals_CU.png')
end

figure(2)
semilogy(10:40,svals(:,10:40)','linewidth',2)
yticks([10^(-20),10^(-16),10^(-12),10^(-8),10^(-4),1])
set(gca,'fontsize',20);
legend(label_names,'Location','northeastoutside')
axis tight
grid on
set(gcf,'position',[100,100,722,500])
if saveresults
saveas(gcf,'figures/svals_CU_zoom.png')
end



