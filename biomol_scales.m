function biomol_scales(saveresults)
d=3; tol = 1e-7;
Plist = 100:100:600;
mlist = [9,17,33,65,129,257,513,1025];
num_P = length(Plist);
num_m = length(mlist);

svals = zeros(num_P,129);
TuckerRanks = zeros(num_P, num_m);
CompressionRatio = zeros(num_P, num_m);
CompressionPct = zeros(num_P, num_m);
CPParams = zeros(num_P, num_m);
ChebTuckParams = zeros(num_P, num_m);
label_names = cell(1,num_P);
for ii = 1:num_P
    n1 = 2048;Pn = Plist(ii);
    label_names{ii} = ['N = ' num2str(Pn)];
    for jj = 1:num_m
        m = mlist(jj) * ones(1,d);
        filename = ['data/Data_256-2048/n' num2str(n1) '_Pn' num2str(Pn) '.mat'];

        load(filename,"LRed","Hunif")

        % get the CP tensor
        xi = LRed.LAM1C./Hunif^3;
        U = {LRed.CU1, LRed.CU2, LRed.CU3};
        n = size(U{1}, 1);  % grid size
        R = length(xi);      % CP rank

        % compute ChebTuck
        [f,ff] = ChebTuck({xi,U},m,[],tol);

        % Extract Tucker ranks from core tensor
        core_size = size(f.core);
        r1 = core_size(1); r2 = core_size(2); r3 = core_size(3);
        m_val = m(1);  % all m_i are the same

        TuckerRanks(ii,jj) = max(core_size);

        % Compute parameter counts
        % CP format: 3*n*R (three factor matrices, each n × R)
        cp_params = 3*n*R;

        % ChebTuck format: core tensor + factor matrices
        % core: r1*r2*r3, factors: m*(r1+r2+r3)
        chebtuck_params = r1*r2*r3 + m_val*(r1 + r2 + r3);

        CPParams(ii,jj) = cp_params;
        ChebTuckParams(ii,jj) = chebtuck_params;
        CompressionRatio(ii,jj) = chebtuck_params / cp_params;
        CompressionPct(ii,jj) = 100 * CompressionRatio(ii,jj);

        if m_val == 129
            svals(ii,:) = svd(ff{3});
        end
    end
end


figure(11)
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

figure(12)
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

% Display the TuckerRanks table
fprintf('\nTable 4: Tucker ranks\n');
fprintf('================================================================\n');
fprintf('  N\\m');
fprintf('%6d', mlist);
fprintf('\n');
fprintf('----------------------------------------------------------------\n');
for ii = 1:num_P
    fprintf('%5d', Plist(ii));
    fprintf('%6d', TuckerRanks(ii,:));
    fprintf('\n');
end
fprintf('================================================================\n');

% Display the Compression Percentage table
fprintf('\nCompression percentage (ChebTuck params / CP params) * 100\n');
fprintf('================================================================\n');
fprintf('  N\\m');
fprintf('%8d', mlist);
fprintf('\n');
fprintf('----------------------------------------------------------------\n');
for ii = 1:num_P
    fprintf('%5d', Plist(ii));
    for jj = 1:num_m
        fprintf('%8.1f', CompressionPct(ii,jj));
    end
    fprintf('\n');
end
fprintf('================================================================\n');

% Display parameter counts for m=129
m_idx = find(mlist == 129);
if ~isempty(m_idx)
    fprintf('\nTable 5: Parameter counts at m=129:\n');
    fprintf('================================================================\n');
    fprintf('    N     CP params     ChebTuck params   Compression(%%)   Savings(%%)\n');
    fprintf('----------------------------------------------------------------\n');
    for ii = 1:num_P
        pct = CompressionPct(ii, m_idx);
        savings = 100 - pct;
        fprintf('%5d %13.0f %17.0f %15.1f %12.1f\n', ...
            Plist(ii), CPParams(ii,m_idx), ChebTuckParams(ii,m_idx), pct, savings);
    end
    fprintf('================================================================\n');
end

% Show savings for m=129
if ~isempty(m_idx)
    fprintf('\nCompression savings at m=129:\n');
    for ii = 1:num_P
        pct = CompressionPct(ii, m_idx);
        savings = 100 - pct;
        fprintf('N=%3d: %.1f%% parameters saved (ChebTuck = %.1f%% of CP)\n', Plist(ii), savings, pct);
    end
end
end

