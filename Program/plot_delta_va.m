close all
rng('default')
% fileName = fullfile('Results','delta_va','deltaVariability_e_coli_lambda_1e-05.mat');
% fileName = fullfile('Results','delta_va','deltaVariability_e_coli_lambda_1e-11.mat');
fileName = fullfile('Results','delta_va','deltaVariability_s_cerevisiae_lambda_1e-07.mat');
% fileName = fullfile('Results','delta_va','deltaVariability_s_cerevisiae_lambda_1e-10.mat');
lambda = str2double(regexp(fileName,'1e[-]*\d+','match'));

load(fileName)

validIdx = true(size(deltaVal));

colors = {'#6999e6','#8c9e4d','#cc5a71'};

[~,ia] = sort(deltaVal(validIdx),'ascend');
origK = corrKcats(validIdx)-deltaVal(validIdx);
minDeltaValid = minDelta(validIdx)/3600;
maxDeltaValid = maxDelta(validIdx)/3600;
deltaValid = deltaVal(validIdx);
sortedDelta = deltaVal(ia);
sortedMin = minDeltaValid(ia);
sortedMax = maxDeltaValid(ia);

figure
hold on

nValidSamples = sum(any(deltaSamplingMat));
if exist('deltaSamplingMat','var') && nValidSamples>5000
    deltaSamplingMat = deltaSamplingMat(ia,:) / 3600;
    
    for i=1:size(deltaSamplingMat,1)
        scatter(i+0.2*rand(size(deltaSamplingMat,2),1)-0.1,deltaSamplingMat(i,:),...
            2,'MarkerEdgeColor',[.7 .7 .7])
    end
    
    count_delta_25_27 = 0;
    count_delta_5_95 = 0;
    for i=1:size(deltaSamplingMat,1)
        p5 = prctile(deltaSamplingMat(i,:),5);
        p25 = prctile(deltaSamplingMat(i,:),25);
        p75 = prctile(deltaSamplingMat(i,:),75);
        p95 = prctile(deltaSamplingMat(i,:),95);
        
        if sortedDelta(i)>=p25 && sortedDelta(i)<=p75
            count_delta_25_27 = count_delta_25_27 + 1;
        end
        
        if sortedDelta(i)>=p5 && sortedDelta(i)<=p95
            count_delta_5_95 = count_delta_5_95 + 1;
        end
        
        if abs(p25-p75)>0
            %             line([i i],[p25 p75],...
            %                 'linewidth',4,'color','b')%colors{1})
            rectangle('Position',[i-.3 p5 .6 p95-p5],'LineWidth',1,...
                'EdgeColor',[.3 .8 0.1])
            rectangle('Position',[i-.2 p25 .4 p75-p25],'LineWidth',1,...
                'EdgeColor','b')
            
        else
            scatter(i,p25, 'Marker', '_', 'LineWidth', 1, 'MarkerEdgeColor', 'b')%colors{1})
        end
    end
    
    scatter(1:sum(validIdx),sortedMin,5,'MarkerEdgeColor',colors{3},...
        'Marker', '*')
else
    fprintf('Number of valid sampled points: %d\n',nValidSamples)
    sortedMin(sortedMin==0) = 1e-6;
    for i=1:numel(deltaVal)
        line([i i],[sortedMin(i) sortedMax(i)],'LineWidth',3,'Color',[.7 .7 .7])
    end
    scatter(1:sum(validIdx),sortedMin,15,'filled','MarkerEdgeColor',colors{3},...
        'MarkerFaceColor',colors{3},'Marker', '^')
    scatter(1:sum(validIdx),sortedMax,15,'filled','MarkerEdgeColor',colors{1},...
        'MarkerFaceColor',colors{1},'Marker', 'v')
    scatter(1:sum(validIdx),sortedDelta,'MarkerEdgeColor',colors{2},...
        'Marker', '_', 'LineWidth', 2)
end


ylabel('\delta [s^{-1}]')
xlabel(sprintf('Enzymes (n = %d)',numel(deltaVal)))
text(.05,.95,sprintf('\\lambda = 10^{%d}',log10(lambda)),...
    'Units','normalized','FontSize',14)
set(gca,'YScale','log','FontSize',14,'LineWidth',1.3,'box','on')

xticks([])

if exist('deltaSamplingMat','var') && nValidSamples > 5000
    h = zeros(3, 1);
    h(1) = plot(NaN,NaN,'*','MarkerFaceColor',colors{3},'MarkerEdgeColor',colors{3},'MarkerSize',10);
    h(2) = plot(NaN,NaN,'s','MarkerFaceColor','none','MarkerEdgeColor','b',...
        'MarkerSize',10,'LineWidth',1.5);
    h(3) = plot(NaN,NaN,'s','MarkerFaceColor','none','MarkerEdgeColor',[.3 .8 0.1],...
        'MarkerSize',10,'LineWidth',1.5);
    h(4) = plot(NaN,NaN,'.','Color',[.7 .7 .7],'MarkerSize',15);
    legend(h,{'min \delta','25% to 75% percentile','5% to 95% percentile',sprintf('sampled \\delta (n=%d)',nValidSamples)},...
        'Box','off','location','southeast')
else
    h = zeros(4, 1);
    h(1) = plot(NaN,NaN,'^','MarkerFaceColor',colors{3},'MarkerEdgeColor',colors{3},'MarkerSize',10);
    h(2) = plot(NaN,NaN,'v','MarkerFaceColor',colors{1},'MarkerEdgeColor',colors{1},'MarkerSize',10);
    h(3) = plot(NaN,NaN,'s','MarkerFaceColor',[.7 .7 .7],'MarkerEdgeColor','none','MarkerSize',10);
    h(4) = plot(NaN,NaN,'_','Color',colors{2},'LineWidth',3);
    legend(h,{'min \delta','max \delta','min/max interval','predicted \delta'},...
        'Box','off','location','north','NumColumns',2)
end
set(gcf,'OuterPosition',[353.6667   88.3333  800.0000  610.6667])


% print('delta_variability_e_coli_1_e_minus_5','-dtiff')
% print('delta_variability_e_coli_1_e_minus_11','-dtiff')
% print('delta_variability_s_cerevisiae_1_e_minus_7','-dtiff')
% print('delta_variability_s_cerevisiae_1_e_minus_10','-dtiff')

% find the number of deltas that span at most one order of magnitude
n_delta_within_1_mgn = sum(log10(range(deltaSamplingMat'))<=1);
fprintf('%d deltas have range of at most one order of magnitude\n',n_delta_within_1_mgn)
n_delta_within_2_mgn = sum(log10(range(deltaSamplingMat'))<=2);
fprintf('%d deltas have range of at most two orders of magnitude\n',n_delta_within_2_mgn)

log_ranges_presto = log10(range(deltaSamplingMat'));
log_ranges_presto(isinf(log_ranges_presto)) = 0;
av_range_small_deltas = mean(log_ranges_presto(sortedDelta<prctile(sortedDelta,50)));
av_range_large_deltas = mean(log_ranges_presto(sortedDelta>=prctile(sortedDelta,50)));
cv_presto = std(deltaSamplingMat,0,2,'omitnan') ./ mean(deltaSamplingMat,2,'omitnan');

residues_presto = sqrt((deltaSamplingMat - mean(deltaSamplingMat,2,'omitnan')).^2);

% determine the average range of kcats in BRENDA
kcat_tab = readtable(fullfile('..','GECKO', 'databases', 'max_KCAT.txt'));
uniq_ec = unique(kcat_tab.Var1);
ranges = cellfun(@(x)range(kcat_tab.Var4(ismember(kcat_tab.Var1,x))),uniq_ec);
mean_brenda = cellfun(@(x)mean(kcat_tab.Var4(ismember(kcat_tab.Var1,x))),uniq_ec);
std_brenda = cellfun(@(x)std(kcat_tab.Var4(ismember(kcat_tab.Var1,x))),uniq_ec);
cv_brenda = mean_brenda ./ std_brenda;
n_per_ec = cellfun(@(x)sum(ismember(kcat_tab.Var1,x)),uniq_ec);
log_ranges_brenda = log10(ranges(n_per_ec>2));
cv_brenda = cv_brenda(n_per_ec>2);
average_range_magn = mean(log_ranges_brenda(~isinf(log_ranges_brenda)));

residues_brenda = cellfun(@(x)...
    sqrt(...
    (kcat_tab.Var4(ismember(kcat_tab.Var1,x))-...
    mean(kcat_tab.Var4(ismember(kcat_tab.Var1,x)),'omitnan')).^2),...
    uniq_ec,'un',0);
residues_brenda = residues_brenda(n_per_ec>2);
residues_brenda = cell2mat(residues_brenda);

fprintf('Average range of kcats per EC number in BRENDA: %.2f orders of magnitude\n',average_range_magn)

figure
h_brenda = histogram(log10(residues_brenda),'Normalization','probability','BinLimits',[-4 6],'NumBins',30);
val_brenda = h_brenda.Values;
h_presto = histogram(log10(residues_presto),'Normalization','probability','BinLimits',[-4 6],'NumBins',30);
val_presto = h_presto.Values;
xval = h_presto.BinEdges(1:end-1) + h_presto.BinWidth / 2;

% brenda
F = fill([xval xval(end:-1:1)],[zeros(size(xval)),val_brenda(end:-1:1)],'k');
F.FaceColor=colors{1};
F.FaceAlpha=0.5;
F.EdgeColor = 'k';
F.LineWidth = 1;
hold on
% presto
F = fill([xval xval(end:-1:1)],[zeros(size(xval)),val_presto(end:-1:1)],'k');
F.FaceColor=colors{2};
F.FaceAlpha=0.5;
F.EdgeColor = 'k';
F.LineWidth = 1;
hold off

% annotation('textbox','Position',[0.2 0.8 .15 0.08], 'String', '$\sqrt{(x-\bar{x})^2}$',...
%     'BackgroundColor', 'white', 'Interpreter', 'latex', 'FitBoxToText', 'on')
xlabel('log_{10}(root squared error)')
ylabel('frequency')
ylim([0 0.2])
legend({'BRENDA','PRESTO'},'box','off')
set(gca,'LineWidth',1.3,'FontName','Arial','FontSize',14)
set(gcf,'OuterPosition',[353.6667  280.3333  487.3333  418.6667])
% exportgraphics(gca,fullfile('Results','delta_va',...
%     'residues_brenda_presto_s_cerevisiae_1_e_minus_7.png'),...
%     'Resolution',300)
% exportgraphics(gca,fullfile('Results','delta_va',...
%     'residues_brenda_presto_e_coli_1_e_minus_5.png'),...
%     'Resolution',300)

fprintf('Average log. eucl. dist. to mean (BRENDA): %.3g\n',...
    mean(log10(residues_brenda(residues_brenda>0)),'omitnan'))
fprintf('Average log. eucl. dist. to mean (PRESTO): %.3g\n',...
    mean(log10(residues_presto),'all','omitnan'))

%% boxplot for EC classes
kcat_tab = readtable(fullfile('..','GECKO', 'databases', 'max_KCAT.txt'));
% kcat_tab = kcat_tab(contains(kcat_tab.Var3,'saccharomyces cerevisiae'),:);
% kcat_tab = kcat_tab(contains(kcat_tab.Var3,'fungi'),:);
% kcat_tab = kcat_tab(contains(kcat_tab.Var3,'escherichia coli'),:);

close
% ec_class = regexp(kcat_tab.Var1,'EC\d.\d+.\d+','match');
ec_class = regexp(kcat_tab.Var1,'EC\d.\d+','match');
% ec_class = regexp(kcat_tab.Var1,'EC\d','match');
ec_class = [ec_class{:}]';
u_ec_class = unique(ec_class);
top_level_prev = 1;
y_limits = [log10(min(kcat_tab.Var4))-1 log10(max(kcat_tab.Var4))+0.5];

for i=1:numel(u_ec_class)
    X = kcat_tab.Var4(strcmp(ec_class,u_ec_class{i}));
    
    u_ec_class{i} = [u_ec_class{i} ' (n=' num2str(numel(X)) ')'];
    
    % get interquartile range
    p = prctile(X,[25 75]);
    
    % plot whiskers
    w_lower = prctile(X,0.35);
    w_upper = prctile(X,99.65);
    line(...
        [i i],...
        log10([w_lower w_upper]),...
        'linewidth',1,...
        'Color','k'...
        )
    
    line(...
        [i-.1 i+.1],...
        log10([w_lower w_lower]),...
        'linewidth',1,...
        'Color','k'...
        )
    line(...
        [i-.1 i+.1],...
        log10([w_upper w_upper]),...
        'linewidth',1,...
        'Color','k'...
        )
    
    % plot interquartile range
    rectangle('Position',[i-.2 log10(p(1)) .4 log10(p(2))-log10(p(1))],'LineWidth',1,...
        'EdgeColor',colors{1}, 'FaceColor', [.7 .7 .7])
    hold on
    
    % plot median
    plot(...
        i,...
        log10(median(X,'omitnan')),...
        'o',...
        'MarkerSize',5,...
        'MarkerFaceColor',colors{3},...
        'MarkerEdgeColor',colors{3}...
        )
    
    
    % plot outliers
    outlier_idx = X < w_lower | X > w_upper;
    plot(...
        repelem(i,sum(outlier_idx),1),...
        log10(X(outlier_idx)),...
        '*',...
        'MarkerEdgeColor','r',...
        'MarkerSize',1 ...
        )
    
    % separate main EC classes
    top_level_curr = str2double(erase(strtok(u_ec_class{i},'.'),'EC'));
    if top_level_curr > top_level_prev
        line(...
            [i-.5 i-.5],...
            y_limits,...
            'linewidth',2 ...
            )
    end
    top_level_prev = top_level_curr;
    
%     text(...
%         i,...
%         y_limits(1)+.5,...
%         num2str(numel(X)),...
%         'HorizontalAlignment','center',...
%         'FontSize',8 ...
%         )
    
    fprintf(...
        'Median (+/- sd) %s: %g (+/- %g) (n=%d)/s\n',...
        u_ec_class{i},...
        median(X,'omitnan'),...
        std(X,'omitnan'),...
        numel(X)...
        )
end
hold off
xlim([0 numel(u_ec_class)+1])
xticks(1:numel(u_ec_class));
xticklabels(u_ec_class)
xtickangle(90)
ylim(y_limits)
ylabel(...
    'log_{10} k_{cat} [s^{-1}]',...
    'FontSize',14 ...
    )
set(gca,...
    'FontName','Arial',...
    'box','on',...
    'linewidth',1.3,...
    'ticklength',[0.005 0.005]...
    )

set(gcf,...
    'OuterPosition',[0.0123    0.1897    1.1987    0.4887]*1000 ...
    )

% exportgraphics(gca,'kcat_ec_class.png','Resolution',300)
% exportgraphics(gca,'kcat_ec_class_sce.png','Resolution',300)
% exportgraphics(gca,'kcat_ec_class_eco.png','Resolution',300)