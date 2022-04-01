clear;clc;close all

% load('Results/lambda_cv/lambda_fitting_s_cerevisiae_batch.mat')
load('Results/lambda_cv/lambda_fitting_e_coli_batch.mat')

sumsDelta(sumsDelta==0) = NaN;

lwd = 3;
fsz = 14;
colors = {'#6999e6','#8c9e4d','#cc5a71'};

avRelErr = mean(relErr,2,'omitnan');
avSumsDelta = mean(sumsDelta,2,'omitnan');
avErrVar = mean(errVar,2,'omitnan');

figure
tiledlayout(3,2)
nexttile([3 1])
xlabel('log_{10} \lambda','FontSize',fsz)
xlim([min(lambdaParams) max(lambdaParams)])

yyaxis left

err_standard = (relErr-min(relErr)) ./ range(relErr);
log_delta_standard = (log10(sumsDelta)-min(log10(sumsDelta))) ./ range(log10(sumsDelta));

lambdaScores = mean(err_standard.*log_delta_standard,2,'omitnan');
plot(lambdaParams, lambdaScores,'LineWidth',lwd,'Color',colors{3},...
    'LineStyle',':')
hold on

plot(lambdaParams,avRelErr,'LineWidth',lwd,'Color',colors{1},'LineStyle', '-')
hold on
for i=1:numel(lambdaParams)
        plot(lambdaParams(i),avRelErr(i),'.','Color',colors{1},'MarkerSize',25)
end

% plot for relative error sd
avErrStd = sqrt(avErrVar);
x = [lambdaParams(1);lambdaParams(1);lambdaParams(2:end)';...
    lambdaParams(end:-1:2)'];
y = [avRelErr(1)+avErrStd(1);avRelErr(1)-avErrStd(1);avRelErr(2:end)-avErrStd(2:end);...
    avRelErr(end:-1:2)+avErrStd(end:-1:2)];
x(isnan(y)) = [];
y(isnan(y)) = [];
F = fill(x,y,'k');
F.FaceColor=colors{1};
F.FaceAlpha=0.3;
ylabel('Relative error to experimental \mu','FontSize',fsz)

yyaxis right
plot(lambdaParams,log10(avSumsDelta),'LineWidth',lwd,'Color',colors{2});
plot(lambdaParams,log10(avSumsDelta),'.','MarkerSize',25,'Color',colors{2});
ylabel('log_{10} \Sigma \delta','FontSize',fsz)

ax = gca;
ax.YAxis(1).Color = colors{1};
ax.YAxis(1).FontWeight = 'bold';
ax.YAxis(2).Color = colors{2};
ax.YAxis(2).FontWeight = 'bold';
ax.LineWidth = 1.3;
ax.XScale = 'log';

% draw vertical line at optimal lambda parameter
d1 = gradient(lambdaScores,log10(lambdaParams));
d2 = gradient(d1,log10(lambdaParams));
disp([log10(lambdaParams);sign(d2')])
ifp_pos_1 = strfind(sign(d2'),[-1 1])+1;

% line([1e-7 1e-7],get(gca,'YLim'),'color','k','linestyle','--') % yeast
line([1e-5 1e-5],get(gca,'YLim'),'color','k','linestyle','--') % E. coli

% draw vertical line at second lambda parameter
d1 = gradient(avSumsDelta,log10(lambdaParams));
d2 = gradient(d1,log10(lambdaParams));
disp([avSumsDelta';sign(d2')])
ifp_pos_2 = strfind(sign(d2'),[-1 1])+1;

% line([1e-10 1e-10],get(gca,'YLim'),'color','k','linestyle','--') % yeast
line([1e-11 1e-11],get(gca,'YLim'),'color','k','linestyle','--') % E. coli

h = zeros(2,1);
h(1) = plot(NaN,NaN,'LineWidth',lwd,'Color',colors{3});
h(2) = plot(NaN,NaN,'Marker','s','MarkerSize',15,'MarkerFaceColor',colors{1},...
    'MarkerEdgeColor',colors{1},'LineStyle','none');
legend(h,{'e_r^* \cdot log_{10}(\Sigma\delta)^*',...
    'e_r \pm sd(CV)'},'FontSize',14,'Box','off',...
    'location','best')

text(-.13,.97,'a','Units','normalized','FontSize',14,'FontName','Arial',...
    'FontWeight', 'bold')
% calculate the number of corrected kcats per lambda and iteration
av_n_corr = mean(cellfun(@numel,corrKcatProts),2);

% calculate Jaccard distance between corrected kcat sets per lambda and iteration
JD_within_iterations = zeros(1,size(corrKcatProts,1));
for i=1:size(corrKcatProts,1)
    jd = nan(1,size(corrKcatProts,2));
    for j=1:size(corrKcatProts,2)-1
        for k=j+1:size(corrKcatProts,2)
            jd(j,k) = 1 - numel(intersect(corrKcatProts{i,j},corrKcatProts{i,k})) / ...
            numel(union(corrKcatProts{i,j},corrKcatProts{i,k}));
        end
    end
    JD_within_iterations(i) = mean(jd,'all','omitnan');
end

uniq_sets = cell(1,size(corrKcatProts,1));
for i=1:size(corrKcatProts,1)
    uniq_sets{i} = unique(vertcat(corrKcatProts{i,:}));
end

JD_between_uniq_sets = zeros(numel(uniq_sets));
for i=1:numel(uniq_sets)-1
    for j=i+1:numel(uniq_sets)
        JD_between_uniq_sets(i,j) = 1 - numel(intersect(uniq_sets{i},uniq_sets{j})) / ...
            numel(union(uniq_sets{i},uniq_sets{j}));
        JD_between_uniq_sets(j,i) = JD_between_uniq_sets(i,j);
    end
end

nexttile([1 1])
plot(mean(avJD,2),'Color',colors{1},'LineWidth',3,'XData',-14:-1)
xticks(-14:-1)
ylabel('Jaccard distance','FontSize',12,'FontName','Arial')
text(-.2,.97,'b','Units','normalized','FontSize',14,'FontName','Arial',...
    'FontWeight', 'bold')

nexttile([2 1])
% colormapeditor
% save('presto_colormap','presto_cm')
load('presto_colormap')
imagesc(log10(lambdaParams),log10(lambdaParams),JD_between_uniq_sets)
xlabel('log_{10} \lambda', 'FontName', 'Arial', 'FontSize', 14)
ylabel('log_{10} \lambda', 'FontName', 'Arial', 'FontSize', 14)
set(gca,'Colormap',presto_cm)
colorbar
text(-.2,.97,'c','Units','normalized','FontSize',14,'FontName','Arial',...
    'FontWeight', 'bold')
text(1.18,0.5,'Jaccard distance','HorizontalAlignment','center',...
    'Rotation',90,'units','normalized','FontSize', 14,...
    'FontName','Arial','Color',[.2 .2 .2])

set(gcf,'OuterPosition',1000*[-1.1017    0.0310    1.0720    0.5633])

% exportgraphics(gcf,'Results/lambda_cv/Fig_S1.tiff','resolution',300)
% exportgraphics(gcf,'Results/lambda_cv/Fig_S6.tiff','resolution',300)