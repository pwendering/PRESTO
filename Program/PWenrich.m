function PWenrich (KegginfTable, pIds, figprefix)
%Function to analyze the adapted kcats for pathway enrichment
%and plot a barplot
%INPUT:
% - table KegginfTable: A table of Kegg information linked to uniprot
%                       protein IDs of all proteins considered for
%                       correction in PRESTO (universe)
% - cell pIds:  A character cell array giving the uniprot IDs of proteins 
%               whose kcats where changed (sample)
% - char figprefix: A character vector giving the path to save the output
%                   figure
%remove eco01100 for the sake of this analysis
KegginfTable.PathwayInfo(~cellfun(@isempty, KegginfTable.PathwayInfo))= ...
    cellfun(@(x) x(~contains(x,{'eco01100', 'sce01100'})), KegginfTable.PathwayInfo(~cellfun(@isempty, KegginfTable.PathwayInfo)), 'UniformOutput', false);
pIds=pIds(ismember(pIds, KegginfTable.UniprotID));
sample=KegginfTable(ismember(KegginfTable.UniprotID, pIds), :);

% get frequencies of sample pathways
samplePW = horzcat(sample.PathwayInfo{:});
uniqsPW = unique(samplePW);
samplenumPerPW = cellfun(@(x)sum(ismember(samplePW,x)),uniqsPW);
[~,ia] = sort(samplenumPerPW,'descend');

%number of terms for which pvalues should be calculated
n=sum(samplenumPerPW>2);
%
wholePW=horzcat(KegginfTable.PathwayInfo{:});
uniqwPW = unique(wholePW);
numPerPW = cellfun(@(x)sum(ismember(wholePW,x)),uniqwPW);

% Bar graph
if ~isempty(figprefix)

figure
data = samplenumPerPW(ia);
labels = uniqsPW(ia);
pval=nan(n,1);
for i=1:n
    pval(i)=hygecdf(data(i), size(KegginfTable,1), numPerPW(ismember(uniqwPW, labels(i))), length(pIds),'upper');
end
[adj_pval, ~, h]=fdr_BH(pval, 0.05);


b= bar(data(h),'Horizontal',true, 'FaceColor', 'flat');
pv_idx=find(h);
b.CData=b.CData(:,1);
for i=1:sum(h)
    b.CData(i) = adj_pval(pv_idx(i));
    %b.CData(i) = i;
end
yticks(1:sum(h))
yticklabels(regexprep(labels(h), {'sce\d+  ', 'eco\d+  '}, ''))
caxis([0, 0.05]);
xlabel( {'Number of enzymes with', 'corrected kcats value'})
cb=colorbar;
cb.Label.String='adjusted p-value';
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0 0 6.6 3]);
saveas(gcf, [figprefix, 'PRESTO_pathways.svg'])
end
