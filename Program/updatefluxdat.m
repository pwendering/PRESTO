function updatefluxdat
%Function to compare and update Chen et al. dat with publication
%supplements
chendat=readtable('Data/Chen/ProteomicsFlux.xlsx', 'Sheet', 'Flux');
lahtveedat=readtable('Data/Chen/Lahtvee2017sup.xlsx');
lahtveedat=lahtveedat(:, [1,2, 9:14]);
%add labels from chendat
lahtveedat.Properties.VariableNames(3:end)=chendat.note(2:7);
%rename conditions
lahtveedat(ismember(lahtveedat.(2), {'Ref'}),2)={'REF'};
lahtveedat(contains(lahtveedat.(2), 'g/L'), 2)=cellfun(@(x) ['EtOH' x(1:2)], lahtveedat.(2)(contains(lahtveedat.(2), 'g/L')), 'UniformOutput', false);
lahtveedat(contains(lahtveedat.(2), 'mM'), 2)=cellfun(@(x) ['Osmo0' x(1)], lahtveedat.(2)(contains(lahtveedat.(2), 'mM')), 'UniformOutput', false);
lahtveedat(contains(lahtveedat.(2), 'oC'), 2)=cellfun(@(x) ['Temp' x(1:2)], lahtveedat.(2)(contains(lahtveedat.(2), 'oC')), 'UniformOutput', false);
lahtveedat=lahtveedat(:,2:end);
lahtveedat.(1)=strcat('Lahtvee2017_',lahtveedat.(1));
lahtveedat.Properties.VariableNames(1)={'CondNames'};
if any(~ismember(lahtveedat.CondNames, chendat.Properties.VariableNames), 'all')
    error('Parsing of lahtvee condition names failed. Aborting...')
end
%calculate condition means
nlahtveedat=nan(sum(ismember(chendat.Properties.VariableNames, lahtveedat.CondNames)),size(lahtveedat, 2)-1);
incCondNames=chendat.Properties.VariableNames(ismember(chendat.Properties.VariableNames, lahtveedat.CondNames));
for i=1:length(incCondNames)
    for j=1:size(nlahtveedat,2)
    nlahtveedat(i,j)=mean(table2array(lahtveedat(ismember(lahtveedat.CondNames, incCondNames(i)), j+1)));
    end
end
%Parse into chendat format
nlahtveedat=array2table(nlahtveedat'*-1, 'VariableNames', incCondNames, 'RowNames', lahtveedat.Properties.VariableNames(2:end));

%Parse Yu data
yudat=readtable('Data/Chen/Yu2021sup.xlsx', 'Sheet', 'Supplementary File 1e');
%remove duplicated nitrogen row
yudat(:,7)=[];
nyudat=table2array(yudat(:,2:end))';
%Adapt Chen et al labelling
if any(table2array(chendat(2,contains(chendat.Properties.VariableNames, 'Yu2021')))+nyudat(1,:)>0.001, 'all')
    error('Supplemental table of Yu2021 has different sample order as chen et al. Aborting...')
end
nyudat=array2table(nyudat, 'VariableNames', chendat.Properties.VariableNames(contains(chendat.Properties.VariableNames, 'Yu2021')), 'RowNames', yudat.Properties.VariableNames(2:end));

%adapt chen et al ordering
nyudat.Properties.RowNames=cellfun(@(x) x(2:4), nyudat.Properties.RowNames, 'UniformOutput', false);
nyudat.Properties.RowNames(contains(nyudat.Properties.RowNames, 'O2_'))={'O2'};
nyudat.Properties.RowNames(contains(nyudat.Properties.RowNames, 'EtO'))={'Eth'};
bait_chendat=cell(size(chendat.note));
for i=1:length(bait_chendat)
    if length(chendat.note{i})>2
        bait_chendat(i)=extractBefore(chendat.note(i), 4);
    else
        bait_chendat(i)=chendat.note(i);
    end
end
[match, idx]=ismember(bait_chendat,nyudat.Properties.RowNames);
tmpdat=nyudat(idx(match), :);
tmpdat=[tmpdat; nyudat(~ismember(nyudat.Properties.RowNames, bait_chendat),:)];
nyudat=tmpdat;
%adapt signs
%Flux sign from chen et al and positive for succinat
yu_sign=[-1, 1, 1, 1, -1, 1, -1, -1, -1, -1, 1, 1];
tmpdat=table2array(nyudat);
for i=1:size(tmpdat, 1)
    tmpdat(i,:)=tmpdat(i,:)*yu_sign(i);
end
nyudat=array2table(tmpdat, 'VariableNames', nyudat.Properties.VariableNames, 'RowNames', nyudat.Properties.RowNames);
%check differences
tmpdat(isnan(tmpdat))=0;
tmpchendat=table2array(chendat(2:end,contains(chendat.Properties.VariableNames, 'Yu2021')));
tmpchendat(isnan(tmpchendat))=0;
test=tmpchendat-tmpdat(1:(end-1), :);
disp('Molecules which have different exchange fluxes in Yu et al 2021 and Chen et al:')
disp(nyudat.Properties.RowNames(sum(test,2)>0.01))

%parse Yu2020 data
yu20dat=readtable('Data/Chen/Yu2020sup.xlsx');
yu20dat=array2table(table2array(yu20dat(:, 2:end))', 'VariableNames', yu20dat.Var1, 'RowNames', yu20dat.Properties.VariableNames(2:end));
yu20dat.Properties.VariableNames=regexprep(yu20dat.Properties.VariableNames, '-|/|=| |ep', '');
yu20dat.Properties.VariableNames=strcat('Yu2020_',regexprep(yu20dat.Properties.VariableNames, 'r', 'R'));
if any(~ismember(yu20dat.Properties.VariableNames, chendat.Properties.VariableNames), 'all')
    error('Conditions names from Yu2020 et al could not be matched with Chen et al. Aborting...')
end
%adapt chendat format
yu20dat.Properties.RowNames=extractBetween(yu20dat.Properties.RowNames, 'r', '_');
yu20dat.Properties.RowNames(cellfun(@length, yu20dat.Properties.RowNames)>2)=extractBefore(yu20dat.Properties.RowNames(cellfun(@length, yu20dat.Properties.RowNames)>2), 4);
match=[];
idx=[];

for i=1:length(bait_chendat)
    TF=strcmpi(bait_chendat(i), yu20dat.Properties.RowNames);
    hits=sum(TF);
    if sum(TF)>1
        error('ambigous naming between yu2020 and chendat detected. Aborting...')
    elseif hits==1
        idx=[idx find(TF)];
    else
        idx=[idx hits];
    end
    match=[match hits];
end
tmpdat=yu20dat(idx(logical(match)),:);
tmpdat=[tmpdat; yu20dat(setdiff(1:size(yu20dat,1), idx(logical(match))), :)];
yu20dat=tmpdat;
tmpdat=table2array(tmpdat);
%add flux sign from chendat 
yu20_sign=[-1, 1, 1,1,-1, 1, 1, 1];
for i = 1:size(tmpdat,1)
    tmpdat(i,:)=tmpdat(i,:)*yu20_sign(i);
end
yu20dat=array2table(tmpdat, 'VariableNames', yu20dat.Properties.VariableNames, 'RowNames', yu20dat.Properties.RowNames);
%compare levels
tmpchendat=table2array(chendat(logical(match), contains(chendat.Properties.VariableNames, 'Yu2020')));
tmpchendat(isnan(tmpchendat))=0;
test=tmpdat(1:sum(match),:)-tmpchendat;
disp('Molecules which have different exchange fluxes in Yu et al 2020 and Chen et al:')
disp(yu20dat.Properties.RowNames(sum(test,2)>0.01))


%% Assemble a new Table
%create identical reaction indices
nlahtveedat.Properties.RowNames(cellfun(@length, nlahtveedat.Properties.RowNames)>3)=extractBefore(nlahtveedat.Properties.RowNames(cellfun(@length, nlahtveedat.Properties.RowNames)>3), 4);
yu20dat.Properties.RowNames=cellfun(@(x) [upper(x(1)), x(2:end)], yu20dat.Properties.RowNames, 'UniformOutput', false);
supdat=struct('lahtvee', nlahtveedat, 'yu2021', nyudat, 'yu2020', yu20dat);
new_rxnn=bait_chendat;
fn=fieldnames(supdat);
for i=1:length(fieldnames(supdat))
    new_rxnn=[new_rxnn; supdat.(fn{i}).Properties.RowNames(~ismember(supdat.(fn{i}).Properties.RowNames, new_rxnn))'];
end
%Assemble table
%add DiBartolomeo from chendat
[match, idx]=ismember(new_rxnn, bait_chendat);
tmpdat=nan(length(new_rxnn), sum(contains(chendat.Properties.VariableNames, 'DiBartolomeo')));
tmpdat(match,:)=table2array(chendat(idx(logical(match)), contains(chendat.Properties.VariableNames, 'DiBartolomeo')));
newcomp=[table(new_rxnn), array2table(tmpdat, 'VariableNames', chendat.Properties.VariableNames(contains(chendat.Properties.VariableNames, 'DiBartolomeo')))];
%Add data from other supplemental tables
for i=1:length(fieldnames(supdat))
    tmptab=supdat.(fn{i});
    tmpdat=nan(length(new_rxnn), size(tmptab, 2));
    [match, idx]=ismember(new_rxnn, tmptab.Properties.RowNames);
    tmpdat(match,:)=table2array(tmptab(idx(logical(match)), :));
    newcomp=[newcomp, array2table(tmpdat, 'VariableNames', tmptab.Properties.VariableNames)];
end
%add growth rates
[match, idx]=ismember(newcomp.Properties.VariableNames, chendat.Properties.VariableNames);
newcomp(1, logical(match))=chendat(1, idx(logical(match)));
%compare full tables
tmpchendat=table2array(chendat(:, idx(logical(match))));
tmpdat=table2array(newcomp(ismember(newcomp.new_rxnn, bait_chendat), 2:end));
tmpchendat(isnan(tmpchendat))=0;
tmpdat(isnan(tmpdat))=0;
test=tmpchendat-tmpdat;
disp('Molecules which have different exchange fluxes in final table and Chen et al:')
disp(newcomp.new_rxnn(sum(test,2)>0.01))
%Change to default for unknowns
tmpdat=table2array(newcomp(:,2:end));
for i=find(sum(isnan(tmpdat),2))'
    switch newcomp.new_rxnn{i}
        case 'O2'
            tmpdat(i,isnan(tmpdat(i,:)))=-1000;
        case {'CO2', 'Pyr', 'Suc'}
            tmpdat(i,isnan(tmpdat(i,:)))=+1000;
        case 'NH4'
            for j=find(isnan(tmpdat(i,:)))
                if contains(newcomp.Properties.VariableNames(j), 'Yu2021')
                    tmpdat(i,j)=0;
                else
                    tmpdat(i,j)=-1000;
                end
            end
        case {'Gln', 'Phe', 'Ile'}
            tmpdat(i,isnan(tmpdat(i,:)))=0;
    end
end
%convert back to chendat note and exchangeRxn fiels
exchangeRxn=cell(size(newcomp,1),1);
note=cell(size(newcomp,1),1);
[match,idx]=ismember(newcomp.new_rxnn, bait_chendat);
exchangeRxn(match)=chendat.exchangeRxn(idx(logical(match)));
exchangeRxn(ismember(newcomp.new_rxnn, {'Suc'}))={'r_2056'};
note(match)=chendat.note(idx(logical(match)));
note(ismember(newcomp.new_rxnn, {'Suc'}))={'Succinate'};
newcomp=[table(exchangeRxn, note), array2table(tmpdat, 'VariableNames', newcomp.Properties.VariableNames(2:end))];

writetable(newcomp, 'Data/Chen/ProteomicsFlux.xlsx', 'Sheet', 'Flux', 'WriteMode', 'overwritesheet')
end
