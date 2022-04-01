function[loginf, kcat_mat] = getloginf(logfile, org_name)
%Function to extract model fitting data created by get_GKOmod.m 
%Input: 
% - char logfile: absolute path to the log file written by get_GKOmod.m
% - char org: model organism
%Output:
% - tab loginf: A table containing fitted GAM, the Error of GAM fitting and
%                    fitted Sigma per model
%- tab kcat_mat: A table where the first three columns give information
%           about a kcat changed during automatic relaxation and the other
%           columns indicate if the given kcat where changed in the
%           respective model

switch org_name
    case 'saccharomyces cerevisiae'
        configuration_yeast
        condNames=readChenetal(topDir);
        [loginf, kcat_tabs]=readlogfileyeast(logfile, orgBasename, condNames);
    case 'escherichia coli'
        configuration_ecoli
        condNames=readDavidi2016([], topDir);
        [loginf,kcat_tabs]=readlogfileecoli(logfile, orgBasename, condNames);
end


%comile matrix of kcat adaptions - n= rxn with adaptions m= model, if
%adaption contained in model (n,m)= 1 else 0
kcat_mat=zeros(0, length(loginf.condNames));
%adpation names
kcat_mat_row=strings(0,3);
for i=1:length(kcat_tabs)
    kcat_tab=string(kcat_tabs{i});
    for j=1:size(kcat_tab,1)
        %add new edge to rows
        if ~ismember(kcat_tab(j,:), kcat_mat_row, 'rows')
            kcat_mat_row=[kcat_mat_row; kcat_tab(j,:)];
            temp_row=zeros(1, length(loginf.condNames));
            temp_row(i)=1;
            kcat_mat=[kcat_mat; temp_row];
        else
            %add entry in existing row
            idx=find(ismember(kcat_mat_row, kcat_tab(j,:), 'rows'));
            kcat_mat(idx, i)=1;
        end
    end
end
disp('Model fitting information:')
disp(loginf)
%disp('Kcat adaption info:')
kcat_mat=[array2table(kcat_mat_row, 'VariableNames', {'Protein', 'Reaction#', 'Name'}), array2table(kcat_mat, 'VariableNames', loginf.condNames)];
%disp(kcat_mat)
end

function[loginf, kcat_tabs]= readlogfileyeast(logfile, orgBasename, condNames)
%% Import Info from log of GECKO run for yeast
%INPUT:
% - char logfile: absolute path to the log file written by get_GKOmod.m
% - char orgBasename: prefix for GECKO models build (used to detect data in
%                     logfile
% - cell condNames: cell characters vector of condition names for which
%                   models where build (qsanity check to make sure correct
%                   log file is used)
lf=fopen(logfile);
%logtxt=textscan(lf, '%s' );
%fclose(lf);
PAR_M=nan(0, 3); %GAM fitting parameters
Mod_n=cell(0); %model names
kcat_tabs=cell(0); %kcat adpations;
logtxt=fgetl(lf);
% m=0
% o=0
while ischar(logtxt)
    if contains(logtxt, 'Fitted GAM')
        %Initialize parameter field in this conditions since GAM fitting
        %takes place before sigma fitting
        if exist('temp_PAR', 'var')&&any(isnan(temp_PAR))
            pause 
            error('NaN values in temp_PAR detected. Check log file format')
        end
        temp_PAR=nan(1,3);
        temp_PAR(1:2)=sscanf(logtxt, 'Fitted GAM = %f -> Error = %f', [1,2]);
        %since this is the first entry found for a new model intialize kcat
        %adaption table
        if exist('kcat_tab', 'var')
            kcat_tabs=[kcat_tabs {kcat_tab}];
        end
        kcat_tab=cell(0,3);
    elseif contains(logtxt, 'Protein:')
        %Import information of relaxed kcats
        kcat_entry=cell(1,3);
        kcat_entry(1)=extractBetween(logtxt, 'Protein:', ' Rxn#');
        kcat_entry(2)=extractBetween(logtxt, 'Rxn#:', ' name:');
        kcat_entry{3}=extractAfter(logtxt, 'name: ');
        %add to table
        kcat_tab=[kcat_tab; kcat_entry];
    elseif contains(logtxt, 'Saving adp')
        %import corresponding condition names
        temp_n=extractBetween(logtxt, [orgBasename '_'], '_batch');
        Mod_n=[Mod_n , temp_n];
        if ~ismember(temp_n, condNames)
            error(['Found model name not matching with conditions:' ...
               temp_n])
        end
        %m=m+1;
    elseif contains(logtxt, 'Sigma')
        temp_PAR(3)=sscanf(logtxt, 'Sigma factor (fitted for growth on glucose): %f');
        PAR_M=[PAR_M;temp_PAR];
        %o=o+1;
    end
    logtxt=fgetl(lf);
end
fclose(lf);
%compile table with condition names GAM and sigma
loginf=[cell2table(Mod_n', 'VariableNames', {'condNames'}), array2table(PAR_M, 'VariableNames', {'GAM', 'Error', 'Sigma'})] ;
end

function[loginf, kcat_tabs] = readlogfileecoli(logfile, orgBasename, condNames)
%% Import Info from log of GECKO run for e.coli
lf=fopen(logfile);
%logtxt=textscan(lf, '%s' );
%fclose(lf);
PAR_M=nan(0, 1); %GAM fitting parameters
Mod_n=cell(0); %model names
kcat_tabs=cell(0); %kcat adpations;
logtxt=fgetl(lf);
% m=0
% o=0
while ischar(logtxt)
    if contains(logtxt, 'Protein:')
        %Import information of relaxed kcats
        kcat_entry=cell(1,3);
        kcat_entry(1)=extractBetween(logtxt, 'Protein:', ' Rxn#');
        kcat_entry(2)=extractBetween(logtxt, 'Rxn#:', ' name:');
        kcat_entry{3}=extractAfter(logtxt, 'name: ');
        %if kcat tab and temp_Par does not exist yet
        if ~(exist('kcat_tab', 'var'))
            %initialize temp_PAR array since kcat refitting is the first log
            %entry for new model
            temp_PAR=nan(1,1);
            kcat_tab=kcat_entry;
        else
            kcat_tab=[kcat_tab; kcat_entry];
        end
    elseif contains(logtxt, 'Saving adp')
        %check if info is parsed in correct order, this is last point
        if exist('temp_PAR', 'var')&&any(isnan(temp_PAR))
            pause 
            error('NaN values in temp_PAR detected. Check log file format')
        end
        %import corresponding condition names
        temp_n=extractBetween(logtxt, [orgBasename '_'], '_batch');
        Mod_n=[Mod_n , temp_n];
        %check if condition name is in experimental data
        if ~ismember(temp_n, condNames)
            error(['Found model name not matching with conditions:' ...
               temp_n])
        end
        %if kcats were adapted save adaptions
        if exist('kcat_tab', 'var')
            kcat_tabs=[kcat_tabs, {kcat_tab}];
        else 
            kcat_tabs=[kcat_tabs, cell(1,1)];
        end
        clear kcat_tab temp_PAR
        %m=m+1;
    elseif contains(logtxt, 'Sigma')
        %in case there is no kcat fitting
        if ~(exist('temp_PAR', 'var'))
            warning('Could not find kcat adaption log. If kcats were adapted the log is corrupted')
            temp_PAR=nan(1,1);
        end
        temp_PAR(1)=sscanf(logtxt, 'Sigma factor (fitted for growth on glucose): %f');
        PAR_M=[PAR_M;temp_PAR];
        %o=o+1;
    end
    logtxt=fgetl(lf);
end
fclose(lf);
%compile table with condition names GAM and sigma
loginf=[cell2table(Mod_n', 'VariableNames', {'condNames'}), array2table(PAR_M, 'VariableNames', {'Sigma'})] ;
end