function[davidi_mod]=create_pFBAmod(raw_batchmod, enzMetPfx)
%function to create a batchmodell from pFBA data using davidi data for
%E.coli and Chen data for yeast
%INPUT:
% - struc raw_batchmod: A unmodified edModel obtained from GECKO
% - char enzMetPfx: Character vector marking enzyme pseudometabolites

  %read Davidi data 
if contains(raw_batchmod.name, 'Escherichia coli', 'IgnoreCase', true)
    kapp=readtable('Data/Davidi/davidi_suppl.xlsx', 'Sheet', 'kmax 1s');
    %adapt naming
    kapp.Properties.VariableNames([1 4])={'ReactionID', 'kappmax'};
    %get reaction mapping information
    map=readtable('Data/Davidi/davidi_mapping_manual.xlsx');
    %filter unmatched
    map=map(~cellfun(@isempty, map.final_match),[1,4]);
    %duplicate isozyme entries. 
elseif contains(raw_batchmod.name, 'Yeast', 'IgnoreCase', true)
    kapp=readtable('Data/Chen/chen_kmax.xlsx', 'Sheet', 'kmax 1s');
    %adapting naming 
    kapp.Properties.VariableNames([1 6])={'ReactionID', 'kappmax'};
    %get reaction mapping information
    map=readtable('Data/Chen/chen_mapping.xlsx');
    %filter unmatched
    map=map(~cellfun(@isempty, map.final_match),[1,4]);
    %duplicate isozyme entries. 
end

    i=1;
    while i<=size(map,1)
        if contains(map.final_match(i), ';')
            isozms=strsplit(map.final_match{i}, '; ')';
            tmp_map=map(repelem(i, length(isozms)), :);
            tmp_map.final_match=isozms;
            map=[map(1:(i-1), :); tmp_map; map((i+1):end,:)];
            i=i+length(isozms);
        else
            i=i+1;
        end
    end
     enzMetIdx = find(contains(raw_batchmod.mets, enzMetPfx) &~ismember(raw_batchmod.mets, {'prot_pool'}));

davidi_mod=raw_batchmod;
    for i=1:size(map,1)
%         %retrieve davidi kapp%                 disp('oldkcat:')
         new_kapp=-1/(kapp.kappmax(ismember(kapp.ReactionID, map.reactionName(i)))*3600);
        %retrieve max presto kcat 
        enz_idx=find(raw_batchmod.S(enzMetIdx,ismember(raw_batchmod.rxns, map.final_match(i))));
        if length(enz_idx)~=1
            if length(unique(raw_batchmod.S(enzMetIdx(enz_idx),ismember(raw_batchmod.rxns, map.final_match(i)))))==1
                %complexes should not be found but accept them if kcat is
                %identical

                davidi_mod.S(enzMetIdx(enz_idx), ismember(raw_batchmod.rxns, map.final_match(i)))=new_kapp;
            else
                warning(['reaction with missing enzyme or enzyme complex detected comparison impossible. skipping reaction ', ...
                    raw_batchmod.rxns{ismember(raw_batchmod.rxns, map.final_match(i))}])
            end
        else

            davidi_mod.S(enzMetIdx(enz_idx), ismember(raw_batchmod.rxns, map.final_match(i)))=new_kapp;
        end
    end
    end