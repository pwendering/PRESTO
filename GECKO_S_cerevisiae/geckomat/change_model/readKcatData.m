%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [eModel, out_tab] = readKcatData(model_data,kcats)
% Reads the output of the getBRENDAdata module, and creates the modified
% enzyme model, with additional metabolites (the enzymes) and reactions
% (for all isoenzymes and also the enzyme exchange reactions).
%
% INPUT:
% model_data        model and EC numbers and substrates/products from each
%                   reaction (output from "getECnumbers.m")
% kcats             kcats for each reaction/enzyme (output from
%                   "matchKcats.m")
%
% OUTPUT:
% eModel            modified model accounting for enzymes
% out_tab:    A table of size (forward+revers reaxtions)x20 that gives the
%              origin for the kcat linked to each reaction. 

% Cheng Zhang               2015-12-03
% Benjamin J. Sanchez       2018-08-11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [eModel, out_tab] = readKcatData(model_data,kcats)

%create a table that tracks kcat origins
 n_r=size(kcats.forw.org_s, 1);
 out_tab=zeros(n_r+sum(model_data.model.rev), size(kcats.forw.org_s,2));
 for i=1:size(out_tab, 1)
     for j=1:size(out_tab, 2)
         
         %index of reverse reactions 
         rev_idx=find(model_data.model.rev);
         %for forward reactions
         if i<=n_r
             orig=find([kcats.forw.org_s(i,j)  kcats.forw.rest_s(i,j)  kcats.forw.org_ns(i,j) ...
                 kcats.forw.org_sa(i,j) kcats.forw.rest_ns(i,j) kcats.forw.rest_sa(i,j)], 1, 'first');
             if ~isempty(orig)
             out_tab(i,j)=orig;
             end
         else %for reverse reactions
             orig=find([kcats.back.org_s(rev_idx(i-n_r),j)  kcats.back.rest_s(rev_idx(i-n_r),j) ...
                 kcats.back.org_ns(rev_idx(i-n_r),j) kcats.back.org_sa(rev_idx(i-n_r),j) kcats.back.rest_ns(rev_idx(i-n_r),j) ...
                 kcats.back.rest_sa(rev_idx(i-n_r),j)], 1, 'first');
             if ~isempty(orig)
                 out_tab(i,j)=orig;
             end
         end
     end
 end

%Get kcat value for both directions:
Fkcat = kcats.forw.kcats;
Bkcat = kcats.back.kcats;
rev   = logical(model_data.model.rev);
kcats = [Fkcat;Bkcat(rev,:)];



%Update uniprots with both directions:
uniprots = [model_data.uniprots; model_data.uniprots(rev,:)];

%Update matched genes with both directions:
matchedGenes = [model_data.matchedGenes; model_data.matchedGenes(rev,:)];

%Convert to irreversible model with RAVEN function (will split in 2 any reversible rxn):
model = convertToIrrev(model_data.model);
out_tab = [table(model.rxns), array2table(out_tab)];

%Convert original model to enzyme model according to uniprots and kcats:
eModel = convertToEnzymeModel(model,matchedGenes,uniprots,kcats);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
