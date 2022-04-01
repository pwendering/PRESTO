%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% model = saveECmodel(model,toolbox, root_name,name,version)
% If toolbox=='COBRA' a the unconverted RAVEN model with csense field is saved as .mat file
% INPUT:
% - struc model:    model structure generated from readKcatData.m or
%                           getConstrainedModel.m
% - char toolbox:   name of the toolbox the model should be optimized for
%                           'COBRA' or 'RAVEN'
% - char root_name:     path of the directory (eg. ('ecYeast')
% - char name:  name of the model file (eg. 'Yeast_batch')
% - char version: version of the model
% Benjamin J. Sanchez. Last edited: 2018-10-25
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = saveECmodel(model,toolbox, root_name, name,version)

fprintf(['Saving ' name version ':\n'])

%Define file path for storage:
struct_name = 'ecModel';
if endsWith(name,'_batch')
    struct_name = [struct_name '_batch'];
end

if ~isdir(root_name)
    mkdir(root_name)
end
file_name = [root_name '/' name];

%Model description:
model.description = [struct_name ' of ' lower(name(3)) name(4:end)];
model.id          = [name '_v' version];

%Format S matrix: avoid long decimals
for i = 1:length(model.mets)
    for j = 1:length(model.rxns)
        if model.S(i,j) ~= 0
            orderMagn    = ceil(log10(abs(model.S(i,j))));
            model.S(i,j) = round(model.S(i,j),6-orderMagn);
        end
    end
end

%For functional models, save upper bounds as +1000:
model.ub(isinf(model.ub)) = 1000;

%Remove model.rules (added by COBRA functions)
model = takeOutField(model,'rules');

if strcmp(toolbox,'COBRA')
    %add csense field
    model.csense=repmat('E', length(model.mets), 1);
    %Save GECKO version as mat file
    save([file_name '.mat'], 'model')
    
    %Save the unconverted model as mat file
    
    %Transform model back to COBRA for saving purposes:
    model_cobra = ravenCobraWrapper(model);
    %Remove fields from COBRA model (temporal):
    model_cobra = takeOutField(model_cobra,'metCharges');
    model_cobra = takeOutField(model_cobra,'metChEBIID');
    model_cobra = takeOutField(model_cobra,'metKEGGID');
    model_cobra = takeOutField(model_cobra,'metNotes');
    model_cobra = takeOutField(model_cobra,'metSBOTerms');
    model_cobra = takeOutField(model_cobra,'rxnConfidenceScores');
    model_cobra = takeOutField(model_cobra,'rxnECNumbers');
    model_cobra = takeOutField(model_cobra,'rxnKEGGID');
    model_cobra = takeOutField(model_cobra,'rxnReferences');
    model_cobra.subSystems = cell(size(model_cobra.rxns));
    model_cobra = takeOutField(model_cobra,'rxnSBOTerms');
    %Save model as sbml and text:
    writeCbModel(model_cobra,'sbml',[file_name '.xml']);
    writeCbModel(model_cobra,'text',[file_name '.txt']);
else
    exportForGit(model,name,root_name,{'xml','yml','txt','mat'});
end

%Convert notation "e-005" to "e-05 " in stoich. coeffs. to avoid
%inconsistencies between Windows and MAC:
copyfile([file_name '.xml'],'backup.xml')
fin  = fopen('backup.xml', 'r');
fout = fopen([file_name '.xml'], 'w');
still_reading = true;
while still_reading
    inline = fgets(fin);
    if ~ischar(inline)
        still_reading = false;
    else
        if ~isempty(regexp(inline,'-00[0-9]','once'))
            inline = strrep(inline,'-00','-0');
        elseif ~isempty(regexp(inline,'-01[0-9]','once'))
            inline = strrep(inline,'-01','-1');
        end
        fwrite(fout, inline);
    end
end
fclose('all');
delete('backup.xml');

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = takeOutField(model,field)

if isfield(model,field)
    model = rmfield(model,field);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
