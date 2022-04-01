function[model] = preprocessYeast(model)
% Function to remove wrong GPR annotations from Yeast GEM

if ~contains(model.modelName, 'Yeast', 'IgnoreCase', true)
    warning('Could not detect Yeast in model name, make sure changes are applied to a yeast GEM')
end

%Remove r_0227 assignment of P40009/YER005W
gID=find(strcmp(model.genes, 'YER005W'));
patt=['\| x\(' num2str(gID) '\) '];
model.rules(findRxnIDs(model, 'r_0227'))= regexprep(model.rules(findRxnIDs(model, 'r_0227')), patt, '');
if contains(model.rules(findRxnIDs(model, 'r_0227')), num2str(gID))
    error('Failed to remove YER005W from r_0227. check rules')
end
    
%Remove r_0322', 'r_0990' assignment of YKL060C
gID=find(strcmp(model.genes, 'YKL060C'));
patt=['x\(' num2str(gID) '\)'];
model.rules(findRxnIDs(model, {'r_0322', 'r_0990'}))=regexprep(model.rules(findRxnIDs(model, {'r_0322', 'r_0990'})), patt, '');
if any(contains(model.rules(findRxnIDs(model, {'r_0322', 'r_0990'})), num2str(gID)))
    error('Failed to remove YKL060C from r_0322 and r_0990. check rules')
end
end

