function model  = readGKOmodel(file)
% Function to read a mat file containing a ecMOdel generated with GECKO
% INPUT:
%  - char file:    The filepath to the .mat file 
% OUTPUT:
%     - struc model: A ecModel in GECKO/RAVEN fromat
model=load(file);
fn=fieldnames(model);
if length(fn)==1
    model=model.(fn{1});
else
    error('.mat file contains more than 1 variable')
end
end
