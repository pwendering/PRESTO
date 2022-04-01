function get_rawecmod()
%Script to obtain the different EC models from GECKO branch
hdir=pwd();
%Load GEM
if ~isfolder('Logs')
    mkdir Logs
end

diary(['Logs/' regexprep(date, '-', '') '_getrawecmods.log'])
%Load GEM
yeastmod=readGKOmodel('Data/yeast-GEM.mat');
ecolimod=readGKOmodel('Data/ecoli-GEM.mat');
cd('GECKO_S_cerevisiae/geckomat')

enhanceGEM(yeastmod, 'COBRA', true, 'ecYeast', '8')

cd('../../GECKO_E_coli/geckomat')

enhanceGEM(ecolimod, 'COBRA', true, 'ecEcoli', 'iML1515')






%copyfile('ecYeast/*.mat' , [topDir '/Data/'])
cd('../../')
diary off
end


