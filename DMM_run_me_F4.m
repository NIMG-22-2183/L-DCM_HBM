% This code was used as part of the L-DCM developments. 

clc
clear all
close all
spm('defaults','EEG');
D_F = '/imaging/...';
%--------------------Addres to my Imaging space folders--------------------
 aj08_root       =  '/imaging/...';
%--------------------------------------------------------------------------
AF_sub = {'subject list'}; % Wite your subject'names here

subjects_name = AF_sub; 
GCM = {} ;

 parfor i = 1:length(subjects_name)
    cd(D_F)              % change dir
    counter              = subjects_name(i);
    subject              = subjects_name{1,i}   ;
    meg_data             = spm_select('Fplist',fullfile(D_F,subject), '^.*\open_DMN.mat$');
    sub                  = append('DCM','_', subject);
    D1{i,1}              = DCM_CSD_F4(sub,meg_data,i);  
end
 
save F4_test2  D1
%%
