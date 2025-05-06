function DCM   = DCM_CSD_F4(sub,meg_data,i)
% A demo of DCM for CSD. This code can be used for one channel LFP data
% that are extracted from MEG. 
%==========================================================================
 % not always a good practice to close all ! 

% Add your SPM to Matlab path
%==========================================================================

% Spesify (un-estiamted) DCM
%==========================================================================

DCM.options.analysis = 'CSD';       % DCM for cross-spectral density 
DCM.options.spatial  = 'LFP';       % virtual electrode data as data 
DCM.options.model    = 'CMM_NMDA_AJ';  % generative model is  spm_fx_cmc
DCM.options.trials   = [1 2];        % index of BL vs AF within file// we donot need that we have trail=1:Dconditon
tim_ax(1)            = 1;           % first sample of MEG data (time domain)
tim_ax(2)            = 501;        % last sample of MEG data (time domain)
DCM.options.Tdcm     = [tim_ax(1) tim_ax(2)];  % time axis in milli seconds

DCM.options.Fdcm     = [1   64];    % frequency range  
DCM.M.Hz             = (1 : 64)';
DCM.options.D        = 1;            % downsampling
DCM.Sname            = {'MPFC', 'PRC', 'LAG','RAG'}; % name(s) of the source

% Adress to data
%==========================================================================
 
DCM.xY.Dfile         = meg_data; 

% if model one LFP channel 
%==========================================================================

DCM.xY.modality      = 'LFP';
DCM.xY.Ic            =  1:4;
DCM.xY.name          = {'MPFC', 'PRC', 'LAG','RAG'}; 

% A, B ,C parameters spesficiation in the generative model  spm_fx_cmc
% for one region source  A =[0]. For one condtion B = [0] 
%==========================================================================
Ns                  = 4; 
MPFC                = 1 ; 
PRC                 = 2; 
LAG                 = 3; 
RAG                 = 4;  
%-------Example------------------------------------------------------
DCM.A{1}             = ones(Ns,Ns);
DCM.A{1}             = DCM.A{1} - diag(diag(DCM.A{1}));

DCM.A{2}             = ones(Ns,Ns);
DCM.A{2}             = DCM.A{2} - diag(diag(DCM.A{2}));

% DCM.A{3}             = zeros(Ns,Ns);

% DCM.A{3}                      = [1   0   0   0
%                                 0   1    0   0
%                                 0   0    1   0
%                                 0   0    0   1] ;
DCM.G                    = ones(Ns,1);
DCM.B                    = zeros(Ns,Ns); % For between condtion case this needs to be defined

DCM.C                   = sparse(length(DCM.A{1}), 0);
DCM.options.D           = 1;
DCM.M.dipfit.Nm         = 1;
DCM.M.dipfit.model      = DCM.options.model;
DCM.M.dipfit.type       = DCM.options.spatial;
DCM.M.dipfit.Nc         = Ns;
DCM.M.dipfit.Ns         = Ns;
DCM.M.U                 = eye(Ns,Ns); 
DCM.name                = append(sub, '_F4');     %'sub_P1002';  % sprintf('DCM_%s_%s',sub{1,1},BT);

DCM.xU.X = [0; 1]; % how from one trial to next we modualte by the B matrix
    DCM.xU.name = {'F4'};
  
    
    %==========================================================================
    
[pE,pC]     = spm_dcm_neural_priors_AJAJ(DCM.A,DCM.B,DCM.C,DCM.options.model);
[pE,pC]     = spm_L_priors_AJAJ(DCM.M.dipfit,pE,pC); % 2024
[pE,pC]     = spm_ssr_priors(pE,pC);
DCM.M.pE    = pE;
DCM.M.pC    = pC;
    % spm_dcm_csd_data needs to be called for one channel data within many lfps.
    %==========================================================================
    
    DCM.options.DATA     = 0;
    DCM             = spm_dcm_csd_data(DCM);
% Model inversion/estiamtion
%==========================================================================

 DCM                  = spm_dcm_csd_AJAJ(DCM);

% return
% %**************************************************************************
% % Results
% %==========================================================================
% 
% %  spm_dcm_csd_results(DCM,'spectral data');%na
% %  spm_dcm_csd_results(DCM,'Coupling (A)'); % not useful for 1 region
% %  spm_dcm_csd_results(DCM,'Coupling (B)'); % not useful for 1 condtion
% %  spm_dcm_csd_results1(DCM,'Coupling (C)'); %na
% %  spm_dcm_csd_results(DCM,'trial-specific effects');%na
% %  spm_dcm_csd_results(DCM,'Input'); % Spectrum of 1/f b and c ch noise 
%   spm_dcm_csd_results(DCM,'Transfer functions'); % Transfer function of CMC
% %  spm_dcm_csd_results(DCM,'Cross-spectra (sources)'); 
% %  spm_dcm_csd_results(DCM,'Cross-spectra (channels)'); 
% %  spm_dcm_csd_results(DCM,'Coherence (sources)')% na
% %  spm_dcm_csd_results(DCM,'Coherence (channels)') %na
% %  spm_dcm_csd_results(DCM,'Covariance (sources)');
% %  spm_dcm_csd_results(DCM,'Covariance (channels)');
% % FORMAT spm_dcm_csd_results(DCM,'Dipoles');  % na only for ECD/IMG
% 
% %Population specific cross spectra
% %==========================================================================
% 
% M             = rmfield(DCM.M,'U'); 
% M.dipfit.type = 'LFP';
% M             = DCM.M;
% M.U           = 1; 
% M.l           = DCM.M.m;
% qp            = DCM.Ep;
% qp.L          = ones(1,M.l);              % set electrode gain to unity
% qp.b          = qp.b - 32;                % and suppress non-specific and
% qp.c          = qp.c - 32;                % specific channel noise
% 
% % response of j-th population in the i-th source
% %--------------------------------------------------------------------------
% 
% i           = 1; % we only have one source 
% j           = 3; % please see spm_fx_cmc for index
% % qp.J{i}     = spm_zeros(qp.J{i}); % if you have many sources use {}
% % qp.J{i}(j)  = 1; % if you have many sources use {}
% 
% qp.J        = zeros(size(qp.J));
% qp.J(j)     = 1; % virtual electrod unity gain (selector of a population)
% 
% [Hs ,Hz, dtf] = spm_csd_mtf(qp,M,DCM.xU); % conditional cross spectra
% [ccf]     = spm_csd2ccf(Hs,DCM.M.Hz); % conditional correlation functions
% [coh]     = spm_csd2coh(Hs,DCM.M.Hz); % conditional covariance
% figure
% plot (Hz, abs(dtf{1, 1}));
% box off 
% axis tight
% end

