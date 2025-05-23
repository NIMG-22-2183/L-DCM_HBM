function [x,M] = spm_x_cmm_NMDA_AJAJ(P)
% initialises a state structure for a mean field model
% FORMAT [x,M] = spm_x_cmm(P)
%
% P - parameter structure
% M - model structure
%
% x        - array of states
% x(i,j,k) - k-th state of j-th population on i-th source
%
%   population: 1 - excitatory spiny stellate cells (input cells)
%               2 - superficial pyramidal cells     (forward output cells)
%               3 - inhibitory interneurons         (intrisic interneuons)
%               4 - deep pyramidal cells            (backward output cells)
%
%        state: 1 V  - voltage
%               2 gE - conductance (excitatory)
%               3 gI - conductance (inhibitory)
%               4 gN - conductance (slow, voltage dependent excitatory)
%
% M - model structure
%
% see also: spm_x_mfm
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_x_cmm_NMDA.m 5741 2013-11-13 12:10:48Z guillaume $




% dimensions
%--------------------------------------------------------------------------
ns   = size(P.A{1},1);                           % number of sources
np   = 4;                                        % number of populations
 
% create (initialise voltage at -50mV)
%--------------------------------------------------------------------------
x        = zeros(ns,np,4) + 1/8;
x(:,:,1) = -50;

% steady-state solution 
%==========================================================================

% create model
%--------------------------------------------------------------------------
M.g   = {};
M.f   = 'spm_fx_cmm_NMDA_AJAJ';
M.x   = x;
M.pE  = P;
M.n   = length(spm_vec(x));
M.u   = sparse(ns,1);

% solve for steady state
%--------------------------------------------------------------------------
x     = spm_dcm_neural_x(P,M);



