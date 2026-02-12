% This function launches comparision cases for RadiativeTransferMonteCarlo
%
% when type='all' or empty, all validation cases are launched
%
% list of possible comparison cases with reference to literature
%  - '2dIsotropicAcoustic' (1,2,5)
%  - '3dIsotropicAcoustic' (1,2,5)
%  - '2dIsotropicElastic'  (3,5)
%  - '3dIsotropicElastic'  (4,5)
%  - '2dIsotropicAcoustic'
%  - '3dIsotropicAcoustic'
%  - '3dIsotropicElastic'
% Our code is compared to the results described in the following papers:
% (1) J. C. J. Paasschens. Solution of the time-dependent Boltzmann equation,
%     Phys. Rev. E 56(1), pp. 1135-1141 (1997).
% (2) M. Hoshiba. Simulation of multiple-scattered coda wave excitation
%     based on the energy conservation law. Phys. Earth Planet. Int. 67,
%     pp. 123-136 (1991).
% (3) H. Nakahara, K. Yoshimoto. Radiative transfer of elastic waves in
%     two-dimensional isotropic scattering media: semi-analytical approach
%     for isotropic source radiation. Earth Planets Space 63, pp. 459-468
%     (2011).
% (4) H. Sato. Multiple isotropic scattering model including P-S conversions
%     for the seismogram envelope formation. Geophys. J. Int 117,
%     pp. 487-494 (1994).
% (5) K. Yoshimoto, Monte Carlo simulation of seismogram envelopes in
%     scattering media. Journal of Geophysical Research: Solid Earth (2000)
% 
% 
clc
close all
clear all
%% 
% do not run ''by parts''
% this script is a a kind of non regression tests

% procedure:
% 1) You should run this script with F5
% 2) For each new method created you should add it in the example of the
% class
% 3) For a new class a new example script should be generated and inserted
% in the TestClass Method
%%
cd(fileparts(mfilename('fullpath')))

tags = {'Elastic','Acoustic','2D','3D','Material','All'};
tic
results = runtests({'NonRegTestClass'},'Tag',tags{end});
table(results)
close all
time = toc;
% check if there are some error in the execution
save('RunAllResults.mat')
%% To read the diagnosys from the results that went wrong
clc 
I_err = find([results.Failed]');
I_inc = find([results.Incomplete]');

if ~isempty(I_err)
    index = 1;
    results(I_err(index)).Details.DiagnosticRecord.Report
    results(I_err(index))
end
%% To Re RUN only the failed tests saved in the file 'RunAllResults.mat'
if 0
    %%
    clear all
    close all
    clc
    
    load('RunAllResults.mat');
    I_err = find([results.Failed]');
    for i_redo = 1 : numel(I_err)
        results(I_err(i_redo)) = runtests(results(I_err(i_redo)).Name);
    end
    save('RunAllResults.mat');
    I_err = find([results.Failed]');
    if ~isempty(results(I_err))
        table(results(I_err))
    else
        table(results)
    end
end