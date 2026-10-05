%% MaterialClass_example.m

close all; clc;

thisDir = fileparts(mfilename('fullpath'));
addpath(fileparts(thisDir));

nPass = 0;
nFail = 0;

function [np, nf] = runTest(name, np, nf, testFcn)
    try
        testFcn();
        fprintf('[PASS] %s\n', name);
        np = np + 1;
    catch ME
        fprintf('[FAIL] %s  —  %s\n', name, ME.message);
        nf = nf + 1;
    end
end

%% 1. CONSTRUCTOR
fprintf('\n=== 1. CONSTRUCTOR ===\n');

[nPass, nFail] = runTest('Empty constructor', nPass, nFail, ...
    @() assert(MaterialClass().d == 3));

geo.dimension = 3;
freq = 50;
lc = 0.1;

vAc = 2000;
covAc = [0.05, 0.05];
ccAc = 0;

mAc = MaterialClass(geo, freq, true, vAc, covAc, ccAc, 'exp', lc);
[nPass, nFail] = runTest('Acoustic constructor', nPass, nFail, ...
    @() assert(mAc.acoustics && mAc.v == vAc));

vp = 6000;
vs = vp / sqrt(3);
covEl = [0.05 0.05 0.05];
ccEl = [0 0 0];

mEl = MaterialClass(geo, freq, false, [vp vs], covEl, ccEl, 'exp', lc);
[nPass, nFail] = runTest('Elastic constructor', nPass, nFail, ...
    @() assert(~mEl.acoustics && mEl.vp == vp && mEl.vs == vs));

%% 2. COPYOBJ
fprintf('\n=== 2. COPYOBJ ===\n');

mCopy = mAc.copyobj();
mCopy.v = 9999;

[nPass, nFail] = runTest('copyobj original unchanged', nPass, nFail, ...
    @() assert(mAc.v == vAc));
[nPass, nFail] = runTest('copyobj modified', nPass, nFail, ...
    @() assert(mCopy.v == 9999));

%% 3. PRESET
fprintf('\n=== 3. PRESET ===\n');

[nPass, nFail] = runTest('preset acoustic', nPass, nFail, ...
    @() assert(MaterialClass.preset(1).acoustics));
[nPass, nFail] = runTest('preset elastic', nPass, nFail, ...
    @() assert(~MaterialClass.preset(3).acoustics));

%% 4. PSDF METHODS
fprintf('\n=== 4. PSDF METHODS ===\n');

mP = MaterialClass();
mP.d = 3;
mP.acoustics = true;

[nPass, nFail] = runTest('exp', nPass, nFail, ...
    @() testExp(mP, lc));
[nPass, nFail] = runTest('power_law', nPass, nFail, ...
    @() testPower(mP, lc));
[nPass, nFail] = runTest('gaussian', nPass, nFail, ...
    @() testGauss(mP, lc));
[nPass, nFail] = runTest('triangular', nPass, nFail, ...
    @() testTri(mP, lc));
[nPass, nFail] = runTest('low_pass', nPass, nFail, ...
    @() testLow(mP, lc));
[nPass, nFail] = runTest('VonKarman', nPass, nFail, ...
    @() testVK(mP, lc));

%% 5. getPSDF
fprintf('\n=== 5. getPSDF ===\n');

mG = MaterialClass(geo, freq, true, vAc, covAc, ccAc, 'exp', lc);
[nPass, nFail] = runTest('getPSDF', nPass, nFail, ...
    @() testGetPSDF(mG));

%% 6. CalcSigma
fprintf('\n=== 6. CalcSigma ===\n');

mAc2 = MaterialClass(geo, freq, true, vAc, covAc, ccAc, 'exp', lc);
mAc2.CalcSigma();

[nPass, nFail] = runTest('sigma acoustic handle', nPass, nFail, ...
    @() assert(isa(mAc2.sigma{1}, 'function_handle')));

mEl2 = MaterialClass(geo, freq, false, [vp vs], covEl, ccEl, 'exp', lc);
mEl2.CalcSigma();

[nPass, nFail] = runTest('sigma elastic handle', nPass, nFail, ...
    @() assert(isa(mEl2.sigma{1,1}, 'function_handle')));

%% 7. prepareSigmaOne
fprintf('\n=== 7. prepareSigmaOne ===\n');

[Sig, ~, invcdfFn] = MaterialClass.prepareSigmaOne(mAc2.sigma{1}, 3);

[nPass, nFail] = runTest('Sigma positive', nPass, nFail, ...
    @() assert(Sig > 0));
[nPass, nFail] = runTest('invcdf valid', nPass, nFail, ...
    @() assert(invcdfFn(0.5) >= 0 && invcdfFn(0.5) <= pi));

%% 8. prepareSigma
fprintf('\n=== 8. prepareSigma ===\n');

mAc3 = MaterialClass.prepareSigma(mAc2.copyobj(), 3);

[nPass, nFail] = runTest('Diffusivity positive', nPass, nFail, ...
    @() assert(mAc3.Diffusivity > 0));

mEl3 = MaterialClass.prepareSigma(mEl2.copyobj(), 3);
[nPass, nFail] = runTest('elastic transport mean free times', ...
    nPass, nFail, @() testElasticTransportMeanFreeTimes(mEl3));

%% 9. HOMOGENEOUS MEDIA
fprintf('\n=== 9. HOMOGENEOUS MEDIA ===\n');

zeroSigma = @(angle) zeros(size(angle));
[zeroTotal, zeroFirstMoment, zeroInverseCDF] = ...
    MaterialClass.prepareSigmaOne(zeroSigma,3);
[nPass, nFail] = runTest('zero scattering operator', nPass, nFail, ...
    @() testZeroScatteringOperator( ...
        zeroTotal,zeroFirstMoment,zeroInverseCDF));

homogeneousAcoustic = MaterialClass( ...
    geo,freq,true,vAc,[0 0],0,'exp',lc);
homogeneousAcoustic = MaterialClass.prepareSigma(homogeneousAcoustic,3);
[nPass, nFail] = runTest('homogeneous acoustic material', nPass, nFail, ...
    @() testHomogeneousAcoustic(homogeneousAcoustic));

homogeneousElastic = MaterialClass( ...
    geo,freq,false,[vp vs],[0 0 0],[0 0 0],'exp',lc);
homogeneousElastic = MaterialClass.prepareSigma(homogeneousElastic,3);
[nPass, nFail] = runTest('homogeneous elastic material', nPass, nFail, ...
    @() testHomogeneousElastic(homogeneousElastic));

% Exact zero mode-conversion channels are also valid when same-mode
% scattering remains active.
noConversionElastic = MaterialClass();
noConversionElastic.d = 3;
noConversionElastic.acoustics = false;
noConversionElastic.vp = vp;
noConversionElastic.vs = vs;
isotropicSigma = @(angle) ones(size(angle));
noConversionElastic.sigma = {isotropicSigma,zeroSigma; ...
                             zeroSigma,isotropicSigma};
noConversionElastic = MaterialClass.prepareSigma(noConversionElastic,3);
[nPass, nFail] = runTest('zero mode-conversion channels', nPass, nFail, ...
    @() testZeroModeConversion(noConversionElastic));

%% 10. CalcLc
fprintf('\n=== 10. CalcLc ===\n');

mLc = MaterialClass(); mLc.d = 3;
mLc.Exponential(lc);

[nPass, nFail] = runTest('Lc positive', nPass, nFail, ...
    @() assert(mLc.CalcLc() > 0));

%% SUMMARY
fprintf('\n============================\n');
fprintf('TOTAL %d | PASS %d | FAIL %d\n', nPass+nFail, nPass, nFail);
fprintf('============================\n');

%% ================= LOCAL TEST FUNCTIONS =================

function testExp(mP, lc)
    m = mP.copyobj();
    m.SpectralLaw = 'exp';
    m.Exponential(lc);
    assert(m.Phi(1) > 0);
end

function testZeroScatteringOperator(Sigma,SigmaPrime,inverseCDF)
    if Sigma ~= 0 || SigmaPrime ~= 0
        error('An exactly zero DSCS must have zero angular integrals.');
    end
    probabilities = [0 0.5 1];
    if any(inverseCDF(probabilities) ~= 0)
        error('The zero-DSCS inverse-CDF placeholder is invalid.');
    end
end

function testElasticTransportMeanFreeTimes(material)
    velocities = [material.vp; material.vs];
    expectedTimes = material.transportMeanFreePath ./ velocities;

    if ~isequal(size(material.transportMeanFreeTime),[2 1])
        error('Elastic transport mean free times must be a 2-by-1 vector.');
    end

    relativeError = norm(material.transportMeanFreeTime-expectedTimes) / ...
        norm(expectedTimes);
    if relativeError >= 1e-12
        error('The transport-time relative error is too large: %e', ...
            relativeError);
    end
end

function testHomogeneousAcoustic(material)
    if material.Sigma ~= 0 || material.Sigmapr ~= 0
        error('A homogeneous acoustic material must have Sigma = 0.');
    end
    if ~isinf(material.meanFreeTime) || ~isinf(material.meanFreePath)
        error('A homogeneous acoustic material must have infinite free scales.');
    end
    if ~isnan(material.Diffusivity) || ~isnan(material.g)
        error('Diffusivity and g must be undefined without scattering.');
    end
end

function testHomogeneousElastic(material)
    if any(material.Sigma ~= 0,'all') || any(material.Sigmapr ~= 0,'all')
        error('A homogeneous elastic material must have Sigma = 0.');
    end
    if any(~isinf(material.meanFreeTime)) || ...
            any(~isinf(material.meanFreePath))
        error('A homogeneous elastic material must have infinite free scales.');
    end
    if material.P2P ~= 1 || material.S2S ~= 1
        error('Safe same-mode placeholders must be used without scattering.');
    end
    if ~isnan(material.Diffusivity)
        error('Diffusivity must be undefined without scattering.');
    end
end

function testZeroModeConversion(material)
    if material.Sigma(1,2) ~= 0 || material.Sigma(2,1) ~= 0
        error('The zero mode-conversion rates were not preserved.');
    end
    if material.P2P ~= 1 || material.S2S ~= 1
        error('Particles must retain their modes when conversion rates vanish.');
    end
    if any(~isfinite(material.meanFreeTime))
        error('Same-mode scattering must retain finite mean free times.');
    end
end

function testPower(mP, lc)
    m = mP.copyobj();
    m.SpectralLaw = 'power_law';
    m.PowerLaw(lc);
    assert(m.Phi(1) > 0);
end

function testGauss(mP, lc)
    m = mP.copyobj();
    m.SpectralLaw = 'gaussian';
    m.Gaussian(lc);
    assert(m.Phi(1) > 0);
end

function testTri(mP, lc)
    m = mP.copyobj();
    m.SpectralLaw = 'triangular';
    m.Triangular(lc);
    assert(m.Phi(0.5) >= 0);
end

function testLow(mP, lc)
    m = mP.copyobj();
    m.SpectralLaw = 'low_pass';
    m.LowPass(lc);
    assert(m.Phi(0.5) >= 0);
end

function testVK(mP, lc)
    m = mP.copyobj();
    m.SpectralLaw = 'VonKarman';
    m.VonKarman(lc, 0.5);
    assert(m.Phi(1) > 0);
end

function testGetPSDF(mG)
    mG.getPSDF();
    assert(~isempty(mG.Phi));
end
