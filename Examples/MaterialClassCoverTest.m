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

%% 9. CalcLc
fprintf('\n=== 9. CalcLc ===\n');

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