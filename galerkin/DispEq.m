clear all

%=================================Constants================================
e0 = 8.85418782 * 10^ - 12;
m0 = 12.5663706144 * 10^ - 7;
C_LIGHT = 1 / sqrt(e0 * m0);
%==================================UserDefs================================
e1 = 2.2;
e2 = 1;

n = 1;

M = 5;
K = 5;

BETA_BEGIN = 0.01;
BETA_STEP = 0.01;
BETA_END = sqrt(e1) + BETA_STEP;

Norm = 10e0;
intTol = 1.0e-8;
begOfIntReg = 0;

FREQ_BEGIN = 0.1;
FREQ_STEP = 0.1;
FREQ_END = 50;
FREQ_SCALE = 10^9;

radius = 0.02;
%==========================================================================

Index = 1;
Ind = 1;

freq = FREQ_BEGIN;

de = e0 * (e2 - e1) / radius;

while freq <= FREQ_END
    w = 2 * pi * freq * FREQ_SCALE;
    k0 = w / C_LIGHT;

    for beta = BETA_BEGIN:BETA_STEP:BETA_END
        fprintf('freq = %f;  beta = %f\n', freq, beta);

        matrixDets(Index) = CalcDet(beta, w, n, e1, e2, radius, M, K, intTol, begOfIntReg);
        betas(Index) = beta;

        Index = Index + 1;
    end

    ansIndexes = zerosIndexes(matrixDets);
    ansBeta = betas(ansIndexes);

    for i = 1:1:size(ansBeta, 2)
        wAns(Ind) = freq * FREQ_SCALE;
        betaAns(Ind) = ansBeta(i) - BETA_STEP / 2;
        Ind = Ind + 1;
    end

    matrixDets = [];
    betas = [];
    Index = 1;

    freq = freq + FREQ_STEP;

end

plot(wAns, betaAns, '.');
