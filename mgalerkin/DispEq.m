clear all

%=================================Constants================================
e0 = 8.85418782 * 10^ - 12;
m0 = 12.5663706144 * 10^ - 7;
C_LIGHT = 1 / sqrt(e0 * m0);

e1 = 2.2;
e2 = 1;

n = 1;

M = 5;
K = 5;

BETA_BEGIN = 0.1;
BETA_STEP = 0.1;
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
wAns = [];
betaAns = [];

freq = FREQ_BEGIN;

w = 2 * pi * freq * FREQ_SCALE;

betaRoots = [];
wRoots = [];

while freq <= FREQ_END
    w = 2 * pi * freq * FREQ_SCALE;

    tempRoots = [];

    for beta = BETA_BEGIN:BETA_STEP:BETA_END
        fprintf('freq = %f;  beta = %f\n', freq, beta);

        betaRoot = NaN;

        try
            [betaRoot, fval] = fzero(@(betaX) CalcDet(betaX, w, n, e1, e2, radius, M, K, intTol, begOfIntReg), [beta, beta + BETA_STEP]);
            fprintf('\t betaRoot = %f; function value = %f\n', betaRoot, fval);

            if (abs(fval) > 0.001)
                continue;
            end

        catch
            fprintf('\t cant find a root on interval [%f, %f]\n', beta, beta + BETA_STEP);
            continue;
        end

        tempRoots = [tempRoots, betaRoot];

    end

    tempRoots(tempRoots <= 0) = NaN;
    tempRoots(tempRoots >= sqrt(e1)) = NaN;
    tempRoots(diff(tempRoots) < 10^ - 5) = NaN;
    tempRoots = tempRoots(~isnan(tempRoots));
    w_ = [];
    w_(1:size(tempRoots, 2)) = w;
    wRoots = [wRoots, w_];
    betaRoots = [betaRoots, tempRoots];

    freq = freq + FREQ_STEP;

end

plot(wRoots, betaRoots, '.');
