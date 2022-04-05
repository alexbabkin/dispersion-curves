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
BETA_END = sqrt(e1) - BETA_STEP;

Norm = 10e0;
intTol = 1.0e-8;
begOfIntReg = 0;

FREQ_BEGIN = 5;
FREQ_STEP = 1;
FREQ_END = 50;

radius = 0.02;
%==========================================================================

Index = 1;
Ind = 1;

freq = FREQ_BEGIN;

de = e0 * (e2 - e1) / radius;

while freq <= FREQ_END
    w = 2 * pi * freq * 10^9;
    k0 = w / C_LIGHT;

    for beta = BETA_BEGIN:BETA_STEP:BETA_END
        fprintf('freq = %f;  beta = %f\n', freq, beta);

        BETA = beta * k0;

        A1 = e0 * (e2 - e1) * m0 * w^2 / radius;
        A2 = e1 * e0 * m0 * w^2 - BETA^2;

        subMatrixA = zeros(M, K);
        subMatrixB = zeros(M, K);
        subMatrixC = zeros(M, K);
        subMatrixD = zeros(M, K);

        %================================a(k,m)====================================
        for k = 1:1:K
            kapa_nk = dJ1_RootDivRad(k, radius);

            for m = 1:1:M
                alpha_nm = J1_RootDivRad(m, radius);

                subMatrixA(k, m) = SubMatrixA(BETA, w, n, de, radius, kapa_nk, alpha_nm, A1, A2, intTol, begOfIntReg);

            end

        end

        %================================d(k,m)====================================
        for m = 1:1:M
            alpha_nm = J1_RootDivRad(m, radius);

            for k = 1:1:K
                kapa_nk = dJ1_RootDivRad(k, radius);

                subMatrixD(m, k) = SubMatrixD(BETA, w, n, de, radius, kapa_nk, alpha_nm, A1, A2, intTol, begOfIntReg);

            end

        end

        %================================b(k,m)====================================
        for k_sh = 1:1:K
            kapa_nk_sh = dJ1_RootDivRad(k_sh, radius);

            for k = 1:1:K
                kapa_nk = dJ1_RootDivRad(k, radius);

                subMatrixB(k_sh, k) = SubMatrixB(kapa_nk_sh, kapa_nk, w, radius, n, A1, A2, de, intTol, begOfIntReg);

            end

        end

        %================================c(k,m)====================================
        for m_sh = 1:1:M
            alpha_nm_sh = J1_RootDivRad(m_sh, radius);

            for m = 1:1:M
                alpha_nm = J1_RootDivRad(m, radius);

                subMatrixC(m_sh, m) = SubMatrixC(alpha_nm_sh, alpha_nm, e1, e2, de, w, radius, n, BETA, A1, A2, intTol, begOfIntReg);

            end

        end

        %==========================================================================
        MATRIX = [subMatrixA subMatrixB
            subMatrixC subMatrixD];

        matrixDets(Index) = det(MATRIX);
        betas(Index) = beta;

        Index = Index + 1;
    end

    ansIndexes = zerosIndexes(matrixDets);
    ansBeta = betas(ansIndexes);

    for i = 1:1:size(ansBeta, 2)
        wAns(Ind) = freq * 10^9;
        betaAns(Ind) = ansBeta(i) - BETA_STEP / 2;
        Ind = Ind + 1;
    end

    d = diff(matrixDets);
    %matrixDets = [];
    betas = [];
    Index = 1;
    freq = freq + FREQ_STEP;
end

plot(wAns, betaAns, '.');
