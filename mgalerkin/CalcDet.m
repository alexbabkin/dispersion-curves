function retVal = CalcDet(beta, w, n, e1, e2, radius, M, K, intTol, begOfIntReg)
    e0 = 8.85418782 * 10^ - 12;
    m0 = 12.5663706144 * 10^ - 7;
    C_LIGHT = 1 / sqrt(e0 * m0);

    de = e0 * (e2 - e1) / radius;

    k0 = w / C_LIGHT;
    BETA = beta * k0;

    A1 = e0 * (e2 - e1) * m0 * w^2 / radius;
    A2 = e1 * e0 * m0 * w^2 - BETA^2;

    %================================a(k,m)====================================

    kapa_nk = dJ1_RootDivRad(K, radius);
    alpha_nm = J1_RootDivRad(M, radius);

    subMatrixA = SubMatrixA(BETA, w, n, de, radius, kapa_nk', alpha_nm, A1, A2, intTol, begOfIntReg);

    %================================d(k,m)====================================

    alpha_nm = J1_RootDivRad(M, radius);
    kapa_nk = dJ1_RootDivRad(K, radius);

    subMatrixD = SubMatrixD(BETA, w, n, de, radius, kapa_nk, alpha_nm', A1, A2, intTol, begOfIntReg);

    %================================b(k,m)====================================
    kapa_nk_sh = dJ1_RootDivRad(K, radius);
    kapa_nk = dJ1_RootDivRad(K, radius);

    subMatrixB = SubMatrixB(kapa_nk_sh', kapa_nk, w, radius, n, A1, A2, de, intTol, begOfIntReg);

    %================================c(k,m)====================================
    alpha_nm_sh = J1_RootDivRad(M, radius);
    alpha_nm = J1_RootDivRad(M, radius);

    subMatrixC = SubMatrixC(alpha_nm_sh', alpha_nm, e1, e2, de, w, radius, n, BETA, A1, A2, intTol, begOfIntReg);

    %==========================================================================
    MATRIX = [subMatrixA subMatrixB
        subMatrixC subMatrixD];

    retVal = det(MATRIX);
