function retVal = SubMatrixB(kapa_nk_sh, kapa_nk, w, radius, n, A1, A2, de, intTol, begOfIntReg)

    m0 = 12.5663706144 * 10^ - 7;

    coeff1 = -w * m0 * kapa_nk^2;
    coeff2 = -w^3 * m0^2 * de * kapa_nk;
    coeff3 = m0 * w;

    unInt1 = @(r) r .* besselj(n, kapa_nk .* r) .* besselj(n, kapa_nk_sh .* r) ./ (A1 .* r + A2);
    unInt2 = @(r) r .* dbesselj(n, kapa_nk .* r) .* besselj(n, kapa_nk_sh .* r) ./ (A1 .* r + A2).^2;
    unInt3 = @(r) r .* besselj(n, kapa_nk .* r) .* besselj(n, kapa_nk_sh .* r);

    retVal = coeff1 * integral(unInt1, begOfIntReg, radius, 'AbsTol', intTol) + ...
        coeff2 * integral(unInt2, begOfIntReg, radius, 'AbsTol', intTol) + ...
        coeff3 * integral(unInt3, begOfIntReg, radius, 'AbsTol', intTol);
