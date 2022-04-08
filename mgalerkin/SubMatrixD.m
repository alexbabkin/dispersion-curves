function retVal = SubMatrixD(BETA, w, n, de, radius, kapa_nk, alpha_nm, A1, A2, intTol, begOfIntReg)
    m0 = 12.5663706144 * 10^ - 7;

    coeff = -BETA * m0 * w^2 * n * de;

    unInt = @(r) (besselj(n, alpha_nm .* r) .* besselj(n, kapa_nk .* r)) ./ (A1 .* r + A2).^2;

    retVal = coeff * integral(unInt, begOfIntReg, radius, 'AbsTol', intTol, 'ArrayValued', true);
