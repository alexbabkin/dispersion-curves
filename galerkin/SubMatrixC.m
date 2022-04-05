function retVal = SubMatrixB(alpha_nm_sh, alpha_nm, e1, e2, de, w, radius, n, BETA, A1, A2, intTol, begOfIntReg)

    e0 = 8.85418782 * 10^ - 12;

    coeff1 = -w * alpha_nm^2;
    coeff2 = -w * BETA^2 * de * alpha_nm;
    coeff3 = w;

    %er = e0 * ((e2 - e1) * r / radius + e1);
    unInt1 = @(r) r .* e0 .* ((e2 - e1) .* r ./ radius + e1) .* besselj(n, alpha_nm .* r) .* besselj(n, alpha_nm_sh .* r) ./ (A1 .* r + A2);
    unInt2 = @(r) r .* dbesselj(n, alpha_nm .* r) .* besselj(n, alpha_nm_sh .* r) ./ (A1 .* r + A2).^2;
    unInt3 = @(r) r .* e0 .* ((e2 - e1) .* r ./ radius + e1) .* besselj(n, alpha_nm .* r) .* besselj(n, alpha_nm_sh .* r);

    retVal = coeff1 * integral(unInt1, begOfIntReg, radius, 'AbsTol', intTol) + ...
        coeff2 * integral(unInt2, begOfIntReg, radius, 'AbsTol', intTol) + ...
        coeff3 * integral(unInt3, begOfIntReg, radius, 'AbsTol', intTol);
