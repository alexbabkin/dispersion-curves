function RetVar = dbesselj(n, arg)
    RetVar = (n .* besselj(n, arg)) ./ (arg) - besselj(n + 1, arg);
