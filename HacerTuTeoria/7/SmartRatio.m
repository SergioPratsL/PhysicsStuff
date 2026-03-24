function val = SmartRatio(n, d)


if (norm(n) < 10^-14)
    val = 0;
    return
end

val = n / d;