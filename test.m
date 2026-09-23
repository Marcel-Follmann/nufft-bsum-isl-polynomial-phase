
d = [0, 3, -3, 4, -4, 6, -6, 12, -12];
g = calculate_gcd(d);

function g = calculate_gcd(d)
    g = d(2);
    for i = 4:2:numel(d)
        g = gcd(g, d(i));
    end
end