function h = hann(n)
% Calculate a hann window (raised cosine)

if ~rem(n,2)
    % Even length window
    half = n/2;
    h = calc_window(half,n);
    h = [h; h(end:-1:1)];
else
    % Odd length window
    half = (n+1)/2;
    h = calc_window(half,n);
    h = [h; h(end-1:-1:1)];
end

return

function w = calc_window(m,n)
    x = (0:m-1)'/(n-1);
    w = 0.5 - 0.5*cos(2*pi*x);   
return

