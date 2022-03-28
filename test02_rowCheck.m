function [] = test02_rowCheck(w, ySol, Xr, Yb, n, m1)
%% Zestaw Testowy
% w = @(x)[2, 1, 1, -x];
% ySol = @(x)1 / 28 * (14 * x + 9 * sqrt(7) * exp(-x/2) .* sin(sqrt(7)*x/2) + 35 * exp(-x/2) .* cos(sqrt(7)*x/2) - 7);
% Xr = [0, 6];
% Yb = [1, 1];
% n = 10:1:50;
% m1 = 5;

%% Test
h = (Xr(2) - Xr(1)) ./ n;

errAM = zeros(1, length(n));
errH = zeros(1, length(n));
errHXAM = zeros(1, length(n));

for i = 1:length(n)
    x = linspace(Xr(1), Xr(2), n(i)+1);
    [YH, YAM] = HXAMMain(w, Xr, Yb, n(i), m1, 0);
    [~, YHXAM] = HXAMMain(w, Xr, Yb, n(i), m1, 1e-10);
    errH(1, i) = sum(abs(ySol(x)-YH(2, :))) / n(i);
    errAM(1, i) = sum(abs(ySol(x)-YAM(2, :))) / n(i);
    errHXAM(1, i) = sum(abs(ySol(x)-YHXAM(2, :))) / n(i);
end

figure(1)
semilogy(h, errH(1, 1)./errH(1, :), h, errAM(1, 1)./errAM(1, :), h, errHXAM(1, 1)./errHXAM(1, :))
xlabel('h')
ylabel('1/(err(h)/err(h0))')
title('wykres względnej zmiany dokładności')
legend('Heun', 'AM', 'AM+Broyden')
end