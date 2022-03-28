function [] = test01_workTest(w, ySol, Xr, Yb, n)
%% Zestaw testowy y'' + y' = -2*x^2
% w = @(x)[2, 1, 1, -x];
% ySol = @(x)1 / 28 * (14 * x + 9 * sqrt(7) * exp(-x/2) .* sin(sqrt(7)*x/2) + 35 * exp(-x/2) .* cos(sqrt(7)*x/2) - 7);
% Xr = [0, 6];
% Yb = [1, 1];
% n = 20;

[YH, YHXAM] = HXAMMain(w, Xr, Yb, n, 1, 1e-10);
[~, YAM] = HXAMMain(w, Xr, Yb, n, 1, 0);

xSpan = linspace(Xr(1), Xr(2), 100);

figure(1)
plot(YH(1, :), YH(2, :), YAM(1, :), YAM(2, :), YHXAM(1, :), YHXAM(2, :), xSpan, ySol(xSpan));
legend('Heun', 'AM', 'HXAM', 'Sol')
title('Wykres rozwiązania y(x)')
xlabel('x')
ylabel('y')
end