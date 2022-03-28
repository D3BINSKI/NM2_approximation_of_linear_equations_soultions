function [] = test04_errCheck02(w, ySol, Xr, Yb, n, m1max)
%% Zestaw Testowy 1 y'' + y' +2y = x
% w = @(x)[2, 1,  1, -x];
% ySol = @(x)1/28 * (14*x + 9*sqrt(7)*exp(-x/2).*sin(sqrt(7)*x/2)+35*exp(-x/2).*cos(sqrt(7)*x/2)-7);
% Xr = [0, 8];
% Yb = [1, 1];
% n = 35;
% m1max = 4;

%% Test
m1 = 0:1:m1max;

err = zeros(2, length(m1));

x = linspace(Xr(1), Xr(2), n+1);

for i = 1:length(m1)
    [~, YAM] = HXAMMain(w, Xr, Yb, n, m1(i), 0);

    err(1, i) = sum((ySol(x) - YAM(2, :)).^2/sum(ySol(x).^2));
end

figure(1)
semilogy(m1, err(1, :));
xlabel('m1')
ylabel('err')
title('zintegrowany błąd względny w zależności od m1')
