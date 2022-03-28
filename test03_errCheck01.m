function [] = test03_errCheck01(w, ySol, Xr, Yb, n, m1)
%% Zestaw Testowy y'' + y' +2y = x
% w = @(x)[2, 1, 1, -x];
% ySol = @(x)1 / 28 * (14 * x + 9 * sqrt(7) * exp(-x/2) .* sin(sqrt(7)*x/2) + ...
%     35 * exp(-x/2) .* cos(sqrt(7)*x/2) - 7);
% Xr = [0, 8];
% Yb = [1, 1];
% n = 10;
% m1 = 1;

%% Test
h = (Xr(2) - Xr(1)) / n;

text = ['Ilość punktów aproksymacyjnych = ', num2str(n)];
disp(text)
text = ['Wartość kroku = ', num2str(h)];
disp(text)
disp(' ')

X = linspace(Xr(1), Xr(2), n+1);
xDokl = linspace(Xr(1), Xr(2), (Xr(2) - Xr(1))/0.01);
[YH, YHXAM, itB] = HXAMMain(w, Xr, Yb, n, m1, 1e-10);
disp('Średnia liczba iteracji metodą Broydena: ')
disp(itB/n)
figure(1)
plot(X, YH(2, :), X, YHXAM(2, :), xDokl, ySol(xDokl));
legend('YHeun', 'YHXAM', 'yDokl')

errHXAM = sum(ySol(X)-YHXAM(2, :)).^2 / sum(ySol(X).^2);
disp('Wartość zagregowanego błędu względnego dla funkcji y')
disp('aproksymowanej przez AM oraz Broydena:')
disp(errHXAM)

errH = sum(ySol(X)-YH(2, :)).^2 / sum(ySol(X).^2);
disp('Wartość zagregowanego błędu względnego dla funkcji y')
disp('aproksymowanej metodą Heuna:')
disp(errH)
end