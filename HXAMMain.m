function [YH, YHXAM, itD] = HXAMMain(w, Xr, Yb, n, m1, eps)
% autor: Dominik Dębińśki 304 246
% HXAMMAIN   Przybliżanie rozwiązań liniowych równań różniczkowych 2-go rzędu.
%   HXAMMain(w, Xr, Yb, n, m1, m2) przybliża rozwiązanie liniowych
%   równań różniczkowych 2-go rzędu przez zastosowanie metody Heun'a oraz niejawnej
%   metody Eulera 2-go rzędu (metoda A-M 2-go rzędu). Wartości kolejnych
%   przybliżeń dla metody niejawnej obliczane są poprzez wykonanie
%   m1 (>= 0) kroków przez korektor (metoda iteracji prostej) oraz wykonanie
%   dodatkowej korekcji (metoda Broydena) w przypadku, gdy eps ~= 0. Macierzą
%   startową dla obliczenia przybliżeń jest macierz jednostkowa.
%
%   WEJŚCIE
%       w - Wektor 1x4 zawierający w każdej z kolumn uchwyty do funkcji według
%       poniższego schematu:
%           w11 - współczynnik stojący przy y,
%           w12 - współczynnik stojący przy pierwszej pochodnej y,
%           w13 - współczynnik stojący przy drugiej pochodnej y,
%           w14 - funkcja opisująca wyraz wolny.
%       Wszystkie elementy wektora W powinny być zależne od jendej zmiennej x.
%       Xr - wektor 1x2 zawierający wartości określające początek oraz
%       koniec przedziału, na którym przybliżane jest rozwiązanie.
%           Xr = [a, b] => a - początek przedziału,
%                          b - koniec przedziału
%       Yb - wektor 1x2, którego wartości orkeślają warunki początkowe
%       zadania.
%           Yb = [c, d] => c - wartość y0,
%                          d - wartość y'0
%       n - parametr definiujący liczbę przedziałów, na jaki zostanie
%       podzielony obszar określony wektorem Xr.
%       m1 - ilość kroków, które zostaną wykonane przez metodę iteracji
%       prostej w procesie obliczania kolejnego przybliżenia. Wartość
%       minimalna tego parametru wynosi 1.
%       eps - wartość maksymalnej różnicy pomiędzy kolejnymi korekcjami
%       wektora w metodzie Broydena (warunek zakończenia działania metody)
%
%   WYJŚCIE
%       YH - wektor przybliżonych rozwiązań liniowego równania różniczkowego,
%       oblicznoych w punktach leżących w przedziale wyznaczonym przez Xr
%       przy pomocy metody Heuna. Ilość punktów zależna od n.
%       YHXAM - wektor przybliżonych rozwiązań liniowego równania różniczkowego,
%       oblicznoych w punktach leżących w przedziale wyznaczonym przez Xr
%       przy pomocy metody predyktor-korektor (predytkor - metoda Heuna, 
%       korektor1 - metoda AM, korektor2 - Broyden). Ilość punktów zależna od n.
%       itD - sumaryczna liczba iteracji metody Broydena.

%% Sprawdzenie przekazanych parametrów
if nargin < 3
    help HXAMMain
    error("Zbyt mała lizcba argumentów")
end
if nargin < 4
    n = 20;
    m1 = 1;
    eps = 1e-10;
end
%% HEUN
Y = zeros(3, n);

Y(1, 1) = Xr(1, 1);
Y(2, 1) = Yb(1, 1);
Y(3, 1) = Yb(1, 2);

% Obliczanie wartości kroku
h = (Xr(1, 2) - Xr(1, 1)) / n;

for i = 2:n + 1
    Y(:, i) = mHeun(w, Y(:, i-1), h);
end
YH = Y;

%% AM + Broyden
%Inicjacja wektora pomocniczego zerami
Y = zeros(3, n);

% Zapisanie wartości początkowych do wektora pomocniczego
Y(1, 1) = Xr(1, 1); %           |  x0   |
Y(2, 1) = Yb(1, 1); %  Y(:,1) = | y'(0) |
Y(3, 1) = Yb(1, 2); %           | y''(0)|

% iterator pojedynczy (zawiera ilość iteracji pojedyńczego
% wywoładnia metody broyden)
itP = 0;
% iterator docelowy (zawiera sumaryczną wartość ilości iteracji
% wszystkich wywołań metody Broydena
itD = 0;

for i = 2:n + 1
    % Heun (PREDYKTOR)
    %  Wyznacz aproksymację wartości funkcji y w punkcie x_i
    Y(:, i) = mHeun(w, Y(:, i-1), h);
    % AM (KOREKTOR I)
    %  Popraw wstępne przybliżenie m1 razy
    for j = 1:m1
        Y(:, i) = Y(:, i-1) + h / 2 * (Fval(w, Y(:, i), Y(1, i)) + Fval(w, Y(:, i-1), Y(1, i-1)));
    end
    % Broyden (KOREKTOR II)
    %  Jeżeli użytkownik podał wartość maksymalnej różnicy pomiędzy
    %  kolejnymi przybliżeniami różną od zera, to wykonaj dodatkową korekcję
    %  metodą Broyden'a (m1 może również wynosić 0, wtedy jedynym
    %  zastosowanym korektorem będzie metoda Broydena)
    if eps ~= 0
        [Y(:, i), itP] = broyden(w, Y(:, i-1), Y(:, i), eps);
    end
    itD = itP + itD;
end
YHXAM = Y;
end