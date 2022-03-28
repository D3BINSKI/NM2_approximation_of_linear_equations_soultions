function [Y] = mHeun(w, Y0, h)
% MHEUN obliczanie przybliżonej wartości równania różniczkowego 2-go rzędu
% w kolejnym węźle
%
% WEJŚCIE
%       w - Wektor 1x4 zawierający w każdej z kolumn uchwyty do funkcji 
%           według poniższego schematu:
%           w11 - współczynnik stojący przy y,
%           w12 - współczynnik stojący przy pierwszej pochodnej y,
%           w13 - współczynnik stojący przy drugiej pochodnej y,
%           w14 - funkcja opisująca wyraz wolny.
%       Y0 - wektor zawierający wartości poprzednio aproksymowanego
%       punktu
%       h - wartość kroku 
%
% WYJŚCIE
%       Y - wektora zawierający przybliżone wartości równania różniczkowego
%       w kolejnym węźle

F = Fval(w, Y0, Y0(1));
K = Fval(w, Y0+F*h, Y0(1)+h);
Y = Y0 + h / 2 * (F + K);
end