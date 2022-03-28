function [F] = Fval(w, Y, x)
% FVAL funkcja służąca do obliczania wektora pochodnych Y w punkcie x
%
% WEJSCIE 
%       w - Wektor 1x4 zawierający w każdej z kolumn uchwyty do funkcji według
%       poniższego schematu:
%           w11 - współczynnik stojący przy y,
%           w12 - współczynnik stojący przy pierwszej pochodnej y,
%           w13 - współczynnik stojący przy drugiej pochodnej y,
%           w14 - funkcja opisująca wyraz wolny.
%       Wszystkie elementy wektora W powinny być zależne od jendej zmiennej x.
%       Y - wektor wartości funkcji w punkcie x
%       x - punkt w którym chcemy obliczyć wartości pochodnych Y
%
% WYJŚCIE
%       F - wektor zawierająacy obliczone wartości pochodnych Y w punkcie x

    wmat = w(x);
    %      y1' ,  y2'                               y3' = f(y1, y2, y3)                
    F = [   1  , Y(3), (-wmat(1, 1)*Y(2) - wmat(1, 2)*Y(3) - wmat(1, 4))/wmat(1, 3)]';
end