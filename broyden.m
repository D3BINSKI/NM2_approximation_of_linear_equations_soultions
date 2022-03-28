function [Y, it] = broyden(w, Y0, Y1, eps)
% BROYDEN funkcja pomocnicza implementująca metodę korekcji aproksymacji
% kolejnego wyrazu Y1
%
% WEJŚCIE 
%       w - Wektor 1x4 zawierający w każdej z kolumn uchwyty do funkcji według
%       poniższego schematu:
%           w11 - współczynnik stojący przy y,
%           w12 - współczynnik stojący przy pierwszej pochodnej y,
%           w13 - współczynnik stojący przy drugiej pochodnej y,
%           w14 - funkcja opisująca wyraz wolny.
%       Wszystkie elementy wektora 'w' powinny być zależne od jendej zmiennej x.
%       Y0 - wektor zawierający wartości poprzednio aproksymowanego
%       punktu
%       Y1 - wketor zawierający wstępne przybliżenie, które chcemy poprawić
%       eps - wartość maksymalnej różnicy pomiędzy kolejnymi przybliżeniami
%       Y1 (warunek zakończenia działania metody)
%
% WYJŚCIE
%       Y - wektor zawierający wyjściowe wartości przybliżenia równiania różniczkowego
%       w aktualnym węźle
%       it - liczba iteracji metodą Broyden'a, która była konieczna do
%       spełnienia warunku abs(Y1(2) - Y(2)) < abs(eps)

B = eye(3);
it = 0;
h = Y1(1,1)-Y0(1,1);
F = @(Y0, Y) Y - h / 2 * (Fval(w, Y, Y(1)) + Fval(w, Y0, Y0(1))) - Y0;

while(1)
   it = it + 1;
   f1old = F(Y0, Y1);
   Y = Y1 - B*f1old;
   s = Y - Y1;
   z = F(Y0, Y) - f1old;
   B = B+ (s-B*z) * (s'*B / (s'*B*z));
   if abs(Y1(2) - Y(2)) < abs(eps)
       break
   end
   Y1 = Y;
end
end