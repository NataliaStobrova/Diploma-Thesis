%vlastna zaokruhlovacia funkcia, lebo MATLAB taku nema, ze cislo 9.4111 alebo cislo 9.49999 zaokruhli na 9.5
function b = zaokruh(x)
 y = zao(x,4)*100;
 a = int2str(x*100);
 b = a(length(a));
 if str2num(b) > 0 & str2num(b) < 5 %| 
  b = zao(x,3)+0.1;
 else
  b = zao(x,3);
 end
end