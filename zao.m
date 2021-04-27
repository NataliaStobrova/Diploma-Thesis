function [zaokruhlene] = zao(c,presnost)
if c==0
    l=0;
else
l=log10(abs(c));
end
r=floor(l)+1;
c=c*10^(-r);
v=round(c,presnost);
zaokruhlene=v*10^(r);
end