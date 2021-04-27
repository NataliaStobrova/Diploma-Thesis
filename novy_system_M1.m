%NOVY system (s jednotnou davkou po roku (dochodkovy_vek_m - vstup_na_trh_prace)) - 1. a 2. pilier (M1)
%pozmeneny narast_nov_dodat(i) aby deficit bol 1.5% do buducnosti alebo minimalna podmienka narast_nov_dodat(i) = 0.01 
clear all
clc
format long
%Spustenie zaokruhlovacich funkcii zaokruh a zao 
%Zadefinovanie vstupnych parametrov
predpoved = 70; %= input('Zadajte poèet rokov predpovede: '); 
vstup_na_trh_prace = 22; %= input('Zadajte vek vstupujúcich na trh práce: '); %minimalne od 16 rokov
dochodkovy_vek_m0 = 63; %= input('Zadajte vek odchodu do dochodku pre mužov: ');
dochodkovy_vek_z0 = 63; %= input('Zadajte vek odchodu do dochodku pre ženy: ');
posun_doch_vek_m = 2; %= input('Zadajte posun veku odchodu do dochodku pre mužov (v mesiacoch): ');
posun_doch_vek_z = 2; %= input('Zadajte posun veku odchodu do dochodku pre ženy (v mesiacoch): ');
limit_2_pilier = 35; %= input('Zadajte vek, do kedy môže poistenec vstúpi do 2. piliera: ');
miera_nezamestnanosti_m = 0.1053; %= input('Zadajte mieru zamestnanosti pre mužov: '); %priemer 
miera_nezamestnanosti_z = 0.1176; %= input('Zadajte mieru zamestnanosti pre ženy: '); %priemer
valorizacia = 0.026; %= input('Zadajte roènú mieru valorizácie dôchodkov: '); %za rok 2020  
narast_mzdy = 0.0562; %= input('Zadajte roènú mieru narastu mzdy: '); %priemer za roky 2017 - 2020 indexu nominalnej mzdy
priemerna_mzda_mes = 1133; %= input('Zadajte priemernu mesaènú mzdu: '); %priemerna mzda za 06/2020
priemerny_dochodok_mes = 484.94; %= input('Zadajte priemerny mesaèný dôchodok: '); %priemerny dochodok za 30.06.2020
zlom = dochodkovy_vek_m0 - vstup_na_trh_prace; %rok prechodu z aktualneho na novy system
pomer_zakl_dod = 0.7; %cast z predoslej davky (zo sucasneho modelu) 
HDP(1) = 91105000000; %v roku 2020
hdp = HDP(1);

for i=1:predpoved
 HDP(i+1) = HDP(i) * (1 + narast_mzdy);
end

for k=1:predpoved-zlom
 deficit(k) = 1.6;
 narast_nov_dodat(k) = 0.035;
 while ((narast_nov_dodat(k) > 0.01) && (deficit(k) > 1.5))  %pokym deficit > 2 alebo narast_noc_dodat >0.01, ked deficit = 1% hdp, tak doratame ake musi byt narast_nov_dodat(i)
  narast_nov_dodat(k) = narast_nov_dodat(k) - 0.005;
  narast_novy_dodat = narast_nov_dodat(k);
  [celk_m, celk_z, doSP] = def(narast_novy_dodat, k, predpoved, posun_doch_vek_m, posun_doch_vek_z, valorizacia, narast_mzdy, vstup_na_trh_prace, ...
  dochodkovy_vek_m0, dochodkovy_vek_z0, limit_2_pilier, miera_nezamestnanosti_m, miera_nezamestnanosti_z, priemerna_mzda_mes, priemerny_dochodok_mes, pomer_zakl_dod, hdp);
  HDP(zlom+k) = HDP(zlom+k-1) * (1 + narast_mzdy);
  celko_m = celk_m;
  celko_z = celk_z;
  deficit(k) = (celko_m + celko_z - doSP) / HDP(zlom+k) * 100;
  dodatok(k) = narast_novy_dodat;
 end
 dodatok(k) = narast_novy_dodat;
end

% % % % % % % MODELOVANIE SPORITELOV V 1. A 2. PILIERI % % % % % % % 
[Vek_strukt_m3, Vek_strukt_z3, Pp_umrtia_m, Pp_umrtia_z, dv_m, dv_z, spor2_m, spor2_z, pocet_rokov_len1_pilier_m_predpoved, pocet_rokov_len1_pilier_z_predpoved,...
 pocet_rokov_2_pilier_m_predpoved, pocet_rokov_2_pilier_z_predpoved, Vek_strukt_m, Vek_strukt_z, spor1_m, spor1_z, odvody_do_1_piliera_len_1_pilier, ...
 priemerna_mzda_roc, odvody_do_1_piliera_aj_2_pilier] = sporitelia(predpoved, vstup_na_trh_prace, dochodkovy_vek_m0, dochodkovy_vek_z0, posun_doch_vek_m, posun_doch_vek_z,...
 limit_2_pilier, miera_nezamestnanosti_m, miera_nezamestnanosti_z, valorizacia, narast_mzdy, priemerna_mzda_mes, priemerny_dochodok_mes)

% % % % % % % MODELOVANIE vyplacania dochodkov v buducnosti % % % % % % %
%Novopriznane dochodky sa rataju podla vzorca D=POMB*ADH*ODP, s POMB = 1 a ADH_0 = 14.2107
%a ODP (sucet obdobia dochodkoveho poistenia ziskaného ku dnu vzniku naroku na dochodok)
adh(1) = 14.2107;  %v roku 2021
vyplacanie_doch_SP_valor_m(1,2:predpoved) = 0;  %zadefinovanie pociatocnych nul 
vyplacanie_doch_SP_valor_z(1,2:predpoved) = 0;

% % % % % % muzi % % % % % % 
dochodok_mes_m = zeros(99-dv_m(1)+1,predpoved);
dodatok_len1_pilier_mes_m = zeros(99-dv_m(1)+1,predpoved);
for j=1:predpoved
%Posun Vekovej struktury kvoli poctu vyplacanych dochodkov
 for i=1:length(Vek_strukt_m3)
  Vek_strukt_m2(i+1) = Vek_strukt_m3(i) * (1 - Pp_umrtia_m(i));
 end
 if (j < zlom) %po zlom. roku sa budu vyplacat uz jednotne davky (bez kratenia dochodku z 1. piliera pre druhopilieristov)
  %buduci penzisti
  for i=1:j
   adh(j+1) = adh(j) * (1 + narast_mzdy);  %predpoklad, ze adh rastie ako narast_mzdy, kedze adh zavisi od priemernych miezd
   pocet_penzistov_m(dv_m(j)-dv_m(1)+1,j) = Vek_strukt_m2(dv_m(j)+1); 
   if ((j == 1) || (j > 1 & (dv_m(j) == dv_m(j-1))))
    vyplacanie_doch_SP_nove_m(j) = adh(j) * 12 * ((Vek_strukt_m2(dv_m(j)+1)-spor2_m(dv_m(j)-vstup_na_trh_prace+1,j)) * pocet_rokov_len1_pilier_m_predpoved(j,3) +...
                                   spor2_m(dv_m(j)-vstup_na_trh_prace+1,j) * pocet_rokov_2_pilier_m_predpoved(j,4) * (1 - pocet_rokov_2_pilier_m_predpoved(j,5)));
    vyplacanie_doch_SP_valor_m(dv_m(j)-dv_m(1)+1,j) = vyplacanie_doch_SP_nove_m(j); %novopriznane dochodky
    dochodok_len1_pilier_mes_m(dv_m(j)-dv_m(1)+1,j) = adh(j) * pocet_rokov_len1_pilier_m_predpoved(j,3);
    dochodok_2_pilier_mes_m(dv_m(j)-dv_m(1)+1,j) = adh(j) * pocet_rokov_2_pilier_m_predpoved(j,4) * (1 - pocet_rokov_2_pilier_m_predpoved(j,5)); %krateny dochodok z 1. piliera pre druhopilieristov
    dochodok_len1_pilier_roc_m(dv_m(j)-dv_m(1)+1,j) = dochodok_len1_pilier_mes_m(dv_m(j)-dv_m(1)+1,j) * 12;
    dochodok_2_pilier_roc_m(dv_m(j)-dv_m(1)+1,j) = dochodok_2_pilier_mes_m(dv_m(j)-dv_m(1)+1,j) * 12;
   end
   %Valorizovanie priznanych dochodkov
   if (i > 1 & i <= 99 - dv_m(1) + 1)
    if (dochodok_len1_pilier_mes_m(i-1,j-1) == 0) 
     dochodok_len1_pilier_mes_m(i,j) = 0;
     dochodok_2_pilier_mes_m(i,j) = 0;
     dochodok_len1_pilier_roc_m(i,j) = 0;
     dochodok_2_pilier_roc_m(i,j) = 0;
    elseif (dochodok_len1_pilier_mes_m(i-1,j-1) < 361.6)  
     dochodok_len1_pilier_mes_m(i,j) = dochodok_len1_pilier_mes_m(i-1,j-1) + 9.4;
     dochodok_2_pilier_mes_m(i,j) = dochodok_2_pilier_mes_m(i-1,j-1) + 9.4;
     dochodok_len1_pilier_roc_m(i,j) = dochodok_len1_pilier_mes_m(i,j) * 12;
     dochodok_2_pilier_roc_m(i,j) = dochodok_2_pilier_mes_m(i,j) * 12;
    else
    %Zaokruhlovacia funkcia zaokruh
     dochodok_len1_pilier_mes_m(i,j) = dochodok_len1_pilier_mes_m(i-1,j-1) + zaokruh(dochodok_len1_pilier_mes_m(i-1,j-1)*valorizacia); 
     dochodok_2_pilier_mes_m(i,j) = dochodok_2_pilier_mes_m(i-1,j-1) + zaokruh(dochodok_2_pilier_mes_m(i-1,j-1)*valorizacia);
     dochodok_len1_pilier_roc_m(i,j) = dochodok_len1_pilier_mes_m(i,j) * 12;
     dochodok_2_pilier_roc_m(i,j) = dochodok_2_pilier_mes_m(i,j) * 12;
    end
   end
   if (i > 1 & dv_m(j) - vstup_na_trh_prace + 1 + (i - 1) <= length(spor2_m)) %pre i > 1 sa valorizuju priznane dochodky
    pocet_penzistov_m(dv_m(j)-dv_m(1)+i,j) = Vek_strukt_m2(dv_m(j)+i); 
   end
   if (i > 1 & i <= 99 - dv_m(1) + 1)
    vyplacanie_doch_SP_valor_m(i,j) = dochodok_len1_pilier_roc_m(i,j) * (pocet_penzistov_m(i,j) - spor2_m(dv_m(1)-vstup_na_trh_prace+(i-1),j)) +...
                                      dochodok_2_pilier_roc_m(i,j) * spor2_m(dv_m(1)-vstup_na_trh_prace+(i-1),j);
   end
  end
  %dnesni penzisti
  for i=(j+1):99-dv_m(1)+1
   %Valorizovanie priznanych dochodkov
   if (priemerny_dochodok_mes(j) < 361.6)  
    priemerny_dochodok_mes(j+1) = priemerny_dochodok_mes(j) + 9.4;
    priemerny_dochodok_roc(j+1) = priemerny_dochodok_mes(j+1) * 12;
   else
    %Zaokruhlovacia funkcia zaokruh  
    priemerny_dochodok_mes(j+1) = priemerny_dochodok_mes(j) + zaokruh(priemerny_dochodok_mes(j)*valorizacia);
    priemerny_dochodok_roc(j+1) = priemerny_dochodok_mes(j+1) * 12;
   end  
   if (j == 1) %nevedomost o pocte predoslych dochodcov, kolkym sa valorizuju dochodky, tak sa zobrali pre 1. rok prognozy z Vekovej struktury 
    pocet_penzistov_m(i,j) = Vek_strukt_m(dv_m(1)+(i-1)) * (1 - Pp_umrtia_m(dv_m(1)+(i-1)));
   else
    pocet_penzistov_m(i,j) = pocet_penzistov_m(i-1,j-1) * (1 - Pp_umrtia_m(dv_m(1)+(i-1))); 
   end
   vyplacanie_doch_SP_valor_m(i,j) = Vek_strukt_m(dv_m(1)+(i-1)+1) * (1 - Pp_umrtia_m(dv_m(1)+(i-1)+1)) * priemerny_dochodok_roc(j+1); 
  end
 else (j >= zlom)
  if (j == zlom)
   %buduci penzisti v roku zlom
   jednotna_davka_mes(j-zlom+1) = pomer_zakl_dod * adh(zlom-1) * (1 + narast_mzdy) * pocet_rokov_len1_pilier_m_predpoved(j,3); %v 1. roku prognozy noveho
   jednotna_davka_roc(j-zlom+1) = jednotna_davka_mes(1) * 12; %v 1. roku prognozy
   dodatok_mes(j-zlom+1) = adh(zlom-1) * (1 + narast_mzdy) * pocet_rokov_len1_pilier_m_predpoved(j,3) - jednotna_davka_mes(j-zlom+1); %v 1. roku prognozy noveho
   dodatok_roc(j-zlom+1) = dodatok_mes(j-zlom+1) * 12; %v 1. roku prognozy
   dochodok_mes_m(dv_m(j)-dv_m(1)+1,j) = jednotna_davka_mes(j-zlom+1);
   dochodok_roc_m(dv_m(j)-dv_m(1)+1,j) = jednotna_davka_roc(j-zlom+1);
   dodatok_len1_pilier_mes_m(dv_m(j)-dv_m(1)+1,j) = dodatok_mes(j-zlom+1);
   dodatok_len1_pilier_roc_m(dv_m(j)-dv_m(1)+1,j) = dodatok_roc(j-zlom+1);
   vyplacanie_doch_SP_nove_m(j) = jednotna_davka_mes(j-zlom+1) * 12 * Vek_strukt_m2(dv_m(j)+1) +...
                                  (Vek_strukt_m2(dv_m(j)+1)-spor2_m(dv_m(j)-vstup_na_trh_prace+1,j)) * 12 * dodatok_mes(j-zlom+1); %niektori dlhsie niektori kratsie takze v priemere koeficient 1
   vyplacanie_doch_SP_valor_m(dv_m(j)-dv_m(1)+1,j) = vyplacanie_doch_SP_nove_m(j);  %do vektora sa priradili novopriznane dochodky
  end
  pocet_penzistov_m(dv_m(j)-dv_m(1)+1,j) = Vek_strukt_m2(dv_m(j)+1); %pocet vsetkych, ktorych sa vyplacaju dochodky
  if ((j > zlom) & (dv_m(j) == dv_m(j-1)))
   %buduci penzisti v dalsich rokoch (po roku zlom)
   jednotna_davka_mes(j-zlom+1) = jednotna_davka_mes(j-zlom) + zaokruh(jednotna_davka_mes(j-zlom)*narast_mzdy); %*narast_nov_doch(1));
   dodatok_mes(j-zlom+1) = dodatok_mes(j-zlom) + zaokruh(dodatok_mes(j-zlom)*dodatok(j-zlom));
   vyplacanie_doch_SP_nove_m(j) = jednotna_davka_mes(j-zlom+1) * 12 * Vek_strukt_m2(dv_m(j)+1) +...
                                  (Vek_strukt_m2(dv_m(j)+1)-spor2_m(dv_m(j)-vstup_na_trh_prace+1,j)) * 12 * dodatok_mes(j-zlom+1); %niektori dlhsie niektori kratsie takze v priemere koeficient 1
   vyplacanie_doch_SP_valor_m(dv_m(j)-dv_m(1)+1,j) = vyplacanie_doch_SP_nove_m(j); %do vektora sa priradili novopriznane dochodky
   dochodok_mes_m(dv_m(j)-dv_m(1)+1,j) = jednotna_davka_mes(j-zlom+1); %rovnaka zakladna davka pre vsetkych
   dochodok_roc_m(dv_m(j)-dv_m(1)+1,j) = jednotna_davka_mes(j-zlom+1) * 12;
   dodatok_len1_pilier_mes_m(dv_m(j)-dv_m(1)+1,j) = dodatok_mes(j-zlom+1);
   dodatok_len1_pilier_roc_m(dv_m(j)-dv_m(1)+1,j) = dodatok_mes(j-zlom+1) * 12;
  end
  %dnesni penzisti s valorizovanou jednotnou davkou   
  for i=dv_m(j)-dv_m(1)+2:dv_m(j)-dv_m(1)+1+j-zlom 
   if (i > dv_m(j) - dv_m(1) + 1 & i <= 99 - dv_m(1) + 1)
    if (dochodok_mes_m(i-1,j-1) == 0) 
     dochodok_mes_m(i,j) = 0;
     dodatok_len1_pilier_mes_m(i,j) = 0;
     dochodok_roc_m(i,j) = 0;
     dodatok_len1_pilier_roc_m(i,j) = 0;
    elseif (dochodok_mes_m(i-1,j-1) < 361.6)  
     dochodok_mes_m(i,j) = dochodok_mes_m(i-1,j-1) + 9.4;
     dodatok_len1_pilier_mes_m(i,j) = dodatok_len1_pilier_mes_m(i-1,j-1) + 9.4;
     dochodok_roc_m(i,j) = dochodok_mes_m(i,j) * 12;
     dodatok_len1_pilier_roc_m(i,j) = dodatok_len1_pilier_mes_m(i,j) * 12;
    else
    %Zaokruhlovacia funkcia zaokruh
     dochodok_mes_m(i,j) = dochodok_mes_m(i-1,j-1) + zaokruh(dochodok_mes_m(i-1,j-1)*valorizacia); 
     dodatok_len1_pilier_mes_m(i,j) = dodatok_len1_pilier_mes_m(i-1,j-1) + zaokruh(dodatok_len1_pilier_mes_m(i-1,j-1) * dodatok(j-zlom));
     dochodok_roc_m(i,j) = dochodok_mes_m(i,j) * 12;
     dodatok_len1_pilier_roc_m(i,j) = dodatok_len1_pilier_mes_m(i,j) * 12;
    end
   end
   if (i > dv_m(j) - dv_m(1) + 1 & dv_m(1) - vstup_na_trh_prace + 1 + (i - 1) <= length(spor2_m)) %pre i > 1 sa valorizuju priznane dochodky
    pocet_penzistov_m(i,j) = Vek_strukt_m2(dv_m(j)+i-(dv_m(j)-dv_m(1)+1)+1); 
   end
   if (i > dv_m(j) - dv_m(1) + 1 & i <= 99 - dv_m(1) + 1)
    vyplacanie_doch_SP_valor_m(i,j) = dodatok_len1_pilier_roc_m(i,j) * (pocet_penzistov_m(i,j) - spor2_m(dv_m(1)-vstup_na_trh_prace+i,j)) +...
                                      pocet_penzistov_m(i,j) * dochodok_roc_m(i,j);
   end
  end
  for i=dv_m(j)-dv_m(1)+1+j-zlom+1:99-dv_m(1)+1
   %Valorizovanie priznanych dochodkov %Zaokruhlovacia funkcia zaokruh  
   if (dochodok_len1_pilier_mes_m(i-1,j-1) == 0)
    dochodok_len1_pilier_mes_m(i,j) = 0;
    dochodok_2_pilier_mes_m(i,j) = 0;
    dochodok_len1_pilier_roc_m(i,j) = 0;
    dochodok_2_pilier_roc_m(i,j) = 0;
   elseif (dochodok_len1_pilier_mes_m(i-1,j-1) < 361.6)  
    dochodok_len1_pilier_mes_m(i,j) = dochodok_len1_pilier_mes_m(i-1,j-1) + 9.4;
    dochodok_2_pilier_mes_m(i,j) = dochodok_2_pilier_mes_m(i-1,j-1) + 9.4;
    dochodok_len1_pilier_roc_m(i,j) = dochodok_len1_pilier_mes_m(i,j) * 12;
    dochodok_2_pilier_roc_m(i,j) = dochodok_2_pilier_mes_m(i,j) * 12;
   else
   %Zaokruhlovacia funkcia zaokruh
    dochodok_len1_pilier_mes_m(i,j) = dochodok_len1_pilier_mes_m(i-1,j-1) + zaokruh(dochodok_len1_pilier_mes_m(i-1,j-1)*valorizacia); 
    dochodok_2_pilier_mes_m(i,j) = dochodok_2_pilier_mes_m(i-1,j-1) + zaokruh(dochodok_2_pilier_mes_m(i-1,j-1)*valorizacia);
    dochodok_len1_pilier_roc_m(i,j) = dochodok_len1_pilier_mes_m(i,j) * 12;
    dochodok_2_pilier_roc_m(i,j) = dochodok_2_pilier_mes_m(i,j) * 12;
   end   
   if (i > dv_m(j) - dv_m(1) + 1 & dv_m(1) - vstup_na_trh_prace + 1 + (i - 1) <= length(spor2_m)) %pre i > 1 sa valorizuju priznane dochodky
    pocet_penzistov_m(i,j) = Vek_strukt_m2(dv_m(1)+i); 
   end
   if (i > dv_m(j) - dv_m(1) + 1 & i <= 99 - dv_m(1) + 1)
    vyplacanie_doch_SP_valor_m(i,j) = dochodok_len1_pilier_roc_m(i,j) * (pocet_penzistov_m(i,j) - spor2_m(dv_m(1)-vstup_na_trh_prace+i,j)) +...
                                      spor2_m(dv_m(1)-vstup_na_trh_prace+i,j) * dochodok_2_pilier_roc_m(i,j);
   end
  end
 end
 vyplacanie_doch_m(j) = sum(vyplacanie_doch_SP_valor_m(:,j));
 celkovy_pocet_penzistov_m(j) = sum(pocet_penzistov_m(:,j));
 Vek_strukt_m3 = [29824; Vek_strukt_m2(2:(end-1))'];
end

% % % % % % zeny % % % % % % (analogicky ako pre muzov)
dochodok_mes_z = zeros(99-dv_z(1)+1,predpoved);
dodatok_len1_pilier_mes_z = zeros(99-dv_z(1)+1,predpoved);
for j=1:predpoved
%Posun Vekovej struktury kvoli poctu vyplacanych dochodkov
 for i=1:length(Vek_strukt_z3)
  Vek_strukt_z2(i+1) = Vek_strukt_z3(i) * (1 - Pp_umrtia_z(i));
 end
 if (j < zlom) %po zlom. roku sa budu vyplacat uz jednotne davky (bez kratenia dochodku z 1. piliera pre druhopilieristov)
 %buduci penzisti
  for i=1:j
   adh(j+1) = adh(j) * (1 + narast_mzdy); %predpoklad, ze adh rastie ako narast_mzdy, kedze adh zavisi od priemernych miezd
   pocet_penzistov_z(dv_z(j)-dv_z(1)+1,j) = Vek_strukt_z2(dv_z(j)+1); 
   if ((j == 1) || (j > 1 & (dv_z(j) == dv_z(j-1))))
    vyplacanie_doch_SP_nove_z(j) = adh(j) * 12 * ((Vek_strukt_z2(dv_z(j)+1)-spor2_z(dv_z(j)-vstup_na_trh_prace+1,j)) * pocet_rokov_len1_pilier_z_predpoved(j,3) +...
                                   spor2_z(dv_z(j)-vstup_na_trh_prace+1,j) * pocet_rokov_2_pilier_z_predpoved(j,4) * (1-pocet_rokov_2_pilier_z_predpoved(j,5)));
    vyplacanie_doch_SP_valor_z(dv_z(j)-dv_z(1)+1,j) = vyplacanie_doch_SP_nove_z(j);  %novopriznane dochodky
    dochodok_len1_pilier_mes_z(dv_z(j)-dv_z(1)+1,j) = adh(j) * pocet_rokov_len1_pilier_z_predpoved(j,3);
    dochodok_2_pilier_mes_z(dv_z(j)-dv_z(1)+1,j) = adh(j) * pocet_rokov_2_pilier_z_predpoved(j,4) * (1 - pocet_rokov_2_pilier_z_predpoved(j,5));  %krateny dochodok z 1. piliera pre druhopilieristov
    dochodok_len1_pilier_roc_z(dv_z(j)-dv_z(1)+1,j) = dochodok_len1_pilier_mes_z(dv_z(j)-dv_z(1)+1,j) * 12;
    dochodok_2_pilier_roc_z(dv_z(j)-dv_z(1)+1,j) = dochodok_2_pilier_mes_z(dv_z(j)-dv_z(1)+1,j) * 12;
   end
   if (i > 1 & i <= 99 - dv_m(1) + 1)
   %Valorizovanie priznanych dochodkov
    if (dochodok_len1_pilier_mes_z(i-1,j-1) == 0)
     dochodok_len1_pilier_mes_z(i,j) = 0;
     dochodok_2_pilier_mes_z(i,j) = 0;
     dochodok_len1_pilier_roc_z(i,j) = 0;
     dochodok_2_pilier_roc_z(i,j) = 0;
    elseif (dochodok_len1_pilier_mes_z(i-1,j-1) < 361.6)
     dochodok_len1_pilier_mes_z(i,j) = dochodok_len1_pilier_mes_z(i-1,j-1) + 9.4;
     dochodok_2_pilier_mes_z(i,j) = dochodok_2_pilier_mes_z(i-1,j-1) + 9.4;
     dochodok_len1_pilier_roc_z(i,j) = 12*dochodok_len1_pilier_mes_z(i,j);
     dochodok_2_pilier_roc_z(i,j) = 12*dochodok_2_pilier_mes_z(i,j);
    else
    %Zaokruhlovacia funkcia zaokruh   
     dochodok_len1_pilier_mes_z(i,j) = dochodok_len1_pilier_mes_z(i-1,j-1) + zaokruh(dochodok_len1_pilier_mes_z(i-1,j-1)*valorizacia); 
     dochodok_2_pilier_mes_z(i,j) = dochodok_2_pilier_mes_z(i-1,j-1) + zaokruh(dochodok_2_pilier_mes_z(i-1,j-1)*valorizacia);
     dochodok_len1_pilier_roc_z(i,j) = dochodok_len1_pilier_mes_z(i,j) *12;
     dochodok_2_pilier_roc_z(i,j) = dochodok_2_pilier_mes_z(i,j) *12;
    end
   end
   if (i > 1 & dv_z(j) - vstup_na_trh_prace + 1 + (i - 1) <= length(spor2_z)) %pre i > 1 sa valorizuju priznane dochodky
    pocet_penzistov_z(dv_z(j)-dv_z(1)+i,j) = Vek_strukt_z2(dv_z(j)+i); 
   end
   if (i > 1 & i <= 99 - dv_z(1) + 1)
    vyplacanie_doch_SP_valor_z(i,j) = dochodok_len1_pilier_roc_z(i,j) * (pocet_penzistov_z(i,j) - spor2_z(dv_z(1)-vstup_na_trh_prace+(i-1),j)) +...
                                      spor2_z(dv_z(1)-vstup_na_trh_prace+(i-1),j) * dochodok_2_pilier_roc_z(i,j);
   end
  end
  %dnesni penzisti
  for i=(j+1):99-dv_z(1)+1
   %Valorizovanie priznanych dochodkov
   if (priemerny_dochodok_mes(j) < 361.6)
    priemerny_dochodok_mes(j+1) = priemerny_dochodok_mes(j) + 9.4;
    priemerny_dochodok_roc(j+1) = priemerny_dochodok_mes(j+1) * 12;
   else
   %Zaokruhlovacia funkcia zaokruh  
    priemerny_dochodok_mes(j+1) = priemerny_dochodok_mes(j) + zaokruh(priemerny_dochodok_mes(j)*valorizacia); 
    priemerny_dochodok_roc(j+1) = priemerny_dochodok_mes(j+1) * 12;
   end        
   if (j == 1)  %nevedomost o pocte predoslych dochodcov, kolkym sa valorizuju dochodky, tak sa zobrali pre 1. rok prognozy z Vekovej struktury
    pocet_penzistov_z(i,j) = Vek_strukt_z(dv_z(1)+(i-1)) * (1 - Pp_umrtia_z(dv_z(1)+(i-1)));
   else
    pocet_penzistov_z(i,j) = pocet_penzistov_z(i-1,j-1) * (1 - Pp_umrtia_z(dv_z(1)+(i-1))); 
   end
   vyplacanie_doch_SP_valor_z(i,j) = Vek_strukt_z(dv_z(1)+(i-1)) * (1 - Pp_umrtia_z(dv_z(1)+(i-1))) * priemerny_dochodok_roc(j+1);
  end
 else (j >= zlom)
  if (j == zlom)
   %buduci penzisti v roku zlom
   jednotna_davka_mes(j-zlom+1) = pomer_zakl_dod * adh(zlom-1) * (1 + narast_mzdy) * pocet_rokov_len1_pilier_z_predpoved(j,3); %v 1. roku prognozy novom
   jednotna_davka_roc(j-zlom+1) = jednotna_davka_mes(1) * 12; %v 1. roku prognozy
   dodatok_mes(j-zlom+1) = adh(zlom-1) * (1 + narast_mzdy) * pocet_rokov_len1_pilier_z_predpoved(j,3) - jednotna_davka_mes(j-zlom+1); %v 1. roku prognozy novom
   dodatok_roc(j-zlom+1) = dodatok_mes(j-zlom+1) * 12; %v 1. roku prognozy
   dochodok_mes_z(dv_z(j)-dv_z(1)+1,j) = jednotna_davka_mes(j-zlom+1);
   dochodok_roc_z(dv_z(j)-dv_z(1)+1,j) = jednotna_davka_roc(j-zlom+1);
   dodatok_len1_pilier_mes_z(dv_z(j)-dv_z(1)+1,j) = dodatok_mes(j-zlom+1);
   dodatok_len1_pilier_roc_z(dv_z(j)-dv_z(1)+1,j) = dodatok_roc(j-zlom+1);
   vyplacanie_doch_SP_nove_z(j) = jednotna_davka_mes(j-zlom+1) * 12 * Vek_strukt_z2(dv_z(j)+1) +...
                                  (Vek_strukt_z2(dv_z(j)+1)-spor2_z(dv_z(j)-vstup_na_trh_prace+1,j)) * 12 * dodatok_mes(j-zlom+1);
   vyplacanie_doch_SP_valor_z(dv_z(j)-dv_z(1)+1,j) = vyplacanie_doch_SP_nove_z(j);  
  end
  pocet_penzistov_z(dv_z(j)-dv_z(1)+1,j) = Vek_strukt_z2(dv_z(j)+1); %pocet vsetkych, ktorych sa vyplacaju dochodky
  if ((j > zlom) & (dv_z(j) == dv_z(j-1)))
   %buduci penzisti v dalsich rokoch (po roku zlom)
   jednotna_davka_mes(j-zlom+1) = jednotna_davka_mes(j-zlom) + zaokruh(jednotna_davka_mes(j-zlom)*narast_mzdy);
   dodatok_mes(j-zlom+1) = dodatok_mes(j-zlom) + zaokruh(dodatok_mes(j-zlom)*dodatok(j-zlom));
   vyplacanie_doch_SP_nove_z(j) = jednotna_davka_mes(j-zlom+1) * 12 * Vek_strukt_z2(dv_z(j)+1) +...
                                  (Vek_strukt_z2(dv_z(j)+1)-spor2_z(dv_z(j)-vstup_na_trh_prace+1,j)) * 12 * dodatok_mes(j-zlom+1); %niektori dlhsie niektori kratsie takze v priemere koeficient 1
   vyplacanie_doch_SP_valor_z(dv_z(j)-dv_z(1)+1,j) = vyplacanie_doch_SP_nove_z(j);  %do vektora sa priradili novopriznane dochodky
   dochodok_mes_z(dv_z(j)-dv_z(1)+1,j) = jednotna_davka_mes(j-zlom+1); %rovnaka zakladna davka pre vsetkych
   dochodok_roc_z(dv_z(j)-dv_z(1)+1,j) = jednotna_davka_mes(j-zlom+1) * 12;
   dodatok_len1_pilier_mes_z(dv_z(j)-dv_z(1)+1,j) = dodatok_mes(j-zlom+1);
   dodatok_len1_pilier_roc_z(dv_z(j)-dv_z(1)+1,j) = dodatok_mes(j-zlom+1) * 12; 
  end
  for i=dv_z(j)-dv_z(1)+2:dv_z(j)-dv_z(1)+1+j-zlom 
   %dnesni penzisti s valorizovanou jednotnou davkou 
   if (i > dv_z(j) - dv_z(1) + 1 & i <= 99 - dv_z(1) + 1)
    if (dochodok_mes_z(i-1,j-1) == 0) 
     dochodok_mes_z(i,j) = 0;
     dodatok_len1_pilier_mes_z(i,j) = 0;
     dochodok_roc_z(i,j) = 0;
     dodatok_len1_pilier_roc_z(i,j) = 0;
    elseif (dodatok_len1_pilier_mes_z(i-1,j-1) < 361.6)  
     dochodok_mes_z(i,j) = dochodok_mes_z(i-1,j-1)+ 9.4;
     dodatok_len1_pilier_mes_z(i,j) = dodatok_len1_pilier_mes_z(i-1,j-1) + 9.4;
     dochodok_roc_z(i,j) = dochodok_mes_z(i,j) * 12;
     dodatok_len1_pilier_roc_z(i,j) = dodatok_len1_pilier_mes_z(i,j) * 12;
    else
     %Zaokruhlovacia funkcia zaokruh
     dochodok_mes_z(i,j) = dochodok_mes_z(i-1,j-1) + zaokruh(dochodok_mes_z(i-1,j-1)*valorizacia); 
     dodatok_len1_pilier_mes_z(i,j) = dodatok_len1_pilier_mes_z(i-1,j-1) + zaokruh(dodatok_len1_pilier_mes_z(i-1,j-1)*dodatok(j-zlom)); 
     dochodok_roc_z(i,j) = dochodok_mes_z(i,j) * 12;
     dodatok_len1_pilier_roc_z(i,j) = dodatok_len1_pilier_mes_z(i,j) * 12;
    end
   end
   if (i > dv_z(j) - dv_z(1) + 1 & dv_z(1) - vstup_na_trh_prace + 1 + (i - 1) <= length(spor2_z)) %pre i > 1 sa valorizuju priznane dochodky
    pocet_penzistov_z(i,j) = Vek_strukt_z2(dv_z(j)+i-(dv_z(j)-dv_z(1)+1)+1); 
   end
   if (i > dv_z(j) - dv_z(1) + 1 & i <= 99 - dv_z(1) + 1)
    vyplacanie_doch_SP_valor_z(i,j) = dodatok_len1_pilier_roc_z(i,j) * (pocet_penzistov_z(i,j) - spor2_z(dv_z(1)-vstup_na_trh_prace+i,j)) +...
                                      pocet_penzistov_z(i,j) * dochodok_roc_z(i,j);
   end
  end
 %dnesni penzisti
  for i=dv_z(j)-dv_z(1)+1+j-zlom+1:99-dv_z(1)+1
   %Valorizovanie priznanych dochodkov %Zaokruhlovacia funkcia zaokruh  
   if (dochodok_len1_pilier_mes_z(i-1,j-1) == 0)
    dochodok_len1_pilier_mes_z(i,j) = 0;
    dochodok_2_pilier_mes_z(i,j) = 0;
    dochodok_len1_pilier_roc_z(i,j) = 0;
    dochodok_2_pilier_roc_z(i,j) = 0;
   elseif (dochodok_len1_pilier_mes_z(i-1,j-1) < 361.6)  
    dochodok_len1_pilier_mes_z(i,j) = dochodok_len1_pilier_mes_z(i-1,j-1) + 9.4;
    dochodok_2_pilier_mes_z(i,j) = dochodok_2_pilier_mes_z(i-1,j-1) + 9.4;
    dochodok_len1_pilier_roc_z(i,j) = dochodok_len1_pilier_mes_z(i,j) * 12;
    dochodok_2_pilier_roc_z(i,j) = dochodok_2_pilier_mes_z(i,j) * 12;
   else
   %Zaokruhlovacia funkcia zaokruh
    dochodok_len1_pilier_mes_z(i,j) = dochodok_len1_pilier_mes_z(i-1,j-1) + zaokruh(dochodok_len1_pilier_mes_z(i-1,j-1)*valorizacia); 
    dochodok_2_pilier_mes_z(i,j) = dochodok_2_pilier_mes_z(i-1,j-1) + zaokruh(dochodok_2_pilier_mes_z(i-1,j-1)*valorizacia);
    dochodok_len1_pilier_roc_z(i,j) = dochodok_len1_pilier_mes_z(i,j) * 12;
    dochodok_2_pilier_roc_z(i,j) = dochodok_2_pilier_mes_z(i,j) * 12;
   end   
   if (i > dv_z(j) - dv_z(1) + 1 & dv_z(1) - vstup_na_trh_prace + 1 + (i - 1) <= length(spor2_z)) %pre i > 1 sa valorizuju priznane dochodky
    pocet_penzistov_z(i,j) = Vek_strukt_z2(dv_z(1)+i); 
   end
   if (i > dv_z(j) - dv_z(1) + 1 & i <= 99 - dv_z(1) + 1)
    vyplacanie_doch_SP_valor_z(i,j) = dochodok_len1_pilier_roc_z(i,j) * (pocet_penzistov_z(i,j) - spor2_z(dv_z(1)-vstup_na_trh_prace+i,j)) +...
                                      spor2_z(dv_z(1)-vstup_na_trh_prace+i,j) * dochodok_2_pilier_roc_z(i,j);
   end
  end
 end
 vyplacanie_doch_z(j) = sum(vyplacanie_doch_SP_valor_z(:,j));
 celkovy_pocet_penzistov_z(j) = sum(pocet_penzistov_z(:,j));
 Vek_strukt_z3 = [28330; Vek_strukt_z2(2:(end-1))'];
end

% MODELOVANIE prijmov Socialnej Poistovne do buducnosti % % % % % % % 
% Produktivny vek (zadane parametre)
for j=1:predpoved
 doSocPoist(j) = (sum(spor1_m(1:dv_m(j)-vstup_na_trh_prace,j)) - sum(spor2_m(1:dv_m(j)-vstup_na_trh_prace,j))) * (1 - miera_nezamestnanosti_m) * odvody_do_1_piliera_len_1_pilier * priemerna_mzda_roc(j+1)+...
                 (sum(spor1_z(1:dv_z(j)-vstup_na_trh_prace,j)) - sum(spor2_z(1:dv_z(j)-vstup_na_trh_prace,j))) * (1 - miera_nezamestnanosti_z) * odvody_do_1_piliera_len_1_pilier * priemerna_mzda_roc(j+1)+...
                 odvody_do_1_piliera_aj_2_pilier(j+1) * priemerna_mzda_roc(j+1) * (sum(spor2_m(1:dv_m(j)-vstup_na_trh_prace,j)) * (1 - miera_nezamestnanosti_m) + sum(spor2_z(1:dv_z(j)-vstup_na_trh_prace,j)) * (1 - miera_nezamestnanosti_z));  
 pocet_pracujucich(j) = sum(spor1_m(1:dv_m(j)-vstup_na_trh_prace,j)) * (1 - miera_nezamestnanosti_m) + sum(spor1_z(1:dv_z(j)-vstup_na_trh_prace,j)) * (1 - miera_nezamestnanosti_z);
end

%Dolezite vystupy 1 pre cele predikovane obdobie:
POCET_PRACUJUCICH = pocet_pracujucich'
POCET_DOCHODCOV_M = celkovy_pocet_penzistov_m' 
POCET_DOCHODCOV_Z = celkovy_pocet_penzistov_z';
PRIJMY_SOC_POISTOVNE_Z_ODVODOV_PRACUJUCICH = doSocPoist' 
VYDAVKY_SOC_POISTOVNE_NA_DOCHODKY = vyplacanie_doch_m' + vyplacanie_doch_z'

%Grafy prijmov / vydavkov Soc. Poistovne a poctu sporitelov / dochodcov
% figure
% scatter(pocet_rokov_2_pilier_m_predpoved(1:end,1), doSocPoist, 'filled', 'g') % 'b' pre posun_doch_vek_m/z = 2
% hold on
% scatter(pocet_rokov_2_pilier_m_predpoved(1:end,1), vyplacanie_doch_m + vyplacanie_doch_z, 25, 'filled', 'MarkerFaceColor', [1 0. 0.5]) %[0.6350 0.0780 0.1840]) %[1 0. 0.5] pre posun_doch_vek_m/z = 2 
% title({'Príjmy a výdavky Sociálnej poisovne'})
% xlabel('Rok prognózy')
% ylabel('Príjmy / Výdavky Sociálnej poisovne')
% ylim([0 200000000000])
% legend({'príjmy ','výdavky'},'Location','northwest')
  
% figure
% scatter(pocet_rokov_2_pilier_m_predpoved(1:end,1), pocet_pracujucich, 'filled', 'g')
% hold on
% scatter(pocet_rokov_2_pilier_m_predpoved(1:end,1), celkovy_pocet_penzistov_m + celkovy_pocet_penzistov_z, 'filled', 'r')
% title({'Poèet prispievate¾ov vs. poèet poberate¾ov'})
% xlabel('Rok prognózy')
% ylabel('Poèet prispievate¾ov / poèet poberate¾ov')
% ylim([500000 2000000])
% legend({'poèet prispievate¾ov', 'poèet poberate¾ov'},'Location','northeast')

for i=1:predpoved
 HDP(i+1) = HDP(i) * (1 + narast_mzdy);
 deficit(i) = ((vyplacanie_doch_m(i) + vyplacanie_doch_z(i)) - doSocPoist(i)) / HDP(i+1) * 100;
end

% figure
% scatter(pocet_rokov_2_pilier_m_predpoved(1:end,1), deficit, 25, 'filled', 'MarkerFaceColor', [1 0. 0.5]) %[0.6350 0.0780 0.1840]) %[1 0. 0.5] pre posun_doch_vek_m/z = 2 
% title({'Deficit starobného fondu Sociálnej poisovne [% HDP]'})
% xlabel('Rok prognózy')
% ylabel('Deficit starobného fondu ')
% ylim([0 5])