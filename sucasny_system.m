%%%%%%% SUCASNY system - 1. a 2. pilier %%%%%%%
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
narast_mzdy = 0.0562;  %= input('Zadajte roènú mieru narastu mzdy: '); %priemer za roky 2017 - 2020 indexu nominalnej mzdy
priemerna_mzda_mes = 1133; %= input('Zadajte priemernu mesaènú mzdu: '); %priemerna mzda za 06/2020
priemerny_dochodok_mes = 484.94; %= input('Zadajte priemerny mesaèný dôchodok: '); %priemerny dochodok za 30.06.2020
HDP(1) = 91105000000; %v roku 2020

% % % % % % % MODELOVANIE SPORITELOV V 1. A 2. PILIERI % % % % % % % 
[Vek_strukt_m3, Vek_strukt_z3, Pp_umrtia_m, Pp_umrtia_z, dv_m, dv_z, spor2_m, spor2_z, pocet_rokov_len1_pilier_m_predpoved, pocet_rokov_len1_pilier_z_predpoved,...
 pocet_rokov_2_pilier_m_predpoved, pocet_rokov_2_pilier_z_predpoved, Vek_strukt_m, Vek_strukt_z, spor1_m, spor1_z, odvody_do_1_piliera_len_1_pilier, ...
 priemerna_mzda_roc, odvody_do_1_piliera_aj_2_pilier] = sporitelia(predpoved, vstup_na_trh_prace, dochodkovy_vek_m0, dochodkovy_vek_z0, posun_doch_vek_m, posun_doch_vek_z,...
 limit_2_pilier, miera_nezamestnanosti_m, miera_nezamestnanosti_z, valorizacia, narast_mzdy, priemerna_mzda_mes, priemerny_dochodok_mes)

% % % % % % % MODELOVANIE vyplacania dochodkov v buducnosti % % % % % % %
%Novopriznane dochodky sa rataju podla vzorca D=POMB*ADH*ODP, s POMB = 1 a ADH_0 = 14.2107
%a ODP (sucet obdobia dochodkoveho poistenia ziskaného ku dnu vzniku naroku na dochodok)
adh(1) = 14.2107; %v roku 2021
vyplacanie_doch_SP_valor_m(1,2:predpoved) = 0;  %zadefinovanie pociatocnych nul 
vyplacanie_doch_SP_valor_z(1,2:predpoved) = 0;
%muzi
for j=1:predpoved
%Posun Vekovej struktury kvoli poctu vyplacanych dochodkov
 for i=1:length(Vek_strukt_m3)
  Vek_strukt_m2(i+1) = Vek_strukt_m3(i) * (1 - Pp_umrtia_m(i));
 end
 %buduci penzisti
 for i=1:j
  adh(j+1) = adh(j) * (1 + narast_mzdy); %predpoklad, ze adh rastie ako narast_mzdy, kedze adh zavisi od priemernych miezd
  pocet_penzistov_m(dv_m(j)-dv_m(1)+1,j) = Vek_strukt_m2(dv_m(j)+1); 
  if ((j == 1) || (j > 1 & (dv_m(j) == dv_m(j-1))))
   vyplacanie_doch_SP_nove_m(j) = adh(j) * 12 * ((Vek_strukt_m2(dv_m(j)+1) - spor2_m(dv_m(j)-vstup_na_trh_prace+1,j)) * pocet_rokov_len1_pilier_m_predpoved(j,3) +...
                                  spor2_m(dv_m(j)-vstup_na_trh_prace+1,j) * pocet_rokov_2_pilier_m_predpoved(j,4) * (1 - pocet_rokov_2_pilier_m_predpoved(j,5)));    vyplacanie_doch_SP_valor_m(dv_m(j)-dv_m(1)+1,j) = vyplacanie_doch_SP_nove_m(j); %novopriznane dochodky
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
  if (j == 1)  %nevedomost o pocte predoslych dochodcov, kolkym sa valorizuju dochodky, tak sa zobrali pre 1. rok prognozy z Vekovej struktury 
   pocet_penzistov_m(i,j) = Vek_strukt_m(dv_m(1)+(i-1)) * (1 - Pp_umrtia_m(dv_m(1)+(i-1)));
  else
   pocet_penzistov_m(i,j) = pocet_penzistov_m(i-1,j-1) * (1 - Pp_umrtia_m(dv_m(1)+(i-1))); 
  end
  vyplacanie_doch_SP_valor_m(i,j) = Vek_strukt_m(dv_m(1)+(i-1)+1) * (1 - Pp_umrtia_m(dv_m(1)+(i-1)+1)) * priemerny_dochodok_roc(j+1); 
 end
 vyplacanie_doch_m(j) = sum(vyplacanie_doch_SP_valor_m(:,j));
 celkovy_pocet_penzistov_m(j) = sum(pocet_penzistov_m(:,j));
 Vek_strukt_m3 = [29824; Vek_strukt_m2(2:(end-1))'];
end

%zeny
for j=1:predpoved
%Posun Vekovej struktury kvoli poctu vyplacanych dochodkov
 for i=1:length(Vek_strukt_z3)
  Vek_strukt_z2(i+1) = Vek_strukt_z3(i) * (1 - Pp_umrtia_z(i));
 end
 %buduci penzisti
 for i=1:j
 adh(j+1) = adh(j) * (1 + narast_mzdy);  %predpoklad, ze adh rastie ako narast_mzdy, kedze adh zavisi od priemernych miezd
 pocet_penzistov_z(dv_z(j)-dv_z(1)+1,j) = Vek_strukt_z2(dv_z(j)+1); 
 if ((j == 1) || (j > 1 & (dv_z(j) == dv_z(j-1))))
  vyplacanie_doch_SP_nove_z(j) = adh(j) * 12 * ((Vek_strukt_z2(dv_z(j)+1) - spor2_z(dv_z(j)-vstup_na_trh_prace+1,j)) * pocet_rokov_len1_pilier_z_predpoved(j,3) +...
                                 spor2_z(dv_z(j)-vstup_na_trh_prace+1,j) * pocet_rokov_2_pilier_z_predpoved(j,4) * (1 - pocet_rokov_2_pilier_z_predpoved(j,5)));
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
    dochodok_len1_pilier_roc_z(i,j) = dochodok_len1_pilier_mes_z(i,j) * 12;
    dochodok_2_pilier_roc_z(i,j) = dochodok_2_pilier_mes_z(i,j) * 12;
   else
    %Zaokruhlovacia funkcia zaokruh   
    dochodok_len1_pilier_mes_z(i,j) = dochodok_len1_pilier_mes_z(i-1,j-1) + zaokruh(dochodok_len1_pilier_mes_z(i-1,j-1)*valorizacia); 
    dochodok_2_pilier_mes_z(i,j) = dochodok_2_pilier_mes_z(i-1,j-1) + zaokruh(dochodok_2_pilier_mes_z(i-1,j-1)*valorizacia);
    dochodok_len1_pilier_roc_z(i,j) = dochodok_len1_pilier_mes_z(i,j) * 12;
    dochodok_2_pilier_roc_z(i,j) = dochodok_2_pilier_mes_z(i,j) * 12;
  end
 end
  if (i > 1 & dv_z(j) - vstup_na_trh_prace + 1 + (i - 1) <= length(spor2_z))  %pre i > 1 sa valorizuju priznane dochodky
   pocet_penzistov_z(dv_z(j)-dv_z(1)+i,j) = Vek_strukt_z2(dv_z(j)+i); 
  end
  if (i > 1 & i <= 99 - dv_z(1) + 1)
   vyplacanie_doch_SP_valor_z(i,j) = dochodok_len1_pilier_roc_z(i,j) * (pocet_penzistov_z(i,j) - spor2_z(dv_z(1)-vstup_na_trh_prace+(i-1),j)) +...
                                     dochodok_2_pilier_roc_z(i,j) * spor2_z(dv_z(1)-vstup_na_trh_prace+(i-1),j);
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
 vyplacanie_doch_z(j) = sum(vyplacanie_doch_SP_valor_z(:,j));
 celkovy_pocet_penzistov_z(j) = sum(pocet_penzistov_z(:,j));
 Vek_strukt_z3 = [28330; Vek_strukt_z2(2:(end-1))'];
end
% % % % % % % MODELOVANIE prijmov Socialnej Poistovne do buducnosti % % % % % % % 
%Produktivny vek (zadane parametre)
for j=1:predpoved
 doSocPoist(j) = (sum(spor1_m(1:dv_m(j)-vstup_na_trh_prace,j)) - sum(spor2_m(1:dv_m(j)-vstup_na_trh_prace,j))) * (1 - miera_nezamestnanosti_m) * odvody_do_1_piliera_len_1_pilier * priemerna_mzda_roc(j+1)+...
                 (sum(spor1_z(1:dv_z(j)-vstup_na_trh_prace,j)) - sum(spor2_z(1:dv_z(j)-vstup_na_trh_prace,j))) * (1 - miera_nezamestnanosti_z) * odvody_do_1_piliera_len_1_pilier * priemerna_mzda_roc(j+1)+...
                  odvody_do_1_piliera_aj_2_pilier(j+1) * priemerna_mzda_roc(j+1) * ( sum(spor2_m(1:dv_m(j)-vstup_na_trh_prace,j)) * (1 - miera_nezamestnanosti_m) + sum(spor2_z(1:dv_z(j)-vstup_na_trh_prace,j)) * (1 - miera_nezamestnanosti_z) );  
 pocet_pracujucich(j) = sum(spor1_m(1:dv_m(j)-vstup_na_trh_prace,j)) * (1 - miera_nezamestnanosti_m) + sum(spor1_z(1:dv_z(j)-vstup_na_trh_prace,j)) * (1 - miera_nezamestnanosti_z);
end
%Dolezite vystupy 1 pre cele predikovane obdobie:
POCET_PRACUJUCICH = pocet_pracujucich'
POCET_DOCHODCOV_M = celkovy_pocet_penzistov_m' 
POCET_DOCHODCOV_Z = celkovy_pocet_penzistov_z'
PRIJMY_SOC_POISTOVNE_Z_ODVODOV_PRACUJUCICH = doSocPoist' 
VYDAVKY_SOC_POISTOVNE_NA_DOCHODKY = vyplacanie_doch_m' + vyplacanie_doch_z'

%Grafy prijmov / vydavkov Soc. Poistovne a poctu sporitelov / dochodcov
% figure
% scatter(pocet_rokov_2_pilier_m_predpoved(1:end,1), doSocPoist, 25, 'filled', 'g') % 'b' pre posun_doch_vek_m/z = 2
% hold on
% scatter(pocet_rokov_2_pilier_m_predpoved(1:end,1), vyplacanie_doch_m + vyplacanie_doch_z, 25, 'filled', 'r') % 'm' pre posun_doch_vek_m/z = 2
% title({'Príjmy a výdavky Sociálnej poisovne'})
% xlabel('Rok prognózy')
% ylabel('Príjmy / Výdavky Sociálnej poisovne')
% ylim([0 200000000000])
% legend({'príjmy ','výdavky'},'Location','northwest')
% legend({'príjmy (posun\_doch\_vek_{m/z} = 0)','výdavky (posun\_doch\_vek_{m/z} = 0)','príjmy (posun\_doch\_vek_{m/z} = 2)','výdavky (posun\_doch\_vek_{m/z} = 2)'},'Location','northwest')
 
% figure
% scatter(pocet_rokov_2_pilier_m_predpoved(1:end,1), pocet_pracujucich, 'filled', 'g') % 'b' pre posun_doch_vek_m/z = 2
% hold on
% scatter(pocet_rokov_2_pilier_m_predpoved(1:end,1), celkovy_pocet_penzistov_m + celkovy_pocet_penzistov_z, 'filled', 'r') % 'm' pre posun_doch_vek_m/z = 2
% title({'Poèet prispievate¾ov vs. poèet poberate¾ov'})
% xlabel('Rok prognózy')
% ylabel('Poèet prispievate¾ov / poèet poberate¾ov')
% ylim([500000 2000000])
% legend({'poèet prispievate¾ov (posun\_doch\_vek_{m/z} = 0)' ,'poèet poberate¾ov (posun\_doch\_vek_{m/z} = 0)', 'poèet prispievate¾ov (posun\_doch\_vek_{m/z} = 2)' ,'poèet poberate¾ov (posun\_doch\_vek_{m/z} = 2)'},'Location','northeast')

for i=1:predpoved
 HDP(i+1) = HDP(i) * (1 + narast_mzdy);
 deficit(i) = ((vyplacanie_doch_m(i) + vyplacanie_doch_z(i)) - doSocPoist(i)) / HDP(i+1) * 100;
end

% figure
% scatter(pocet_rokov_2_pilier_m_predpoved(1:end,1), deficit, 25, 'filled', 'r')  % 'm' pre posun_doch_vek_m/z = 2
% title({'Deficit starobného fondu Sociálnej poisovne [% HDP]'})
% xlabel('Rok prognózy')
% ylabel('Deficit starobného fondu ')
% ylim([0 5])
% legend({'deficit (posun\_doch\_vek_{m/z} = 0)', 'deficit (posun\_doch\_vek_{m/z} = 2)'},'Location','southwest')