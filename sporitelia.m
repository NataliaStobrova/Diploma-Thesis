function [Vek_strukt_m3, Vek_strukt_z3, Pp_umrtia_m, Pp_umrtia_z, dv_m, dv_z, spor2_m, spor2_z, pocet_rokov_len1_pilier_m_predpoved, pocet_rokov_len1_pilier_z_predpoved,...
 pocet_rokov_2_pilier_m_predpoved, pocet_rokov_2_pilier_z_predpoved, Vek_strukt_m, Vek_strukt_z, spor1_m, spor1_z, odvody_do_1_piliera_len_1_pilier,...
 priemerna_mzda_roc, odvody_do_1_piliera_aj_2_pilier]  = sporitelia(predpoved, vstup_na_trh_prace, dochodkovy_vek_m0, dochodkovy_vek_z0, posun_doch_vek_m, posun_doch_vek_z,...
 limit_2_pilier, miera_nezamestnanosti_m, miera_nezamestnanosti_z, valorizacia, narast_mzdy, priemerna_mzda_mes, priemerny_dochodok_mes)
priemerna_mzda_mes(1) = priemerna_mzda_mes;
priemerna_mzda_roc(1) = priemerna_mzda_mes(1) * 12;
priemerny_dochodok_mes(1) = priemerny_dochodok_mes;
priemerny_dochodok_roc(1) = priemerny_dochodok_mes(1) * 12;

%Vlozenie dat z excelu
Vek_strukt_m0 = xlsread('F:\diplomová práca\Vekova_struktura_2019.xlsx',1, 'A2:A101');  %v kazdej zlozke je od veku 0 do veku 99 pocet ludi (100 zlozkovy vektor)
Vek_strukt_z0 = xlsread('F:\diplomová práca\Vekova_struktura_2019.xlsx',1, 'B2:B101');  %na 100. zlozke je pre pocet ludi vo veku 99 rokov
Pp_umrtia_m = xlsread('F:\diplomová práca\data_DP.xlsx',1, 'D3:D102');
Pp_umrtia_z = xlsread('F:\diplomová práca\data_DP.xlsx',1, 'E3:E102');
Vek_strukt_m = Vek_strukt_m0;
Vek_strukt_z = Vek_strukt_z0;
Vek_strukt_19_20_m = 29824; %znama najnovsia porodnost v roku 2019
Vek_strukt_19_20_z = 28330; %pre rok 2020 sme pouzili porodnost z roku 2019

%Posuv Vekovej struktury o rok z 2019 na 2020   
for i=1:length(Vek_strukt_m0)-1
 Vek_strukt_m2(i+1) = Vek_strukt_m(i) * (1 - Pp_umrtia_m(i));
 Vek_strukt_z2(i+1) = Vek_strukt_z(i) * (1 - Pp_umrtia_z(i));
end
Vek_strukt_m = Vek_strukt_m2;
Vek_strukt_m(1) = Vek_strukt_19_20_m;
Vek_strukt_z = Vek_strukt_z2;
Vek_strukt_z(1) = Vek_strukt_19_20_z;
Vek_strukt_m = Vek_strukt_m';
Vek_strukt_z = Vek_strukt_z';
Vek_strukt_m3 = Vek_strukt_m;
Vek_strukt_z3 = Vek_strukt_z;

%Zadefinovanie tabulky krátenia dochodku z 1. piliera pre druhopilieristov
kratenie_dochodku = xlsread('F:\diplomová práca\data_DP.xlsx',2, 'C2:C11'); %koeficient kratenia
kratenie_dochodku = [kratenie_dochodku, [7 5 1 1 1 1 1 1 1 1]']; %pocet rokov z povodnej tabulky od 2005-2012, 2012-2016, ...

poistenci_prvy_pilier_m = xlsread('F:\diplomová práca\Poistenci_31.12.2019.xlsx',1, 'I3:I85'); %od 15 do 97 rokov
poistenci_prvy_pilier_z = xlsread('F:\diplomová práca\Poistenci_31.12.2019.xlsx',1, 'J3:J85'); %od 15 do 97 rokov
pocet_poistencov_v_1_pilieri_m = sum(poistenci_prvy_pilier_m(vstup_na_trh_prace-15+1:83)); %suma od vstup_na_trh_prace po 97 rokov (83 = 97 - 15 + 1)
pocet_poistencov_v_1_pilieri_z = sum(poistenci_prvy_pilier_z(vstup_na_trh_prace-15+1:83)); %suma od vstup_na_trh_prace po 97 rokov

%Pociatocny stav poistencov v 1. pilieri (od vstup_na_trh_prace po 97 rokov) vypocitany ako:  
%sumarny kumulativny priemer (ako faktor vhodny na skalovanie) * vekova struktura obyvatelstva
poistenci_v_1_pilieri_m = Vek_strukt_m(vstup_na_trh_prace+1:98) * pocet_poistencov_v_1_pilieri_m / sum(Vek_strukt_m(vstup_na_trh_prace+1:98)); %od vstup_na_trh_prace po 97
poistenci_v_1_pilieri_z = Vek_strukt_z(vstup_na_trh_prace+1:98) * pocet_poistencov_v_1_pilieri_z / sum(Vek_strukt_z(vstup_na_trh_prace+1:98)); %od vstup_na_trh_prace po 97
%Poistenci z 1. piliera v 0. roku a pocet poistencov z 1. piliera v 0. roku
poistenci_v_1_pilieri_m;
poistenci_v_1_pilieri_z;
sum(poistenci_v_1_pilieri_m) + sum(poistenci_v_1_pilieri_z);

%Pre mladsich ako vstup_na_trh_prace, sa da z Vekovej struktury vyratat pocet vstup_na_trh_prace - rocnych prvopilieristov pre kazdy rok predpovede
for j=1:vstup_na_trh_prace
 for i=j:-1:1
  if (j == i)
   mladsi_v_1_pilieri_m(j) = Vek_strukt_m(vstup_na_trh_prace-i+1) * (1 - Pp_umrtia_m(vstup_na_trh_prace-i+1));
   mladsi_v_1_pilieri_z(j) = Vek_strukt_z(vstup_na_trh_prace-i+1) * (1 - Pp_umrtia_z(vstup_na_trh_prace-i+1));  
  else
   mladsi_v_1_pilieri_m(j) = mladsi_v_1_pilieri_m(j) * (1 - Pp_umrtia_m(vstup_na_trh_prace-i+1));
   mladsi_v_1_pilieri_z(j) = mladsi_v_1_pilieri_z(j) * (1 - Pp_umrtia_z(vstup_na_trh_prace-i+1));
  end
 end
end
mladsi_v_1_pilieri_m = mladsi_v_1_pilieri_m * pocet_poistencov_v_1_pilieri_m / sum(Vek_strukt_m(vstup_na_trh_prace+1:98));
mladsi_v_1_pilieri_z = mladsi_v_1_pilieri_z * pocet_poistencov_v_1_pilieri_z / sum(Vek_strukt_z(vstup_na_trh_prace+1:98));

%Grafy pociatocneho stavu prvopilieristov
% figure
% x = [vstup_na_trh_prace:1:97]
% barh(x,poistenci_v_1_pilieri_m, 'b')
% title({'Veková štruktúra poistencov v 1. pilieri (muži)'})
% xlabel('Poèet sporite¾ov v danom veku')
% ylabel('Vek')
% legend({'poistenci\_v\_1\_pilieri_m (muži)'})
  
% figure
% x = [vstup_na_trh_prace:1:97]
% barh(x,poistenci_v_1_pilieri_z, 'm')
% title({'Veková štruktúra poistencov v 1. pilieri (ženy)'})
% xlabel('Poèet sporite¾ov v danom veku')
% ylabel('Vek')
% legend({'poistenci\_v\_1\_pilieri_z (ženy)'})

sporitelia_druhy_pilier_m = xlsread('F:\diplomová práca\Vekova struktura II. pilier_31.12.2019.xlsx',1, 'I3:I75'); %od 16 do 88 rokov
sporitelia_druhy_pilier_z = xlsread('F:\diplomová práca\Vekova struktura II. pilier_31.12.2019.xlsx',1, 'J3:J75'); %od 16 do 88 rokov
pocet_sporitelov_v_2_pilieri_m = sum(sporitelia_druhy_pilier_m(vstup_na_trh_prace-16+1:73)); %od vstup_na_trh_prace po 88 rokov (73 = 88 - 16 + 1)
pocet_sporitelov_v_2_pilieri_z = sum(sporitelia_druhy_pilier_z(vstup_na_trh_prace-16+1:73)); %od vstup_na_trh_prace po 88 rokov

%Pociatocny stav sporitelov v 2. pilieri (od vstup_na_trh_prace po 88 rokov) vypocitany ako:                                    
%sumarny kumulativny priemer (ako faktor vhodny na skalovanie) * poistenci_v_1_pilieri 
sporitelia_v_2_pilieri_m = poistenci_v_1_pilieri_m(1:88-vstup_na_trh_prace+1) * pocet_sporitelov_v_2_pilieri_m / sum(poistenci_v_1_pilieri_m(1:88-vstup_na_trh_prace+1)); %1.zlozka vstup_na_trh_prace
sporitelia_v_2_pilieri_z = poistenci_v_1_pilieri_z(1:88-vstup_na_trh_prace+1) * pocet_sporitelov_v_2_pilieri_z / sum(poistenci_v_1_pilieri_z(1:88-vstup_na_trh_prace+1)); %posledna zlozka 88 - vstup_na_trh_prace + 1

%Pre mladsich ako vstup_na_trh_prace, sa da z Vekovej struktury vyratat pocet vstup_na_trh_prace - rocnych druhopilieristov 
%pre kazdy rok predpovede do buducnosti
mladsi_v_2_pilieri_m = mladsi_v_1_pilieri_m * pocet_sporitelov_v_2_pilieri_m / sum(poistenci_v_1_pilieri_m(1:88-vstup_na_trh_prace+1));
mladsi_v_2_pilieri_z = mladsi_v_1_pilieri_z * pocet_sporitelov_v_2_pilieri_z / sum(poistenci_v_1_pilieri_z(1:88-vstup_na_trh_prace+1));
mladsi_v_2_pilieri_m0 = mladsi_v_2_pilieri_m;
mladsi_v_2_pilieri_z0 = mladsi_v_2_pilieri_z;

%Grafy pociatocneho stavu druhopilieristov 
% figure
% x = [vstup_na_trh_prace:1:88]
% barh(x,sporitelia_v_2_pilieri_m, 'b')
% title({'Veková štruktúra sporite¾ov v 2. pilieri (muži)'})
% xlabel('Poèet sporite¾ov v danom veku')
% ylabel('Vek')
% legend({'sporitelia\_v\_2\_pilieri_m (muži)'})
  
% figure
% x = [vstup_na_trh_prace:1:88]
% barh(x,sporitelia_v_2_pilieri_z, 'm')
% title({'Veková štruktúra sporite¾ov v 2. pilieri (ženy)'})
% xlabel('Poèet sporite¾ov v danom veku')
% ylabel('Vek')
% legend({'sporitelia\_v\_2\_pilieri_z (ženy)'})

% % % % % % % MODELOVANIE NARASTU SPORITELOV DO 2. PILIERA % % % % % % % 
%Narast poistencov je modelovany tak, aby vzdy populacny rocnik, v case dovrsenia veku 35 rokov bol na 90% v 2. pilieri
dochodkovy_vek_m = dochodkovy_vek_m0;
dochodkovy_vek_z = dochodkovy_vek_z0;
for j=1:predpoved  %j = 1 znamena posunut sa dalej o 1 rok / 1. rok prognozy
 poistenci_v_1_pilieri_m_2 = zeros(99-vstup_na_trh_prace+1,1);  
 poistenci_v_1_pilieri_z_2 = zeros(99-vstup_na_trh_prace+1,1);  
 for i=1:length(poistenci_v_1_pilieri_m)
  poistenci_v_1_pilieri_m_2(i+1) = poistenci_v_1_pilieri_m(i) * (1 - Pp_umrtia_m(vstup_na_trh_prace+i)); %lebo Pp_umrtia zacina od 0-rocnych 
  poistenci_v_1_pilieri_z_2(i+1) = poistenci_v_1_pilieri_z(i) * (1 - Pp_umrtia_z(vstup_na_trh_prace+i));
  %Podmienka pri dovrseni limit_2_pilier rokov na 90% v 2. pilieri   
  if (i == limit_2_pilier - vstup_na_trh_prace + 1) 
   sporitelia_2_pilier_90_percent_m(j) = 0.9 * poistenci_v_1_pilieri_m_2(i); %v 1. zlozke je (limit_druhy_pilier-1)-rocny a aky treba narast, ked v case dovrsenia 35 ma byt 90% v 2. pilieri
   sporitelia_2_pilier_90_percent_z(j) = 0.9 * poistenci_v_1_pilieri_z_2(i);
  end
 end
 if (j <= vstup_na_trh_prace) %od 0 po vstup_na_trh_prace je znamna Vekova struktura
  poistenci_v_1_pilieri_m = [mladsi_v_1_pilieri_m(j); poistenci_v_1_pilieri_m_2(2:99-vstup_na_trh_prace+1)]; %po 99, lebo ti 100-rocni zomru
  poistenci_v_1_pilieri_z = [mladsi_v_1_pilieri_z(j); poistenci_v_1_pilieri_z_2(2:99-vstup_na_trh_prace+1)]; %v 2. zlozke su vstup_na_trh_prace-rocni
 else %dalej do budunosti je Vekova struktura neznama, tak predpoklad ze rovnaky pocet ako minuly rok
  poistenci_v_1_pilieri_m = [mladsi_v_1_pilieri_m(length(mladsi_v_1_pilieri_m)); poistenci_v_1_pilieri_m_2(2:99-vstup_na_trh_prace+1)]; %po 99, lebo ti 100-rocni zomru
  poistenci_v_1_pilieri_z = [mladsi_v_1_pilieri_z(length(mladsi_v_1_pilieri_z)); poistenci_v_1_pilieri_z_2(2:99-vstup_na_trh_prace+1)]; 
 end 
 %Vypocet narastu sporitelov v 2. pilieri 
 if (j <= limit_2_pilier - vstup_na_trh_prace)
  for w=1:j
   if (j ~= w)
    sporitelia_v_2_pilieri_m(limit_2_pilier-vstup_na_trh_prace+1-j) = sporitelia_v_2_pilieri_m(limit_2_pilier-vstup_na_trh_prace+1-j) * (1 - Pp_umrtia_m(limit_2_pilier-(j-w)));
    sporitelia_v_2_pilieri_z(limit_2_pilier-vstup_na_trh_prace+1-j) = sporitelia_v_2_pilieri_z(limit_2_pilier-vstup_na_trh_prace+1-j) * (1 - Pp_umrtia_z(limit_2_pilier-(j-w)));
   else
    narast_sporitelov_2_pilier_m0(j) = (sporitelia_2_pilier_90_percent_m(j) - sporitelia_v_2_pilieri_m(limit_2_pilier-vstup_na_trh_prace+1-j) * (1 - Pp_umrtia_m(limit_2_pilier))) / (j);
    narast_sporitelov_2_pilier_z0(j) = (sporitelia_2_pilier_90_percent_z(j) - sporitelia_v_2_pilieri_z(limit_2_pilier-vstup_na_trh_prace+1-j) * (1 - Pp_umrtia_z(limit_2_pilier))) / (j);
   end
  end 
 elseif (j > limit_2_pilier - vstup_na_trh_prace & j <= limit_2_pilier)
  mladsi_v_2_pilieri_m = mladsi_v_2_pilieri_m0;
  mladsi_v_2_pilieri_z = mladsi_v_2_pilieri_z0;
  for w=vstup_na_trh_prace+1:1:limit_2_pilier-1
   mladsi_v_2_pilieri_m(j-(limit_2_pilier-vstup_na_trh_prace)) = mladsi_v_2_pilieri_m(j-(limit_2_pilier-vstup_na_trh_prace)) * (1 - Pp_umrtia_m(w));
   mladsi_v_2_pilieri_z(j-(limit_2_pilier-vstup_na_trh_prace)) = mladsi_v_2_pilieri_z(j-(limit_2_pilier-vstup_na_trh_prace)) * (1 - Pp_umrtia_z(w));
  end
  narast_sporitelov_2_pilier_m0(j) = (sporitelia_2_pilier_90_percent_m(j) - mladsi_v_2_pilieri_m(j-(limit_2_pilier-vstup_na_trh_prace)) * (1 - Pp_umrtia_m(limit_2_pilier))) / (j);
  narast_sporitelov_2_pilier_z0(j) = (sporitelia_2_pilier_90_percent_z(j) - mladsi_v_2_pilieri_z(j-(limit_2_pilier-vstup_na_trh_prace)) * (1 - Pp_umrtia_z(limit_2_pilier))) / (j);
 else
  narast_sporitelov_2_pilier_m0(j) = narast_sporitelov_2_pilier_m0(limit_2_pilier);
  narast_sporitelov_2_pilier_z0(j) = narast_sporitelov_2_pilier_z0(limit_2_pilier);
 end 
end
narast_sporitelov_2_pilier_m0 = narast_sporitelov_2_pilier_m0(end:-1:1)'; %posledna zlozka je pocet ludi, ktori pribudli do 2. piliera v 1. roku prognozy
narast_sporitelov_2_pilier_z0 = narast_sporitelov_2_pilier_z0(end:-1:1)'; %tychto ludi pripocitame, aby sme dosiahli predpoklad o 90% v 2. pilieri

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Opatovne zadefinovanie pociatocneho stavu prvopilieristov i druhopilieristov
poistenci_v_1_pilieri_m_0 = Vek_strukt_m(vstup_na_trh_prace+1:98) * pocet_poistencov_v_1_pilieri_m / sum(Vek_strukt_m(vstup_na_trh_prace+1:98)); %od vstup_na_trh_prace po 97
poistenci_v_1_pilieri_z_0 = Vek_strukt_z(vstup_na_trh_prace+1:98) * pocet_poistencov_v_1_pilieri_z / sum(Vek_strukt_z(vstup_na_trh_prace+1:98)); %od vstup_na_trh_prace po 97
poistenci_v_1_pilieri_m = [poistenci_v_1_pilieri_m_0; zeros(99-97,1)]; %2 zlozky nulove, aby sme mali obmedzenie po 99 rokov, kedze predpoklad, ze v 100 vsetci zomru
poistenci_v_1_pilieri_z = [poistenci_v_1_pilieri_z_0; zeros(99-97,1)];
sporitelia_v_2_pilieri_m = poistenci_v_1_pilieri_m(1:88-vstup_na_trh_prace+1) * pocet_sporitelov_v_2_pilieri_m / sum(poistenci_v_1_pilieri_m(1:88-vstup_na_trh_prace+1)); %1. zlozka je vstup_na_trh_prace
sporitelia_v_2_pilieri_z = poistenci_v_1_pilieri_z(1:88-vstup_na_trh_prace+1) * pocet_sporitelov_v_2_pilieri_z / sum(poistenci_v_1_pilieri_z(1:88-vstup_na_trh_prace+1)); %posledna zlozka je 88 - vstup_na_trh_prace + 1
posun_doch_vek_m = posun_doch_vek_m/12; %prehodenie mesiacov na roky
posun_doch_vek_z = posun_doch_vek_z/12; %prehodenie mesiacov na roky
odvody_do_1_piliera_len_1_pilier = 0.18; %ti, co maju len 1. pilier tak odvod do Soc. Poistovne v 0. roku bol 0.18
odvody_do_1_piliera_aj_2_pilier(1) = 0.13; %ti, co maju aj 2. pilier tak odvod do Soc. Poistovne v 0. roku 0.13
odvody_do_2_piliera_aj_2_pilier(1) = odvody_do_1_piliera_len_1_pilier - odvody_do_1_piliera_aj_2_pilier(1); %ti, co maju aj 2. pilier tak odvod do DSS v 0. roku 0.05

% % % % % % % MODELOVANIE 1. A 2. PILIERA % % % % % % % 
for j=1:predpoved  %j = 0 znamena, ze sme v roku 2020
%Osoby, ktore boli aspon raz dochodkovo poistene (1. pilier), sa mozu rozhodnut pre vstup do 2. piliera do dovrsenia veku 35 rokov 
%tento parameter limit_druhy_pilier sa da menit a dalsi predpoklad ze poistenci do 2. piliera neodchadzaju ale len prichadzaju 
%Odvody do Socialnej poistovne
 if (j < 4)
  odvody_do_1_piliera_aj_2_pilier(j+1) = odvody_do_1_piliera_aj_2_pilier(j) - 0.0025; %j = 1 znamena rok 2021 / o 0.25% rocne sa znizuju odvody do 1. piliera
  odvody_do_2_piliera_aj_2_pilier(j+1) = odvody_do_1_piliera_len_1_pilier - odvody_do_1_piliera_aj_2_pilier(j+1);
 elseif (j >= 4)
  odvody_do_1_piliera_aj_2_pilier(j+1) = 0.12;
  odvody_do_2_piliera_aj_2_pilier(j+1) = 0.06; 
 end
 poistenci_v_1_pilieri_m_2 = zeros(99-vstup_na_trh_prace+1,1);  %zadefinovanie pomocnej premennej
 poistenci_v_1_pilieri_z_2 = zeros(99-vstup_na_trh_prace+1,1); 
 %Podmienka na narast sporitelov
 if (j <= limit_2_pilier - vstup_na_trh_prace)
  narast_sporitelov_2_pilier_m = narast_sporitelov_2_pilier_m0(((predpoved-j+2)-(limit_2_pilier-vstup_na_trh_prace)):(predpoved-j+1)); 
  narast_sporitelov_2_pilier_z = narast_sporitelov_2_pilier_z0(((predpoved-j+2)-(limit_2_pilier-vstup_na_trh_prace)):(predpoved-j+1));
 else
  narast_sporitelov_2_pilier_m = [narast_sporitelov_2_pilier_m(1); narast_sporitelov_2_pilier_m(1:end)];
  narast_sporitelov_2_pilier_z = [narast_sporitelov_2_pilier_z(1); narast_sporitelov_2_pilier_z(1:end)];
 end
 %Vyvoj prvopilieristov
 for i=1:99-vstup_na_trh_prace+1  
  poistenci_v_1_pilieri_m_2(i+1) = poistenci_v_1_pilieri_m(i) * (1 - Pp_umrtia_m(vstup_na_trh_prace+i)); %lebo Pp_umrtia zacina od 0-rocnych  
  poistenci_v_1_pilieri_z_2(i+1) = poistenci_v_1_pilieri_z(i) * (1 - Pp_umrtia_z(vstup_na_trh_prace+i));
 end
 %Vyvoj druhopilieristov
 for i=1:length(sporitelia_v_2_pilieri_m)
  if (i <= limit_2_pilier - vstup_na_trh_prace) %sporitelia do 2. piliera pribudaju od veku vstup_na_trh_prace do veku limit_druhy_pilier rokov
   sporitelia_v_2_pilieri_m_2(i+1) = sporitelia_v_2_pilieri_m(i) * (1 - Pp_umrtia_m(vstup_na_trh_prace+i)) + narast_sporitelov_2_pilier_m(i);
   sporitelia_v_2_pilieri_z_2(i+1) = sporitelia_v_2_pilieri_z(i) * (1 - Pp_umrtia_z(vstup_na_trh_prace+i)) + narast_sporitelov_2_pilier_z(i);  
  elseif (i > limit_2_pilier - vstup_na_trh_prace && i < 99 - vstup_na_trh_prace + 1) %predpoklad, ze vsetci 100-rocni zomru
   sporitelia_v_2_pilieri_m_2(i+1) = sporitelia_v_2_pilieri_m(i) * (1 - Pp_umrtia_m(vstup_na_trh_prace+i));
   sporitelia_v_2_pilieri_z_2(i+1) = sporitelia_v_2_pilieri_z(i) * (1 - Pp_umrtia_z(vstup_na_trh_prace+i));
  end
  %Pomocna matica na vypocet kratenia dochodkov - priemerny pocet rokov v 2. pilieri (neskor potrebne)
  if (i <= limit_2_pilier - vstup_na_trh_prace)
   matica_pomoc_m(i,j) = sporitelia_v_2_pilieri_m_2(i+1) - sporitelia_v_2_pilieri_m(i);
   matica_pomoc_z(i,j) = sporitelia_v_2_pilieri_z_2(i+1) - sporitelia_v_2_pilieri_z(i);
  end
 end
 if (j <= vstup_na_trh_prace)
  poistenci_v_1_pilieri_m = [mladsi_v_1_pilieri_m(j); poistenci_v_1_pilieri_m_2(2:99-vstup_na_trh_prace+1)];  %po 99, lebo predpoklad ze vsetci 100-rocni zomru
  poistenci_v_1_pilieri_z = [mladsi_v_1_pilieri_z(j); poistenci_v_1_pilieri_z_2(2:99-vstup_na_trh_prace+1)];  %v 2. zlozke su vstup_na_trh_prace-rocni
  sporitelia_v_2_pilieri_m = [mladsi_v_2_pilieri_m(j); sporitelia_v_2_pilieri_m_2(2:end)'];
  sporitelia_v_2_pilieri_z = [mladsi_v_2_pilieri_z(j); sporitelia_v_2_pilieri_z_2(2:end)'];   
 else 
  poistenci_v_1_pilieri_m = [mladsi_v_1_pilieri_m(length(mladsi_v_1_pilieri_m)); poistenci_v_1_pilieri_m_2(2:99-vstup_na_trh_prace+1)];  %po 99, lebo predpoklad ze vsetci 100-rocni zomru
  poistenci_v_1_pilieri_z = [mladsi_v_1_pilieri_z(length(mladsi_v_1_pilieri_z)); poistenci_v_1_pilieri_z_2(2:99-vstup_na_trh_prace+1)];  %v 2. zlozke su vstup_na_trh_prace-rocni 
  sporitelia_v_2_pilieri_m = [mladsi_v_2_pilieri_m0(length(mladsi_v_2_pilieri_m0)); sporitelia_v_2_pilieri_m_2(2:end)'];
  sporitelia_v_2_pilieri_z = [mladsi_v_2_pilieri_z0(length(mladsi_v_2_pilieri_z0)); sporitelia_v_2_pilieri_z_2(2:end)']; 
 end
 %Uprava dochodkoveho veku - napr. ak dochodkovy_vek_m2 < 62.5 tak bude 62 a ak >=62.5 bude 63
 if (round(dochodkovy_vek_m,2) < floor(dochodkovy_vek_m) + 0.5)
  dochodkovy_vek_m2 = floor(dochodkovy_vek_m);
 else
  dochodkovy_vek_m2 = ceil(dochodkovy_vek_m);
 end
 if (round(dochodkovy_vek_z,2) < floor(dochodkovy_vek_z) + 0.5)
  dochodkovy_vek_z2 = floor(dochodkovy_vek_z);
 else
  dochodkovy_vek_z2 = ceil(dochodkovy_vek_z);
 end
 dochodkovy_vek_m = dochodkovy_vek_m + posun_doch_vek_m;
 dochodkovy_vek_z = dochodkovy_vek_z + posun_doch_vek_z; 
 %Tabulky prvopilieristov a druhopilieristov
 for i=1:99-vstup_na_trh_prace+1 
  spor1_m(i,j) = poistenci_v_1_pilieri_m(i);
  spor1_z(i,j) = poistenci_v_1_pilieri_z(i);
  if (i <= length(sporitelia_v_2_pilieri_m))
   spor2_m(i,j) = sporitelia_v_2_pilieri_m(i); 
   spor2_z(i,j) = sporitelia_v_2_pilieri_z(i);
  else
   spor2_m(i,j) = 0; 
   spor2_z(i,j) = 0;  
  end
  if (i <= dochodkovy_vek_m2 - vstup_na_trh_prace+1)
   sp2_m(i,j) = sporitelia_v_2_pilieri_m(i);
   sp1_m(i,j) = poistenci_v_1_pilieri_m(i);
  end
  if (i <= dochodkovy_vek_z2 - vstup_na_trh_prace+1)
   sp2_z(i,j) = sporitelia_v_2_pilieri_z(i);
   sp1_z(i,j) = poistenci_v_1_pilieri_z(i);
  end   
 end 
 %Zvysovanie priemernej mzdy o parameter narast_mzdy
 priemerna_mzda_roc(j+1) = priemerna_mzda_roc(j) * (1 + narast_mzdy);
end

% % % % % % % VYPOCET udajov pre druhopilieristov % % % % % % % 
%Vysvetlenie jednotlivych stlpcov vo vektore pocet_rokov_2_pilier_m_predpoved / pocet_rokov_2_pilier_z_predpoved poporadi:
%1 - rok v ktorom dovrsia doch_vek                    %4 - pocet rokov v 1. pilieri
%2 - pocet rokov v 2. pilieri                         %5 - kratenie dochodku z 1. piliera
%3 - pocet sporitelov v 2.pilieri v doch_veku 

%Vypocet poctu rokov v 2. pilieri
%Pre vstup_na_trh_prace rocnych a starsich v 0. roku predpovede 
dochodkovy_vek_m = dochodkovy_vek_m0;
dochodkovy_vek_z = dochodkovy_vek_z0;
for j=1:limit_2_pilier-vstup_na_trh_prace 
 j = j-1;
 suma_m = 0;
 sumaa_m = 0;
 suma_z = 0;
 sumaa_z = 0;
 for i=1:limit_2_pilier-vstup_na_trh_prace  
  suma_m = suma_m + (dochodkovy_vek_m - vstup_na_trh_prace - i + 1) * matica_pomoc_m(i,i+j); %ti co pribudli krat pocet rokov ostavajucich do dochodkoveho veku
  sumaa_m = sumaa_m + matica_pomoc_m(i,i+j); %pocet ludi kolko celkovo pribudlo do veku limit_2_pilier
  suma_z = suma_z + (dochodkovy_vek_z - vstup_na_trh_prace - i + 1) * matica_pomoc_z(i,i+j);
  sumaa_z = sumaa_z + matica_pomoc_z(i,i+j);
 end
 stlpec_poctu_rok_pil_m(j+1) = round(suma_m/sumaa_m); %pocet rokov v 2. pilieri
 stlpec_poctu_rok_pil_z(j+1) = round(suma_z/sumaa_z);
 dochodkovy_vek_m = dochodkovy_vek_m + posun_doch_vek_m;
 dochodkovy_vek_z = dochodkovy_vek_z + posun_doch_vek_z;
end
%Pre mladsich a rovnych ako vstup_na_trh_prace rokov je predpoklad, ze budu rovnaky pocet rokov v 2. pilieri ako predosli sporitelia
%Pre starsich ako limit_2_pilier v 0. roku predpovede 
for i=1:dochodkovy_vek_m0-42 %preto 42, lebo ti mali v roku 2013 35 rokov, pokial limit_2_pilier
 od_15_m(i) = 16 + (i - 1);
end
for i=1:dochodkovy_vek_z0-42 %preto 42, lebo ti mali v roku 2013 35 rokov, pokial limit_2_pilier
 od_15_z(i) = 16 + (i - 1);
end

%Pridanie poctu rokov v 2. pilieri do vektora pocet_rokov_2_pilier_m_predpoved
ind_m = 0;
ind_z = 0;
for j=1:predpoved
 %muzi
 if (j <= length(od_15_m)) 
  pocet_rokov_2_pilier_m_predpoved(j) = od_15_m(j); %od 15 do 35 rokov v 2. pilieri pre 27 az 47 rocnych v 2005 (42 - 62 r. v 2020)
 elseif (j > length(od_15_m) & j <= dochodkovy_vek_m0 - limit_2_pilier + 1)
  pocet_rokov_2_pilier_m_predpoved(j) = round((stlpec_poctu_rok_pil_m(length(stlpec_poctu_rok_pil_m))+od_15_m(length(od_15_m)))/2);   
 elseif (j > dochodkovy_vek_m0 - limit_2_pilier + 1 & j <= dochodkovy_vek_m0 - vstup_na_trh_prace + 1) %preto 42, lebo ti mali v roku 2013 35 rokov, pokial limit_2_pilier
  pocet_rokov_2_pilier_m_predpoved(j) = stlpec_poctu_rok_pil_m(length(stlpec_poctu_rok_pil_m)-ind_m); %staci 1. zlozka, lebo ostatne su rovnake
  ind_m = ind_m + 1; 
 else
  pocet_rokov_2_pilier_m_predpoved(j) = stlpec_poctu_rok_pil_m(1);
 end
 %zeny
 if (j <= length(od_15_z)) 
  pocet_rokov_2_pilier_z_predpoved(j) = od_15_z(j); %od 15 do 35 rokov v 2. pilieri pre 27 až 47 rocnych v 2005 (42 - 62 r. v 2020)
 elseif (j > length(od_15_z) & j <= dochodkovy_vek_z0 - limit_2_pilier + 1)
  pocet_rokov_2_pilier_z_predpoved(j) = round((stlpec_poctu_rok_pil_z(length(stlpec_poctu_rok_pil_z))+od_15_z(length(od_15_z)))/2);   
 elseif (j > dochodkovy_vek_z0 - limit_2_pilier + 1 & j <= dochodkovy_vek_z0 - vstup_na_trh_prace + 1) %preto 42, lebo ti mali v roku 2013 35 rokov, pokial limit_2_pilier
  pocet_rokov_2_pilier_z_predpoved(j) = stlpec_poctu_rok_pil_z(length(stlpec_poctu_rok_pil_z)-ind_z); %staci 1. zlozka, lebo ostatne su rovnake
  ind_z = ind_z + 1; 
 else
  pocet_rokov_2_pilier_z_predpoved(j) = stlpec_poctu_rok_pil_z(1);
 end
end
pocet_rokov_2_pilier_m_predpoved = pocet_rokov_2_pilier_m_predpoved;
pocet_rokov_2_pilier_z_predpoved = pocet_rokov_2_pilier_z_predpoved;

%Uprava dochodkoveho veku - napr. ak dochodkovy_vek_m2 < 62.5 tak bude 62 a
%ak >62.5 bude 63 + zadefinovanie premennych dv_m / dv_z
dochodkovy_vek_m = dochodkovy_vek_m0;
dochodkovy_vek_z = dochodkovy_vek_z0;
for j=1:predpoved
 if (round(dochodkovy_vek_m,2) < floor(dochodkovy_vek_m)+0.5)
  pocet_m(j) = sp2_m(floor(dochodkovy_vek_m) - vstup_na_trh_prace+1,j);
  dv_m(j) = floor(dochodkovy_vek_m);
 else
  pocet_m(j) = sp2_m(ceil(dochodkovy_vek_m) - vstup_na_trh_prace+1,j);
  dv_m(j) = ceil(dochodkovy_vek_m);    
 end
 if (dv_m(j) >= 67)
  dv_m(j) = 67;
 end
 if (round(dochodkovy_vek_z,2) < floor(dochodkovy_vek_z)+0.5)
  pocet_z(j) = sp2_z(floor(dochodkovy_vek_z) - vstup_na_trh_prace+1,j);
  dv_z(j) = floor(dochodkovy_vek_z);
 else
  pocet_z(j) = sp2_z(ceil(dochodkovy_vek_z) - vstup_na_trh_prace+1,j);
  dv_z(j) = ceil(dochodkovy_vek_z); 
 end
 if (dv_z(j) >= 67)
  dv_z(j) = 67;
 end
 dochodkovy_vek_m = dochodkovy_vek_m + posun_doch_vek_m;
 dochodkovy_vek_z = dochodkovy_vek_z + posun_doch_vek_z;
end
%Pridanie udajov pre druhopilieristov
pocet_rokov_2_pilier_m_predpoved = [2020+(1:length(pocet_rokov_2_pilier_m_predpoved))', pocet_rokov_2_pilier_m_predpoved', pocet_m'];
pocet_rokov_2_pilier_z_predpoved = [2020+(1:length(pocet_rokov_2_pilier_z_predpoved))', pocet_rokov_2_pilier_z_predpoved', pocet_z'];

%Grafy poctu sporitelov v 2. pilier pre muzov i zeny
% figure
% scatter(pocet_rokov_2_pilier_m_predpoved(1:end,1),pocet_rokov_2_pilier_m_predpoved(1:end,3),'filled', 'b')
% title({'Poèet sporite¾ov (muži) v 2. pilieri v dôchodkovom veku'})
% xlabel('Rok prognózy')
% ylabel('Poèet sporite¾ov')
% legend({'sporitelia\_v\_2\_pilieri_m (muži)'},'Location','northwest')
 
% figure
% scatter(pocet_rokov_2_pilier_z_predpoved(1:end,1),pocet_rokov_2_pilier_z_predpoved(1:end,3),'filled', 'm')
% title({'Poèet sporite¾ov (ženy) v 2. pilieri v dôchodkovom veku'})
% xlabel('Rok prognózy')
% ylabel('Poèet sporite¾ov')
% legend({'sporitelia\_v\_2\_pilieri_z (ženy)'},'Location','northwest')

%Pocet rokov prace v 1. pilieri
for i=1:length(pocet_rokov_2_pilier_m_predpoved)
 pocet_rokov_2_pilier_m_predpoved(i,4) = dv_m(i) - vstup_na_trh_prace; %dochodkovy_vek - vstup na trh prace -  predpoklad, ze vsetci vstupili a budu vstupovat na trh prace vo veku vstup_na_trh_prace
 pocet_rokov_2_pilier_z_predpoved(i,4) = dv_z(i) - vstup_na_trh_prace;
end

%Vypocet pomeru kratenia dochodku z 1. piliera
%muzi
for i=1:length(pocet_rokov_2_pilier_m_predpoved)
 suma = 0;
 w = 0;
 pocet_rokov_2_pilier_m_predpoved(i,5) = 0;
 po = pocet_rokov_2_pilier_m_predpoved(i,1) - pocet_rokov_2_pilier_m_predpoved(i,2); 
 %od ktoreho roku treba ratat pomery kratenia / dava sa +1 lebo do pocet_rokov_2_pilier_m_predpoved(i,1)+1 sa krati a v pocet_rokov_2_pilier_m_predpoved(i,1) ide do dochodku 
 if (po <= 2012 & pocet_rokov_2_pilier_m_predpoved(i,1) >= 2012)
  suma = suma + (2012 - po);
  pocet_rokov_2_pilier_m_predpoved(i,5) = pocet_rokov_2_pilier_m_predpoved(i,5) + kratenie_dochodku(1,1) * (2012 - po) / pocet_rokov_2_pilier_m_predpoved(i,4);
 end
 for j=2016:1:2024  %lebo od 2016 do 2024 sa koeficienty kratenia lisia rok co rok
  if (po + suma <= j & pocet_rokov_2_pilier_m_predpoved(i,1) > j)
   if (j == 2016)
    pocet_rokov_2_pilier_m_predpoved(i,5) = pocet_rokov_2_pilier_m_predpoved(i,5) + kratenie_dochodku(j-2014,1) * (j - (po + suma - 1)) / pocet_rokov_2_pilier_m_predpoved(i,4);
    suma = suma + (j - (po + suma));
   else 
    pocet_rokov_2_pilier_m_predpoved(i,5) = pocet_rokov_2_pilier_m_predpoved(i,5) + kratenie_dochodku(j-2014,1) * (j - (po + suma)) / pocet_rokov_2_pilier_m_predpoved(i,4);
    suma = suma + (j - (po + suma));   
   end
  end
 end
 if (po + suma >= 2024 & pocet_rokov_2_pilier_m_predpoved(i,1) > 2025)
  pocet_rokov_2_pilier_m_predpoved(i,5) = pocet_rokov_2_pilier_m_predpoved(i,5) + kratenie_dochodku(10,1) * (pocet_rokov_2_pilier_m_predpoved(i,1) - 1 - (po + suma)) / pocet_rokov_2_pilier_m_predpoved(i,4);
 end
end
%zeny
for i=1:length(pocet_rokov_2_pilier_z_predpoved)
 suma = 0;
 w = 0;
 pocet_rokov_2_pilier_z_predpoved(i,5) = 0;
 po = pocet_rokov_2_pilier_z_predpoved(i,1) - pocet_rokov_2_pilier_z_predpoved(i,2);  
 %od ktoreho roku treba ratat pomery kratenia / dava sa +1 lebo do pocet_rokov_2_pilier_m_predpoved(i,1)+1 sa krati a v pocet_rokov_2_pilier_m_predpoved(i,1) ide do dochodku 
 if (po <= 2012 & pocet_rokov_2_pilier_z_predpoved(i,1) >= 2012)
  suma = suma + (2012 - po);
  pocet_rokov_2_pilier_z_predpoved(i,5) = pocet_rokov_2_pilier_z_predpoved(i,5) + kratenie_dochodku(1,1) * (2012 - po) / pocet_rokov_2_pilier_z_predpoved(i,4);
 end
 for j=2016:1:2024 %lebo od 2016 do 2024 sa koeficienty kratenia lisia rok co rok
  if (po + suma <= j & pocet_rokov_2_pilier_z_predpoved(i,1) > j)
   if (j == 2016)
    pocet_rokov_2_pilier_z_predpoved(i,5) = pocet_rokov_2_pilier_z_predpoved(i,5) + kratenie_dochodku(j-2014,1)*(j - (po + suma - 1)) / pocet_rokov_2_pilier_z_predpoved(i,4);
    suma = suma + (j - (po + suma));
   else 
    pocet_rokov_2_pilier_z_predpoved(i,5) = pocet_rokov_2_pilier_z_predpoved(i,5) + kratenie_dochodku(j-2014,1)*(j - (po + suma)) / pocet_rokov_2_pilier_z_predpoved(i,4);
    suma = suma + (j - (po + suma));   
   end
  end
 end
 if (po + suma >= 2024 & pocet_rokov_2_pilier_z_predpoved(i,1) > 2025)
  pocet_rokov_2_pilier_z_predpoved(i,5) = pocet_rokov_2_pilier_z_predpoved(i,5) + kratenie_dochodku(10,1) * (pocet_rokov_2_pilier_z_predpoved(i,1) - 1 - (po + suma)) / pocet_rokov_2_pilier_z_predpoved(i,4);
 end
end
% Graf pomerov kratenia dochodkov z 1. piliera pre druhopilieristov (muzi - zeny)
% figure
% scatter(pocet_rokov_2_pilier_m_predpoved(1:end,1),pocet_rokov_2_pilier_m_predpoved(1:end,5),'filled', 'g')
% title({'Pomer krátenia dôchodkov z 1. piliera pre sporite¾ov (muži - ženy) v 2. pilieri v dôchodkovom veku'})
% xlabel('Rok prognózy')
% ylabel('Pomer krátenia')
% legend({'sporitelia\_v\_2\_pilieri_{m/z} (muži - ženy)'},'Location','northwest')

% % % % % % % VYPOCET udajov pre prvopilieristov % % % % % % % 
%Vysvetlenie jednotlivych stlpcov vo vektore pocet_rokov_len1_pilier_m_predpoved / pocet_rokov_len1_pilier_z_predpoved poporadi:
%1 - rok v ktorom dovrsia doch_vek 
%2 - pocet poistencov len v 1.pilieri v doch_veku
%3 - pocet rokov v 1. pilieri

%Vypocet poctu poistencov len v 1. pilieri v dochodkovom veku
for j=1:predpoved
 pocet_rokov_len1_pilier_m_predpoved(j) = sp1_m(dv_m(j)-vstup_na_trh_prace+1,j) - pocet_rokov_2_pilier_m_predpoved(j,3);
 pocet_rokov_len1_pilier_z_predpoved(j) = sp1_z(dv_z(j)-vstup_na_trh_prace+1,j) - pocet_rokov_2_pilier_z_predpoved(j,3);
end
pocet_rokov_len1_pilier_m_predpoved = [2020+(1:length(pocet_rokov_len1_pilier_m_predpoved)); pocet_rokov_len1_pilier_m_predpoved; dv_m-vstup_na_trh_prace]';
pocet_rokov_len1_pilier_z_predpoved = [2020+(1:length(pocet_rokov_len1_pilier_z_predpoved)); pocet_rokov_len1_pilier_z_predpoved; dv_z-vstup_na_trh_prace]';

%Graf pocet rokov v 1. pilieri
% figure
% scatter(pocet_rokov_2_pilier_m_predpoved(1:end,1), pocet_rokov_len1_pilier_m_predpoved(1:end,3),'filled', 'g')
% ylim([40 75]);
% hold on 
% scatter(pocet_rokov_2_pilier_m_predpoved(1:end,1), dv_m(1:end),'filled', 'r')
% title({'Poèet rokov v 1. pilieri a dôchodkový vek (muži - ženy)'})
% xlabel('Rok prognózy')
% ylabel('Poèet rokov v 1. pilieri / Dôchodkový vek')
% legend({'pocet\_rokov\_v\_1\_pilieri_{m/z}','dochodkovy\_vek_{m/z}'},'Location','northwest')
% hold off
end