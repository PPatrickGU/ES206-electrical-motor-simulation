clear all;

%Valeur des constants
time = 7e-3;
E = 400;
Wv = 220;
Pp = 4;
Ld = 1.15e-3; %H
Lq = 3.31e-3; %H
J = 800e-6;
Phim = 200e-3;
Rs = 0.18;
dt = 1e-7;
tp = 100e-6;
Gp = 5000;
Gi = 0;
RCmax = 50; %(50N.m)Couple maximum
pas = time/dt;

i = 1:1:pas;
t = i * dt;

%Pseudo-Initialisation des variables
RVd = zeros(1,pas); 
RVq = zeros(1,pas); 
RV1 = zeros(1,pas); 
RV2 = zeros(1,pas); 
RV3 = zeros(1,pas);
RVn = zeros(1,pas);
RVo1 = zeros(1,pas); 
RVo2 = zeros(1,pas); 
RVo3 = zeros(1,pas); 
Vo1 = zeros(1,pas); 
Vo2 = zeros(1,pas);
Vo3 = zeros(1,pas); 
Vd = zeros(1,pas); 
Vq = zeros(1,pas);
Phid = zeros(1, pas);
dPhid = zeros(1, pas);
Phiq = zeros(1, pas);
dPhiq = zeros(1, pas);
Id = zeros(1, pas);
Iq = zeros(1, pas);
I1 = zeros(1, pas);
I2 = zeros(1, pas);
I3 = zeros(1, pas);
Itot = zeros(1, pas);
MI1 = zeros(1, pas);
MI2 = zeros(1, pas);
MI3 = zeros(1, pas);
C = zeros(1, pas);
Wm = zeros(1, pas);
dWm = zeros(1, pas);
We = zeros(1, pas);
A = zeros(1, pas);
MA = zeros(1, pas);
tint = zeros(1, pas);
tint(pas) = 0;
tr1 = zeros(1, pas);
tr2 = zeros(1, pas);
tr3 = zeros(1, pas);
tf1 = zeros(1, pas);
tf2 = zeros(1, pas);
tf3 = zeros(1, pas);
MId = zeros(1, pas);
MIq = zeros(1, pas);
MWe = zeros(1, pas);
MPhid = zeros(1, pas);
MPhiq = zeros(1, pas);
Srd = zeros(1, pas);
Srq = zeros(1, pas);
EPhid = zeros(1, pas);
EPhiq = zeros(1, pas);
RPhid = zeros(1, pas); 
RPhiq = zeros(1, pas);
RC = zeros(1, pas);

% etablir la couple RC
seuil1 = (500e-6)/dt;
seuil2 = (1500e-6)/dt;
for j = 1:seuil1
    RC(j) = 0;
end
for j = fix((seuil1+1)):seuil2
    RC(j) = RC(j-1) + RCmax/(1000e-6/dt); % to determine
end
for j = fix((seuil2+1)):pas
    RC(j) = RCmax;
end
for i = 2:1:(pas-1)
    tint(i) = mod(i * dt, tp) - dt;
    if (-0.5< (tint(i)/dt) < 0.5) 
      continue
    else
      RC(i+1) = RC(i);
    end
end

% Initialisation des Boucles
Vo1(1) = -E;
Vo2(1) = -E;
Vo3(1) = -E;
Phid(1) = Phim;
Phiq(1) = 0;
A(1) = 0;
MA(1) = 0;
Wm(1) = 0;
Id(1) = (Phid(1) - Phim)/Ld;
Iq(1) =  Phiq(1)/Lq;
C(1) = Pp * (Phid(1) * Iq(1) - Phiq(1) * Id(1));
Vd(1) = (2/3)^(1/2) * (cos(A(1))* Vo1(1) + cos(A(1) -2*pi/3)*Vo2(1) + cos(A(1) - 4*pi/3)*Vo3(1));
Vq(1) = (2/3)^(1/2) * (-sin(A(1))* Vo1(1) - sin(A(1) -2*pi/3)*Vo2(1) - sin(A(1) - 4*pi/3)*Vo3(1));
RV1(1) = (2/3)^(1/2) * (cos(MA(1))*RVd(1) - sin(MA(1))*RVq(1));
RV2(1) = (2/3)^(1/2) * (cos(MA(1) -2/3*pi)*RVd(1) - sin(MA(1)-2/3*pi)*RVq(1));
RV3(1) = (2/3)^(1/2) * (cos( MA(1) -4/3*pi)*RVd(1) - sin(MA(1) -4/3*pi)*RVq(1));
RVMAX = max([RV1(1),RV2(1),RV3(1)]);
RVMIN = min([RV1(1),RV2(1),RV3(1)]);
RVn(1) = -(RVMAX +RVMIN)/2;
RVo1(1) = RV1(1) + RVn(1);
RVo2(1) = RV2(1) + RVn(1); 
RVo3(1) = RV3(1) + RVn(1);
tr1(1) = tp * (1 - RVo1(1)/E)/4; 
tf1(1) = tp * (3 + RVo1(1)/E)/4;
tr2(1) = tp * (1 - RVo2(1)/E)/4; 
tf2(1) = tp * (3 + RVo2(1)/E)/4;
tr3(1) = tp * (1 - RVo3(1)/E)/4; 
tf3(1) = tp * (3 + RVo3(1)/E)/4;
MI1(1) = I1(1);
MI2(1) = I2(1);
MI3(1) = I3(1);
MId(1) = (2/3)^(1/2)*(cos(MA(1))*MI1(1) + cos(MA(1) - 2*pi/3)*MI2(1) + cos(MA(1) - 4*pi/3)*MI3(1));  
MIq(1) = (2/3)^(1/2)*(-sin(MA(1))*MI1(1) - sin(MA(1) - 2*pi/3)*MI2(1) - sin(MA(1) - 4*pi/3)*MI3(1));
MWe(1) = 0;
MPhid(1) = MId(1)*Ld + Phim;
MPhiq(1) = MIq(1)*Lq;
RPhid(1) = Phim;
Rphiq(1) = (Lq*RC(1))/(Pp*Phim);
EPhid(1) = RPhid(1) - MPhid(1);
EPhiq(1) = RPhiq(1) - MPhiq(1);
Srd(1) = EPhid(1)*Gp;
Srq(1) = EPhiq(1)*Gp;
RVd(1) = Srd(1) - MPhiq(1)*MWe(1);
RVq(1) = Srq(1) + MPhid(1)*MWe(1);



tint(pas) = 0;

% boucle
for i = 1:1:(pas-1)

  % Calcul de Id et Iq (Courbes 3) 
    dPhid(i) = Phiq(i) * We(i) + Vd(i) - Id(i) * Rs;
    dPhiq(i) = Vq(i) - Phid(i) * We(i) - Iq(i) * Rs;
    Phid(i+1) = Phid(i) +  dPhid(i) * dt;
    Phiq(i+1) = Phiq(i) +  dPhiq(i) * dt;
    Id(i+1) = (Phid(i+1) - Phim)/Ld;
    Iq(i+1) = Phiq(i+1)/Lq;
  
    %Grandeur mecanique (Courbe 6)
    C(i+1) = Pp * (Phid(i+1) * Iq(i+1) - Phiq(i+1) * Id(i+1));
    dWm(i) = (C(i) )/J;
    Wm(i+1) = Wm(i) + dWm(i) * dt;
    We(i+1) = Wm(i+1) * Pp;
    A(i+1) = A(i) + We(i) * dt;
    
    % Calcul de I1 I2 I3
    I1(i+1) = (2/3)^(1/2) * (cos(A(i+1))*Id(i+1) - sin(A(i+1))*Iq(i+1));
    I2(i+1) = (2/3)^(1/2) * (cos(A(i+1) - pi*2/3)*Id(i+1) - sin(A(i+1) - pi*2/3)*Iq(i+1));
    I3(i+1) = (2/3)^(1/2) * (cos(A(i+1) - pi*4/3)*Id(i+1) - sin(A(i+1) - pi*4/3)*Iq(i+1));
    
    % Renouvellement de Mesure ou pas  
    tint(i) = mod(i * dt, tp) - dt ;    
    if  (tint(i) > tr1(i)) && (tint(i) < tf1(i)) 
        Vo1(i+1) = E;
    else
        Vo1(i+1) = -E;
    end
    
    if  (tint(i) > tr2(i)) && (tint(i) < tf2(i))
        Vo2(i+1) = E;
    else
        Vo2(i+1) = -E;
    end
    
    if  (tint(i) > tr3(i)) && (tint(i) < tf3(i))
        Vo3(i+1) = E;
    else
        Vo3(i+1) = -E;
    end
    
    if  (-0.5< (tint(i)/dt) < 0.5) 
        MA(i+1) = A(i+1);
        MI1(i+1) = I1(i+1);
        MI2(i+1) = I2(i+1);
        MI3(i+1) = I3(i+1);
        MId(i+1) = (2/3)^(1/2)*(cos(MA(i+1))*MI1(i+1) + cos(MA(i+1) - 2*pi/3)*MI2(i+1) + cos(MA(i+1) - 4*pi/3)*MI3(i+1));  
        MIq(i+1) = (2/3)^(1/2)*(-sin(MA(i+1))*MI1(i+1) - sin(MA(i+1) - 2*pi/3)*MI2(i+1) - sin(MA(i+1) - 4*pi/3)*MI3(i+1));
        MPhid(i+1) = MId(i+1)*Ld + Phim;
        MPhiq(i+1) = MIq(i+1)*Lq;
        MWe(i+1) = (MA(i+1) - MA(i))/tp;

    else
        MA(i+1) = MA(i);
        MI1(i+1) = MI1(i);
        MI2(i+1) = MI2(i);
        MI3(i+1) = MI3(i);
        MId(i+1) = MId(i); 
        MIq(i+1) = MIq(i);
        MPhiq(i+1) = MPhiq(i);
        MPhid(i+1) = MPhid(i);
        MWe(i+1) = MWe(i);
    end
    
    RPhid(i+1) = Phim;
    RPhiq(i+1) = (Lq*RC(i+1))/(Pp*Phim);
    
    EPhid(i+1) = RPhid(i+1) - MPhid(i+1);
    EPhiq(i+1) = RPhiq(i+1) - MPhiq(i+1);
    Srd(i+1) = EPhid(i+1)*Gp;
    Srq(i+1) = EPhiq(i+1)*Gp;
    RVd(i+1) = Srd(i+1) - MPhiq(i+1)*MWe(i+1);
    RVq(i+1) = Srq(i+1) + MPhid(i+1)*MWe(i+1);

    
    % Calcul de RVi (Courbes 8) 
    RV1(i+1) = (2/3)^(1/2) * (cos(MA(i+1))*RVd(i+1) - sin(MA(i+1))*RVq(i+1));
    RV2(i+1) = (2/3)^(1/2) * (cos(MA(i+1) -2/3*pi)*RVd(i+1) - sin(MA(i+1)-2/3*pi)*RVq(i+1));
    RV3(i+1) = (2/3)^(1/2) * (cos( MA(i+1) -4/3*pi)*RVd(i+1) - sin(MA(i+1) -4/3*pi)*RVq(i+1));
 
    RVMAX = max([RV1(i+1),RV2(i+1),RV3(i+1)]);
    RVMIN = min([RV1(i+1),RV2(i+1),RV3(i+1)]);
    RVn(i+1) = -(RVMAX +RVMIN)/2;
    
    RVo1(i+1) = RV1(i+1) + RVn(i+1);
    RVo2(i+1) = RV2(i+1) + RVn(i+1); 
    RVo3(i+1) = RV3(i+1) + RVn(i+1);
    
    %Renouvellement de tr ou pas
    if (-0.5< (tint(i)/dt) < 0.5) 
        tr1(i+1) = tp * (1 - RVo1(i+1)/E)/4; 
        tf1(i+1) = tp * (3 + RVo1(i+1)/E)/4;
        tr2(i+1) = tp * (1 - RVo2(i+1)/E)/4; 
        tf2(i+1) = tp * (3 + RVo2(i+1)/E)/4;
        tr3(i+1) = tp * (1 - RVo3(i+1)/E)/4; 
        tf3(i+1) = tp * (3 + RVo3(i+1)/E)/4;
    else
        tr1(i+1) = tr1(i);
        tf1(i+1) = tf1(i);
        tr2(i+1) = tr2(i);
        tf2(i+1) = tf2(i);
        tr3(i+1) = tr3(i);
        tf3(i+1) = tf3(i);
    end
    
    % Calcul de Vd et Vq (Courbes 2)    
    Vd(i+1) = (2/3)^(1/2) * (cos(A(i+1))* Vo1(i+1) + cos(A(i+1)-2*pi/3)*Vo2(i+1) + cos(A(i+1) - 4*pi/3)*Vo3(i+1));  
    Vq(i+1) = (2/3)^(1/2) * (-sin(A(i+1))* Vo1(i+1) - sin(A(i+1) -2*pi/3)*Vo2(i+1) - sin(A(i+1) - 4*pi/3)*Vo3(i+1)); 
end


%Courbes



subplot(3,3,1)
plot(t, Vo1+50, t, Vo2, t, Vo3-50,'LineWidth',2);
xlim([0 2*tp]);
hold on;
plot(t,tint*1e6,'--',t,tr1*1e6,'--',t,tf1*1e6,'--','linewidth',2);
title('Tensions triphas¨¦es de l¡¯onduleur')
xlabel('Temps (s)')
ylabel('Tension (volt)')
legend({'Vo1+50V','Vo2','Vo3-50V','tint*10^6','tr1*10^6','tf1*10^6'},'Location','northeast');

subplot(3,3,2)
plot(t, RVd, t, RVq, t, -RVn,'LineWidth',2);
title('Tensions direct et quadratique')
xlabel('Temps (s)')
ylabel('Tension (Volt)')
legend('RVd','RVq', '-RVn')


subplot(3,3,3)
plot(t, Phid, t, Phiq, t, RPhid, t, RPhiq, t, MPhid, t, MPhiq, 'LineWidth',2);
title('Flux direct et quadratique')
xlabel('Temps (s)')
ylabel('Phi')
legend('Phid','Phiq','RPhid','RPhiq','MPhid','MPhiq')


subplot(3,3,4)
plot(t, Id, t, Iq, t, MId, t, MIq, 'LineWidth',2);
title('Courants direct et quadratique')
xlabel('Temps (s)')
ylabel('Courant (A)')
legend('Id','Iq','MId','MIq')


subplot(3,3,5)
plot(t, I1, t, I2, t, I3, t, Itot, t , MI1, t, MI2, t, MI3, 'LineWidth',2);
title('Courants triphas¨¦es')
xlabel('Temps (s)')
ylabel('Courant (A)')
legend('I1','I2','I3 ','I1+I2+I3', 'MI1','MI2','MI3 ')


subplot(3,3,6)
plot(t, C*100, t, A*1000, t, We*10, t, MA*1000, t, RC*100, t, MWe*10,t, RVd.*Id + RVq .* Iq, 'LineWidth',2);
title('Grandeurs m¨¦caniques')
xlabel('Temps (s)')
ylabel('')
legend('C*100','A*1000','We*10','MA*1000','RC*100','MWe*10','P')

% subplot(3,3,6)
% plot(t, C*10, t, A*100, t, We, t, MA*100, t, RC*10, t, MWe, 'LineWidth',2);
% title('Grandeurs m¨¦caniques')
% xlabel('Temps (s)')
% ylabel('')
% legend('C*10','A*100','We','MA*100','RC*10','MWe')

subplot(3,3,7)
plot(t, RVo1, t, RVo2, t, RVo3,'LineWidth',2);
title('Tensions triphas¨¦es')
xlabel('Temps (s)')
ylabel('Tension (volt)')
legend('RVo1','RVo2','RVo3')

subplot(3,3,8)
plot(t, RV1, t, RV2, t, RV3, t,-RVn, 'LineWidth',2);
%xlim([0 10*tp]);
title('RV1, RV2, RV3, -RVn')
xlabel('Temps (s)')
ylabel('Tension (volt)')
legend('RV1', 'RV2', 'RV3', '-RVn')









