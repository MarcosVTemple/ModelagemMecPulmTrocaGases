# Marcos Temple
# 25/05/2021
close all
clear all
clc
# Mecânica pulmonar

# 1- DEFINIÇÕES INICIAIS

  # Parâmetros
  f = 0.25;                #[rad/s]
  N = 60000;
  dt = 0.001;                        #[s]
  t = 0:dt:((N-1)*dt);              #[s]
  Pfis = 1.5;                      #Amplitude de Pw
  phi = 0;                       #Fase
  
  Q_O2 = 0.004083;                        # vazão do gás O2 [L/s] 200-300 ml/min
  Q_CO2 = 1.416e-5*Q_O2;                  # vazão do gás CO2 [L/s]

  RK_Q=0;  # inicial
  RK_dQ=0;

  x0=[RK_Q;RK_dQ];
  x=zeros(2,N);
  x(:,1)=x0;
  y=zeros(5,N);
  u_plot = zeros(2,N);
# 2- INTEGRAÇÃO NUMÉRICA - Runge-Kutta 4ªOrdem
  for i=1:(N-1)
    
    #f = f+0.002;
    phi = 2*pi*f*t(i);
    u = entrada_mecanica_pulmonar_RK(f, phi, Pfis);         #dPw
    u_plot(:,i) = u; 
    
    t_umeio=(t(i)+(t(i)+dt))/2;
    phi_umeio = phi+2*pi*f*(dt/2);
    u_meio=entrada_mecanica_pulmonar_RK(f, phi_umeio, Pfis);
    
    t_udt = t(i)+dt;
    phi_udt = phi+2*pi*f*dt;
    u_dt = entrada_mecanica_pulmonar_RK(f, phi_udt, Pfis);
    
    k1 = derivada_mecanica_pulmonar_RK(x(:,i),t(i),u(2)); 
    k2 = derivada_mecanica_pulmonar_RK(x(:,i)+k1*dt/2,t_umeio,u_meio(2));
    k3 = derivada_mecanica_pulmonar_RK(x(:,i)+k2*dt/2,t_umeio,u_meio(2));
    k4 = derivada_mecanica_pulmonar_RK(x(:,i)+k3*dt, (t(i)+dt),u_dt(2));

     % Calcula estado em t+dt
    x(:,i+1) = x(:,i) + dt*(k1+2*k2+2*k3+k4)/6;

     % Calcula saida em t+dt
    y(:,i+1) = saida_mecanica_pulmonar_RK(x(:,i+1),t(:,i+1));
    
    
    phi = phi_udt;
    #fprintf("\n fase: %.8f",phi)
    

  end #for rk

  
# 3- PLOTAGEM DOS RESULTADOS
Q_RK = x(1,:);
Q_SAIDA = y(1,:);
QA_SAIDA = y(2,:);
diff_Q_QA = y(3,:);

Pw = u_plot(1,:)/1.36;  # [mmHg]
dPw = u_plot(2,:)/1.36;
Paw = y(4,:)/1.36;
PA = y(5,:)/1.36;
#t=t(:);

figure(1)
subplot(2,1,1)
plot(t, Q_SAIDA,'r', 'linewidth', 3);
  set(gca, 'FontName', 'Times New Roman', 'FontSize', 10) 
xlabel("Tempo [s]");
ylabel("Fluxo de ar [L/s]");
legend("Q");
title("Ventilação Total (Q) e Alveolar (QA) - Aumento da Resistência Periférica");
#title("Ventilação Total (Q) e Alveolar (QA) - Condições Normais");

figure(1)
subplot(2,1,2)
plot(t, QA_SAIDA,'k', 'linewidth', 3);
  set(gca, 'FontName', 'Times New Roman', 'FontSize', 10) 
xlabel("Tempo [s]");
ylabel("Fluxo de ar [L/s]");
legend("QA");
print -djpg Ventilacao_aum.jpg;

# sum(Q_RK(5200:7200))/1000
V = trapz(Q_SAIDA(7300:9300));
VolumeMinutoCalc = V*f;
fprintf("\nVolume corrente pulmonar: %.3f",V);
fprintf("\nVolume minuto calculado: %.3f",VolumeMinutoCalc);


fprintf("\nAmplitude da pressão alveolar: %.3f",max(PA));

figure(2)
subplot(3,1,1)
plot(t, Pw,'r', 'linewidth', 3);
  set(gca, 'FontName', 'Times New Roman', 'FontSize', 20) 
xlabel("Tempo [s]");
ylabel("Pressão [mmHg]");
legend("Pw");
#title("RK");

figure(2)
subplot(3,1,2)
plot(t, Paw,'b', 'linewidth', 3);
  set(gca, 'FontName', 'Times New Roman', 'FontSize', 20) 
xlabel("Tempo [s]");
ylabel("Pressão [mmHg]");
legend("Paw");

figure(2)
subplot(3,1,3)
plot(t, PA, 'g', 'linewidth', 3);
  set(gca, 'FontName', 'Times New Roman', 'FontSize', 20) 
xlabel("Tempo [s]");
ylabel("Pressão [mmHg]");
legend("PA");

#######################
figure(3)
plot(t, Q_SAIDA,'r', 'linewidth', 3, t, QA_SAIDA,'k', 'linewidth', 3);
  set(gca, 'FontName', 'Times New Roman', 'FontSize', 20) 
xlabel("Tempo [s]");
ylabel("Fluxo de ar [L/s]");
legend("Q", "QA");
#title("Ventilação Total e Alveolar");


figure(4)
plot(t, Pw,'r', 'linewidth', 3, t, Paw,'b', 'linewidth', 3, t, PA, 'g', 'linewidth', 3);
  set(gca, 'FontName', 'Times New Roman', 'FontSize', 10) 
xlabel("Tempo [s]");
ylabel("Pressão [mmHg]");
legend("Pw","Paw","PA");
title("Pressão na Cavidade Torácica (Pw), Vias Aéreas Centrais (Paw) e Alvéolos (PA) - Aumento da Resistência Periférica");
print -djpg Pressao_aum.jpg;
