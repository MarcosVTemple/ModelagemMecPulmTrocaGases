function y = saida_mecanica_pulmonar_RK(x, t)
  
  # Constantes
  Cs = 0.005;          #[L/cmH2O]
  Rc = 1;           #[cmH2O.s/L]
  Rp = 0.5;          #[cmH2O.s/L]
  Rp = 5;          #[cmH2O.s/L]
  Cl = 0.2;         #[L/cmH2O]
  
  # Variaveis de interesse
  y1 =  x(1);                       #Q
  y2 =  x(1) + Cs*Rc*x(2);     #QA
  y3 =  y1-y2; #Q-QA
  y4 =  -x(1)*(Rc); #Paw
  y5 = -(Rc*y1+Rp*y2); #PA = -Rc*Q - Rp*QA
  y=[y1;y2;y3;y4;y5];
endfunction