%% första korr. 

hbar = 1.0545718e-34; % Plancks reducerade konstant Js
m_e = 9.10938356e-31; % Elektronens massa kg
c = 3e8;              % Ljusets hastighet m/s
a = 1e-9;             %bredden på lådan m

%Diskretisering
N = 50;  % Antal diskreta punkter
delta=1/(N+1);
x = linspace(-1/2+delta, 1/2-delta, N);  % Diskret position

% Startvärde för optimering
start_value = 28;  % Startvärde för f

%fminsearch (optimeringsverktyg) tillsammans med error(felfunktion)
error_function = @(f) calculate_error(f, x, delta, hbar, m_e, a);
options = optimset('Display', 'iter', 'TolFun', 1e-9, 'TolX', 1e-9);  % Optimeringsinställningar
optimal_f = fminsearch(error_function, start_value, options);

%Utskrivet optimalt f-värde
fprintf('Optimalt f: %.9f\n', optimal_f);

% Funktion för att beräkna felet
function error = calculate_error(f, x, delta, hbar, m_e, a)
  
%Dimensionslös potenital
   vk = (f * x).^2 / 2;

 %Tridiagonala matrisen
 main_diag = 1./(delta.^2) ;  % Huvuddiagonal
 off_diag = -1./(2.*delta^2) ;  % Super- och subdiagonal
 main_diag1 = main_diag + vk;  % Justerar huvuddiagonalen med potentialen
 H = diag(main_diag1.*ones(1,length(x))) + diag(off_diag.*ones(1,length(x)-1),1) + diag(off_diag.*ones(1,length(x)-1),-1);

 % Egenvärden och egenvektorer
  [wavefcn, eigenvalue] = eig(H);
  eigenvalue_num1 = eigenvalue(1,1)  % Grundtillståndets numeriska egenvärde

 % Grundtillståndets energi utan störningsteori
   energi_ground = (pi.^2) / 2;  
  
  % Störningsteori 1'a ordn.
   corr = (hbar.^2 * f.^2 * (pi - 6)) / (24 * m_e * pi.^2) * (m_e * a.^2) / (hbar.^2);
   energi_corr = energi_ground + corr  % Korrigerad grundtillståndsenergi mha störningsteori
  

   % Beräkning utav felet mellan korrigerat och numeriskt egenvärde
   error = abs(energi_corr - eigenvalue_num1);  % Absolutvärdet av felet

end
