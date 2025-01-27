%%ana./num. N=10

clear;clc;clf

%given data
N=10; %diskreta pkter
delta=1/(N+1); %diskret steglängd
vk=0; %pot.energi

% Bygg tridiagonala matrisen
main_diag = 1./(delta.^2) * ones(N,1); % Huvuddiagonal
off_diag = -1./(2.*delta^2) * ones(N-1,1); % Super- och subdiagonalen

% Justera huvuddiagonalen för att inkludera potentialen
main_diag1= main_diag + vk;

% Konstruera tridiagonala matrisen
H = diag(main_diag1) + diag(off_diag, 1) + diag(off_diag, -1);

%Hitta egenvärden och egenvektorer
[wavefcn, eigenvalue] = eig(H);
disp(eigenvalue); 

%Beräkna de analytiska egenvärdena för att kunna jämföra dem med de
%numeriska i en plot
n=1:N; %vektor med enerinivåindex där n=1,2,3,....,N
ana_eigenvalue=(pi^2*n.^2)/2; %se analytisk lösn. för framtagning av ekv.

disp('Analytiska egenvärden för de 10 första energinivåerna');
disp(ana_eigenvalue);

%Normalisera:
norm_factors=zeros(1, N); 

for n=1:N
    psi_n=wavefcn(:,n);
    norm_factors(n)=sqrt(sum(abs(psi_n).^2)*delta); %normaliseringsfaktor
    wavefcn(:,n)=wavefcn(:,n)./norm_factors(n); %vågfkn normaliserad
    % 
    % %Definerar gränsvärdena =0
    % wavefcn(1,n)=0; %x=-1/2
    % wavefcn(N,n)=0; %x=1/2
end 

%Normeringsfaktorer
disp('Normeringsfaktorer för varje kolonn i wavefcn');
disp(norm_factors);

%plotta och jämför de analytiska och numeriska funktionerna:
for n=1:N
    figure(n);
    hold on
    x=linspace(-1/2+delta, 1/2-delta, N); %definera intervallet för x-vektorn
    plot(x, wavefcn(:,n), 'm', 'LineWidth', 1.5, 'DisplayName', ['\psi_' num2str(n)]);
    
     x=linspace(-1/2, 1/2);
    if mod(n,2)==1 %udda n
        psi_ana=sqrt(2)*cos(n*pi*x);
    else %jämn n
        psi_ana=sqrt(2)*sin(n*pi*x);
    end 

 %analytisk funktion:
 plot(x, psi_ana, 'c', 'LineWidth', 2.5, 'DisplayName', ['Analytisk \psi_ana_' num2str(n)]);
 
 xlabel('Position x')
 ylabel(['\psi_' num2str(n) '(x)']);
 title(['Numerisk och analytisk lösning för n = ' num2str(n)]);
 legend show;

 hold off;

end 

%%N=50
format long g

%given data:
N=50; %diskreta punkter,antalet kan ändras 
delta=1/(N+1); %ska det vara minus här också??
vk=0; %pot.energi

% Bygg tridiagonala matrisen
main_diag = 1./(delta.^2) * ones(N,1); % Huvuddiagonal
off_diag = -1./(2.*delta^2) * ones(N-1,1); % Super- och subdiagonalen

% Justera huvuddiagonalen för att inkludera potentialen
main_diag1= main_diag + vk;

% Konstruera tridiagonala matrisen
H = diag(main_diag1) + diag(off_diag, 1) + diag(off_diag, -1);

%Hitta egenvärden och egenvektorer
[wavefcn, eigenvalue] = eig(H);
eigenvalues=diag(eigenvalue);
disp('egenvärden num.');
disp(eigenvalues); 

n=1:N; %vektor med enerinivåindex där n=1,2,3,....,N
ana_eigenvalues=(pi^2*n.^2)/2; %se analytisk lösn. för framtagning av ekv.
ana_eigenvalues=ana_eigenvalues.'; %transporterar eigenvektorn
disp('egenvärden ana. för 10 första energinivåerna')
disp(ana_eigenvalues);


%%a,fmin,fmax
clear; clc; clf;

% Konstanter
hbar = 1.0545718e-34; % Plancks reducerande konstant J*s
m_e = 9.10938356e-31; % Elektronens massa kg
c =3e8;              % Ljusets hastighet m/s

% Givna våglängder 
L_max=790e-9;    %ger ett fmin
L_mid=593e-9;    %ger ett fmid (medelvärdet av fmin och fmax)
L_min=431e-9;    %ger ett fmax
L_target=[L_min, L_max, L_mid];

% Startvärden för optimering
a_start = 1e-9;        % Startvärde för a
fmin_start = 10;       % Startvärde för fmin
fmax_start = 35;       % Startvärde för fmax
start_values = [a_start, fmin_start, fmax_start];


%fminsearch (optimeringsverktyg) tillsammans med error(felfunktion)
options = optimset('Display', 'iter', 'TolFun', 1e-9, 'TolX', 1e-9);
optimal_values = fminsearch(@(parameters) calculate_error(parameters, L_target, hbar, m_e, c), start_values, options);

% Framtagning av optimala värden
[optimal_a, optimal_fmin, optimal_fmax]= deal(optimal_values(1), optimal_values(2), optimal_values(3));

% Utskrivna optimala värden 
fprintf('Optimalt a (nm): %.9f\n', optimal_a * 1e9);
fprintf('Optimalt fmin: %.9f\n', optimal_fmin);
fprintf('Optimalt fmax: %.9f\n', optimal_fmax);
fprintf('Optimalt fmedel: %.9f\n', (optimal_fmin + optimal_fmax) / 2);

%Beräkning av tot.felet mellan numerisk och givna våglängder
function error = calculate_error(parameters, L_target, hbar, m_e, c)
   
    a = parameters(1);
    fmin = parameters(2);
    fmax = parameters(3);
    fmid = (fmax + fmin) / 2; 

    % Diskretisering
    N = 100; 
    x = linspace(-1/2, 1/2, N); 
    delta = x(2) - x(1);    

    %Dimensionslös potential tillhörande varje f
     vk_max = (fmax * x).^2 / 2;
     vk_min = (fmin * x).^2 / 2;
     vk_mid = (fmid * x).^2 / 2;
 
    %Tridiagonala Hamiltion-matrisen tillhörande varje f
     main_diag = 1./(delta.^2);  % Huvuddiagonal
     off_diag = -1./(2.*delta^2) * ones(N-1,1);  % Super- och subdiagonal
     main_diag1 = main_diag + vk_max;  % Justerar huvuddiagonalen med potentialen
     H_max = diag(main_diag1) + diag(off_diag, 1) + diag(off_diag, -1);

     main_diag = 1./(delta.^2);  % Huvuddiagonal
     off_diag = -1./(2.*delta^2) * ones(N-1,1);  % Super- och subdiagonal
     main_diag1 = main_diag + vk_min;  % Justerar huvuddiagonalen med potentialen
     H_min = diag(main_diag1) + diag(off_diag, 1) + diag(off_diag, -1);

     main_diag = 1./(delta.^2);  % Huvuddiagonal
     off_diag = -1./(2.*delta^2) * ones(N-1,1);  % Super- och subdiagonal
     main_diag1 = main_diag + vk_mid;  % Justerar huvuddiagonalen med potentialen
     H_mid = diag(main_diag1) + diag(off_diag, 1) + diag(off_diag, -1);

     %Egenvärden (energier) tillhörande varje f
     [~, max_eigenvalues] = eig(H_max);
     max_energies = sort(diag(max_eigenvalues)) * (hbar^2 / (m_e * a^2));

     [~, min_eigenvalues] = eig(H_min);
     min_energies = sort(diag(min_eigenvalues)) * (hbar^2 / (m_e * a^2));

     [~, mid_eigenvalues] = eig(H_mid);
     mid_energies = sort(diag(mid_eigenvalues)) * (hbar^2 / (m_e * a^2));

     %Beräkna energidiffereansen för varje f
     deltaE12_max = abs( max_energies(2) -  max_energies(1));

     deltaE12_min = abs(min_energies(2) - min_energies(1));

     deltaE12_mid = abs(mid_energies(2) - mid_energies(1));

     %Beräkna numerisk våglängd
     L_num_mid = (hbar * 2 * pi * c) / deltaE12_mid;

     L_num_max = (hbar * 2 * pi * c) / deltaE12_min;

     L_num_min = (hbar * 2 * pi * c) / deltaE12_max;

     L_num = [L_num_min, L_num_max, L_num_mid];

     % Beräkna kvadrerat fel mellan beräknade och givna våglängder
    error = sum((L_target - L_num).^2);
end

