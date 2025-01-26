% Konstanter
hbar = 1.0545718e-34; % Plancks reducerade konstant Js
m_e = 9.10938356e-31; % Elektronens massa kg
c = 3e8;              % Ljusets hastighet m/s
a=1e-9;               %lådans bredd m

%Diskretisering
N = 50;  
delta=1/(N+1);
f=100;
x=linspace(-1/2+delta,1/2-delta,N);

%Dimensionlös potential
vk = (f*x).^2/2; 

%Tridiagonala matrisen
main_diag = 1./(delta.^2) * ones(1,N); % Huvuddiagonal
off_diag = -1./(2.*delta^2) * ones(1,N-1); % Super- och subdiagonalen
main_diag1= main_diag + vk;
H = diag(main_diag1) + diag(off_diag, 1) + diag(off_diag, -1);

%Egenvärden och egenvektorer
 [wavefcn, eigenvalue] = eig(H);
fprintf('numeriska egenvärden vid grundtillstånd(dim.lös)');
eigenvalue_num1=eigenvalue(1:1)

%Grundtillståndet energi, med störningsteori
energi_ground=(pi.^2)/2;
corr=(hbar.^2*f.^2*(pi-6))/(24*m_e*pi.^2)*(m_e*a.^2)./(hbar.^2); %korrigering 1'a ordn. dimensionslös
energi_corr=energi_ground+corr;
fprintf('1a ordn. korrigering energinivå vid grundtillstånd (dim.lös)');
disp(energi_corr);

