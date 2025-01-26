% 2. b
% Konstanter
hbar = 1.0545718e-34;  % Plancks reducerade konstant J·s
m_e = 9.10938356e-31;  % Elektronens massa kg
c = 3e8;               % Ljusets hastighet m/s

% Bredd på lådan
a = 1e-9;  % För enkelhets skull väljs a = 1e-9 m

% Diskretisering
N = 100;
x = linspace(-1/2, 1/2, N);
delta = x(2) - x(1);  % Steglängd

% Valda värden på f
f_values = [10 30 50 100 120 140 200 220 240];

% Skapar en figur och håller den öppen för multipla plottar
figure;
hold on;

% Loopar över olika f-värden
for i = 1:length(f_values)
    f = f_values(i);

    % Dimensionslös potential
    vk = (f * x).^2 / 2;

    % Tridiagonal Hamilton-matris
    main_diag = 1./(delta.^2) * ones(N,1);  % Huvuddiagonal
    off_diag = -1./(2.*delta.^2) * ones(N-1,1);  % Super- och subdiagonal
    main_diag1 = main_diag + vk(:);  % Justerar huvuddiagonalen med potentialen

    % Skapar tridiagonala H-matrisen
    H = diag(main_diag1) + diag(off_diag, 1) + diag(off_diag, -1);

    % Egenvärden och egenvektorer för H
    [wavefcn, eigenvalue] = eig(H);

    % Vågfunktion för grundtillståndet (första egenvektorn)
    ground_state_wavefcn = wavefcn(:, 1);

    % Normalisering av vågfunktionen
    ground_state_wavefcn = ground_state_wavefcn / sqrt(sum(abs(ground_state_wavefcn).^2) * delta);

    % Plottar vågfunktionen för grundtillståndet
    plot(x, ground_state_wavefcn , 'DisplayName', ['f = ', num2str(f)], 'LineWidth', 1.5);
end

%Plot design
xlabel('x');
ylabel('Vågfunktion \psi(x)');
title('Grundtillståndets vågfunktion för olika värden på f');
legend('show');
grid on;
hold off;




