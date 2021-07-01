tic

grid = 1e-3;        % Gitterweite [m]
l = 1000;           % Länge der Schiene [m]

mu_0 = 1.25664e-6;  % magnetische Feldkonstante [Vs/Am]
mu_r = 50;          % relative Permeabilität
rho = 2e-7;         % spezifischer Widerstand [Ωm] Werte aus: DIN EN 50641, F. Kießling, R. Puschmann, A. Schmieder, Fahrleitungen elektrischer Bahnen: Planung, Berechnung, Ausführung, Betrieb , 3. Aufl. Erlangen: Publicis Publishing,

f = 50;             % Frequenz [Hz]
I_0 = 100;          % Schienenstrom [A]

profile = [74, 75, 75, 75, 75, 75, 75, 75, 75, 74, 73, 71, 58, 50, 47, 44, 41, 39, 36, 33, 31, 28, 25, 22, 20, 17, 15, 14, 13, 13, 12, 12, 12, 12, 11, 11, 11, 11, 11, 10, 10, 10, 10, 10, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 10, 10, 10, 10, 10, 10, 11, 11, 11, 11, 11, 12, 12, 12, 13, 13, 14, 15, 17, 20, 22, 24, 27, 30, 32, 35, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 37, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 35, 35, 34, 34, 33, 32, 30, 28, 25, 18];
n = 2 * sum(profile);

% Bestimmen der Positionen der Leiter
position = zeros(n, 2); % Spalte 1: y-Koordinate, Spalte 2: z-Koordinate
ip = 1;
for i1 = 1:length(profile)
    for i2 = 1:profile(i1)
        position(ip, 1) = (i2 - 0.5) * grid;
        position(ip, 2) = (i1 - 0.5) * grid;
        position(ip + 1, 1) = - (i2 - 0.5) * grid;
        position(ip + 1, 2) = (i1 - 0.5) * grid;
        ip = ip + 2;
    end
end
Zii = rho * l / grid ^ 2 + 1j * f * mu_0 * mu_r * l * (log(l / grid) - 0.2235 * log(2 * grid / l) + 0.5);

A = [[diag(Zii * ones(n, 1)) (-ones(n, 1))];[ones(1, n) 0]];

for i1 = 1:(n - 1)
   for i2 = (i1 + 1):n
       dij = sqrt((position(i1, 1) - position(i2, 1))^2 + (position(i1, 2) - position(i2, 2))^2);
       Xij = 1j * f * mu_0 * mu_r * l * (asinh(l / dij) - sqrt(1 + (dij / l)^2) + dij / l);
       A(i1, i2) = Xij;
       A(i2, i1) = Xij;
   end
end
b = [zeros(n, 1); I_0];

x = linsolve(A, b);

toc