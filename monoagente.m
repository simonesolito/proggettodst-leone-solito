clear
close all

%% inizializzazione componenti del sistema

% Parametri orbitali
mu = 3.986e14;                                  % parametro gravitazionale standard (Terra)
a = 7000e3;                                     % lunghezza del semiasse maggiore (7000 km)
e = 0.1;                                        % Orbita ellittica
Ts = 10;                                        % Campionamento ogni 10 secondi
Tmax = 2 * pi * sqrt(a^3/mu);                   % durata dell'orbita

T = 0 : Ts : Tmax;                              % asse dei tempi
Nsamp = length(T);

omega_0 = sqrt(mu/a^3);                         % velocità angolare media del satellite
M = omega_0 * T;                                % anomalia media
m = 500e3;                                      % massa satellite

% calcolo dell'anomalia eccentrica E partendo dalla equazione M = E - e * sin(E)
% Newton-Raphson

E = M; 
for i = 1:10 
    E = E - (E - e*sin(E) - M) / (1 - e*cos(E));
end

% calcolo di posizione, velocità e accelerazione angolare
f = 2 * atan2(sqrt(1 + e) * sin(E/2) , sqrt(1 - e) * cos(E/2));
df = omega_0 * (1 + e * cos(f)).^2/(1 - e^2)^(3/2);
d2f = -2* omega_0^2 * e * sin(f) .* (1 + e * cos(f)).^3/(1 - e^2)^3;

rt = a * (1 - e^2)./(1 + e * cos(f));               % distanza satellite terra 

A1 = zeros(3,3,Nsamp);                          % LA PARTE CONINUA LA POTREI MEMORIZZARE MOMENTANEAMENTE
A2 = A1;                                            % RIDUCENDO LE MATRICI A DUE DIMENSIONI
Acont = zeros(6,6,Nsamp);
Bcont = 1/m * [zeros(3,3) ; eye(3)];

Adisc = Acont;
Bdisc = Bcont;


for i = 1:Nsamp 
    % calcolo delle matrici del sistema continuo
    A1(:,:,i) = [df(i)^2 + 2*mu/rt(i)^3, d2f(i), 0;
                 -d2f(i), df(i)^2 - mu/rt(i)^3, 0;
                 0,     0,      -mu/rt(i)^3         ];

    A2(1,2,i) = 2*df(i);
    A2(2,1,i) = -A2(1,2,i);

    Acont(:,:,i) = [zeros(3,3) , eye(3) ; A1(:,:,i) , A2(:,:,i)];

    % discretizzazione del sistema
    Maux = [Acont(:,:,i) , Bcont ; zeros(3,9)] * Ts;
    Mexp = expm(Maux);

    Adisc(:,:,i) = Mexp(1:6,1:6);
    Bdisc(:,:,i) = Mexp(1:6,7:9);

end

%% Calcolo della traiettoria

Q = eye(6);                             % Peso della precisione
R = eye(3);                             % Peso del consumo

P = zeros(6,6,Nsamp);
P(:,:,Nsamp) = Q;
NiterRiccati = 5;

for i = 1 : NiterRiccati
    for k = Nsamp-1:-1:1
        Ak = Adisc(:,:,k);
        Bk = Bdisc(:,:,k);
        Pk_next = P(:,:,k+1);
            
        % Termini intermedi per l'Eq. 32 (per leggibilità e velocità)
        temp_inv = (R + Bk' * Pk_next * Bk);
            
        % Equazione di Riccati Discreta (Eq. 32 del paper)
        P(:,:,k) = Q + Ak' * Pk_next * Ak - ...
                (Ak' * Pk_next * Bk) * (temp_inv \ (Bk' * Pk_next * Ak));
    end
    % Per la periodicità: l'inizio del giro precedente diventa la fine del nuovo
    P(:,:,Nsamp) = P(:,:,1);
end

flag = zeros(1,3);
Kk = zeros(3,6,Nsamp);

for k = 1:Nsamp-1
    Ak = Adisc(:,:,k);
    Bk = Bdisc(:,:,k);
    Pk = P(:,:,k);
    Pk_next = P(:,:,k+1);
    
    
    Kk(:,:,k) = (R + Bk' * Pk_next * Bk) \ (Bk' * Pk_next * Ak);

    % Analisi di stabiliità
    Acl = Ak - Bk * Kk(:,:,k);

    [~, flag(1)] = chol(Pk);
    Q = Pk - transpose(Acl)*Pk_next*Acl;
    [~, flag(2)] = chol(Q);

   
    if flag(1) > 0 || flag(2) > 0
        fprintf('La matrice del campione numero %d non è definita positiva\n', k);
        flag(3)=1;
    end
end

if flag(3) == 0
    disp('Il sistema è asintoticamente stabile');
end


Kk(:,:,Nsamp) = Kk(:,:,1); % Continuità periodica



%% Definizione del sistema
zeta = zeros(6, Nsamp);
zeta(:,1) = [100; -50; 20; 0; 0; 0]; % Errore iniziale: 100m radiale, -50m lungo-traccia...

for k = 1:Nsamp-1
    % Calcolo dell'ingresso di controllo u = -K*zeta
    uk = -Kk(:,:,k) * zeta(:,k);
    
    % Evoluzione dello stato (Eq. 11)
    zeta(:,k+1) = Adisc(:,:,k) * zeta(:,k) + Bdisc(:,:,k) * uk;
end

% Plot della posizione 3D
figure;
plot3(zeta(1,:), zeta(2,:), zeta(3,:), 'LineWidth', 2);
grid on; hold on;
plot3(0,0,0, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r'); % Il Target
xlabel('x [m] (Radiale)'); ylabel('y [m] (Lungo-traccia)'); zlabel('z [m]');
title('Traiettoria di Rendezvous monoagente');

