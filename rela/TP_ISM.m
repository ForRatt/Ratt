

% TP ISM

% Filtrage et calcul du Seeing.
% ----------------------------------

% 2) Initilisation et commentaires.
clear all              
close all              
Nx = 640;               % Largeur des images en pixels.
Ny = 480;               % Longeur des images en pixels.
Nb = 128;
Ng = 2;                % Nom d'images � traiter.
bp = 25;
X = zeros(Nx,Ny,3);     % Initialisation de la matrice de l'image.
                        % La matrice X est tridimensionnel.
A = zeros(1,Ng); B = A; Dg_maxI = A;                       

% 3) Boucle qui lit les images l'une apr�s l'autre.
for i = 1 : Ng          % It�ration sur chaque image (de 1 � 50).
    i = num2str(i);     % Convertir la valeur enti�re de "i" en caract�re.
    nomFichier = strcat(i,'.jpeg'); % Nom du fichier image dans le dossier.
    X = imread(nomFichier); % Mettre les donn�es de l'image dans une 
                            % matrice de taille Nx*Ny*3
% 4) Visualisation des images successivement.
    figure(1);          % Visualisation des images sur la m�me figure.
    image(X);           % Afficher l'image dont les donn�es sont stock�es
                        % dans la matrice X.

% 5) Extraire l'image sous forme 'double'.    
    X0(:,:) = double(X(:,:,1)); % La matrice X0 est bidimensionnel.
                        % La taille de la matrice X0 est % La taille de la
                        % matrice C est [480,640].
    C = max(X,[],2);    % Chercher le maximum dans la matrice X.
                        % La taille de la matrice C est [480,1,3].
                     
% 6) Recherche de la position du maximum de l'image, ensuite la centrer.
% Ici on cherche le centre de la figure de diffraction en utilisant
% la commande Matlab 'max'.
    [max0, ind0] = max(C);
% Renvoie les indices des valeurs maximales dans le vecteur 'ind0'.
% Si les valeurs le long de la premi�re dimension non-singleton contiennent
% plus qu'un �l�ment maximal, l'indice du premier est renvoy�.
    I = X0(ind0(:,:,2),:);
    
% 7) Visualiser la coupe de l'intensit� 'I(x)' et enregistrer un exemple
% sur une figure.
    figure(2);          % Visualiser la courbe de l'intensit� 'I(x)'
                        % de chaque image succissivement sur la figure(2)
    plot((-Nx/2 : Nx/2 - 1), I, 'k--');
    
% 8) Recentrer encore une autre fois.
    [max0, ind0] = max(I);
    I2 = zeros(1,Nx);% centrer l'image 
    i0 = ind0 - Nx/2 - 1;
    for kk = 2 + abs(i0) : Nx - abs(i0) - 1
        I2(kk) = I(kk + i0);
    end
    
% 9) Dans une zone ou il n'y a pas de signal, calculer la valeur moyenne
% 'mean' du bruit dans l'image.
    NoiseI = mean(I2(1 : Nx/2 - Nb));
    
% 10) Retrancher cette valeur de l'image 'I(x)'.
    I = I2 - NoiseI;
    
% 11) Cr�er un tableau (Ib) de taille [128,1] qui contient la coupe de
% l'intensit� I(x).
    Ib = I(Nx/2 + 1 - Nb/2 : Nx/2 + 1 + Nb/2 - 1);
    
% 12) Calculer la transform�e de fourier inverse de l'intensit� (Ib) qui 
% prendre les deux parties : r�elle 'IIr' et imaginaire 'IIi';
    IIi = imag((fftshift(ifft(fftshift(Ib))))).*Nb;
    IIr = real((fftshift(ifft(fftshift(Ib))))).*Nb;

% 13) Visualiser les courbes 'IIi' et 'IIr'.
    figure(4);
    M = max(max(IIr),max(IIi));
    plot((-Nb/2 : Nb/2 - 1), IIr./M, 'r', (-Nb/2 : Nb/2 - 1), IIi./M, 'b');
    grid on;
    
% 14) Calculer des quantit�s A et B.
    x1 = 25;            % On prend la position x1 = 25.
    A(Ng) = IIi(Nb/2 + 1 - x1) - IIi(Nb/2 + 1 + x1);
    B(Ng) = IIr(Nb/2 + 1 - x1) + IIr(Nb/2 + 1 + x1);

% 15) Calculer la valeur moyenne (sur les 50 images) du maximum de
% l'intensit� Imax.
    Dg_maxI(Ng) = max(I); 
end;

maxI0 = mean(Dg_maxI);

% 16) Normalisation des deux quantit�s A et B.
A = A./maxI0;
B = B./maxI0;

% 17) Calculer les variances Va et Vb des deux quantit�s A et B.
Va = mean(var(A));
Vb = mean(var(B));

% 18) Calcul de la fonction de structure D de la phase.
D = Va + 0.8 * Vb;

%19) En d�duire la valeur du sseing.

seeing = 0.6 * D^(3/5);
%disp(['La valeur du Seeing est : ', num2str(seeing)]);

    

