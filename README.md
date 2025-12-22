## Calcul erreur Runge-Kutta
**RK4 error (dt = 0.005):**
0.04839186301086329 %

**RK4 error (dt = 0.01):**
1.1230357354165932 %

**RK4 error (dt = 0.001):**
1.848640525281553e-5 %

## Idée
- Utiliser Python pour avoir la positions des pendules.
  - OpenCV
  - masquage des couleurs avec couleur rouge pendule comme base
- Rapprochement avec Gradient descent (descente de gradient) -> travail avec paramètres

## Changé
- Tracker App pour avoir les 70 premières frames.
- Trouver par rapport à ces 70 frames d'abord.
  - Check RMSE (avec premières valeurs mises au bol (0.022 et 0.0018) j'ai trouvé RMSE total de ~0.05)
  un peu de changement à la main et arrivé sur RMSE de 0.015479 (0.020283 et 0.002083)
- un peu de recherches et trouvé manière d'optimiser avec packages "LsqFit" ou "Optim".
  - utilisé Optim car bonne doc https://julianlsolvers.github.io/Optim.jl/stable/user/minimization/ avec descente de gradient (indice du prof)
  - nouveau fichier -> optimisation_pendule.jl