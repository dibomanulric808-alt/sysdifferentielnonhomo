# Systèmes Différentiels Non Homogènes (R Implementation)

Ce dépôt contient l'implémentation en **R** des exemples théoriques et des méthodes présentés dans l'exposé de mathématiques sur les **Systèmes Différentiels Non Homogènes** ($X'(t) = A(t)X(t) + B(t)$).

## Forme générale
$$\mathbf{X}'(t) = \mathbf{A}(t)\mathbf{X}(t) + \mathbf{B}(t)$$

- **Solution générale** : $X(t) = X_h(t) + X_p(t)$
- **Variation des constantes** : $X(t) = e^{At}C + e^{At} \int e^{-As}B(s) ds$

## Contenu du projet
- `R/main_systems_ode.R` : Script R principal avec résolutions analytiques et numériques.
- `output/` : Emplacement pour les graphiques générés.

## Prérequis
Installez les packages R suivants avant d'exécuter le script :
```R
install.packages(c("deSolve", "expm", "ggplot2", "gridExtra"))
```

## Structure du dépôt GitHub
```text
ode-non-homogeneous-systems/
├── README.md
├── LICENSE
├── R/
│   └── main_systems_ode.R
└── output/
```
