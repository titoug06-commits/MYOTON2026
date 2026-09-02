        Simulation LSDYNA de l'impact MyotonPRO sur silicone multi-couche 
        PILOTAGE EN FORCE (Résultats préliminaires, instabilités observées)

Le fichier .zip contient le modèle lsdyna (fichier .k une fois décompréssé).

1. Objectif de la simulation:
Modélisation numérique du comportement mécanique d'un empilement de
silicone, imitant la superposition peau/gras/fascia/muscle, soumis à une force 
imposée d'une sonde rigide représentant le MyotonPRO. Cette approche vise à tester 
une stratégie alternative au pilotage en déplacement, plus proche des conditions 
réelles où la sonde a une inertie propre.

 
2. Géométrie et Matériaux
Identique au modèle pilotage en déplacement. Quatre couches 
déformables (peau Dragon Skin 30, gras Ecoflex 00-10 3,36 mm, fascia Dragon Skin 30, 
muscle Ecoflex 00-50 26,03 mm) + sonde rigide.

Les paramètres matériaux (Ogden + Prony) sont ceux du tableau du mémoire.

 
3. Conditions aux limites et Chargements

  Mouvement imposé : Force dynamique appliquée à la sonde rigide selon l'axe Z 
  via 'LOAD_RIGID_BODY', piloté par la courbe 'LCID 3'. Cette courbe provient de 
  mesures expérimentales du capteur de force placé directement sous la sonde MyotonPRO 
  (fichier : capteur_force_experimentale.csv, unités : N, temps en ms).
  
  La courbe comporte deux phases :
    - Montée progressive (précharge + impulsion) : 0 → ~0,81 N sur ~20 ms
    - Relaxation et retour : ~0,81 N → 0 N sur ~200 ms
  
  Blocages (SPC) : Plateau de base encastré en translation complète.
  Gravité : Appliquée selon Z ('LOAD_BODY_Z').
  Frontières : Non-réfléchissantes ('BOUNDARY_NON_REFLECTING').

 
4. Contacts et Paramètres de calcul

  Contacts : Surface-surface automatique ('CONTACT_AUTOMATIC_SURFACE_TO_SURFACE') 
  avec option 'SOFT=2'.
  
  Temps : 210 ms total.
  Sorties : D3PLOT (dt=2 ms), GLSTAT, RCFORC, NODEFOR, et ELOUT (dt=0,1 ms).

 
5. Différences par rapport au pilotage en déplacement

  Le pilotage en force impose une **force mesurée expérimentalement** plutôt qu'un 
  déplacement prescrit. Le déplacement de la sonde émerge de la réponse mécanique 
  dynamique du système (inertie + raideur + amortissement du matériau).
  
  Avantage théorique : Plus proche de la physique réelle.
  
  Inconvénient observé : Instabilités numériques significatives et résultats loins de la réalité(voir section 6).

 
6. Résultats et limitations observées

  Les courbes de déformation principale maximale (E1) en fonction du temps, issues 
  du fichier ELOUT pour plusieurs points de profondeur croissante, révèlent une 
  **forte instabilité numérique après le pic d'impulsion** (~t > 20 ms).
  
  Observation détaillée :
    - Pic principal : ~0,30–0,31 (amplitude beaucoup plus élevée que le pilotage 
      en déplacement, x10)
    - Phase de relaxation : oscillations sur-amorties et non-physiques, amplitude 
      restant autour de 0,10–0,12 indéfiniment
    - Pas de convergence vers zéro malgré l'absence de nouvelle excitation
  
  Ces oscillations suggèrent un problème de **convergence du contact dynamique** 
  ou une **sensibilité de l'algorithme aux gradients abrupts** de la courbe de force 
  expérimentale brute.
  
  **Conclusion** : Le pilotage en force, bien que numériquement stable en termes de 
  crash (pas d'arrêt prématuré de la simulation), produit des résultats non exploitables 
  à ce stade. Le pilotage en déplacement (voir dossier voisin) est privilégié pour 
  l'analyse.

 
7. Pistes de correction futures

  - Lissage / filtrage passe-bas de la courbe de force expérimentale
  - Ajustement du paramètre SOFT ou passage à un contact de type FORM=1 (penalty)
  - Augmentation de la finesse du maillage au contact sonde-gras
  - Réduction du timestep initial (DTINIT)
  - Modification des conditions limites et contact (friction, GLUE contact)
