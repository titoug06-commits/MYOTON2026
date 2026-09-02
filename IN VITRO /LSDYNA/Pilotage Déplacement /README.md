        Simulation LSDYNA de l'impact MyotonPro sur silicone multi-couche 
        PILOTAGE EN DÉPLACEMENT
  
1. Objectif de la simulation:
Modélisation numérique du comportement mécanique d'un empilement de caoutchouc de silicone, 
imitant la superposition de la peau, le gras, le fascia et le muscle, soumis à 
l'enfoncement (déplacement imposé) d'une sonde rigide.

 
2. Géométrie et Matériaux

Le modèle intègre un corps rigide et 4 couches hyperélastiques/viscoélastiques 
(loi d'Ogden et loi de Prony 'MAT_OGDEN_RUBBER_TITLE (MAT_077_0)') :

  Sonde (PID 1) : Corps rigide ('MAT_RIGID').
  Peau (PID 2 & 21) : Dragon Skin 30 (1.56 mm).
  Gras (PID 3 & 31) : Ecoflex 00-10 (3.36 mm).
  Fascia (PID 4 & 41) : Dragon Skin 30 (2.70 mm).
  Muscle (PID 5 & 51) : Ecoflex 00-50 (26.03 mm).

Notes : Les paramètres matériaux sont ceux obtenus via les fits expérimentaux 
(cf. dossier Fits_Ogden_Prony, tableau dans le mémoire). En LS-DYNA, la normalisation 
Prony utilise VFLAG=1, où g_∞ = 1 - Σg_i.

 
3. Conditions aux limites et Chargements

  Mouvement imposé : Déplacement dynamique de la sonde selon l'axe Z piloté par la 
  courbe 'LCID 3' ('BOUNDARY_PRESCRIBED_MOTION_RIGID_ID'), obtenue à partir des 
  données de stéréocorrélation (DIC) lors d'une impulsion isolée mesurée sur 
  l'empilement à gras faible (3,36 mm). Durée d'impulsion ~7 ms, suivie d'une 
  relaxation sur ~200 ms (total 210 ms).
  
  Blocages (SPC) : Le plateau de base est totalement bloqué en translation 
  (dX, dY, dZ) via 'BOUNDARY_SPC_SET_ID'.
  
  Gravité : Champ de force appliqué suivant l'axe Z ('LOAD_BODY_Z').
  
  Frontières non-réfléchissantes : Absorption des ondes de bord ('BOUNDARY_NON_REFLECTING').

 
4. Contacts et Paramètres de calcul

  Contacts : Interactions de type surface-sur-surface automatiques 
  ('CONTACT_AUTOMATIC_SURFACE_TO_SURFACE') avec option 'SOFT=2' segment-to-segment.
  
  Temps de simulation : Durée totale fixée à t = 210.0 ms ('CONTROL_TERMINATION').
  
  Sorties : Génération des fichiers binaires de résultats (D3PLOT toutes les 2.0 ms) 
  et des historiques globaux (GLSTAT, RCFORC, NODEFOR, fréquence 0,1 ms).
