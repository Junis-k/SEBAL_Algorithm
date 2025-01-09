# Modèle SEBAL implémenté dans GRASS GIS
Modèle automatisé de récupération de l'évapotranspiration réelle avec des données climatiques au sol utilisant des images satellites Landsat 8 OLI/TIRS

## Exigences

* GRASS GIS 7.X [https://grasswiki.osgeo.org/wiki/GRASS-Wiki](https://grasswiki.osgeo.org/wiki/GRASS-Wiki)
* Python 3.1.X [https://www.python.org/](https://www.python.org/)
* Images Landsat 8 OLI/TIRS [Earth Explorer](http://earthexplorer.usgs.gov/)
* Données climatiques de la zone d'étude (vitesse du vent, évapotranspiration de référence)

## Documentation

Bastiaanssen, W. G. M., Menenti, M., Feddes, R. A., & Holtslag, A. A. M. (1998). 
A remote sensing surface energy balance algorithm for land (SEBAL). 1. Formulation. 
Journal of Hydrology, 212‑213, 198‑212. <https://doi.org/10.1016/S0022-1694(98)00253-4>

Bastiaanssen, W. G. M., Molden, D. J., & Makin, I. W. (2000). Remote sensing for irrigated agriculture : 
Examples from research and possible applications. 
Agricultural Water Management, 46(2), 137‑155. <https://doi.org/10.1016/S0378-3774(00)00080-9>

Hamimed, A., & Khaldi, A. (2000). Utilisation des données satellitaires TM de Landsat 
pour le suivi de l’état hydrique d’un couvert végétal dans les conditions semi -arides en Algérie. Teledetection, 2(29‑38), 10.
