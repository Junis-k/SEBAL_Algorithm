# -*- coding: utf-8

import numpy
import math 
import os
import grass.script as grass
from os import system
from grass.script import core as g
  
system('clear')

MTLfile = [i for i in os.listdir('.') if i.endswith('MTL.txt')]

files = [i for i in os.listdir('.') if i.endswith(('.TIF','.tif'))]

u_2m = float(input('Placez la valeur de la vitesse du vent pour la hauteur de 2 m mesurée dans la station météorologique - (m/s) : '))

EToi = float(input('Placez la valeur instantanée de l évapotranspiration de référence (EToi) de la station météorologique, pour le moment du passage du satellite - (mm) : '))

ETo = float(input('Placez la valeur quotidienne de l évapotranspiration de référence (ETo) de la station météorologique - (mm) :'))

g.parse_command('g.region',flags='p',rast='MDT_Sebal@PERMANENT',quiet=True)

runCC = g.parse_command('g.list',type='raster', pattern='CC_432')
runRLo = g.parse_command('g.list',type='raster', pattern='RLo')

if runCC == {}:
        print ('Importation d images Landsat 8, soyez patient...')
        for i in range(len(files)):
                g.parse_command('r.in.gdal', input=files[i], output=os.path.splitext(files[i])[0], overwrite=True)
        print ('Terminé!')
 
        print ('Calcule la réflectance et la température du sommet de l atmosphère pour Landsat 8, soyez patient...')
        g.parse_command('i.landsat.toar', input=files[0].split('_B')[0] + '_B', output= 'LS8_corre', metfile=MTLfile, sensor='oli8',overwrite=True)
        print ('Terminé!')     

        g.parse_command('g.remove', type='raster', pattern=files[0].split('_B')[0] + '*', flags = 'f')

        print ('Composition colorée R=B4 G=B3 B=B2 Landsat 8, soyez patient...')
        grass.run_command('i.colors.enhance', red='LS8_corre4',green='LS8_corre3',blue='LS8_corre2',quiet=True)
        grass.run_command('r.composite',red='LS8_corre4', green='LS8_corre3', blue='LS8_corre2', output='CC_432',quiet=True,overwrite=True)
        print ('Terminé!')
if runRLo == {}:
        print ('Calcul NDVI...')
        grass.mapcalc('NDVI=($LS8_corre5-$LS8_corre4)/($LS8_corre5+$LS8_corre4)',
                LS8_corre5='LS8_corre5',
                LS8_corre4='LS8_corre4',
                overwrite='true',
                quiet='true')
        print ('Terminé!')
        
        print ('Calcul de SAVI avec L égal à 0,5...')
        Lsavi = 0.5 #float(input('Entre com o valor de L: '))
        grass.mapcalc('SAVI= (($LS8_corre5-$LS8_corre4)/($LS8_corre5+$LS8_corre4+$Lsavi))*(1+$Lsavi)',
                LS8_corre4='LS8_corre4',
                LS8_corre5='LS8_corre5',
                Lsavi=Lsavi,
                overwrite='true',
                quiet='true')
        print ('Terminé!')
        
        print ('Calcul de l indice de surface foliaire (LAI)...')
        grass.mapcalc('LAI=if(SAVI < 0.1, 0.00001,(if(0.1 < SAVI && SAVI < 0.687,-log((0.69-SAVI)/0.59)/0.91,if(SAVI > 0.687,6,0))))',
                overwrite='true',
                quiet='true')
        print ('Terminé!')
        
        print ('Calcul de l émissivité de surface à bande étroite (eNBf)...')
        grass.mapcalc('eNB=if(LAI<3 && NDVI>0.,0.97+0.0033*LAI,if(LAI>=3 && NDVI>0.,0.98))',
                overwrite='true',
                quiet='true')
        grass.mapcalc('eNBf=if(eNB == 0, 0.99, eNB)',
                overwrite='true',
                quiet='true')
        print ('Terminé!')
        
        print ('Calcul de l émissivité de surface à large bande (e0f)...')
        grass.mapcalc('e0=if(LAI<3 && NDVI>0.,0.95+0.01*LAI,if(LAI>=3 && NDVI>0.,0.98))',
                overwrite='true',
                quiet='true') 
        grass.mapcalc('e0f=if(e0 == 0, 0.985, e0)',
                overwrite='true',
                quiet='true')
        print ('Terminé!')
        
        print ('Calcul la temperature de surface (Ts) - K...')
        grass.mapcalc('Ts = $LS8_corre10/(1+((10.8*$LS8_corre10)/14380)*log($eNBf))',
              LS8_corre10='LS8_corre10', 
              eNBf='eNBf',               
              overwrite='true',
              quiet='true')
        print ('Terminé!')
        
        Ts_median=g.parse_command('r.univar', flags='ge', map='Ts', quiet = True)['median']        
        print ('Temperature moyenne de surface:', Ts_median,'K')
        
        for line in open(MTLfile[0]):
                if 'EARTH_SUN_DISTANCE' in line:
                        d=float(line.split('=')[-1])
                elif'RADIANCE_MAXIMUM_BAND_1'in line:
                        RADIANCE_MAXIMUM_BAND_1=float(line.split('=')[-1])
                elif 'RADIANCE_MAXIMUM_BAND_2' in line:
                        RADIANCE_MAXIMUM_BAND_2=float(line.split('=')[-1])
                elif 'RADIANCE_MAXIMUM_BAND_3' in line:
                        RADIANCE_MAXIMUM_BAND_3=float(line.split('=')[-1])
                elif 'RADIANCE_MAXIMUM_BAND_4' in line:
                        RADIANCE_MAXIMUM_BAND_4=float(line.split('=')[-1])
                elif 'RADIANCE_MAXIMUM_BAND_5' in line:
                        RADIANCE_MAXIMUM_BAND_5=float(line.split('=')[-1])
                elif 'RADIANCE_MAXIMUM_BAND_6' in line:
                        RADIANCE_MAXIMUM_BAND_6=float(line.split('=')[-1])
                elif 'RADIANCE_MAXIMUM_BAND_7' in line:
                        RADIANCE_MAXIMUM_BAND_7=float(line.split('=')[-1])
                elif 'REFLECTANCE_MAXIMUM_BAND_1' in line:
                        REFLECTANCE_MAXIMUM_BAND_1=float(line.split('=')[-1])
                elif 'REFLECTANCE_MAXIMUM_BAND_2' in line:
                        REFLECTANCE_MAXIMUM_BAND_2=float(line.split('=')[-1])
                elif 'REFLECTANCE_MAXIMUM_BAND_3' in line:
                        REFLECTANCE_MAXIMUM_BAND_3=float(line.split('=')[-1])
                elif 'REFLECTANCE_MAXIMUM_BAND_4' in line:
                        REFLECTANCE_MAXIMUM_BAND_4=float(line.split('=')[-1])
                elif 'REFLECTANCE_MAXIMUM_BAND_5' in line:
                        REFLECTANCE_MAXIMUM_BAND_5=float(line.split('=')[-1])
                elif 'REFLECTANCE_MAXIMUM_BAND_6' in line:
                        REFLECTANCE_MAXIMUM_BAND_6=float(line.split('=')[-1])
                elif 'REFLECTANCE_MAXIMUM_BAND_7' in line:
                        REFLECTANCE_MAXIMUM_BAND_7=float(line.split('=')[-1])
                elif 'SUN_ELEVATION' in line:
                        SUN_ELEVATION=float(line.split('=')[-1])

        ESUN_B1=(math.pi*d*d)*(RADIANCE_MAXIMUM_BAND_1/REFLECTANCE_MAXIMUM_BAND_1) 
        ESUN_B2=(math.pi*d*d)*(RADIANCE_MAXIMUM_BAND_2/REFLECTANCE_MAXIMUM_BAND_2) 
        ESUN_B3=(math.pi*d*d)*(RADIANCE_MAXIMUM_BAND_3/REFLECTANCE_MAXIMUM_BAND_3) 
        ESUN_B4=(math.pi*d*d)*(RADIANCE_MAXIMUM_BAND_4/REFLECTANCE_MAXIMUM_BAND_4) 
        ESUN_B5=(math.pi*d*d)*(RADIANCE_MAXIMUM_BAND_5/REFLECTANCE_MAXIMUM_BAND_5) 
        ESUN_B6=(math.pi*d*d)*(RADIANCE_MAXIMUM_BAND_6/REFLECTANCE_MAXIMUM_BAND_6) 
        ESUN_B7=(math.pi*d*d)*(RADIANCE_MAXIMUM_BAND_7/REFLECTANCE_MAXIMUM_BAND_7) 
        ESUN=[ESUN_B1,ESUN_B2,ESUN_B3,ESUN_B4,ESUN_B5,ESUN_B6,ESUN_B7]
        print (ESUN)
      
        W = []
        for i in range(len(ESUN)):
                W += [ESUN[i]/sum(ESUN)]
        W1=W[0]        
        W2=W[1]
        W3=W[2]
        W4=W[3]
        W5=W[4]
        W6=W[5]
        W7=W[6]
        print (W)
        
        print ('Calcul de l albédo au sommet de l atmosphère (aTOA)...')
        grass.mapcalc('aTOA=$LS8_corre1*$W1+$LS8_corre2*$W2+$LS8_corre3*$W3+$LS8_corre4*$W4+$LS8_corre5*$W5+$LS8_corre6*$W6+$LS8_corre7*$W7',
                      LS8_corre1='LS8_corre1',W1=W1,
                      LS8_corre2='LS8_corre2',W2=W2, 
                      LS8_corre3='LS8_corre3',W3=W3, 
                      LS8_corre4='LS8_corre4',W4=W4, 
                      LS8_corre5='LS8_corre5',W5=W5, 
                      LS8_corre6='LS8_corre6',W6=W6, 
                      LS8_corre7='LS8_corre7',W7=W7,         
                      overwrite='true',
                      quiet='true')
        print ('Terminé!')              
        print ('Calcul de la transmissivité de l air en ondes courtes (Tsw)...')
        grass.mapcalc('Tsw=0.75+0.00002*$MDT_Sebal',
                      MDT_Sebal='MDT_Sebal',
                      overwrite='true',
                      quiet='true')
        print ('Terminé!')
                   
        print ('Calcul l albédo de surface (aS)...')
        grass.mapcalc('aS=($aTOA-0.03)/($Tsw^2)',
                      aTOA='aTOA',
                      Tsw='Tsw',
                      overwrite='true',
                      quiet='true')
        print ('Terminé!')
        
        print ('Calcul du rayonnement de courte longueur d onde entrant (Rsi) - W/m2...')
        grass.mapcalc('Rsi=1367*cos(90-$SUN_ELEVATION)*(1/($d^2))*$Tsw',
                      SUN_ELEVATION=SUN_ELEVATION,
                      Tsw='Tsw',
                      d=d,
                      overwrite='true',
                      quiet='true')
        print ('Terminé!')
        
        print ('Calcul du rayonnement sortant à grande longueur d onde (RLo) - W/m2...')
        grass.mapcalc('RLo=$e0f*5.67e-8*$Ts^4',
                      Ts='Ts',
                      e0f='e0f',
                      overwrite='true',
                      quiet='true')
        print ('Terminé!')

print ('Fabrication du masque de pixels froids...')
Ts_median=g.parse_command('r.univar', flags='ge', map='Ts', quiet = True)['median']        
grass.mapcalc('Pcold=if($NDVI>0.4 && $Ts<$Ts_median,$Ts,null())',
              NDVI='NDVI',
              aS='aS',
              Ts='Ts',
              Ts_median=Ts_median,
              overwrite='true',
              quiet='true')
print ('Choisissez les coordonnées des pixels froids dans les zones d irrigation. Utilisez les matrices Pcold et CC_432 pour vous aider...')
xy_Pcold = str(input('Placez les coordonnées (est,nord): ')).strip('()')
print (xy_Pcold)
TsPcold_z=g.parse_command('r.what', map='Ts', coordinates=xy_Pcold)

print("####################################" , TsPcold_z, list(TsPcold_z.keys()), "*********************")

z_TsPcold = float(list(TsPcold_z.keys())[0].split('|')[3])
##### z_TsPcold=float(dict.keys(TsPcold_z)[0].split('|')[3])

print ('Température du pixel froid:',z_TsPcold,'K')
grass.write_command('v.in.ascii', input='-', output='Pcold', sep=',', stdin=xy_Pcold, overwrite='true', quiet='true')

print ('Calcul du rayonnement de grande longueur d onde entrant (RLi) - W/m2...')
grass.mapcalc('RLi=0.85*((-log($Tsw))^0.09)*5.67e-8*$z_TsPcold^4',
              z_TsPcold=z_TsPcold,
              Tsw='Tsw',
              overwrite='true',
              quiet='true')
print ('Terminé!')
print ('Calcul le rayonnement net (Rn) - W/m2...')
grass.mapcalc('Rn=(1-$aS)*$Rsi+$RLi-$RLo-(1-$e0f)*$RLi',
              aS='aS',
              RLi='RLi',
              RLo='RLo',
              Rsi='Rsi',
              e0f='e0f',
              overwrite='true',
              quiet='true')
print ('Terminé!')
print ('Calcul le flux de chaleur du sol (G) - W/m2...')
grass.mapcalc('G_Rn=if(NDVI<0,0.5,(($Ts-273.15)/$aS)*(0.0038*$aS+0.0074*$aS^2)*(1-0.98*$NDVI^4))',
              aS='aS',
              Ts='Ts',
              #Rn='Rn',
              NDVI='NDVI',
              overwrite='true',
              quiet='true')
grass.mapcalc('G=$G_Rn*$Rn',
              G_Rn='G_Rn',
              Rn='Rn',
              overwrite='true',
              quiet='true')
print ('Terminé!')
print ('Fabrication du masque de pixels chauds...')
grass.mapcalc('Phot=if($SAVI>0.18 && $SAVI<0.3, $Ts,null())',
              SAVI='SAVI',
              aS='aS',
              Ts='Ts',
              Ts_median=Ts_median,
              #Ts_max=Ts_max,
              overwrite='true',
              quiet='true')
print ('Choisissez les coordonnées du pixel chaud dans les zones de sol nu. Utilisez le raster Phot et CC_432 pour vous aider...')
xy_Phot = str(input('Placez les coordonnées (est,nord): ')).strip('()')
print (xy_Phot)
print ('Terminé!')
TsPhot_z=g.parse_command('r.what', map='Ts', coordinates=xy_Phot)

z_TsPhot=float(list(TsPhot_z.keys())[0].split('|')[3])
###### z_TsPhot=float(dict.keys(TsPhot_z)[0].split('|')[3])

print ('Temperature du pixel chaud:',z_TsPhot,'K')
grass.write_command('v.in.ascii', input='-', output='Phot', sep=',', stdin=xy_Phot, overwrite='true', quiet='true')
print ('Calcul de la vitesse de frottement (u*) pour la station météorologique - m/s...')
h=0.15
Zom=0.123*h
u_ast=0.41*u_2m/(math.log(2/Zom))
u_200m=u_ast*(math.log(200/Zom))/0.41
print ('Terminé!')
print ('Vitesse de frottement', u_ast, 'm/s')
print ('Calcul de la carte de longueur de rugosité de l élan (Z0map) - m...')
grass.mapcalc('Z0map=exp(-5.809+5.62*$SAVI)',
              SAVI='SAVI',
              overwrite='true',
              quiet='true')
print ('Terminé!')
print ('Calcul de la carte des vitesses de frottement (u*map) - m/s...')
grass.mapcalc('u_astmap=0.41*$u_200m/log(200/$Z0map)',
              Z0map='Z0map',
              u_200m=u_200m,
              overwrite='true',
              quiet='true')
print ('Terminé!')
print ('Calcul de la carte de résistance aérodynamique au transport de chaleur en termes de stabilité neutre (rah) - s/m...')
grass.mapcalc('rah=log(2/0.1)/($u_astmap*0.41)',
              u_astmap='u_astmap',
              overwrite='true',
              quiet='true')
print ('Terminé!')

GPhot_z=g.parse_command('r.what', map='G', coordinates=xy_Phot)

z_GPhot=float(list(GPhot_z.keys())[0].split('|')[3])
##### z_GPhot=float(dict.keys(GPhot_z)[0].split('|')[3])

rahPhot_z=g.parse_command('r.what', map='rah', coordinates=xy_Phot)

z_rahPhot=float(list(rahPhot_z.keys())[0].split('|')[3])
##### z_rahPhot=float(dict.keys(rahPhot_z)[0].split('|')[3])

RnPhot_z=g.parse_command('r.what', map='Rn', coordinates=xy_Phot)

z_RnPhot=float(list(RnPhot_z.keys())[0].split('|')[3])
##### z_RnPhot=float(dict.keys(RnPhot_z)[0].split('|')[3])

print ('Correction de la résistance aérodynamique, soyez patient...')

z_rahPhot_i=0
i=1
while (abs(z_rahPhot_i - z_rahPhot) > 0.00001):
        print ('Nombre d itérations:',i)
        i = i+1
        z_rahPhot_i=z_rahPhot
        print ('Estimation des coefficients d équation >> dT = a.Ts + b')
        a=(z_RnPhot-z_GPhot)* z_rahPhot_i/((z_TsPhot-z_TsPcold)*1.25*1004)
        b=-a*z_TsPcold
        print ('a:',a, 'b:', b )
        print (' - K...')
        grass.mapcalc('dT=$a*$Ts+$b',
                      Ts='Ts',
                      a=a,
                      b=b,
                      overwrite='true',
                      quiet='true')
        print ('Terminé!')
        print ('Calcul le flux de chaleur sensible (H) - W/m2...')
        grass.mapcalc('H=($dT/$rah)*1.25*1004',
                      dT='dT',
                      rah='rah',
                      overwrite='true',
                      quiet='true')
        print ('Terminé!')
        print ('Calcul de la carte de longueur Monin-Obukhov (L) - m...')
        grass.mapcalc('L=-(1.25*1004*$Ts*$u_astmap^3)/(0.41*9.81*$H)',
                      Ts='Ts',
                      u_astmap='u_astmap',
                      H='H',
                      overwrite='true',
                      quiet='true')
        print ('Terminé!')
        print ('Calcul de la correction de la stabilité atmosphérique (L200m, L2m, L01m)...')
        grass.mapcalc('L200m=if($L<0,2*log((1+(1-16*(200/$L))^0.25)/2)+log((1+(1-16*(200/$L))^0.5)/2)-2*atan((1-16*(200/$L))^0.25)+0.5*3.14159265,if($L>0,-5*(2/$L),0))',
                      L='L',
                      overwrite='true',
                      quiet='true')
        grass.mapcalc('L2m=if($L<0,2*log((1+(1-16*(2/$L))^0.5)/2),if($L>0,-5*(2/$L),0))',
                      L='L',
                      overwrite='true',
                      quiet='true')
        grass.mapcalc('L01m=if($L<0,2*log((1+(1-16*(0.1/$L))^0.5)/2),if($L>0,-5*(0.1/$L),0))',
                      L='L',
                      overwrite='true',
                      quiet='true')
        print ('Terminé!')
        print ('Calcul de la carte de vitesse de friction corrigée (u*map) - m/s...')
        grass.mapcalc('u_astmap=0.41*$u_200m/(log(200/$Z0map)-$L200m)',
                      Z0map='Z0map',
                      u_200m=u_200m,
                      L200m='L200m',
                      overwrite='true',
                      quiet='true')
        print ('Terminé!')
        print ('Calcul de la résistance aérodynamique corrigée au transport de chaleur (rah) - s/m...')
        grass.mapcalc('rah=(log(2/0.1)-$L2m+$L01m)/($u_astmap*0.41)',
                      u_astmap='u_astmap',
                      L2m='L2m',
                      L01m='L01m',
                      overwrite='true',
                      quiet='true')
        print ('Terminé!')
        rahPhot_z=g.parse_command('r.what', map='rah', coordinates=xy_Phot)

        z_rahPhot=float(list(rahPhot_z.keys())[0].split('|')[3])
        ##### z_rahPhot=float(dict.keys(rahPhot_z)[0].split('|')[3])

HPhot_z=g.parse_command('r.what', map='H', coordinates=xy_Phot)

z_HPhot=float(list(HPhot_z.keys())[0].split('|')[3])
##### z_HPhot=float(dict.keys(HPhot_z)[0].split('|')[3])

dTPhot_z=g.parse_command('r.what', map='dT', coordinates=xy_Phot)

z_dTPhot=float(list(dTPhot_z.keys())[0].split('|')[3])
##### z_dTPhot=float(dict.keys(dTPhot_z)[0].split('|')[3])

print ('Résultats du pixel chaud, vérifiez si Rhhot - Ghot = Hhot')
print ('Hhot:', z_HPhot,)
print ('rahhot:', z_rahPhot,)
print ('Ghot:', z_GPhot,)
print ('Rnhot:', z_RnPhot, 'dT:',z_dTPhot)
print ('Calcul du flux de chaleur latente (LET) - W/m2...')
grass.mapcalc('LET=Rn-G-H',
              H='H',
              Rn='Rn',
              G='G',
              overwrite='true',
              quiet='true')
print ('Terminé!')
print ('Calcul de l évapotranspiration instantanée (ETi) - mm/h...')
grass.mapcalc('ETi=if(3600*(LET/(2.45*10^6))<0,0,3600*(LET/(2.45*10^6)))',
              LET='LET',
              overwrite='true',
              quiet='true')
print ('Terminé!')
print ('Calcul de la fraction d évapotranspiration de référence (ETof)...')
grass.mapcalc('ETof=$ETi/$EToi',
              ETi='ETi',
              EToi=EToi,
              overwrite='true',
              quiet='true')
print ('Terminé!')
print ('Calcul de l évapotranspiration quotidienne (ETday) - mm/day...')
grass.mapcalc('ETday=$ETof*$ETo',
              ETof='ETof',
              ETo=ETo,
              overwrite='true',
              quiet='true')
print ('Terminé!')
