#!/bin/python3
import math
import matplotlib.pyplot as plt
import numpy as np

MASSE_VOL_AIR = 1.225 # kg/m3
MASSE_VOL_HYDROGENE = 0.082 # kg/m3

COEF_PENETRATION_AEROSTAT = 0.04 #

# Valeurs par défaut pour les paramètres variables
_charge_utile = 1.5 #215000. # kg
_ratio = 3

# Valeurs globales
_surface_penetrante = 0

def demander_utilisateur():
    charge_utile = float(input("Charge utile de votre aérostat(kg): "))
    ratio = float(input("Ratio Longueur/Hauteur: "))


def get_resistance_air(vitesse, masse_vol_air = MASSE_VOL_AIR):
    return 0.5 * masse_vol_air * vitesse**2 * COEF_PENETRATION_AEROSTAT * _surface_penetrante


def get_dimensions_aerostat(volume, ratio):
    # Calcul des dimensions de l'aérostat
    diametre = (12*volume/(_ratio*1.3**2*math.pi)) ** (1./3.)
    longueur = diametre * ratio
    surface_ballon = 2./3.*longueur**(1/2)*2*math.pi*1.3*diametre
    return diametre, longueur, surface_ballon


def main():
    pouvoir_porteur_hydrogene = MASSE_VOL_AIR - MASSE_VOL_HYDROGENE
    #demander_utilisateur()
    volume = _charge_utile/pouvoir_porteur_hydrogene

    # Calcul des dimensions de l'aérostat
    diametre, longueur, surface_ballon = get_dimensions_aerostat(volume, _ratio)

    # Calcul de la force nécessaire pour se mouvoir à 130 km/h
    vitesse = 130 * 1000.0 / 3600.0 # m/s
    global _surface_penetrante
    _surface_penetrante = (diametre/2.)**2*math.pi
    resistance_aerodynamique = get_resistance_air(vitesse)

    print("------  Resultats  -------")
    print("Volume d'Hydrogène nécessaire: "+str(volume)+" m3")
    print("Dimensions de l'aerostat: "+str(round(longueur,2))+" m de longueur sur "+str(round(diametre,2))+" m de hauteur")
    print("Surface de l'aérostat: "+str(round(surface_ballon,2))+" m²")
    print("Resistance aérodynamique: "+str(round(resistance_aerodynamique))+" N")
    print("Resistance aérodynamique pifométrique: "+str(round(resistance_aerodynamique*1.2))+" N")


def get_pression_densite(altitude):
    p0 = 101325 # Pa
    T0 = 288.15 # K
    g = 9.80665 # m/s²
    L = 0.0065 # K/m
    R = 8.31447 # J/(mol*K)
    M = 0.0289654 # kg/mol   masse molaire de l'air sec

    pression_air = p0*(1 - L*altitude/T0)**(g*M/(R*L))
    densite_air = pression_air*M/(R*T0*(1 - L*altitude/T0))
    temperature = T0 - L*altitude
    return pression_air, densite_air, temperature


def courbes_altitude_energie():
    vitesses = [100, 110, 121, 133, 146, 161, 177, 195]

    datas = {} # {100:[[],[]], ...}
    max_energie = 0
    for vitesse in vitesses:
        datas[vitesse] = [[],[]]
        vitesse_m_par_sec = vitesse / 3.6
        for altitude in range(0, 10000, 100):
            pression_air, densite_air, temperature = get_pression_densite(altitude)
            energie = get_resistance_air(vitesse_m_par_sec, densite_air)
            if max_energie < energie:
                max_energie = energie
            datas[vitesse][0].append(altitude)
            datas[vitesse][1].append(energie)

    plt.plot(datas[vitesses[0]][0], datas[vitesses[0]][1], 'b-', label=str(vitesses[0])+" km/h")
    vitesses.remove(vitesses[0])
    for vitesse in vitesses:
        plt.plot(datas[vitesse][0], datas[vitesse][1], 'r-', label=str(vitesse)+" km/h")

    plt.ylabel('Energie (N)')
    plt.xlabel('Altitude (m)')

    plt.xlim(0, 10000)
    plt.xticks(np.arange(0, 10000, step=500))
    graduation_energie = np.arange(0, max_energie, step=max_energie/20)
    for i,val in enumerate(graduation_energie):
        graduation_energie[i] = round(val/100)*100
    plt.yticks(graduation_energie)
    plt.legend()
    plt.grid(color='#AAAAAA', linestyle='-', linewidth=1)
    plt.show()

def courbes_pression_densite():
    courbe_pression_densite = [[],[],[],[]]
    for altitude in range(0, 10000, 100):
        pression_air, densite_air, temperature = get_pression_densite(altitude)
        courbe_pression_densite[0].append(altitude)
        courbe_pression_densite[1].append(pression_air/101325)
        courbe_pression_densite[2].append(densite_air/MASSE_VOL_AIR)
        courbe_pression_densite[3].append(temperature/288.15)

    plt.plot(courbe_pression_densite[1], courbe_pression_densite[0], 'g-', label='pression')
    plt.plot(courbe_pression_densite[2], courbe_pression_densite[0], 'r-', label='densite')
    plt.plot(courbe_pression_densite[3], courbe_pression_densite[0], 'b-', label='temperature')
    plt.legend()
    plt.xlim(.2, 1)
    plt.ylim(0, 10000)
    plt.yticks(np.arange(0, 10000, step=500))
    plt.xticks(np.arange(0.2, 1.0, step=.1))
    plt.grid(color='#AAAAAA', linestyle='-', linewidth=1)
    plt.xlabel('Ratio')
    plt.ylabel('Altitude (m)')
    plt.show()


def get_densite_hy(pression, temperature):
    M = .00201565 # kg/mol
    R = 8.31447 # J/(mol*K)
    densite = pression * M / (R * temperature)
    return densite

def get_portance(altitude):
    pression_air, densite_air, temperature = get_pression_densite(altitude)
    densite_hy = get_densite_hy(pression_air, temperature)
    return densite_air - densite_hy


def courbes_portance():
    courbe_portance = [[],[]]
    for altitude in range(0, 10000, 100):
        courbe_portance[0].append(altitude)
        pression_air, densite_air, temperature = get_pression_densite(altitude)
        densite_hy = get_densite_hy(pression_air, temperature)
        portance = densite_air - densite_hy
        nouveau_vol = MASSE_VOL_HYDROGENE / densite_hy
        courbe_portance[1].append(portance*nouveau_vol)

    plot1 = plt.figure(1)
    plt.plot(courbe_portance[0], courbe_portance[1], 'g-')
    #plt.xlim(.2, 1)
    #plt.ylim(0, 10000)
    plt.xticks(np.arange(0, 10000, step=500))
    #plt.xticks(np.arange(0.2, 1.0, step=.1))
    plt.grid(color='#AAAAAA', linestyle='-', linewidth=1)
    plt.ylabel('Portance (kg/m3)')
    plt.xlabel('Altitude (m)')

    courbe_portance = [[],[]]
    for altitude in range(0, 25000, 100):
        courbe_portance[0].append(altitude)
        pression_air, densite_air, temperature = get_pression_densite(altitude)
        densite_hy = get_densite_hy(pression_air, temperature)
        portance = densite_air - densite_hy
        nouveau_vol = MASSE_VOL_HYDROGENE / densite_hy
        courbe_portance[1].append(nouveau_vol)

    plot1 = plt.figure(2)
    plt.plot(courbe_portance[0], courbe_portance[1], 'g-')
    #plt.xlim(.2, 1)
    plt.xlim(0, 10000)
    plt.xticks(np.arange(0, 10000, step=500))
    #plt.xticks(np.arange(0.2, 1.0, step=.1))
    plt.grid(color='#AAAAAA', linestyle='-', linewidth=1)
    plt.ylabel('Volume (m3)')
    plt.xlabel('Altitude (m)')
    plt.show()


def courbes_ratio_volume_diametre():
    datas = [[],[]]
    for volume in range (1, 100000):
        diametre, longueur, surface_ballon = get_dimensions_aerostat(volume, _ratio)
        datas[0].append(volume)
        datas[1].append(surface_ballon)

    plt.plot(datas[0], datas[1], 'g-')
    #plt.xlim(.2, 1)
    #plt.ylim(0, 10000)
    #plt.xticks(np.arange(0, 10000, step=500))
    #plt.xticks(np.arange(0.2, 1.0, step=.1))
    plt.grid(color='#AAAAAA', linestyle='-', linewidth=1)
    plt.ylabel('Surface (m²)')
    plt.xlabel('Volume (m3)')
    plt.show()



def show_courbes():
    courbes_altitude_energie()
    courbes_pression_densite()
    courbes_portance()
    courbes_ratio_volume_diametre()


main()
#show_courbes()
