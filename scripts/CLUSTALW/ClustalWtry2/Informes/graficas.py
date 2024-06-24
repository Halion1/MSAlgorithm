import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import json
import os


carpeta_alineamientos = os.path.join(os.getcwd(), "Informes/Alineamientos")

D = []

for archivo in os.listdir(carpeta_alineamientos):
    if archivo.endswith(".json"):
        # Construir la ruta completa al archivo
        ruta_completa = os.path.join(carpeta_alineamientos, archivo)

    datos_json = {}
    with open(ruta_completa, "r") as f:
            datos_json = json.load(f)
            datos_json["Alineamiento"]["Grupo"] = len(datos_json["Cadenas"])
            datos_json["Alineamiento"]["Cadenas"] = datos_json["Cadenas"]
    D.append( datos_json["Alineamiento"] )
D = pd.DataFrame(D)

## Grafica de Longitud segun Grupo

df = D.drop_duplicates(subset=['Grupo'])[['Grupo', 'Cadenas']]

def grupbarplot(df,set_ylabel,set_title,width = 0.3,distancia = 0.1):
    fig, ax = plt.subplots(layout='constrained')

    for columna in df.values:
        
        long = len(columna[1])
        
        if (long %2 !=0):
            long -=1 
        long= long//2
        x  = [a for a in range(-long,long+1,1)]
        if (long %2 !=0):
            x.remove(0)
        x = np.array(x)*(distancia+width)+columna[0]
        rects = ax.bar(x, columna[1], width)
        ax.bar_label(rects, padding=3)

    x = df['Grupo']
    ax.set_xticks(x, df['Grupo'])

    ax.set_ylabel(set_ylabel)
    ax.set_title(set_title)
    plt.show()

grupbarplot(df,'Longitud Cadenas','Longitudes de grupos')

## Grafico de Varianza de cada tiempo

variance_by_group = D[['Grupo', 'Tiempo']].groupby('Grupo').var()['Tiempo']

variance_df = variance_by_group.to_frame().reset_index()
variance_df.rename(columns={'Tiempo': 'Tiempo_Variance'}, inplace=True)


plt.bar(variance_df['Grupo'],variance_df['Tiempo_Variance'])
plt.xticks(variance_df['Grupo'])
plt.show()

## Grafico de Error estandar

variance_by_group = D[['Grupo', 'Tiempo']].groupby('Grupo').std()['Tiempo']
variance_by_group_count = np.sqrt(D[['Grupo', 'Tiempo']].groupby('Grupo').count())['Tiempo']

result = variance_by_group / variance_by_group_count

variance_df = result.to_frame().reset_index()
variance_df.rename(columns={'Tiempo': 'Tiempo_Variance'}, inplace=True)



plt.bar(variance_df['Grupo'],variance_df['Tiempo_Variance'])
plt.xticks(variance_df['Grupo'])
plt.show()


## Grafico de Score y Tiempo

variance_by_group = D[['Grupo', 'Tiempo', "Score"]].groupby('Grupo').mean()

rects = plt.bar(variance_by_group['Score'],variance_by_group['Tiempo'],width=0.1)
plt.xticks(variance_by_group['Score'])
plt.show()

## Grafico de Score y Longitud

variance_by_group = D[['Grupo', 'Long', "Score"]].groupby('Grupo').mean()

rects = plt.bar(variance_by_group['Score'],variance_by_group['Long'],width=0.1)
plt.xticks(variance_by_group['Score'])
plt.show()

## Grafico de Score y Varianza

## Grafico de Varianza de cada tiempo

variance_by_group = D[['Grupo', 'Tiempo']].groupby('Grupo').var()['Tiempo']
Score_by_group = D[['Grupo', 'Score']].groupby('Grupo')

variance_df = variance_by_group.to_frame().reset_index()
variance_df.rename(columns={'Tiempo': 'Tiempo_Variance'}, inplace=True)

variance_df ["Score"] = Score_by_group

print(variance_df)


