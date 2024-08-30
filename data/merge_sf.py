import pandas as pd

# Cargar los dos archivos CSV
try:
    data = pd.read_csv("/home/bioxaxi/PycharmProjects/protein-metamorphisms-is/data/SF.71pairs.regions.checked.ago.2024.csv")
    data_2 = pd.read_csv("/home/bioxaxi/PycharmProjects/protein-metamorphisms-is/data/metamorphic.regions.checked.ago.2024.csv")
    print("Archivos CSV cargados correctamente.")
except Exception as e:
    print(f"Error al cargar los archivos CSV: {e}")

# Descomponer pairA en pdb_id y chain en 'data'
try:
    data[['pdb_id', 'chain']] = data['pairA'].str.split('_', expand=True)
    data[['pdb_idB', 'chainB']] = data['pairB'].str.split('_', expand=True)
    data.drop(columns=['pairA', 'pairB', 'pdb_idB', 'chainB'], inplace=True)  # Eliminamos columnas innecesarias
    print("Columnas 'pairA' y 'pairB' descompuestas y eliminadas correctamente en 'data'.")
except Exception as e:
    print(f"Error al procesar 'pairA' y 'pairB' en 'data': {e}")

# Asignar correctamente pairA a pdb_id y pairB a pdb_idB en data_2 y luego descomponer
try:
    data_2['pdb_id'], data_2['chain'] = data_2['pairA'].str[:4], data_2['pairA'].str[4:]
    data_2['pdb_idB'], data_2['chainB'] = data_2['pairB'].str[:4], data_2['pairB'].str[4:]
    data_2.drop(columns=['pairA', 'pairB'], inplace=True)  # Eliminamos columnas innecesarias
    print("Columnas 'pairA' y 'pairB' descompuestas correctamente en 'data_2'.")
except Exception as e:
    print(f"Error al descomponer 'pairA' y 'pairB' en 'data_2': {e}")

# Realizar el merge de los datasets usando 'pdb_id' y 'chain'
try:
    merged_data = pd.merge(
        data,
        data_2,
        how='inner',
        on=['pdb_id', 'chain']  # Especifica las columnas clave comunes entre los dos datasets
    )
    print("Merge realizado correctamente.")
except Exception as e:
    print(f"Error al realizar el merge: {e}")

# Guardar el dataset mergeado en un nuevo archivo CSV
try:
    merged_data.to_csv("/home/bioxaxi/PycharmProjects/protein-metamorphisms-is/data/merged_sf_metamorphic_regions.csv", index=False)
    print("Dataset mergeado guardado correctamente.")
except Exception as e:
    print(f"Error al guardar el dataset mergeado: {e}")

# Mostrar las primeras filas del dataset mergeado
try:
    print("Primeras filas de 'merged_data':")
    print(merged_data.head())
except Exception as e:
    print(f"Error al mostrar 'merged_data': {e}")
