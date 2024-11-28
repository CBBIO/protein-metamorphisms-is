import os


def concatenate_files(file_list, output_file):
    """
    Concatena el contenido de los archivos en file_list en un solo archivo de salida.

    :param file_list: Lista de rutas de los archivos a concatenar.
    :param output_file: Ruta del archivo de salida.
    """
    with open(output_file, 'w') as outfile:
        for file_path in file_list:
            if os.path.isfile(file_path):
                with open(file_path, 'r') as infile:
                    outfile.write(infile.read())
                    outfile.write("\n")  # Añade una nueva línea entre archivos
            else:
                print(f"El archivo {file_path} no existe o no es un archivo válido.")


if __name__ == "__main__":
    # Lista de archivos a concatenar
    file_list = [
        '/home/bioxaxi/PycharmProjects/protein-metamorphisms-is2/protein_metamorphisms_is/sql/model/entities/sequence/sequence.py',
        '/home/bioxaxi/PycharmProjects/protein-metamorphisms-is2/protein_metamorphisms_is/sql/model/operational/clustering/cluster.py',
        '/home/bioxaxi/PycharmProjects/protein-metamorphisms-is2/protein_metamorphisms_is/sql/model/entities/protein/protein.py',
        '/home/bioxaxi/PycharmProjects/protein-metamorphisms-is2/protein_metamorphisms_is/sql/model/entities/protein/accesion.py',
        '/home/bioxaxi/PycharmProjects/protein-metamorphisms-is2/protein_metamorphisms_is/sql/model/entities/go_annotation/go_term.py',
        '/home/bioxaxi/PycharmProjects/protein-metamorphisms-is2/protein_metamorphisms_is/sql/model/entities/go_annotation/go_annotation.py',
        '/home/bioxaxi/PycharmProjects/protein-metamorphisms-is2/protein_metamorphisms_is/sql/model/entities/embedding/sequence_embedding.py',
        '/home/bioxaxi/PycharmProjects/protein-metamorphisms-is2/protein_metamorphisms_is/sql/model/entities/embedding/structure_3di.py',
        '/home/bioxaxi/PycharmProjects/protein-metamorphisms-is2/protein_metamorphisms_is/sql/model/entities/structure/chain.py',
        '/home/bioxaxi/PycharmProjects/protein-metamorphisms-is2/protein_metamorphisms_is/sql/model/entities/structure/state.py',
        '/home/bioxaxi/PycharmProjects/protein-metamorphisms-is2/protein_metamorphisms_is/sql/model/entities/structure/structure.py'



        # Añade más archivos según sea necesario
    ]

    # Ruta del archivo de salida
    output_file = ('/home/bioxaxi/PycharmProjects/protein-metamorphisms-is2/protein_metamorphisms_is/schema.txt')

    # Llama a la función para concatenar archivos
    concatenate_files(file_list, output_file)
    print(f"Archivos concatenados en {output_file}")
