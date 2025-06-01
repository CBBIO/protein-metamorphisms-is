import yaml


def read_yaml_config(filepath):
    """
    Lee un archivo YAML y devuelve su contenido como un diccionario.

    :param filepath: Ruta al archivo YAML.
    :return: Diccionario con la configuración.
    """
    with open(filepath, "r") as file:
        config = yaml.safe_load(file)
    return config
