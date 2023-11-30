import re


def extract_float(s):
    """
    Extrae el primer número flotante de una cadena de texto.

    Esta función busca en la cadena proporcionada el primer patrón que
    corresponda a un número flotante. Si lo encuentra, lo devuelve como
    un valor de tipo float. Si no encuentra ningún número flotante,
    devuelve None.

    :param s: La cadena de texto de la cual extraer el número flotante.
    :type s: str
    :return: El primer número flotante encontrado en la cadena, o None si no
             se encuentra ninguno.
    :rtype: float or None

    :Example:

    >>> extract_float("1.75 A")
    1.75
    >>> extract_float("El valor es -3.5")
    -3.5
    >>> extract_float("-")
    None

    .. note::
        Esta función solo extrae el primer número flotante encontrado
        en la cadena.
    """

    match = re.search(r"(\d+\.\d+)", s)
    if match:
        return float(match.group(1))
    else:
        return None


def safe_convert_to_int(value):
    """
    Intenta convertir un valor a un entero.

    Esta función intenta convertir el valor dado a un entero. Si la conversión
    falla debido a un ValueError (por ejemplo, si el valor es una cadena vacía
    o una cadena que no representa un entero), maneja la excepción retornando
    un valor por defecto.

    :param value: El valor que se intentará convertir a entero.
        Puede ser de cualquier tipo que sea convertible a entero,
        típicamente una cadena o un número.
    :type value: str o int
    :return: El valor convertido a entero, o None si la conversión falla.
    :rtype: int o None

    """
    try:
        return int(value)
    except ValueError:
        # Manejar el error o retornar un valor por defecto como None o 0
        # dependiendo de tus requerimientos
        return None


def procesar_chain_string(chain):
    """
    Procesa una cadena que representa un objeto de cadena y su rango.

    Esta función toma una cadena en el formato
    'nombreCadena=rangoInicio-rangoFin',
    la divide y convierte los rangos de inicio y fin a enteros.
    Retorna el nombre de la cadena y los valores de inicio y fin
    del rango como enteros.

    :param chain_obj: La cadena de objeto de cadena a procesar.
    :type chain_obj: str
    :return: El nombre de la cadena, el valor de inicio y
    el valor de fin del rango.

    """

    chain = chain.replace(" ", "")
    if '=' not in chain:
        return None, None, None

    chain_name, seq_range = chain.split('=')
    try:
        start, end = seq_range.split('-')
        start = safe_convert_to_int(start)
        end = safe_convert_to_int(end)
    except Exception:
        start, end = None, None

    return chain_name, start, end
