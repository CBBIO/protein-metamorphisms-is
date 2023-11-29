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
    >>> extract_float("-")
    None

    .. note:: Esta función solo extrae el primer número flotante encontrado
    en la cadena.
    """
    match = re.search(r"(\d+\.\d+)", s)
    if match:
        return float(match.group(1))
    else:
        return None
