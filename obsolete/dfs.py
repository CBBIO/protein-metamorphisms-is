import requests
from requests.auth import HTTPBasicAuth
import json

def check_messages_in_memory(rabbitmq_host, rabbitmq_user, rabbitmq_password, queue_name):
    url = f'http://{rabbitmq_host}:15672/api/queues/%2F/{queue_name}'
    auth = HTTPBasicAuth(rabbitmq_user, rabbitmq_password)
    response = requests.get(url, auth=auth)
    if response.status_code == 200:
        queue_info = response.json()
        print(queue_info)
        # Aquí asumimos que la API devuelve una clave que indica cuántos mensajes están en memoria, esto es un ejemplo y deberás ajustarlo según la respuesta real de tu API.
        messages_in_memory = queue_info.get('messages_ram', 'No disponible')
        return messages_in_memory
    else:
        return f"Error al acceder a la API de RabbitMQ: {response.status_code}"

# Configuración
rabbitmq_host = "localhost"
rabbitmq_user = "guest"
rabbitmq_password = "guest"
queue_name = "pdbextractor_extractor"

# Uso de la función
print(check_messages_in_memory(rabbitmq_host, rabbitmq_user, rabbitmq_password, queue_name))
