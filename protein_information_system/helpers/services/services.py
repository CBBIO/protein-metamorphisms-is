import pika
from psycopg2 import OperationalError
from sqlalchemy import text

from protein_information_system.sql.base.database_manager import DatabaseManager


def check_services(conf, logger):
    try:
        db_manager = DatabaseManager(conf)
        session = db_manager.get_session()
        session.execute(text("SELECT 1"))  # ✅ Correcto
        session.close()
        logger.info("PostgreSQL connection OK.")
    except OperationalError as e:
        raise ConnectionError(
            f"Could not connect to PostgreSQL at {conf['DB_HOST']}:{conf['DB_PORT']}. "
            f"Verify your DB settings.\nDocs: https://fantasia.readthedocs.io"
        ) from e

    # Comprobación de RabbitMQ usando Pika (como haces tú mismo)
    try:
        connection_params = pika.ConnectionParameters(
            host=conf["rabbitmq_host"],
            port=conf.get("rabbitmq_port", 5672),
            credentials=pika.PlainCredentials(
                conf["rabbitmq_user"], conf["rabbitmq_password"]
            ),
            blocked_connection_timeout=3,
            heartbeat=600,
        )
        connection = pika.BlockingConnection(connection_params)
        connection.close()
        logger.info("RabbitMQ connection OK.")
    except Exception as e:
        raise ConnectionError(
            f"Could not connect to RabbitMQ at {conf['rabbitmq_host']}:{conf.get('rabbitmq_port', 5672)}. "
            f"Verify your MQ settings.\nDocs: https://fantasia.readthedocs.io"
        ) from e
