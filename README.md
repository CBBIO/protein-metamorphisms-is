[![codecov](https://codecov.io/gh/frapercan/python-poetry-template/graph/badge.svg?token=hqzrADVeRy)](https://codecov.io/gh/frapercan/protein-data-handler)
[![PyPI version](https://badge.fury.io/py/protein-data-handler.svg)](https://badge.fury.io/py/protein-data-handler)
[![Documentation Status](https://readthedocs.org/projects/protein-data-handler/badge/?version=latest)](https://protein-data-handler.readthedocs.io/es/latest/?badge=latest)
![Linting Status](https://github.com/frapercan/protein-data-handler/actions/workflows/test-lint.yml/badge.svg?branch=main)

# Manejador de Datos de Proteínas

## Descripción
Este proyecto, **Manejador de Datos de Proteínas**, es una biblioteca de Python diseñada para facilitar el manejo y análisis de datos relacionados con proteínas. Ofrece herramientas para interactuar con bases de datos de proteínas, procesar datos y realizar diversas operaciones específicas del dominio.

## Características
- **Interacción con UniProt**: Módulos para consultar y procesar datos desde UniProt.
- **Herramientas de Base de Datos**: Funcionalidades para manejar bases de datos relacionadas con proteínas.
- **Configuración Flexible**: Soporte para configuraciones personalizadas a través de archivos YAML.

## Instalación
Para instalar esta biblioteca, ejecuta el siguiente comando en tu terminal:
```bash
pip install protein-data-handler
```

## Ejemplos
Podemos encontrar ejemplos de muestra bajo el directorio **`examples`**, configurar el directorio de trabajo donde se encuentra el script de muestra.

```python
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from protein_data_handler.helpers.config.yaml import read_yaml_config
from protein_data_handler.models.uniprot import Base
from protein_data_handler.uniprot import cargar_codigos_acceso, extraer_entradas
import logging
logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s - %(levelname)s - %(message)s')

if __name__ == "__main__":
    config = read_yaml_config("./config.yaml")
    DATABASE_URI = \
        (f"postgresql+psycopg2://{config['DB_USERNAME']}:"
         f"{config['DB_PASSWORD']}"
         f"@{config['DB_HOST']}:"
         f"{config['DB_PORT']}/"
         f"{config['DB_NAME']}")
    engine = create_engine(DATABASE_URI)
    Base.metadata.drop_all(engine)

    Base.metadata.create_all(engine)
    Session = sessionmaker(bind=engine)
    session = Session()

    cargar_codigos_acceso(criterio_busqueda=config['criterio_busqueda'],limite=config['limit'], session=session)
    extraer_entradas(session=session)
```

