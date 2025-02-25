Setup Instructions
==================

0. NVIDIA Requirements
-----------------------

Before proceeding with the installation, ensure your system meets these NVIDIA requirements:

Recommended Configuration
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **NVIDIA Driver**: Version 550.120 or later
- **CUDA Version**: 12.4 or later
- **CUDA Compiler (nvcc)**: 12.0 or later

Verify installation by running:

.. code-block:: bash

    nvidia-smi  # Check driver and CUDA version
    nvcc --version  # Check CUDA compiler version

If these commands fail, install or update:

1. NVIDIA drivers

2. CUDA Toolkit (includes nvcc)

For installation guide, visit `NVIDIA CUDA Installation Guide <https://docs.nvidia.com/cuda/cuda-installation-guide-linux/index.html>`_.

1. Install Docker
-----------------

Ensure Docker is installed on your system. If it’s not, you can download it from `here <https://docs.docker.com/get-docker/>`_.

2. Set Up PostgreSQL with pgvector
----------------------------------

Run the following command to start a PostgreSQL container with the pgvector extension:

.. code-block:: bash

    docker run -d --name pgvectorsql \
        -e POSTGRES_USER=usuario \
        -e POSTGRES_PASSWORD=clave \
        -e POSTGRES_DB=BioData \
        -p 5432:5432



Once the container is running, connect to the database and enable the `vector` extension:


.. code-block:: bash

        docker exec -it pgvectorsql psql -U usuario -d BioData -c "CREATE EXTENSION IF NOT EXISTS vector;"

3. (Optional) Connect to the Database
--------------------------------------

You can use **pgAdmin 4**, a graphical interface for managing and interacting with PostgreSQL databases, or any other SQL client. Use the connection details defined during the `docker run` command (`name`, `user`, `password`, `database`, `port`) to access the database.

For more information or to download pgAdmin 4, visit the official website: `pgAdmin 4 <https://www.pgadmin.org/>`_.

After connecting to the database, you can:

- **Verify the database status**:
  Check that the database is running and accepting connections.

- **Inspect tables and schema**:
  Explore the database structure, including tables, columns, and relationships.

- **Run queries**:
  Perform SQL operations such as inserting, updating, or retrieving data.

- **Monitor activity**:
  Trace database activity by reviewing logs and active processes.

3. Set Up RabbitMQ
------------------

Start a RabbitMQ container using the command below:

.. code-block:: bash

    docker run -d --name rabbitmq \
        -p 15672:15672 \
        -p 5672:5672 \
        rabbitmq:management

4. (Optional) Manage RabbitMQ
------------------------------

Once RabbitMQ is running, you can access its management interface to monitor and manage queues, exchanges, and messages.

Open your browser and go to: `RabbitMQ Management Interface <http://localhost:15672/#/queues>`_.

From this interface, you can:

- **Inspect queues**:
  View the list of queues, their message rates, and other details.

- **Monitor message flow**:
  Track the rate of incoming and outgoing messages in real-time.


Ensure that the RabbitMQ container is running and accessible at `localhost:15672`. Use the default credentials (`guest`/`guest`) unless you have configured different ones.

