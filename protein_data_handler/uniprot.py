from concurrent.futures import ThreadPoolExecutor, as_completed
from http.client import HTTPException
from urllib.parse import quote
import requests

import logging

from Bio import ExPASy, SwissProt
from sqlalchemy import func
from sqlalchemy.exc import NoResultFound

from protein_data_handler.helpers.parser.parser import (extract_float,
                                                        procesar_chain_string)

from protein_data_handler.models.uniprot import (
    Accession,
    GOTerm,
    PDBReference,
    Proteina, Chain
)

logger = logging.getLogger(__name__)


def cargar_codigos_acceso(criterio_busqueda, limite, session):
    """
    Busca en UniProtKB y carga códigos de acceso en la base de datos.

    Realiza una petición HTTP a UniProtKB usando criterios de búsqueda
    y un límite especificado. Los códigos de acceso obtenidos son
    almacenados en la base de datos mediante la sesión proporcionada.

    :param criterio_busqueda: Criterio para la consulta a UniProtKB.
    :type criterio_busqueda: str
    :param limite: Máximo de resultados de UniProtKB.
    :type limite: int
    :param session: Sesión de SQLAlchemy para la base de datos.
    :type session: sqlalchemy.orm.session.Session

    :raises Exception: Por errores en petición HTTP o al interactuar
                       con la base de datos.

    :return: None
    """

    criterio_busqueda_codificado = quote(criterio_busqueda)
    url = (
        f"https://rest.uniprot.org/uniprotkb/stream?"
        f"query={criterio_busqueda_codificado}"
        f"&format=list&size={limite}"
    )
    logger.info(f"URL solicitada: {url}")
    try:
        respuesta = requests.get(url)
        respuesta.raise_for_status()
        accessions = respuesta.text.strip().split("\n")
        logger.info(f"Número de Accessions en UniProt: {len(accessions)}")

        accessions_encontradas = []
        accessions_nuevas = []

        for accession_code in accessions:
            try:
                accession = (session.query(Accession)
                             .filter_by(accession_code=accession_code)
                             .one())

                accessions_encontradas.append(accession_code)

                # Actualiza la fecha de actualización
                accession.updated_at = func.now()
                accession.disappeared = False

            except NoResultFound:
                # Si no existe, crea una nueva instancia
                accession = Accession(accession_code=accession_code)
                session.add(accession)
                accessions_nuevas.append(accession)

        consulta_no_encontradas = (
            session.query(Accession)
            .filter(~Accession.accession_code.in_(accessions))
            .filter(not Accession.disappeared)
        )

        accessions_no_encontradas = consulta_no_encontradas.all()
        logger.info(f"Accesions no encontrados: "
                    f"{len(accessions_no_encontradas)}")

        consulta_no_encontradas.update({Accession.disappeared: True},
                                       synchronize_session=False)

        session.commit()

        # Loguear la información
        logger.info(f"Accessions ya existentes: {len(accessions_encontradas)}")
        logger.info(f"Accessions nuevas: {len(accessions_nuevas)}")

    except Exception as e:
        session.rollback()
        logger.error(f"Error: {e}")


def descargar_registro(accession_code):
    """
    Descarga información de una proteína desde ExPASy/SwissProt.

    :param accession_code: Código de acceso para la descarga.
    :type accession_code: str

    :raises Exception: Si hay errores en la descarga o procesamiento.

    :return: Registro de la proteína o None en caso de error.
    :rtype: SwissProt.Record o None
    """

    try:
        handle = ExPASy.get_sprot_raw(accession_code)
        record = SwissProt.read(handle)
        return record

    except Exception as e:
        logger.error(f"Error al descargar la entrada {accession_code}: {e}")
        return None


def extraer_entradas(session, max_workers=10):
    """
    Descarga y procesa entradas de UniProt concurrentemente con multihilos.

    Recupera códigos de acceso desde la base de datos y descarga información
    detallada de UniProt concurrentemente. Utiliza `ThreadPoolExecutor` para
    manejar múltiples descargas simultáneas.

    :param session: Sesión de SQLAlchemy para la base de datos.
    :type session: sqlalchemy.orm.session.Session
    :param max_workers: Número máximo de hilos para descargas.
    :type max_workers: int, opcional
    :raises Exception: Por fallos en operaciones de descarga o procesamiento.
    """

    logger.info("Iniciando la descarga de entradas de UniProt.")

    accessions = session.query(Accession).all()
    logger.info(f"Total de proteínas a descargar: {len(accessions)}")

    # Usar ThreadPoolExecutor para manejar múltiples descargas simultáneamente
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Crear un futuro para cada descarga de proteína
        future_to_uniprot_id = {
            executor.submit(
                descargar_registro, accession.accession_code
            ): accession
            for accession in accessions
        }

        # Procesar los resultados a medida que se completan
        for future in as_completed(future_to_uniprot_id):
            uniprot_id = future_to_uniprot_id[future]
            try:
                data = future.result()
                if data:
                    almacenar_entrada(data, session)
            except Exception as e:
                # Obtener la información de seguimiento de la excepción
                logger.error(f"Error al procesar la entrada {uniprot_id}: {e}")


def almacenar_entrada(data, session):
    """
    Almacena datos de UniProt en la base de datos.

    Toma un registro de UniProt y lo almacena en la base de datos.
    Actualiza entradas existentes con el mismo nombre o crea una nueva
    si no existe. También gestiona referencias cruzadas de cada
    entrada de UniProt.

    :param data: Datos de la entrada de UniProt.
    :type data: SwissProt.Record
    :param session: Sesión de SQLAlchemy para operaciones de base de datos.
    :type session: sqlalchemy.orm.session.Session
    :raises Exception: Por errores en operaciones de base de datos.
    """
    try:
        # Buscar si la proteína ya existe
        proteina = (
            session.query(Proteina)
            .filter_by(entry_name=data.entry_name)
            .first()
        )

        # Si no existe, crear una nueva instancia de Proteina
        if not proteina:
            proteina = Proteina(entry_name=data.entry_name)

        # Actualizar campos de la proteína con los datos proporcionados
        proteina.data_class = data.data_class
        proteina.molecule_type = data.molecule_type
        proteina.sequence_length = data.sequence_length
        proteina.created_date = data.created[0]
        proteina.sequence = data.sequence
        proteina.sequence_update_date = data.sequence_update[0]
        proteina.annotation_update_date = data.annotation_update[0]
        proteina.description = data.description
        proteina.gene_name = str(data.gene_name)
        proteina.organism = data.organism
        proteina.organelle = data.organelle
        proteina.organism_classification = ",".join(
            data.organism_classification
        )
        proteina.taxonomy_id = ",".join(data.taxonomy_id)
        proteina.host_organism = ",".join(data.host_organism)
        proteina.host_taxonomy_id = ",".join(data.host_taxonomy_id)
        proteina.comments = "; ".join(data.comments)
        proteina.keywords = data.keywords
        proteina.protein_existence = data.protein_existence
        proteina.seqinfo = data.seqinfo

        # Añadir la proteína a la sesión
        session.add(proteina)

        # Procesar cada código de acceso
        for accession_code in data.accessions:
            accession = (
                session.query(Accession)
                .filter_by(accession_code=accession_code)
                .first()
            )
            if accession is None:
                accession = Accession(accession_code=accession_code)
            accession.proteina_entry_name = proteina.entry_name
            session.add(accession)

        # Procesar cada referencia cruzada
        for reference in data.cross_references:
            if reference[0] == "PDB":
                pdb_ref = (
                    session.query(PDBReference)
                    .filter_by(pdb_id=reference[1])
                    .first()
                )

                if pdb_ref is None:
                    pdb_ref = PDBReference(
                        pdb_id=reference[1],
                        method=reference[2],
                        resolution=extract_float(reference[3]),
                        proteina=proteina,
                    )
                    session.add(pdb_ref)
                chains = reference[4].split(',')
                for chain_obj in chains:
                    chain_name, start, end = procesar_chain_string(chain_obj)

                    pdb_id = (session.query(PDBReference)
                              .filter_by(pdb_id=reference[1]).first().id)

                    chain = (
                        session.query(Chain)
                        .filter_by(pdb_reference_id=pdb_id, chain=chain_name)
                        .first()
                    )
                    if chain is None:
                        chain = Chain(pdb_reference_id=pdb_id,
                                      chain=chain_name,
                                      seq_start=start,
                                      seq_end=end)
                        session.add(chain)
            elif reference[0] == "GO":
                go_term = (
                    session.query(GOTerm).filter_by(go_id=reference[1]).first()
                )
                if go_term is None:
                    go_term = GOTerm(
                        go_id=reference[1],
                        category=reference[2].split(":")[0],
                        description=reference[2].split(":")[1],
                        proteina=proteina,
                    )
                    session.add(go_term)

        # Guarda todos los cambios en la base de datos
        session.commit()

    except HTTPException as e:
        logger.error(f"Error al volcar la entrada: {e}")
        session.rollback()
    finally:
        session.close()
