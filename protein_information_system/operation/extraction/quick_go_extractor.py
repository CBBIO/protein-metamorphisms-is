import traceback
from typing import Any, Optional

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

from protein_information_system.sql.model.entities.go_annotation.go_annotation import (ProteinGOTermAnnotation,
)
from protein_information_system.sql.model.entities.go_annotation.go_term import GOTerm
from protein_information_system.sql.model.entities.protein.accesion import Accession
from protein_information_system.tasks.queue import QueueTaskInitializer


class QuickGOAnnotationExtractor(QueueTaskInitializer):
    """
    Extracts Gene Ontology annotations from the QuickGO REST API.

    The extractor obtains, for every UniProt accession:

    - UniProt gene product identifier.
    - GO term identifier.
    - GO term name.
    - GO aspect/category.
    - GO evidence code, such as EXP, IDA or IEA.
    - Supporting PubMed identifier, when available.

    Accessions must already exist and be linked to proteins before
    executing this extractor.
    """

    BASE_URL = "https://www.ebi.ac.uk/QuickGO/services"
    ANNOTATION_ENDPOINT = f"{BASE_URL}/annotation/search"

    def __init__(self, conf):
        """
        Initializes the QuickGO extractor.

        Args:
            conf: Configuration dictionary.
        """
        super().__init__(conf, session_required=True)

        self.reference_attribute = "Accession"

        self.page_size = conf.get("quickgo_page_size", 100)
        self.request_timeout = conf.get("quickgo_timeout", 60)

        self.http_session = self._create_http_session()

        self.logger.info(
            "QuickGOAnnotationExtractor initialized with page size %s.",
            self.page_size,
        )

    def _create_http_session(self) -> requests.Session:
        """
        Creates an HTTP session with automatic retry handling.
        """
        retry_strategy = Retry(
            total=5,
            connect=5,
            read=5,
            status=5,
            backoff_factor=1.0,
            status_forcelist=(429, 500, 502, 503, 504),
            allowed_methods=frozenset({"GET"}),
            respect_retry_after_header=True,
        )

        adapter = HTTPAdapter(
            max_retries=retry_strategy,
            pool_connections=10,
            pool_maxsize=10,
        )

        session = requests.Session()
        session.mount("https://", adapter)

        session.headers.update(
            {
                "Accept": "application/json",
                "User-Agent": "protein-information-system/quickgo-extractor",
            }
        )

        return session

    def enqueue(self):
        """
        Enqueues one task for every accession stored in the database.

        Accessions must already be linked to proteins by UniProtExtractor.
        """
        try:
            query = (
                self.session.query(Accession)
                .filter(Accession.protein_id.isnot(None))
                .order_by(Accession.code)
            )

            limit_execution = self.conf.get("limit_execution")

            if isinstance(limit_execution, int) and limit_execution > 0:
                query = query.limit(limit_execution)

            accessions = query.all()

            if not accessions:
                self.logger.info(
                    "No accessions linked to proteins were found. "
                    "QuickGO enqueue operation skipped."
                )
                return

            for accession in accessions:
                self.logger.debug(
                    "Publishing QuickGO task for accession '%s'.",
                    accession.code,
                )

                self.publish_task(accession.code)

            self.logger.info(
                "Successfully enqueued %s QuickGO tasks.",
                len(accessions),
            )

        except Exception as exc:
            self.logger.error(
                "Failed to enqueue QuickGO tasks: %s\n%s",
                exc,
                traceback.format_exc(),
            )
            raise

    def process(self, accession_code: str) -> dict[str, Any]:
        """
        Downloads all QuickGO annotations for one UniProt accession.

        Args:
            accession_code:
                UniProt accession, for example P04637 or UniProtKB:P04637.

        Returns:
            Dictionary containing the accession and its QuickGO annotations.
        """
        if not accession_code:
            raise ValueError("No accession code provided.")

        normalized_accession = self._normalize_accession(accession_code)
        gene_product_id = f"UniProtKB:{normalized_accession}"

        try:
            annotations = self._fetch_all_annotations(gene_product_id)

            self.logger.info(
                "Retrieved %s QuickGO annotations for accession '%s'.",
                len(annotations),
                normalized_accession,
            )

            return {
                "accession_code": normalized_accession,
                "gene_product_id": gene_product_id,
                "annotations": annotations,
            }

        except requests.RequestException as exc:
            self.logger.error(
                "QuickGO request failed for accession '%s': %s\n%s",
                normalized_accession,
                exc,
                traceback.format_exc(),
            )
            raise

        except Exception as exc:
            self.logger.error(
                "Unexpected QuickGO error for accession '%s': %s\n%s",
                normalized_accession,
                exc,
                traceback.format_exc(),
            )
            raise

    def _fetch_all_annotations(
        self,
        gene_product_id: str,
    ) -> list[dict[str, Any]]:
        """
        Retrieves every annotation page for one gene product.
        """
        annotations = []
        page = 1
        number_of_hits = None

        while True:
            params = {
                "geneProductId": gene_product_id,
                "includeFields": "goName",
                "limit": self.page_size,
                "page": page,
            }

            response = self.http_session.get(
                self.ANNOTATION_ENDPOINT,
                params=params,
                timeout=self.request_timeout,
            )

            response.raise_for_status()

            payload = response.json()
            page_results = payload.get("results", [])

            if number_of_hits is None:
                number_of_hits = payload.get("numberOfHits", 0)

            if not page_results:
                break

            annotations.extend(page_results)

            self.logger.debug(
                "QuickGO page %s retrieved for %s: %s records.",
                page,
                gene_product_id,
                len(page_results),
            )

            if len(annotations) >= number_of_hits:
                break

            page += 1

        return annotations

    def store_entry(self, data: dict[str, Any]) -> None:
        """
        Stores QuickGO annotations supported by a PubMed reference.

        Annotations whose reference is not a PMID are skipped and are not
        inserted into ProteinGOTermAnnotation.
        """
        accession_code = data["accession_code"]
        annotations = data.get("annotations", [])

        try:
            accession = (
                self.session.query(Accession)
                .filter(Accession.code == accession_code)
                .first()
            )

            if accession is None:
                raise ValueError(
                    f"Accession '{accession_code}' does not exist "
                    "in the database."
                )

            if accession.protein_id is None:
                raise ValueError(
                    f"Accession '{accession_code}' is not linked "
                    "to a protein."
                )

            protein_id = accession.protein_id

            created = 0
            existing = 0
            skipped_invalid = 0
            skipped_filtered = 0
            skipped_duplicate = 0

            processed_keys: set[tuple[str, str, Optional]] = set()

            for annotation_data in annotations:
                if not self._is_allowed_annotation(annotation_data):
                    skipped_filtered += 1
                    continue

                go_id = annotation_data.get("goId")
                evidence_code = annotation_data.get("goEvidence")

                if go_id is not None:
                    go_id = str(go_id).strip()

                if evidence_code is not None:
                    evidence_code = str(evidence_code).strip()

                if not go_id or not evidence_code:
                    self.logger.warning(
                        "Skipping incomplete QuickGO annotation: %s.",
                        annotation_data,
                    )
                    skipped_invalid += 1
                    continue

                publication_id = self.extract_publication_id(
                    annotation_data.get("reference")
                )

                association_key = (
                    go_id,
                    evidence_code,
                    publication_id,
                )

                if association_key in processed_keys:
                    skipped_duplicate += 1
                    continue

                processed_keys.add(association_key)

                go_term = self.get_or_create_go_term(
                    annotation_data
                )

                status = self.get_or_create_association(
                    protein_id=protein_id,
                    go_id=go_term.go_id,
                    annotation_data=annotation_data,
                )

                if status == "created":
                    created += 1
                elif status == "existing":
                    existing += 1
                else:
                    skipped_invalid += 1

            self.session.commit()

            self.logger.info(
                "Stored QuickGO annotations for '%s': "
                "%s created, %s existing, %s duplicated, "
                "%s filtered and %s invalid.",
                accession_code,
                created,
                existing,
                skipped_duplicate,
                skipped_filtered,
                skipped_invalid,
            )

        except Exception as exc:
            self.session.rollback()

            self.logger.error(
                "Failed to store QuickGO annotations for '%s': %s\n%s",
                accession_code,
                exc,
                traceback.format_exc(),
            )
            raise

    def _is_allowed_annotation(
        self,
        annotation_data: dict[str, Any],
    ) -> bool:
        """
        Checks whether a QuickGO annotation passes the configured
        evidence filters.

        QuickGO has two related fields:

        - goEvidence: EXP, IDA, IEA, IMP...
        - evidenceCode: ECO:0000269, ECO:0000314...
        """
        allowed_evidences = self.conf.get("allowed_evidences")
        allowed_eco_codes = self.conf.get("allowed_eco_codes")

        go_evidence = annotation_data.get("goEvidence")
        eco_code = annotation_data.get("evidenceCode")

        if allowed_evidences:
            if go_evidence not in allowed_evidences:
                return False

        if allowed_eco_codes:
            if eco_code not in allowed_eco_codes:
                return False

        return True

    def get_or_create_go_term(
        self,
        annotation_data: dict[str, Any],
    ) -> GOTerm:
        """
        Retrieves or creates a GO term from QuickGO annotation data.
        """
        go_id = annotation_data.get("goId")

        if not go_id:
            raise ValueError(
                f"QuickGO annotation does not contain goId: "
                f"{annotation_data}"
            )

        go_term = (
            self.session.query(GOTerm)
            .filter(GOTerm.go_id == go_id)
            .first()
        )

        go_name = annotation_data.get("goName")
        go_aspect = annotation_data.get("goAspect")

        if go_term is None:
            go_term = GOTerm(
                go_id=go_id,
                category=go_aspect,
                description=go_name,
            )

            self.session.add(go_term)

            self.logger.debug(
                "Created GO term '%s'.",
                go_id,
            )

        else:
            if not go_term.description and go_name:
                go_term.description = go_name

            if not go_term.category and go_aspect:
                go_term.category = go_aspect

        return go_term

    

    @staticmethod
    def extract_publication_id(reference: Any) -> Optional[str]:
        if reference is None:
            return None

        if isinstance(reference, list):
            for item in reference:
                publication_id = (
                    QuickGOAnnotationExtractor.extract_publication_id(item)
                )
                if publication_id is not None:
                    return publication_id
            return None

        reference = str(reference).strip()

        if not reference.upper().startswith("PMID:"):
            return None

        publication_id = reference.split(":", maxsplit=1)[1].strip()
        return publication_id or None

    @staticmethod
    def _normalize_accession(accession_code: str) -> str:
        """
        Removes the UniProtKB prefix when it is present.
        """
        accession_code = accession_code.strip()

        if accession_code.startswith("UniProtKB:"):
            return accession_code.split(":", maxsplit=1)[1]

        return accession_code






    def get_or_create_association(
            self,
            protein_id: str,
            go_id: str,
            annotation_data: dict[str, Any],
            ) -> str:
        
        evidence_code = annotation_data.get("goEvidence")

        if evidence_code is not None:
            evidence_code = str(evidence_code).strip()

        if not evidence_code:
            self.logger.warning(
                "Skipping annotation without GO evidence: "
                "protein=%s, GO=%s.",
                protein_id,
                go_id,
            )
            return "skipped"

        publication_id = self.extract_publication_id(
            annotation_data.get("reference")
        )

        try:
            existing_association = (
                self.session.query(ProteinGOTermAnnotation)
                .filter_by(
                    protein_id=protein_id,
                    go_id=go_id,
                    evidence_code=evidence_code,
                    publication_id=publication_id,
                )
                .first()
            )

            if existing_association is not None:
                self.logger.debug(
                    "QuickGO association already exists: "
                    "protein=%s, GO=%s, evidence=%s, publication=%s.",
                    protein_id,
                    go_id,
                    evidence_code,
                    publication_id,
                )
                return "existing"

            association = ProteinGOTermAnnotation(
                protein_id=protein_id,
                go_id=go_id,
                evidence_code=evidence_code,
                publication_id=publication_id,
            )

            self.session.add(association)

            self.logger.debug(
                "Created QuickGO association: "
                "protein=%s, GO=%s, evidence=%s, publication=%s.",
                protein_id,
                go_id,
                evidence_code,
                publication_id,
            )

            return "created"

        except Exception as exc:
            self.logger.error(
                "Failed to create QuickGO association: "
                "protein=%s, GO=%s, evidence=%s, "
                "publication=%s, error=%s.",
                protein_id,
                go_id,
                evidence_code,
                publication_id,
                exc,
                exc_info=True,
            )
            raise