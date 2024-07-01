from Bio import ExPASy, SwissProt


def download_record(accession_code):
    """
    Download detailed protein information from UniProt using ExPASy and SwissProt.
    ExPASy is a Bioinformatics Resource Portal which provides access to scientific databases and software tools,
    while SwissProt is a manually annotated and reviewed protein sequence database part of UniProt.
    Args:
        accession_code (str): The accession code of the protein record to download.
    Returns:
        record: A SwissProt record object containing detailed protein information.
    Raises:
        Exception: Raises an exception with a message indicating what went wrong during the download.
    """
    try:
        handle = ExPASy.get_sprot_raw(accession_code)
        record = SwissProt.read(handle)
        return record
    except Exception as e:
        raise Exception(f"Error downloading the entry {accession_code}: {e}")