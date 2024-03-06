from protein_metamorphisms_is.sql.model import StructuralComplexityLevel, StructuralAlignmentType


def handle_structural_complexity_levels(session, constants):
    # Cargar los niveles de complejidad desde el archivo YAML
    structural_complexity_levels = constants['structural_complexity_levels']

    for level_data in structural_complexity_levels:
        # Comprobar si el nivel de complejidad ya existe por nombre
        exists = session.query(StructuralComplexityLevel).filter_by(name=level_data['name']).first()
        if not exists:
            # Si no existe, crear y añadir el nuevo nivel de complejidad
            complexity_level = StructuralComplexityLevel(**level_data)
            session.add(complexity_level)
    session.commit()


def handle_structural_alignment_types(session, constants):
    # Cargar los tipos de alineamiento desde el archivo YAML
    structural_alignment_types = constants['structural_alignment_types']

    for level_data in structural_alignment_types:
        # Comprobar si el tipo de alineamiento ya existe por nombre
        exists = session.query(StructuralAlignmentType).filter_by(name=level_data['name']).first()
        if not exists:
            # Si no existe, crear y añadir el nuevo tipo de alineamiento
            alignment_type = StructuralAlignmentType(**level_data)
            session.add(alignment_type)

    # Comprometer los cambios en la base de datos
    session.commit()
