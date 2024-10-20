from goatools.obo_parser import GODag
from goatools.godag.go_tasks import get_go2parents
from goatools.semantic import min_branch_length

# Cargar el archivo OBO con la estructura de la ontología GO
go_dag = GODag("/home/bioxaxi/PycharmProjects/protein-metamorphisms-is2/data/go-basic.obo")

# Definir los dos términos GO para los que se quiere calcular el MBL
term_1 = "GO:0006397"  # biological_process
term_2 = "GO:0048255"  # cellular process


print(min_branch_length(term_1, term_2,go_dag, branch_dist=None))