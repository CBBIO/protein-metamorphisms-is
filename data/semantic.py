from goatools.obo_parser import GODag
from collections import deque

def calculate_shortest_path(go_dag, term1, term2):
    """Calcula la distancia mínima entre dos términos GO en el DAG."""
    def bfs_shortest_path(start, goal):
        """Encuentra la distancia más corta entre dos términos en un grafo usando BFS."""
        queue = deque([(start, 0)])  # Cola para BFS con (nodo_actual, distancia_actual)
        visited = set()  # Para evitar ciclos

        while queue:
            current, dist = queue.popleft()
            if current == goal:
                return dist
            if current in visited:
                continue
            visited.add(current)

            # Añadir padres (movimiento hacia ancestros en el DAG)
            for parent in go_dag[current].parents:
                queue.append((parent.id, dist + 1))

        return float('inf')  # Si no hay camino entre los nodos

    # Ejecutar BFS desde term1 a term2 y viceversa (bidireccionalidad implícita en el DAG)
    distance1 = bfs_shortest_path(term1, term2)
    distance2 = bfs_shortest_path(term2, term1)

    return min(distance1, distance2)

# Cargar la ontología GO
go_dag = GODag("go-basic.obo")

# Términos GO
term1 = "GO:0045944"
term2 = "GO:0042789"

# Calcular y mostrar la distancia mínima
distance = calculate_shortest_path(go_dag, term1, term2)
if distance == float('inf'):
    print(f"No hay un camino entre {term1} y {term2}.")
else:
    print(f"La distancia mínima entre {term1} y {term2} es: {distance}")
