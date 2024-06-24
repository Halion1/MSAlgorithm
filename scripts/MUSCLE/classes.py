# Import .py files

class Sequence:
    def __init__(self, id, data):
        self.id = id
        self.data = data


class DistanceMatrixOwn:
    def __init__(self, size):
        self.matrix = [[0.0 for _ in range(size)] for _ in range(size)]

    def set_distance(self, i, j, distance):
        self.matrix[i][j] = distance
        self.matrix[j][i] = distance  # symmetric matrix

    def get_distance(self, i, j):
        return self.matrix[i][j]
