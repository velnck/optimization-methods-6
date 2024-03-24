# multiply matrix by vector
def multiply_matrix_by_vector(matrix, vector):
    size = len(matrix)
    result_vector = []
    for i in range(size):
        result_vector.append(sum(matrix[i][j] * vector[j] for j in range(size)))
    return result_vector


# calculate q matrix
def calculate_q_matrix(inverted_a, vector_x, index):
    vector_l = multiply_matrix_by_vector(inverted_a, vector_x)
    try:
        coefficient = (-1) / vector_l[index]
    except ZeroDivisionError:
        print('Matrix cannot be inverted.')
        return None
    vector_l[index] = -1
    vector_l = [value * coefficient for value in vector_l]

    q_matrix = []
    size = len(inverted_a)
    for row_idx in range(size):
        q_row = []
        for col_idx in range(size):
            if col_idx == index:
                q_row.append(vector_l[row_idx])
            elif row_idx == col_idx:
                q_row.append(1)
            else:
                q_row.append(0)
        q_matrix.append(q_row)
    return q_matrix


# multiply q matrix by inverted a matrix
def multiply_q_by_inverted_a(matrix_q, matrix_inverted_a, index):
    result_matrix = []

    size = len(matrix_q)
    for row_idx in range(size):
        result_matrix_row = []
        if row_idx == index:
            for col_idx in range(size):
                result_matrix_row.append(
                        matrix_q[index][index] * matrix_inverted_a[index][col_idx])
        else:
            for col_idx in range(size):
                result_matrix_row.append((
                        matrix_inverted_a[row_idx][col_idx] +
                        matrix_q[row_idx][index] * matrix_inverted_a[index][col_idx]))
        result_matrix.append(result_matrix_row)
    return result_matrix


# compound function that calculates final inverted matrix
def calculate_inverted_matrix(inverted_matrix_a, vector_x, position):
    q_matrix = calculate_q_matrix(inverted_matrix_a, vector_x, position)
    if q_matrix is not None:
        result_matrix = multiply_q_by_inverted_a(q_matrix, inverted_matrix_a, position)
        return result_matrix
    else:
        return None


# test if calculated matrix is correct
# by multiplying it by the initial matrix
def test_inverted_matrix(inverted_matrix, initial_matrix):
    size = len(inverted_matrix)
    for i in range(size):
        for j in range(size):
            el = 0
            for k in range(size):
                el += inverted_matrix[i][k] * initial_matrix[k][j]
            if i == j:
                if (abs(el) - 1) > 1e-10:
                    return False
            elif abs(el) > 1e-10:
                return False
    return True
