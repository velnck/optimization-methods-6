import numpy as np
from invert_matrix import calculate_inverted_matrix


# поиск оптимального плана симплекс-методом
def simplex_method(a_matrix, b_basis_indexes, c, x):
    print("A matrix\n", a_matrix)
    print("basis indexes\t", b_basis_indexes)
    print("c vector\t", c)
    print("x vector\t", x)

    iteration = 1
    while iteration < 1000:
        print("\nITERATION ", iteration)
        # построить матрицу a_b_matrix и найти обратную ей
        a_b_matrix = np.array([
            [
                a_matrix[i][j] for j in b_basis_indexes
            ] for i in range(len(b_basis_indexes))
        ])
        print("Ab matrix\n", str(a_b_matrix))
        if iteration > 1:
            a_b_inv = np.array(calculate_inverted_matrix(prev_a_b_inv, a_matrix[:, j0], k))
            if a_b_inv is None:
                return None
        else:
            try:
                a_b_inv = np.linalg.inv(a_b_matrix)
            except np.linalg.LinAlgError:
                print("Матрица не обратима.")
                return None
            prev_a_b_inv = a_b_inv
        print("Ab inverse matrix\n", a_b_inv)

        # построить вектор c_basis_vector
        c_basis_vector = np.array([c[i] for i in b_basis_indexes])
        print("c_basis_vector\t", c_basis_vector)

        # найти вектор потенциалов potentials_vector
        potentials_vector = c_basis_vector.dot(a_b_inv)
        print("potentials_vector\t", potentials_vector)

        # найти вектор оценок delta
        delta = potentials_vector.dot(a_matrix) - c
        print("delta vector\t", delta)

        # проверить условие оптимальности
        if (delta >= 0).all():
            return x

        # найти индекс первого отрицательного значения в векторе оценок j0
        j0 = -1
        for i in range(len(delta)):
            if delta[i] < 0:
                j0 = i
                break

        print("j0 = ", j0)

        # найти вектор z
        z = a_b_inv.dot(a_matrix[:, j0])
        print("z vector\t", z)

        # построить вектор theta и найти в нем минимальное значение
        theta = []
        for i in range(len(a_matrix)):
            if z[i] > 0:
                theta.append(x[b_basis_indexes[i]] / z[i])
            else:
                theta.append(np.Infinity)
        print("theta vector\t", theta)
        theta_min = np.min(theta)
        print("theta min = ", theta_min)

        # проверить условие ограниченности целевой функции
        if theta_min == np.Infinity:
            print("Целевой функционал неограничен сверху на множестве допустимых планов.")
            return None

        # найти первый индекс на котором достигается минимум theta j*
        k = np.argmin(theta)
        j_star = b_basis_indexes[k]
        print("k = ", k)
        print("j* = ", j_star)

        # в b заменить индекс j* на j0
        b_basis_indexes[k] = j0
        print("b vector\t", b_basis_indexes)

        # обновить компоненты плана
        x[j_star] = 0

        for i in range(len(b_basis_indexes)):
            x[b_basis_indexes[i]] -= theta_min * z[i]

        x[j0] = theta_min
        print("x vector\t", x)

        iteration += 1

    print("Превышено количество итераций.")
    return None


def main():
    # входные данные
    # c = [1, 1, 0, 0, 0]
    # a = np.array([
    #         [-1, 1, 1, 0, 0],
    #         [1, 0, 0, 1, 0],
    #         [0, 1, 0, 0, 1]
    #      ])
    #
    # # начальный базисный допустимый план
    # b_basis_indexes = [2, 3, 4]
    # x = [0, 0, 1, 3, 2]

    c = [3, 5]
    a = np.array([
        [-1, 1, 1, 0, 0, 1, 0, 0],
        [1, 0, 0, 1, 0, 0, 1, 0],
        [0, 1, 0, 0, 1, 0, 0, 1]
    ])

    # начальный базисный допустимый план
    b_basis_indexes = [5, 6, 7]
    x = [0, 0, 0, 0, 0, 1, 3, 2]

    # ------------------------------------

    # c = np.random.randint(-5, high=5, size=5)
    # a = np.random.randint(-5, high=5, size=(3, 5))
    # b_basis_indexes = np.random.randint(0, high=5, size=3)
    # x = np.random.randint(-5, high=5, size=5)

    optimal_plan = simplex_method(a, b_basis_indexes, c, x)
    print(optimal_plan)


if __name__ == '__main__':
    main()
