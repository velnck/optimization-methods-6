import copy

import numpy as np
from invert_matrix import calculate_inverted_matrix


# поиск оптимального плана симплекс-методом
def simplex_main_phase(a_matrix, b_basis_indexes, c, x):
    iteration = 1
    while iteration < 1000:
        # построить матрицу a_b_matrix и найти обратную ей
        a_b_matrix = a_matrix[:, b_basis_indexes]
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
        prev_a_b_inv = copy.deepcopy(a_b_inv)

        # построить вектор c_basis_vector
        c_basis_vector = np.array([c[i] for i in b_basis_indexes])

        # найти вектор потенциалов potentials_vector
        potentials_vector = c_basis_vector.dot(a_b_inv)

        # найти вектор оценок delta
        delta = potentials_vector.dot(a_matrix) - c

        # проверить условие оптимальности
        if (delta >= 0).all():
            return x

        # найти индекс первого отрицательного значения в векторе оценок j0
        j0 = -1
        for i in range(len(delta)):
            if delta[i] < 0:
                j0 = i
                break

        # найти вектор z
        z = a_b_inv.dot(a_matrix[:, j0])

        # построить вектор theta и найти в нем минимальное значение
        theta = []
        for i in range(len(a_matrix)):
            if z[i] > 0:
                theta.append(x[b_basis_indexes[i]] / z[i])
            else:
                theta.append(np.Infinity)
        theta_min = np.min(theta)

        # проверить условие ограниченности целевой функции
        if theta_min == np.Infinity:
            print("Целевой функционал неограничен сверху на множестве допустимых планов.")
            return None

        # найти первый индекс на котором достигается минимум theta j*
        k = np.argmin(theta)
        j_star = b_basis_indexes[k]

        # в b заменить индекс j* на j0
        b_basis_indexes[k] = j0

        # обновить компоненты плана
        x[j_star] = 0
        for i in range(len(b_basis_indexes)):
            x[b_basis_indexes[i]] -= theta_min * z[i]
        x[j0] = theta_min

        iteration += 1

    print("Превышено количество итераций.")
    return None
