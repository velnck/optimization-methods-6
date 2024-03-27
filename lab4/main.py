import copy

import numpy as np

from invert_matrix import calculate_inverted_matrix
from simplex_first_phase import simplex_first_phase


def dual_simplex(a_matrix, b, c, basis_indexes):
    n = len(a_matrix[0])
    m = len(a_matrix)

    iteration = 1
    while iteration < 1000:
        # print("\nITERATION", iteration)
        # построить матрицу a_basis_matrix и найти обратную ей
        a_basis_matrix = a_matrix[:, basis_indexes]
        if iteration > 1:
            # ???
            a_basis_inv = np.array(calculate_inverted_matrix(prev_a_basis_inv, a_matrix[:, j0], k))
            if a_basis_inv is None:
                return None
        else:
            try:
                a_basis_inv = np.linalg.inv(a_basis_matrix)
            except np.linalg.LinAlgError:
                print("Матрица не обратима.")
                return None
        prev_a_basis_inv = copy.deepcopy(a_basis_inv)
        # print("A_basis =\n", a_basis_matrix)
        # print("A_basis_inv =\n", a_basis_inv)

        # построить вектор c_basis_vector
        c_basis_vector = np.array([c[i] for i in basis_indexes])
        # print("c_basis_vector =", c_basis_vector)

        # найти базисный допустимый план двойственной задачи dual_basis_plan
        dual_basis_plan = c_basis_vector.dot(a_basis_inv)
        # print("dual_basis_plan =", dual_basis_plan)

        # найти псевдоплан pseudoplan
        pseudoplan = np.array([
            a_basis_inv[np.where(basis_indexes == i)[0]].dot(b)[0]
            if i in basis_indexes else 0 for i in range(n)
        ])
        # print("pseudoplan =", pseudoplan)

        # если псевдоплан >= 0, то он является оптимальным планом прямой задачи
        if np.all(pseudoplan >= 0):
            return pseudoplan

        # выделим отрицательную компоненту псевдоплана и сохраним её индекс в j_k
        j_k = np.where(pseudoplan < 0)[0][-1]
        k = np.where(basis_indexes == j_k)[0][0]
        # print("j_k =", j_k)
        # print("k =", k)

        # для каждого небазисного индекса вычислим mu
        is_infeasible = True
        mu = [None for _ in range(n)]
        for idx in range(n):
            if idx not in basis_indexes:
                mu[idx] = a_basis_inv[k].dot(a_matrix[:, idx])
                if mu[idx] < 0:
                    is_infeasible = False

        # print("mu =", mu)

        # если для всех небазисных индексов mu >= 0, то прямая задача не совместна
        if is_infeasible:
            print("Задача не совместна.")
            return None

        # для каждого небазисного индекса такого, что mu < 0, вычислим sigma
        # найдем минимальное значение sigma и сохраним его индекс в j0
        min_sigma = None
        min_sigma_index = None
        for i in range(n):
            if i not in basis_indexes:
                current_sigma = (c[i] - a_matrix[:, i].dot(dual_basis_plan)) / mu[i]
                # print("sigma[", i, "] =", current_sigma)
                if min_sigma is None or current_sigma < min_sigma:
                    min_sigma = current_sigma
                    min_sigma_index = i
        # print("min_sigma =", min_sigma)
        # print("min_sigma_index =", min_sigma_index)

        # в множестве базисных индексов заменим k-ый индекс на j0
        j0 = min_sigma_index
        basis_indexes[k] = min_sigma_index
        # print("basis_indexes =", basis_indexes)

        iteration += 1

    print("Превышено количество итераций.")
    return None



def main():
    c_vector = np.array([-4, -3, -7, 0, 0])
    A_matrix = np.array([
        [-2,-1, -4, 1, 0],
        [-2, -2, -2, 0, 1]
    ])
    b_vector = np.array([-1, -1.5])
    # ------------------------------------

    # c_vector = np.random.randint(-5, high=5, size=5)
    # A_matrix = np.random.randint(-5, high=5, size=(3, 5))
    # b_vector = np.random.randint(-5, high=5, size=3)
    # ------------------------------------

    print("Исходные данные:\nA =\n", A_matrix)
    print("b =", b_vector)
    print("c =", c_vector)
    initial_plan = simplex_first_phase(copy.deepcopy(A_matrix), copy.deepcopy(b_vector))
    if initial_plan is not None:
        print("Начальный базисный план:\nx = ", initial_plan[0], "\nB = ", initial_plan[1])
        optimal_plan = dual_simplex(A_matrix, b_vector, c_vector, initial_plan[1])
        # optimal_plan = dual_simplex(A_matrix, b_vector, c_vector, np.array([3, 4]))
        if optimal_plan is not None:
            print("Оптимальный план:\nx =", optimal_plan)


if __name__ == '__main__':
    main()
