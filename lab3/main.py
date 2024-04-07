import numpy as np
from simplex_main_phase import simplex_main_phase


def simplex_first_phase(A_matrix, b_vector):
    n = len(A_matrix[0])
    m = len(A_matrix)

    # 1. Преобразовать задачу таким образом, чтобы b_vector был неотрицательным.
    for i in range(len(b_vector)):
        if b_vector[i] < 0:
            b_vector[i] *= -1
            A_matrix[i] = A_matrix[i] * (-1)

    # 2. Составить вспомогательную задачу.
    c_additional = np.array([0 for _ in range(n)] + [-1 for _ in range(m)])

    A_matrix_additional = np.concatenate((A_matrix, np.eye(m)), axis=1)

    x_additional = np.concatenate((np.zeros(shape=n), b_vector))

    # 3. Построить начальный базисный допустимый план (x_additional, basis_indexes) вспомогательной задачи.
    basis_indexes = np.arange(m) + n

    # 4. Решить вспомогательную задачу основной фазой симплекс-метода. Решение: (x_additional, basis_indexes).
    if simplex_main_phase(A_matrix_additional, basis_indexes, c_additional, x_additional) is None:
        return None

    # 5. Проверить условие совместности задачи.
    if not np.all(x_additional[n:] == x_additional[n]):
        print("Задача несовместна.")
        return None

    # 6. Сформировать допустимый план исходной задачи.
    x_vector = x_additional[:n]

    iteration = 1
    while iteration < 1000:
        # 7. Проверить, является ли полученный результат решением задачи.
        if np.all(basis_indexes < n):
            return x_vector, basis_indexes

        # 8. Выбрать в наборе базисных индексов максимальный индекс искусственной переменной.
        max_index = np.max(basis_indexes)
        k = np.argmax(basis_indexes)
        i = max_index - n + 1

        # 9. Для каждого индекса index in {1, 2, ..., n} \ basis_indexes вычислить вектор l(index).
        A_matrix_additional_basis = A_matrix_additional[:, basis_indexes]
        A_matrix_additional_basis_inv = np.linalg.inv(A_matrix_additional_basis)
        constraints_are_linearly_dependent = True
        for index in range(n):
            if index not in basis_indexes:
                l = A_matrix_additional_basis_inv.dot(A_matrix_additional[:, index])
                # 10. Если найдется такой индекс, что l(index)_k != 0,
                # то заменим в наборе basis_indexes значение basis_indexes[k] на index.
                if l[k] != 0:
                    basis_indexes[k] = index
                    constraints_are_linearly_dependent = False
                    break

        # 11. Если таких индексов нет, то i-ое основное ограничение задачи
        # линейно выражается через остальные и его необходимо удалить.
        if constraints_are_linearly_dependent:
            A_matrix = np.delete(A_matrix, obj=i - 1, axis=0)
            b_vector = np.delete(b_vector, obj=i - 1, axis=0)
            basis_indexes = np.delete(basis_indexes, obj=k, axis=0)
            A_matrix_additional = np.delete(A_matrix_additional, obj=i - 1, axis=0)

        iteration += 1

    print("Превышено количество итераций.")
    return None


def main():
    c_vector = np.array([1, 0, 0])
    A_matrix = np.array([
        [1, 1, 1],
        [2, 2, 2]
    ])
    b_vector = np.array([0, 0])


    print("Исходные данные:\nc =", c_vector)
    print("A = \n", A_matrix)
    print("b =", b_vector)
    result = simplex_first_phase(A_matrix, b_vector)
    if result is not None:
        print("Решение:\nx = ", result[0], "\nB = ", result[1])


if __name__ == '__main__':
    main()
