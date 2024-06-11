import numpy as np


def solve(A, c, D, basis_plan, basis, extended_basis):
    iteration = 0
    while True:
        iteration += 1
        # находим векторы c(x), u(x), delta(x)
        c_x = c + np.dot(D, basis_plan)
        c_basis = c_x[basis]
        A_basis = A[:, basis]
        u_x = np.dot(-1 * c_basis, A_basis)
        delta_x = np.dot(u_x, A) + c_x

        # проверим условие оптимальности
        if np.all(delta_x >= 0):
            return basis_plan

        # выбираем отрицательную компоненту delta_x
        j0 = np.where(delta_x < 0)[0][0]

        # найдем вектор l_vector
        l_vector = np.zeros(len(D))
        l_vector[j0] = 1

        # для этого составим матрицу H
        H = np.array(D[np.ix_(extended_basis, extended_basis)])
        H = np.concatenate([H, A[:, extended_basis].transpose()], axis=1)
        H = np.concatenate([
            H,
            np.concatenate([
                A[:, extended_basis],
                np.zeros(shape=(len(A), len(A)))],
                axis=1)],
            axis=0)
        H_inv = np.linalg.inv(H)

        # строим вектор b*
        b_star = np.array(D[extended_basis, j0])
        b_star = np.concatenate([b_star, A[:, j0]])

        # получаем вектор x
        x_vector = np.dot(-1 * H_inv, b_star)
        l_vector[extended_basis] = x_vector[:len(extended_basis)]

        # для каждого индекса j из extended_basis найдем theta[j] и вычислим theta_j0
        delta = np.dot(np.dot(l_vector, D), l_vector)
        if delta == 0:
            theta_j0 = np.Infinity
        else:
            theta_j0 = abs(delta_x[j0]) / delta
        theta = np.empty(len(extended_basis))
        for i in range(len(extended_basis)):
            if l_vector[extended_basis[i]] < 0:
                theta[i] = -1 * basis_plan[extended_basis[i]] / l_vector[extended_basis[i]]
            else:
                theta[i] = np.Infinity

        # находим минимальное значение theta_min и индекс j_star, на котором оно достигается
        theta_min = min(theta)
        if theta_j0 < theta_min:
            theta_min = theta_j0
            j_star = j0
        else:
            j_star = np.argmin(theta)
        if theta_min == np.Infinity:
            print("Целевая функция задачи не ограничена снизу на множестве допустимых планов.")
            return None

        # обновим текущий план
        basis_plan = basis_plan + theta_min * l_vector

        # обновим опору ограничений basis и расширенную опору ограничений extended_basis
        if j_star == j0:
            extended_basis = np.concatenate([extended_basis, [j_star]])
        elif j_star in extended_basis and j_star not in basis:
            extended_basis = np.delete(extended_basis, np.where(extended_basis == j_star))
        elif j_star in basis:
            s = np.where(basis == j_star)[0][0]
            all_zeros = True
            for j in extended_basis:
                if j not in basis:
                    vect = np.dot(np.linalg.inv(A_basis), A[:, j])
                    if vect[s] != 0:
                        all_zeros = False
                        # в опоре ограничений заменяем j_star на j
                        basis[np.where(basis == j_star)[0][0]] = j
                        extended_basis = np.delete(extended_basis, np.where(extended_basis == j_star))
                        break
            if all_zeros:
                basis[np.where(basis == j_star)[0][0]] = j0
                extended_basis[np.where(extended_basis == j_star)[0][0]] = j0


if __name__ == '__main__':
    c = np.array([-8, -6, -4, -6])
    A = np.array([
        [1, 0, 2, 1],
        [0, 1, -1, 2]
    ])
    D = np.array([
        [2, 1, 1, 0],
        [1, 1, 0, 0],
        [1, 0, 1, 0],
        [0, 0, 0, 0]
    ])
    b = np.array([2, 3])
    basis = np.array([0, 1])
    extended_basis = np.array([0, 1])
    basis_plan = np.array([2, 3, 0, 0])

    print("Исходные данные:")
    print("Матрица A:\n", A)
    print("Вектор b: ", b)

    print("Вектор c:", c)
    print("Матрица D:\n", D)

    print("Правильный опорный план x:", basis_plan)

    optimal_plan = solve(A, c, D, basis_plan, basis, extended_basis)
    print("\nОптимальный план =\n", optimal_plan)
