import copy

import numpy as np


def get_edges(vertices):
    # по заданному списку вершин графов найдем все ребра
    edges = dict()
    for i in range(len(vertices)):
        edges[vertices[i]] = []
        for j in range(len(vertices)):
            if i == j:
                continue
            if vertices[i][0] == vertices[j][0]:
                no_vertices_between = True
                for col_num in range(min(vertices[i][1], vertices[j][1]) + 1, max(vertices[i][1], vertices[j][1])):
                    if (vertices[i][0], col_num) in vertices:
                        no_vertices_between = False
                        break
                if no_vertices_between:
                    edges[vertices[i]].append(vertices[j])
            elif vertices[i][1] == vertices[j][1]:
                no_vertices_between = True
                for row_num in range(min(vertices[i][0], vertices[j][0]) + 1, max(vertices[i][0], vertices[j][0])):
                    if (row_num, vertices[i][1]) in vertices:
                        no_vertices_between = False
                        break
                if no_vertices_between:
                    edges[vertices[i]].append(vertices[j])
    return edges


def find_cycle(vertices):
    edges = get_edges(vertices)
    current_vertex = vertices[0]
    path = []
    while True:
        # добавляем вершину к текущему пути
        path.append(current_vertex)
        # выбираем новую текущую вершину
        current_vertex = edges[current_vertex][0]
        # удаляем из списка ребро, соединяющее текущую вершину с предыдущей
        if path[-1] in edges[current_vertex]:
            edges[current_vertex].remove(path[-1])
        if current_vertex in edges[path[-1]]:
            edges[path[-1]].remove(current_vertex)

        if current_vertex in path[:-1]:
            # если добавленная вершина уже встречалась ранее, то цикл найден
            return path[path.index(current_vertex):]

        while len(edges[current_vertex]) == 0:
            # если в текущей вершине не осталось ребер, то возвращаемся к предыдущей вершине
            current_vertex = path.pop()


def potentials_method(production_vector, consumption_vector, cost_matrix):
    # проверим условие баланса
    total_production = sum(production_vector)
    total_consumption = sum(consumption_vector)
    if total_production > total_consumption:
        # если производится больше, чем потребляется, то вводим новый пункт потребления
        consumption_vector = np.append(consumption_vector, [total_production - total_consumption])
        cost_matrix = np.concatenate((cost_matrix, np.zeros(shape=(len(cost_matrix), 1))), axis=1)
    elif total_production < total_consumption:
        # если потребляется больше, чем производится, то создадим новый пункт производства
        production_vector = np.append(production_vector, [total_consumption - total_production])
        cost_matrix = np.concatenate((cost_matrix, np.zeros(shape=(1, len(cost_matrix[0])))), axis=0)

    production_points_count = len(production_vector)
    consumption_points_count = len(consumption_vector)
    transport_plan = np.zeros((production_points_count, consumption_points_count))   # план перевозок
    basis = []

    # сформируем начальное множество базисных позиций
    production_point = 0
    consumption_point = 0
    while not (np.all(production_vector == 0) and np.all(consumption_vector == 0)):
        basis.append((production_point, consumption_point))
        if production_vector[production_point] > consumption_vector[consumption_point]:
            # количество продукции в текущем пункте производства больше чем заявка в пункте потребления
            transport_plan[production_point][consumption_point] = consumption_vector[consumption_point]
            production_vector[production_point] -= consumption_vector[consumption_point]
            consumption_vector[consumption_point] = 0
            consumption_point += 1   # переход к следующему пункту потребления
        elif production_vector[production_point] <= consumption_vector[consumption_point]:
            # количество продукции в текущем пункте производства меньше чем заявка в пункте потребления
            transport_plan[production_point][consumption_point] = production_vector[production_point]
            consumption_vector[consumption_point] -= production_vector[production_point]
            production_vector[production_point] = 0
            production_point += 1   # переход к следующему пункту производства

    # начало итерации
    iteration = 0
    while True:
        # составим систему уравнений для потенциалов
        potentials_equations_count = len(basis) + 1
        potentials_equations_vector = np.zeros(shape=potentials_equations_count)

        production_potentials_equations_matrix = np.zeros((potentials_equations_count, production_points_count))
        consumption_potentials_equations_matrix = np.zeros((potentials_equations_count, consumption_points_count))
        for i in range(len(basis)):
            production_potentials_equations_matrix[i][basis[i][0]] = 1
            consumption_potentials_equations_matrix[i][basis[i][1]] = 1
            potentials_equations_vector[i] = cost_matrix[basis[i][0]][basis[i][1]]
        production_potentials_equations_matrix[-1][0] = 1
        potentials_equations_vector[-1] = 0
        # найдем решение полученной системы
        potentials = np.linalg.solve(
            np.concatenate(
                (production_potentials_equations_matrix, consumption_potentials_equations_matrix),
                axis=1
            ),
            potentials_equations_vector
        )
        production_potentials = potentials[:production_points_count]
        consumption_potentials = potentials[production_points_count:]

        # проверим условие оптимальности текущего плана
        is_optimal = True
        new_basis_position = (-1, -1)
        for i in range(production_points_count):
            for j in range(consumption_points_count):
                if (i, j) not in basis:
                    is_optimal = production_potentials[i] + consumption_potentials[j] <= cost_matrix[i][j]
                    if not is_optimal:
                        new_basis_position = (i, j)
                        break
            if not is_optimal:
                break
        if is_optimal:
            return transport_plan

        # если план не является оптимальным
        # добавляем в базис позицию, для которой не выполняется условие оптимальности
        basis.append(new_basis_position)
        # находим цикл в графе базиса
        cycle = find_cycle(copy.deepcopy(basis))
        cycle = cycle[cycle.index(new_basis_position):] + cycle[:cycle.index(new_basis_position)]

        # находим угловые вершины
        corner_positions = []
        for i in range(len(cycle)):
            if cycle[i - 1][0] != cycle[(i + 1) % len(cycle)][0] and cycle[i - 1][1] != cycle[(i + 1) % len(cycle)][1]:
                corner_positions.append(cycle[i])

        # находим минимальное значение в угловых вершинах со знаком '-'
        min_corner_value = np.Infinity
        for i in range(1, len(corner_positions), 2):
            min_corner_value = min(min_corner_value, transport_plan[corner_positions[i][0]][corner_positions[i][1]])
        # удаляем из базиса позицию, на которой достигается найденное минимальное значение
        basis.remove(next((pos for pos in basis if transport_plan[pos[0]][pos[1]] == min_corner_value)))

        # обновляем план перевозок
        for i in range(len(corner_positions)):
            if i % 2 == 0:
                transport_plan[corner_positions[i][0]][corner_positions[i][1]] += min_corner_value
            else:
                transport_plan[corner_positions[i][0]][corner_positions[i][1]] -= min_corner_value


if __name__ == '__main__':
    production = np.array([100, 300, 300])      # вектор поставок
    consumption = np.array([300, 200, 200])     # вектор заявок
    cost = np.array([
        [8, 4, 1],
        [8, 4, 3],
        [9, 7, 5]
    ])                                          # матрица стоимостей

    print("Исходные данные:")
    print("Вектор поставок: ", production)
    print("Вектор заявок:", consumption)
    print("Матрица стоимостей:\n", cost)

    optimal_plan = potentials_method(production, consumption, cost)
    print("\nОптимальный план =\n", optimal_plan)
