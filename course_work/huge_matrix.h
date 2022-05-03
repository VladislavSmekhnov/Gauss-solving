#pragma once

#include <vector>
#include <iostream>
#include <random>

template <typename T>
inline void CreateHugeMatrix(std::vector<std::vector<T>>& matrix, int lower_limit, int size)
{
	std::random_device random_device; // Источник энтропии
	std::mt19937 generator(random_device()); // Генератор случайных чисел
	std::uniform_int_distribution<> distribution(lower_limit, lower_limit * 2);

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (i == j) matrix[i].push_back(distribution(generator));
			else matrix[i].push_back(rand() % 2);
		}
	}
}

template <typename T>
inline void CreateFreeMatrixColumn(std::vector<T> &free_matrix_column, int size)
{
	for (int i = 0; i < size; i++)
		free_matrix_column.push_back((rand() % 100) / 10);
}
