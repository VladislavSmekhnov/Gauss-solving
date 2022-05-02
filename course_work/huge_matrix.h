#pragma once

#include <vector>
#include <iostream>
#include <random>

template <typename T>
inline void create_huge_matrix(std::vector<std::vector<T>>& matrix, int lowerLimit, int size)
{
	//std::vector<std::vector<T>> huge_matrix;
	std::random_device random_device; // Источник энтропии
	std::mt19937 generator(random_device()); // Генератор случайных чисел
	std::uniform_int_distribution<> distribution(lowerLimit, lowerLimit * 2);

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			if (i == j)
				matrix[i].push_back(distribution(generator));
			else
				matrix[i].push_back(rand() % 2);
		}
	}
}

template <typename T>
inline void create_freeMatrixColumn(std::vector<T> &freeMatrixColumn, int size)
{
	for (int i = 0; i < size; i++)
		freeMatrixColumn.push_back((rand() % 100) / 10);
}
