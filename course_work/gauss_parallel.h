#pragma once
/*Здесь реализован параллельный способ решения СЛАУ методом Гаусса*/
#include <vector>
#include <iostream>
#include <omp.h>
//=============================================================================================================
int parallel_no_solution()
{
    std::cerr << "Нет решения!" << std::endl;
    exit(-1);
    return 0;
}
template <typename T>
int parallel_find_max_inCol(const std::vector<std::vector<T>>& matrix, int col, int rowsAmount);
template <typename T>
int parallel_triangulate_matrix(std::vector<std::vector<T>>& matrix, int rowsAmount);
template <typename T>
std::vector<T> parallel_gauss_solving(std::vector<std::vector<T>>& matrix, std::vector<T>& freeMatrixColumn, int rowsAmount);
//=============================================================================================================
template<typename T>
inline int parallel_find_max_inCol(const std::vector<std::vector<T>>& matrix, int col, int rowsAmount)
{ // Здесь ищем максимальный элемент в столбце
    /*Суть алгоритма:
    Для «обнуления» элементов i-го столбца матрицы достаточно
    ко всем строкам с номерами j = i+1, ... n прибавить i-тую строку, домноженную на ( -a[j][i]/a[i][i] ).
    При выполнении такой операции может возникать деление на 0, если элемент на главной диагонали
    окажется равным нулю — в этом случае выполняют перестановку строк матрицы.
    Наиболее эффективным подходом к перестановке строк является перестановка i-той строки со строкой,
    имеющей максимальный по модулю элемент в i-том столбце.*/
    
    T max = std::abs(matrix[col][col]); // модуль значения элемента
    int maxPos = col; // номер строки, на которой находится максимальный элемент
    #pragma omp parallel for shared(max, maxPos) num_threads(10)
    for (int i = col + 1; i < rowsAmount; i++)
    {
        T element = std::abs(matrix[i][col]); // Берём новый элемент из столбца
        if (element > max) // Если значение по модулю нового элемента больше уже имеющегося максимума, то
        {
            max = element; // обновляем максимум
            maxPos = i; // и запоминаем номер строки, на которой находится максимальный элемент
        }
    }

    return maxPos;
}

template<typename T>
inline int parallel_triangulate_matrix(std::vector<std::vector<T>>& matrix, int rowsAmount)
{ // Здесь триангулируем матрицу (приводим матрицу к треугольному виду)
    unsigned int swapCounter = 0; // количество перестановок строк

    if (rowsAmount == 0) return swapCounter; // Если размерность матрицы = 0, то ничего не делаем
    // Полагаем, что все строки матрицы имеют одинаковую длину.
    const int num_cols = matrix[0].size(); // количество столбцов в матрице
    unsigned int imax; // индекс (номер) максимального элемента в текущем столбце
    
    for (int i = 0; i < rowsAmount - 1; i++) // 'matrixDimension - 1', чтобы не рассматривать
    {                                               // столбец свободных членов матрицы
        imax = parallel_find_max_inCol(matrix, i, rowsAmount); // Нашли индекс максимального элемента в солбце

        if (i != imax) // Если текущий элемент не является максимальным, то
        {
            swap(matrix[i], matrix[imax]); // переставляем строки
            swapCounter++; // считаем количество таких операций, т. к. такое действие меняет знак определителя
        }
        #pragma omp parallel for shared(num_cols, imax, swapCounter) collapse(2) num_threads(10)
        for (int j = i + 1; j < rowsAmount; j++)
        {
            T coeff = -matrix[j][i] / matrix[i][i]; // вычислили коэффициент для обнуления элемента i-й строки
            for (int k = i; k < num_cols; k++)
                matrix[j][k] += matrix[i][k] * coeff; // обнуляем элемент и складываем строки
        }
    }

    return swapCounter;
}

template<typename T>
inline std::vector<T> parallel_gauss_solving(std::vector<std::vector<T>>& matrix, std::vector<T>& freeMatrixColumn, int rowsAmount)
{ // Здесь решаем СЛАУ методом Гаусса
    std::vector<T> solution(rowsAmount); // результирующий вектор
    #pragma omp parallel for num_threads(10)
    for (int i = 0; i < rowsAmount; i++)
        matrix[i].push_back(freeMatrixColumn[i]); // добавляем столбец свободных членов

    parallel_triangulate_matrix(matrix, rowsAmount); // приводим матрицу к треугольному виду

    /* При использовании метода Гаусса для решения
    систем линейных алгебраических уравнений следует избегать приближенных вычислений */
    #pragma omp parallel for collapse(2) num_threads(10)
    for (int i = rowsAmount - 1; i >= 0; i--)
    { // Обходим уравнения с конца матрицы
        if (std::abs(matrix[i][i]) < 0.0001)
            throw parallel_no_solution();
        // Уравнения перебираются начиная с последнего
        for (int j = i + 1; j < rowsAmount; j++) // для каждого из них вычисляется сумма matrix[i][j] * solution[j]
        {
            #pragma omp critical
            matrix[i][rowsAmount] -= matrix[i][j] * solution[j];
        }
        //#pragma omp critical
        solution[i] = matrix[i][rowsAmount] / matrix[i][i]; // находим значение корня текущего уравнения
    }

    return solution;
}
