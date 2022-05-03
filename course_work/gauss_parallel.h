#pragma once
/*Здесь реализован параллельный способ решения СЛАУ методом Гаусса*/
#include <vector>
#include <iostream>
#include <omp.h>
//=============================================================================================================
int ParallelNoSolution()
{
    std::cerr << "Нет решения!" << std::endl;
    exit(-1);
    return 0;
}
template <typename T>
void PrintVectorP(std::vector<T>& x);
template <typename T>
int ParallelFindMaxInCol(const std::vector<std::vector<T>>& matrix, int col, int matrix_dimension);
template <typename T>
int ParallelTriangulateMatrix(std::vector<std::vector<T>>& matrix, int matrix_dimension);
template <typename T>
std::vector<T> SolveInParallel(std::vector<std::vector<T>>& matrix, std::vector<T>& freeMatrixColumn, int matrix_dimension);
//=============================================================================================================
template<typename T>
inline int ParallelFindMaxInCol(const std::vector<std::vector<T>>& matrix, int col, int matrix_dimension)
{ // Здесь ищем максимальный элемент в столбце
    /*Суть алгоритма:
    Для «обнуления» элементов i-го столбца матрицы достаточно
    ко всем строкам с номерами j = i+1, ... n прибавить i-тую строку, домноженную на ( -a[j][i]/a[i][i] ).
    При выполнении такой операции может возникать деление на 0, если элемент на главной диагонали
    окажется равным нулю — в этом случае выполняют перестановку строк матрицы.
    Наиболее эффективным подходом к перестановке строк является перестановка i-той строки со строкой,
    имеющей максимальный по модулю элемент в i-том столбце.*/

    int max_pos = col; // номер строки, на которой находится максимальный элемент
    #pragma omp parallel shared(max_pos) num_threads(8)
    {
        T max = std::abs(matrix[col][col]); // модуль значения элемента
        #pragma omp for
        for (int i = col + 1; i < matrix_dimension; i++)
        {
            T element = std::abs(matrix[i][col]); // Берём новый элемент из столбца
            #pragma omp critical
            {
                if (element > max) // Если значение по модулю нового элемента больше уже имеющегося максимума, то
                {
                    max = element; // обновляем максимум
                    max_pos = i; // и запоминаем номер строки, на которой находится максимальный элемент
                }
            }
        }
    }

    return max_pos;
}

template<typename T>
inline int ParallelTriangulateMatrix(std::vector<std::vector<T>>& matrix, int matrix_dimension)
{ // Здесь триангулируем матрицу (приводим матрицу к треугольному виду)
    unsigned int swap_counter = 0; // количество перестановок строк

    if (matrix_dimension == 0) return swap_counter; // Если размерность матрицы = 0, то ничего не делаем
    // Полагаем, что все строки матрицы имеют одинаковую длину.
    const int num_cols = matrix[0].size(); // количество столбцов в матрице
    unsigned int i_max; // индекс (номер) максимального элемента в текущем столбце
    
    for (int i = 0; i < matrix_dimension - 1; i++) // 'matrixDimension - 1', чтобы не рассматривать
    {                                               // столбец свободных членов матрицы
        i_max = ParallelFindMaxInCol(matrix, i, matrix_dimension); // Нашли индекс максимального элемента в солбце

        if (i != i_max) // Если текущий элемент не является максимальным, то
        {
            std::swap(matrix[i], matrix[i_max]); // переставляем строки
            swap_counter++; // считаем количество таких операций, т. к. такое действие меняет знак определителя
        }
        
        for (int j = i + 1; j < matrix_dimension; j++)
        {
            T coeff = -matrix[j][i] / matrix[i][i]; // вычислили коэффициент для обнуления элемента i-й строки
            
            #pragma omp parallel for shared(num_cols) num_threads(8)
            for (int k = i; k < num_cols; k++)
                matrix[j][k] += matrix[i][k] * coeff; // обнуляем элемент и складываем строки
        }
    }

    return swap_counter;
}

template<typename T>
inline std::vector<T> SolveInParallel(std::vector<std::vector<T>>& matrix, std::vector<T>& freeMatrixColumn, int matrix_dimension)
{ // Здесь решаем СЛАУ методом Гаусса
    std::vector<T> solution(matrix_dimension); // результирующий вектор

    #pragma omp parallel for num_threads(8)
    for (int i = 0; i < matrix_dimension; i++)
        matrix[i].push_back(freeMatrixColumn[i]); // добавляем столбец свободных членов

    ParallelTriangulateMatrix(matrix, matrix_dimension); // приводим матрицу к треугольному виду

    /* При использовании метода Гаусса для решения
    систем линейных алгебраических уравнений следует избегать приближенных вычислений */

    for (int i = matrix_dimension - 1; i >= 0; i--)
    { // Обходим уравнения с конца матрицы
        if (std::abs(matrix[i][i]) < 0.0001)
            throw ParallelNoSolution();
        // Уравнения перебираются начиная с последнего
        for (int j = i + 1; j < matrix_dimension; j++) // для каждого из них вычисляется сумма matrix[i][j] * solution[j]
            matrix[i][matrix_dimension] -= matrix[i][j] * solution[j];

        solution[i] = matrix[i][matrix_dimension] / matrix[i][i]; // находим значение корня текущего уравнения
    }

    return solution;
}

template <typename T>
inline void PrintVectorP(std::vector<T>& x)
{
    unsigned int y, length;

    std::cout << "Количество решений СЛАУ = " << x.size() << std::endl;

    if (x.size() > 50)
    {
        std::cout << "Вывести полностью или частично? (1 - полностью, 2 - частично): ";
        std::cin >> y;
        if (y == 1)
        {
            for (int i = 0; i < x.size(); i++)
            {
                std::cout << x[i] << " ";
                if (i % 12 == 0 && i > 0) printf("\n");
            }
            printf("\n");
        }
        else if (y == 2)
        {
            printf("Количество элементов для вывода за раз: ");
            std::cin >> length;
            int step = length;
            for (int i = 0; i < x.size(); )
            {
                for (; i < length; i++)
                {
                    std::cout << x[i] << " ";
                }
                printf("\n");
                length = (length + step) % (x.size() + 1);
                if (length == x.size()) break;
                printf("0 - выход, 1 - продолжить вывод корней СЛАУ: ");
                std::cin >> y;
                if (y == 0) break;
            }
        }
    }
    else
    {
        for (int i = 0; i < x.size(); i++)
        {
            std::cout << x[i] << " ";
            if (i % 12 == 0 && 1 > 0) printf("\n");
        }
        printf("\n");
    }

}
