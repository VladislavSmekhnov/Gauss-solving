#include <iostream>
//#include <omp.h>
#include <vector>
#include <cmath>
#include "gauss_serial.h"
#include "gauss_parallel.h"

int main()
{
    setlocale(LC_ALL, "utf-8");

    int m; // количество строк
    int n; // количество столбцов
    double value; // значение элемента массива
    std::vector<std::vector<double>> matrix_s, matrix_p; // матрица элементов без столбца свободных членов
    std::vector<double> freeMatrixColumn_s, freeMatrixColumn_p; // столбец свободных членов

    printf("Решение СЛАУ методом Гаусса\n");
    // Задаём размерность матрицы
    printf("Укажите размерность матрицы\n");
    printf("Количество строк: ");
    std::cin >> m;
    //freeMatrixColumn_s.resize(m);
    matrix_s.resize(m);
    //freeMatrixColumn_p.resize(m);
    matrix_p.resize(m);
    printf("Количество столбцов: ");
    std::cin >> n;
    // Заполняем матрицу
    printf("Заполните матрицу (без столбца свободных членов):\n");
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cin >> value;
            matrix_s[i].push_back(value);
            matrix_p[i].push_back(value);
        }
    }
    printf("Теперь заполните столбец свободных членов:\n");
    for (int i = 0; i < m; i++)
    {
        std::cin >> value;
        freeMatrixColumn_s.push_back(value);
        freeMatrixColumn_p.push_back(value);
    }

    /* Последовательное решение : */
    printf("Решаю...\n\n");
    printf("Ответ (последовательное решение):\n");
    std::vector<double> answer1 = gauss_solving(matrix_s, freeMatrixColumn_s, m);
    // Ответ:
    for (int i = 0; i < answer1.size(); i++)
        std::cout << answer1[i] << " ";

    printf("\n");
    /* Конец последовательного решения */

    /* Параллельное решение : */
    printf("Ответ (параллельное решение):\n");
    std::vector<double> answer2 = parallel_gauss_solving(matrix_p, freeMatrixColumn_p, m);
    for (int i = 0; i < answer2.size(); i++)
        std::cout << answer2[i] << " ";

    printf("\n");
    /* Конец параллельного решения */

    system("pause");
    return 0;
}
