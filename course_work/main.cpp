#include <iostream>
#include <vector>
#include <cmath>
#include <omp.h>
#include "gauss_serial.h"
#include "gauss_parallel.h"
#include "huge_matrix.h"

bool check_answer(std::vector<double> answer, std::vector<double> freeMatrixCloumn, std::vector<std::vector<double>> matrix);
void create_new_freeColumn(std::vector<double>& v, std::vector<double>& answer, std::vector<std::vector<double>> matrix);
void keep_4_digits_after_point(std::vector<double>& v);

int main()
{
    setlocale(LC_ALL, "utf-8");

    double timeStart_s, timeEnd_s, tick_s;
    double timeStart_p, timeEnd_p, tick_p;
    std::vector<double> answer_s;
    std::vector<double> answer_p;
    int m; // количество строк
    int n; // количество столбцов
    double value; // значение элемента массива
    std::vector<std::vector<double>> matrix_s, matrix_p, matrixForChecking; // матрица элементов без столбца свободных членов
    std::vector<double> freeMatrixColumn_s, freeMatrixColumn_p, freeMatrixColumnForChecking; // столбец свободных членов

    printf("Решение СЛАУ методом Гаусса\n");
    // Задаём размерность матрицы
    printf("Укажите размерность матрицы\n");
    printf("Количество строк: ");
    std::cin >> m;
    matrix_s.resize(m);
    matrix_p.resize(m);
    matrixForChecking.resize(m);
    printf("Количество столбцов: ");
    std::cin >> n;

    create_huge_matrix(matrix_s, 5, n);
    matrixForChecking = matrix_p = matrix_s;

    create_freeMatrixColumn(freeMatrixColumn_s, n);
    freeMatrixColumnForChecking = freeMatrixColumn_p = freeMatrixColumn_s;

    //// Заполняем матрицу
    //printf("Заполните матрицу (без столбца свободных членов):\n");
    //for (int i = 0; i < m; i++)
    //{
    //    for (int j = 0; j < n; j++)
    //    {
    //        std::cin >> value;
    //        matrix_s[i].push_back(value);
    //        matrix_p[i].push_back(value);
    //        matrixForChecking[i].push_back(value);
    //    }
    //}

    //printf("Теперь заполните столбец свободных членов:\n");
    //for (int i = 0; i < m; i++)
    //{
    //    std::cin >> value;
    //    freeMatrixColumn_s.push_back(value);
    //    freeMatrixColumn_p.push_back(value);
    //    freeMatrixColumnForChecking.push_back(value);
    //}

    /* Последовательное решение : */
    printf("\nРешаю...\n\n");
    printf("Последовательно...\n");
    timeStart_s = omp_get_wtime();
    answer_s = gauss_solving(matrix_s, freeMatrixColumn_s, m);
    timeEnd_s = omp_get_wtime();
    tick_s = omp_get_wtick();
    printf("Параллельно...\n");
    timeStart_p = omp_get_wtime();
    answer_p = parallel_gauss_solving(matrix_p, freeMatrixColumn_p, m);
    timeEnd_p = omp_get_wtime();
    tick_p = omp_get_wtick();

    //printf("Ответ (последовательное решение):\n");
    // Ответ:
    //for (int i = 0; i < answer_s.size(); i++)
    //    std::cout << answer_s[i] << " ";

    //printf("\n");
    ///* Конец последовательного решения */

    ///* Параллельное решение : */
    //printf("Ответ (параллельное решение):\n");
    // Ответ:
    //for (int i = 0; i < answer_p.size(); i++)
    //    std::cout << answer_p[i] << " ";

    //printf("\n");
    /* Конец параллельного решения */
    
    // Измеряем время, затраченное на решение СЛАУ обоими методами:
    printf("\nВремя, затраченное в последовательной области:\t%1f\n", timeEnd_s - timeStart_s);
    printf("Время, затраченное в параллельной области:\t%1f\n\n", timeEnd_p - timeStart_p);
    printf("Точность таймера для последовательной области:\t%1f\n", tick_s);
    printf("Точность таймера для параллельной области:\t%1f\n\n", tick_p);

    printf("\nПроверка решения:\n");
    /* Проверка последовательного решения */
    std::cout << std::boolalpha << "Последовательное решение верно:\t" <<
        check_answer(answer_s, freeMatrixColumnForChecking, matrixForChecking) << std::endl;

    /* Проверка параллельного решения */
    std::cout << std::boolalpha << "Параллельное решение верно:\t" <<
        check_answer(answer_p, freeMatrixColumnForChecking, matrixForChecking) << std::endl;

    system("pause");
    return 0;
}

bool check_answer(std::vector<double> answer, std::vector<double> freeMatrixCloumn, std::vector<std::vector<double>> matrix)
{
    std::vector<double> myFreeColumn(freeMatrixCloumn.size(), 0.0);

    create_new_freeColumn(myFreeColumn, answer, matrix);
    for (int i = 0; i < freeMatrixCloumn.size(); i++)
        if (myFreeColumn[i] != freeMatrixCloumn[i]) return false;

    return true;
}

void create_new_freeColumn(std::vector<double> &v, std::vector<double>& answer, std::vector<std::vector<double>> matrix)
{
    for (int i = 0; i < matrix.size(); i++)
    {
        for (int j = 0; j < matrix.size(); j++)
            v[i] += matrix[i][j] * answer[j];
    }
    keep_4_digits_after_point(v);
}

void keep_4_digits_after_point(std::vector<double>& v)
{
    for (int i = 0; i < v.size(); i++)
    {
        v[i] *= 10000;
        v[i] = nearbyint(v[i]);
        v[i] /= 10000;
    }
}
