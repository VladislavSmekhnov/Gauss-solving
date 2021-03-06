#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <omp.h>
#include "gauss_sequential.h"
#include "gauss_parallel.h"
#include "huge_matrix.h"

void Keep4DigitsAfterPoint(std::vector<double>& v);
void CreateNewFreeColumn(std::vector<double>& v, std::vector<double>& answer, std::vector<std::vector<double>> matrix);
bool CheckAnswer(std::vector<double> answer, std::vector<double> freeMatrixCloumn, std::vector<std::vector<double>> matrix);
template <typename T>
std::string CompareComputationTime(T sequential_t, T parallel_t);

int main()
{
    setlocale(LC_ALL, "utf-8");

    double timeStart_s, timeEnd_s, tick_s; // время работы последовательной области
    double timeStart_p, timeEnd_p, tick_p; // время работы параллельной области
    int matrix_dimension; // размерность матрицы
    unsigned int number_of_threads;
    double value; // значение элемента массива
    std::vector<double> answer_s; // ответ последовательного решения
    std::vector<double> answer_p; // ответ параллельного решения
    std::vector<std::vector<double>> matrix_s, matrix_p, matrix_for_checking; // матрица элементов без столбца свободных членов
    std::vector<double> freeMatrixColumn_s, freeMatrixColumn_p, free_matrix_column_for_checking; // столбец свободных членов

    printf("Решение СЛАУ методом Гаусса\n\n");
    printf("Задайте количество потоков: ");
    std::cin >> number_of_threads;
    omp_set_num_threads(number_of_threads);
    printf("\n");
    // Задаём размерность матрицы
    printf("Укажите размерность матрицы: ");
    std::cin >> matrix_dimension;
    while (matrix_dimension < 2)
    {
        std::cerr << "Минимальный размер матрицы: 2×2! Попробуйте ещё раз." << std::endl;
        printf("Укажите размерность матрицы: ");
        std::cin >> matrix_dimension;
    }
    matrix_s.resize(matrix_dimension);
    matrix_p.resize(matrix_dimension);
    matrix_for_checking.resize(matrix_dimension);

    printf("\n1 - ввести значения для матрицы коэффициентов \'A\' и вектора свободных членов \'b\' вручную.\n");
    printf("2 - сгенерировать все значения автоматически\n\n");
    do
    {
        printf("Введите номер выбранного вами варианта: ");
        std::cin >> value;
        if (value != 1 && value != 2) std::cerr << "Введено неправильное значение! Попробуйте ещё раз." << std::endl;
    } while (value != 1 && value != 2);
    if (value == 1)
    {
        printf("Заполните матрицу (без столбца свободных членов):\n");
        for (int i = 0; i < matrix_dimension; i++)
        {
            for (int j = 0; j < matrix_dimension; j++)
            {
                printf("\t");
                std::cin >> value;
                matrix_s[i].push_back(value);
                matrix_p[i].push_back(value);
                matrix_for_checking[i].push_back(value);
            }
        }

        printf("Теперь заполните столбец свободных членов:\n");
        for (int i = 0; i < matrix_dimension; i++)
        {
            printf("\t");
            std::cin >> value;
            freeMatrixColumn_s.push_back(value);
            freeMatrixColumn_p.push_back(value);
            free_matrix_column_for_checking.push_back(value);
        }
    }
    else if (value == 2)
    {
        printf("\nСоздаю матрицу размером: %d×%d", matrix_dimension, matrix_dimension);
        CreateHugeMatrix(matrix_s, 5, matrix_dimension);
        matrix_for_checking = matrix_p = matrix_s;
        CreateFreeMatrixColumn(freeMatrixColumn_s, matrix_dimension);
        free_matrix_column_for_checking = freeMatrixColumn_p = freeMatrixColumn_s;
    }
    else
        std::cerr << "Введено неправильное значение!" << std::endl;

    printf("\n\nРешаю...\n");

    printf("Последовательно...\n");
    /* Последовательное решение : */
    timeStart_s = omp_get_wtime();
    answer_s = SolveSequentially(matrix_s, freeMatrixColumn_s, matrix_dimension);
    timeEnd_s = omp_get_wtime();
    tick_s = omp_get_wtick();

    printf("Параллельно...\n");
    /* Параллельное решение : */
    timeStart_p = omp_get_wtime();
    answer_p = SolveInParallel(matrix_p, freeMatrixColumn_p, matrix_dimension);
    timeEnd_p = omp_get_wtime();
    tick_p = omp_get_wtick();
    
    // Измеряем время, затраченное на решение СЛАУ обоими способами:
    std::cout << "\nВремя, затраченное в последовательной области:\t" << std::round((timeEnd_s - timeStart_s) * 10000) / 10000 << std::endl;
    std::cout << "Время, затраченное в параллельной области:\t" << std::round((timeEnd_p - timeStart_p) * 10000) / 10000 << std::endl;
    std::cout << CompareComputationTime(timeEnd_s - timeStart_s, timeEnd_p - timeStart_p) << std::endl;
    printf("Точность таймера для последовательной области:\t%f\n", tick_s);
    printf("Точность таймера для параллельной области:\t%f\n\n", tick_p);

    printf("Проверка решения:\n");
    /* Проверка последовательного решения */
    std::cout << std::boolalpha << "Последовательное решение верно:\t" <<
        CheckAnswer(answer_s, free_matrix_column_for_checking, matrix_for_checking) << std::endl;
    /* Проверка параллельного решения */
    std::cout << std::boolalpha << "Параллельное решение верно:\t" <<
        CheckAnswer(answer_p, free_matrix_column_for_checking, matrix_for_checking) << std::endl;

    printf("Вывести ли ответы? (0 - нет, 1 - да): ");
    std::cin >> value;
    if (value == 1) PrintVectorS(answer_s);

    system("pause");
    return 0;
}

template <typename T>
std::string CompareComputationTime(T sequential_t, T parallel_t)
{
    char buffer[256] = {};
    sequential_t = std::round(sequential_t * 100) / 100;
    parallel_t = std::round(parallel_t * 100) / 100;
    
    if (sequential_t < parallel_t)
    {
        _gcvt_s(buffer, 255, std::round((parallel_t / sequential_t) * 100) / 100, 5);
        std::string result(buffer, sizeof(buffer));
        return "Последовательные вычисления быстрее параллельных в " + result + " раза";
    }
    else if (parallel_t < sequential_t)
    {
        _gcvt_s(buffer, 255, std::round((sequential_t / parallel_t) * 100) / 100, 5);
        std::string result(buffer, sizeof(buffer));
        return "Параллельные вычисления быстрее последовательных в " + result + " раза";
    }
        
    else
        return "Время последовательных и параллельных вычислений одинаково";
}

bool CheckAnswer(std::vector<double> answer, std::vector<double> freeMatrixCloumn, std::vector<std::vector<double>> matrix)
{
    std::vector<double> myFreeColumn(freeMatrixCloumn.size(), 0.0);

    CreateNewFreeColumn(myFreeColumn, answer, matrix);
    for (int i = 0; i < freeMatrixCloumn.size(); i++)
        if (myFreeColumn[i] != freeMatrixCloumn[i]) return false;

    return true;
}

void CreateNewFreeColumn(std::vector<double> &v, std::vector<double>& answer, std::vector<std::vector<double>> matrix)
{
    for (int i = 0; i < matrix.size(); i++)
    {
        for (int j = 0; j < matrix.size(); j++)
            v[i] += matrix[i][j] * answer[j];
    }
    Keep4DigitsAfterPoint(v);
}

void Keep4DigitsAfterPoint(std::vector<double>& v)
{
    for (int i = 0; i < v.size(); i++)
    {
        v[i] *= 10000;
        v[i] = nearbyint(v[i]);
        v[i] /= 10000;
    }
}
