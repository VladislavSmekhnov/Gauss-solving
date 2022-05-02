#include <iostream>
#include <vector>
#include <cmath>
#include "gauss_serial.h"
#include "gauss_parallel.h"

bool check_answer(std::vector<double> answer, std::vector<double> freeMatrixCloumn, std::vector<std::vector<double>> matrix);
void create_new_freeColumn(std::vector<double>& v, std::vector<double>& answer, std::vector<std::vector<double>> matrix);
void keep_4_digits_after_point(std::vector<double>& v);

int main()
{
    setlocale(LC_ALL, "utf-8");

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
    // Заполняем матрицу
    printf("Заполните матрицу (без столбца свободных членов):\n");
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cin >> value;
            matrix_s[i].push_back(value);
            matrix_p[i].push_back(value);
            matrixForChecking[i].push_back(value);
        }
    }
    printf("Теперь заполните столбец свободных членов:\n");
    for (int i = 0; i < m; i++)
    {
        std::cin >> value;
        freeMatrixColumn_s.push_back(value);
        freeMatrixColumn_p.push_back(value);
        freeMatrixColumnForChecking.push_back(value);
    }

    /* Последовательное решение : */
    printf("\nРешаю...\n\n");
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
    printf("\nПроверка решения\n\n");
    /* Проверка последовательного решения */
    std::cout << std::boolalpha << "Последовательное решение верно:\t" <<
        check_answer(answer1, freeMatrixColumnForChecking, matrixForChecking) << std::endl;

    /* Проверка параллельного решения */
    std::cout << std::boolalpha << "Параллельное решение верно:\t" <<
        check_answer(answer2, freeMatrixColumnForChecking, matrixForChecking) << std::endl;

    system("pause");
    return 0;
}

bool check_answer(std::vector<double> answer, std::vector<double> freeMatrixCloumn, std::vector<std::vector<double>> matrix)
{
    std::vector<double> myFreeColumn(freeMatrixCloumn.size(), 0.0);

    create_new_freeColumn(myFreeColumn, answer, matrix);
    
    for (int i = 0; i < freeMatrixCloumn.size(); i++)
    {
        if (myFreeColumn[i] != freeMatrixCloumn[i])
            return false;
    }

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
        //printf("Умножили на 10000\n");
        v[i] = nearbyint(v[i]);
        //printf("Округляем умноженное число (избавляемся от знаков после точки\n");
        v[i] /= 10000;
        //printf("Возвращаем число в исходную форму\n");
    }
}
