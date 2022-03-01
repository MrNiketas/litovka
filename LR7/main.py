import numpy as np
import matplotlib.pyplot as plot
import math

# Коэффициенты
CoefficientsA = [0.0, -1, 1.5, 0.5, 0.5]
CheckCoefficients = 4

# Вариант 6
FirstX = 0
SecondX = 1 / math.e
FirstT = 0
SecondT = 1

# Погрешность
E = 0.01
H = 0.2


# Функция экстремали
def functionX(t):
    return (math.pow(math.e, -t) + (1 + math.e) * t * math.pow(math.e, -t) - 1) / 2


# Пробная функция w0
def functionW0(t):
    return t / math.e


# Производная пробной функции
def functionW0Derivative():
    return 1 / math.e


# Переключатель между пробными функциями
def functionWX(Coefficients, t):
    if CheckCoefficients == 2:
        return functionW2(Coefficients, t)
    elif CheckCoefficients == 3:
        return functionW3(Coefficients, t)
    else:
        return functionW4(Coefficients, t)


# Переключатель между производными пробных функций
def functionWXDerivative(Coefficients, t):
    if CheckCoefficients == 2:
        return functionW2Derivative(Coefficients, t)
    elif CheckCoefficients == 3:
        return functionW3Derivative(Coefficients, t)
    else:
        return functionW4Derivative(Coefficients, t)


# Пробная функция с двумя коэффициентами
def functionW2(Coefficients, t):
    return Coefficients[1] * (t ** 1) * (t / math.e - 1 / math.e) + Coefficients[2] * (t ** 2) * (
            t / math.e - 1 / math.e) + functionW0(t)


# Пробная функция с тремя коэффициентами
def functionW3(Coefficients, t):
    return CoefficientsA[1] * (t ** 1) * (t / math.e - 1 / math.e) + \
           Coefficients[2] * (t ** 2) * (t / math.e - 1 / math.e) + \
           Coefficients[3] * (t ** 3) * (t / math.e - 1 / math.e) + functionW0(t)


# Пробная функция с четыремя коэффициентами
def functionW4(Coefficients, t):
    return Coefficients[1] * (t ** 1) * (t / math.e - 1 / math.e) + \
           Coefficients[2] * (t ** 2) * (t / math.e - 1 / math.e) + \
           Coefficients[3] * (t ** 3) * (t / math.e - 1 / math.e) + \
           Coefficients[4] * (t ** 4) * (t / math.e - 1 / math.e) + functionW0(t)


# Производная пробной функции с двумя коэффициентами
def functionW2Derivative(Coefficients, t):
    return Coefficients[1] * ((1 + 1) * (t ** 1) / math.e - 1 * (t ** (1 - 1) / math.e)) + \
           Coefficients[2] * ((2 + 1) * (t ** 2) / math.e - 2 * (t ** (2 - 1) / math.e)) + \
           functionW0Derivative()


# Производная пробной функции с тремя коэффициентами
def functionW3Derivative(Coefficients, t):
    return Coefficients[1] * ((1 + 1) * (t ** 1) / math.e - 1 * (t ** (1 - 1) / math.e)) + \
           Coefficients[2] * ((2 + 1) * (t ** 2) / math.e - 2 * (t ** (2 - 1) / math.e)) + \
           Coefficients[3] * ((3 + 1) * (t ** 3) / math.e - 3 * (t ** (3 - 1) / math.e)) + functionW0Derivative()


# Производная пробной функции с четыремя коэффициентами
def functionW4Derivative(Coefficients, t):
    return Coefficients[1] * ((1 + 1) * (t ** 1) / math.e - 1 * (t ** (1 - 1) / math.e)) + \
           Coefficients[2] * ((2 + 1) * (t ** 2) / math.e - 2 * (t ** (2 - 1) / math.e)) + \
           Coefficients[3] * ((3 + 1) * (t ** 3) / math.e - 3 * (t ** (3 - 1) / math.e)) + \
           Coefficients[4] * ((4 + 1) * (t ** 4) / math.e - 4 * (t ** (4 - 1) / math.e)) + functionW0Derivative()


# Метод прямоугольников
def MethodRectangle(C):
    h = H
    a = FirstT
    n = (SecondT - FirstT) / h
    sumIntegral = 0.0
    for i in range(int(n)):
        x = a + i * h
        if (i == 0) | (i == n):
            sumIntegral = sumIntegral + functionWX(C, x)
        elif i % 2 == 0:
            sumIntegral = sumIntegral + (2 * functionWX(C, x))
        else:
            sumIntegral = sumIntegral + (4 * functionWX(C, x))
    return sumIntegral * h / 3


def MethodRitz():
    # Запоминаем начальные приближения
    coefficients = CoefficientsA

    # Считаем интеграл с начальными приближениями коэффициентов методом прямоугольников
    integralPrev = MethodRectangle(coefficients)

    # Указатель на коэффициент
    i = 1

    # Направления изменения коэффициентов
    move = -1

    while True:

        # Шаг изменения коэффициента
        h = 0.27

        # Вспомогательная переменная для отслеживания интеграла
        distancePrev = 0

        while True:

            coefficients[i] = coefficients[i] + move * h
            integralNext = MethodRectangle(coefficients)
            distanceNext = math.fabs(integralNext - integralPrev)

            # Условия выхода из 2го цикла
            if math.fabs(integralNext - integralPrev) < E:
                break

            # Если предыдущий интеграл больше найденного,то остаемся в этом направлении изменения коэффициентов
            elif integralPrev > integralNext:
                move = 1
            # Иначе меняем направление
            else:
                move = -1

            # Отслеживание интеграла
            if distancePrev == 0:
                distanceNext = distancePrev
            if distanceNext > distancePrev:
                move = move * -1

            # Изменения необходимые для того чтобы не зациклилось
            h = h * 0.92
            integralPrev = integralNext

        # Увеличиваем указатель на 1
        i = i + 1
        # Условия выхода из 1го цикла
        if i == CheckCoefficients + 1:
            break
    return [coefficients]


# Отрисовщик
def plotCreator(Coefficients):
    # Массивы для отрисовки экстремали из уравнения Эйлера
    tArray = np.linspace(FirstT, SecondT, 100)
    xArray = [functionX(i) for i in tArray]

    # Массив для отрисовки экстремали полученной методом Ритца
    xResultArray = [functionWX(Coefficients, i) for i in tArray]

    plot.plot(tArray, xArray, tArray, xResultArray, "-o")
    plot.axhline(0, color='black', linestyle='--')
    plot.axvline(0, color='black', linestyle='--')
    plot.grid(True)
    plot.xlabel('t')
    plot.ylabel('x')
    plot.legend(['Метод сопряженных градиентов', "Точка минимума"])
    plot.show()


if __name__ == '__main__':
    list_out = MethodRitz()  # Запуск алгоритма метода Ритца, list_out - выходные результаты метода
    coefficients = list_out[0]  # Так как мы ищем коэффициенты, то нас интересуют именно они

    # Для красоты вывода
    print("Количество коэффициентов:", CheckCoefficients)
    if CheckCoefficients >= 2:
        print("Коэффициент 1 :", coefficients[1])
        print("Коэффициент 2 :", coefficients[2])
    if CheckCoefficients >= 3:
        print("Коэффициент 3 :", coefficients[3])
    if CheckCoefficients >= 4:
        print("Коэффициент 4 :", coefficients[4])

    # Запуск отрисовщика
    plotCreator(coefficients)
