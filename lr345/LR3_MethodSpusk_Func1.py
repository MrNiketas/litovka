import numpy as np
import matplotlib.pyplot as plt
import math as mth

A = -2
B = 4
C = 2
D = 5
ALPHA = 35
E = 0.00001


def Function(X):  # ФУНКЦИЯ
    return (((X[0] - A) * mth.cos(ALPHA) + (X[1] - B) * mth.sin(ALPHA)) ** 2) / (C ** 2) + \
           (((X[1] - B) * mth.cos(ALPHA) - (X[0] - A) * mth.sin(ALPHA)) ** 2) / (D ** 2)


def GradientFunction(x):  # ГРАДИЕНТ ФУНКЦИИ
    return np.array([-C * (-(x[0] - A) * mth.sin(ALPHA) + (x[1] - B) * mth.cos(ALPHA)) * mth.sin(ALPHA) / D ** 2 + (
            (x[0] - A) * mth.cos(ALPHA) + (x[1] - B) * mth.sin(ALPHA)) * mth.cos(ALPHA) / C,
                     C * (-(x[0] - A) * mth.sin(ALPHA) + (x[1] - B) * mth.cos(ALPHA)) * mth.cos(ALPHA) / D ** 2 + (
                             (x[0] - A) * mth.cos(ALPHA) + (x[1] - B) * mth.sin(ALPHA)) * mth.sin(ALPHA) / C])


def LineSearch(x, beginF, H):
    flag = 0
    a = 0
    b = 1
    fX = Function(x)  # Значение функции в точке x
    gF = GradientFunction(x)  # Значение градиента функции в точке x
    rgF = -GradientFunction(x)  # Значение анти-градиента функции в точке x
    dfX = np.dot(gF, rgF)
    jump = b * 0.0001
    while flag == 0:
        newfX = Function(x + jump * rgF)
        if (newfX - fX) <= (H * jump * dfX):
            if (newfX - fX) >= ((1 - H) * jump * dfX):
                flag = 1
            else:
                a = jump
                if b < beginF:
                    jump = (a + b) / 2
                else:
                    jump = 2 * jump
    return jump


def MethodSpusk(xFirst):  # МЕТОД НАИСКОРЕЙШЕГО СКУСКА
    x0 = xFirst  # Начальная точка с координатами (X1;X2)
    iterationsMax = 10000
    X = np.zeros((2, iterationsMax))
    X[:, 0] = x0
    iteration = 1
    gradientFunction = GradientFunction(x0)
    xNext = x0
    xPrev = x0
    delta = sum(gradientFunction ** 2)
    while iteration < iterationsMax and delta > E:
        S = -GradientFunction(xNext)
        h = LineSearch(xNext, 1, 0.1)
        xNext = xPrev + h*S
        X[:, iteration] = xNext
        iteration = iteration + 1
        while Function(xNext) < Function(xPrev):
            xPrev = xNext
            xNext = xPrev + h * S

            if Function(xNext) > Function(xPrev):
                xNext = xPrev
                break
            X[:, iteration] = xNext
            iteration = iteration + 1
        gradientNext = GradientFunction(xNext)
        delta = sum(gradientNext ** 2)
        X[:, iteration] = xNext
        iteration = iteration + 1

    print("Результаты: ")
    print("///Точка минимума :", xNext)
    print("///Количество итераций :", iteration)
    X = X[:, 0: iteration]
    return [X]


def Paint(fx, fy):
    xArrayPlot = np.arange(-10, 10, 0.05)
    yArrayPlot = np.arange(-10, 10, 0.05)
    [X, Y] = np.meshgrid(xArrayPlot, yArrayPlot)
    plt.contour(X, Y, Function([X, Y]), 250)
    plt.plot(fx, fy, 'o-',A, B, 'o-')
    plt.axhline(0, color='black', linestyle='--')
    plt.axvline(0, color='black', linestyle='--')
    plt.grid(True)
    plt.xlabel('x1')
    plt.ylabel('x2')
    plt.legend(['Метод наискорейшего спуска', "Точка минимума"])
    plt.show()


if __name__ == "__main__":
    print("Введите X-координату: ")
    x = float(input())
    print("Введите Y-координату: ")
    y = float(input())
    result = MethodSpusk([x, y])
    X = result[0]
    Paint(X[0, :], X[1, :])