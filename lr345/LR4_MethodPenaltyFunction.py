import numpy as np
import matplotlib.pyplot as plt
import math as mth

A = 2
B = 1
C = 1
D = 4
ALPHA = 20
E = 0.001
KEY = 0.4


def Function(X):  # НАЧАЛЬНАЯ ФУНКЦИЯ
    return (((X[0] - A) * mth.cos(ALPHA) + (X[1] - B) * mth.sin(ALPHA)) ** 2) / (C ** 2) + \
           (((X[1] - B) * mth.cos(ALPHA) - (X[0] - A) * mth.sin(ALPHA)) ** 2) / (D ** 2)


def FunctionPenalty(X, K):  # ФУНКЦИЯ ШТРАФА
    return Function(X) + FunctionUnequalOUTSIDE(X, K)


def FunctionUnequalOUTSIDE(X, K):  # ШТРАФНАЯ ФУНКЦИЯ
    return K*(max(FirstUnequal(X), 0)+min(SecondUnequal(X), 0))


def FirstUnequal(X): # ОГРАНИЧЕНИЕ 1
    return mth.e**(X[0]/X[1])-X[0]


def SecondUnequal(X): # ОГРАНИЧЕНИЕ 2
    return X[1]


def CheckPoint(X): # ОГРАНИЧЕНИЯ В ВИДЕ НЕРАВЕНСТВА
    if(FirstUnequal(X) <= 0) and (SecondUnequal(X) >= 0):
        return True
    else:
        return False


def MethodSimplex(x0, K):
    X = np.zeros([2, 10000])
    x1 = x0
    x2 = np.array([x0[0] + KEY, x0[1]])
    x3 = np.array([x0[0], x0[1] + KEY])
    f1 = FunctionPenalty(x1, K)
    f2 = FunctionPenalty(x2, K)
    f3 = FunctionPenalty(x3, K)
    last = x0
    i = 0
    while True:
        alpha = 2
        fMax = max(f1, f2, f3)
        if fMax == f1:
            dx = (x2[0] + x3[0]) / 2 - x1[0]
            dy = (x2[1] + x3[1]) / 2 - x1[1]
            while True:
                x = np.array([x1[0]+alpha*dx, x1[1]+alpha*dy])
                fdX = FunctionPenalty(x, K)
                if fdX < f1:
                    x1 = x
                    f1 = fdX
                    last = x1

                    break
                else:
                    alpha = alpha*0.92
                if alpha < E:
                    break
        elif fMax == f2:
            dx = (x1[0] + x3[0]) / 2 - x2[0]
            dy = (x1[1] + x3[1]) / 2 - x2[1]
            while True:
                x = np.array([x2[0]+alpha*dx, x2[1]+alpha*dy])
                fdX = FunctionPenalty(x, K)
                if fdX < f2:
                    x2 = x
                    f2 = fdX
                    last = x2
                    break
                else:
                    alpha = alpha*0.92
                if alpha < E:
                    break
        else:
            dx = (x2[0] + x1[0]) / 2 - x3[0]
            dy = (x2[1] + x1[1]) / 2 - x3[1]
            while True:
                x = np.array([x3[0]+alpha*dx, x3[1]+alpha*dy])
                fdX = FunctionPenalty(x, K)
                if fdX < f3:
                    x3 = x
                    f3 = fdX
                    last = x3
                    break
                else:
                    alpha = alpha*0.92
                if alpha < E:
                    break
        i=i+1
        if alpha < E:
            break
    print("Количество итераций симплекса:", i)
    return last


def MethodPenaltyFunction(X):
    XArray = np.zeros([2,100])
    XArray[:, 0] = X
    K = 1
    xPrev = X
    xNext = X
    iteration = 1
    print("Значения первого ограничения (f<=0) в начальной точке:", FirstUnequal(xPrev))
    print("Значения второго ограничения (f>=0) в начальной точке:", SecondUnequal(xPrev))
    if not (CheckPoint(X)):
        print("Начальная точка не лежит в допустимой области")
    print()
    while True:
        xNext = MethodSimplex(xPrev, K)
        if not (CheckPoint(xNext)):
            print("Мы ушли из допустимой области")
            K = K * 10
            continue
        print(xNext)
        print()
        XArray[:, iteration] = xNext
        if mth.fabs(distance(xPrev, xNext)) < E:
            break
        xPrev = xNext
        iteration = iteration + 1
        K = K*10
    XArray = XArray[:, 0: iteration]
    print()
    print("Количество итераций для Метода штрафных функций: ", iteration)
    print("Точка минимума:", xNext)
    print("Значения первого ограничения (f<=0) в конечной точке:", FirstUnequal(xNext))
    print("Значения второго ограничения (f>=0) в конечной точке:", SecondUnequal(xNext))
    return [XArray]


def Paint(fx, fy):
    xArrayPlot = np.arange(-10, 10, 0.05)
    yArrayPlot = np.arange(-10, 10, 0.05)
    [X, Y] = np.meshgrid(xArrayPlot, yArrayPlot)
    plt.contour(X, Y, Function([X, Y]), 200)
    plt.contour(X, Y, FirstUnequal([X, Y]), 1)
    plt.contour(X, Y, SecondUnequal([X, Y]), 1)
    plt.plot(fx, fy, 'o-',A, B, 'o-')
    plt.axhline(0, color='black', linestyle='--')
    plt.axvline(0, color='black', linestyle='--')
    plt.grid(True)
    plt.xlabel('x1')
    plt.ylabel('x2')
    plt.legend(['Метод штрафных функций', "Точка минимума"])
    plt.show()


def distance(X0, X1):
    return np.sqrt((X1[0] - X0[0]) ** 2 + (X1[1] - X0[1]) ** 2)


def FunctionUnequalINSIDE(X, K):  # ФУНКЦИЯ ШТРАФА ВНУТРЕННЯЯ ФУНКЦИЯ
    if(FirstUnequal(X) == 0) or (SecondUnequal(X) == 0):
        return 0
    return K * ((FirstUnequal(X)+SecondUnequal(X))/(FirstUnequal(X)*SecondUnequal(X)))


if __name__ == "__main__":
    print("Введите X-координату: ")
    x = float(input())
    print("Введите Y-координату: ")
    y = float(input())
    result = MethodPenaltyFunction([x, y])
    X = result[0]
    Paint(X[0, :], X[1, :])