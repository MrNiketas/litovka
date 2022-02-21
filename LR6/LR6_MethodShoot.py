import numpy as np
import matplotlib.pyplot as plot
import math

# Z0 - угловой коэффициенты, DeltaT - перемещение по T
Z0 = 0.5
DeltaT = 0.05

# Вариант 6
FirstX = 0
SecondX = 1 / math.e
FirstT = 0
SecondT = 1

# Погрешность
E1 = 0.0019
E2 = 0.001


# Функция экстремали
def functionX(t):
    return (math.pow(math.e, -t) + (1 + math.e) * t * math.pow(math.e, -t) - 1) / 2


# Функция из ур-ия Эйлера
def functionFI(x, zPrev):
    return -zPrev - x - 0.5


# Метод пристрелки
def MethodShoot():
    globalT = 0
    deltaT = DeltaT
    deltaX = 0

    iMax = 10000
    xArray = np.zeros(iMax)
    tArray = np.zeros(iMax)

    z0 = Z0
    zPrev = 0
    zNext = 0

    xPrev = 0
    xNext = 0

    while True:
        globalT = FirstT
        delT = deltaT

        zPrev = z0
        zNext = z0
        xPrev = FirstX
        xNext = FirstX

        xArray[0] = FirstX
        tArray[0] = FirstT
        i = 1
        while (SecondT - globalT) > E1:
            xNext = xPrev + zPrev * delT
            xPrev = xNext
            xArray[i] = xNext

            if (i == 1):
                zNext = zPrev + functionFI(xNext, 0) * delT
                zPrev = zNext
            else:
                zNext = zPrev + functionFI(xNext, zPrev) * delT
                zPrev = zNext

            globalT = globalT + delT
            tArray[i] = globalT

            i = i + 1

        if math.fabs(xPrev - SecondX) < E1:
            print(z0)
            break
        else:
            if deltaX == 0:
                deltaX = math.fabs(xPrev - SecondX)
                z0 = z0 * 1.2
            elif deltaX > math.fabs(xPrev - SecondX):
                deltaX = math.fabs(xPrev - SecondX)
                z0 = z0 * 1.2
            else:
                if z0 < E2:
                    print(z0)
                    break
                z0 = z0 * 0.97
    print("Метод пристрелки:")
    xArray = xArray[0:i]
    tArray = tArray[0:i]
    return [xArray, tArray]


# Matplotlib
def plotCreator(xResult, tResult):
    tArray = np.linspace(FirstT, SecondT, 100)
    xArray = [functionX(i) for i in tArray]
    plot.plot(tArray, xArray, tResult, xResult, "-o")
    plot.axhline(0, color='black', linestyle='--')
    plot.axvline(0, color='black', linestyle='--')
    plot.grid(True)
    plot.xlabel('t')
    plot.ylabel('x')
    plot.legend(['Метод сопряженных градиентов', "Точка минимума"])
    plot.show()


if __name__ == '__main__':
    list_out = MethodShoot()
    xArray = list_out[0]
    tArray = list_out[1]
    plotCreator(xArray, tArray)
