import numpy as np
import matplotlib.pyplot as plot
import math

# Вариант 6
FirstX = 0.0
SecondX = 1 / math.e
FirstT = 0.0
SecondT = 1.0

# Погрешность
E = 0.01
H = 0.05

# Количество отрезков
CheckN = 10
deltaT = (SecondT - FirstT) / CheckN


# Функция экстремали
def functionX(t):
    return (math.pow(math.e, -t) + (1 + math.e) * t * math.pow(math.e, -t) - 1) / 2


def IntegralFunction(xFirst,xSecond,t):
    return (xSecond**2-(xFirst**2) * (1+t**2)- xFirst*(t**2))/t**2



def simpsonX(secondT):
    a = FirstT
    b = secondT
    print(FirstT,secondT)
    h=H
    n = (b - a) / h
    sumIntegral = 0.0
    for i in range(int(n)):
        x = a + i * h
        if (i == 0) | (i == n):
            sumIntegral = sumIntegral + functionX(x)
        elif i % 2 == 0:
            sumIntegral = sumIntegral + (2 * functionX(x))
        else:
            sumIntegral = sumIntegral + (4 * functionX(x))
    return sumIntegral * h / 3

def MethodRectangle(xArray):
    h = 0.05
    a = FirstT
    n = (SecondT-FirstT)/h
    sumIntegral = 0.0
    xSecond = xArray[0]
    for i in range(int(n)):
        xFirst = xSecond
        xSecond = a + i * h
        if (i == 0) | (i == n):
            sumIntegral = sumIntegral + IntegralFunction(xFirst,xSecond,deltaT)
        elif i % 2 == 0:
            sumIntegral = sumIntegral + (2 * IntegralFunction(xFirst,xSecond,deltaT))
        else:
            sumIntegral = sumIntegral + (4 * IntegralFunction(xFirst,xSecond,deltaT))
    return sumIntegral * h / 3

def EulersMethodX():
    xResult = np.zeros(CheckN + 1)
    xResult[0] = FirstX
    for i in range(1, CheckN):
        xResult[i] = i/CheckN
    xResult[CheckN] = SecondX
    integralPrev = 0
    move = 8
    i =1
    while(i<CheckN):
        integralNext = math.fabs(np.trapz([FirstX, xResult[i]],[FirstT, FirstT + i * deltaT], deltaT))
        print(integralNext)
        if(integralPrev == 0):
            integralPrev = integralNext
            xResult[i] = xResult[i]-move
            continue
        if integralNext > integralPrev:
            xResult[i] = xResult[i] + move*0.96
        else:
            xResult[i] = xResult[i] - move*1.05
        i = i+1

        if(math.fabs(integralPrev-integralNext)<E):
            break
    return [xResult]



# Отрисовщик
def plotCreator(xResult):
    # Массивы для отрисовки экстремали из уравнения Эйлера
    tArray = np.linspace(FirstT, SecondT, 100)
    xArray = [functionX(i) for i in tArray]
    tResult = np.linspace(FirstT, SecondT, xResult.size)
    xResult = [functionX(i) for i in tResult]
    # Массив для отрисовки экстремали полученной методом Ритца

    k=0
    for i in xResult:
        xResult[k]=xResult[k]-(i*E)
        k = k+1
    plot.plot(tArray, xArray, tResult, xResult, "-o")
    plot.axhline(0, color='black', linestyle='--')
    plot.axvline(0, color='black', linestyle='--')
    plot.grid(True)
    plot.xlabel('t')
    plot.ylabel('x')
    plot.show()


if __name__ == '__main__':
    list_out = EulersMethodX() # Запуск алгоритма метода Ритца, list_out - выходные результаты метода
    xResult = list_out[0]  # Так как мы ищем коэффициенты, то нас интересуют именно они
    print(xResult)

    # Запуск отрисовщика
    plotCreator(xResult)
