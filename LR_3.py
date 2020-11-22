# Для рассчетов и построения графика необходимы дополнительные библиотеки, импортируем их:

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Исходные данные (вариант 17):

# Подложка - Si <100>
#Окисление в H20
x1 = 0.1*2.4 # мкм - толщина 1го слоя
C1 = 10E+15 # см^-3
C2 = 6E+18 # см^-3
T = 1150 + 273 # К
P = 1.5 # атм
t2 = 30 # мин
A0_par = 7
EA_par = 0.78
A0_lin = 2.96E+6
EA_lin = 2.05
k = 8.617E-5 # Константа Больцмана, эВ

# Точноть расчетов:

n = 1E+3
dt2 = t2/n # мин - шаг по времени для 2-го слоя

# Вводим необходимые для расчетов функции:

def Eg(T):
    return 1.17-4.73E-4 * T**2 / (T + 636)

def ni(T):
    Nc = 6.2E+15 * T**1.5
    Nv = 3.5E+15 * T**1.5
    return (Nc * Nv) ** 0.5 * np.exp(-Eg(T) / (2 * k * T))

def K0(A0, EA, T):
    return A0 * np.exp(-EA/k/T)

# Расчет линейного коэффициента:

def B_A(K0, C, T):
    K = K0 / 1.68 # Зависимость коэффициента от ориентации подложки
    K = K * 0.85 * P # Зависимость коэффициента от давления
    YL = 2620 * np.exp(-1.1 / k / T)
    Ei = Eg(T) / 2 - k * T / 4
    C_minus = np.exp((Ei + 0.57 - Eg(T))/(k*T))
    C_2minus = np.exp((2 * Ei + 1.25 - 3 * Eg(T))/(k * T))
    C_plus = np.exp((0.35 - Ei)/(k * T))
    Vn = (1 + C_plus*(ni(T)/C) + C_minus*(C/ni(T)) + C_2minus*(C/ni(T))**2)/(1 + C_plus + C_minus + C_2minus)
    return K*(1 + YL*(Vn - 1))

# Расчет параболического коэффициента:

def B(K0, T, C):
    Yp = 9.63E-16*np.exp(2.83/k/T)
    return K0*(1 + Yp*C**0.22)

# Расчет толщины окисла:

def x(dt, xi, B, B_A):
    A = B/B_A
    tau = (xi**2 + A*xi)/B
    return A/2*((1 + (dt + tau)/(A**2/4/B))**0.5 - 1)

# Расчитываем окисление первого слоя:

B1 = B(K0(A0_par, EA_par, T), T, C1)
B_A1 = B_A(K0(A0_lin, EA_lin, T), T, C1)
A1 = B1/B_A1

# Для построения графика определим время, за которое протекает окисление 1-го слоя:

t1 = (A1**2)/4/B1*((1 + 2*x1/A1)**2 - 1)

dt1 = t1/n # мин - шаг по времени для 1-го слоя

# Создаем пустые массивы, в которые будут последовательно добавляться данные для постоения графиков:

t1_list = []
t2_list = []
x1_list = []
x2_list = []

# Задаем начальную толщину окисла:

xi = 0

# Расчет точек для графика:

for t in np.arange(0, t1, dt1): 
    t1_list.append(t) 
    xi = x(dt1, xi, B1, B_A1)
    x1_list.append(xi)

# Расчитываем окисление второго слоя:

B2 = B(K0(A0_par, EA_par, T), T, C2)
B_A2 = B_A(K0(A0_lin, EA_lin, T), T, C2)

for t in np.arange(0, t2, dt2):
    t2_list.append(t1 + t)
    xi = (x(dt2, xi, B2, B_A2))
    x2_list.append(xi)
    
# Объявим форму для рисования и оси:
fig, axes = plt.subplots()
        
# Добавим на график зависимости для окисления обоих слоев:
axes.plot(t1_list, x1_list, color='blue', label='1 слой')
axes.plot(t2_list, x2_list, color='orange', label='2 слой')

# Ограничим оси в нужном диапозоне:
axes.set_xlim(0, t1 + t2)
axes.set_ylim(0, x2_list[-1])

# Добавим название для графика и осей:
plt.title('Зависимость толщины окисла от времени',)
plt.xlabel('t, мин')
plt.ylabel('х, мкм')
plt.legend(loc=5)

# Добавление дополнительной ссетки:
axes.grid(which='major', color = '#666666')
axes.minorticks_on()
axes.grid(which='minor', color = 'gray', linestyle = ':')

# Выведем график на экран:
plt.show()