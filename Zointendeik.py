import numdifftools as nd
from scipy.optimize import linprog
import math

# g = lambda x: (x ** 4) + x + 1
# grad1 = nd.Gradient(g)([1])
# print("Gradient of x ^ 4 + x+1 at x = 1 is ", grad1)

tochka = [0,0]
array_ogranicheniy = [1,1,2,1,5,5,-1,0,0,0,-1,0]
def rosen(x):
    return 2 * (x[0] ** 2) + 2 * x[1] ** 2 - 2 * x[1] * x[0] - 4 * x[0] - 6 * x[1]
   # return  (1-x[0]) ** 2 + 100 * (x[1] - x[0] ** 2) ** 2

array_new_ogranicheniy = []
def active_lin_ogranich(tochka,array_ogranicheniy):
    counter = 0
    for i in range(int(len(array_ogranicheniy)/3)):
        if array_ogranicheniy[counter] * tochka[0] + array_ogranicheniy[counter + 1] * tochka[1] == array_ogranicheniy[counter + 2]:
            print("Активное ограничение номер", i + 1)
            array_new_ogranicheniy.append([array_ogranicheniy[counter],array_ogranicheniy[counter + 1],array_ogranicheniy[counter + 2]])
        elif array_ogranicheniy[counter] * tochka[0] + array_ogranicheniy[counter + 1] * tochka[1] >= array_ogranicheniy[counter + 2]:
            print("Точка находится вне допустимой области")
        counter += 3

    # active = []
    # if tochka[0] + tochka[1] <= 4:
    #     active.append(1)
    # elif 3 * tochka[1] + 2 * tochka[2] <= 0:
    #     active.append(2)



# active_lin_ogranich(tochka,array_ogranicheniy)
gradient = nd.Gradient(rosen)(tochka)
print("Gradient ", gradient)



# Составим новую задачу лп
print("Новіе ограничения", array_new_ogranicheniy)


lhs_ineq = []
for i in array_new_ogranicheniy:
    lhs_ineq.append([i[0],i[1]])
print(lhs_ineq)
# lhs_ineq = [[-1,0],  # левая сторона неравенства
#              [0, -1]  # левая сторона второго неравенства
#             ]
print("Левая сторона",lhs_ineq)


rhs_ineq = []
for i in array_new_ogranicheniy:
    rhs_ineq.append(0)
# rhs_ineq = [0,  # правая сторона неравенства
#            0]  # правая сторона неравенства
print("Правая сторона",rhs_ineq)


# В методе Зойнтендека всегда присутствуют ограничения на переменные от -1 до 1
n = 2
ogranichenie = (-1, 1)
bnd = [ogranichenie for i in range(n)]  # закидываем ограничениями по количеству мерности нашей функции


opt = linprog(c=gradient, A_ub=lhs_ineq, b_ub=rhs_ineq, bounds=bnd,
              method="revised simplex")
print(opt.fun)
epsil = 0.01
if opt.fun < epsil:
    print("Мы дошли до нужной точки")
print(opt)
a = list(opt.x)
print(a, "наш вектор направления")
# print("Этап 2 Добро пожаловать в расчет длины шага. Вас встречает Свен.")

def sven(S, tekyachya_tocka):
    xtochki = []
    fvalue = []
    # print(tekyachya_tocka)
    y0 = opt.fun
    # print(y0,"  значение функциии в начальной точке")
    xtochki.append(tekyachya_tocka)
    fvalue.append(y0)
    dl = 0.1 * (math.sqrt(tekyachya_tocka[0] ** 2 + tekyachya_tocka[1] ** 2)) / (math.sqrt(S[0] ** 2 + S[1] ** 2))
    if dl == 0:
        dl = 0.01
    print(dl)
    f1 = rosen([tekyachya_tocka[0] + dl * S[0], tekyachya_tocka[1] + dl * S[1]])
    f2 = rosen([tekyachya_tocka[0] - dl * S[0], tekyachya_tocka[1] - dl * S[1]])
    print(f1,f2)
    if f1 <= y0:
        s = 1
        # x.append([x0[0] + dl * S[0], x0[1] + dl * S[1]])
        # print("Додатне лямбда")
        xtochki.append([tekyachya_tocka[0] + dl * S[0],tekyachya_tocka[1] + dl * S[1]])
        fvalue.append(f1)
    elif f2 >= y0:
        s = -1
        # print("Отрицательное лямбда")
        xtochki.append([tekyachya_tocka[0] - dl * S[0], tekyachya_tocka[1] - dl * S[1]])
        fvalue.append(f2)

    dl = s * dl
    while fvalue[-2] >= fvalue[-1]:
        dl = dl * 2
        xtochki.append([xtochki[-1][0] + dl * S[0], xtochki[-1][1] + dl * S[1]])
        fvalue.append(rosen(xtochki[-1]))

    xtochki.append([(xtochki[-1][0] + xtochki[-2][0]) / 2, (xtochki[-1][1] + xtochki[-2][1]) / 2])
    fvalue.append(rosen(xtochki[-1]))

    n = fvalue.index(min(fvalue))
    if len(xtochki) >= n + 2:
        # print(xtochki)
        return [xtochki[n - 1], xtochki[n + 1]]
    else:
        # print(xtochki)
        return [xtochki[n - 1], xtochki[n]]


# sven(a,tekyachya_tocka=tochka)
# aa = sven(a, tochka)[0]
# bb = sven(a, tochka)[1]
# print(aa,bb)
epsilon = 0.1

def golden_func(x0,S):
    a = sven(S, x0)[0]
    b = sven(S, x0)[1]

    L = [abs(b[0] - a[0]), abs(b[1] - a[1])]
    l1 = [a[0] + 0.382 * L[0], a[1] + 0.382 * L[1]]
    l2 = [a[0] + 0.618 * L[0], a[1] + 0.618 * L[1]]
    f1 = rosen([x0[0] + l1[0] * S[0], x0[1] + l1[1] * S[1]])
    f2 = rosen([x0[0] + l2[0] * S[0], x0[1] + l2[1] * S[1]])

    while L[0] > epsilon and L[1] > epsilon:
        if f1 > f2:
            a = l1
            b = b
        elif f1 < f2:
            a = a
            b = l2

        L = [abs(b[0] - a[0]), abs(b[1] - a[1])]
        l1 = [a[0] + 0.382 * L[0], a[1] + 0.382 * L[1]]
        l2 = [a[0] + 0.618 * L[0], a[1] + 0.618 * L[1]]
        f1 = rosen([x0[0] + l1[0]*S[0], x0[1] + l1[1] * S[1]])
        f2 = rosen([x0[0] + l2[0]*S[0], x0[1] + l2[1] * S[1]])
    print(L)
    return L



def pauel_func(x0, S):
    l1 = sven(S, x0)[0]
    l3 = sven(S, x0)[1]
    l2 = [(l1[0] + l3[0]) / 2, (l1[1] + l3[1]) / 2]

    l = [l2]
    f1 = rosen(l1)
    f2 = rosen(l2)
    f3 = rosen(l3)
    y = [f2]

    if f1 <= f3:
        l.append(l1)
        y.append(f1)
        dl = [-abs(l2[0] - l1[0]), -abs(l2[1] - l1[1])]
    elif f1 >= f3:
        l.append(l3)
        y.append(f3)
        dl = [abs(l2[0] - l3[0]), abs(l2[1] - l3[1])]
    l1 = l[0]
    l2 = l[-1]

    while y[-2] > y[-1]:
        if abs(y[-2] - y[-1]) > epsilon:
            dl = [2*dl[0], 2*dl[1]]
            l.append([l2[0] + dl[0], l2[1] + dl[1]])
            y.append(rosen(l[-1]))
        else:
            break

    l.append([(l[-1][0] + l[-2][0]) / 2, (l[-1][1] + l[-2][1])/2])
    y.append(rosen(l[-1]))
    n = y.index(min(y))

    return [abs(l[n-1][0] - l[n][0]), abs(l[n-1][1] - l[n][1])]


shag = golden_func([0,0],[1,1]) # интервал для шага меньший епсилон
shag2 = pauel_func([0,0],[1,1]) # интервал для шага меньший епсилон
print(shag)
print(shag2)
new_tochka = [tochka[0] + min(shag2) * a[0] + tochka[1] + min(shag2)] # добавляем шаг по направлению к точке
# шаг берем меньший из возможного интервала.
print(new_tochka)








