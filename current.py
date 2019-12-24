#! /usr/bin/env python
# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import tkinter
from math import acos
import linalg
import time

#Это вообще не надо будет
#По двум точкам на плоскости "рисует" дугу
def misha_arc(a, b):
    c = a  / (a[0] ** 2 + a[1] ** 2)
    len_b = (a[0] - c[0])** 2 + (a[1] - c[1])** 2 
    len_c = (a[0] - b[0])** 2 + (a[1] - b[1])** 2 
    len_a = (b[0] - c[0])** 2 + (b[1] - c[1])** 2
    if(np.fabs(-len_a ** 2 - len_b ** 2 - len_c ** 2 + 2 * (len_a * len_b + len_b * len_c + len_c * len_a)) < 1e-20):
        return np.array([0, 0, 0, 0, 0, 2])
    x_O = (len_a * (len_b + len_c - len_a) * a[0] + len_b * (len_a + len_c - len_b) * b[0] + len_c * (len_b + len_a - len_c) * c[0]) / (-len_a ** 2 - len_b ** 2 - len_c ** 2 + 2 * (len_a * len_b + len_b * len_c + len_c * len_a))
    y_O = (len_a * (len_b + len_c - len_a) * a[1] + len_b * (len_a + len_c - len_b) * b[1] + len_c * (len_b + len_a - len_c) * c[1]) / (-len_a ** 2 - len_b ** 2 - len_c ** 2 + 2 * (len_a * len_b + len_b * len_c + len_c * len_a))
    if(-1 + x_O ** 2 + y_O ** 2 < 0 or np.fabs(-1 + x_O ** 2 + y_O ** 2) < 1e-20):
        return np.array([0, 0, 0, 0, 0, 2])
    R = (-1 + x_O ** 2 + y_O ** 2) ** 0.5
    flag = 1
    if(np.fabs((a[0] - x_O) / R) <= 1):
        angle_1 = acos((a[0] - x_O) / R)
    elif(a[0] - x_O > 0):
        angle_1 = 0
    else:
        angle_1 = np.pi 
    if(np.fabs((b[0] - x_O) / R) <= 1):
        angle_2 = acos((b[0] - x_O) / R)
    elif(b[0] - x_O > 0):
        angle_2 = 0
    else:
        angle_2 = np.pi
    if(a[1] < y_O):
        angle_1 = 2 * np.pi - angle_1
    if(b[1] < y_O):
        angle_2 = 2 * np.pi - angle_2
    if(angle_1 > angle_2):
        end = angle_1
        begin = angle_2
        flag = -1
    else:
        begin = angle_1
        end = angle_2
    if(end - begin > np.pi):
        d = begin
        begin = end - 2 * np.pi
        end = d
        flag = -flag
#Возвращает центр, аргумент начала дуги и аргумент конца дуги, а также радиус и флаг - в каком порядке ее надо рисовать
    return np.array([x_O, y_O, begin, end, R, flag])

# def Draw(p, q, T, word):
#     ar = []
#     xf = []
#     yf = []
#     for i in range(p + 1):
#         angle = 2 * i * np.pi / p
#         x = np.cos(angle)
#         y = np.sin(angle)
#         r = np.cos(np.pi / q + np.pi / p) * (np.cos(np.pi / q) ** 2 - np.sin(np.pi / p) ** 2) ** (-0.5)
#         x *= r
#         y *= r
#         pt = change_to_circle(mult_mat_vec(T, change_to_r3([x, y])))
#         if i != 0:
#             tmp = misha_arc(pt, old)
#             if(tmp[5] == -1):
#                 theta = np.arange(tmp[2], tmp[3], (tmp[3] - tmp[2]) / 100)
#                 X = tmp[0] + tmp[4] * np.cos(theta)
#                 Y = tmp[1] + tmp[4] * np.sin(theta)
#                 xf = np.hstack((xf, X))
#                 yf = np.hstack((yf, Y))
#             elif(tmp[5] == 1):                
#                 theta = np.arange(tmp[3], tmp[2], (tmp[2] - tmp[3]) / 100)
#                 X = tmp[0] + tmp[4] * np.cos(theta)
#                 Y = tmp[1] + tmp[4] * np.sin(theta)
#                 xf = np.hstack((xf, X))
#                 yf = np.hstack((yf, Y))
#             else:
#                 theta = np.arange(0, 1, 0.01)
#                 X = old[0] + (-old[0] + pt[0]) * theta
#                 Y = old[1] + (-old[1] + pt[1]) * theta
#                 xf = np.hstack((xf, X))
#                 yf = np.hstack((yf, Y))               
#         old = pt
#     plt.fill(xf, yf, color = decode(word, p, q))

  
#перемножает матрицы. типа очень быстро пермножает
def mult_mat3(a, b):
    return np.array(a) @ np.array(b)

#матрица на вектор
def mult_mat_vec(a, b):
    c = [0, 0, 0]
    for i in range(3):
        for j in range(3):
            c[i] += a[i][j] * b[j]
    return c

#замена координат
def change_to_circle(P):
    a = 1 / (1 + P[2])
    return [a * P[0],  a * P[1]]

#замена координат
def change_to_r3(P):
    a = 1 / (1 - P[0] ** 2 - P[1] ** 2)
    return [2 * a * P[0], 2 * a * P[1], (1 + P[0] ** 2 + P[1] ** 2) * a]

#Делает матрицу гиперболического поворота на тетта в одной какой-то(вроде ХЗет) плоскости, и альфа в плоскости ХУ
def make_Q(alpha, theta):
    A = [[np.cos(alpha), -np.sin(alpha), 0],
        [np.sin(alpha), np.cos(alpha), 0],
        [0, 0, 1]]
    B = [[np.cosh(theta), 0, np.sinh(theta)],
        [0, 1, 0],
        [np.sinh(theta), 0, np.cosh(theta)]]
    return np.array(A) @ np.array(B)


def E():
    return np.identity(3, dtype=float)

#Просто поворот на 2 пи / р
def RotateP(p, q):
    return np.array([[np.cos(2 * np.pi / p), -np.sin(2 * np.pi / p), 0],
                    [np.sin(2 * np.pi / p), np.cos(2 * np.pi / p), 0],
                    [0, 0, 1]
                    ])
#То же самое на удвоенный угол
def Rotate2P(p, q):
    return np.array([[np.cos(4 * np.pi / p), -np.sin(4 * np.pi / p), 0],
                    [np.sin(4 * np.pi / p), np.cos(4 * np.pi / p), 0],
                    [0, 0, 1]])
#Утроенный
def Rotate3P(p, q):
    return np.array([[np.cos(6 * np.pi / p), -np.sin(6 * np.pi / p), 0],
                    [np.sin(6 * np.pi / p), np.cos(6 * np.pi / p), 0],
                    [0, 0, 1]])
#Поворот относительно вершины многоугольника. Образующая группы симметрий
def RotateQ(p, q):
    r = np.cos(np.pi / q + np.pi / p) * (np.cos(np.pi / q) ** 2 - np.sin(np.pi / p) ** 2) ** (-0.5)
    alpha = 0.
    betha = 2 * np.pi / q
    theta = np.log((r + 1) / (1 - r))
    R = [
        [np.cos(betha), -np.sin(betha), 0],
        [np.sin(betha), np.cos(betha), 0],
        [0, 0, 1]
    ]
    return np.array(make_Q(alpha, theta)) @ np.array(R) @ np.array(np.linalg.inv(make_Q(alpha, theta)))


#Adjacancey = True Edge, False - vertex
#Алгоритм из книжки
#Входные: матрица предыдущего шага, кол-во слоев, По ребру пришли или по вершине, пэ, ку, ворд - задумывался как способ кодировки цвета
def replicate(initial_tran, LayersToDo, AdjacancyType, p, q, word, OMAE):
    if(LayersToDo <= 0):
        return
    Draw(p, q, initial_tran, word, OMAE, LayersToDo)
    if(AdjacancyType):
        ExposedEdges = p - 3
        RotateCenter = mult_mat3(initial_tran, Rotate3P(p, q))
        word[0] += 3
    else:
        ExposedEdges = p - 2
        RotateCenter = mult_mat3(initial_tran, Rotate2P(p, q))
        word[0] += 2
    for i in range(1, ExposedEdges + 1):
        RotateVertex = mult_mat3(RotateCenter, RotateQ(p, q))
        word [1] += 1
        replicate(RotateVertex, LayersToDo - 1, True, p, q, word, OMAE)
        if i < ExposedEdges:
            VertexPgons = q - 3
        else:
            VertexPgons = q - 4
        for j in range(1, VertexPgons + 1):
            RotateVertex = mult_mat3(RotateVertex, RotateQ(p, q))
            word[1] += 1
            replicate(RotateVertex, LayersToDo - 1, False, p, q, word, OMAE)
        RotateCenter = mult_mat3(RotateCenter, RotateP(p, q))
        word[0] += 1

def RotateQ_2(p, q):
    angle = 2 * np.pi / p
    x = np.cos(angle)
    y = np.sin(angle)
    r = np.cos(np.pi / q + np.pi / p) * (np.cos(np.pi / q) ** 2 - np.sin(np.pi / p) ** 2) ** (-0.5)
    x *= r
    y *= r
    B = [r, 0.]
    tmp = misha_arc([x,y], B)
    alpha = np.pi / p
    betha = np.pi
    # print(tmp)
    # print((tmp[0] ** 2 + tmp[1] ** 2 - tmp[4] ** 2) ** 0.5)
    Rad = (tmp[0] ** 2 + tmp[1] ** 2) ** 0.5 - tmp[4]
    theta = np.log((Rad + 1) / (1 - Rad))
    R = [
        [np.cos(betha), -np.sin(betha), 0],
        [np.sin(betha), np.cos(betha), 0],
        [0, 0, 1]
    ]
    return np.array(make_Q(alpha, theta)) @ np.array(R) @ np.array(np.linalg.inv(make_Q(alpha, theta)))



def Curve(pts, pts_2, p, q):
    theta = np.arange(0, 1, 0.01)
    fig = plt.figure()
    ax = fig.gca()
    cur_x = []
    cur_y = []
    cur_x_2 = []
    cur_y_2 = []
    for i in range(len(pts) - 1):
        cur_x.append(pts[i][0] + (pts[i + 1][0] - pts[i][0]) * theta)
        cur_y.append(pts[i][1] + (pts[i + 1][1] - pts[i][1]) * theta)
    X = np.concatenate(cur_x)
    Y = np.concatenate(cur_y)
    a = 2 * np.pi / p
    X_2 = X * np.cos(a) - Y * np.sin(a)
    Y_2 = X * np.sin(a) + Y * np.cos(a)

    for i in range(len(pts_2) - 1):
        cur_x_2.append(pts_2[i][0] + (pts_2[i + 1][0] - pts_2[i][0]) * theta)
        cur_y_2.append(pts_2[i][1] + (pts_2[i + 1][1] - pts_2[i][1]) * theta)
    X_arc = np.concatenate(cur_x_2)
    Y_arc = np.concatenate(cur_y_2)
    X_kek = []
    Y_kek = []
    
    for i in range(len(X_arc)):
        tmp = change_to_circle(mult_mat_vec(RotateQ_2(p, q), change_to_r3([X_arc[i], Y_arc[i]])))
        X_kek.append(tmp[0])#ТУТ УСКОРИТЬ
        Y_kek.append(tmp[1])
    X_draw = np.concatenate([X, X_arc, X_kek[::-1], X_2[::-1]])
    Y_draw = np.concatenate([Y, Y_arc, Y_kek[::-1], Y_2[::-1]])

    X_res = []
    Y_res = []

    for i in range(p):
        a = i * 2 * np.pi / p
        X_res.append(X_draw * np.cos(a) - Y_draw * np.sin(a))
        Y_res.append(X_draw * np.sin(a) + Y_draw * np.cos(a))
        # plt.fill(X_draw * np.cos(a) - Y_draw * np.sin(a), X_draw * np.sin(a) + Y_draw * np.cos(a))
        # print(i) 
        # plt.show()

    return [X_res, Y_res]


def Draw(p, q, T, word, OMAE, lw_width):
    x_draw = []
    y_draw = []
    for i in range(p):
        for j  in range(len(OMAE[1][i])):
            tmp = change_to_circle(mult_mat_vec(T, change_to_r3([OMAE[0][i][j], OMAE[1][i][j]])))
            x_draw.append(tmp[0])
            y_draw.append(tmp[1])
        if (lw_width < 3):
            plt.fill(x_draw, y_draw, color = decode(word, p, q), lw=0.1)
        else:
            plt.fill(x_draw, y_draw, color = decode(word, p, q), lw=0.4)
        
        word[0] += 1
        x_draw = []
        y_draw = []


def Draw_low(p, q, T, word):
    ar = []
    xf = []
    yf = []
    for i in range(p + 1):
        angle = 2 * i * np.pi / p
        x = np.cos(angle)
        y = np.sin(angle)
        r = np.cos(np.pi / q + np.pi / p) * (np.cos(np.pi / q) ** 2 - np.sin(np.pi / p) ** 2) ** (-0.5)
        x *= r
        y *= r
        pt = change_to_circle(mult_mat_vec(T, change_to_r3([x, y])))
        if i != 0:
            tmp = misha_arc(pt, old)
            if(tmp[5] == -1):
                theta = np.arange(tmp[2], tmp[3], (tmp[3] - tmp[2]) / 100)
                X = tmp[0] + tmp[4] * np.cos(theta)
                Y = tmp[1] + tmp[4] * np.sin(theta)
                xf = np.hstack((xf, X))
                yf = np.hstack((yf, Y))
            elif(tmp[5] == 1):                
                theta = np.arange(tmp[3], tmp[2], (tmp[2] - tmp[3]) / 100)
                X = tmp[0] + tmp[4] * np.cos(theta)
                Y = tmp[1] + tmp[4] * np.sin(theta)
                xf = np.hstack((xf, X))
                yf = np.hstack((yf, Y))
            else:
                theta = np.arange(0, 1, 0.01)
                X = old[0] + (-old[0] + pt[0]) * theta
                Y = old[1] + (-old[1] + pt[1]) * theta
                xf = np.hstack((xf, X))
                yf = np.hstack((yf, Y))               
        old = pt
    plt.fill(xf, yf, color = decode(word, p, q))



def make_pt_q(p, q):
    angle = 2 * np.pi / p
    x = np.cos(angle)
    y = np.sin(angle)
    r = np.cos(np.pi / q + np.pi / p) * (np.cos(np.pi / q) ** 2 - np.sin(np.pi / p) ** 2) ** (-0.5)
    x *= r
    y *= r
    B = [r, 0.]
    tmp = misha_arc([x,y], B)  
    angle = np.pi / p
    x = np.cos(angle)
    y = np.sin(angle)
    x *= (tmp[0] ** 2 + tmp[1] ** 2) ** 0.5 - tmp[4]
    y *= (tmp[0] ** 2 + tmp[1] ** 2) ** 0.5 - tmp[4]
    return (x,y)
        
#P = A
#Q = B

#Слово по которому пришли в заданный многоугольник отображает в цвет
def decode(word, p, q):
    number_A = word[0] % p
    number_B = word[1] % q
    return [0.1 * number_B, 0.3 * number_A / p, 0.1 + (0.8 * number_B * number_A) / (p * q), 0.3]


#Тут будет такая функция: рисует базовую область. а ты в ней точки тыкаешь и круто все
if __name__ == "__main__":  
    an = np.linspace(0, 2 * np.pi, 100)
#С этими цифрами можно играться
    p = int(input('Enter p:'))
    q = int(input('Enter q:'))
    n = int(input('Enter number of layers:')) + 1
    #n = 3 #Количество слоев прорисовки, ну или ему пропорциональное число
    curve_length = int(input('Введите количество звеньев ломанной_1 и тыкните столько точек. Она соединит точку (0, 0) и точку (1, 0):'))
    Draw_low(p, q, E(), [1, 0])
    plt.plot(0 +  0.001 * np.cos(an), 0 + 0.001 * np.sin(an))
    plt.plot(np.cos(np.pi / q + np.pi / p) * (np.cos(np.pi / q) ** 2 - np.sin(np.pi / p) ** 2) ** (-0.5) +  0.001 * np.cos(an), 0 + 0.001 * np.sin(an))
    input_pts = plt.ginput(curve_length, timeout = -1)
    input_pts.insert(0, (0., 0.))
    input_pts.insert(curve_length + 1, (np.cos(np.pi / q + np.pi / p) * (np.cos(np.pi / q) ** 2 - np.sin(np.pi / p) ** 2) ** (-0.5), 0.))
    pts_draw = np.asarray(input_pts)
    plt.plot(pts_draw[:,0], pts_draw[:,1])
    pt_q = make_pt_q(p, q)
    plt.plot(pt_q[0] + 0.001 * np.cos(an), pt_q[1]+ 0.001 * np.sin(an))
    plt.draw()
    curve_length_2 = int(input('Введите количество звеньев ломанной_2. Она соединит точку (1, 0) и сердину дуги:'))
    input_pts_2 = plt.ginput(curve_length_2, timeout = -1)
    input_pts_2.insert(0, (np.cos(np.pi / q + np.pi / p) * (np.cos(np.pi / q) ** 2 - np.sin(np.pi / p) ** 2) ** (-0.5), 0.))
    
    input_pts_2.insert(curve_length_2 + 1,  pt_q)


    OMAE = Curve(input_pts, input_pts_2, p, q) #массив содержащий ломанную
    RotateCenter = E()
    word = [1, 0]
    Draw(p, q, E(), word,  OMAE, 3)
    start_time = time.time()
    for i in range(1, p + 1):
        RotateVertex = mult_mat3(RotateCenter, RotateQ(p, q))
        word[1] += 1
        replicate(RotateVertex, n - 2, False, p, q, word, OMAE)
        for j in range(1, q - 2):
            RotateVertex = mult_mat3(RotateVertex, RotateQ(p, q))
            word[1] += 1
            replicate(RotateVertex, n - 2, False, p, q, word, OMAE)
        RotateCenter = mult_mat3(RotateCenter, RotateP(p, q))
        word[0] += 1
    plt.xlim(-1.2, 1.2)
    plt.ylim(-1.2, 1.2)
    
    plt.plot( np.cos(an), np.sin(an))
    plt.show()
