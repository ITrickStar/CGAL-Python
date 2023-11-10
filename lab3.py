# bspline3.py
from __future__ import print_function
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQGLViewer import *
from OpenGL.GL import *
import random


class BSpline:
    def __init__(self, reference_points, discrete_num=10, degree=False, closed=False):
        self.points = reference_points
        self.d_num = int(discrete_num)
        self.closed = closed
        self.degree = degree

        # Генерация коэффициентов для сгенеренных вершин B-сплайна 3 порядка
        self.coefs = []
        for i in range(self.d_num):
            if degree:
                spline_segm_coef = self.calc_spline3_coef(i / self.d_num)
            else:
                spline_segm_coef = self.calc_spline2_coef(i / self.d_num)

            self.coefs.append(spline_segm_coef)

    @staticmethod
    def clamp(value, minval, maxval):
        return max(minval, min(value, maxval))

    @staticmethod
    def calc_spline2_coef(t):
        coefs = [0, 0, 0]
        coefs[0] = (1.0 - t) * (1.0 - t) / 2.0
        coefs[1] = (-2.0 * t * t + 2.0 * t + 1) / 2
        coefs[2] = t * t / 2.0
        return coefs

    @staticmethod
    def calc_spline3_coef(t):
        coefs = [0, 0, 0, 0]
        coefs[0] = (1.0 - t) * (1.0 - t) * (1.0 - t) / 6.0
        coefs[1] = (3.0 * t * t * t - 6.0 * t * t + 4) / 6.0
        coefs[2] = (-3.0 * t * t * t + 3 * t * t + 3 * t + 1) / 6.0
        coefs[3] = t * t * t / 6.0
        return coefs

    def draw_spline_curve(self):
        if not self.closed:
            segmentsCount = len(self.points) - 1
            glBegin(GL_LINE_STRIP)
        else:
            segmentsCount = len(
                self.points
            )  # Сегмент между первой и последней вершиной
            glBegin(GL_LINE_LOOP)
        glColor3f(1.0, 1.0, 0.0)
        for i in range(segmentsCount):
            self.draw_glvertex_for_one_segment_of_spline(i)
        glEnd()

    def draw_glvertex_for_one_segment_of_spline(self, segment_id):
        pNum = len(self.points)
        # Вычисление номеров вершин в списке вершин для построения сплайна
        if not self.closed:
            p0 = self.clamp(segment_id - 1, 0, pNum - 1)
            p1 = self.clamp(segment_id, 0, pNum - 1)
            p2 = self.clamp(segment_id + 1, 0, pNum - 1)
            p3 = self.clamp(segment_id + 2, 0, pNum - 1)
        else:
            p0 = (segment_id - 1 + pNum) % pNum
            p1 = (segment_id + pNum) % pNum
            p2 = (segment_id + 1 + pNum) % pNum
            p3 = (segment_id + 2 + pNum) % pNum
        # По заранее вычисленным коэффициентам
        # вычисляем промежуточные точки сплайна
        # и выводим их в OpenGL
        if self.degree:
            for i in range(self.d_num):
                x = (
                    self.coefs[i][0] * self.points[p0][0]
                    + self.coefs[i][1] * self.points[p1][0]
                    + self.coefs[i][2] * self.points[p2][0]
                    + self.coefs[i][3] * self.points[p3][0]
                )
                y = (
                    self.coefs[i][0] * self.points[p0][1]
                    + self.coefs[i][1] * self.points[p1][1]
                    + self.coefs[i][2] * self.points[p2][1]
                    + self.coefs[i][3] * self.points[p3][1]
                )
                z = (
                    self.coefs[i][0] * self.points[p0][2]
                    + self.coefs[i][1] * self.points[p1][2]
                    + self.coefs[i][2] * self.points[p2][2]
                    + self.coefs[i][3] * self.points[p3][2]
                )
                glVertex3f(x, y, z)

        else:
            for i in range(self.d_num):
                x = (
                    self.coefs[i][0] * self.points[p0][0]
                    + self.coefs[i][1] * self.points[p1][0]
                    + self.coefs[i][2] * self.points[p2][0]
                )
                y = (
                    self.coefs[i][0] * self.points[p0][1]
                    + self.coefs[i][1] * self.points[p1][1]
                    + self.coefs[i][2] * self.points[p2][1]
                )
                z = (
                    self.coefs[i][0] * self.points[p0][2]
                    + self.coefs[i][1] * self.points[p1][2]
                    + self.coefs[i][2] * self.points[p2][2]
                )
                glVertex3f(x, y, z)


class BSpline2:
    def __init__(self, reference_points, discrete_num=10, closed=False):
        self.points = reference_points
        self.d_num = int(discrete_num)
        self.closed = closed

        # Генерация коэффициентов для сгенеренных вершин B-сплайна 3 порядка
        self.coefs = []
        for i in range(self.d_num):
            spline_segm_coef = self.calc_spline3_coef(i / self.d_num)
            self.coefs.append(spline_segm_coef)


class Edges:
    def __init__(self, reference_points):
        self.points = reference_points

    def draw(self):
        for ind in range(1, len(self.points)):
            glBegin(GL_LINES)
            glColor3f(*(1, 0, 1))
            glVertex3f(*self.points[ind - 1])
            glVertex3f(*self.points[ind])
            glEnd()


class Viewer(QGLViewer):
    def __init__(self, parent=None):
        QGLViewer.__init__(self, parent)

    def draw(self):
        # points = ((0, 0, 0), (0, 3, 0), (1, 3, 0), (1, 1, 0), (2, 1, 0), (3, 2, 0), (3, 0, 0))
        points = []
        for i in range(7):
            point = (random.random(), random.random(), random.random())
            points.append(point)

        edges = Edges(points)
        edges.draw()
        spline2 = BSpline(points, 20, False, False)
        spline2.draw_spline_curve()
        spline3 = BSpline(points, 20, True, False)
        spline3.draw_spline_curve()


def main():
    qapp = QApplication([])
    viewer = Viewer()
    viewer.show()
    qapp.exec_()


if __name__ == "__main__":
    main()
