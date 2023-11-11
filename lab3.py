from __future__ import print_function
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQGLViewer import *
from OpenGL.GL import *
import random as rn
from OpenGL.GLU import gluPerspective, gluProject
from CGAL.CGAL_Kernel import Point_3
from CGAL.CGAL_Kernel import Vector_3


def clamp(value, minval, maxval):
    return max(minval, min(value, maxval))


class BSpline:
    def __init__(
        self,
        reference_points,
        power=3,
        discrete_num=10,
        closed=False,
        color=[1.0, 0.0, 1.0],
    ):
        self.points = reference_points
        self.power = power
        self.d_num = int(discrete_num)
        self.closed = False
        self.color = color

        # Генерация коэффициентов для сгенеренных вершин B-сплайна 3 порядка
        self.coefs = []
        if power == 2:
            for i in range(self.d_num):
                spline_segm_coef = self.calc_spline2_coef(i / self.d_num)
                self.coefs.append(spline_segm_coef)
        if power == 3:
            for i in range(self.d_num):
                spline_segm_coef = self.calc_spline3_coef(i / self.d_num)
                self.coefs.append(spline_segm_coef)

    def calc_spline2_coef(self, t):
        coefs = [0] * 3
        coefs[0] = (1.0 - t) * (1.0 - t) / 2.0
        coefs[1] = (-2.0 * t * t + 2.0 * t + 1.0) / 2.0
        coefs[2] = t * t / 2.0
        return coefs

    def calc_spline3_coef(self, t):
        coefs = [0] * 4
        coefs[0] = (1.0 - t) * (1.0 - t) * (1.0 - t) / 6.0
        coefs[1] = (3.0 * t * t * t - 6.0 * t * t + 4) / 6.0
        coefs[2] = (-3.0 * t * t * t + 3 * t * t + 3 * t + 1) / 6.0
        coefs[3] = t * t * t / 6.0
        return coefs

    def draw_spline_curve(self):
        glPointSize(2.5)
        glBegin(GL_POINTS)
        glColor3f(0.0, 1.0, 0.0)
        for point in self.points:
            glVertex3f(point.x(), point.y(), point.z())
        glEnd()
        glBegin(GL_LINE_STRIP)
        glColor3f(1.0, 0.0, 0.0)
        for point in self.points:
            glVertex3f(point.x(), point.y(), point.z())
        glEnd()

        if not self.closed:
            segmentsCount = len(self.points) - 1
            glBegin(GL_LINE_STRIP)
        else:
            segmentsCount = len(
                self.points
            )  # Сегмент между первой и последней вершиной
            glBegin(GL_LINE_LOOP)
        glColor3f(self.color[0], self.color[1], self.color[2])
        for i in range(segmentsCount):
            self.draw_glvertex_for_one_segment_of_spline(i)
        glEnd()

    def draw_glvertex_for_one_segment_of_spline(self, segment_id):
        pNum = len(self.points)
        # Вычисление номеров вершин в списке вершин для построения сплайна

        if self.power == 2:
            if not self.closed:
                p0 = clamp(segment_id - 1, 0, pNum - 1)
                p1 = clamp(segment_id, 0, pNum - 1)
                p2 = clamp(segment_id + 1, 0, pNum - 1)
            else:
                p0 = (segment_id - 1 + pNum) % pNum
                p1 = (segment_id + pNum) % pNum
                p2 = (segment_id + 1 + pNum) % pNum
        if self.power == 3:
            if not self.closed:
                p0 = clamp(segment_id - 1, 0, pNum - 1)
                p1 = clamp(segment_id, 0, pNum - 1)
                p2 = clamp(segment_id + 1, 0, pNum - 1)
                p3 = clamp(segment_id + 2, 0, pNum - 1)
            else:
                p0 = (segment_id - 1 + pNum) % pNum
                p1 = (segment_id + pNum) % pNum
                p2 = (segment_id + 1 + pNum) % pNum
                p3 = (segment_id + 2 + pNum) % pNum
        # По заранее вычисленным коэффициентам
        # вычисляем промежуточные точки сплайна
        # и выводим их в OpenGL
        if self.power == 2:
            for i in range(self.d_num):
                x = (
                    self.coefs[i][0] * self.points[p0].x()
                    + self.coefs[i][1] * self.points[p1].x()
                    + self.coefs[i][2] * self.points[p2].x()
                )
                y = (
                    self.coefs[i][0] * self.points[p0].y()
                    + self.coefs[i][1] * self.points[p1].y()
                    + self.coefs[i][2] * self.points[p2].y()
                )
                z = (
                    self.coefs[i][0] * self.points[p0].z()
                    + self.coefs[i][1] * self.points[p1].z()
                    + self.coefs[i][2] * self.points[p2].z()
                )

                glVertex3f(x, y, z)
        if self.power == 3:
            for i in range(self.d_num):
                x = (
                    self.coefs[i][0] * self.points[p0].x()
                    + self.coefs[i][1] * self.points[p1].x()
                    + self.coefs[i][2] * self.points[p2].x()
                    + self.coefs[i][3] * self.points[p3].x()
                )
                y = (
                    self.coefs[i][0] * self.points[p0].y()
                    + self.coefs[i][1] * self.points[p1].y()
                    + self.coefs[i][2] * self.points[p2].y()
                    + self.coefs[i][3] * self.points[p3].y()
                )
                z = (
                    self.coefs[i][0] * self.points[p0].z()
                    + self.coefs[i][1] * self.points[p1].z()
                    + self.coefs[i][2] * self.points[p2].z()
                    + self.coefs[i][3] * self.points[p3].z()
                )
                glVertex3f(x, y, z)


# Make spline
# points = ((0,0,0),(0,3,0),(1,3,0),(1,1,0),(2,1,0),(3,2,0),(3,0,0))


class Viewer(QGLViewer):
    def __init__(self, parent=None):
        QGLViewer.__init__(self, parent)
        self.spline2d = BSpline([])
        self.spline3d = BSpline([])
        self.points = []
        self.n = 3
        self.selected_vertex_index = None
        self.mouse_x, self.mouse_y = 0, 0
        self.closed = False

    def draw(self):
        self.spline2d.draw_spline_curve()
        self.spline3d.draw_spline_curve()

    def keyPressEvent(self, e):
        modifiers = e.modifiers()
        if e.nativeVirtualKey() == Qt.Key_C:
            self.points.clear()
            print("points cleared")
        elif e.nativeVirtualKey() == Qt.Key_Q:
            self.points.append(Point_3(rn.uniform(-1.0, 1.0), rn.uniform(-1.0, 1.0), 0))
            print("new point generated")
        elif e.nativeVirtualKey() == Qt.Key_W:
            i = rn.randrange(0, len(self.points))
            if len(self.points) > 0:
                self.points[i] = Point_3(
                    rn.uniform(-1.0, 1.0), rn.uniform(-1.0, 1.0), 0
                )
            print("point moved")
        elif e.nativeVirtualKey() == Qt.Key_1:
            self.points = []
            for i in range(self.n):
                self.points.append(
                    Point_3(rn.uniform(-1.0, 1.0), rn.uniform(-1.0, 1.0), 0)
                )
            self.spline2d = BSpline(
                self.points, 2, 10, self.closed, color=[0.5, 1, 0.5]
            )
            self.spline3d = BSpline(
                self.points, 3, 10, self.closed, color=[0.5, 0, 0.5]
            )
            print("spline2d and spline3d regenerated")
        elif e.nativeVirtualKey() == Qt.Key_2:
            self.points = []
            for i in range(self.n):
                self.points.append(
                    Point_3(rn.uniform(-1.0, 1.0), rn.uniform(-1.0, 1.0), 0)
                )
            self.spline2d = BSpline(
                self.points, 2, 10, self.closed, color=[0.5, 1, 0.5]
            )
            print("spline2d regenerated")
        elif e.nativeVirtualKey() == Qt.Key_3:
            self.points = []
            for i in range(self.n):
                self.points.append(
                    Point_3(rn.uniform(-1.0, 1.0), rn.uniform(-1.0, 1.0), 0)
                )
            self.spline3d = BSpline(
                self.points, 3, 10, self.closed, color=[0.5, 0, 0.5]
            )
            print("spline3d regenerated")
        self.updateGL()

    def project_to_screen(self, vertex):
        model_view = glGetDoublev(GL_MODELVIEW_MATRIX)
        projection = glGetDoublev(GL_PROJECTION_MATRIX)
        viewport = glGetIntegerv(GL_VIEWPORT)

        screen_pos = gluProject(
            float(vertex.x()),
            float(vertex.y()),
            float(vertex.z()),
            model_view,
            projection,
            viewport,
        )
        return screen_pos

    def mousePressEvent(self, event):
        if event.button() == Qt.LeftButton:
            # При нажатии левой кнопки мыши, проверяем, близко ли курсор к какой-либо вершине
            if len(self.points) > 0:
                for i, vertex in enumerate(self.points):
                    screen_pos = self.project_to_screen(vertex)
                    if (
                        abs(screen_pos[0] - event.x()) < 10
                        and abs(screen_pos[1] - (self.height() - event.y())) < 10
                    ):
                        self.selected_vertex_index = i
                        self.mouse_x = event.x()
                        self.mouse_y = event.y()
                        print(
                            f"Point selected event: ({self.points[self.selected_vertex_index]})"
                        )
                        break

    def mouseMoveEvent(self, event):
        if self.selected_vertex_index is not None:
            # print(f"Mouse move event: ({event.x()}, {event.y()})")
            # If a vertex is selected, update its coordinates based on mouse movement
            delta_x = event.x() - self.mouse_x
            delta_y = self.mouse_y - event.y()
            print("delta_x ", delta_x)
            print("delta_y ", delta_y)

            #     # Create a Vector_3 with the movement
            movement_vector = Vector_3(
                delta_x / self.width() * 3, delta_y / self.height() * 2, 0
            )

            #     # Update the vertex by adding the movement vector
            self.points[self.selected_vertex_index] = (
                self.points[self.selected_vertex_index] + movement_vector
            )
            # print(f"Point move event: ({self.points[self.selected_vertex_index]})")

            self.mouse_x = event.x()
            self.mouse_y = event.y()

            # Update the screen
            self.update()

    def mouseReleaseEvent(self, event):
        if event.button() == Qt.LeftButton:
            self.selected_vertex_index = None
        self.updateGL()


def main():
    qapp = QApplication([])
    viewer = Viewer()
    viewer.show()
    qapp.exec_()


if __name__ == "__main__":
    main()
