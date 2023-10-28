# viewer3d_triangle.py
import random
import numpy as np
from PyQt5.QtCore import *
from PyQt5.QtWidgets import *
from PyQGLViewer import *
from OpenGL.GL import *
from OpenGL.GLUT import *

objects_to_draw = []


class MySphere:
    def __init__(self, r, points, color, rand_points=0, xyz=(0, 0, 0)):
        self.vertices = []
        for i in np.linspace(0.0, np.pi, points):
            for j in np.linspace(0.0, 2 * np.pi, points):
                x = r * np.sin(i) * np.cos(j) + xyz[0]
                y = r * np.sin(i) * np.sin(j) + xyz[1]
                z = r * np.cos(i) + xyz[2]
                self.vertices.append((x, y, z))

        for i in np.linspace(0.0, np.pi, rand_points):
            for j in np.linspace(0.0, 2 * np.pi, rand_points):
                phi = np.random.uniform(0, np.pi)
                theta = np.random.uniform(0, 2 * np.pi)
                x = r * np.sin(phi) * np.cos(theta) + xyz[0]
                y = r * np.sin(phi) * np.sin(theta) + xyz[1]
                z = r * np.cos(phi) + xyz[2]
                self.vertices.append((x, y, z))

        self.color = color


class MyCube:
    def __init__(self, size, edge_size, color, rand_points=0, xyz=(0, 0, 0)):
        self.vertices = []
        for x in np.linspace(0.0, size, edge_size):
            for y in np.linspace(0.0, size, edge_size):
                for z in np.linspace(0.0, size, edge_size):
                    self.vertices.append((x + xyz[0], y + xyz[1], z + xyz[2]))
        for x in range(rand_points):
            self.vertices.append(
                (
                    random.uniform(xyz[0], xyz[0] + size),
                    random.uniform(xyz[1], xyz[1] + size),
                    random.uniform(xyz[2], xyz[2] + size),
                )
            )
        self.color = color


class MyObject:
    def __init__(self, vertices, color):
        self.vertices = vertices
        self.color = color


class Viewer(QGLViewer):
    def __init__(self, parent=None):
        QGLViewer.__init__(self, parent)

    # def initializeGL(self):
    #     glEnable(GL_DEPTH_TEST)

    def draw(self):
        glPointSize(5.0)
        glBegin(GL_POINTS)
        for obj in objects_to_draw:
            glColor3f(*obj.color)
            for vertex in obj.vertices:
                glVertex3f(*vertex)
        glEnd()

    def clear_vertices(self):
        # Очистить список вершин
        objects_to_draw.clear()
        self.updateGL()  # Перерисовать сцену

    def keyPressEvent(self, e):
        modifiers = e.modifiers()
        if e.nativeVirtualKey() == Qt.Key_W:
            glPolygonMode(GL_FRONT_AND_BACK, GL_LINE)
        elif e.nativeVirtualKey() == Qt.Key_F:
            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL)
        elif e.nativeVirtualKey() == Qt.Key_C:
            self.clear_vertices()

        self.updateGL()


def main():
    # object = MyObject(
    #     vertices=[(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.5, 1.0, 0.0)],
    #     color=(1.0, 0.0, 0.0),
    # )
    # objects_to_draw.append(object)

    cube = MyCube(
        size=1,
        edge_size=10,
        rand_points=100,
        color=(1.0, 0.0, 0.0),
    )
    objects_to_draw.append(cube)

    # sphere = MySphere(r=0.1, points=100, rand_points=0, color=(0.0, 1.0, 0.0))
    # objects_to_draw.append(sphere)

    qapp = QApplication([])
    viewer = Viewer()
    viewer.show()
    qapp.exec_()


if __name__ == "__main__":
    main()
