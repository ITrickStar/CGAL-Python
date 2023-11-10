from PyQt5.QtWidgets import QApplication, QMainWindow, QOpenGLWidget
from PyQt5.QtGui import QMouseEvent
from PyQt5.QtCore import Qt
from OpenGL.GL import *
from OpenGL.GLU import gluPerspective, gluProject
from CGAL.CGAL_Kernel import Point_3
from CGAL.CGAL_Kernel import Vector_3


class MyOpenGLWidget(QOpenGLWidget):
    def __init__(self, parent=None):
        super(MyOpenGLWidget, self).__init__(parent)
        self.vertices = [Point_3(0, 0, 0), Point_3(1, 1, 1), Point_3(2, 2, 2)]
        self.selected_vertex_index = None
        self.mouse_x, self.mouse_y = 0, 0

    def initializeGL(self):
        glEnable(GL_DEPTH_TEST)

    def resizeGL(self, w, h):
        glViewport(0, 0, w, h)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(45, w / h, 1, 100)
        glMatrixMode(GL_MODELVIEW)

    def paintGL(self):
        #glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glLoadIdentity()

        glTranslatef(0, 0, -5)

        glColor3f(1, 1, 1)
        glPointSize(5.0)
        glBegin(GL_POINTS)

        for vertex in self.vertices:
            glVertex3f(float(vertex.x()), float(vertex.y()), float(vertex.z()))

        glEnd()

        # Обновляем экран
        self.update()

    def mousePressEvent(self, event: QMouseEvent):
        if event.button() == Qt.LeftButton:
            # При нажатии левой кнопки мыши, проверяем, близко ли курсор к какой-либо вершине
            for i, vertex in enumerate(self.vertices):
                screen_pos = self.project_to_screen(vertex)
                if abs(screen_pos[0] - event.x()) < 5 and abs(screen_pos[1] - (self.height() - event.y())) < 5:
                    self.selected_vertex_index = i
                    print(f"Point selected event: ({self.vertices[self.selected_vertex_index]})")
                    break

    def mouseReleaseEvent(self, event: QMouseEvent):
        if event.button() == Qt.LeftButton:
            self.selected_vertex_index = None

    def mouseMoveEvent(self, event: QMouseEvent):
        print(f"Mouse move event: ({event.x()}, {event.y()})")
        
        if self.selected_vertex_index is not None:
            # If a vertex is selected, update its coordinates based on mouse movement
            delta_x = event.x() - self.mouse_x
            delta_y = (self.height() - event.y()) - self.mouse_y

            # Create a Vector_3 with the movement
            movement_vector = Vector_3(delta_x / 100, delta_y / 100, 0)

            # Update the vertex by adding the movement vector
            self.vertices[self.selected_vertex_index] = self.vertices[self.selected_vertex_index] + movement_vector
            print(f"Point move event: ({self.vertices[self.selected_vertex_index]})")

        self.mouse_x = event.x()
        self.mouse_y = self.height() - event.y()

        # Update the screen
        self.update()

    def project_to_screen(self, vertex):
        model_view = glGetDoublev(GL_MODELVIEW_MATRIX)
        projection = glGetDoublev(GL_PROJECTION_MATRIX)
        viewport = glGetIntegerv(GL_VIEWPORT)

        screen_pos = gluProject(float(vertex.x()), float(vertex.y()), float(vertex.z()), model_view, projection, viewport)
        return screen_pos

class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.setGeometry(100, 100, 800, 600)
        self.central_widget = MyOpenGLWidget(self)
        self.setCentralWidget(self.central_widget)

if __name__ == '__main__':
    app = QApplication([])
    main_window = MainWindow()
    main_window.show()
    app.exec_()
