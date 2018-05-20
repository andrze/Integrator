#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    ~MainWindow();

private:
    Ui::MainWindow *ui;
    void update();

private slots:

    void on_fileButton_clicked();
    void on_runButton_clicked();
    void on_velocityBox_valueChanged(double);
    void on_positionBox_valueChanged(double);
    void on_startTimeBox_valueChanged(double);
    void on_endTimeBox_valueChanged(double);
    void on_deltaTBox_valueChanged(double);
};

#endif // MAINWINDOW_H
