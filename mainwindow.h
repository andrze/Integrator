#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "integrator.h"

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
    Integrator integrator;

    void update();
    PlotSet integrate();
    PlotSet results;

private slots:

    void on_fileButton_clicked();
    void on_saveButton_clicked();
    void on_kappaBox_valueChanged(double);
    void on_lambdaBox_valueChanged(double);
    void on_startTimeBox_valueChanged(double);
    void on_endTimeBox_valueChanged(double);
    void on_deltaTBox_valueChanged(double);
    void on_etaBox_valueChanged(double);
    void on_uBox_valueChanged(double);
    void on_ksenseBox_valueChanged(int arg1);
    void on_usenseBox_valueChanged(int arg1);
    void on_lsenseBox_valueChanged(int arg1);
};

#endif // MAINWINDOW_H
