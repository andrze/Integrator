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
    void on_alphaBox_valueChanged(double);
    void on_uBox_valueChanged(double);
    void on_tauBox_valueChanged(double);
    void on_ZSpinBox_valueChanged(double);
    void on_YSpinBox_valueChanged(double);
    void on_startTimeBox_valueChanged(double);
    void on_endTimeBox_valueChanged(double);
    void on_deltaTBox_valueChanged(double);
    void on_asenseBox_valueChanged(int arg1);
    void on_usenseBox_valueChanged(int arg1);
    void on_tsenseBox_valueChanged(int arg1);
    //void on_fpButton_clicked();
    void on_pushButton_clicked();
};

#endif // MAINWINDOW_H
