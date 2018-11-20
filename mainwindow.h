#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "qcustomplot.h"
#include "integrator.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private:
    Ui::MainWindow *ui;
    Integrator integrator;

    void update();
    std::vector<Plot> integrate();
    std::vector<Plot>  results;
    void reset_xAxis();
    std::vector<QCustomPlot*> plots;

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
    void on_clearButton_clicked();
    void on_xMinSpinBox_valueChanged(double);
    void on_xMaxSpinBox_valueChanged(double);
};

#endif // MAINWINDOW_H
