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

    std::vector<Plot> integrate();
    std::vector<Plot>  results;
    void reset_xAxis();
    std::vector<QCustomPlot*> plots;

private slots:

    void on_fileButton_clicked();
    void on_saveButton_clicked();
    void on_tsenseBox_valueChanged(int arg1);
    //void on_fpButton_clicked();
    void on_pushButton_clicked();
    void on_clearButton_clicked();
    void on_xMinSpinBox_valueChanged(double);
    void on_xMaxSpinBox_valueChanged(double);
};

#endif // MAINWINDOW_H
