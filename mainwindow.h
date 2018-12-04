#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include "qcustomplot.h"
#include "integrator.h"
#include <thread>
#include <mutex>

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
    //std::unique_ptr<std::thread> worker;

    std::vector<Plot> integrate();
    std::vector<Plot> integrate(RealVector start_point);
    std::vector<std::vector<Plot> > results;
    void reset_xAxis();
    std::vector<QCustomPlot*> plots;
    //std::vector<std::pair<double,double> > critical_line;
    std::atomic<bool> stop;
    std::mutex data_mutex;

    double find_fp(RealVector start_point);

    void plot_result(std::vector<Plot> p);

signals:
    void calculated(std::vector<Plot>* res);

private slots:

    void plot_result(std::vector<Plot>* res);
    void on_fileButton_clicked();
    void on_saveButton_clicked();
    void on_tsenseBox_valueChanged(int arg1);
    void on_pushButton_clicked();
    void on_clearButton_clicked();
    void on_xMinSpinBox_valueChanged(double);
    void on_xMaxSpinBox_valueChanged(double);
    void on_cubicRadioButton_clicked();
    void on_hexRadioButton_clicked();
    void on_fpButton_clicked();
    void on_CLButton_clicked();
    void on_stopButton_clicked();
};

#endif // MAINWINDOW_H
