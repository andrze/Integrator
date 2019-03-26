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

    std::pair<PlotSet,int> integrate();
    std::pair<PlotSet,int> integrate(RealVector start_point);
    std::vector<PlotSet > results;
    void reset_axes();
    std::vector<QCustomPlot*> plots;
    std::vector<bool> is_plot_log;
    //std::vector<std::pair<double,double> > critical_line;
    std::atomic<bool> stop;
    std::mutex data_mutex;

    double find_fp(RealVector start_point);

    void plot_result(std::vector<Plot> p);
    void restart_integrator();

signals:
    void calculated(PlotSet* res);

private slots:

    void plot_result(PlotSet* res);
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
    void on_dimensionBox_valueChanged(double arg1);
};

#endif // MAINWINDOW_H
