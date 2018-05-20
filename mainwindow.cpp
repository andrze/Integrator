#include <vector>
#include <functional>
#include <fstream>
#include <QString>
#include <QFileDialog>
#include <algorithm>
#include <cmath>
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "integrator.h"
#include "realvector.h"
#include "physics.cpp"


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::update(){
    std::vector<double> coords({ui->positionBox->value(), ui->velocityBox->value()});
    RealVector start(coords);

    EquationSet eq;
    eq.derivatives = std::vector<std::function<double(RealVector)> >{k_der, l_der};

    Integrator i(eq);

    PlotSet p = i.integrate(ui->startTimeBox->value(),
                            ui->endTimeBox->value(),
                            ui->deltaTBox->value(), start);

    //std::ofstream file(ui->filenameEdit->text().toStdString());
    //file << p;

    // generate some data:
    std::vector<std::vector<double> > results = p.transpose();

    std::vector<QCustomPlot*> plots{ui->posPlot, ui->velPlot};

    QVector<double> x=QVector<double>::fromStdVector(p.time_exp());

    for(size_t i=0; i<plots.size(); i++){
        auto plot = plots[i];
        std::vector<double> vals = results[i];
        QVector<double> y = QVector<double>::fromStdVector(vals);
        plot->addGraph();
        plot->graph(0)->setData(x, y);
        plot->xAxis->setRange(0.,1.5);
        auto min = std::min_element(vals.begin(), vals.end());
        auto max = std::max_element(vals.begin(), vals.end());
        plot->yAxis->setRange((*min)-0.25, (*max)+0.25),
        plot->replot();
    }
}

void MainWindow::on_fileButton_clicked()
{
    {
        std::string filename = QFileDialog::getSaveFileName(this, tr("Save file"), "", tr("CSV files (*.csv)")).toStdString();
        if( filename.size() == 0 ){
            return;
        }
        if( filename.size() < 4 || filename.substr(filename.size()-4,4) != ".csv" ){
            filename += ".csv";
        }
        ui->runButton->setEnabled(true);
        ui->filenameEdit->setText(tr(&filename[0]));
    }
}


void MainWindow::on_runButton_clicked()
{
    this->update();
}

void MainWindow::on_velocityBox_valueChanged(double)
{
    this->update();
}


void MainWindow::on_positionBox_valueChanged(double)
{
    this->update();
}

void MainWindow::on_startTimeBox_valueChanged(double)
{
    this->update();
}

void MainWindow::on_endTimeBox_valueChanged(double)
{
    this->update();
}

void MainWindow::on_deltaTBox_valueChanged(double)
{
    this->update();
}
