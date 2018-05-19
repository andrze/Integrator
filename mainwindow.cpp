#include <vector>
#include <functional>
#include <fstream>
#include <QString>
#include <QFileDialog>
#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "integrator.h"
#include "realvector.h"

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

double pos_der(RealVector x){
    return x.coords[1];
}

double vel_der(RealVector x){
    return -x.coords[0];
}

void MainWindow::on_runButton_clicked()
{
    std::vector<double> coords({ui->positionBox->value(), ui->velocityBox->value()});
    RealVector start(coords);

    EquationSet eq;
    eq.derivatives = std::vector<std::function<double(RealVector)> >{pos_der, vel_der};

    Integrator i(eq);

    PlotSet p = i.integrate(ui->startTimeBox->value(),
                            ui->endTimeBox->value(),
                            ui->deltaTBox->value(), start);

    std::ofstream file(ui->filenameEdit->text().toStdString());
    file << p;

}
