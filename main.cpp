#include <QApplication>
#include "integrator.h"
#include "realvector.h"
#include "mainwindow.h"
#include <fenv.h>

int main(int argc, char *argv[])
{
    feraiseexcept(FE_INVALID | FE_OVERFLOW);
    QApplication a(argc, argv);
    MainWindow w;
    w.show();
    return a.exec();
}
