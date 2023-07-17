#ifndef MAINWINDOW_H
#define MAINWINDOW_H


#include "filterdesignwidget.h"
#include "coefficientswidget.h"

#include <QMainWindow>


class MainWindow : public QMainWindow
{
	Q_OBJECT

public:
	MainWindow(QWidget *parent = nullptr);
	~MainWindow();


private:
	FilterDesignWidget *filterdesignWidget;
	CoefficientsWidget *coefficientsWidget;
	QLabel *impulsePlot;
	QLabel *fftPlot;

};
#endif // MAINWINDOW_H
