#include "mainwindow.h"


#include <QDockWidget>


MainWindow::MainWindow(QWidget *parent)
	: QMainWindow(parent)
{
	coefficientsWidget = new CoefficientsWidget(this);
	impulsePlot = new QLabel;
	fftPlot = new QLabel;
	filterdesignWidget = new FilterDesignWidget(this);
	setCentralWidget(filterdesignWidget);

	auto fftPlotDock = new QDockWidget("Frequency Response", this);
	auto impulsePlotDock = new QDockWidget("Impulse Response", this);
	auto coeffsDock = new QDockWidget("Coefficients", this);

	fftPlotDock->setWidget(fftPlot);
	fftPlotDock->setAllowedAreas(Qt::AllDockWidgetAreas);
	addDockWidget(Qt::TopDockWidgetArea, fftPlotDock);

	impulsePlotDock->setWidget(impulsePlot);
	impulsePlotDock->setAllowedAreas(Qt::AllDockWidgetAreas);
	addDockWidget(Qt::TopDockWidgetArea, impulsePlotDock);

	coeffsDock->setWidget(coefficientsWidget);
	coeffsDock->setAllowedAreas(Qt::AllDockWidgetAreas);
	addDockWidget(Qt::RightDockWidgetArea, coeffsDock);

	filterdesignWidget->setFftPlot(fftPlot);
	filterdesignWidget->setImpulsePlot(impulsePlot);

	connect(filterdesignWidget, &FilterDesignWidget::createdCoefficients, coefficientsWidget, &CoefficientsWidget::setCoefficients);
}

MainWindow::~MainWindow()
{
}

