#ifndef FILTERDESIGNWIDGET_H
#define FILTERDESIGNWIDGET_H

#include <QLabel>
#include <QSlider>
#include "slider.h"
#include <QPushButton>
#include <QLineEdit>
#include <QComboBox>
#include <QWidget>
#include <QCheckBox>
#include <QEvent>
#include <QResizeEvent>
#include <QDebug>

#include <functional>
#include <optional>

struct PlotInfo
{
	bool connectDots{true};
	QColor plotColor{0,128,0,255};
	QColor gridColor{192, 192, 192, 192};
	QColor txtColor{QColor{64,0,0,255}};
	QBrush backgroundColor{QColor{255,255,255,255}};

	std::optional<std::pair<double, double>> range;
	std::optional<size_t> size;

	std::pair<bool, bool> gridlines{true, true};
	std::pair<int, int> precision{0,1};

	bool showCrossings{false};
};

class PixmapResizer : public QObject
{
	Q_OBJECT

public:
	PixmapResizer(QPixmap& p, QObject *parent = nullptr) : _pm(std::ref(p)), QObject(parent){}

protected:
	bool eventFilter(QObject *obj, QEvent *event) override {
		if(event->type() == QEvent::Resize) {

			QResizeEvent *resizeEvent = static_cast<QResizeEvent*>(event);
			qDebug() << resizeEvent->size();
			QPixmap &p = _pm.get();
			p = std::move(p.scaled(resizeEvent->size()));
			return true;
		}
		return QObject::eventFilter(obj, event);
	}

	// store a reference to pixmap
	std::reference_wrapper<QPixmap> _pm;
};


class FilterDesignWidget : public QWidget
{
	Q_OBJECT

	enum FilterType {
		LowPass,
		HighPass
	};

public:
	explicit FilterDesignWidget(QWidget *parent = nullptr);

	template<typename FloatType>
	std::vector<FloatType> makeFilterCoefficients();
	void setImpulsePlot(QLabel *newImpulsePlot);

	void setFftPlot(QLabel *newFftPlot);

public slots:
	void makeFilter();

signals:
	void createdCoefficients(const std::vector<double> &coefficients);

private:
	QLabel *impulsePlot{nullptr};
	QLabel *fftPlot{nullptr};
//	QPixmap fftPlotpixmap{1024, 768};
//	QPixmap impulsePlotPixmap{1024, 768};

	QPixmap fftPlotpixmap{512, 384};
	QPixmap impulsePlotPixmap{512, 384};

	QPushButton *generateButton{nullptr};
	QComboBox *filterType{nullptr};
	Slider *startFrequency{nullptr};
	Slider *transitionBand{nullptr};
	Slider *dynamicRange{nullptr};
	QCheckBox *minPhaseCheckbox{nullptr};

	double kaiserBeta;

	void plotResponse(const QVector<double> &impulse);
	void plotVector(QPixmap *pm, const std::vector<double> &data, const PlotInfo& plotInfo);

	static bool convertToLinear(std::vector<double> &s, bool fromMagSquared);
	static bool convertToDb(std::vector<double> &s, bool fromMagSquared);
};

#endif // FILTERDESIGNWIDGET_H
