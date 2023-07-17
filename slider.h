#ifndef SLIDER_H
#define SLIDER_H

#include <QSlider>
#include <QWidget>
#include <QLineEdit>
#include <QLabel>

#include <QHBoxLayout>
#include <QVBoxLayout>

class Slider : public QWidget
{
	Q_OBJECT
public:
	Slider(Qt::Orientation orientation, const QString &title, QWidget *parent = nullptr, bool titleAbove = true);

public:
	double value() const;
	bool getLogarithmic() const;
	void setLogarithmic(bool newLogarithmic);

	bool getTitleAbove() const;
	void setTitleAbove(bool newTitleAbove);

public slots:
	void setOrientation(Qt::Orientation orientation);
	void setRange(double min, double max);
	void setSteps(int numSteps);
	void setValue(double v);

signals:
	void valueChanged(double value);

private:
	QLabel *label;
	QSlider *m_slider;
	QLineEdit *valueEdit;
	QHBoxLayout *h;
	QVBoxLayout *v;

	double m;
	double m_1;
	double a{0.0};
	double b{1.0};

	bool titleAbove{true};
	bool logarithmic;
	int steps{1000};
};

#endif // SLIDER_H
