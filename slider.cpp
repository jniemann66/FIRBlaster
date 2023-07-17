#include "slider.h"

#include <QDoubleValidator>

Slider::Slider(Qt::Orientation orientation, const QString& title, QWidget *parent, bool titleAbove)
	: QWidget(parent), titleAbove(titleAbove)
{
	label = new QLabel(title);
	m_slider = new QSlider(orientation, this);
	valueEdit = new QLineEdit(this);
	auto dv = new QDoubleValidator(this);
	valueEdit->setValidator(dv);

	h = new QHBoxLayout;
	v = new QVBoxLayout;

	setOrientation(orientation);

	connect(m_slider, &QSlider::valueChanged, this, [this](int i){
		double v = value();
		valueEdit->setText(QString::number(v, 'f', 4));
		emit valueChanged(v);
	});

	connect(valueEdit, &QLineEdit::editingFinished, this, [this](){
		QString s = valueEdit->text();
		bool ok = false;
		double v = s.toDouble(&ok);
		if(ok) {
			setValue(v);
		}
	});
}

bool Slider::getLogarithmic() const
{
	return logarithmic;
}

void Slider::setLogarithmic(bool newLogarithmic)
{
	logarithmic = newLogarithmic;
}

void Slider::setOrientation(Qt::Orientation orientation)
{
	if(orientation == Qt::Horizontal) {
		if(titleAbove) {
			v->addWidget(label);
			h->addWidget(m_slider, 10);
			h->addWidget(valueEdit, 1);
			v->addLayout(h);
		} else {
			h->addWidget(label, 2);
			h->addWidget(m_slider, 10);
			h->addWidget(valueEdit, 1);
			v->addLayout(h);
		}
		setLayout(v);
	} else {
		v->addWidget(label, 2);
		v->addWidget(m_slider, 10);
		v->addWidget(valueEdit, 1);
		h->addStretch(2);
		h->addLayout(v,1);
		h->addStretch(2);
		setLayout(h);
	}
}

void Slider::setRange(double min, double max)
{
	a = min;
	b = max;
	m = (b - a) / steps;
	m_1 = 1.0 / m;

	m_slider->setRange(1, steps);
}

void Slider::setSteps(int numSteps)
{
	steps = numSteps;
	m = (b - a) / steps;
	m_1 = 1.0 / m;
}

void Slider::setValue(double v)
{
	m_slider->setValue(m_1 * (v -a));
}

double Slider::value() const
{
	return m_slider->value() * m + a;
}

bool Slider::getTitleAbove() const
{
	return titleAbove;
}

void Slider::setTitleAbove(bool newTitleAbove)
{
	titleAbove = newTitleAbove;
}
