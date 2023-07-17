#ifndef COEFFICIENTSWIDGET_H
#define COEFFICIENTSWIDGET_H

#include <QWidget>
#include <QTextEdit>
#include <QLabel>

class CoefficientsWidget : public QWidget
{
	Q_OBJECT
public:
	explicit CoefficientsWidget(QWidget *parent = nullptr);
	QString getCoefficients() const;

public slots:
	void setCoefficients(const std::vector<double> &newCoefficients);

signals:

private:
	QLabel *coefficientCountLabel{nullptr};
	QTextEdit *coeffs;


};

#endif // COEFFICIENTSWIDGET_H
