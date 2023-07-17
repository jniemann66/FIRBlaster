#include "coefficientswidget.h"

#include <QVBoxLayout>

CoefficientsWidget::CoefficientsWidget(QWidget *parent)
	: QWidget{parent}
{
	coeffs = new QTextEdit;
	coefficientCountLabel = new QLabel;
	coeffs->setReadOnly(true);
	auto mainLayout = new QVBoxLayout;
	mainLayout->addWidget(coefficientCountLabel);
	mainLayout->addWidget(coeffs);
	setLayout(mainLayout);
}

QString CoefficientsWidget::getCoefficients() const
{
	return coeffs->toPlainText();
}

void CoefficientsWidget::setCoefficients(const std::vector<double> &newCoefficients)
{
	coefficientCountLabel->setText(QStringLiteral("<h1>%1</h1>").arg(newCoefficients.size()));

	QString coeffsText;

	QStringList t;
	for(double v : newCoefficients) {
		t.append(QStringLiteral("\t%1")
				 .arg(QString::number(v, 'g', 16))
				 );
	}

	coeffsText = QStringLiteral("const double coefficients[] {\n%1\n};")
				 .arg(t.join(",\n"));

	coeffs->setText(coeffsText);

}
