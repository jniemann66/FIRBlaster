#include "filterdesignwidget.h"
#include "firgenerator.h"

#include <algorithm>

#include <QHBoxLayout>
#include <QVBoxLayout>
#include <QString>
#include <QStringLiteral>
#include <QDebug>
#include <QPainter>
#include <QGroupBox>

#define FILTERSIZE_LIMIT 131071
#define FILTERSIZE_BASE 103

FilterDesignWidget::FilterDesignWidget(QWidget *parent)
	: QWidget{parent}
{
//	fftPlot = new QLabel;
//	impulsePlot = new QLabel;
	generateButton = new QPushButton("Generate");
	filterType = new QComboBox;
	startFrequency = new Slider(Qt::Horizontal, "Start Frequency");
	transitionBand = new Slider(Qt::Horizontal, "Transition Band Width");
	dynamicRange = new Slider(Qt::Vertical, "Sidelobe Attenuation");
	minPhaseCheckbox = new QCheckBox("Minimum Phase");	

//	fftPlot->installEventFilter(new PixmapResizer(fftPlotpixmap, this));
//	impulsePlot->installEventFilter(new PixmapResizer(impulsePlotPixmap, this));

	//installEventFilter


	filterType->addItem("Low Pass", static_cast<int>(FilterType::LowPass));
	filterType->addItem("High Pass", static_cast<int>(FilterType::HighPass));
	filterType->setCurrentIndex(0);

	startFrequency->setSteps(1000);
	startFrequency->setRange(0, 1.0);
	startFrequency->setValue(0.95);

	transitionBand->setSteps(1000);
	transitionBand->setRange(0, 1.0);
	transitionBand->setValue(0.05);

	dynamicRange->setSteps(1000);
	dynamicRange->setRange(-200, 0);
	dynamicRange->setValue(-72);

	auto filterTypeLabel = new QLabel("Filter Type");
	auto filterTypeLayout = new QHBoxLayout;
	auto startFrequencyLayout = new QHBoxLayout;
	auto transitionBandLayout = new QHBoxLayout;
	auto dynamicRangeLayout = new QVBoxLayout;
	auto specsLeftLayout = new QVBoxLayout;
	auto calcLayout = new QHBoxLayout;
	auto specsLayout = new QHBoxLayout;
	auto specsBox = new QGroupBox("Specifications");
	auto bottomLayout = new QHBoxLayout;
	auto mainlayout = new QVBoxLayout;

	specsBox->setLayout(specsLayout);

	filterTypeLayout->addWidget(filterTypeLabel);
	filterTypeLayout->addWidget(filterType);
	filterTypeLayout->addWidget(minPhaseCheckbox);
	filterTypeLayout->addStretch();

	startFrequencyLayout->addWidget(startFrequency);
	transitionBandLayout->addWidget(transitionBand);

	dynamicRangeLayout->addWidget(dynamicRange);

	calcLayout->addWidget(generateButton, 1);
	calcLayout->addStretch(3);

	specsLeftLayout->addLayout(filterTypeLayout);
	specsLeftLayout->addLayout(startFrequencyLayout);
	specsLeftLayout->addLayout(transitionBandLayout);
	specsLeftLayout->addStretch();
	specsLeftLayout->addLayout(calcLayout);

	specsLayout->addLayout(specsLeftLayout, 8);
	specsLayout->addLayout(dynamicRangeLayout, 1);

	bottomLayout->addWidget(specsBox);

	mainlayout->addLayout(bottomLayout);

	setLayout(mainlayout);

	connect(filterType, QOverload<int>::of(&QComboBox::currentIndexChanged), this, &FilterDesignWidget::makeFilter);
	connect(startFrequency, &Slider::valueChanged, this,  &FilterDesignWidget::makeFilter);
	connect(transitionBand, &Slider::valueChanged, this,  &FilterDesignWidget::makeFilter);
	connect(dynamicRange, &Slider::valueChanged, this,  &FilterDesignWidget::makeFilter);
	connect(generateButton, &QPushButton::pressed, this, &FilterDesignWidget::makeFilter);
	connect(minPhaseCheckbox, &QCheckBox::toggled, this, &FilterDesignWidget::makeFilter);

	QMetaObject::invokeMethod(this, &FilterDesignWidget::makeFilter, Qt::QueuedConnection);
}

void FilterDesignWidget::makeFilter()
{
	setCursor(QCursor(Qt::BusyCursor));
	auto results = makeFilterCoefficients<double>();
	plotResponse(QVector<double>{results.cbegin(), results.cend()});
	emit createdCoefficients(results);
	unsetCursor();
}

void FilterDesignWidget::setFftPlot(QLabel *newFftPlot)
{
	fftPlot = newFftPlot;
}

void FilterDesignWidget::setImpulsePlot(QLabel *newImpulsePlot)
{
	impulsePlot = newImpulsePlot;
}

template<typename FloatType>
std::vector<FloatType> FilterDesignWidget::makeFilterCoefficients() {

	// determine cutoff frequency and steepness
	double ft = 0.5 * startFrequency->value();
	double steepness = 0.090909091 / (0.5 * transitionBand->value());

	// determine filtersize
	int filterSize = static_cast<int>(
		std::min<int>(FILTERSIZE_BASE *  steepness, FILTERSIZE_LIMIT)
		| 1 // ensure that filter length is always odd
	);

	const bool minPhase = minPhaseCheckbox->isChecked();

	FloatType* pFilterTaps;
	std::vector<FloatType> filterTaps;

	static constexpr bool padIR = false;
	if constexpr(padIR) {
		filterTaps.resize(2 * filterSize, 0.0);
		pFilterTaps = filterTaps.data() + filterSize / 2;
	} else {
		filterTaps.resize(filterSize, 0.0);
		pFilterTaps = filterTaps.data();
	}

	FilterType f = static_cast<FilterType>(filterType->currentData().toInt());
	switch(f) {
	case LowPass:
		makeLPF<FloatType>(pFilterTaps, filterSize, ft);
		break;
	case HighPass:
		makeHPF<FloatType>(pFilterTaps, filterSize, ft);
		break;
	}

	kaiserBeta = calcKaiserBeta(std::abs<double>(dynamicRange->value()));
	dynamicRange->setToolTip(QStringLiteral("Kaiser β=%1").arg(kaiserBeta,0, 'f', 2));
	applyKaiserWindow2<FloatType>(pFilterTaps, filterSize, kaiserBeta);

	// conditionally convert filter coefficients to minimum-phase:
	if (minPhaseCheckbox->isChecked()) {
		makeMinPhase<FloatType>(pFilterTaps, filterSize);
	}

	return filterTaps;
}

void FilterDesignWidget::plotVector(QPixmap* pm, const std::vector<double>& data, const PlotInfo& plotInfo)
{
	QPainter p(pm);

	double w = fftPlotpixmap.width();
	double h = fftPlotpixmap.height();

	std::pair<double, double> range
			= plotInfo.range.value_or(std::make_pair(
										  1.1 * *min_element(data.begin(), data.end()),
										  1.1 * *max_element(data.begin(), data.end())
										  ));

	const auto& [a,b] = range;
	const double sy = h / (b - a);

	double penThickness = 0.0; // 0.0 = cosmetic pen
	QPen plotPen{plotInfo.plotColor, penThickness, Qt::SolidLine, Qt::RoundCap,};
	QPen gridPen{plotInfo.gridColor, penThickness};
	QPen txtPen{plotInfo.txtColor};
	QBrush bg{plotInfo.backgroundColor};

	p.setBrush(bg);
	p.setBackground(bg);
	p.fillRect(fftPlotpixmap.rect(), bg);
	p.setCompositionMode(QPainter::CompositionMode_SourceOver);
	p.setPen(plotPen);
	p.setRenderHint(QPainter::Antialiasing);

	size_t sz = plotInfo.size.value_or(data.size());

	p.translate(0, h + a * sy);
	p.scale(1.0, -sy);

	double xs = 1.0 * w / sz;
	QVector<QPointF> points;
	for(size_t i = 0; i < sz; i++) {
		double x = i * xs;
		QPointF pt{x,  data.at(i)};
		points.append(pt);
	}

	if(plotInfo.connectDots) {
		p.drawPolyline(points);
	} else {
		p.drawPoints(points);
	}

	double o = std::pow(10.0, std::floor(log10(b-a)) - 1.0);

	if(plotInfo.gridlines.first) {
		// horizontal gridlines
		double gy = a - std::fmod(a, o);
		while (gy < b) {
			p.setPen(gridPen);
			p.drawLine(0, gy, w, gy);

			QPointF t = p.transform().map(QPointF{10, gy});
			p.setWorldMatrixEnabled(false);
			p.setPen(txtPen);
			p.drawText(t, QString::number(gy, 'f', plotInfo.precision.first));
			p.setWorldMatrixEnabled(true);
			gy += o;
		}
	}

	if(plotInfo.gridlines.second) {
		// vertical gridlines
		double w_10 = w / 10.0;

		const double tyPos = 0.5 * o;

		for(int r = 0; r <= 10; r++) {
			p.setPen(gridPen);
			double rw_10 = r * w_10;
			p.drawLine(rw_10, a, rw_10, b);

			QPointF t = p.transform().map(QPointF{rw_10, tyPos});
			p.setWorldMatrixEnabled(false);
			p.setPen(txtPen);
			p.drawText(t, QString::number(0.1 * r, 'f', plotInfo.precision.second));
			p.setWorldMatrixEnabled(true);
		}
	}
}

void FilterDesignWidget::plotResponse(const QVector<double>& impulse)
{
	if(fftPlot == nullptr && impulsePlot == nullptr)
		return;

	// plot impulse response
	if(impulsePlot != nullptr) {

		PlotInfo plotInfo;
		plotInfo.range = {-0.5, 1.1};
		plotInfo.gridlines = {true, false};
		plotInfo.precision = {2, 1};
		plotInfo.plotColor = QColor{0, 0, 128, 255};
		plotVector(&impulsePlotPixmap, std::vector<double>{impulse.begin(), impulse.end()}, plotInfo);
		impulsePlot->setPixmap(impulsePlotPixmap);
	}

	// plot fft
	if(fftPlot != nullptr) {
		size_t length = impulse.length();

		static constexpr size_t minFFTsize = 1024;

		auto pow2length =  std::max(minFFTsize, static_cast<size_t>(pow(2, 2 + floor(log2(length)))));
		alignas(16) std::vector <std::complex<double>> complexInput(pow2length, {0.0, 0.0});

		//qDebug() << "FFT size" << pow2length;

		size_t a = (pow2length - length + 1) / 2;

		for (int n = 0; n < length; n++) {
			complexInput[n + a] = {impulse[n], 0.0};
		}

		alignas(16) const std::vector <std::complex<double>> complexOutput = fftV(complexInput);

		std::vector<double> mag2Response;

		for (auto m  = 0; m < pow2length; m++) {
			const std::complex<double> z = complexOutput.at(m);
			mag2Response.push_back(z.real() * z.real() + z.imag() * z.imag());
		}

		convertToDb(mag2Response, /* magSquared= */ true);
		PlotInfo plotInfo;
		plotInfo.range = {-190, 10};
		plotInfo.precision.first = 0;
		plotInfo.size = mag2Response.size() / 2;

//		convertToLinear(mag2Response, /* magSquared= */ true);
//		PlotInfo plotInfo;
//		plotInfo.range = {-0.01, 1.01};
//		plotInfo.precision.first = 1;
//		plotInfo.size = mag2Response.size() / 2;

		plotVector(&fftPlotpixmap, mag2Response, plotInfo);
		fftPlot->setPixmap(fftPlotpixmap);
	}
}

bool FilterDesignWidget::convertToLinear(std::vector<double> &s, bool fromMagSquared)
{
	bool hasSignal{false};
	int numBins = s.size();

	// find peak
	double peak{0.0};
	if(fromMagSquared) {
		for(int b = 0; b < numBins; b++) {
			peak = std::max(peak, std::sqrt (s[b]));
		}
	} else {
		for(int b = 0; b < numBins; b++) {
			peak = std::max(peak, s[b]);
		}
	}

	if(std::fpclassify(peak) != FP_ZERO) {
		hasSignal = true;
		const double scale = 1.0 / peak;

		// function to convert to percentage of fullScale
		if(fromMagSquared) {
			auto scaleFunc = [scale] (double v) -> double {
				return scale * std::sqrt(v);
			};

			// scale the data
			std::transform (s.begin(), s.end(), s.begin(), scaleFunc);

		} else {
			auto scaleFunc = [scale] (double v) -> double {
				return scale * v;
			};

			// scale the data
			std::transform (s.begin(), s.end(), s.begin(), scaleFunc);
		}
	}

	return hasSignal;
}

bool FilterDesignWidget::convertToDb(std::vector<double> &s, bool fromMagSquared)
{
	bool hasSignal{false};
	int numBins = s.size();

	// find peak
	double peak{0.0};

	for(int b = 0; b < numBins; b++) {
		peak = std::max(peak, s[b]);
	}


	if(std::fpclassify(peak) != FP_ZERO) { // scale the data

		hasSignal = true;
		double dBMult = fromMagSquared ? 10.0 : 20.0;

		// set a floor to avoid log(0) problems
		double floor = std::max(std::numeric_limits<double>::min(), peak * pow(10.0, -300.0 / dBMult)); // 300dB below peak or smallest normal number

		// function to convert to dB
		auto scaleFunc = [scale = 1.0 / peak, dBMult, floor] (double v) -> double {
			return dBMult * std::log10(std::max(scale * v, floor));
		};


		std::transform (s.begin(), s.end(), s.begin(), scaleFunc);

	}

	return hasSignal;
}
