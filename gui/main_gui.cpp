#include <QApplication>
#include <QMainWindow>
#include <QWidget>
#include <QVBoxLayout>
#include <QFormLayout>
#include <QLineEdit>
#include <QPushButton>
#include <QTextEdit>
#include <QFile>
#include <QTextStream>
#include <QMessageBox>
#include "../APlug.h"

class MainWindow : public QMainWindow {
public:
    MainWindow() {
        QWidget *central = new QWidget;
        QVBoxLayout *layout = new QVBoxLayout;
        QFormLayout *form = new QFormLayout;

        tasksEdit = new QLineEdit("0");
        sampleSizeEdit = new QLineEdit("0");
        sampleFileEdit = new QLineEdit;
        minWorkersEdit = new QLineEdit("1");
        maxWorkersEdit = new QLineEdit(QString::number(MAX_NUM_WORKERS));
        stepEdit = new QLineEdit("15");

        form->addRow("Number of tasks", tasksEdit);
        form->addRow("Sample size", sampleSizeEdit);
        form->addRow("Sample file path", sampleFileEdit);
        form->addRow("Minimum workers", minWorkersEdit);
        form->addRow("Maximum workers", maxWorkersEdit);
        form->addRow("Step size", stepEdit);
        layout->addLayout(form);

        QPushButton *runBtn = new QPushButton("Run");
        layout->addWidget(runBtn);

        outputEdit = new QTextEdit;
        outputEdit->setReadOnly(true);
        layout->addWidget(outputEdit);

        central->setLayout(layout);
        setCentralWidget(central);
        setWindowTitle("APlug GUI");

        connect(runBtn, &QPushButton::clicked, this, [this]() { runAPlug(); });
    }

private:
    void runAPlug() {
        unsigned int numTasks = tasksEdit->text().toUInt();
        unsigned int sampleSize = sampleSizeEdit->text().toUInt();
        QString sampleFile = sampleFileEdit->text();
        int minWorkers = minWorkersEdit->text().toInt();
        int maxWorkers = maxWorkersEdit->text().toInt();
        int step = stepEdit->text().toInt();

        if (sampleFile.isEmpty()) {
            QMessageBox::warning(this, "Input error", "Sample file path is empty");
            return;
        }

        QFile file(sampleFile);
        if (!file.open(QIODevice::ReadOnly | QIODevice::Text)) {
            QMessageBox::warning(this, "File error", "Cannot open sample file");
            return;
        }

        QTextStream in(&file);
        APlug plug;
        while (!in.atEnd()) {
            QString line = in.readLine();
            plug.pushTime(line.toDouble());
        }
        file.close();

        if (sampleSize == 0)
            sampleSize = plug.getSampleSize();
        if (numTasks == 0 || numTasks < sampleSize)
            numTasks = sampleSize;

        plug.setSampleSizeAndNumOfTasks(sampleSize, numTasks);
        plug.getRecommendation(minWorkers, maxWorkers, step);

        QString output;
        output += QString("Number of tasks:  %1\n").arg(plug.getNumOfTasks());
        output += QString("Sample size:      %1\n").arg(plug.getSampleSize());
        output += QString("Sample sum:       %1\n").arg(plug.getSum());
        output += QString("Sample min:       %1\n").arg(plug.getMin());
        output += QString("Sample max:       %1\n").arg(plug.getMax());
        output += QString("Sample mean:      %1\n").arg(plug.getMean());
        output += QString("Sample stdev:     %1\n").arg(plug.getStdev());
        output += QString("Estimated serial time (histogram): %1 sec\n").arg(plug.getEstSerialTimeSecHistogram());
        output += QString("Estimated serial time (gamma): %1 sec\n").arg(plug.getEstSerialTimeSecGamma());
        output += QString("Suggested max workers: %1\n").arg(plug.getMaxWorkers());
        output += QString("Recommended number of workers: %1\n").arg(plug.getRecommendedNumOfWorkers());

        outputEdit->setPlainText(output);
    }

    QLineEdit *tasksEdit;
    QLineEdit *sampleSizeEdit;
    QLineEdit *sampleFileEdit;
    QLineEdit *minWorkersEdit;
    QLineEdit *maxWorkersEdit;
    QLineEdit *stepEdit;
    QTextEdit *outputEdit;
};

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);
    MainWindow w;
    w.show();
    return app.exec();
}

