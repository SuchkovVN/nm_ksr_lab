#include "helpform.h"
#include "ui_helpform.h"
#include <QPixmap>

#ifndef IMAGE_PATH
#define IMAGE_PATH ""
#endif

HelpForm::HelpForm(QWidget *parent) :
    QDialog(parent),
    ui(new Ui::HelpForm)
{
    ui->setupUi(this);
}

HelpForm::~HelpForm()
{
    delete ui;
}
