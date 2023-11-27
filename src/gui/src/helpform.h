#ifndef HELPFORM_H
#define HELPFORM_H

#include <QDialog>

namespace Ui {
class HelpForm;
}

class HelpForm : public QDialog
{
    Q_OBJECT

public:
    explicit HelpForm(QWidget *parent = nullptr);
    ~HelpForm();

private:
    Ui::HelpForm *ui;
};

#endif // HELPFORM_H
