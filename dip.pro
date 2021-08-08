QT += widgets
requires(qtConfig(filedialog))
qtHaveModule(printsupport): QT += printsupport
QMAKE_CXXFLAGS+= -openmp

HEADERS       = mainwindow.h \
                fftw3.h \
                imageprocess.h \
                mdichild.h \
                padding.h \
                transform.h
SOURCES       = main.cpp \
                imagepocess.cpp \
                mainwindow.cpp \
                mdichild.cpp \
                padding.cpp \
                transform.cpp

RESOURCES     = \
    dip.qrc

# install
target.path = $$[QT_INSTALL_EXAMPLES]/widgets/mainwindows/mdi
INSTALLS += target

win32: LIBS += -L$$PWD/./ -llibfftw3-3 -llibfftw3f-3 -llibfftw3l-3

INCLUDEPATH += F:/opencv/build/include


INCLUDEPATH += $$PWD/.
DEPENDPATH += $$PWD/.

TRANSLATIONS += $$PWD/languages/zh_CN.ts \
               $$PWD/languages/en_US.ts


win32:CONFIG(release, debug|release): LIBS += -LF:/opencv/build/x64/vc14/lib/ -lopencv_world440
else:win32:CONFIG(debug, debug|release): LIBS += -LF:/opencv/build/x64/vc14/lib/ -lopencv_world440d

INCLUDEPATH += F:/opencv/build/x64/vc14
DEPENDPATH += F:/opencv/build/x64/vc14
