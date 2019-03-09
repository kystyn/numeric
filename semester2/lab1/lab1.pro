TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp \
    distribution.cpp

HEADERS += \
    controller.h \
    distribution.h \
    interpolator.h
