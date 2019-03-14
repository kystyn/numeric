TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        main.cpp \
    matr.cpp \
    relax.cpp \
    simple.cpp \
    ldr.cpp

HEADERS += \
    distribution.h \
    controller.h \
    interpolator.h \
    matr.h \
    relax.h \
    def.h \
    tabulated_basis_func.h \
    simple.h \
    ldr.h
