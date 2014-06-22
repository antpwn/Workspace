TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

# include SDL
LIBS += -L/usr/lib -lSDL2

INCLUDEPATH = usr/include

SOURCES += main.cpp

