TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt
CONFIG += c++11

SOURCES += main.cpp \
    Game.cpp

unix|win32: LIBS += -lSDL2

HEADERS += \
    Game.h
