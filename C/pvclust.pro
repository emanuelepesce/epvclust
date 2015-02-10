TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

LIBS += -lm

SOURCES +=

SOURCES += pvclust.c
SOURCES += pvclust.h

SOURCES += cluster.c
SOURCES += cluster.h

SOURCES += driverCluster.c
SOURCES += driverCluster.h

SOURCES += libwwd.c
SOURCES += libwwd.h

SOURCES += libsort.c
SOURCES += libsort.h

SOURCES += libClust.c
SOURCES += libClust.h


