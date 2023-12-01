// stoffdaten.h : Definition der Datenstruktur
// Autor: majo (majo02@avt.rwth-aachen.de), thbe (thbe02@avt.rwth-aachen.de)
// UML-Diagramm:
// -
// Quellen:
// -

#ifndef STOFFDATEN
#define STOFFDATEN

typedef struct tpropdata{ //fasst grundlegende Informationen des Gemisches zusammen
   void *libdata; //Grundlage für tplibdata
   int  nstoff; //gibt die Anzahl der Stoffe an
   char **stname; //gibt die Stoffbezeichnungen an
   int refined_enrtl; //sysmetrisches(true?) oder unsymetrisches(false?) Modell
 } tpropdata;

#endif //STOFFDATEN