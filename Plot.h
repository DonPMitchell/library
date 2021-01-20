#pragma once
//
//  plot - basic scientific plotting (similar to my 1982 SRL plotting procedures)
//  D. P. Mitchell  2004/07/04.
//
#include "Image.h"

class ScientificPlot {
public:
                    ScientificPlot() {}

    int             NewPlot(int nWidth, int nHeight);       // size, default viewport and window (0, nWidth-1)
    int             NewPlot(char *szColorImage);            // start with a background image
    int             SetViewPort(int iLowX, int jLowY, int iHighX, int jHighY); //subrectangle of nWidth x nHeight
    int             SetWindow(double xLow, double yLow, double xHigh, double yHigh); // mapping to view port
    int             Move(double x, double y);
    int             Draw(double x, double y);
    double          Print(char *sz);
    double          Plot(int nSymbol);
    DisplayRGB      SetPlotColor(DisplayRGB rgb);
    double          SetLineThickness(double fThick);
    double          SetTextSize(double fFont);
    double          SetTextDirection(double fDir);
    int             SetTextJustification(int nJust);
    int             XAxis(double xOrg, double yOrg, double xLength, int nSegments); // drawn on viewpoint boundary
    int             YAxis(double xOrg, double yOrg, double yLength, int nSegments);
    int             Label(char *szLabelFormat, char *szTitle, int nLabelFont, int nTitleFont);
    int             SetTics(int nTicLength, int nTicType);
    int             Fill(DisplayRGB rgb);
    int             WriteBMP(char *szFileName);

//private:
    DisplayRGB      m_rgbPlotColor;
    Image           m_imPlot;                               // plot image
    int             m_iLowX, m_jLowY, m_iHighX, m_jHighY;   // viewport
    double          m_xLow, m_yLow, m_xHigh, m_yHigh;       // window
    double          m_xA, m_xB, m_yA, m_yB;                 // window -> viewport transform
    double          m_xCursor, m_yCursor;                   // current drawing origin
    double          m_fThick;
    double          m_nTextSize;
    double          m_nTextDirection;
    int             m_nTextJustification;
    int             m_nTicType;
    int             m_nTicLength, m_nTicScale;
    double          m_xOrg, m_yOrg, m_fLength;              // axis info used for labelling
    int             m_bXAxis;
    int             m_nSegments;
};

#define ML_TICS_LEFT    0
#define ML_TICS_UP      0
#define ML_TICS_RIGHT   1
#define ML_TICS_DOWN    1
#define ML_TICS_BOTH    2
#define ML_TICS_IN      3
#define ML_TICS_OUT     4

#define ML_TICS_SHORT   4
#define ML_TICS_MEDIUM  2
#define ML_TICS_LONG    1

#define FONT_14         14.0
