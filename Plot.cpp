#include "stdafx.h"
#include "Image.h"
#include "Plot.h"
#include <stdio.h>
#define D_PI  3.14159265358979323846264338327950288419716939937510

int
ScientificPlot::NewPlot(int nWidth, int nHeight)
{
    if (m_imPlot.NewImage(nWidth, nHeight, 3) == 0)
        return 0;
    m_imPlot.m_bWrapX = m_imPlot.m_bWrapY = WRAP_CLAMP;
    SetViewPort(0, 0, nWidth-1, nHeight-1);
    //
    //  NOTE: y axis usually points upwards, unlike image pixel j axis.
    //
    SetWindow(0.0, 0.0, double(nWidth-1), double(nHeight-1));
    Fill(ML_Black);
    SetPlotColor(ML_White);                     // reset all defaults
    SetLineThickness(1.0);
    SetTextSize(FONT_14);
    SetTextDirection(TEXT_RIGHT);
    SetTextJustification(JUSTIFIED_LEFT);
    SetTics(ML_TICS_MEDIUM, ML_TICS_IN);
    return 1;
}

int
ScientificPlot::NewPlot(char *szColorImage)
{
    Image imTmp;
    int nWidth, nHeight, i, j;
    float f;

    if (m_imPlot.ReadBMP(szColorImage) == 0)
        return 0;
    nWidth = m_imPlot.m_nWidth;
    nHeight = m_imPlot.m_nHeight;
    if (m_imPlot.m_nChannels != 3) {            // convert B/W image into RGB
        m_imPlot.NewImage(nWidth, nHeight, 3);
        imTmp.ReadBMP(szColorImage);
        for (j = 0; j < nHeight; j++) {
            for (i = 0; i < nWidth; i++) {
                f = imTmp.Get(i, j);
                m_imPlot.Set(f, i, j, 0);
                m_imPlot.Set(f, i, j, 1);
                m_imPlot.Set(f, i, j, 2);
            }
        }
    }
    m_imPlot.m_bWrapX = m_imPlot.m_bWrapY = WRAP_CLAMP;
    SetViewPort(0, 0, nWidth-1, nHeight-1);
    SetWindow(0.0, 0.0, double(nWidth-1), double(nHeight-1));
    SetPlotColor(ML_White);                     // reset all defaults
    SetLineThickness(1.0);
    SetTextSize(FONT_14);
    SetTextDirection(TEXT_RIGHT);
    SetTextJustification(JUSTIFIED_LEFT);
    SetTics(ML_TICS_MEDIUM, ML_TICS_IN);
    return 1;
}

int
ScientificPlot::WriteBMP(char *szFileName)
{
    return m_imPlot.WriteBMP(szFileName);
}

int
ScientificPlot::Fill(DisplayRGB rgb)
{
    m_imPlot.FillRGB(rgb);
    return 1;
}

int
ScientificPlot::SetViewPort(int iLowX, int jLowY, int iHighX, int jHighY)
{
    int iTmp;

    if (iLowX > m_imPlot.m_nWidth || iLowX < 0 || iHighX > m_imPlot.m_nWidth || iHighX < 0)
        return 0;
    if (jLowY > m_imPlot.m_nHeight || jLowY < 0 || jHighY > m_imPlot.m_nHeight || jHighY < 0)
        return 0;
    //
    //  Put into canonical form.  Note image viewport has origin at the upper left, but plot
    //  window has origin in lower left.
    //
    if (iLowX > iHighX) {
        iTmp = iLowX;
        iLowX = iHighX;
        iHighX = iTmp;
    }
    if (jHighY > jLowY) {
        iTmp = jLowY;
        jLowY = jHighY;
        jHighY = iTmp;
    }
    m_iLowX = iLowX;
    m_jLowY = jLowY;
    m_iHighX = iHighX;
    m_jHighY = jHighY;
    return 1;
}

int
ScientificPlot::SetWindow(double xLow, double yLow, double xHigh, double yHigh)
{
    if (xLow < xHigh) {     // force into canonical order for easier clipping
        m_xLow = xLow;
        m_xHigh = xHigh;
    } else {
        m_xLow = xHigh;
        m_xHigh = xLow;
    }
    if (yLow < yHigh) {
        m_yLow = yLow;
        m_yHigh = yHigh;
    } else {
        m_yLow = yHigh;
        m_yHigh = yLow;
    }
    m_xA = double(m_iHighX - m_iLowX)/(m_xHigh - m_xLow);   // window -> viewport transform
    m_xB = -m_xLow * m_xA + INDEX_TO_SAMPLE(m_iLowX);
    m_yA = double(m_jHighY - m_jLowY)/(m_yHigh - m_yLow);
    m_yB = -m_yLow * m_yA + INDEX_TO_SAMPLE(m_jLowY);
    //if (m_iHighX - m_iLowX > m_jHighY - m_jLowY)
    //    m_nTicScale = (m_iHighX - m_iLowX)/40;              // long tic size
    //else
        m_nTicScale = (m_jHighY - m_jLowY)/40;
    return 1;
}

DisplayRGB
ScientificPlot::SetPlotColor(DisplayRGB rgb)
{
    DisplayRGB rgbOld;

    rgbOld = m_rgbPlotColor;
    m_rgbPlotColor = rgb;
    return rgbOld;
}

double
ScientificPlot::SetLineThickness(double fThick)
{
    double fOldThick;

    fOldThick = m_fThick;
    m_fThick = fThick;
    return fOldThick;
}

int
ScientificPlot::SetTics(int nTicLength, int nTicType)
{
    m_nTicLength = nTicLength;
    m_nTicType = nTicType;
    return 1;
}

double
ScientificPlot::SetTextSize(double fFont)
{
    double f;

    f = m_nTextSize;
    m_nTextSize = 0.6 * fFont/12.0;       // font size divided by vector font size.  0.6 is a fudge factor to match Caltech plots
    return f * 12.0 / 0.6;
}

double
ScientificPlot::SetTextDirection(double fDir)
{
    double f;

    f = m_nTextDirection;
    m_nTextDirection = fDir;
    return f;
}

int
ScientificPlot::SetTextJustification(int nJust)
{
    int n;

    n = m_nTextJustification;
    m_nTextJustification = nJust;
    return n;
}

int
ScientificPlot::Move(double x, double y)
{
    m_xCursor = m_xA * x + m_xB;    // cursor in viewport coordinates
    m_yCursor = m_yA * y + m_yB;
    return 1;
}

int
ScientificPlot::Draw(double x, double y)
{
    x = m_xA * x + m_xB;
    y = m_yA * y + m_yB;
    m_imPlot.DrawLineRGB(m_rgbPlotColor, m_xCursor, m_yCursor, x, y, m_fThick);
    m_xCursor = x;
    m_yCursor = y;
    return 1;
}

int
ScientificPlot::XAxis(double xOrg, double yOrg, double xLength, int nSegments)
{
    double xStart, xStop, yStartTic, yStopTic, fTicLength, xTic;
    int i, nTicType;

    m_xOrg = xOrg;
    m_yOrg = yOrg;
    m_fLength = xLength;
    m_bXAxis = 1;
    m_nSegments = nSegments;
    nTicType = m_nTicType;
    if (xOrg < m_xLow)              // clip axis to window
        xStart = m_xLow;
    else
        xStart = xOrg;
    if (xOrg + xLength > m_xHigh)
        xStop = m_xHigh;
    else
        xStop = xOrg + xLength;
    Move(xStart, yOrg);
    Draw(xStop, yOrg);
    fTicLength = fabs(double(m_nTicScale)/m_yA/double(m_nTicLength));
    if (nTicType == ML_TICS_IN) {
        if (yOrg < 0.5*(m_yLow + m_yHigh))
            nTicType = ML_TICS_UP;
        else
            nTicType = ML_TICS_DOWN;
    }
    if (nTicType == ML_TICS_OUT) {
        if (yOrg < 0.5*(m_yLow + m_yHigh))
            nTicType = ML_TICS_DOWN;
        else
            nTicType = ML_TICS_UP;
    }
    if (nTicType == ML_TICS_DOWN)
        yStartTic = yOrg;
    else
        yStartTic = yOrg + fTicLength;
    if (nTicType == ML_TICS_UP)
        yStopTic = yOrg;
    else
        yStopTic = yOrg - fTicLength;
    for (i = 1; i < nSegments; i++) {
        xTic = xOrg + (i * xLength)/double(nSegments);
        if (xTic >= m_xLow && xTic <= m_xHigh) {
            Move(xTic, yStartTic);
            Draw(xTic, yStopTic);
        }
    }
    return 1;
}

int
ScientificPlot::YAxis(double xOrg, double yOrg, double yLength, int nSegments)
{
    double yStart, yStop, xStartTic, xStopTic, fTicLength, yTic;
    int i, nTicType;

    m_xOrg = xOrg;
    m_yOrg = yOrg;
    m_fLength = yLength;
    m_bXAxis = 0;
    m_nSegments = nSegments;
    nTicType = m_nTicType;
    if (yOrg < m_yLow)              // clip axis to window
        yStart = m_yLow;
    else
        yStart = yOrg;
    if (yOrg + yLength > m_yHigh)
        yStop = m_yHigh;
    else
        yStop = yOrg + yLength;
    Move(xOrg, yStart);
    Draw(xOrg, yStop);
    fTicLength = fabs(double(m_nTicScale)/m_xA/double(m_nTicLength));
    if (nTicType == ML_TICS_IN) {
        if (xOrg < 0.5*(m_xLow + m_xHigh))
            nTicType = ML_TICS_RIGHT;
        else
            nTicType = ML_TICS_LEFT;
    }
    if (nTicType == ML_TICS_OUT) {
        if (xOrg < 0.5*(m_xLow + m_xHigh))
            nTicType = ML_TICS_LEFT;
        else
            nTicType = ML_TICS_RIGHT;
    }
    if (nTicType == ML_TICS_LEFT)
        xStartTic = xOrg;
    else
        xStartTic = xOrg + fTicLength;
    if (nTicType == ML_TICS_RIGHT)
        xStopTic = xOrg;
    else
        xStopTic = xOrg - fTicLength;
    for (i = 1; i < nSegments; i++) {
        yTic = yOrg + (i * yLength)/double(nSegments);
        if (yTic >= m_yLow && yTic <= m_yHigh) {
            Move(xStartTic, yTic);
            Draw(xStopTic, yTic);
        }
    }
    return 1;
}

double
ScientificPlot::Print(char *sz)
{
    double fWide, dx, dy, y;

    dx = dy = 0.0;
    if (m_nTextDirection == TEXT_RIGHT)         // in plotting, text centered on path, not above it.
        dy -= m_nTextSize*6.0;
    else if (m_nTextDirection == TEXT_UP)
        dx += m_nTextSize*6.0;
    else if (m_nTextDirection == TEXT_DOWN)
        dx -= m_nTextSize*6.0;
    y = double(m_imPlot.m_nHeight) - m_yCursor + dy;
    fWide = m_imPlot.DrawTextRGB(m_rgbPlotColor, m_xCursor + dx, y, sz, m_nTextSize, m_nTextDirection, m_nTextJustification);
    if (fWide == 0)
        return 0;
    if (m_nTextJustification == JUSTIFIED_CENTER)
        fWide = fWide/2.0;                  // print followed by left-justified plots works
    else if (m_nTextJustification == JUSTIFIED_RIGHT)
        fWide = -fWide;                     // prints and plots work if all right-justified
    if (m_nTextDirection == TEXT_RIGHT)
        m_xCursor += fWide;
    else if (m_nTextDirection == TEXT_UP)
        m_yCursor -= fWide;
    else if (m_nTextDirection == TEXT_DOWN)
        m_yCursor += fWide;
    return fWide;
}

double
ScientificPlot::Plot(int nSymbol)
{
    double fWide, y, dx, dy;
    char sz[2];

    sz[0] = nSymbol;
    sz[1] = 0;
    dx = dy = 0.0;
    if (m_nTextDirection == TEXT_RIGHT) {           // in plotting, text centered on path, not above it.
        dy -= m_nTextSize*6.0;
        dx -= m_nTextSize*4.0;
    } else if (m_nTextDirection == TEXT_UP) {
        dx += m_nTextSize*6.0;
        dy -= m_nTextSize*4.0;
    } else if (m_nTextDirection == TEXT_DOWN) {
        dx -= m_nTextSize*6.0;
        dy += m_nTextSize*4.0;
    }
    y = double(m_imPlot.m_nHeight) - m_yCursor + dy;
    fWide = m_imPlot.DrawTextRGB(m_rgbPlotColor, m_xCursor + dx, y, sz, m_nTextSize, m_nTextDirection, m_nTextJustification);
    if (fWide == 0)
        return 0;
    if (m_nTextJustification == JUSTIFIED_CENTER)
        fWide = fWide/2;                // print followed by left-justified plots works
    else if (m_nTextJustification == JUSTIFIED_RIGHT)
        fWide = -fWide;                 // prints and plots work if all right-justified
    if (m_nTextDirection == TEXT_RIGHT)
        m_xCursor += fWide;
    else if (m_nTextDirection == TEXT_UP)
        m_yCursor -= fWide;
    else if (m_nTextDirection == TEXT_DOWN)
        m_yCursor += fWide;
    return fWide;
}

int
ScientificPlot::Label(char *szLabelFormat, char *szTitle, int nLabelFont, int nTitleFont)
{
    int nTextJust, i, j;
    double xLabel, yLabel, xTitle, yTitle, nFont, nTextDir;
    char szLabel[128];

    nTextJust = SetTextJustification(JUSTIFIED_CENTER);
    if (m_bXAxis) {
        //
        //  Horizontal labels and titles (X-Axis)
        //
        nTextDir = SetTextDirection(TEXT_RIGHT);
        xLabel = xTitle = 0.5*(m_xLow + m_xHigh);
        if (m_yOrg > 0.5*(m_yLow + m_yHigh)) {  // top
            if (szLabelFormat == 0) {
                yLabel = yTitle = m_yHigh;
            } else {
                yLabel = m_yHigh + fabs((nLabelFont/2 + 3)/m_yA);
                yTitle = m_yHigh + fabs((nLabelFont + 3)/m_yA);
            }
            if (m_nTicType == ML_TICS_UP || m_nTicType == ML_TICS_OUT || m_nTicType == ML_TICS_BOTH) {
                yLabel += fabs(m_nTicScale/m_nTicLength/m_yA);
                yTitle += fabs(m_nTicScale/m_nTicLength/m_yA);
            }
            if (szTitle)
                yTitle += fabs((3*nTitleFont/4 + 1)/m_yA);
        } else {
            if (szLabelFormat == 0) {
                yLabel = yTitle = m_yLow;
            } else {
                yLabel = m_yLow - fabs((nLabelFont/2 + 3)/m_yA);
                yTitle = m_yLow - fabs((nLabelFont + 3)/m_yA);
            }
            if (m_nTicType == ML_TICS_DOWN || m_nTicType == ML_TICS_OUT || m_nTicType == ML_TICS_BOTH) {
                yLabel -= fabs(m_nTicScale/m_nTicLength/m_yA);
                yTitle -= fabs(m_nTicScale/m_nTicLength/m_yA);
            }
            if (szTitle)
                yTitle -= fabs((3*nTitleFont/4 + 1)/m_yA);
        }
        if (szLabelFormat) {
            nFont = SetTextSize(nLabelFont);
            for (i = 1; i < m_nSegments; i++) {   // nSegment-1 tics are drawn
                xLabel = double(i * m_fLength)/m_nSegments + m_xOrg;
                if (xLabel >= m_xLow && xLabel <= m_xHigh) {
                    if (fabs(xLabel) < fabs(1.0e-10*m_fLength))
                        xLabel = 0.0;
                    //sprintf(szLabel, szLabelFormat, xLabel);
                    sprintf_s(szLabel, sizeof(szLabel), szLabelFormat, xLabel);
                    Move(xLabel, yLabel);
                    Print(szLabel);
                }
            }
            SetTextSize(nFont);
        }
        if (szTitle) {
            nFont = SetTextSize(nTitleFont);
            Move(xTitle, yTitle);
            Print(szTitle);
            SetTextSize(nFont);
        }
    } else {
        //
        //  Vertical labels and titles (Y-Axis)
        //
        yLabel = yTitle = 0.5*(m_yLow + m_yHigh);
        if (m_xOrg < 0.5*(m_xLow + m_xHigh)) {      // left
            nTextDir = SetTextDirection(TEXT_UP);
            if (szLabelFormat == 0) {
                xLabel = m_xLow;
            } else {
                xLabel = m_xLow - fabs((nLabelFont/2 + 3)/m_xA);
                xTitle = m_xLow - fabs((nLabelFont + 3)/m_xA);
            }
            if (m_nTicType == ML_TICS_LEFT || m_nTicType == ML_TICS_OUT || m_nTicType == ML_TICS_BOTH) {
                xLabel -= fabs(m_nTicScale/m_nTicLength/m_xA);
                xTitle -= fabs(m_nTicScale/m_nTicLength/m_xA);
            }
            if (szTitle)
                xTitle -= fabs((3*nTitleFont/4 + 1)/m_xA);
        } else {
            nTextDir = SetTextDirection(TEXT_DOWN);
            if (szLabelFormat == 0) {
                xLabel = m_xHigh;
            } else {
                xLabel = m_xHigh + fabs((nLabelFont/2 + 3)/m_xA);
                xTitle = m_xHigh + fabs((nLabelFont + 3)/m_xA);
            }
            if (m_nTicType == ML_TICS_RIGHT || m_nTicType == ML_TICS_OUT || m_nTicType == ML_TICS_BOTH) {
                xLabel += fabs(m_nTicScale/m_nTicLength/m_xA);
                xTitle += fabs(m_nTicScale/m_nTicLength/m_xA);
            }
            if (szTitle)
                xTitle += fabs((3*nTitleFont/4 + 1)/m_xA);
        }
        if (szLabelFormat) {
            nFont = SetTextSize(nLabelFont);
            for (j = 1; j < m_nSegments; j++) {   // nSegment-1 tics are drawn
                yLabel = (j * m_fLength)/m_nSegments + m_yOrg;
                if (fabs(yLabel) < fabs(1.0e-10*m_fLength))
                    yLabel = 0.0;
                sprintf_s(szLabel, sizeof(szLabel), szLabelFormat, yLabel);
                Move(xLabel, yLabel);
                Print(szLabel);
            }
            SetTextSize(nFont);
        }
        if (szTitle) {
            nFont = SetTextSize(nTitleFont);
            Move(xTitle, yTitle);
            Print(szTitle);
            SetTextSize(nFont);
        }
    }
    SetTextJustification(nTextJust);
    SetTextDirection(nTextDir);
    return 1;
}

static double
SimpleRand()
{
    static unsigned x = 1234567;

    x = 1099087573*x + 2654435761;
    return double(x) / 4294967296.0;
}

static void
Septagon()
{
    int i;
    double theta, x, y;

    for (i = 0; i <= 7; i++) {
        theta = double(i) * 360.0 / 7.0;
        theta = D_PI * theta / 180.0;
        x = floor(4.0 + 4.0*sin(theta) + 0.5);
        y = floor(6.0 + 4.0*cos(theta) + 0.5);
        printf("P(%d,%d), ", int(x), int(y));
    }
    printf("\n");
    for (i = 0; i <= 7; i++) {
        theta = double(i) * 360.0 / 7.0;
        theta = D_PI * theta / 180.0;
        x = floor(2.0 + 2.0*sin(theta) + 0.5);
        y = floor(10.0 + 2.0*cos(theta) + 0.5);
        printf("P(%d,%d), ", int(x), int(y));
    }
    printf("\n");
}

void
TestPlot()
{
    ScientificPlot pl;
    Image im;
    int i;
    double x, y, width, fSize;
    DisplayRGB rgb;
    static unsigned char sz[8] = "SYM: ";

    Septagon();

    printf("Testing ML_Plot class\n");
    pl.NewPlot(640, 480);
    pl.WriteBMP("Plot_Blank.bmp");

    im.NewImage(640, 640, 3);
    width = im.DrawTextRGB(ML_White, 100.0, 100.0, "Hello ", 1.5);
    im.DrawTextRGB(ML_Red, 100.0+width, 100.0, "World", 1.5);
    im.WriteBMP("TestImageText.bmp");

    pl.NewPlot(640, 400);
    pl.SetPlotColor(ML_White);
    pl.SetViewPort(40, 40, 600, 360);
    pl.SetWindow(0.0, 0.0, 2.0, 1.0);
    pl.XAxis(0.0, 0.0, 2.0, 4);
    pl.XAxis(0.0, 1.0, 2.0, 8);
    pl.Label(0, "Test of Plotting Routines, 5", 18, 24);
    fSize = pl.SetTextSize(24);
    pl.Plot('~' + 1);
    pl.Plot('~' + 2);
    pl.Plot('~' + 3);
    pl.Plot('~' + 4);
    pl.Plot('~' + 5);
    pl.Plot('~' + 6);
    pl.Move(1.0, 0.5);
    pl.Plot(ML_SYM_SQUARE);
    pl.SetTextSize(fSize);
    //pl.Plot(ML_GREEK_LOW+7);
    //pl.Plot('=');
    //pl.Plot('7');
    //pl.Plot(ML_SYM_DEG);
    pl.XAxis(0.0, 0.0, 2.0, 8);
    pl.Label("%g", 0, 18, 24);
    pl.YAxis(0.0, 0.0, 1.0, 5);
    pl.Label("%0.3f", "Axis Y", 14, 18);
    //pl.Plot(ML_SYM_SUP3);
    pl.YAxis(2.0, 0.0, 1.0, 3);
    pl.Label("%g", "The Right Side", 14, 18);
    pl.Move(0.5, 0.5);
    pl.Print("Now is the ");
    pl.Print("time for all ");
    //pl.Plot(ML_GREEK_LOW+7);
    //pl.Plot('=');
    //pl.Plot('7');
    pl.Move(1.1, 0.2);
    pl.Draw(1.6, 0.6);
    pl.Draw(1.2, 0.8);
    pl.WriteBMP("Plot_Test.bmp");

    pl.NewPlot(640, 396);
    pl.SetViewPort(40, 40, 600, 396-40);
    pl.SetWindow(0.0, -0.5, 20.0, 1.0);
    pl.XAxis(0.0, -0.5, 20.0, 10);
    pl.Label("%g", 0, 14, 14);
    pl.XAxis(0.0, 1.0, 20.0, 10);
    pl.Label(0, "Bessel Functions", 14, 24);
    pl.YAxis(0.0, -0.5, 1.5, 6);
    pl.Label("%g", 0, 14, 14);
    pl.YAxis(20.0, -0.5, 1.5, 6);
    pl.Label("%g", 0, 14, 14);

    for (i = 0; i < 5; i++) {
        rgb = DisplayRGB(SimpleRand(), SimpleRand(), SimpleRand());
        rgb = rgb/rgb.MaxRGB();
        pl.SetPlotColor(rgb);
        for (x = 0.0; x < 20.0; x += 0.1) {
            y = _jn(i, x);
            if (x == 0.0)
                pl.Move(x, y);
            pl.Draw(x, y);
        }
    }
    pl.Move(10.0, 0.9);
    pl.SetPlotColor(DisplayRGB(1.0, 1.0, 0.3));
    //pl.Print("Bessel Functions");
    pl.WriteBMP("Plot_Bessel.bmp");

}