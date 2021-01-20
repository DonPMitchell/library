#include "stdafx.h"
#include "Video.h"
#pragma intrinsic(sqrt)

//static unsigned char rgchFrame[3][720][1280];
//static unsigned char rgchHeap[3*720*1280];

int
Video::NewVideo(char *szFileName, Image &im)
{
    char rgcHeader[512];
    unsigned long nBytesWritten;

    if (im.m_nChannels != 3 || im.m_nHeight == 0 || im.m_nWidth == 0 || im.m_nHeight % 4 || im.m_nWidth % 4)      // friendly to 8x8 DCT coding eventually
        return 0;
    m_nHeight = im.m_nHeight;
    m_nWidth  = im.m_nWidth;
    hFile = CreateFile(szFileName, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, 
                                   FILE_ATTRIBUTE_NORMAL, NULL);
    if (hFile == INVALID_HANDLE_VALUE)
        return 0;
    pchBuffer = new unsigned char[3*m_nWidth*m_nHeight];
    rgcHeader[0] = 0;
    //
    // frame dimensions, 30 fps, progressive, square pixels, YCbCr
    //
    sprintf_s(rgcHeader, sizeof(rgcHeader), "YUV4MPEG2 W%d H%d F30:1 Ip A1:1 C444\n", m_nWidth, m_nHeight);  
    printf("Test Video |%s|\n", rgcHeader);
	if (WriteFile(hFile, rgcHeader, (DWORD)strlen(rgcHeader), &nBytesWritten, NULL) == 0)
        return 0;
    return 1;
}

int
Video::Close()
{
    CloseHandle(hFile);
    m_nHeight = m_nWidth = 0;
    delete pchBuffer;
    pchBuffer = 0;
    return 1;
}

int
Video::WriteFrame(Image &im)
{
    int i, j;
    unsigned long nBytesWritten;
    double Y, Cb, Cr;
    DisplayRGB rgb;

	if (WriteFile(hFile, "FRAME\n", (DWORD)strlen("FRAME\n"), &nBytesWritten, NULL) == 0)
        return 0;
    for (j = 0; j < m_nHeight; j++) {
        for (i = 0; i < m_nWidth; i++) {
            rgb = im.GetRGB(i, j);
            if (rgb.red < 0.0 || rgb.red > 1.0 
             || rgb.grn < 0.0 || rgb.grn > 1.0
             || rgb.blu < 0.0 || rgb.blu > 1.0)
                rgb = ClampToGamut(rgb);
            rgb.red = pow(rgb.red, 1.0f/m_fGamma);    // simple gamma 2.0
            rgb.grn = pow(rgb.grn, 1.0f/m_fGamma);
            rgb.blu = pow(rgb.blu, 1.0f/m_fGamma);
            Y  =  16 +  65.481*rgb.red + 128.553*rgb.grn +  24.966*rgb.blu;
            Cb = 128 -  37.797*rgb.red -  74.203*rgb.grn + 112.000*rgb.blu;
            Cr = 128 + 112.000*rgb.red -  93.786*rgb.grn -  18.214*rgb.blu;
            pchBuffer[0*m_nHeight*m_nWidth + j*m_nWidth + i] = int(Y);
            pchBuffer[1*m_nHeight*m_nWidth + j*m_nWidth + i] = int(Cb);
            pchBuffer[2*m_nHeight*m_nWidth + j*m_nWidth + i] = int(Cr);
        }
    }
	if (WriteFile(hFile, pchBuffer, 3*m_nWidth*m_nHeight, &nBytesWritten, NULL) == 0)
        return 0;
    return 1;
}

static double
RandByte()
{
    static unsigned n = 1;
    double f;

    n = 2654435761 + n * 1099087573;
    f= double(n)/4294967296.0;
    return 0.25 + 0.75*f;
}

#define D_PI  3.14159265358979323846264338327950288419716939937510

void
TestVideoClass()
{
    Image im;
    Video vd;
    static DisplayRGB rgColors[360];
    DisplayRGB rgb;
    double fTheta, fRadians, fTurn, xS, yS, xE, yE;
    int i;

    for (i = 0; i < 360; i++)
        rgColors[i] = DisplayRGB(RandByte(), RandByte(), RandByte());
    im.NewImage(1280, 720, 3);
    vd.NewVideo("TestVideo.mpg", im);
    for (fTurn = 0.0; fTurn < 360.0; fTurn += 1.0) {
        im.FillRGB(ML_Black);
        for (fTheta = 0.0; fTheta < 360.0; fTheta += 10.0) {
            fRadians = (fTheta + fTurn) * D_PI/180.0;
            rgb = rgColors[int(fTheta)];
            xS = 640.0 + 25.0*sin(fRadians);
            yS = 360.0 + 25.0*cos(fRadians);
            xE = 640.0 + 300.0*sin(fRadians);
            yE = 360.0 + 300.0*cos(fRadians);
            im.DrawLineRGB(rgb, xS, yS, xE, yE, 1.0);
        }
        vd.WriteFrame(im);
        if (fTurn == 0.0)
            im.WriteBMP("VideoFrame.bmp");
    }
    vd.Close();
}