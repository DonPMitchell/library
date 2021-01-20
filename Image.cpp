#include "stdafx.h"
#include "Image.h"
//#include "Montecarlo.h"
#pragma intrinsic(log, exp, sqrt, atan2, cos, sin, fabs)
#define D_PI  3.14159265358979323846264338327950288419716939937510

int
Image::NewImage(int nWidth, int nHeight, int nChannels, int bWrapX, int bWrapY)
{
    Index64 nSize, nOldSize;

    nSize = nWidth*Index64(nHeight*nChannels);
    nOldSize = m_nWidth*Index64(m_nHeight*m_nChannels);
    if (m_pfImage == 0 || nSize != nOldSize) {
        delete [] m_pfImage;
        m_pfImage = new float[nSize];
    }
    m_nWidth = nWidth;
    m_nHeight = nHeight;
    m_nChannels = nChannels;
    m_bWrapX = bWrapX;
    m_bWrapY = bWrapY;
    return m_pfImage != 0;
}

float
Image::Fill(double f, int kChan)
{
    Index64 i, iStart, iEnd, nArea;


    if (kChan >= m_nChannels)
        return float(f);
    nArea = m_nHeight*Index64(m_nWidth);
    iStart = kChan*nArea;
    iEnd = iStart + nArea;
    for (i = iStart; i < iEnd; i++)
        m_pfImage[i] = float(f);
    return float(f);
}

DisplayRGB
Image::FillRGB(const DisplayRGB &rgb)
{
    Fill(rgb.red, 0);
    Fill(rgb.grn, 1);
    Fill(rgb.blu, 2);
    return rgb;
}

inline int
ML_Modulus(int n, int m)
{
    n = n % m;
    if (n < 0)
        n += m;
    return n;
}

inline int
ML_Wrap(int &i, int nLength, int nWrap)
{
    int bDraw;

    bDraw = 1;
    if (nWrap == WRAP_AROUND) {
        i = ML_Modulus(i, nLength);
    } else if (nWrap == WRAP_MIRROR) {   //REVIEW: seems to fail on large offsets
        i = ML_Modulus(i, 2*nLength);
        if (i >= nLength)
            i = 2*nLength - i - 1;
    } else {
        bDraw = 0;                      // don't draw when clipped
        if (i < 0)
            i = 0;
        else if (i >= nLength)
            i = nLength - 1;
        else
            bDraw = 1;
    }
    //if (i < 0 || i >= nLength)
    //    printf("ML_Wrap error, i == %d\n", i);
    return bDraw;
}

int
Image::Wrap(int &iCol, int &jRow) const
{
    return ML_Wrap(iCol, m_nWidth, m_bWrapX) & ML_Wrap(jRow, m_nHeight, m_bWrapY);      //beware not to use &&, early out
}

//
//  I/O for version 3 windows bmp image format.
//  OS/2 format not supported.  RLE not supported.
//
int
Image::ReadBMP(char *szFileName, double fGamma)
{
    HANDLE hFile;
	unsigned char	bfType[2];
	struct {                            // Version 3 Header
		unsigned int	bfSize;
		unsigned short	bfReserved1;
		unsigned short	bfReserved2;
		unsigned int	bfOffBits;

		unsigned int	biSize;
		unsigned int	biWidth;
		unsigned int	biHeight;
		unsigned short	biPlanes;
		unsigned short	biBitCount;
		unsigned int	biCompression;
		unsigned int	biSizeImage;
		unsigned int	biXPelsPerMeter;
		unsigned int	biYPelsPerMeter;
		unsigned int	biClrUsed;
		unsigned int	biClrImportant;
	} H;
    struct {
        unsigned nRed;
        unsigned nGrn;
        unsigned nBlu;
    } rgMasks;
    struct {
        unsigned char   b, g, r, a;
    } rgColors[256];
	int i, j, nLineSize, nOffset, nColors, iRed, iGrn, iBlu, nRedShift, nGrnShift, nBluShift;
    int nWidth, nHeight, nChannels;
    unsigned nPixel;
    unsigned long nBytesRead;
    unsigned char *pchLine = 0;
    double fScale;

    //
    //  Open file and read header
    //
    hFile = CreateFile(szFileName, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING,
                                   FILE_ATTRIBUTE_NORMAL | FILE_FLAG_SEQUENTIAL_SCAN,
                                   NULL);
    if (hFile == INVALID_HANDLE_VALUE)
        goto Error;
	if (ReadFile(hFile, bfType, 2, &nBytesRead, NULL) == 0)
        goto Error;
	if (bfType[0] != 'B' || bfType[1] != 'M')
        goto Error;
    if (ReadFile(hFile, &H, sizeof(H), &nBytesRead, NULL) == 0)
        goto Error;
    nHeight = H.biHeight;
    nWidth = H.biWidth;
    //
    //  Build masks and shifts for 16 and 32-bit modes.
    //
    nOffset = 14 + 40;
    if (H.biCompression == 0) {
        if (H.biBitCount == 16) {
            rgMasks.nRed = 0x7C00;
            rgMasks.nGrn = 0x03E0;
            rgMasks.nBlu = 0x001F;
        } else if (H.biBitCount == 32) {
            rgMasks.nRed = 0x00FF0000;
            rgMasks.nGrn = 0x0000FF00;
            rgMasks.nBlu = 0x000000FF;
        }
    } else if (H.biCompression == 3) {
	    if (ReadFile(hFile, &rgMasks, 12, &nBytesRead, NULL) == 0)
            goto Error;
        nOffset += 12;
    } else
        goto Error;     // only RGB and BITFIELD compression supported
    if (H.biBitCount == 16 || H.biBitCount == 32) {
        for (nRedShift = 0; nRedShift < 32; nRedShift++)
            if (rgMasks.nRed & (1 << nRedShift))
                break;
        for (nGrnShift = 0; nGrnShift < 32; nGrnShift++)
            if (rgMasks.nGrn & (1 << nGrnShift))
                break;
        for (nBluShift = 0; nBluShift < 32; nBluShift++)
            if (rgMasks.nBlu & (1 << nBluShift))
                break;
    }
    //
    //  Read color map, if one exists.  Note that 16, 24 and 32 bit images can
    //  have an optional colormap, recommended in case the image is quantized.
    //
    if (H.biClrUsed) {
        nColors = H.biClrUsed;
    } else if (H.biBitCount == 1)
        nColors = 2;
    else if (H.biBitCount == 4)
        nColors = 16;
    else if (H.biBitCount == 8)
        nColors = 256;
    else if (H.biBitCount >= 16)
        nColors = 0;
    else
        goto Error;
    if (nColors && ReadFile(hFile, rgColors, 4*nColors, &nBytesRead, NULL) == 0)
        goto Error;
    nOffset += 4*nColors;
    if (nColors == 0) {
        if (H.biBitCount >= 16)
            nChannels = 3;    // Color image
        else
            nChannels = 1;    // grey or bitmap
    } else {                    // inspect colormap for non-zero hue
        nChannels = 1 + 2*(H.biBitCount >= 16);
        for (i = 0; i < nColors; i++)
            if (rgColors[i].r != rgColors[i].g || rgColors[i].r != rgColors[i].b) {
                nChannels = 3;    // map contains color
                break;
            }
    }
    //
    //  Read the raster data, convert to linear-luminance floating-point.
    //
	nLineSize = 4*((H.biWidth*H.biBitCount + 31) >> 5);   // pad to 32-bit boundary
    pchLine = new unsigned char[nLineSize];
    if (pchLine == 0)
        goto Error;
    if (NewImage(nWidth, nHeight, nChannels) == 0)
        goto Error;
    if (H.biBitCount == 16)
        fScale = 1.0/31.0;      // 5-5-5
    else
        fScale = 1.0/255.0;     // color map or RGB have 8-bit channels'
    SetGamma(fGamma);
    for (j = m_nHeight - 1; j >= 0; --j) {
	    if (ReadFile(hFile, pchLine, nLineSize, &nBytesRead, NULL) == 0)
            goto Error;
        for (i = 0; i < m_nWidth; i++) {
            switch (H.biBitCount) {

            case 1:
                iRed = (pchLine[i >> 3] >> (7 - (i & 7))) & 1;
                break;
            case 4:
                iRed = (pchLine[i >> 1] >> (4 - 4*(i & 1))) & 15;
                break;
            case 8:
                iRed = pchLine[i];
                break;
            case 16:
                nPixel = pchLine[2*i + 0] | (pchLine[2*i + 1] << 8);
                iBlu = (nPixel & rgMasks.nBlu) >> nBluShift;
                iGrn = (nPixel & rgMasks.nGrn) >> nGrnShift;
                iRed = (nPixel & rgMasks.nRed) >> nRedShift;
                iBlu = (255*iBlu) / (rgMasks.nBlu >> nBluShift);    //REVIEW: Marginal for 555 and 565 RGB
                iGrn = (255*iGrn) / (rgMasks.nGrn >> nGrnShift);
                iRed = (255*iRed) / (rgMasks.nRed >> nRedShift);
                break;
            case 24:
                iRed = pchLine[3*i + 2];    // BGR, not RGB
                iGrn = pchLine[3*i + 1];
                iBlu = pchLine[3*i + 0];
                break;
            case 32:
                nPixel =  pchLine[4*i + 0]        | (pchLine[4*i + 1] <<  8)
                       | (pchLine[4*i + 2] << 16) | (pchLine[4*i + 3] << 24);
                iBlu = (nPixel & rgMasks.nBlu) >> nBluShift;        //REVIEW: Only works for 888 RGB
                iGrn = (nPixel & rgMasks.nGrn) >> nGrnShift;
                iRed = (nPixel & rgMasks.nRed) >> nRedShift;
                break;
            }
            //
            //  Look up color if indexed.  Now iRed is 8-bits.
            //
            if (H.biBitCount <= 8) {
                iBlu = rgColors[iRed].b;
                iGrn = rgColors[iRed].g;
                iRed = rgColors[iRed].r;
            }
            Set(GammaDecode(iRed), i, j, 0);
            if (m_nChannels >= 3) {
                Set(GammaDecode(iGrn), i, j, 1);
                Set(GammaDecode(iBlu), i, j, 2);
            }
        }
    }
    delete [] pchLine;
    CloseHandle(hFile);
    return 1;
Error:
    printf("%s Error return\n", szFileName);
    delete [] pchLine;
    delete [] m_pfImage;
    m_pfImage = 0;
    CloseHandle(hFile);
    return 0;
}

//
//  clip bright in-gamut color to unit RGB cube, preservig luminance if possible.
//
DisplayRGB
ClampToGamut(DisplayRGB &rgb)
{
    float fLum, t;
    DisplayRGB rgbClamped;

    fLum = Luminance(rgb);
    if (fLum < 0.0)
        return DisplayRGB(0.0);
    else if (fLum > 1.0)
        return DisplayRGB(1.0);
    else {
        //
        //  Use clipping to find the intersection of the RGB cube
        //  with the ray from (fLum, fLum, fLum) to (rgb.red, rgb.grn, rgb.blu).
        //
        t = 1.0;
        if (rgb.red < 0.0 && -fLum > (rgb.red - fLum)*t)
            t = -fLum/(rgb.red - fLum);
        if (rgb.red > 1.0 && 1.0 - fLum < (rgb.red - fLum)*t)
            t = (1.0f - fLum)/(rgb.red - fLum);
        if (rgb.grn < 0.0 && -fLum > (rgb.grn - fLum)*t)
            t = -fLum/(rgb.grn - fLum);
        if (rgb.grn > 1.0 && 1.0 - fLum < (rgb.grn - fLum)*t)
            t = (1.0f - fLum)/(rgb.grn - fLum);
        if (rgb.blu < 0.0 && -fLum > (rgb.blu - fLum)*t)
            t = -fLum/(rgb.blu - fLum);
        if (rgb.blu > 1.0 && 1.0 - fLum < (rgb.blu - fLum)*t)
            t = (1.0f - fLum)/(rgb.blu - fLum);
        rgbClamped.red = fLum + t*(rgb.red - fLum);
        rgbClamped.grn = fLum + t*(rgb.grn - fLum);
        rgbClamped.blu = fLum + t*(rgb.blu - fLum);
        return rgbClamped;
    }
}

static HANDLE
WriteHeaderBMP(const char *szFile, int nWidth, int nHeight, int nChannels, int nDPI)
{
	unsigned char	bfType[2];
	struct {
		unsigned int	bfSize;
		unsigned short	bfReserved1;
		unsigned short	bfReserved2;
		unsigned int	bfOffBits;

		unsigned int	biSize;
		unsigned int	biWidth;
		unsigned int	biHeight;
		unsigned short	biPlanes;
		unsigned short	biBitCount;
		unsigned int	biCompression;
		unsigned int	biSizeImage;
		unsigned int	biXPelsPerMeter;
		unsigned int	biYPelsPerMeter;
		unsigned int	biClrUsed;
		unsigned int	biClrImportant;
	} H;
    struct {
        unsigned char   b, g, r, a;
    } rgColors[256];
    unsigned long nBytesWritten, nLineSize;
    int i;
    HANDLE hFile;

    hFile = CreateFile(szFile, GENERIC_WRITE, 0, NULL, CREATE_ALWAYS, 
                                    FILE_ATTRIBUTE_NORMAL, NULL);
    if (hFile == INVALID_HANDLE_VALUE)
        goto Error;

    //
    //  Build header structure
    //
    H.biSize = 40;
    H.biWidth = nWidth;
    H.biHeight = nHeight;
    H.biPlanes = 1;
    if (nChannels <= 2)
        H.biBitCount = 8;       // Grayscale image, needs a colormap
    else
        H.biBitCount = 24;      // RGB color image, or first 3 channels of deeper image
    H.biCompression = 0;
	H.biSizeImage = 0;
	H.biXPelsPerMeter = int(double(nDPI)/0.0254 + 0.5);
	H.biYPelsPerMeter = int(double(nDPI)/0.0254 + 0.5);
	H.biClrUsed = 0;
	H.biClrImportant = 0;
	bfType[0] = 'B';
	bfType[1] = 'M';
	nLineSize = 4*((H.biWidth*H.biBitCount + 31) >> 5);   // pad to 32-bit boundary
	H.bfSize = 40 + 14 + nHeight*nLineSize;
	H.bfReserved1 = 0;
	H.bfReserved2 = 0;
	H.bfOffBits = 40 + 14;
    //
    //  2 channels would be an FFT magnitude/phase image.
    //  The magnitude image will be written
    //
    if (nChannels <= 2) {
        H.bfSize += sizeof(rgColors);
        H.bfOffBits += sizeof(rgColors);
    }
    //
    //  Write header and colormap
    //
	if (WriteFile(hFile, bfType, 2, &nBytesWritten, NULL) == 0)
        goto Error;
    if (WriteFile(hFile, &H, sizeof(H), &nBytesWritten, NULL) == 0)
        goto Error;
    if (nChannels <= 2) {
        for (i = 0; i < 256; i++)
            rgColors[i].r = rgColors[i].g = rgColors[i].b = i;
        if (WriteFile(hFile, rgColors, sizeof(rgColors), &nBytesWritten, NULL) == 0)
            goto Error;
    }
    return hFile;

Error:
    CloseHandle(hFile);
    return INVALID_HANDLE_VALUE;
}

int
Image::WriteBMP(char *szFileName, double fGamma)
{
    HANDLE hFile;
	int i, j, nLineSize;
    unsigned long nBytesWritten;
    unsigned char *pchLine = 0;
    DisplayRGB rgb;
    float *pfRed, *pfGrn, *pfBlu;

    if (m_pfImage == 0)
        return 0;
    SetGamma(fGamma);
    hFile = WriteHeaderBMP(szFileName, m_nWidth, m_nHeight, m_nChannels, 300);
    if (hFile == INVALID_HANDLE_VALUE)
        goto Error;
    if (m_nChannels <= 2)
	    nLineSize = 4*((8 * m_nWidth + 31) >> 5);   // pad to 32-bit boundary
    else
	    nLineSize = 4*((24 * m_nWidth + 31) >> 5);   // pad to 32-bit boundary
    pchLine = new unsigned char[nLineSize];
    if (pchLine == 0)
        goto Error;
    //
    //  Quantize the floating point values to 8-bit channels (grey or RGB).
    //  If color, properly clamp the color to the sRGB gamut.
    //
    for (j = m_nHeight - 1; j >= 0; --j) {
        if (m_nChannels <= 2) {
            for (i = 0; i < m_nWidth; i++) {
                pchLine[i] = GammaEncode(Get(i, j, 0));
            }
        } else {
            pfRed = m_pfImage + m_nWidth*j;
            pfGrn = pfRed + m_nWidth*m_nHeight;
            pfBlu = pfGrn + m_nWidth*m_nHeight;
            for (i = 0; i < m_nWidth; i++) {
                rgb.red = pfRed[i];
                rgb.grn = pfGrn[i];
                rgb.blu = pfBlu[i];
                if (rgb.red < 0.0 || rgb.red > 1.0 
                 || rgb.grn < 0.0 || rgb.grn > 1.0
                 || rgb.blu < 0.0 || rgb.blu > 1.0)
                    rgb = ClampToGamut(rgb);
                pchLine[3*i + 0] = GammaEncode(rgb.blu); // 80 percent of CPU time
                pchLine[3*i + 1] = GammaEncode(rgb.grn);
                pchLine[3*i + 2] = GammaEncode(rgb.red);
            }
        }
	    if (WriteFile(hFile, pchLine, nLineSize, &nBytesWritten, NULL) == 0)
            goto Error;
    }
    delete [] pchLine;
    CloseHandle(hFile);
    return 1;
Error:
    delete [] pchLine;
    CloseHandle(hFile);
    return 0;
}

int
Image::WriteBMP(char *szFileName, char *szFile2, double fGamma)
{
    char sz[256];

    strcpy_s(sz, sizeof(sz), szFileName);
    strcat_s(sz, sizeof(sz), szFile2);
    return WriteBMP(sz, fGamma);
}

int
Image::WriteBMP(char *szFileName, char *szFile2, char *szFile3, double fGamma)
{
    char sz[256];

    strcpy_s(sz, sizeof(sz), szFileName);
    strcat_s(sz, sizeof(sz), szFile2);
    strcat_s(sz, sizeof(sz), szFile3);
    return WriteBMP(sz, fGamma);
}

void
Image::SetGamma(double fGamma)
{
    int n;
    double f;

    if (fGamma == GAMMA_SRGB) {
        for (n = 0; n < 256; n++) {
            f = double(n)*(1.0/255.0);
            if (f <= 0.03928)
                f = f/12.92;
            else
                f = float(pow((f + 0.055)/1.055, 2.4));
            m_rgfGammaT[n] = float(f);
        }
    } else
        for (n = 0; n < 256; n++)
            m_rgfGammaT[n] = float(pow(double(n)*(1.0/255.0), fGamma));
    for (n = 1; n < 256; n++)
        m_rgfGammaA[n] = (m_rgfGammaT[n - 1] + m_rgfGammaT[n])/2;
}

//
//  Read NASA VICAR format.  .LBL files contain the label, .IMG contains the data.
//  8 and 16-bit integer (bytes swapped) and 32-bit VAX floating point are typical formats.
//  Note, some image files contain the label data at the start of that file, so the same
//  name can be given as both LBL and IMG files.
//
static int
ReadLine(HANDLE hFile, char *sz, int nSize)
{
    unsigned long nBytesRead;
    char rgch[1];
    int i;

    for (i = 0; i < nSize-1; i++) {
        if (ReadFile(hFile, rgch, 1, &nBytesRead, NULL) == 0 || nBytesRead == 0) {
            *sz = 0;
            return 0;
        }
        if (rgch[0] == '\n' || rgch[0] == '\r')
            break;
        if (isspace(rgch[0]))
            continue;           // eat all white space
        *sz++ = rgch[0];
    }
    *sz = 0;
    return 1;
}

#define VICAR_UCHAR     1
#define VICAR_USHORT    2
#define VICAR_VAXFLOAT  3
#define VICAR_UNKNOWN   0

static int
VICARlabel(char *szFileName, int &nHeader, int &nWidth, int &nHeight, int &nType,
           int &iChannel, char szName[48], float &fExposure, char szDate[48])
{
    HANDLE hFile;
    char szLine[128], szPixelType[32];
    int nState, nHeadRecSize, nHeadRecords, nPixelBits, nAuxWide, nAuxLong;

    hFile = CreateFile(szFileName, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING,
                                   FILE_ATTRIBUTE_NORMAL | FILE_FLAG_SEQUENTIAL_SCAN,
                                   NULL);
    if (hFile == INVALID_HANDLE_VALUE) {
        return 0;
    }
    nHeadRecSize = nHeadRecords = nPixelBits = 0;
    nAuxWide = nAuxLong = 0;
    nWidth = nHeight = 0;
    iChannel = -1;
    nState = 0;
    while (ReadLine(hFile, szLine, 128)) {
        //
        //  State-machine parser, set states
        //
        if (strcmp(szLine, "END") == 0)
            break;
        if (strcmp(szLine, "OBJECT=IMAGE_HEADER") == 0)     // Looking at image header
            nState = 1;
        if (strcmp(szLine, "OBJECT=FOREIGN_LABEL") == 0)
            nState = 1;
        if (strcmp(szLine, "OBJECT=IMAGE") == 0)            // Looking at image object
            nState = 2;
        if (strlen(szLine) > 11 && strcmp(szLine + strlen(szLine) - 11, "=SFDU_LABEL") == 0)
            nState = 3;                                     // Looking at PDS header
        if (strcmp(szLine, "OBJECT=SPECTROMETER_RECORD") == 0)
            nState = 4;                                     // Fobos-2 auxiliary header
        if (strncmp(szLine, "END_OBJECT", 10) == 0)
            nState = 0;
        //
        //  Look or the various indications of filter # or channel # in
        //  images from spacecrafts.  Also check for a target name.
        //
        if (nState == 3 && strncmp(szLine, "IMAGE_CHANNEL=", 14) == 0)  // e.g., vega 
            iChannel = atoi(szLine + 14);
        if (nState == 3 && strncmp(szLine, "FILTER_NUMBER=", 14) == 0)  // e.g., fobos, galileo
            iChannel = atoi(szLine + 14);
        if (nState == 3 && strncmp(szLine, "TARGET_NAME=", 12) == 0) {
            strncpy_s(szName, 48, szLine+12, 47);
            szName[47] = 0;
        }
        if (nState == 3 && strncmp(szLine, "EXPOSURE_DURATION=", 18) == 0) {
            fExposure = float(atof(szLine + 18));
        }
        //
        //  There is no consistant way to express image time in VICAR/PDS formats
        //
        szDate[0] = 0;
        if (nState == 3 && strncmp(szLine, "IMAGE_TIME=", 11) == 0) {
            strncpy_s(szDate, 48, szLine+11, 47);
            szName[47] = 0;
        }
        if (nState == 3 && strncmp(szLine, "OBSERVATION_TIME=", 17) == 0) {
            strncpy_s(szDate, 48, szLine+17, 47);
            szName[47] = 0;
        }
        //
        //  Look for various object and header entries, for image size
        //
        if (nState == 1 && strncmp(szLine, "BYTES=", 6) == 0)
            nHeadRecSize = atoi(szLine + 6);
        if (nState == 1 && strncmp(szLine, "RECORDS=", 8) == 0)
            nHeadRecords = atoi(szLine + 8);
        if (nState == 2 && strncmp(szLine, "LINES=", 6) == 0)
            nHeight = atoi(szLine + 6);
        if (nState == 2 && strncmp(szLine, "LINE_SAMPLES=", 13) == 0)
            nWidth = atoi(szLine + 13);
        if (nState == 3 && strncmp(szLine, "RECORD_BYTES=", 13) == 0)
            nHeadRecSize = atoi(szLine + 13);
        if (nState == 3 && strncmp(szLine, "LABEL_RECORDS=", 14) == 0)
            nHeadRecords = atoi(szLine + 14);
        if (nState == 4 && strncmp(szLine, "BYTES=", 6) == 0)
            nAuxWide = atoi(szLine + 6);
        if (nState == 4 && strncmp(szLine, "RECORDS=", 8) == 0)
            nAuxLong = atoi(szLine + 8);
        //
        //include the prefix, it is often needed to find sync
        //
        if (nState == 2 && strncmp(szLine, "LINE_PREFIX_BYTES=", 18) == 0)
            nWidth += atoi(szLine + 18);
        if (nState == 2 && strncmp(szLine, "SAMPLE_BITS=", 12) == 0)
            nPixelBits = atoi(szLine + 12);
        if (nState == 2 && strncmp(szLine, "SAMPLE_TYPE=", 12) == 0)
            strcpy_s(szPixelType, 32, szLine + 12);
    }
    if (strncmp(szFileName, "C:\\Don\\Venera15_SAR", 19) == 0 && nWidth == 1050) {
        nWidth = 1052;
        //printf("Fixing error in Venera-15 label file\n");
    }
    CloseHandle(hFile);
    nHeader = nHeadRecSize * nHeadRecords + nAuxWide * nAuxLong;
    //UNSIGNED_INTEGER, MSB_INTEGER (signed) usually all means "unsigned".
    if (nPixelBits == 8)
        nType = VICAR_UCHAR;
    else if (nPixelBits == 16)
        nType = VICAR_USHORT;
    else if (nPixelBits == 32 && strcmp(szPixelType, "VAX_REAL") == 0)
        nType = VICAR_VAXFLOAT;
    else {
        printf("Unknown VICAR pixel format: %d bits, %s\n", nPixelBits, szPixelType);
        nType = VICAR_UNKNOWN;
        return 0;
    }
    if (nWidth && nHeight)
        return 1;
    else
        return 0;
}

//
//  Convert NASA IBM tapes written with VAX floating point binary data
//  Note this entails byte ordering of VAX and also that imposed by IBM 9-track tape
//
static double
ML_VaxFloat(unsigned nVaxFloat)
{
    union {
        unsigned        nVax;
        float           fIEEE;
        unsigned char   rgnBytes[4];
    } u;
    unsigned char nTmp;

    u.nVax = nVaxFloat;
    nTmp = u.rgnBytes[0];
    u.rgnBytes[0] = u.rgnBytes[2];
    u.rgnBytes[2] = nTmp;
    nTmp = u.rgnBytes[3];
    u.rgnBytes[3] = u.rgnBytes[1];
    u.rgnBytes[1] = nTmp;
    return double(u.fIEEE)/4.0;     // avoid underflow by going to double
};

static union {
    unsigned char rgch[40960];
    unsigned int  rgn[10240];
} s_uBuffer;

int
Image::ReadVICAR(char *szName)
{
    char szFileName[256], szTarget[48], szDate[48];
    float fExposure;
    int nWidth, nHeight, nHeader, nType;
    int i, iPixel, iExpected, iChannel;
    unsigned long nBytesRead, nTotalBytes;
    HANDLE hFile;

    strcpy_s(szFileName, 256, szName);
    strcat_s(szFileName, 256, ".lbl");
    iChannel = -1;
    szTarget[0] = 0;
    szDate[0] = 0;
    fExposure = 0.0;
    if (VICARlabel(szFileName, nHeader, nWidth, nHeight, nType,
                   iChannel, szTarget, fExposure, szDate) == 0) {
        //
        //  Header is sometimes written into beginning of .img file
        //
        strcpy_s(szFileName, 256, szName);
        strcat_s(szFileName, 256, ".img");
        if (VICARlabel(szFileName, nHeader, nWidth, nHeight, nType,
                       iChannel, szTarget, fExposure, szDate) == 0)
            return 0;
    }
    if (NewImage(nWidth, nHeight) == 0)
        return 0;
    m_iChannel = iChannel;
    strncpy_s(m_szName, 48, szTarget, 47);
    m_fExposure = fExposure;
    strncpy_s(m_szDate, 48, szDate, 47);
    strcpy_s(szFileName, 256, szName);
    strcat_s(szFileName, 256, ".img");
    hFile = CreateFile(szFileName, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING,
                                   FILE_ATTRIBUTE_NORMAL | FILE_FLAG_SEQUENTIAL_SCAN,
                                   NULL);
    if (hFile == INVALID_HANDLE_VALUE) {
        printf("Cannot open image-data file, %s\n", szFileName);
        return 0;
    }
    if (nHeader) {
        ReadFile(hFile, s_uBuffer.rgch, nHeader, &nBytesRead, NULL);
        //ML_TestAssert(nBytesRead == nHeader);
    }
    iPixel = 0;
    nTotalBytes = 0;
    iExpected = nWidth*nHeight;
    while (ReadFile(hFile, s_uBuffer.rgch, 4096, &nBytesRead, NULL) && nBytesRead > 0) {
        nTotalBytes += nBytesRead;
        switch (nType) {
            case VICAR_UCHAR:
                for (i = 0; i < int(nBytesRead); i++) {
                    if (iPixel < iExpected)
                        m_pfImage[iPixel++] = float(s_uBuffer.rgch[i]) * 1.0f/255.0f;
                    else
                        iPixel++;
                }                    
                break;
            case VICAR_USHORT:
                for (i = 0; i < int(nBytesRead); i += 2) {
                    if (iPixel < iExpected)
                        m_pfImage[iPixel++] = float(s_uBuffer.rgch[i+1] + (s_uBuffer.rgch[i]<<8)) * 1.0f/65535.0f;
                    else
                        iPixel++;
                }                    
                break;
            case VICAR_VAXFLOAT:
                for (i = 0; i < int(nBytesRead); i += 4) {
                    if (iPixel < iExpected)
                        m_pfImage[iPixel++] = float(ML_VaxFloat(s_uBuffer.rgn[i/4]));
                    else
                        iPixel++;
                }
                break;
        }
    }
    //if (iExpected != iPixel)
    //    printf("%d pixels expected, %d read\n", iExpected, iPixel);
    CloseHandle(hFile);
    return 1;
}

static unsigned char s_rgnEBCDIC[256] = {
      0,  1,  2,  3,156,  9,134,127,151,141,142, 11, 12, 13, 14, 15,
     16, 17, 18, 19,157,133,  8,135, 24, 25,146,143, 28, 29, 30, 31,
    128,129,130,131,132, 10, 23, 27,136,137,138,139,140,  5,  6,  7,
    144,145, 22,147,148,149,150,  4,152,153,154,155, 20, 21,158, 26,
     32,160,161,162,163,164,165,166,167,168, 91, 46, 60, 40, 43, 33,
     38,169,170,171,172,173,174,175,176,177, 93, 36, 42, 41, 59, 94,
     45, 47,178,179,180,181,182,183,184,185,124, 44, 37, 95, 62, 63,
    186,187,188,189,190,191,192,193,194, 96, 58, 35, 64, 39, 61, 34,
    195, 97, 98, 99,100,101,102,103,104,105,196,197,198,199,200,201,
    202,106,107,108,109,110,111,112,113,114,203,204,205,206,207,208,
    209,126,115,116,117,118,119,120,121,122,210,211,212,213,214,215,
    216,217,218,219,220,221,222,223,224,225,226,227,228,229,230,231,
    123, 65, 66, 67, 68, 69, 70, 71, 72, 73,232,233,234,235,236,237,
    125, 74, 75, 76, 77, 78, 79, 80, 81, 82,238,239,240,241,242,243,
     92,159, 83, 84, 85, 86, 87, 88, 89, 90,244,245,246,247,248,249,
     48, 49, 50, 51, 52, 53, 54, 55, 56, 57,250,251,252,253,254,255
};
/*
77                   700     968 700 968 L 1                          
SCMVM73 FDS=0059891    IM-1       
 ERT  YR=74 DAY=037 GMT=17/27/46     CCAMERA=A       
 FILTER=6(UV)        NOMINAL EXPOSURE=   65.3 MSEC       LXXXXXXX X       
 SCET  YR=XX DAY=XXX GMT=XX/XX/XX                       
 CCENTER  SL.RANGE=XXXXXXXX KM   KM/LINE=XXX.XX   KM/SAMPLE=XXX.XX
*/
unsigned char g_MarinerLine[968];

int
Image::ReadMariner10(char *szName, char *sz)
{
    unsigned char rgch[968];
    char szHeader[969];
    unsigned long nBytesRead;
    int i, j;
    HANDLE hFile;

    hFile = CreateFile(szName, GENERIC_READ, FILE_SHARE_READ, NULL, OPEN_EXISTING,
                               FILE_ATTRIBUTE_NORMAL | FILE_FLAG_SEQUENTIAL_SCAN,
                               NULL);
    if (hFile == INVALID_HANDLE_VALUE) {
        // printf("Cannot open Mariner-10 file, %s\n", szName);
        return 0;
    }
    ReadFile(hFile, rgch, 968, &nBytesRead, NULL);
    for (i = 0; i < 968; i++) {
        szHeader[i] = s_rgnEBCDIC[rgch[i]];   // VICAR header in old IBM character set
        if (szHeader[i] < 0)
            szHeader[i] = '_';
    }
    szHeader[356] = 0;
    //
    //  Filters are Clear, Minus UV, Orange, UV Polarizing, UV, Blue
    //
    NewImage(832, 700);
    for (i = 0; i < 340; i++) {
        if (strncmp(szHeader + i, "FILTER=", 7) == 0) {
            m_iChannel = atoi(szHeader + i + 7);
            if (m_iChannel == 0)
                m_iChannel = 8;
            break;
        }
    }
    for (i = 0; i < 340; i++) {
        if (strncmp(szHeader + i, "EXPOSURE=", 9) == 0) {
            m_fExposure = float(atof(szHeader + i + 9));
            break;
        }
    }
    for (i = 0; i < 340; i++) {
        if (strncmp(szHeader + i, "FDS=", 4) == 0) {
            m_nFDS = atoi(szHeader + i + 4);
        }
    }
    for (i = 0; i < 340; i++) {
        if (strncmp(szHeader + i, "CAMERA=", 7) == 0) {
            m_nCamera = szHeader[i + 7];
        }
    }
    for (i = 0; i < 340; i++) {
        if (strncmp(szHeader + i, "ERT  YR=", 7) == 0) {
            strncpy_s(m_szDate, 48, szHeader + i + 5, 26);
            break;
        }
    }
    for (j = 0; j < 700; j++) {
        ReadFile(hFile, rgch, 968, &nBytesRead, NULL);
        for (i = 0; i < 832; i++)
            Set(double(rgch[i])/255.0, i, j);
        if (j == 1)
            for (i = 0; i < 968; i++)
                g_MarinerLine[i] = rgch[i];
    }
    strcpy_s(m_szName, 48, "VENUS");
    CloseHandle(hFile);
    if (sz)
        for (i = 0; i < 969; i++)
            sz[i] = szHeader[i];
    return 1;
}

//
//  Kaiser-Windowed Sinc Filter
//
static double gs_fI0Alpha = 1.0/11.3019219521363;
static double gs_fAlphaSquared = 4.0*4.0;
#define F_PI    3.1415926535897932384626433832795028842f

static float
Kaiser4Filter(float x)
{
    float x1, x2, xTerm, fN, fSinc, fKaiser;

    x1 = x * 1.0f/KAISER4_HALFWIDTH;
    x2 = x1*x1;
    if (x2 > 1.0)
        return 0.0;
    x2 = 0.25f * float(gs_fAlphaSquared) * (1.0f - x2);
    fKaiser = 1.0f + x2;
    fN = 2.0f;
    xTerm = x2;
    while (xTerm > 1.0e-4) {
        xTerm *= x2/(fN*fN);
        fN += 1.0;
        fKaiser += xTerm;
        xTerm *= x2/(fN*fN);
        fN += 1.0;
        fKaiser += xTerm;
    }
    x1 = F_PI*x;
    x2 = x1*x1;
    if (x2 < 0.0001f)
        fSinc = 1.0f + x2*(-1.0f/6.0f + x2*1.0f/120.0f);
    else
        fSinc = sinf(x1)/x1;
    return fSinc * fKaiser * float(gs_fI0Alpha); // * 1.0033f;    // correct normalization
}

static float
Kaiser16Filter(float x)
{
    float x1, x2, xTerm, fN, fSinc, fKaiser;

    x1 = x * 1.0f/KAISER16_HALFWIDTH;
    x2 = x1*x1;
    if (x2 > 1.0)
        return 0.0;
    x2 = 0.25f * float(gs_fAlphaSquared) * (1.0f - x2);
    fKaiser = 1.0f + x2;
    fN = 2.0;
    xTerm = x2;
    while (xTerm > 1.0e-4) {
        xTerm *= x2/(fN*fN);
        fN += 1.0;
        fKaiser += xTerm;
        xTerm *= x2/(fN*fN);
        fN += 1.0;
        fKaiser += xTerm;
    }
    x1 = F_PI*x;
    x2 = x1*x1;
    if (x2 < 0.0001f)
        fSinc = 1.0f + x2*(-1.0f/6.0f + x2*1.0f/120.0f);
    else
        fSinc = sinf(x1)/x1;
    return fSinc * fKaiser * float(gs_fI0Alpha); // * 1.0011f;
}

float
Image::Sample(double xSample, double ySample, double fScale, int kChan) const
{
    double x, y, xLow, xHigh, yLow, yHigh, fWeightX, fWeight, fSumWeight, fRecipScale, fValue, fHalfWidth;
    int i, j, iLow, jLow;
    float rgf[512];

    fHalfWidth = KAISER4_HALFWIDTH;
    fSumWeight = fValue = 0.0;
    xLow  = xSample - fHalfWidth*fScale;
    xHigh = xSample + fHalfWidth*fScale;
    iLow = LOWEST_INDEX(xLow);
    yLow  = ySample - fHalfWidth*fScale;
    yHigh = ySample + fHalfWidth*fScale;
    jLow = LOWEST_INDEX(yLow);
    fRecipScale = 1.0/fScale;
    //
    //  Separable filter, store the y-dimension factors
    //
    for (j = jLow, y = INDEX_TO_SAMPLE(jLow); y <= yHigh; j++, y += 1.0)
        rgf[j - jLow] = Kaiser4Filter(float((y - ySample)*fRecipScale));
    rgf[j - jLow] = 0.0;
    for (i = iLow, x = INDEX_TO_SAMPLE(iLow); x <= xHigh; i++, x += 1.0) {
        fWeightX = Kaiser4Filter(float((x - xSample)*fRecipScale));
        for (j = jLow, y = INDEX_TO_SAMPLE(jLow); y <= yHigh; j++, y += 1.0) {
            fWeight = fWeightX * rgf[j - jLow];
            fSumWeight += fWeight;
            fValue += Get(i, j, kChan) * fWeight;
        }
    }
    return float(fValue/fSumWeight);
}

float
Image::UltraSample(double xSample, double ySample, double fScale, int kChan) const
{
    double x, y, xLow, xHigh, yLow, yHigh, fWeightX, fWeight, fSumWeight, fRecipScale, fValue, fHalfWidth;
    int i, j, iLow, jLow;
    float rgf[512];

    fHalfWidth = KAISER16_HALFWIDTH;
    fSumWeight = fValue = 0.0;
    xLow  = xSample - fHalfWidth*fScale;
    xHigh = xSample + fHalfWidth*fScale;
    iLow = LOWEST_INDEX(xLow);
    yLow  = ySample - fHalfWidth*fScale;
    yHigh = ySample + fHalfWidth*fScale;
    jLow = LOWEST_INDEX(yLow);
    fRecipScale = 1.0/fScale;
    //
    //  Separable filter, store the y-dimension factors
    //
    for (j = jLow, y = INDEX_TO_SAMPLE(jLow); y <= yHigh; j++, y += 1.0)
        rgf[j - jLow] = Kaiser16Filter(float((y - ySample)*fRecipScale));
    rgf[j - jLow] = 0.0;
    for (i = iLow, x = INDEX_TO_SAMPLE(iLow); x <= xHigh; i++, x += 1.0) {
        fWeightX = Kaiser16Filter(float((x - xSample)*fRecipScale));
        for (j = jLow, y = INDEX_TO_SAMPLE(jLow); y <= yHigh; j++, y += 1.0) {
            fWeight = fWeightX * rgf[j - jLow];
            fSumWeight += fWeight;
            fValue += Get(i, j, kChan) * fWeight;
        }
    }
    return float(fValue/fSumWeight);
}

const DisplayRGB   ML_Black(0.0, 0.0, 0.0);
const DisplayRGB   ML_White(1.0, 1.0, 1.0);
const DisplayRGB   ML_Red(1.0, 0.0, 0.0);
const DisplayRGB   ML_Green(0.0, 1.0, 0.0);
const DisplayRGB   ML_Blue(0.0, 0.0, 1.0);
const DisplayRGB   ML_RussianRed( 0.415148, 0.004116, 0.004116 );  

DisplayRGB
Image::SampleRGB(double xSample, double ySample, double fScale) const
{
    double x, y, xLow, xHigh, yLow, yHigh, fWeightX, fWeight, fSumWeight, fRecipScale, fHalfWidth;
    int i, j, iLow, jLow;
    DisplayRGB color;
    float rgf[512];

    fHalfWidth = KAISER4_HALFWIDTH;
    fSumWeight = 0.0;
    color = ML_Black;
    xLow  = xSample - fHalfWidth*fScale;
    xHigh = xSample + fHalfWidth*fScale;
    iLow = LOWEST_INDEX(xLow);
    yLow  = ySample - fHalfWidth*fScale;
    yHigh = ySample + fHalfWidth*fScale;
    jLow = LOWEST_INDEX(yLow);
    fRecipScale = 1.0/fScale;
    //
    //  Separable filter, store the y-dimension factors
    //
    for (j = jLow, y = INDEX_TO_SAMPLE(jLow); y <= yHigh; j++, y += 1.0)
        rgf[j - jLow] = Kaiser4Filter(float((y - ySample)*fRecipScale));
    rgf[j - jLow] = 0.0;
    for (i = iLow, x = INDEX_TO_SAMPLE(iLow); x <= xHigh; i++, x += 1.0) {
        fWeightX = Kaiser4Filter(float((x - xSample)*fRecipScale));
        for (j = jLow, y = INDEX_TO_SAMPLE(jLow); y <= yHigh; j++, y += 1.0) {
            fWeight = fWeightX * rgf[j - jLow];
            fSumWeight += fWeight;
            color.red += float(Get(i, j, 0) * fWeight);
            color.grn += float(Get(i, j, 1) * fWeight);
            color.blu += float(Get(i, j, 2) * fWeight);
        }
    }
    return color/fSumWeight;
}

float
Image::Splat(double f, double xSample, double ySample, Image *pimWeight, double fScale, int kChan)
{
    double x, y, xLow, xHigh, yLow, yHigh, fWeightX, fWeight, fRecipScale, fHalfWidth;
    int i, j, iLow, jLow;
    float rgf[512];

    if (fScale <= 0.0)
        return float(f);
    fHalfWidth = KAISER4_HALFWIDTH;
    xLow  = xSample - fHalfWidth*fScale;
    xHigh = xSample + fHalfWidth*fScale;
    iLow = LOWEST_INDEX(xLow);
    yLow  = ySample - fHalfWidth*fScale;
    yHigh = ySample + fHalfWidth*fScale;
    jLow = LOWEST_INDEX(yLow);
    fRecipScale = 1.0/fScale;
    //
    //  Separable filter, store the y-dimension factors
    //
    for (j = jLow, y = INDEX_TO_SAMPLE(jLow); y <= yHigh; j++, y += 1.0)
        rgf[j - jLow] = Kaiser4Filter(float((y - ySample)*fRecipScale));
    rgf[j - jLow] = 0.0;
    for (i = iLow, x = INDEX_TO_SAMPLE(iLow); x <= xHigh; i++, x += 1.0) {
        fWeightX = Kaiser4Filter(float((x - xSample)*fRecipScale));
        for (j = jLow, y = INDEX_TO_SAMPLE(jLow); y <= yHigh; j++, y += 1.0) {
            fWeight = fWeightX * rgf[j - jLow];
            Splat(f*fWeight, i, j, kChan);
            if (pimWeight)
                pimWeight->Splat(fWeight, i, j);
        }
    }
    return float(f);
}

float
Image::UltraSplat(double f, double xSample, double ySample, Image *pimWeight, double fScale, int kChan)
{
    double x, y, xLow, xHigh, yLow, yHigh, fWeightX, fWeight, fRecipScale, fHalfWidth;
    int i, j, iLow, jLow;
    float rgf[512];

    if (fScale <= 0.0)
        return float(f);
    fHalfWidth = KAISER16_HALFWIDTH;
    xLow  = xSample - fHalfWidth*fScale;
    xHigh = xSample + fHalfWidth*fScale;
    iLow = LOWEST_INDEX(xLow);
    yLow  = ySample - fHalfWidth*fScale;
    yHigh = ySample + fHalfWidth*fScale;
    jLow = LOWEST_INDEX(yLow);
    fRecipScale = 1.0/fScale;
    //
    //  Separable filter, store the y-dimension factors
    //
    for (j = jLow, y = INDEX_TO_SAMPLE(jLow); y <= yHigh; j++, y += 1.0)
        rgf[j - jLow] = Kaiser16Filter(float((y - ySample)*fRecipScale));
    rgf[j - jLow] = 0.0;
    for (i = iLow, x = INDEX_TO_SAMPLE(iLow); x <= xHigh; i++, x += 1.0) {
        fWeightX = Kaiser16Filter(float((x - xSample)*fRecipScale));
        for (j = jLow, y = INDEX_TO_SAMPLE(jLow); y <= yHigh; j++, y += 1.0) {
            fWeight = fWeightX * rgf[j - jLow];
            Splat(f*fWeight, i, j, kChan);
            if (pimWeight)
                pimWeight->Splat(fWeight, i, j);
        }
    }
    return float(f);
}

DisplayRGB
Image::SplatRGB(DisplayRGB &rgb, double xSample, double ySample, Image *pimWeight, double fScale)
{
    double x, y, xLow, xHigh, yLow, yHigh, fWeightX, fWeight, fRecipScale, fHalfWidth;
    int i, j, iLow, jLow;
    float rgf[512];

    if (fScale <= 0.0)
        return rgb;
    fHalfWidth = KAISER4_HALFWIDTH;
    xLow  = xSample - fHalfWidth*fScale;
    xHigh = xSample + fHalfWidth*fScale;
    iLow = LOWEST_INDEX(xLow);
    yLow  = ySample - fHalfWidth*fScale;
    yHigh = ySample + fHalfWidth*fScale;
    jLow = LOWEST_INDEX(yLow);
    fRecipScale = 1.0/fScale;
    //
    //  Separable filter, store the y-dimension factors
    //Kaiser4Filter
    for (j = jLow, y = INDEX_TO_SAMPLE(jLow); y <= yHigh; j++, y += 1.0)
        rgf[j - jLow] = Kaiser4Filter(float((y - ySample)*fRecipScale));
    rgf[j - jLow] = 0.0;
    for (i = iLow, x = INDEX_TO_SAMPLE(iLow); x <= xHigh; i++, x += 1.0) {
        fWeightX = Kaiser4Filter(float((x - xSample)*fRecipScale));
        for (j = jLow, y = INDEX_TO_SAMPLE(jLow); y <= yHigh; j++, y += 1.0) {
            fWeight = fWeightX * rgf[j - jLow];
            SplatRGB(rgb*fWeight, i, j);
            if (pimWeight)
                pimWeight->Splat(fWeight, i, j);
        }
    }
    return rgb;
}

int
Image::Normalize(Image *pimWeight)
{
    int k;
    Index64 i, nSize;
    float *pf, *pfW, fW;

    if (pimWeight->m_nWidth != m_nWidth || pimWeight->m_nHeight!= m_nHeight)
        return 0;
    nSize = m_nWidth*Index64(m_nHeight);   // It's safer to use get/set
    for (k = 0; k < m_nChannels; k++) {
        pf = m_pfImage + k*nSize;
        pfW = pimWeight->m_pfImage;
        for (i = 0; i < nSize; i++) {
            fW = *pfW++;
            if (fW)
                *pf++ /= fW;
            else
                pf++;
        }
    }
    return 1;
}

double
Image::SampleLerp(double x, double y) const
{
    unsigned i1, i2, j1, j2;
    double dx, dy, a, b, c, d;

    x -= 0.5;   // correct phase, pixel samples are on N + 0.5.
    y -= 0.5;
    i1 = unsigned(floor(x));
    i2 = i1 + 1;
    dx = x - floor(x);
    j1 = unsigned(floor(y));
    j2 = j1 + 1;
    dy = y - floor(y);
    a = Get(i1, j1);
    b = Get(i2, j1);
    c = Get(i1, j2);
    d = Get(i2, j2);
    a = a + dx*(b - a);
    c = c + dx*(d - c);
    return double(a + dy*(c - a));
}
//
//  Anti-aliased drawing primatives
//
void
Image::DrawLineRGB(DisplayRGB rgb, double xStart, double yStart, double xStop, double yStop, double fWidth)
{
    double dx, dy, fTan, fSec, fSinCos, fTmp, fDist, fLine;
    double x, y, xLine, yLine, xClosest, yClosest, fSupport;
    int iStart, iStop, i, jStart, jStop, j;
    DisplayRGB rgbBack, rgbBlend;

    fSupport = 0.5*(fWidth - 1.0) + 1.0;    // halfwidth of the filtered line
    dx = xStop - xStart;
    dy = yStop - yStart;
    if (dx)
        fTan = dy/dx;
    else
        fTan = 1.0E+100;
    if (dx == 0.0 && dy == 0.0)
        fTan = 0.0;
    if (fabs(fTan) < 1.0) {
        if (xStart > xStop) {
            fTmp = xStop;
            xStop = xStart;
            xStart = fTmp;
            fTmp = yStop;
            yStop = yStart;
            yStart = fTmp;
        }
        fSec = (1.0 + fTan*fTan);
        fSinCos = fTan/fSec;        // Tangent/Secant**2 = Sin*Cos
        fSec = sqrt(fSec);
        iStart = LOWEST_INDEX(xStart - fSupport);
        iStop  = HIGHEST_INDEX(xStop + fSupport);
        for (i = iStart; i <= iStop; i++) {
            x = INDEX_TO_SAMPLE(i);
            yLine = yStart + (x - xStart)*fTan;
            jStart = LOWEST_INDEX(yLine - fSupport*fSec);
            jStop  = HIGHEST_INDEX(yLine + fSupport*fSec);
            for (j = jStart; j <= jStop; j++) {
                y = INDEX_TO_SAMPLE(j);
                rgbBack = GetRGB(i, j);
                fDist = fabs(yLine - y)/fSec;
                xClosest = x - (yLine - y)*fSinCos;
                if (xClosest <= xStart)
                    fDist = sqrt((x - xStart)*(x - xStart) + (y - yStart)*(y - yStart));
                if (xClosest > xStop)
                    fDist = sqrt((x - xStop)*(x - xStop) + (y - yStop)*(y - yStop));
                if (fDist < 0.5*(fWidth - 1.0)) {
                    fLine = 1.0;
                } else if (fDist < 0.5*(fWidth - 1.0) + 1.0) {
                    fLine = 0.5*(fWidth - 1.0) + 1.0 - fDist;
                    //fLine = 0.5*cos(D_PI*(1.0-fLine)) + 0.5;
                } else
                    fLine = 0.0;
                // if (fBack < fLine*f)
                //     Set(fLine*f, i, j, kChan);
                rgbBlend = rgb*fLine + rgbBack*(1.0 - fLine);
                SetRGB(rgbBlend, i, j);
            }
        }
    } else {
        if (yStart > yStop) {
            fTmp = xStop;
            xStop = xStart;
            xStart = fTmp;
            fTmp = yStop;
            yStop = yStart;
            yStart = fTmp;
        }
        dx = xStop - xStart;
        dy = yStop - yStart;
        //fTan = 1.0/fTan;
        fTan = dx/dy;
        fSec = (1.0 + fTan*fTan);
        fSinCos = fTan/fSec;
        fSec = sqrt(fSec);
        jStart = LOWEST_INDEX(yStart - fSupport);
        jStop  = HIGHEST_INDEX(yStop + fSupport);
        for (j = jStart; j <= jStop; j++) {
            y = INDEX_TO_SAMPLE(j);
            xLine = xStart + (y - yStart)*fTan;
            iStart = LOWEST_INDEX(xLine - fSupport*fSec);
            iStop  = HIGHEST_INDEX(xLine + fSupport*fSec);
            for (i = iStart; i <= iStop; i++) {
                x = INDEX_TO_SAMPLE(i);
                rgbBack = GetRGB(i, j);
                fDist = fabs(xLine - x)/fSec;
                yClosest = y - (xLine - x)*fSinCos;
                if (yClosest <= yStart)
                    fDist = sqrt((x - xStart)*(x - xStart) + (y - yStart)*(y - yStart));
                if (yClosest > yStop)
                    fDist = sqrt((x - xStop)*(x - xStop) + (y - yStop)*(y - yStop));
                if (fDist < 0.5*(fWidth - 1.0)) {
                    fLine = 1.0;
                } else if (fDist < 0.5*(fWidth - 1.0) + 1.0) {
                    fLine = 0.5*(fWidth - 1.0) + 1.0 - fDist;
                    //fLine = 0.5*cos(D_PI*(1.0-fLine)) + 0.5;  // raised-cosine impulse looks jaggy
                } else
                    fLine = 0.0;
                // if (fBack < fLine*f)
                //     Set(fLine*f, i, j, kChan);
                rgbBlend = rgb*fLine + rgbBack*(1.0 - fLine);
                SetRGB(rgbBlend, i, j);
            }
        }
 
    }
}

//
//  Draw text with simple vector font
//
#define FONT_UP 0xFE
#define FONT_LAST 0xFF
#define P(x,y)	((((x) & 0xF) << 4) | (((y) & 0xF) << 0))

typedef struct
{
	unsigned char vertex[8]; 
} VectorFont;

static VectorFont rgFont[ML_NFONT] = {
    { FONT_LAST },                                                                                      // ' '
    { P(4,0), P(3,2), P(5,2), P(4,0), FONT_UP, P(4,4), P(4,12), FONT_LAST },                            // '!'
    { P(2,10), P(2,6), FONT_UP, P(6,10), P(6,6), FONT_LAST },                                           // '"'
    { P(0,4), P(8,4), P(6,2), P(6,10), P(8,8), P(0,8), P(2,10), P(2,2) },                               // '#'
    { P(2,2), P(6,4), P(2,8), P(6,10), FONT_UP, P(4,12), P(4,0), FONT_LAST },                           // '$'
    { P(0,1), P(8,11), FONT_UP, P(2,10), P(2,8), FONT_UP, P(6,4), P(6,2) },                             // '%'
    { P(8,0), P(3,12), P(8,9), P(0,4), P(4,0), P(8,4), FONT_LAST },                                     // '&'
    { P(2,6), P(6,10), FONT_LAST },                                                                     // '''
    { P(6,0), P(2,4), P(2,8), P(6,12), FONT_LAST },                                                     // '('
    { P(2,0), P(6,4), P(6,8), P(2,12), FONT_LAST },                                                     // ')'
    { P(1,3), P(4,11), P(7,3), P(0,8), P(8,8), P(1,3), FONT_LAST },                                     // '*'
    { P(1,6), P(7,6), FONT_UP, P(4,9), P(4,3), FONT_LAST },                                             // '+'
    { P(2,0), P(4,2), FONT_LAST },                                                                      // ','
    { P(1,6), P(7,6), FONT_LAST },                                                                      // '-'
    { P(3,0), P(4,0), FONT_LAST },                                                                      // '.'
    { P(0,0), P(8,12), FONT_LAST },                                                                     // '/'
    { P(0,0), P(8,0), P(8,12), P(0,12), P(0,0), P(8,12), FONT_LAST },                                   // '0' to '9'   
    { P(4,0), P(4,12), P(3,10), FONT_LAST },
    { P(0,12), P(8,12), P(8,7), P(0,5), P(0,0), P(8,0), FONT_LAST },
    { P(0,12), P(8,12), P(8,0), P(0,0), FONT_UP, P(0,6), P(8,6), FONT_LAST },
    { P(0,12), P(0,6), P(8,6), FONT_UP, P(8,12), P(8,0), FONT_LAST },
    { P(0,0), P(8,0), P(8,6), P(0,7), P(0,12), P(8,12), FONT_LAST },
    { P(0,12), P(0,0), P(8,0), P(8,5), P(0,7), FONT_LAST },
    { P(0,12), P(8,12), P(8,6), P(4,0), FONT_LAST },
    { P(0,0), P(8,0), P(8,12), P(0,12), P(0,0), FONT_UP, P(0,6), P(8,6), },
    { P(8,0), P(8,12), P(0,12), P(0,7), P(8,5), FONT_LAST },
    { P(4,9), P(4,7), FONT_UP, P(4,5), P(4,3), FONT_LAST },                                             // ':'
    { P(4,9), P(4,7), FONT_UP, P(4,5), P(1,2), FONT_LAST },                                             // ';'
    { P(6,0), P(2,6), P(6,12), FONT_LAST },                                                             // '<'
    { P(1,4), P(7,4), FONT_UP, P(1,8), P(7,8), FONT_LAST },                                             // '='
    { P(2,0), P(6,6), P(2,12), FONT_LAST },                                                             // '>'
    { P(0,8), P(4,12), P(8,8), P(4,4), FONT_UP, P(4,1), P(4,0), FONT_LAST },                            // '?'
    { P(8,4), P(4,0), P(0,4), P(0,8), P(4,12), P(8,8), P(4,4), P(3,6) },                                // '@'
    { P(0,0), P(0,8), P(4,12), P(8,8), P(8,0), FONT_UP, P(0,4), P(8,4) },                               // 'A' to 'Z'
    { P(0,0), P(0,12), P(4,12), P(8,10), P(4,6), P(8,2), P(4,0), P(0,0) },
    { P(8,0), P(0,0), P(0,12), P(8,12), FONT_LAST },
    { P(0,0), P(0,12), P(4,12), P(8,8), P(8,4), P(4,0), P(0,0), FONT_LAST },
    { P(8,0), P(0,0), P(0,12), P(8,12), FONT_UP, P(0,6), P(6,6), FONT_LAST },
    { P(0,0), P(0,12), P(8,12), FONT_UP, P(0,6), P(6,6), FONT_LAST },
    { P(6,6), P(8,4), P(8,0), P(0,0), P(0,12), P(8,12), FONT_LAST },
    { P(0,0), P(0,12), FONT_UP, P(0,6), P(8,6), FONT_UP, P(8,12), P(8,0) },
    { P(0,0), P(8,0), FONT_UP, P(4,0), P(4,12), FONT_UP, P(0,12), P(8,12) },
    { P(0,4), P(4,0), P(8,0), P(8,12), FONT_LAST },
    { P(0,0), P(0,12), FONT_UP, P(8,12), P(0,6), P(6,0), FONT_LAST },
    { P(8,0), P(0,0), P(0,12), FONT_LAST },
    { P(0,0), P(0,12), P(4,8), P(8,12), P(8,0), FONT_LAST },
    { P(0,0), P(0,12), P(8,0), P(8,12), FONT_LAST },
    { P(0,0), P(0,12), P(8,12), P(8,0), P(0,0), FONT_LAST },
    { P(0,0), P(0,12), P(8,12), P(8,6), P(0,5), FONT_LAST },
    { P(0,0), P(0,12), P(8,12), P(8,4), P(0,0), FONT_UP, P(4,4), P(8,0) },
    { P(0,0), P(0,12), P(8,12), P(8,6), P(0,5), FONT_UP, P(4,5), P(8,0) },
    { P(0,2), P(2,0), P(8,0), P(8,5), P(0,7), P(0,12), P(6,12), P(8,10) },
    { P(0,12), P(8,12), FONT_UP, P(4,12), P(4,0), FONT_LAST },
    { P(0,12), P(0,2), P(4,0), P(8,2), P(8,12), FONT_LAST },
    { P(0,12), P(4,0), P(8,12), FONT_LAST },
    { P(0,12), P(2,0), P(4,4), P(6,0), P(8,12), FONT_LAST },
    { P(0,0), P(8,12), FONT_UP, P(0,12), P(8,0), FONT_LAST },
    { P(0,12), P(4,6), P(8,12), FONT_UP, P(4,6), P(4,0), FONT_LAST },
    { P(0,12), P(8,12), P(0,0), P(8,0), FONT_UP, P(2,6), P(6,6), FONT_LAST },
    { P(6,0), P(2,0), P(2,12), P(6,12), FONT_LAST },                                                    // '['
    { P(0,12), P(8,0), FONT_LAST },                                                                     // '\'
    { P(2,0), P(6,0), P(6,12), P(2,12), FONT_LAST },                                                    // ']'
    { P(2,6), P(4,12), P(6,6), FONT_LAST },                                                             // '^'
    { P(0,0), P(8,0), FONT_LAST },                                                                      // '_'
    { P(2,10), P(6,6), FONT_LAST },                                                                     // '`' grave accent
    { P(0,0), P(0,8), P(4,12), P(8,8), P(8,0), FONT_UP, P(0,4), P(8,4) },                               // lower case drawn as 'A' to 'Z'
    { P(0,0), P(0,12), P(4,12), P(8,10), P(4,6), P(8,2), P(4,0), P(0,0) },
    { P(8,0), P(0,0), P(0,12), P(8,12), FONT_LAST },
    { P(0,0), P(0,12), P(4,12), P(8,8), P(8,4), P(4,0), P(0,0), FONT_LAST },
    { P(8,0), P(0,0), P(0,12), P(8,12), FONT_UP, P(0,6), P(6,6), FONT_LAST },
    { P(0,0), P(0,12), P(8,12), FONT_UP, P(0,6), P(6,6), FONT_LAST },
    { P(6,6), P(8,4), P(8,0), P(0,0), P(0,12), P(8,12), FONT_LAST },
    { P(0,0), P(0,12), FONT_UP, P(0,6), P(8,6), FONT_UP, P(8,12), P(8,0) },
    { P(0,0), P(8,0), FONT_UP, P(4,0), P(4,12), FONT_UP, P(0,12), P(8,12) },
    { P(0,4), P(4,0), P(8,0), P(8,12), FONT_LAST },
    { P(0,0), P(0,12), FONT_UP, P(8,12), P(0,6), P(6,0), FONT_LAST },
    { P(8,0), P(0,0), P(0,12), FONT_LAST },
    { P(0,0), P(0,12), P(4,8), P(8,12), P(8,0), FONT_LAST },
    { P(0,0), P(0,12), P(8,0), P(8,12), FONT_LAST },
    { P(0,0), P(0,12), P(8,12), P(8,0), P(0,0), FONT_LAST },
    { P(0,0), P(0,12), P(8,12), P(8,6), P(0,5), FONT_LAST },
    { P(0,0), P(0,12), P(8,12), P(8,4), P(0,0), FONT_UP, P(4,4), P(8,0) },
    { P(0,0), P(0,12), P(8,12), P(8,6), P(0,5), FONT_UP, P(4,5), P(8,0) },
    { P(0,2), P(2,0), P(8,0), P(8,5), P(0,7), P(0,12), P(6,12), P(8,10) },
    { P(0,12), P(8,12), FONT_UP, P(4,12), P(4,0), FONT_LAST },
    { P(0,12), P(0,2), P(4,0), P(8,2), P(8,12), FONT_LAST },
    { P(0,12), P(4,0), P(8,12), FONT_LAST },
    { P(0,12), P(2,0), P(4,4), P(6,0), P(8,12), FONT_LAST },
    { P(0,0), P(8,12), FONT_UP, P(0,12), P(8,0), FONT_LAST },
    { P(0,12), P(4,6), P(8,12), FONT_UP, P(4,6), P(4,0), FONT_LAST },
    { P(0,12), P(8,12), P(0,0), P(8,0), FONT_UP, P(2,6), P(6,6), FONT_LAST },
    { P(6,0), P(4,2), P(4,10), P(6,12), FONT_UP, P(2,6), P(4,6), FONT_LAST },                           // '{'
    { P(4,0), P(4,5), FONT_UP, P(4,6), P(4,12), FONT_LAST },                                            // '|'
    { P(4,0), P(6,2), P(6,10), P(4,12), FONT_UP, P(6,6), P(8,6), FONT_LAST },                           // '}'
    { P(0,4), P(2,8), P(6,4), P(8,8), FONT_LAST },                                                      // '~'

    //{ P(3,12), P(5,11), P(6,8), P(4,6), P(2,6), P(0,8), P(1,11), P(3,12) },                           // degree sign
    { P(2,12), P(4,11), P(4,10), P(3,8), P(1,8), P(0,10), P(0,11), P(2,12) },
    { P(1,6), P(7,6), FONT_UP, P(4,9), P(4,3), FONT_UP, P(1,2), P(7,2) },                               // plus/minus sign
    { P(1,6), P(2,8), P(6,6), P(7,8), FONT_UP, P(1,4), P(7, 4), FONT_LAST },                            // aproximate equal
    { P(1,3), P(7,9), FONT_UP, P(7,3), P(1,9), FONT_LAST },                                             // times sign
    { P(0,2), P(8,2), P(8,10), P(0,10), P(0,2), FONT_LAST },                                            // square
    { P(4,10), P(7,8), P(8,5), P(6,2), P(2,2), P(0,5), P(1,8), P(4,10) },                               // circle
    { P(0,2), P(8,2), P(4,9), P(0,2), FONT_LAST },                                                      // triangle
    { P(0,9), P(8,9), P(4,2), P(0,9), FONT_LAST },                                                      // del
    { P(4,2), P(8,6), P(4,10), P(0,6), P(4,2), FONT_LAST },                                             // diamond
    { P(7,6), P(4,1), P(2,1), P(0,2), P(0,5), P(2, 6), P(4,6), P(7,1) },                                // alpha
    { P(0,0), P(0,10), P(3,12), P(7,10), P(3,6), P(7,2), P(4,0), P(0,2) },                              // beta
};

double
Image::DrawTextRGB(DisplayRGB rgb, double xStart, double yStart, char *sz, double fSize, double fAngle, int nJust)
{
    double dx, dy, x, y, s, c;
    int i, j, iChar, nVertex, nPenUp;
    static double xPen, yPen, xOrg, yOrg, xJust, yJust;

    s = sin(D_PI*fAngle/180.0);
    c = cos(D_PI*fAngle/180.0);
    if (nJust == JUSTIFIED_LEFT) {
        xJust = 0.0;
        yJust = 0.0;
    } else if (nJust == JUSTIFIED_RIGHT) {
        xJust = -c*12*fSize*strlen(sz);
        yJust = -s*12*fSize*strlen(sz);
    } else {
        xJust = -c*12*fSize*strlen(sz)/2.0;
        yJust = -s*12*fSize*strlen(sz)/2.0;
    }
    for (i = 0; sz[i]; i++) {
        xOrg = xPen = xStart + c*12*i*fSize + xJust;      
        yOrg = yPen = yStart + s*12*i*fSize + yJust;
        nPenUp = 1;
        iChar = unsigned char(sz[i]) - ' ';
        if (iChar < 0 || iChar >= ML_NFONT)
            iChar = 0;
        for (j = 0; j < 8; j++) {
            nVertex = rgFont[iChar].vertex[j];
            if (nVertex == FONT_LAST)
                break;
            if (nVertex == FONT_UP) {
                nPenUp = 1;
                continue;
            }
            x = fSize * ((nVertex >> 4) & 0xF);
            y = fSize * ((nVertex >> 0) & 0xF);
            if (sz[i] >= 'a' && sz[i] <= 'z')               // smallcaps
                y *= 0.75;
            dx = c*x - s*y;
            dy = s*x + c*y;
            if (!nPenUp)
                DrawLineRGB(rgb, xPen, m_nHeight - yPen, xOrg + dx, m_nHeight - yOrg - dy, 1.0*fSize);
            xPen = xOrg + dx;
            yPen = yOrg + dy;
		    nPenUp = 0;
        }
    }
    return 12.0*fSize*strlen(sz);           // return width of text in pixels
}
//
//  Test image class
//
static double
RandByte()
{
    static unsigned n = 1;
    double f;

    n = 2654435761 + n * 1099087573;
    f= double(n)/4294967296.0;
    return 0.25 + 0.75*f;
}

void
TestImageClass()
{
    Image im, im2;
    DisplayRGB rgb;
    double fTheta, fRadians, xS, yS, xE, yE, x, y;
    int i, j;
    char sz[2];

    printf("TestImageClass:\n");
    im.NewImage(1280, 720, 3);
    im.FillRGB(ML_Black*0.5);
    rgb = ML_White;
    im.SplatRGB(rgb, 10.0, 10.0);
    for (fTheta = 0.0; fTheta < 360.0; fTheta += 10.0) {
        fRadians = fTheta * D_PI/180.0;
        rgb = DisplayRGB(RandByte(), RandByte(), RandByte());
        xS = 320.0 + 25.0*sin(fRadians);
        yS = 256.0 + 25.0*cos(fRadians);
        xE = 320.0 + 200.0*sin(fRadians);
        yE = 256.0 + 200.0*cos(fRadians);
        im.DrawLineRGB(rgb, xS, yS, xE, yE, 1.0);
    }
    im.DrawTextRGB(ML_White, 960.0, 100.0, "Hello World!", 2.0, 0.0, JUSTIFIED_RIGHT);
    im.DrawTextRGB(ML_White, 960.0, 200.0, "Hello World!", 2.0, 0.0, JUSTIFIED_LEFT);
    im.DrawTextRGB(ML_White, 960.0, 300.0, "Hello World!", 2.0, 0.0, JUSTIFIED_CENTER);
    im.WriteBMP("TestImage.bmp");

    im.ReadBMP("TestImage.bmp");
    rgb = im.GetRGB(100, 100);
    rgb.Print();
    im.NewImage(1280, 1280, 3);
    im.FillRGB(ML_Black);
    for (fTheta = 0.0; fTheta < 360.0; fTheta += 30.0) {
        fRadians = fTheta * D_PI/180.0;
        xS = 640.0 + 120.0*cos(fRadians);
        yS = 640.0 + 120.0*sin(fRadians);
        im.DrawTextRGB(ML_White, xS, yS, "Hello World!", 3.0, fTheta, JUSTIFIED_LEFT);
    }
    im.WriteBMP("TestText2.bmp");

    im.ReadBMP("Artur.bmp");
    for (fTheta = 0.0; fTheta < 360.0; fTheta += 30.0) {
        fRadians = fTheta * D_PI/180.0;
        xS = 640.0 + 120.0*cos(fRadians);
        yS = 640.0 + 120.0*sin(fRadians);
        im.DrawTextRGB(ML_Black, xS, yS, "Hello World!", 3.0, fTheta, JUSTIFIED_LEFT);
    }
    im.WriteBMP("ArturText.bmp");
    im.FillRGB(ML_Black);
    for (i = 1; i < 12; i++) {
        im.DrawTextRGB(ML_White, 100.0, 90.0*i, "Hello World!", double(i)/2.0);
    }
    im.WriteBMP("TestText.bmp");

    im.FillRGB(ML_Black);
    sz[1] = 0;
    for (i = 0; i < ML_NFONT; i++) {
        x = 100.0 + 100.0*(i%10);
        y = 100.0 + 100.0*(i/10);
        sz[0] = ' ' + i;
        im.DrawTextRGB(ML_White, x, y, sz, 5.0);
    }
    im.WriteBMP("VectorFont.bmp");

    im.ReadBMP("Artur.bmp");
    im.m_bWrapX = WRAP_CLAMP;
    im.m_bWrapY = WRAP_CLAMP;
    im2.NewImage(1600, 1600, 3);
    for (j = 0; j < 1600; j++) {
        for (i = 0; i < 1600; i++) {
            im2.SetRGB(im.GetRGB(i, j), i, j);
        }
    }
    im2.WriteBMP("Artur_Clamp.bmp");
}