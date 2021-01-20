//
//  Image class
//  D.P. Mitchell 2016/06/17.
//  D.P. Mitchell 2003/10/12.
//
//  Pixel coordinate convention:    (0,0) --->  i  ---------+
//                                    |                     |
//                                    |                     |
//                                    V                     |
//                                                          |
//                                    j                     |
//                                                          |
//                                    |                     |
//                                    +-----------------(W-1,H-1)
//
#pragma once
#define WRAP_CLAMP          0
#define WRAP_AROUND         1
#define WRAP_MIRROR         2

#define GAMMA_SRGB          -1.0
#define KAISER4_HALFWIDTH   4.0f
#define KAISER16_HALFWIDTH  16.0f

#define SAMPLE_TO_INDEX(X)   int(floor(X))
#define INDEX_TO_SAMPLE(I)   (double(I) + 0.5)
#define LOWEST_INDEX(X)      int(floor(X + 0.5))
#define HIGHEST_INDEX(X)     int(floor(X - 0.5))

#define JUSTIFIED_RIGHT     0
#define JUSTIFIED_LEFT      1
#define JUSTIFIED_CENTER    2

#define TEXT_RIGHT            0.0
#define TEXT_LEFT           180.0
#define TEXT_UP              90.0
#define TEXT_DOWN           270.0

#define ML_ASCII        95
#define ML_SYMBOLS      (ML_ASCII) 
#define ML_NFONT        (ML_ASCII + 11)

#define ML_SYM_DEG      ('~' + 1)
#define ML_SYM_PLUSMIN  ('~' + 2)
#define ML_SYM_TIMES    ('~' + 3)
#define ML_SYM_SQUARE   ('~' + 4)
#define ML_SYM_CIRCLE   ('~' + 5)
#define ML_SYM_TRIANGLE ('~' + 6)

typedef __int64 Index64;
//
//  sRGB color value
//
struct DisplayRGB {
    union {
        float   rgb[3];
        struct {
            float   red;
            float   grn;
            float   blu;
        };
    };

                DisplayRGB() {}
                DisplayRGB(double r, double g, double b) : red(float(r)), grn(float(g)), blu(float(b)) {}
                DisplayRGB(double f) : red(float(f)), grn(float(f)), blu(float(f)) {}
    DisplayRGB  operator +(const DisplayRGB &v) const { return DisplayRGB(red+v.red, grn+v.grn, blu+v.blu); }
    DisplayRGB  operator -(const DisplayRGB &v) const { return DisplayRGB(red-v.red, grn-v.grn, blu-v.blu); }
    DisplayRGB  operator *(double f)            const { return DisplayRGB(red*f, grn*f, blu*f); }
    DisplayRGB  operator /(double f)            const { return DisplayRGB(red/f, grn/f, blu/f); }
    DisplayRGB  operator *(DisplayRGB rgb)      const { return DisplayRGB(red*rgb.red, grn*rgb.grn, blu*rgb.blu); }
    float       MaxRGB() const { if (red > grn && red > blu) return red;
                                 else if (grn > blu) return grn;
                                 else return blu; }

    void        Print()                         const { printf("(%f %f %f)\n", red, grn, blu); }
};

extern DisplayRGB ClampToGamut(DisplayRGB &rgb);    // clamp color to unit RGB cube, preserving luminance if possible.

static inline DisplayRGB
ClampToMax(DisplayRGB &rgb)
{
    if (rgb.red == 0.0 && rgb.grn == 0.0 && rgb.blu == 0.0)
        return rgb;
    else if (rgb.red >= rgb.grn && rgb.red >= rgb.blu)
        return rgb/rgb.red;
    else if (rgb.grn >= rgb.blu)
        return rgb/rgb.grn;
    else
        return rgb/rgb.blu;
}

extern const DisplayRGB    ML_Black;
extern const DisplayRGB    ML_White;
extern const DisplayRGB    ML_Red;
extern const DisplayRGB    ML_Green;
extern const DisplayRGB    ML_Blue;
extern const DisplayRGB    ML_RussianRed;  // sRGB = 171, 21, 21, 

inline float           
Luminance(const DisplayRGB &rgb) 
{ 
    return 0.2126f*rgb.red        // From ML_RGBtoXYZ matrix, Y(D65)
         + 0.7152f*rgb.grn
         + 0.0722f*rgb.blu; 
}

extern DisplayRGB ClampToGamut(DisplayRGB &rgb);
//
//  Mariner-10 image data extensions
//
struct MarinerData {
    int     nFDS;       // flight data subsystem clock
    double  fJD;        // date and time in image header
    int     nCamera;    // camera A or B (0 or 1)
    int     nFilter;    // color filter
    double  fExposure;  // exposure in milliseconds
    double  fJD2;       // date and time from Mariner ephemeris file
    double  x, y, z;    // spacecraft coordinates
    double  c1, c2, c3; // spacecraft popinting angles
};
//
//  Floating-point image with basic .bmp I/O
//
struct Image {
    float   *m_pfImage;
    int     m_nHeight;
    int     m_nWidth;
    int     m_nChannels;
    int     m_bWrapX;
    int     m_bWrapY;
    float   *m_pfBorder;
    float   m_rgfGammaT[256];
    float   m_rgfGammaA[256];

    int     m_iChannel;         // optional indicator of channel or filter number (VICAR and Mariner 10 images)
    int     m_nFDS;             // spacecraft clock
    int     m_nCamera;          // camera 'A' or 'B'
    int     m_nBitErrors;
    float   m_fExposure;
    char    m_szName[48];       // optional name or target
    char    m_szDate[48];


            Image() : m_pfImage(0), m_pfBorder(0) {};
           ~Image() { delete [] m_pfImage; m_pfImage = 0; delete [] m_pfBorder; m_pfBorder = 0; }
    int     NewImage(int nWidth, int nHeight, int nChannels = 1,
                     int bWrapX = WRAP_MIRROR, int bWrapY = WRAP_MIRROR);
    int     NewImage(Image &im) { return NewImage(im.m_nWidth, im.m_nHeight, im.m_nChannels,
                                                im.m_bWrapX, im.m_bWrapY); }
    int     ReadBMP(char *szFileName, double fGamma = GAMMA_SRGB);   // Microsoft image format 
    int     WriteBMP(char *szFileName, double fGamma = GAMMA_SRGB);
    int     WriteBMP(char *szFileName, char *szFile2, double fGamma = GAMMA_SRGB);
    int     WriteBMP(char *szFileName, char *szFile2, char *szFile3, double fGamma = GAMMA_SRGB);
    int     ReadVICAR(char *szName);                                    // NASA image format
    int     ReadMariner10(char *szFileName, char *szHeader = 0);      // NASA format for Mariner-10 images
    int         Wrap(int &iCol, int &jRow) const;
    float       Fill(double f, int kChan = 0);
    float       Get(int iCol, int jRow, int kChan = 0) const;
    float       Set(double f, int iCol, int jRow, int kChan = 0);
    float       Splat(double f, int iCol, int jRow, int kChan = 0);
    DisplayRGB  FillRGB(const DisplayRGB &rgb);
    DisplayRGB  GetRGB(int iCol, int jRow) const;
    DisplayRGB  SetRGB(const DisplayRGB &rgb, int iCol, int jRow);
    DisplayRGB  SplatRGB(DisplayRGB &rgb, int iCol, int jRow);
    float       Sample(double xSample, double ySample, double fScale = 1.0, int kChan = 0) const;
    float       Splat(double f, double xSample, double ySample, Image *pimWeight = 0,
                      double fScale = 1.0, int kChan = 0);
    float       UltraSample(double xSample, double ySample, double fScale = 1.0, int kChan = 0) const;
    float       UltraSplat(double f, double xSample, double ySample, Image *pimWeight = 0,
                      double fScale = 1.0, int kChan = 0);
    DisplayRGB  SampleRGB(double xSample, double ySample, double fScale = 1.0) const;
    DisplayRGB  SplatRGB(DisplayRGB &rgb, double xSample, double ySample, Image *pimWeight = 0, double fScale = 1.0);
    int         Image::Normalize(Image *pimWeight);
    double      SampleLerp(double x, double y) const;
    void        DrawLineRGB(DisplayRGB rgb, double xStart, double yStart, double xStop, double yStop, double fWidth);
    double      DrawTextRGB(DisplayRGB rgb, double xStart, double yStart, char *sz, double fSize, 
                            double fAngle = 0.0, int nJust = JUSTIFIED_LEFT);

    float       GammaDecode(int n)
                {
                    if (n < 0)
                        return 0.0;
                    else if (n > 255)
                        return 1.0;
                    else
                        return m_rgfGammaT[n];
                }
    int         GammaEncode(double f)
                {
                    int n;

                    n = 0;
                    if (f > m_rgfGammaA[128   ]) n = 128;
                    if (f > m_rgfGammaA[n + 64]) n = n + 64;
                    if (f > m_rgfGammaA[n + 32]) n = n + 32;
                    if (f > m_rgfGammaA[n + 16]) n = n + 16;
                    if (f > m_rgfGammaA[n +  8]) n = n +  8;
                    if (f > m_rgfGammaA[n +  4]) n = n +  4;
                    if (f > m_rgfGammaA[n +  2]) n = n +  2;
                    if (f > m_rgfGammaA[n +  1]) n = n +  1;
                    return n;
                }
    void        SetGamma(double fGamma);
};
//
//  Basic pixel get/set/splat
//
inline float
Image::Get(int iCol, int jRow, int kChan) const
{
    int i, j;

    i = iCol;
    j = jRow;
    if (unsigned(i) >= unsigned(m_nWidth) || unsigned(j) >= unsigned(m_nHeight))
        Wrap(i, j);
    return m_pfImage[i + m_nWidth*Index64(j + m_nHeight*kChan)];
}

inline float
Image::Set(double f, int iCol, int jRow, int kChan)
{
    int i, j;

    i = iCol;
    j = jRow;
    if (unsigned(i) >= unsigned(m_nWidth) || unsigned(j) >= unsigned(m_nHeight)) {
        if (Wrap(i, j) == 0)
            return float(f); // wont draw if wrap-around is clamped
    }
    return m_pfImage[i + m_nWidth*Index64(j + m_nHeight*kChan)] = float(f);
}

inline float
Image::Splat(double f, int iCol, int jRow, int kChan)
{
    int i, j;

    i = iCol;
    j = jRow;
    if (unsigned(i) >= unsigned(m_nWidth) || unsigned(j) >= unsigned(m_nHeight)) {
        if (Wrap(i, j) == 0)
            return float(f); // wont draw if wrap-around is clamped
    }
    return m_pfImage[i + m_nWidth*Index64(j + m_nHeight*kChan)] += float(f);
}

inline DisplayRGB
Image::GetRGB(int iCol, int jRow) const
{
    int i, j;

    i = iCol;
    j = jRow;
    if (unsigned(i) >= unsigned(m_nWidth) || unsigned(j) >= unsigned(m_nHeight))
        Wrap(i, j);
    if (m_nChannels < 3)
        return DisplayRGB(m_pfImage[i + m_nWidth*Index64(j)]);
    else
        return DisplayRGB(m_pfImage[i + m_nWidth*Index64(j + m_nHeight*0)],
                          m_pfImage[i + m_nWidth*Index64(j + m_nHeight*1)],
                          m_pfImage[i + m_nWidth*Index64(j + m_nHeight*2)]);
}

inline DisplayRGB
Image::SetRGB(const DisplayRGB &rgb, int iCol, int jRow)
{
    int i, j;

    i = iCol;
    j = jRow;
    if (unsigned(i) >= unsigned(m_nWidth) || unsigned(j) >= unsigned(m_nHeight)) {
        if (Wrap(i, j) == 0)
            return rgb; // wont draw if wrap-around is clamped
    }
    if (m_nChannels < 3)
        return DisplayRGB(m_pfImage[i + m_nWidth*Index64(j)] = float(Luminance(rgb)));
    else {
        m_pfImage[i + m_nWidth*Index64(j + m_nHeight*0)] = rgb.red;
        m_pfImage[i + m_nWidth*Index64(j + m_nHeight*1)] = rgb.grn;
        m_pfImage[i + m_nWidth*Index64(j + m_nHeight*2)] = rgb.blu;
        return rgb;
    }
}

inline DisplayRGB
Image::SplatRGB(DisplayRGB &rgb, int iCol, int jRow)
{
    int i, j;

    i = iCol;
    j = jRow;
    if (unsigned(i) >= unsigned(m_nWidth) || unsigned(j) >= unsigned(m_nHeight)) {
        if (Wrap(i, j) == 0)
            return rgb; // wont draw if wrap-around is clamped
    }
    if (m_nChannels < 3)
        return DisplayRGB(m_pfImage[i + m_nWidth*Index64(j)] += float(Luminance(rgb)));
    else {
        m_pfImage[i + m_nWidth*Index64(j + m_nHeight*0)] += rgb.red;
        m_pfImage[i + m_nWidth*Index64(j + m_nHeight*1)] += rgb.grn;
        m_pfImage[i + m_nWidth*Index64(j + m_nHeight*2)] += rgb.blu;
        return rgb;
    }
}

extern Image im44;